#!/usr/bin/env python3

"""
TDO-scanner (tracrRNA-dependent off-target scanner)
Developed by Xiaoguo Zheng,mail: zhengxiaoguo@sjtu.edu.cn

A tool to scan Genome and Transcriptome (Mature RNA) for specific sequence motifs 
with variable inserts (Type 1) or fixed inserts with backbone mismatches (Type 2).
"""

import argparse
import sys
import re
import regex
import pysam
from collections import defaultdict

# -----------------------------------------------------------------------------
# Utils
# -----------------------------------------------------------------------------

def reverse_complement(seq):
    complement = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
    return seq.translate(complement)[::-1]a

def parse_input_pattern(pattern_str):
    match = re.match(r"([ACGTNacgtn]+)\(([ACGTNacgtn]+)\)([ACGTNacgtn]+)", pattern_str)
    if not match:
        print(f"Error: Pattern format '{pattern_str}' is invalid. Expected format: LEFT(MIDDLE)RIGHT")
        sys.exit(1)
    return match.groups()

def get_upstream_20bp(fasta, chrom, start_0b, end_0b, strand):
    """
    Get 20bp upstream relative to the feature direction.
    """
    try:
        if strand == '+':
            p_start = max(0, start_0b - 20)
            return fasta.fetch(chrom, p_start, start_0b).upper()
        else:
            p_end = end_0b + 20
            seq = fasta.fetch(chrom, end_0b, p_end).upper()
            return reverse_complement(seq)
    except:
        return "N" * 20

def fmt_loc(chrom, start, end):
    """Format 0-based [start, end) to 1-based chr:start-end"""
    return f"{chrom}:{start+1}-{end}"

# -----------------------------------------------------------------------------
# Core Logic
# -----------------------------------------------------------------------------

class TDOscanner:
    def __init__(self, fasta_path, gtf_path, pattern, var_len_range, mismatch_limit):
        self.fasta_path = fasta_path
        self.gtf_path = gtf_path
        self.left, self.mid, self.right = parse_input_pattern(pattern)
        
        try:
            vl = var_len_range.split('-')
            self.min_var = int(vl[0])
            self.max_var = int(vl[1])
            self.mismatches = int(mismatch_limit)
        except ValueError:
            sys.exit("Error: Invalid numeric parameters.")

        # Type 1: Variable Gap
        self.regex_type1 = regex.compile(f"({self.left})([ACGTN]{{{self.min_var},{self.max_var}}})({self.right})", regex.IGNORECASE)
        
        # Type 2: Fixed Backbone with mismatches
        full_ref = self.left + self.mid + self.right
        # (?e) enables error counting. 
        self.regex_type2 = regex.compile(f"(?e)(({full_ref}){{s<={self.mismatches}}})", regex.IGNORECASE)

        # Storage for aggregation
        self.gene_hits_buffer = defaultdict(list)
        self.rna_hits_buffer = defaultdict(list) 

    def load_gtf_data(self):
        print("Loading GTF annotation...")
        self.genes = {}
        self.transcripts = defaultdict(list)
        self.tx_meta = {} 

        with open(self.gtf_path, 'r') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                
                feat = parts[2]
                chrom = parts[0]
                strand = parts[6]
                start = int(parts[3]) - 1
                end = int(parts[4])
                
                attr = parts[8]
                attr_dict = {}
                for x in attr.split(';'):
                    if not x.strip(): continue
                    try:
                        k, v = x.strip().split(' ', 1)
                        attr_dict[k] = v.strip('"')
                    except:
                        pass

                if feat == 'gene':
                    gid = attr_dict.get('gene_id')
                    gname = attr_dict.get('gene_name', gid)
                    gtype = attr_dict.get('gene_biotype') or attr_dict.get('gene_type', 'NA')
                    
                    if gid:
                        self.genes[gid] = {
                            'chr': chrom, 'start': start, 'end': end, 'strand': strand,
                            'name': gname, 'biotype': gtype, 'id': gid
                        }
                
                elif feat == 'exon':
                    tid = attr_dict.get('transcript_id')
                    gid = attr_dict.get('gene_id')
                    if tid:
                        self.transcripts[tid].append((start, end))
                        if tid not in self.tx_meta:
                            tname = attr_dict.get('transcript_name', tid)
                            ttype = attr_dict.get('transcript_biotype') or attr_dict.get('transcript_type', 'NA')
                            gname = attr_dict.get('gene_name', gid)
                            
                            self.tx_meta[tid] = {
                                'chr': chrom, 'strand': strand, 'gid': gid,
                                'tname': tname, 'ttype': ttype, 'gname': gname
                            }

        for tid in self.transcripts:
            self.transcripts[tid].sort()

    def get_rna_genomic_fragments(self, tid, rna_start_idx, rna_end_idx):
        exons = self.transcripts[tid]
        strand = self.tx_meta[tid]['strand']
        
        genomic_coords = []
        if strand == '+':
            for s, e in exons:
                genomic_coords.extend(range(s, e))
        else:
            for s, e in reversed(exons):
                genomic_coords.extend(range(e-1, s-1, -1))
        
        if rna_start_idx >= len(genomic_coords) or rna_end_idx > len(genomic_coords):
             return "Error_Bounds"

        match_g_indices = genomic_coords[rna_start_idx : rna_end_idx]
        match_g_indices.sort()
        
        if not match_g_indices: return ""
        
        blocks = []
        current_start = match_g_indices[0]
        current_prev = match_g_indices[0]
        
        for pos in match_g_indices[1:]:
            if pos == current_prev + 1:
                current_prev = pos
            else:
                blocks.append((current_start, current_prev + 1))
                current_start = pos
                current_prev = pos
        blocks.append((current_start, current_prev + 1))
        
        chrom = self.tx_meta[tid]['chr']
        block_strs = [f"{chrom}:{b_s+1}-{b_e}" for b_s, b_e in blocks]
        
        return ";".join(block_strs)

    def run(self):
        self.load_gtf_data()
        
        files = {
            'gene_t1': open("output_gene_type1.txt", "w"),
            'gene_t2': open("output_gene_type2.txt", "w"),
            'rna_t1': open("output_matureRNA_type1.txt", "w"),
            'rna_t2': open("output_matureRNA_type2.txt", "w")
        }
        
        h_gene_t1 = "Gene_location\tStrand\tGeneID\tGeneName\tGeneBiotype\tMatchedSequence\tTar_location\tVariablePart\tupstream20bp\n"
        h_gene_t2 = "Gene_location\tStrand\tGeneID\tGeneName\tGeneBiotype\tMatchedSequence\tTar_location\tMismatchCount\tupstream20bp\n"
        h_rna_t1 = "Transcript_location\tStrand\tTranscriptID\tTranscriptBiotype\tGeneID\tGeneName\tMatchedSequence\tVariablePart\tRNA_Location\tGenomic_Location\tupstream20bp\n"
        h_rna_t2 = "Transcript_location\tStrand\tTranscriptID\tTranscriptBiotype\tGeneID\tGeneName\tMatchedSequence\tMismatchCount\tRNA_Location\tGenomic_Location\tupstream20bp\n"

        files['gene_t1'].write(h_gene_t1)
        files['gene_t2'].write(h_gene_t2)
        files['rna_t1'].write(h_rna_t1)
        files['rna_t2'].write(h_rna_t2)

        fasta = pysam.FastaFile(self.fasta_path)
        
        print("Scanning Genes...")
        for gid, info in self.genes.items():
            chrom = info['chr']
            strand = info['strand']
            try:
                seq = fasta.fetch(chrom, info['start'], info['end']).upper()
            except: continue
            
            search_seq = reverse_complement(seq) if strand == '-' else seq
            self._scan_gene(search_seq, gid, info, fasta)

        self._flush_gene_buffer(files, fasta)

        print("Scanning Transcripts...")
        for tid, exons in self.transcripts.items():
            info = self.tx_meta[tid]
            chrom = info['chr']
            strand = info['strand']
            
            parts = []
            for s, e in exons:
                try:
                    parts.append(fasta.fetch(chrom, s, e).upper())
                except:
                    parts.append("N" * (e-s))
            
            full_seq = "".join(parts)
            if strand == '-': full_seq = reverse_complement(full_seq)
            
            self._scan_rna(full_seq, tid, info, files, fasta)

        self._flush_rna_buffer(files)

        for f in files.values(): f.close()
        fasta.close()
        print("Done.")

    def _is_core_intact(self, matched_seq):
        """Check if middle part matches self.mid strictly"""
        len_left = len(self.left)
        len_mid = len(self.mid)
        # Assuming substitution only (s<=2), length is preserved.
        if len(matched_seq) != (len_left + len_mid + len(self.right)):
             return False
        core_seq = matched_seq[len_left : len_left + len_mid]
        return core_seq.upper() == self.mid.upper()

    def _scan_gene(self, seq, gid, info, fasta):
        chrom = info['chr']
        strand = info['strand']
        g_start_offset = info['start']

        def get_genomic_coords(s_idx, e_idx):
            if strand == '+':
                g_s = g_start_offset + s_idx
                g_e = g_start_offset + e_idx 
            else:
                g_s = info['end'] - e_idx
                g_e = info['end'] - s_idx 
            return g_s, g_e

        # Type 1 (Regex already handles consume logic fine as overlaps are rare/handled by distinct pattern)
        for m in self.regex_type1.finditer(seq):
            g_s, g_e = get_genomic_coords(m.start(), m.end())
            up20 = get_upstream_20bp(fasta, chrom, g_s, g_e, strand)
            var_part = m.group(2)
            match_seq = m.group(0)
            key = (chrom, g_s, g_e, strand, match_seq, var_part, 0, up20)
            self.gene_hits_buffer["t1"].append({'loc_key': key, 'info': info})

        # Type 2 - ADDED overlapped=True
        # This matches Perl's window-sliding behavior
        for m in self.regex_type2.finditer(seq, overlapped=True):
            match_seq = m.group(0)
            
            if not self._is_core_intact(match_seq):
                continue

            g_s, g_e = get_genomic_coords(m.start(), m.end())
            up20 = get_upstream_20bp(fasta, chrom, g_s, g_e, strand)
            cnt = sum(m.fuzzy_counts)
            key = (chrom, g_s, g_e, strand, match_seq, "", cnt, up20)
            self.gene_hits_buffer["t2"].append({'loc_key': key, 'info': info})

    def _flush_gene_buffer(self, files, fasta):
        for mtype in ["t1", "t2"]:
            aggregated = defaultdict(list)
            for item in self.gene_hits_buffer[mtype]:
                aggregated[item['loc_key']].append(item['info'])
            
            for key, info_list in aggregated.items():
                chrom, g_s, g_e, strand, match_seq, var_part, mm_count, up20 = key
                
                g_locs = ";".join([fmt_loc(i['chr'], i['start'], i['end']) for i in info_list])
                ids = ";".join([i['id'] for i in info_list])
                names = ";".join([i['name'] for i in info_list])
                types = ";".join([i['biotype'] for i in info_list])
                match_loc = fmt_loc(chrom, g_s, g_e)
                
                if mtype == "t1":
                    line = f"{g_locs}\t{strand}\t{ids}\t{names}\t{types}\t{match_seq}\t{match_loc}\t{var_part}\t{up20}\n"
                    files['gene_t1'].write(line)
                else:
                    line = f"{g_locs}\t{strand}\t{ids}\t{names}\t{types}\t{match_seq}\t{match_loc}\t{mm_count}\t{up20}\n"
                    files['gene_t2'].write(line)

    def _scan_rna(self, seq, tid, info, files, fasta):
        chrom = info['chr']
        strand = info['strand']
        
        exons = self.transcripts[tid]
        min_s = exons[0][0]
        max_e = exons[-1][1]
        tx_loc_str = fmt_loc(chrom, min_s, max_e)

        def buffer_rna_hit(match, mtype, var_part="", mm_count=0):
            s_idx = match.start()
            e_idx = match.end()
            match_seq = match.group(0)
            
            rna_loc_str = f"{s_idx+1}-{e_idx}"
            genomic_loc_str = self.get_rna_genomic_fragments(tid, s_idx, e_idx)
            
            up20 = ""
            if strand == '+':
                first_base_loc_str = self.get_rna_genomic_fragments(tid, s_idx, s_idx+1)
                try:
                    fb_part = first_base_loc_str.split(':')[1].split('-')[0]
                    g_s_point = int(fb_part) - 1
                    up20 = get_upstream_20bp(fasta, chrom, g_s_point, g_s_point, strand)
                except:
                    up20 = "N"*20
            else:
                first_base_loc_str = self.get_rna_genomic_fragments(tid, s_idx, s_idx+1)
                try:
                    fb_part = first_base_loc_str.split(':')[1].split('-')[1]
                    g_e_point = int(fb_part)
                    up20 = get_upstream_20bp(fasta, chrom, 0, g_e_point, strand)
                except:
                    up20 = "N"*20
            
            loc_key = (genomic_loc_str, match_seq, var_part, mm_count, up20, strand)
            
            details = {
                'tx_loc': tx_loc_str,
                'tid': tid,
                'ttype': info['ttype'],
                'gid': info['gid'],
                'gname': info['gname'],
                'rna_loc': rna_loc_str
            }
            
            self.rna_hits_buffer[mtype].append({'key': loc_key, 'val': details})

        # Scan T1
        for m in self.regex_type1.finditer(seq):
            buffer_rna_hit(m, "t1", var_part=m.group(2))
            
        # Scan T2 - ADDED overlapped=True
        for m in self.regex_type2.finditer(seq, overlapped=True):
            match_seq = m.group(0)
            
            if not self._is_core_intact(match_seq):
                continue
            
            cnt = sum(m.fuzzy_counts)
            buffer_rna_hit(m, "t2", mm_count=cnt)

    def _flush_rna_buffer(self, files):
        print("Aggregating RNA results...")
        for mtype in ["t1", "t2"]:
            grouped = defaultdict(list)
            for item in self.rna_hits_buffer[mtype]:
                grouped[item['key']].append(item['val'])
            
            for key, details_list in grouped.items():
                genomic_loc_str, match_seq, var_part, mm_count, up20, strand = key
                
                tx_locs = ";".join([d['tx_loc'] for d in details_list])
                tids    = ";".join([d['tid'] for d in details_list])
                ttypes  = ";".join([d['ttype'] for d in details_list])
                gids    = ";".join([d['gid'] for d in details_list])
                gnames  = ";".join([d['gname'] for d in details_list])
                r_locs  = ";".join([d['rna_loc'] for d in details_list])
                
                base = f"{tx_locs}\t{strand}\t{tids}\t{ttypes}\t{gids}\t{gnames}"
                
                if mtype == "t1":
                    line = f"{base}\t{match_seq}\t{var_part}\t{r_locs}\t{genomic_loc_str}\t{up20}\n"
                    files['rna_t1'].write(line)
                else:
                    line = f"{base}\t{match_seq}\t{mm_count}\t{r_locs}\t{genomic_loc_str}\t{up20}\n"
                    files['rna_t2'].write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TDO scanner with overlapped regex search")
    parser.add_argument("fasta", help="Genomic FASTA file")
    parser.add_argument("gtf", help="Annotation GTF file")
    parser.add_argument("pattern", help="Format: LEFT(MIDDLE)RIGHT, e.g., ACG(T)GCA")
    parser.add_argument("range", help="Variable length range, e.g., 1-5")
    parser.add_argument("mismatch", help="Max mismatches allowed for Type 2")
    
    args = parser.parse_args()
    
    scanner = TDOscanner(args.fasta, args.gtf, args.pattern, args.range, args.mismatch)
    scanner.run()
