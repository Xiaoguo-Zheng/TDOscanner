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
import regex  # Requires check: pip install regex
import pysam
from collections import defaultdict

# -----------------------------------------------------------------------------
# Helper Classes and Functions
# -----------------------------------------------------------------------------

def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def parse_input_pattern(pattern_str):
    """
    Parses the pattern string like "GTTTA(GA)GCTA".
    Returns: (Left_seq, Middle_seq, Right_seq)
    """
    match = re.match(r"([ACGTNacgtn]+)\(([ACGTNacgtn]+)\)([ACGTNacgtn]+)", pattern_str)
    if not match:
        print(f"Error: Pattern format '{pattern_str}' is invalid. Expected format: LEFT(MIDDLE)RIGHT")
        sys.exit(1)
    return match.groups()

def get_upstream_seq(fasta, chrom, start_0b, strand, length=20):
    """
    Retrieves the 20bp upstream sequence based on strand.
    start_0b: 0-based start index of the feature.
    """

    try:
        if strand == '+':
            up_start = max(0, start_0b - length)
            seq = fasta.fetch(chrom, up_start, start_0b)
            return seq.upper()
        else:
            
            up_start = start_0b
            up_end = start_0b + length
            seq = fasta.fetch(chrom, up_start, up_end)
            return reverse_complement(seq.upper())
    except Exception as e:
        return "N" * length

# -----------------------------------------------------------------------------
# Core Scanner Logic
# -----------------------------------------------------------------------------

class TDOscanner:
    def __init__(self, fasta_path, gtf_path, pattern, var_len_range, mismatch_limit):
        self.fasta_path = fasta_path
        self.gtf_path = gtf_path
        
        # Parse Pattern
        self.left, self.mid, self.right = parse_input_pattern(pattern)
        
        # Parse params
        try:
            vl = var_len_range.split('-')
            self.min_var = int(vl[0])
            self.max_var = int(vl[1])
            self.mismatches = int(mismatch_limit)
        except ValueError:
            print("Error: Invalid numeric parameters.")
            sys.exit(1)

        # Prepare Regex Objects
        # Type 1: Left + Variable(min,max) + Right (Exact backbone matches)
        # Note: regex library allows variable length lookbehinds/aheads, but simple composition is easier.
        self.regex_type1 = regex.compile(f"({self.left})([ACGTN]{{{self.min_var},{self.max_var}}})({self.right})", regex.IGNORECASE)
        
        # Type 2: Full sequence with mismatches, but strictly check middle later
        self.full_ref_seq = self.left + self.mid + self.right
        # Use fuzzy regex allowing substitution only (s<=X)
        self.regex_type2 = regex.compile(f"(?e)(({self.full_ref_seq}){{s<={self.mismatches}}})", regex.IGNORECASE)

        # Output handlers
        self.out_files = {}
        
    def open_outputs(self):
        self.out_files['gene_t1'] = open("output_gene_type1.txt", "w")
        self.out_files['gene_t2'] = open("output_gene_type2.txt", "w")
        self.out_files['rna_t1'] = open("output_matureRNA_type1.txt", "w")
        self.out_files['rna_t2'] = open("output_matureRNA_type2.txt", "w")
        
        # Write Headers
        header_base = "GeneID\tChr\tStart\tEnd\tStrand\tUpstream20bp\tMatchSeq"
        self.out_files['gene_t1'].write(f"{header_base}\tVariablePart\n")
        self.out_files['gene_t2'].write(f"{header_base}\tMismatchCount\n")
        
        header_rna = "TranscriptID\tGenomicChr\tGenomicStart\tGenomicEnd\tStrand\tTramscriptPos\tUpstream20bp\tMatchSeq"
        self.out_files['rna_t1'].write(f"{header_rna}\tVariablePart\n")
        self.out_files['rna_t2'].write(f"{header_rna}\tMismatchCount\n")

    def close_outputs(self):
        for f in self.out_files.values():
            f.close()

    def load_gtf_data(self):
        """
        Parses GTF to extract Gene coordinates and Exon structures.
        Only keeps 'gene' and 'exon' features.
        """
        print("Loading GTF annotation...")
        self.genes = {} # gene_id -> {chrom, start, end, strand}
        self.transcripts = defaultdict(list) # transcript_id -> list of (start, end) tuples
        self.tx_to_gene = {} # transcript_id -> gene_id
        self.tx_info = {} # transcript_id -> {chrom, strand}

        with open(self.gtf_path, 'r') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                
                feat_type = parts[2]
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]

                # Simple attribute parser
                attr_dict = {}
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if not attr: continue
                    if " " in attr:
                        key, val = attr.split(' ', 1)
                        val = val.replace('"', '')
                        attr_dict[key] = val

                if feat_type == 'gene':
                    gid = attr_dict.get('gene_id')
                    if gid:
                        self.genes[gid] = {'chr': chrom, 'start': start, 'end': end, 'strand': strand}
                
                elif feat_type == 'exon':
                    tid = attr_dict.get('transcript_id')
                    gid = attr_dict.get('gene_id')
                    if tid:
                        # Store essential info
                        self.transcripts[tid].append((start, end)) # 1-based GTF coords
                        self.tx_to_gene[tid] = gid
                        self.tx_info[tid] = {'chr': chrom, 'strand': strand}

        # Sort exons for transcripts
        for tid in self.transcripts:
            self.transcripts[tid].sort()

    def run(self):
        self.open_outputs()
        self.load_gtf_data()
        
        # Open Fasta
        with pysam.FastaFile(self.fasta_path) as fasta:
            
            # --- Analysis 1: GENE Level ---
            print("Scanning Gene sequences...")
            for gene_id, info in self.genes.items():
                chrom = info['chr']
                start = info['start'] - 1 # Convert to 0-based
                end = info['end'] # pysam end is exclusive
                strand = info['strand']
                
                try:
                    seq_dna = fasta.fetch(chrom, start, end).upper()
                except KeyError:
                    # Contig in GTF not in Fasta
                    continue

                if strand == '-':
                    seq_search = reverse_complement(seq_dna)
                else:
                    seq_search = seq_dna

                self.scan_and_write(seq_search, gene_id, chrom, start, end, strand, fasta, mode="gene")

            # --- Analysis 2: Mature RNA Level ---
            print("Scanning Mature RNA sequences...")
            for tid, exons in self.transcripts.items():
                chrom = self.tx_info[tid]['chr']
                strand = self.tx_info[tid]['strand']
                
                # Construct Mature RNA
                rna_seq = ""
                exon_seqs = []
                
                # Exons are 1-based start, inclusive end
                coords_map = [] # To map RNA index back to Genome index
                
                temp_seq = ""
                coords_list = [] # List of genomic coords for each base in RNA
                
                # Extract sequence for raw exons
                raw_exons_seq = []
                for (es, ee) in exons:
                    try:
                        # 0-based start, ee is exclusive for pysam logic
                        # es is 1-based in GTF -> es-1 for 0-based
                        seq = fasta.fetch(chrom, es-1, ee).upper()
                        raw_exons_seq.append(seq)
                        
                        # Generate genomic coordinates for every base
                        for pos in range(es, ee + 1):
                            coords_list.append(pos)
                    except:
                        pass
                
                if not raw_exons_seq: continue

                full_dnaseq = "".join(raw_exons_seq)
                
                # Handle Strand
                if strand == '-':
                    rna_seq = reverse_complement(full_dnaseq)
                    # Coordinates also need to be reversed to match the 5'->3' RNA sequence
                    coords_list = coords_list[::-1]
                else:
                    rna_seq = full_dnaseq
                
                self.scan_and_write(rna_seq, tid, chrom, 0, 0, strand, fasta, mode="rna", rna_coords_map=coords_list)

        self.close_outputs()
        print("Done. Check output files.")

    def scan_and_write(self, sequence, obj_id, chrom, g_start_0b, g_end, strand, fasta_handle, mode="gene", rna_coords_map=None):
        """
        Generic function to scan a sequence and write results.
        """
        
        # --- Type 1 Scan ---
        for match in self.regex_type1.finditer(sequence):
            # match.groups() -> (Left, VarPart, Right)
            full_match_seq = match.group(0)
            var_part = match.group(2)
            
            start_idx = match.start()
            end_idx = match.end()
            
            self.write_hit(mode=mode, type_idx="t1", 
                           obj_id=obj_id, chrom=chrom, strand=strand, 
                           seq_start_idx=start_idx, seq_end_idx=end_idx, 
                           match_seq=full_match_seq, extra_info=var_part, 
                           fasta=fasta_handle, 
                           gene_g_start=g_start_0b, gene_g_end=g_end,
                           rna_coords_map=rna_coords_map)

        # --- Type 2 Scan ---
        for match in self.regex_type2.finditer(sequence):
            full_match_seq = match.group(0)
            
            
            l_len = len(self.left)
            m_len = len(self.mid)
            r_len = len(self.right)
            
            
            if self.mid not in full_match_seq:
                # If the exact middle sequence isn't found, it's definitely not a Type 2 match
                # (since Middle must be Fixed)
                continue
            
            # Calculate mismatch count
            errors = sum(match.fuzzy_counts)
            
            self.write_hit(mode=mode, type_idx="t2", 
                           obj_id=obj_id, chrom=chrom, strand=strand, 
                           seq_start_idx=match.start(), seq_end_idx=match.end(), 
                           match_seq=full_match_seq, extra_info=str(errors), 
                           fasta=fasta_handle, 
                           gene_g_start=g_start_0b, gene_g_end=g_end,
                           rna_coords_map=rna_coords_map)

    def write_hit(self, mode, type_idx, obj_id, chrom, strand, 
                  seq_start_idx, seq_end_idx, match_seq, extra_info, 
                  fasta, gene_g_start=0, gene_g_end=0, rna_coords_map=None):
        
        # Calculate Genomic Coordinates of the match start/end
        genomic_start = 0
        genomic_end = 0
        upstream20 = ""
        
        if mode == "gene":
            # seq follows strand.
            # gene_g_start is the 0-based lower numeric coordinate of the gene block on chr.
            
            if strand == '+':
                # Sequence is 5'->3' identical to forward genomic
                # relative start + absolute start
                abs_start = gene_g_start + seq_start_idx
                abs_end = gene_g_start + seq_end_idx
                
                genomic_start = abs_start + 1 # 1-based output? standard bioinformatics is 1-based often
                genomic_end = abs_end
                
                upstream20 = get_upstream_seq(fasta, chrom, abs_start, strand, 20)
                
            else: # Strand '-'

                # Length of gene seq
                gene_len = gene_g_end - gene_g_start
                
                
                # Start (higher coord on genome)
                g_high = gene_g_end - seq_start_idx 
                # End (lower coord on genome)
                g_low = gene_g_end - seq_end_idx
                
                genomic_start = g_low + 1 # Convert to likely 1-based for output convention
                genomic_end = g_high
                
                # Upstream for '-' strand is having higher coordinates
                # We need the base physically following g_high
                # Pass the coordinate corresponding to the 5' end of the motif (which is g_high)
                upstream20 = get_upstream_seq(fasta, chrom, g_high, strand, 20)

        elif mode == "rna":
            # Map RNA index to genomic coords
            try:
                # rna_coords_map is a list of genomic positions (1-based usually if derived from GTF directly, let's assume 1-based)
                # But in load_gtf, I stored int(parts[3]).
                
                # Get the genomic coordinate of the FIRST base of the match
                g_start_val = rna_coords_map[seq_start_idx]
                g_end_val = rna_coords_map[seq_end_idx - 1] # inclusive last base
                
                genomic_start = g_start_val
                genomic_end = g_end_val
                
                if strand == '+':
                    # Upstream is lower coordinate. g_start_val is 1-based.
                    # Convert to 0-based for pysam fetch
                    up_ref = g_start_val - 1
                    upstream20 = get_upstream_seq(fasta, chrom, up_ref, strand, 20)
                else:

                    
                    upstream20 = get_upstream_seq(fasta, chrom, g_start_val, strand, 20)
                    
            except IndexError:
                genomic_start = "BoundErr"
                genomic_end = "BoundErr"
                upstream20 = "NNN"

        # Write to file
        file_key = f"{mode}_{type_idx}"
        outfile = self.out_files[file_key]
        
        if mode == "gene":
            # columns: GeneID, Chr, LoopStart, LoopEnd, LoopStrand, Upstream, MatchSeq, Extra
            outfile.write(f"{obj_id}\t{chrom}\t{genomic_start}\t{genomic_end}\t{strand}\t{upstream20}\t{match_seq}\t{extra_info}\n")
        else:
            # columns: TransID, Chr, LoopStart, LoopEnd, LoopStrand, RNA_Index, Upstream, MatchSeq, Extra
            rna_pos_str = f"{seq_start_idx+1}-{seq_end_idx}"
            outfile.write(f"{obj_id}\t{chrom}\t{genomic_start}\t{genomic_end}\t{strand}\t{rna_pos_str}\t{upstream20}\t{match_seq}\t{extra_info}\n")

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TDOscanner: Search gene/RNA motifs. Developed by LazyMan.")
    parser.add_argument("fasta", help="Path to reference genome (hg38.fa)")
    parser.add_argument("gtf", help="Path to annotation file (hg38.gtf)")
    parser.add_argument("pattern", help="Sequence pattern, e.g., GTTTA(GA)GCTA")
    parser.add_argument("range", help="Variable length range for Type 1, e.g., 2-6")
    parser.add_argument("mismatch", help="Max allowed mismatches for Type 2, e.g., 2")
    
    args = parser.parse_args()
    
    scanner = TDOscanner(args.fasta, args.gtf, args.pattern, args.range, args.mismatch)
    scanner.run()
