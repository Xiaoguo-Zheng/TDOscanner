import pysam
import argparse
import math
import statistics
from collections import defaultdict, Counter

# ==========================================
# Core Function: Official Window Clustering
# ==========================================
def build_cluster(chrom, window_reads):
    """
    Aggregate data to match official guideseq.py calculations.
    """
    positions = [item[0] for item in window_reads]
    cluster_start = min(positions)
    cluster_end = max(positions)
    
    plus_mi = 0
    minus_mi = 0
    plus_total = 0
    minus_total = 0
    umi_details_counter = Counter()
    
    for pos, record in window_reads:
        umi_details_counter[record["umi_seq"]] += record["count"]
        if record["strand"] == "+":
            plus_mi += 1                   # Official .mi is the number of unique BAM lines/molecules
            plus_total += record["count"]  # Official .total is the sum of raw PCR counts
        else:
            minus_mi += 1
            minus_total += record["count"]
            
    pos_counts = Counter(positions)
    peak_pos = pos_counts.most_common(1)[0][0]
    
    pos_stdev = statistics.stdev(positions) if len(positions) > 1 else 0.0
    
    bi_sum_mi = plus_mi + minus_mi
    bi_geom_mean = math.sqrt(plus_mi * minus_mi)
    total_sum = plus_total + minus_total
    
    return {
        "chrom": chrom,
        "start": cluster_start,
        "end": cluster_end,
        "peak": peak_pos,
        "plus_mi": plus_mi,
        "minus_mi": minus_mi,
        "bi_sum_mi": bi_sum_mi,
        "bi_geom_mean": bi_geom_mean,
        "plus_total": plus_total,
        "minus_total": minus_total,
        "total_sum": total_sum,
        "pos_stdev": pos_stdev,
        "umi_details_counter": umi_details_counter
    }

def analyze_bam(bam_file, output_file, fasta_file=None, window_size=25, min_mapq=20):
    raw_sites = defaultdict(lambda: defaultdict(list))

    print(f"[INFO] Reading BAM file: {bam_file} ...")
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            # Filter unmapped and low quality reads
            if read.is_unmapped or read.mapping_quality < min_mapq:
                continue

            chrom = read.reference_name
            is_rev = read.is_reverse
            pos = read.reference_start if not is_rev else read.reference_end
            if pos is None:
                continue

            strand = "-" if is_rev else "+"
            
            # Extract Official Count and UMI from read name
            qname = read.query_name
            qname_parts = qname.split('_')
            if len(qname_parts) >= 2:
                try:
                    count = int(qname_parts[-1])
                    umi_seq = "_".join(qname_parts[1:-1])
                except ValueError:
                    count = 1
                    umi_seq = "NoUMI"
            else:
                count = 1
                umi_seq = "NoUMI"
            
            raw_sites[chrom][pos].append({
                "umi_seq": umi_seq,
                "strand": strand,
                "count": count
            })

    print(f"[INFO] Performing official fixed-width clustering (window <= {window_size}bp) ...")
    clusters =[]
    for chrom in sorted(raw_sites.keys()):
        chrom_reads =[]
        for pos in sorted(raw_sites[chrom].keys()):
            for record in raw_sites[chrom][pos]:
                chrom_reads.append((pos, record))
        
        if not chrom_reads:
            continue
            
        curr_window = [chrom_reads[0]]
        
        for item in chrom_reads[1:]:
            pos, record = item
            first_pos = curr_window[0][0]
            # Official logic: distance from FIRST read to current read must be <= window_size
            if pos - first_pos <= window_size:
                curr_window.append(item)
            else:
                clusters.append(build_cluster(chrom, curr_window))
                curr_window = [item]
                
        if curr_window:
            clusters.append(build_cluster(chrom, curr_window))

    fasta = None
    if fasta_file:
        print(f"[INFO] Loading Reference Genome: {fasta_file} ...")
        fasta = pysam.FastaFile(fasta_file)
        
    print("[INFO] Generating official-style dsODN integration report with Strict Filters ...")

    with open(output_file, 'w') as out:
        # Construct exact official headers (Total_Reads removed as requested)
        header =[
            "BED_Chromosome", "BED_Min.Position", "BED_Max.Position", "BED_Name",
            "Position", "WindowSequence",
            "+.mi", "-.mi", "bi.sum.mi", "bi.geometric_mean.mi",
            "+.total", "-.total", "total.sum",
            "position.stdev",
            "Filter_Pass",
            "Reads_per_UMI_Detail"
        ]
        
        out.write("\t".join(header) + "\n")
        
        for cl in clusters:
            # Core filter: Only output sites with total UMI >= 2
            if cl["bi_sum_mi"] < 2:
                continue
                
            chrom = cl["chrom"]
            c_start = cl["start"]
            c_end = cl["end"]
            peak_pos = cl["peak"]
            
            # Match official BED_Name format
            bed_name = f"{chrom}_{peak_pos}_{cl['bi_sum_mi']}"
            
            window_seq = ""
            
            # If FASTA is provided, extract 50bp sequence around the peak as WindowSequence
            if fasta:
                search_start = max(0, peak_pos - 25)
                search_end = peak_pos + 25
                try:
                    window_seq = fasta.fetch(chrom, search_start, search_end).upper()
                except KeyError:
                    pass
            
            # [Core Modification]: Define Filter logic
            total_reads = cl["total_sum"]
            
            # Check conditions: +.mi >= 1 AND -.mi >= 1 AND total_reads >= 3 AND bi.geom_mean >= 2
            # Note: bi.geom_mean >= 2 mathematically requires (+.mi * -.mi) >= 4
            if cl["plus_mi"] >= 1 and cl["minus_mi"] >= 1 and total_reads >= 3 and cl["bi_geom_mean"] >= 2.0:
                filter_pass = "Y"
            else:
                filter_pass = "N"
            
            umi_details = ", ".join([f"{u}:{c}" for u, c in cl["umi_details_counter"].most_common()])
            
            # Assemble row with the new columns (Total_Reads removed)
            row_data =[
                chrom, str(c_start), str(c_end), bed_name,
                str(peak_pos), window_seq,
                str(cl["plus_mi"]), str(cl["minus_mi"]), str(cl["bi_sum_mi"]), f"{cl['bi_geom_mean']:.2f}",
                str(cl["plus_total"]), str(cl["minus_total"]), str(cl["total_sum"]),
                f"{cl['pos_stdev']:.2f}",
                filter_pass,
                umi_details
            ]
            
            out.write("\t".join(row_data) + "\n")

    print(f"[SUCCESS] Analysis complete! Extracted {len(clusters)} dsODN sites. Results saved to: {output_file}")
    if fasta:
        fasta.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Genome-wide dsODN Integration Scanner (Official Format)")
    parser.add_argument("-i", "--input_bam", required=True, help="Input aligned BAM/SAM file")
    parser.add_argument("-o", "--output_txt", required=True, help="Output TSV report file")
    
    # Optional arguments
    parser.add_argument("-f", "--fasta", type=str, help="Reference genome FASTA file (Optional, required for WindowSequence)")
    parser.add_argument("-w", "--window_size", type=int, default=25, help="Maximum span of a cluster (bp), default: 25")
    parser.add_argument("-q", "--min_mapq", type=int, default=20, help="Minimum mapping quality, default: 20")
    
    args = parser.parse_args()
    
    analyze_bam(
        bam_file=args.input_bam, 
        output_file=args.output_txt, 
        fasta_file=args.fasta, 
        window_size=args.window_size, 
        min_mapq=args.min_mapq
    )
