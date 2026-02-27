# TDOscanner
## TracrRNA-Dependent Off-target Scanner

A specialized bioinformatics tool for scanning genome and transcriptome sequences to identify specific DNA/RNA motifs with variable-length spacers (Type 1) or fixed-spacer motifs with backbone mismatches (Type 2).

# Features
### Dual scanning modes:
  - Type 1: Variable-length insertion between conserved flanking sequences
  - Type 2: Fixed middle sequence with allowed mismatches in flanking regions
### Dual-level analysis:
  - Gene-level scanning (genomic DNA)
  - Mature RNA-level scanning (spliced transcripts)

## Installation
### System Requirements
- **Operating System:** Linux (Ubuntu 20.04+, CentOS 7+, or compatible distributions)
- **Memory:** 8 GB RAM minimum, 16+ GB recommended for large genomes
- **Storage:** 50+ GB free disk space for reference genomes

## Create analysis environment
```
conda create -n TDOscanner -c bioconda -c conda-forge python=3.11 conda-forge::regex conda-forge::bzip2 bioconda::pysam bioconda::htslib bioconda::gffread -y
conda activate TDOscanner
```

## Quick Start
```
git clone https://github.com/Xiaoguo-Zheng/TDOscanner.git
cd tdo-scanner
python tdo_scanner.py --help  
python get_candidate_offtarget.py --help   
```
## Usage
### Basic Command
`python tdo_scanner.py <genome.fa> <annotation.gtf> <pattern> <range> <mismatch>`
`python get_candidate_offtarget.py <genome.fa> `

### Arguments
|Argument|Description|Example|
|:-------------:|:--------------------------------------:|:-------------------:|
|**genome.fa**|	Reference genome in FASTA format	|hg38.fa|
|**annotation.gtf**	|Gene annotation in GTF format	|hg38.gtf|
|**pattern**|	Sequence pattern in LEFT(MIDDLE)RIGHT format|	GTTTA(GA)GCTA|
|**range**	|Variable length range for Type 1 (min-max)|	2-6|
|**mismatch**|	Maximum allowed mismatches for Type 2|	2|

##Download reference 
#Homo_sapiens reference  
```
#hg38.fa
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#hg38.gtf
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.113.chr.gtf.gz
```  

#Mus_musculus reference  
```
#mm39.fa
wget https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

#mm39.gtf
wget https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.chr.gtf.gz
gunzip Mus_musculus.GRCm39.113.chr.gtf.gz
```


## Example
```
python tdo_scanner.py hg38.fa hg38.gtf "GTTTTA(GA)GCTA" "2-6" 2
python get_candidate_offtarget.py hg38.fa
```

## Input Pattern Format
### The pattern follows the format: LEFT(MIDDLE)RIGHT
  - LEFT: Conserved 5' flanking sequence
  - MIDDLE: Core sequence (fixed for Type 2, variable for Type 1)
  - RIGHT: Conserved 3' flanking sequence
### Pattern Examples
 - GTTTA(GA)GCTA: LEFT=GTTTA, MIDDLE=GA, RIGHT=GCTA
  - ACGT(NNNN)TTAA: LEFT=ACGT, MIDDLE=NNNN, RIGHT=TTAA
### Output Files
  - The tool generates four output files:

1.output_gene_type1.txt - Type 1 matches at gene level  
2.output_gene_type2.txt - Type 2 matches at gene level  
3.output_matureRNA_type1.txt - Type 1 matches at spliced RNA level  
4.output_matureRNA_type2.txt - Type 2 matches at spliced RNA level  

## Output Columns  
### Gene-level output:
```
#type1 output
GeneID	GeneName	GeneBiotype	MatchedSequence	Tar_location	VariablePart	upstream20bp	MM0_Count	MM1_Count	MM0_Locations	MM1_Locations

#type2 output
Gene_location	Strand	GeneID	GeneName	GeneBiotype	MatchedSequence	Tar_location	MismatchCount	upstream20bp	MM0_Count	MM1_Count	MM0_Locations	MM1_Locations

```
### RNA-level output:
```
#type1 output
Transcript_location	Strand	TranscriptID	TranscriptBiotype	GeneID	GeneName	MatchedSequence	VariablePart	RNA_Location	Genomic_Location	upstream20bp	MM0_Count	MM1_Count	MM0_Locations	MM1_Locations

#type2 output
Transcript_location	Strand	TranscriptID	TranscriptBiotype	GeneID	GeneName	MatchedSequence	MismatchCount	RNA_Location	Genomic_Location	upstream20bp	MM0_Count	MM1_Count	MM0_Locations	MM1_Locations

```
## Algorithm Details
### Type 1 Scanning
- Uses regex pattern: LEFT + [ACGTN]{min,max} + RIGHT
- Matches sequences where the middle region can be any sequence of length within specified range
- Requires exact matches in LEFT and RIGHT regions
### Type 2 Scanning
- Uses fuzzy regex matching with substitution tolerance
- Requires exact match in MIDDLE region
- Allows up to N mismatches in LEFT+RIGHT regions combined
- Uses regex library's fuzzy matching capabilities
### Sequence Processing
1. Gene sequences: Extracted based on GTF gene coordinates
2. RNA sequences: Constructed by concatenating exons from transcripts
3. Reverse strand handling: Automatically reverse-complemented for scanning
4. Upstream sequence: 20bp immediately upstream of match start (genomic coordinates)
## Use Cases
### 1. CRISPR Off-target Detection

- Find potential off-target sites with variable spacers
- Identify sites with partial homology to gRNA sequences
### 2. Regulatory Element Discovery

- Scan for transcription factor binding sites with flexible spacing
- Identify RNA-binding protein targets in transcripts
### 3.Sequence Motif Mining
- Discover novel sequence patterns in genomes/transcriptomes
- Identify conserved motifs with variable intervening sequences
## Limitations
1. Performance: Whole-genome scanning can be time-consuming for large genomes
2. Memory Usage: Loading large FASTA files requires sufficient RAM
3. GTF Compatibility: Requires standard GTF format with gene and exon features
4. Pattern Complexity: Complex patterns with many degenerate bases may reduce specificity

## Troubleshooting
### Common Issues
1. "Error: Pattern format is invalid"

- Ensure pattern follows LEFT(MIDDLE)RIGHT format
- Check for special characters or spaces
2. "KeyError" when fetching sequences

- Verify chromosome names match between FASTA and GTF files
- Check that GTF coordinates are within FASTA sequence boundaries
3. No matches found

- Verify pattern case (tool is case-insensitive but uses uppercase internally)
- Adjust range/mismatch parameters for broader matching
- Check that sequences are being extracted correctly

## Citation
If you use TDO-scanner in your research, please cite:

TDOscanner: A tool for scanning genome and transcriptome for specific sequence motifs.


## License
MIT license

## Support
For questions, bug reports, or feature requests:

- Email: zhengxiaoguo@sjtu.edu.cn

## Version History
v1.0.0 (Current): Initial release with dual scanning modes
Planned features: Parallel processing, BED format input, web interface






