# TRUST4 Benchmark Setup Documentation

## Date
2026-09-18

## TRUST4 Version & Build Information

- **Version**: v1.1.10 (stable)
- **Commit**: 4032f90 (merged PR #410 from liulab-dfci/dev)
- **Repository**: https://github.com/liulab-dfci/TRUST4
- **Paper**: Song, L., Cohen, D., Ouyang, Z. et al. TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data. Nat Methods (2021). https://doi.org/10.1038/s41592-021-01142-2
- **Installation Location**: `/data1/wetlab/TRUST4_benchmark/algorithms/TRUST4/`
- **Build Date**: 2026-09-18
- **Build Status**: Successfully compiled; binaries ready

## Installation & Build

### Build Command
```bash
cd /data1/wetlab/TRUST4_benchmark/algorithms/TRUST4
make
```

### Build Requirements
- gcc (present: `/usr/bin/gcc`)
- make (present: `/usr/bin/make`)
- pthreads (included with system)
- zlib (installed via conda-forge: v1.3.2)
- samtools (installed via bioconda: v1.24)
- perl (conda env trust4_benchmark: v5.38.x)

### Build Outputs
- `run-trust4` - Perl wrapper script (entry point)
- `trust4` - Main C++ binary
- `fastq-extractor` - Read extraction tool
- `bam-extractor` - BAM file read extraction tool
- `annotator` - Gene/CDR annotation tool

## Reference Data

### Human IMGT Reference File
- **Location**: `/data1/wetlab/TRUST4_benchmark/reference/TRUST4/human_IMGT+C.fa`
- **File Size**: 575 KB (compressed from 6.0MB IMGT database download)
- **Number of Gene Sequences**: 1,489
- **Content**: Complete IMGT reference sequences for human BCR/TCR V, D, J, and C genes (constant regions)
- **Source**: IMGT (International ImMunoGeneTics Information System)
- **Download URL**: https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP
- **Generation Method**: 
  ```bash
  perl /data1/wetlab/TRUST4_benchmark/algorithms/TRUST4/BuildImgtAnnot.pl Homo_sapiens > human_IMGT+C.fa
  ```
- **Description**: This file contains:
  - IGH (immunoglobulin heavy chain) genes: V, D, J, and C regions
  - IGK (kappa light chain) genes: V, J, and C regions
  - IGL (lambda light chain) genes: V, J, and C regions
  - This is the recommended reference for FASTQ input

## Test Run Setup: FZ-116 Sample

### Sample Information
- **Sample Name**: FZ-116
- **SRA Accession**: SRR7882936
- **Project**: PRJNA492301
- **Data Type**: Bulk RNA-seq
- **Matched Validation Data**: FZ-116.csv.gz (iRepertoire BCR-seq reference)

### Input Files
- **Read 1**: `/data1/wetlab/TRUST4_benchmark/data/fastq/SRR7882936_1.fastq.gz` (4.4 GB)
- **Read 2**: `/data1/wetlab/TRUST4_benchmark/data/fastq/SRR7882936_2.fastq.gz` (5.3 GB)
- **Total Input**: ~9.7 GB gzipped paired-end RNA-seq reads

### TRUST4 Command
```bash
cd /data1/wetlab/TRUST4_benchmark/results/TRUST4/FZ-116

/data1/wetlab/TRUST4_benchmark/algorithms/TRUST4/run-trust4 \
    -f /data1/wetlab/TRUST4_benchmark/reference/TRUST4/human_IMGT+C.fa \
    --ref /data1/wetlab/TRUST4_benchmark/reference/TRUST4/human_IMGT+C.fa \
    -1 /data1/wetlab/TRUST4_benchmark/data/fastq/SRR7882936_1.fastq.gz \
    -2 /data1/wetlab/TRUST4_benchmark/data/fastq/SRR7882936_2.fastq.gz \
    -t 8 \
    -o TRUST_FZ-116
```

### Command Line Options Explained
- `-f`: Path to FASTA file with V/D/J/C gene genomic coordinates. For FASTQ input, this is the IMGT reference file.
- `--ref`: Detailed reference file from IMGT (recommended for better annotation). Also the IMGT reference.
- `-1`: Path to Read 1 FASTQ file (gzipped)
- `-2`: Path to Read 2 FASTQ file (gzipped)
- `-t 8`: Number of threads (8 threads, matching the published benchmark)
- `-o`: Output prefix for all output files

### Thread Configuration
- **Threads Used**: 8
- **Rationale**: The published TRUST4 benchmark paper (Song et al. 2021, Nat Methods) used 8 threads for their evaluations. Using the same configuration ensures consistent and comparable results.

### Runtime Environment
- **Conda Environment**: `trust4_benchmark`
- **Python Version**: N/A (Perl-based pipeline)
- **Perl**: 5.38.x (from conda)
- **samtools**: v1.24 (from bioconda)

### Job Submission
- **Execution Method**: Direct bash invocation
- **Start Time**: 2026-09-18 14:23:00
- **Script Location**: `/data1/wetlab/TRUST4_benchmark/scripts/run_trust4_test.sh`

## TRUST4 Pipeline Stages

TRUST4 processes data through four sequential stages:

1. **Stage 0 (Extraction)**: Extracts candidate BCR/TCR-bearing reads from RNA-seq FASTQ files
   - Tool: `fastq-extractor`
   - Input: Raw FASTQ files
   - Output: `TRUST_FZ-116_toassemble_1.fq`, `TRUST_FZ-116_toassemble_2.fq`

2. **Stage 1 (Assembly)**: De novo assembly of extracted reads into contigs
   - Assembles reads into consensus contigs for V, D, J, and C genes
   - Constructs full-length BCR/TCR sequences
   - Output: `TRUST_FZ-116_raw.out`, `TRUST_FZ-116_final.out`

3. **Stage 2 (Annotation)**: Aligns contigs to IMGT reference genes
   - Tool: `annotator`
   - Identifies V, D, J, C gene identities
   - Defines CDR1, CDR2, CDR3 boundaries
   - Output: `TRUST_FZ-116_annot.fa`, `TRUST_FZ-116_cdr3.out`

4. **Stage 3 (Report Generation)**: Generates summary reports
   - Creates AIRR-format output
   - Generates CD3-focused report
   - Output: `TRUST_FZ-116_report.tsv`, `TRUST_FZ-116_airr.tsv`

## Expected Output Files

### Primary Output Files

1. **trust_raw.out** (`TRUST_FZ-116_raw.out`)
   - Raw contigs from assembly
   - Fields: contig_id, contig_sequence, nucleotide_abundance_weight
   - Purpose: All assembled consensus sequences before filtering

2. **trust_final.out** (`TRUST_FZ-116_final.out`)
   - Filtered final contigs with read coverage information
   - Fields: contig_id, contig_sequence, read_support_count
   - Purpose: High-confidence assembled sequences

3. **trust_annot.fa** (`TRUST_FZ-116_annot.fa`)
   - FASTA format with detailed annotations in header
   - Header format: `consensus_id consensus_length average_coverage [annotations]`
   - Annotations include: V/D/J/C gene calls, CDR1/2/3 positions and sequences
   - Purpose: Full annotated sequences with gene assignments

4. **trust_cdr3.out** (`TRUST_FZ-116_cdr3.out`)
   - TSV file with per-consensus annotation details
   - Columns:
     - `consensus_id`: Unique contig identifier
     - `index_within_consensus`: Position within consensus (for multi-chain contigs)
     - `V_gene`: Called V gene name(s) (up to 3 ranked by similarity)
     - `D_gene`: Called D gene (heavy chain only)
     - `J_gene`: Called J gene name(s)
     - `C_gene`: Called constant region gene
     - `CDR1`: CDR1 nucleotide sequence
     - `CDR2`: CDR2 nucleotide sequence
     - `CDR3`: CDR3 nucleotide sequence (most important for BCR diversity)
     - `CDR3_score`: Score 1.00 = complete, 0.01 = imputed, others = motif signal strength
     - `read_fragment_count`: Number of reads supporting this consensus
     - `CDR3_germline_similarity`: Alignment similarity to germline sequence
     - `complete_vdj_assembly`: Boolean indicating if full VDJ sequence detected
   - Purpose: Detailed annotation for clonotype analysis

5. **trust_report.tsv** (`TRUST_FZ-116_report.tsv`)
   - Simplified TSV report focused on CDR3 clonotypes
   - Columns:
     - `read_count`: Total reads assigned to this clonotype
     - `frequency`: Proportion within BCR/TCR chain type
     - `CDR3_dna`: CDR3 nucleotide sequence
     - `CDR3_amino_acids`: CDR3 amino acid sequence (for comparison with repertoire databases)
     - `V`: V gene assignment
     - `D`: D gene assignment
     - `J`: J gene assignment
     - `C`: Constant region
     - `consensus_id`: Link to full consensus sequence
     - `consensus_id_complete_vdj`: Link to complete VDJ sequence
   - Purpose: AIRR-compatible clonotype summary for downstream analysis

6. **trust_airr.tsv** (`TRUST_FZ-116_airr.tsv`)
   - Full AIRR (Adaptive Immune Receptor Repertoire) format output
   - Follows: https://docs.airr-community.org/en/latest/datarep/rearrangements.html
   - Purpose: Standard format for sharing with other repertoire analysis tools

### Intermediate Files (Generated During Processing)

- `TRUST_FZ-116_toassemble_1.fq` - Extracted candidate reads (Read 1)
- `TRUST_FZ-116_toassemble_2.fq` - Extracted candidate reads (Read 2)
- Assembly working files (variable k-mer assemblies, alignment results)

### Log Files

- `/data1/wetlab/TRUST4_benchmark/logs/TRUST4_FZ-116.log` - Full TRUST4 execution log

## Output Directory Structure

```
/data1/wetlab/TRUST4_benchmark/results/TRUST4/FZ-116/
├── TRUST_FZ-116_raw.out              # Raw assembly contigs
├── TRUST_FZ-116_final.out            # Final contigs with coverage
├── TRUST_FZ-116_annot.fa             # Annotated sequences
├── TRUST_FZ-116_cdr3.out             # Detailed annotations (CDR3-focused)
├── TRUST_FZ-116_report.tsv           # Simplified clonotype report
├── TRUST_FZ-116_airr.tsv             # AIRR-format output
├── TRUST_FZ-116_toassemble_1.fq      # Extracted reads (R1)
└── TRUST_FZ-116_toassemble_2.fq      # Extracted reads (R2)
```

## Key Points on Output Format

### CDR3 Scoring
- In `trust_cdr3.out`, CDR3_score values indicate:
  - **1.00**: Complete CDR3 with reference-guided imputation (highest confidence)
  - **0.01 to 0.99**: CDR3 with partial motif evidence
  - **0.00**: Partial/uncertain CDR3

### Gene Scoring
- Each gene assignment includes similarity rank (1-3), with coordinates mapped to consensus and reference positions
- Format: `GENE_NAME(ref_length):(cons_start-cons_end):(ref_start-ref_length):similarity_percent`

### Frequency Normalization
- In `trust_report.tsv`, frequencies are normalized separately for:
  - BCR chains (IGH, IGK, IGL)
  - TCR chains (TRA, TRB, TRG, TRD)
  - Sum to 100% within each chain type

### AIRR Compatibility
- All outputs are compatible with AIRR analysis ecosystems
- `trust_airr.tsv` is directly compatible with TCRMatch and other IEDB tools
- Includes provenance (reads, consensus_id, gene identities, CDR3 sequence)

## Validation & Comparison

### Matched iRepertoire Reference
For validation, FZ-116 has matched BCR-seq data:
- **File**: `/data1/wetlab/TRUST4_benchmark/reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz`
- **Format**: CSV with columns: `CDR3(pep),V,D,J,C,CDR3(nuc),copy`
- **Content**: ~100-1000 high-confidence BCR clonotypes from targeted BCR-seq
- **Purpose**: Gold-standard validation; TRUST4 should reconstruct these clonotypes from bulk RNA-seq

## Troubleshooting & Notes

### Known Issues
- TRUST4 requires sufficient disk space during assembly (up to 10x input size for temporary files)
- For very large datasets (>100M reads), extraction can take 1-2 hours on 8 threads
- Extraction phase is the most time-consuming step

### Dependencies in trust4_benchmark Environment
```
samtools v1.24
zlib v1.3.2
perl v5.38.x
gcc (system)
make (system)
```

### Conda Activation
Always activate the environment before running TRUST4:
```bash
conda activate trust4_benchmark
```

## Next Steps

1. **Monitor FZ-116 completion** - Check `/data1/wetlab/TRUST4_benchmark/logs/TRUST4_FZ-116.log` for progress
2. **Inspect output format** - Analyze `TRUST_FZ-116_cdr3.out` and `TRUST_FZ-116_report.tsv` structures
3. **Validate against iRepertoire** - Compare TRUST4 CDR3 calls to FZ-116.csv.gz ground truth
4. **Profile reconstruction accuracy** - Calculate sensitivity and precision for CDR3 detection
5. **Benchmark other samples** - Once validated on FZ-116, run remaining 5 samples
6. **Prepare for TRUST4 modification** - Document baseline performance before implementing graph-based improvements

## References

- Song, L., Cohen, D., Ouyang, Z. et al. (2021). TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data. Nat Methods 18, 627–630. https://doi.org/10.1038/s41592-021-01142-2
- TRUST4 GitHub: https://github.com/liulab-dfci/TRUST4
- IMGT Database: https://www.imgt.org/
- AIRR Data Standards: https://docs.airr-community.org/
