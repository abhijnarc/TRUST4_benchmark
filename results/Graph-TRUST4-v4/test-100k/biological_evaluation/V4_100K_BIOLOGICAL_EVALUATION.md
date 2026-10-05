# Graph-TRUST4-v4 100K Biological Evaluation

## Objective
Determine whether branch-aware, quality-aware graph assembly improves recovery of biologically annotatable BCR sequence relative to the prior graph implementation.

## Input
- 100K paired-end subset; completed v4 contig FASTA.
- FZ-116 iRepertoire reference.
- TRUST4 baseline uses full FZ-116.

## Annotation yield
| metric | all contigs | % all | IGH | % IGH |
|---|---:|---:|---:|---:|
| total_contigs | 131032 | 100.000 | 712 | 100.000 |
| IGH_contigs | 712 | 0.543 | 712 | 100.000 |
| V_assigned | 35384 | 27.004 | 712 | 100.000 |
| J_assigned | 50086 | 38.224 | 0 | 0.000 |
| C_assigned | 107065 | 81.709 | 0 | 0.000 |
| CDR3_assigned | 22966 | 17.527 | 1 | 0.140 |
| CDR3_amino_acid_assigned | 22966 | 17.527 | 1 | 0.140 |
| V+J | 20962 | 15.998 | 0 | 0.000 |
| V+J+CDR3 | 20928 | 15.972 | 0 | 0.000 |
| V+J+C | 20584 | 15.709 | 0 | 0.000 |
| V+J+C+CDR3 | 20552 | 15.685 | 0 | 0.000 |
| productive_IGH | 1 | 0.001 | 1 | 0.140 |
| missing_V | 95648 | 72.996 | 0 | 0.000 |
| missing_J | 80946 | 61.776 | 712 | 100.000 |
| missing_C | 23967 | 18.291 | 712 | 100.000 |
| missing_CDR3 | 108066 | 82.473 | 711 | 99.860 |

## Contig length and completeness
Length-bin annotation counts are in `length_vs_annotation.tsv`; base columns count all chains, and corresponding `_IGH` columns/count percentages restrict to IGH V hits. Length summaries for IGH annotation categories are in `completeness_summary.tsv`.

For 150-199 bp contigs (n=20409, IGH n=178), IGH V assignment=178, CDR3=0, V+J+CDR3=0, V+J+C+CDR3=0.
For >=1500 bp contigs (n=42, IGH n=0), IGH V assignment=0, CDR3=0, V+J+CDR3=0, V+J+C+CDR3=0.
Across bins, IGH V assignments were 150-199: 178, 200-249: 199, 250-299: 146, 300-399: 79, 400-499: 85, 500-749: 25, 750-999: 0, 1000-1499: 0, >=1500: 0; IGH CDR3 assignments were 150-199: 0, 200-249: 0, 250-299: 1, 300-399: 0, 400-499: 0, 500-749: 0, 750-999: 0, 1000-1499: 0, >=1500: 0; matched primary clonotypes were 150-199: 0, 200-249: 0, 250-299: 0, 300-399: 0, 400-499: 0, 500-749: 0, 750-999: 0, 1000-1499: 0, >=1500: 0. No IGH J or C assignments were present, so no length bin yielded a complete primary clonotype.
These are descriptive length strata; they do not establish that length caused annotation or biological recovery.

## Primary benchmark
- v4 unique clonotypes: 0
- Reference clonotypes: 116632
- Matches: 0
- Precision: NOT AVAILABLE
- Sensitivity: 0.0
- Pearson r: NOT AVAILABLE (v4 path abundance cannot be recovered from stored path tables).

## Additional diagnostics
Relaxed matching results are in `matching_breakdown.tsv`. D and exact/collapsed isotype agreement are in `gene_agreement.tsv`.
Abundance-stratified reference counts are in `abundance_stratified_benchmark.tsv`; v4 abundance bins, precision, and Pearson correlation are unavailable because per-contig abundance is not stored. Matches are assigned to strata by iRepertoire reference copy.

## Comparison
TRUST4 is full FZ-116; v3-debug/v4 are 100K subsets. The v3-debug biological output present in the project is full FZ-116, not the 100K subset; no v3 100K biological evaluation exists, so those values are NOT AVAILABLE. Raw count/metric differences across input scales are not directly equivalent. See `TRUST4_V3_V4_COMPARISON.tsv`.

## Interpretation
### 1. Assembly-level observation
The v4 diagnostic assembly contains 131,032 contigs; its length statistics are in the adjacent computational diagnostics. This is an assembly-structure observation only.

### 2. Annotation-level observation
The yield and length-stratified tables show how V, J, CDR3, and productive IGH calls vary with contig length. Increased length is not itself evidence of improved biological recovery.

### 3. Benchmark-level observation
The primary published-key evaluation yielded 0 validated clonotypes from 0 v4 candidate clonotypes against 116632 reference clonotypes. Abundance concordance is unavailable.

The TRUST4 baseline was independently rerun in memory with `scripts/benchmark_trust4_published.py` and reproduced the specified values before v4 metrics were generated. No thresholds were tuned.

## Limitations
- v4 was evaluated on a 100K-read subset.
- iRepertoire is an external matched reference and is not necessarily molecule-identical to RNA-seq.
- Low-abundance reference clonotypes can reduce sensitivity.
- No threshold tuning was performed.
- v4 path tables omit bundle/node membership and per-contig abundance; Pearson abundance correlation and v4 abundance-stratified precision cannot be calculated.
- The available v3-debug biological results are full FZ-116, not a comparable 100K evaluation.
