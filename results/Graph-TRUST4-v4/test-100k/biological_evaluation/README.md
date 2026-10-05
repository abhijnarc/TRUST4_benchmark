# Graph-TRUST4-v4 100K Biological Evaluation

## Inputs
- FASTA: `results/Graph-TRUST4-v4/test-100k/assembled_contigs.fa`
- iRepertoire: `reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz`
- TRUST4: `results/TRUST4/FZ-116/TRUST_FZ-116_report.tsv`

## Reproduction
Run from project root. Annotation writes only into this evaluation folder:

```sh
algorithms/TRUST4/annotator -f reference/TRUST4/human_IMGT+C.fa -a results/Graph-TRUST4-v4/test-100k/assembled_contigs.fa --fasta --needReverseComplement --noImpute -t 64 --outputFormat 1 > results/Graph-TRUST4-v4/test-100k/biological_evaluation/graph_v4_raw_airr.tsv
python3 results/Graph-TRUST4-v4/test-100k/biological_evaluation/evaluate_graph_v4_biology.py
```

The evaluation script imports the functions in `scripts/benchmark_trust4_published.py` and first verifies the TRUST4 baseline. It stops before v4 output generation if any supplied baseline count/metric differs by more than 1e-9.

## Normalization and matching
Use uppercase CDR3 nucleotide sequence; trim 3 nt from both ends of TRUST4-style junctions; remove V/J/D/C allele suffixes; strip iRepertoire's leading `h`; collapse IGHA1/2 to IGHA, IGHG1-4 to IGHG, IGHD1-4 to IGHD, IGHE1/2 to IGHE, and IGHM1/2 to IGHM. Primary key is trimmed CDR3nt + normalized V + normalized J + collapsed C/isotype; D is diagnostic only. Incomplete primary keys are excluded. Identical v4 primary keys count once. Since no v4 abundance is stored, duplicates cannot be ranked by maximum abundance; the first contig is retained only for match-table traceability.

TRUST4 and reference duplicate keys are coalesced by maximum count/copy using the published functions. Precision is matches / candidate unique clonotypes; sensitivity is matches / reference unique clonotypes. Baseline Pearson is the published Pearson correlation between TRUST4 count and iRepertoire copy among matches. V4 Pearson is unavailable because graph path node membership/abundance is not present; contig length is never used as a substitute.

## Outputs
`graph_v4_raw_airr.tsv` preserves raw annotator output; `graph_v4_annotated.tsv` includes standardized fields plus raw fields; `graph_v4_normalized.tsv` contains normalized calls. Other outputs: `annotation_yield.tsv`, `length_vs_annotation.tsv`, `graph_v4_benchmark_metrics.tsv`, `matching_breakdown.tsv`, `gene_agreement.tsv`, `abundance_correlation.tsv`, `graph_v4_vs_iRep_matches.tsv`, `abundance_stratified_benchmark.tsv`, `TRUST4_V3_V4_COMPARISON.tsv`, `completeness_summary.tsv`, and `V4_100K_BIOLOGICAL_EVALUATION.md`.

No assembly, threshold tuning, or biological filtering based on iRepertoire was performed. The existing v3 biological evaluation is full FZ-116; no v3 100K biological output was found.
