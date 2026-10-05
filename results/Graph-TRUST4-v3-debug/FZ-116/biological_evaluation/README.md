# FZ-116 Graph-TRUST4 biological evaluation

Annotation command: `algorithms/TRUST4/annotator -f reference/TRUST4/human_IMGT+C.fa -a results/Graph-TRUST4-v3-debug/FZ-116/graph_fz116/assembled_contigs.fa --fasta --needReverseComplement --noImpute -t 64 --outputFormat 1`.

Evaluation command: `python3 scripts/evaluate_graph_fz116.py`.

Reference: `reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz`. TRUST4 baseline: `results/TRUST4/FZ-116/TRUST_FZ-116_report.tsv`.

Normalization follows `scripts/benchmark_trust4_published.py`: first gene candidate, remove allele suffixes, collapse IGHA1/2 to IGHA and IGHG1-4 to IGHG (plus the existing IGHD/IGHE/IGHM mappings), uppercase sequences, trim 3 nt from each CDR3 end, and coalesce duplicate keys by maximum abundance. Primary key is normalized V + J + collapsed C + trimmed CDR3nt; D is auxiliary.

Graph abundance is the sum of `abundance` values in `bundles.tsv` for the unique bundle nodes on the reconstructed contig path. Paths are reconstructed from `graph_edges.tsv` using the validated path traversal; contigs without a path mapping are excluded from biological clonotype evaluation and counted.

Records are excluded from clonotype evaluation when they are not IGH or lack a valid CDR3 after the published 3-nt trim; missing V/J/C values and unmapped paths are retained and counted in the metrics table. Duplicate normalized keys are coalesced by maximum path abundance. Abundance min/median/mean/max, assembly cycle/repetition/length checks, and all metric values are recorded in `graph_fz116_benchmark_metrics.tsv` and `graph_fz116_normalized.tsv`. No assembly, parameter tuning, other samples, or iRepertoire-guided filtering was performed.
