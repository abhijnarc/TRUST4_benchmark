# PURE Graph FZ-116 Biological Evaluation

## Experiment status

- Assembly reused from the frozen output; Graph-TRUST4-pure was not rerun.
- Annotation reused from `pure_raw_airr.tsv`; no annotation rerun was needed.
- All 209,716 FASTA IDs have exactly one AIRR record. Every AIRR sequence matches its FASTA sequence in the same orientation or as its reverse complement (108,274 same-orientation; 101,442 reverse-complement).
- Matching uses the established published TRUST4 key: trimmed CDR3 nucleotide + normalized V + normalized J + normalized/collapsed C.
- CDR3 nucleotide normalization removes the first and last 3 nt from the TRUST4/AIRR junction. V/J allele suffixes are removed; constant-region subclasses are collapsed as in the published benchmark.
- This does not use the earlier CDR3 amino-acid comparison that produced zero matches.

## COMPUTATIONAL ASSEMBLY RESULT

| Metric | Graph-TRUST4-pure |
|---|---:|
| Contigs | 209,716 |
| N50 | 2,924 bp |
| Assembled bases | 302,589,678 |
| Graph nodes / bundles | 921,684 |
| Accepted edges | 6,027,883 |
| Paths | 209,716 |
| Branches | 71,129 |
| Components | 403,122 |
| Cycles | 0 |
| Raw read records | 2,160,900 |
| Mean overlap | 64.5274 bp |
| Mean identity | 0.993051 |
| Mean QAOS | 0.966499 |

Frozen parameters include k=9, minimum overlap=31, minimum identity=0.9, minimum QAOS=0.9, Q_HIGH=20, candidate cap=500, and 64 threads. Assembly wall time was 43033.6 seconds; peak RSS was 2249 MB.

## Annotation yield and completeness

`pure_raw_airr.tsv` has 209,716 unique annotation records for 209,716 contigs. Every sequence matches its FASTA record either directly or by reverse complement.

| Metric | Count | Percent of contigs |
|---|---:|---:|
| total_contigs | 209,716 | 100.000% |
| IGH_contigs | 26,363 | 12.571% |
| V_assigned | 95,057 | 45.327% |
| D_assigned | 4,234 | 2.019% |
| J_assigned | 13,003 | 6.200% |
| C_assigned | 27,611 | 13.166% |
| CDR3_assigned | 13,446 | 6.412% |
| productive_IGH | 5,223 | 2.491% |

## BIOLOGICAL RECONSTRUCTION RESULT

| Metric | Graph-TRUST4-pure |
|---|---:|
| Pure unique clonotypes | 1,144 |
| iRepertoire reference clonotypes | 116,632 |
| Matched clonotypes | 516 |
| Precision | 0.451049 |
| Sensitivity (primary recall measure) | 0.004424 |
| Pearson abundance correlation | 0.078030 |

The iRepertoire reference has 116,632 unique clonotypes under the frozen CDR3nt + V + J + C benchmark key. No extra CDR3aa-stop filter was applied; this preserves the denominator used by the published TRUST4 comparison.

## Direct comparison with frozen TRUST4 baseline

| Metric | TRUST4 | Graph-TRUST4-pure | Difference (pure - TRUST4) |
|---|---:|---:|---:|
| contigs | 43,720 | 209,716 | 165,996 |
| N50 | 242 | 2,924 | 2,682 |
| IGH_contigs | 11,548 | 26,363 | 14,815 |
| productive_IGH | 11,383 | 5,223 | -6,160 |
| unique_clonotypes | 11,233 | 1,144 | -10,089 |
| matches | 3,815 | 516 | -3,299 |
| precision | 0.339624 | 0.451049 | 0.111425 |
| sensitivity | 0.032710 | 0.004424 | -0.028286 |
| Pearson r | 0.660441 | 0.078030 | -0.582411 |

The TRUST4 assembly figures come from its existing `TRUST_FZ-116_final.out` and AIRR outputs. Its published benchmark metrics were reused without rerunning TRUST4: 11,233 unique clonotypes, 3,815 matches, precision 0.339624, sensitivity 0.032710, and Pearson r 0.660441.

## Abundance and graph provenance

Per-contig graph abundance is the sum of `bundle_abundance` over nodes in that contig's path. Every path-step abundance matches `graph_nodes.tsv`; the sum of graph-node abundances (2,160,900) matches the graph's 2,160,900 raw read records. `node_reads.tsv` has 1,658,216 provenance rows, with at most 64 stored raw read IDs per node. Therefore bundle counts are complete, while the listed raw read IDs are capped samples.

Matched clonotypes link to contigs and path-node IDs in `pure_vs_iRep_matches.tsv`; follow those node IDs in `contig_paths.tsv` and `node_reads.tsv`, with retained alternatives documented in `branch_decisions.tsv`. Alternative paths may share graph nodes and thus read support; graph abundance should be interpreted as path support, not as independent molecule counts.

## Limitations and interpretation

- Only 26,363 of 209,716 assembled contigs have an IGH V assignment; the biological benchmark evaluates the clonotypes that have the required key fields.
- The graph-derived abundance is a sum of constituent bundle read counts. Alternative paths can share nodes, so support across different contigs is not independent.
- Paired-read links are not represented as paired support in this experimental graph output.
- Contig count, N50, and assembled bases are computational assembly measures; they do not establish biological recovery.
- The primary recall comparison is sensitivity: Graph-TRUST4-pure recovers 516 reference clonotypes (0.004424), versus TRUST4's frozen 3,815 matches (0.032710 sensitivity). The graph result therefore does not show improved reference-clonotype recall, despite its larger contig count and N50.

The annotation-yield counts and full benchmark metrics are in `pure_annotation_yield.tsv` and `pure_benchmark_metrics.tsv`; side-by-side values are in `pure_vs_TRUST4.tsv`.
