# Graph-TRUST4-QAOS subset validation report

This report audits the already completed subset run. No assembler source,
threshold, output, or full FZ-116 run was changed or started for this audit.

## A. Test status

Inputs were `test_data/FZ116_test_1.fq` and `test_data/FZ116_test_2.fq`.
The graph run used `-k 9 -m 31 -t 64 -M 32` and exited with status 0.
The existing TRUST4 run on the same pair of input files also exited with
status 0. TRUST4 annotation of the graph FASTA exited with status 0.

The FASTA has 465,135 non-empty records, and `annotated_contigs.tsv` has
465,135 records plus its header. `path_statistics.tsv` also has 465,135
paths. The output files are non-empty and their counts agree. No crash,
signal, or truncation was found in the `/usr/bin/time -v` logs.

The complete logs are [subset.time](../../../logs/Graph-TRUST4-QAOS/subset.time),
[annotation.time](../../../logs/Graph-TRUST4-QAOS/annotation.time), and
[trust4_subset.stdout](../../../logs/Graph-TRUST4-QAOS/trust4_subset.stdout).

## B. Computational performance

| Measurement | Graph-TRUST4-QAOS | TRUST4 annotation |
|---|---:|---:|
| Requested threads | 64 | 64 |
| OpenMP maximum threads | 64 | not emitted by annotator |
| Average CPU utilization | 5,747% | 2,384% |
| Approximate average cores utilized | 57.5 | 23.8 |
| Elapsed time | 29:38.61 | 52:10.15 |
| User CPU time | 102,192.55 s | 74,373.62 s |
| System CPU time | 32.08 s | 256.22 s |
| Peak RSS | 1,009,824 KB (~0.96 GiB) | 37,859,644 KB (~36.1 GiB) |
| Signals / exit status | 0 / 0 | 0 / 0 |

The graph process reported `requested_threads=64 omp_max_threads=64`.
Its average utilization was high but not equivalent to all 64 cores
continuously being busy. A peak CPU percentage was not recorded, so it cannot
be claimed from these artifacts. The graph stage stayed below the 32 GiB
assembler limit. The downstream TRUST4 annotation process peaked at about
36.1 GiB, exceeding a 32 GiB whole-pipeline memory budget.

The machine-level logs are [subset.time](../../../logs/Graph-TRUST4-QAOS/subset.time)
and [annotation.time](../../../logs/Graph-TRUST4-QAOS/annotation.time).

## C. Graph statistics

| Statistic | Value |
|---|---:|
| Input reads counted | 200,000 mate reads (100,000 pairs) |
| Candidate extraction | not performed by QAOS; direct FASTQ loading |
| Graph nodes/bundles | 59,148 |
| Candidate pairs before cap | 233,810,604 |
| Candidate pairs rejected by cap | 209,289,599 |
| Candidate pairs evaluated | 24,521,005 |
| Accepted directed edges | 7,871,713 |
| Connected components | 1,998 |
| Isolated nodes | 620 (1.048%) |
| Largest component | 39,118 nodes (66.136% of nodes) |
| Branch records by source out-degree > 1 | 49,287 |
| Branch events encountered during path enumeration | 158,496 |
| Resolved branches | 0 |
| Ambiguous/retained branches | all reported branches; no biological resolver |
| Cycles | 0 reported; paths use visited-node simple-path traversal |
| Enumerated graph paths/contigs | 465,135 |
| Edge/node ratio | 133.08 |

The high path count is not equivalent to high biological reconstruction:
bounded path enumeration emits many alternatives from a highly connected
graph. The graph itself reports retained alternatives rather than resolving
them with BCR evidence.

## D. QAOS statistics

| Statistic | Value |
|---|---:|
| Mean overlap | 81.4257 bp |
| Median overlap | 76 bp |
| Mean identity | 0.967981 |
| Median identity | 0.972973 |
| Mean QAOS | 0.829817 |
| Median QAOS | 0.866340 |
| Identity/QAOS Pearson correlation | 0.889972 |

Accepted-edge QAOS distribution:

| QAOS range | Edges | Fraction |
|---|---:|---:|
| `< 0.5` | 344,664 | 4.38% |
| `0.5–<0.7` | 1,086,743 | 13.81% |
| `0.7–<0.8` | 1,285,738 | 16.33% |
| `0.8–<0.9` | 2,052,181 | 26.07% |
| `>= 0.9` | 3,102,387 | 39.41% |

All 7,871,713 accepted edges also met the ordinary identity threshold of
0.90. Therefore QAOS was calculated and is distributed broadly, but it did
not change acceptance decisions in this run: every accepted edge would have
passed identity alone. The OR acceptance policy and the absence of any
accepted edge with identity below 0.90 make QAOS non-influential for this
particular edge set.

## E. Biological reconstruction statistics

| Measurement | TRUST4 exact subset | Graph-TRUST4-QAOS |
|---|---:|---:|
| Reconstructed contigs | 707 | 465,135 |
| V assignments | 149 | 211,295 |
| J assignments | 148 | 56,484 |
| V+J assignments | 148 | 5,297 |
| CDR3/junction assignments | 149 | 32,274 |
| Complete V+J+CDR3 records | 148 | 5,050 |
| Complete V+D+J field (`complete_vdj=T`) | not reported here | 0 |

The graph FASTA-to-annotation retention is 100% by record count, but biological
retention is:

- V assignment: 45.42% of graph contigs.
- J assignment: 12.14%.
- V+J: 1.14%.
- V+J+CDR3: 1.09%.

There is no valid candidate-read-to-contig retention ratio for QAOS because
the executable did not inherit TRUST4 candidate extraction; it loaded the
200,000 input mate reads directly. TRUST4's log reports 196,852 reads found
and 75,029 reads assembled after its own candidate filtering. This difference
must not be interpreted as an assembler retention comparison.

The largest graph-to-biology loss is at the V/J compatibility stage:
211,295 V calls fall to 5,297 V+J calls. A further 247 V+J records lack a
CDR3. This confirms the earlier failure mode: many graph paths do not become
complete BCR reconstructions.

## F. Same-subset TRUST4 comparison

| Metric | TRUST4 | Graph-TRUST4-QAOS | Difference (Graph - TRUST4) |
|---|---:|---:|---:|
| Total reconstructed contigs | 707 | 465,135 | 464,428 |
| V assignments | 149 | 211,295 | 211,146 |
| J assignments | 148 | 56,484 | 56,336 |
| CDR3s | 149 | 32,274 | 32,125 |
| Complete V+J+CDR3 BCRs | 148 | 5,050 | 4,902 |
| Unique clonotypes | 127 | 19 | -108 |
| Reference clonotype matches | 0 | 0 | 0 |
| Precision | 0.000000 | 0.000000 | 0.000000 |
| Sensitivity | 0.000000 | 0.000000 | 0.000000 |
| Abundance correlation | NOT AVAILABLE | NOT AVAILABLE | NOT AVAILABLE |
| Runtime | not captured for TRUST4 wrapper | 1,778.61 s | NOT AVAILABLE |
| Peak memory | not captured for TRUST4 wrapper | ~0.96 GiB | NOT AVAILABLE |

The exact-subset clonotype comparison uses the established normalized
CDR3-nucleotide + V + J + C definition. No superiority conclusion is drawn.
Pearson correlation is unavailable because no defensible comparable
clonotype-abundance value was recovered for this subset comparison.

## G. Main failure point

The implementation does form a large graph and emits many paths, but it does
not perform BCR-aware path resolution. The source enumerates paths using
graph topology and a visited-node bound; it does not use V, D, J, C, CDR3, or
locus compatibility while resolving branches. Consequently:

```text
59,148 graph nodes
7,871,713 edges
465,135 paths
211,295 V calls
56,484 J calls
5,297 V+J calls
5,050 complete V+J+CDR3 records
```

The exact stage responsible for the major loss is annotation compatibility
after path generation, compounded by unresolved alternative paths. This is
not evidence that the annotation executable is incomplete: it produced one
record for every graph contig.

## H. Representative examples

The full 15-record table is in
[representative_examples.tsv](./representative_examples.tsv). It contains
five complete V+J+CDR3 examples, five long unannotated contigs, and five
paths marked `AMBIGUOUS_RETAINED`, with contig ID, length, V/D/J/CDR3,
mean identity, geometric-mean QAOS path score, node count, and summed
bundle abundance. The latter is path support, not an independent molecule
count; paths share graph nodes.

Representative categories include complete records such as contigs
70298–70303, long unannotated contigs such as 294093 and 91798, and retained
ambiguous paths beginning at contigs 0–4. No sequences were modified.

## I. Implementation verification

| Intended property | Audit result |
|---|---|
| Candidate extraction inherited from TRUST4 | **Not satisfied by this executable.** It directly loads input FASTQs and does not call TRUST4 preprocessing or `fastq-extractor`. |
| Non-greedy graph construction/path generation | Satisfied structurally: bounded simple-path enumeration retains alternatives and does not choose one outgoing edge. |
| QAOS overlap evaluation | Satisfied: `score_edge` computes quality-aware probabilities and stores QAOS separately. |
| Identity retained separately | Satisfied and reported per edge/path. |
| Abundance as soft evidence | **Not satisfied as a scoring input.** Bundle abundance is stored and emitted, but does not influence edge acceptance, path selection, or consensus. |
| Paired-end information | Partially satisfied: exact sequence-pair links contribute a binary endpoint support value; insert-size/orientation compatibility is not modeled. |
| SHM-tolerant consensus | **Not satisfied as promised.** Consensus chooses the higher-quality base on overlap conflicts; it has no SHM model or abundance-weighted consensus. |
| V/D/J/C in path resolution | **Not satisfied.** Annotation is downstream only; no V/D/J/C compatibility is used for branch/path resolution. |
| Original TRUST4 intact | The QAOS implementation is in a separate directory and the original TRUST4 executable/source was used unchanged for the subset comparison. |

The implementation therefore functions as a quality-scored, paired-link-aware
overlap graph prototype, but it does not yet implement the complete
TRUST4-inherited candidate pipeline or the promised BCR-aware/SHM-aware path
resolution.

## J. Recommendation

Do not scale this implementation to full FZ-116 yet. The subset completed
without a graph crash and within the assembler memory budget, but validation
identified substantive functional gaps: direct raw-read loading rather than
TRUST4 candidate extraction, QAOS not affecting any accepted edge decisions,
no BCR-aware path resolution, no SHM-aware consensus, and downstream
annotation exceeding 32 GiB. The next experimental step should be a debugging
and instrumentation phase addressing these discrepancies on a small fixed
subset before any full-dataset run.
