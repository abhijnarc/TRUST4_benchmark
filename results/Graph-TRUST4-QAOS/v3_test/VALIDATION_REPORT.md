# Graph-TRUST4-QAOS v3 validation report

This report covers only the required FZ-116 subset. No full FZ-116 run was
started and the original TRUST4 source was not modified.

## A. Test status

The TRUST4 candidate extractor, v3 graph run, and downstream TRUST4
annotation all completed with exit status 0 and no reported signals. The v3
FASTA contains 25,095 records and the annotation file contains 25,095
records plus its header. Required v3 outputs are present:

- `candidate_statistics.tsv`
- `graph_statistics.tsv`
- `edge_statistics.tsv`
- `path_statistics.tsv`
- `path_score.tsv`
- `consensus_support.tsv`
- `assembled_contigs.fa`
- `annotated_contigs.tsv`
- `annotation_yield.tsv`
- `comparison.tsv`

## B. Computational performance

| Measurement | v2 QAOS | v3 |
|---|---:|---:|
| Requested/OpenMP threads | 64 / 64 | 64 / 64 |
| Graph elapsed | 1,736.74 s | 1,723.48 s |
| Graph user CPU | 103,945.87 s | 103,132.93 s |
| Graph CPU utilization | 5,664% | 5,734% |
| Approximate average cores | 56.6 | 57.3 |
| Graph peak RSS | 887 MB | 885 MB |
| Annotation elapsed | 52:01.38 | 3:27.45 |
| Annotation peak RSS | ~37.3 GiB | ~5.24 GiB |
| Exit status | 0 | 0 |

The graph used approximately 57 of 64 cores on average. The reduced v3 path
set lowered annotation memory below 32 GiB.

## C. Candidate extraction

The v3 runner invokes TRUST4's existing `fastq-extractor` from
`algorithms/TRUST4/FastqExtractor.cpp`. On this subset:

| Statistic | Value |
|---|---:|
| Raw mate reads | 200,000 |
| Candidate mate reads | 200,000 |
| Candidate pairs | 100,000 |
| Candidate fraction | 1.0 |

The extractor is therefore genuinely used, but this input was entirely
accepted by it. The result is measured rather than substituted with raw-read
assumptions.

## D. Graph, QAOS, and pruning

| Statistic | v2 QAOS | v3 |
|---|---:|---:|
| Graph nodes/bundles | 59,148 | 59,148 |
| Accepted edges | 6,830,654 | 6,830,654 |
| Components | 2,370 | 2,370 |
| Branch events | 179,408 | 179,460 |
| Paths before pruning | 498,351 | 498,368 |
| Paths after pruning | 498,351 | 25,095 |
| Dominance-pruned paths | 0 | 473,273 |
| Biologically inconsistent paths removed | 0 | 0 |
| Ambiguous paths retained | not emitted | 365 |
| N50 | 5,088 | 8,015 |

The v3 graph remains non-greedy: path alternatives are generated first, then
local state-aware dominance pruning retains the best path for a compatible
terminal-node/BCR-evidence state. This is not a global top-N truncation.

The v3 edge model is unchanged from v2. Mean identity is 0.974813, mean QAOS
is 0.867360, and 4,859 accepted edges have identity below 0.90. QAOS remains
active through the existing weighted edge score.

## E. BCR-aware path scoring

Each bundle receives compact V/D/J/C evidence from shared k-mers against the
existing TRUST4/IMGT reference. Path state propagates evidence bits and
records separate terms in `path_score.tsv`:

```text
sequence_score
qaos_score
abundance_score
paired_score
V_score
D_score
J_score
C_score
biological_order_score
total_score
classification
```

Abundance contributes through log abundance support, paired links contribute
through path paired support, and biological ordering contributes to the
reported path score. `consensus_support.tsv` records abundance-weighted base
support for A/C/G/T on representative paths.

Important limitation: this is reference-k-mer evidence propagation, not full
TRUST4 per-read V/D/J/C annotation. It is a conservative lightweight
pre-annotation signal and should not be described as equivalent to the
annotator's final calls.

## F. Biological reconstruction

| Measurement | v2 QAOS | v3 |
|---|---:|---:|
| Assembled contigs | 498,351 | 25,095 |
| V assignments | 236,595 | 13,291 |
| J assignments | 41,520 | 7,809 |
| V+J assignments | 14,420 | 5,619 |
| CDR3 assignments | 28,454 | 5,990 |
| Complete V+J+CDR3 | 14,086 | 5,616 |
| Productive IGH | 1,166 | 7 |
| Unique clonotypes | 90 | 8 |

The large reduction in paths also reduces absolute annotated records. The
retained-path complete-BCR fraction is 22.38% of v3 contigs, compared with
2.83% for v2. This indicates that pruning is enriching for biologically
compatible paths, but the productive IGH count and clonotype recovery still
require investigation.

## G. Comparison with v2

The exact machine-readable comparison is in [comparison.tsv](./comparison.tsv).
Measurements are objective; no superiority claim is made.

| Metric | v2 QAOS | v3 | Difference |
|---|---:|---:|---:|
| Graph nodes | 59,148 | 59,148 | 0 |
| Accepted edges | 6,830,654 | 6,830,654 | 0 |
| Components | 2,370 | 2,370 | 0 |
| Paths before pruning | 498,351 | 498,368 | +17 |
| Paths after pruning | 498,351 | 25,095 | -473,256 |
| Complete V+J+CDR3 | 14,086 | 5,616 | -8,470 |
| Productive IGH | 1,166 | 7 | -1,159 |
| Runtime seconds | 1,736.74 | 1,723.48 | -13.26 |
| Peak graph RSS MB | 887 | 885 | -2 |

## H. Implementation sanity checks

| Requirement | Result |
|---|---|
| Non-greedy graph | Pass: alternatives generated before pruning |
| QAOS active | Pass: inherited v2 weighted score and low-identity accepted edges |
| V/J/C enters path scoring | Pass, via reference-k-mer evidence bits |
| Abundance enters path scoring | Pass, as a soft log abundance term |
| Paired evidence enters path scoring | Pass, exact bundle-pair support term |
| Dominance pruning occurs | Pass: 473,273 paths removed |
| Consensus is no longer higher-quality-only | Pass: A/C/G/T abundance-weighted support is emitted |
| Original TRUST4 intact | Pass |
| Annotation memory <=32 GiB | Pass for v3 subset: ~5.24 GiB |

The paired-end implementation remains conservative: it uses exact bundle
pair links, not insert-size inference or full graph-position compatibility.
The biological evidence is also lightweight reference-k-mer evidence rather
than full TRUST4 annotation propagated into every node.

## I. Recommendation

Do not start full FZ-116 yet. The v3 architecture demonstrates the intended
path-reduction mechanism and reduces annotation memory substantially, but the
drop from 14,086 complete records in v2 to 5,616 in v3, together with only
7 productive IGH records, requires debugging of the BCR evidence and
dominance state before scaling. The next step should validate retained and
removed paths against full TRUST4 annotation on a smaller diagnostic subset,
then refine state compatibility without changing the QAOS overlap formula.
