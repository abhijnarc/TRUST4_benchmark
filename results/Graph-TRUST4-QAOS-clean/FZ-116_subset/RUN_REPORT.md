# Graph-TRUST4-QAOS-clean FZ-116 Subset Validation

## Scope and status

This report evaluates the already completed 100,000-pair FZ-116 subset run.
The graph assembly was not rerun and the full FZ-116 dataset was not started.
The original TRUST4 extractor and annotator were reused; the original TRUST4
source was not modified.

The graph command exited with status 0, delivered no signals, and used no swap.
The clean assembler smoke test passes. The output FASTA has 51,777 records,
matching the 51,777 annotation records. Every annotation sequence matches its
FASTA record either directly (26,839) or as its reverse complement (24,938);
there are no sequence mismatches.

## Computational performance

| Measure | Result |
|---|---:|
| Input | 100,000 raw read pairs |
| Candidate extraction | 100,000 pairs (200,000 reads), fraction 1.0 |
| Requested / OpenMP threads | 64 / 64 |
| Graph-stage wall time (`assembly.time`) | 6:52.65 |
| Graph internal elapsed time | 409.521 s |
| Graph command CPU utilization | 5,763% (about 57.6 CPU cores on average) |
| Internal utilization of requested capacity | 90.7392% |
| Graph-stage user / system CPU time | 23,772.10 s / 10.59 s |
| Peak RSS | 1,777,680 KB (about 1,736 MiB) |
| Full candidate-extraction + graph + annotation workflow | 7:35.09 |
| Full-workflow average CPU utilization | 5,277% |

The graph stage used close to, but not all, 64 cores continuously. The graph
stage and full-workflow timings are distinct; the retained `assembly.time`
file measures the graph command.

## Graph statistics

| Statistic | Result |
|---|---:|
| Raw / candidate reads | 200,000 / 200,000 |
| Graph nodes | 200,000 |
| Candidate k-mer hit events | 11,296,335,957 |
| Oriented candidate pairs scored | 9,999,761 |
| Accepted overlap edges | 6,658,585 |
| Pairs rejected by alignment thresholds | 3,341,202 |
| Candidate-hit events not retained at the cap | 11,286,255,756 |
| Connected components | 23,717 |
| Isolated nodes | 22,598 (11.2995%) |
| Branch nodes | 170,677 |
| Cycles detected | 691 |
| Paths emitted | 51,777 |
| Components reaching the 100-path cap | 157 |
| Maximum nodes per path | 64 |
| Edge/node ratio | 33.293 |
| Largest connected component | 144,616 nodes (72.308% of reads) |

The k-mer-hit counts are repeated seed-hit events, not unique read pairs. The
edge count includes alternative qualifying overlap lengths/orientations
between reads. Branch-resolution and ambiguous-branch counts are not defined
by this implementation: it preserves alternatives and bounds path traversal
rather than resolving branches with biological evidence.

Only 31,986 of 200,000 distinct candidate read nodes occur in any emitted
path (15.993%); 168,014 are absent from emitted paths. Of the 51,777 emitted
paths, 30,297 contain one node. The recorded statistics do not identify how
many omitted nodes are attributable individually to path caps, boundary-state
selection, or cycle termination. The 157 capped components show that the
bounded path limit is reached in a nontrivial number of components.

## QAOS and identity

| Measure on accepted edges | Result |
|---|---:|
| Minimum overlap | 31 bp |
| Mean / median overlap | 32.2515 / 31 bp |
| Mean / median identity | 0.999442 / 1.0 |
| Mean / median QAOS | 0.993150 / 0.997590 |
| Accepted edges with identity exactly 1.0 | 6,527,625 (98.0332%) |
| Accepted edges with identity below 1.0 | 130,960 (1.9668%) |
| Accepted-edge identity/QAOS Pearson correlation | 0.83176 |

Accepted edges pass separate identity >= 0.90 and QAOS >= 0.90 filters. Their
QAOS distribution is 145,629 in [0.90, 0.95), 1,116,403 in [0.95, 0.99), and
5,396,553 >= 0.99; none of the accepted edges is below 0.90. QAOS is therefore
used by the edge evaluator, and identity remains separately recorded.
However, all 3,341,202 rejected pairs are combined into one counter, so the
run outputs do not reveal how many would have passed identity but failed
QAOS. The fraction of edges marginally rejected by QAOS cannot be measured
from this run. The high identity of accepted edges alone does not establish
that QAOS had no effect.

## Biological reconstruction and annotation yield

| Measure | Graph-TRUST4-QAOS-clean |
|---|---:|
| Total contigs / annotation records | 51,777 / 51,777 |
| Contig length min / median / 75th percentile / max | 150 / 150 / 246 / 841 bp |
| N50 / assembled bases | 212 bp / 10,425,308 bp |
| Contigs with IGH V assignment | 1,543 |
| V / D / J / C assignments | 14,753 / 0 / 1,065 / 14,713 |
| CDR3 nucleotide assignments | 681 |
| Productive IGH annotations | 8 |
| V+J assignments, all chains | 170 |
| V+J+CDR3 assignments, all chains | 166 |
| V+J+C+CDR3 assignments, all chains | 132 |
| IGH V+J assignments | 0 |
| Complete IGH V+J+CDR3 reconstructions | 0 |
| Unique all-chain clonotypes under the corrected nucleotide key | 22 |
| Unique IGH clonotypes under the corrected nucleotide key | 0 |

Thus the candidate-read-to-emitted-path retention is 15.993%. The
contig-to-complete-V/J/CDR3 retention across all chains is 166/51,777
(0.3206%); requiring C as well gives 132/51,777 (0.2549%). There are no
complete IGH V+J assignments, despite 1,543 contigs receiving an IGH V call.
The primary biological loss is therefore the combination of poor read-node
coverage by emitted paths and failure to obtain compatible IGH V/J
assignments from the annotated paths. The existing aggregate outputs do not
permit assigning every missing read node to a single traversal cause.

## Corrected reference key and same-subset TRUST4 comparison

The comparison uses the corrected nucleotide-based procedure in
`scripts/evaluate_pure_graph_fz116.py`, not the earlier amino-acid CDR3 key
documented as problematic in `results/benchmark/FZ-116/CRITICAL_FINDINGS.md`.
V/J alleles are removed, the biological C isotype is normalized using the
established subclass-collapse rule, and TRUST4's junction nucleotide
sequence is uppercased and trimmed by three bases at each end before matching
the iRepertoire nucleotide CDR3. The frozen reference normalization produces
116,632 unique reference clonotypes. The subset has no complete IGH
clonotypes to match: Graph-TRUST4-QAOS-clean has 0 matches, and the
same-subset TRUST4 annotation has 0 matches.

| Metric | TRUST4, same subset | Graph-TRUST4-QAOS-clean | Difference |
|---|---:|---:|---:|
| Contigs | 707 | 51,777 | +51,070 |
| N50 (bp) | 277 | 212 | -65 |
| Assembled bases | 180,977 | 10,425,308 | +10,244,331 |
| IGH V-assigned contigs | 6 | 1,543 | +1,537 |
| V / J / C / CDR3 assignments | 387 / 287 / 421 / 227 | 14,753 / 1,065 / 14,713 / 681 | +14,366 / +778 / +14,292 / +454 |
| Productive IGH annotations | 0 | 8 | +8 |
| Complete V+J+CDR3, all chains | 87 | 166 | +79 |
| Complete IGH V+J+CDR3 | 0 | 0 | 0 |
| Unique all-chain clonotypes | 68 | 22 | -46 |
| Unique IGH clonotypes | 0 | 0 | 0 |
| Reference clonotypes | 116,632 | 116,632 | 0 |
| Reference matches | 0 | 0 | 0 |
| Precision | NOT AVAILABLE | NOT AVAILABLE | NOT COMPARABLE |
| Sensitivity | 0.0 | 0.0 | 0.0 |
| Pearson abundance correlation | NOT AVAILABLE | NOT AVAILABLE | NOT COMPARABLE |
| Runtime | NOT MEASURED | 455.09 s, full workflow | NOT COMPARABLE |
| Peak RSS | NOT MEASURED | 1,736 MiB | NOT COMPARABLE |

Precision is not defined because neither same-subset annotation produced an
eligible complete IGH clonotype. Sensitivity is 0/116,632 for both under the
full matched-sample reference denominator; because this is a subset run, it
should not be interpreted as a full-sample sensitivity estimate. The
separately reported frozen full-FZ-116 TRUST4 benchmark (11,233 unique
clonotypes, 3,815 matches, precision 0.339624321196, sensitivity
0.032709719459, Pearson r 0.660440783353) is context only and is not the
same-subset comparator.

No defensible per-clonotype abundance is available for the read-level graph
paths. Path node counts are not used as abundance: reads can be represented
in alternative paths, and the graph-node output has no bundle-abundance
measure. Pearson r is therefore NOT AVAILABLE.

## Representative paths

The examples below are records from the completed run, not manually altered
sequences. The first four are complete light-chain (IGK) annotations, not
IGH reconstructions. No complete IGH example exists. `Path score` is the
implementation's raw `cumulative_QAOS` field (an accumulated score, not a
probability); supporting nodes are the path's read-node count.

| Contig | Length | V | D | J | C | CDR3 nucleotide | Mean identity | Mean QAOS | Path score | Supporting nodes | Example type |
|---|---:|---|---|---|---|---|---:|---:|---:|---:|---|
| contig_17656 | 269 | IGKV1-39*01\|IGKV1D-39*01 | - | IGKJ4*01 | IGKC*01 | TGTCAACAGAGTTACAGTACCCCTCTCACTTTC | 1.000000 | 0.940748 | 0.940748 | 2 | Complete IGK |
| contig_18536 | 269 | IGKV1-16*01 | - | IGKJ4*01 | IGKC*01 | TGCCAACAGTATAATAGTTACCCTCTCACTTTC | 1.000000 | 0.977514 | 0.977514 | 2 | Complete IGK |
| contig_19886 | 265 | IGKV1D-12*01 | - | IGKJ4*01 | IGKC*01 | TGTCAACAGGCTAACAGTTTCCCTCCCACTTTC | 1.000000 | 0.997268 | 0.997268 | 2 | Complete IGK |
| contig_19914 | 263 | IGKV1-17*01 | - | IGKJ1*01 | IGKC*01 | TGTCTACAGCATAATAGTTACCCTCGGACGTTC | 1.000000 | 0.995424 | 0.995424 | 2 | Complete IGK |
| contig_37505 | 592 | - | - | - | TRGC2*06 | - | 0.994565 | 0.957734 | 3.830940 | 5 | Long, no V/J/CDR3 |
| contig_37511 | 592 | - | - | - | TRGC2*07 | - | 0.997101 | 0.985572 | 4.927860 | 6 | Long, no V/J/CDR3 |
| contig_37518 | 592 | - | - | - | TRGC2*06 | - | 1.000000 | 0.993602 | 4.968010 | 6 | Long, no V/J/CDR3 |
| contig_0 | 213 | - | - | - | - | - | 1.000000 | 0.998797 | 0.998797 | 2 | One of 12 paths in component 0 |
| contig_1 | 213 | - | - | - | - | - | 0.988506 | 0.964727 | 0.964727 | 2 | One of 12 paths in component 0 |
| contig_3 | 213 | - | - | - | - | - | 1.000000 | 0.998797 | 0.998797 | 2 | One of 12 paths in component 0 |

## Implementation verification and output limitations

| Intended behavior | Observed implementation |
|---|---|
| Candidate extraction inherited from TRUST4 | Yes; unmodified `fastq-extractor` produced the 100,000 candidate pairs. |
| Non-greedy graph construction | Yes, alternatives are retained; seed postings and oriented candidates are capped. |
| QAOS used in overlap evaluation | Yes; the source applies QAOS and identity as separate edge criteria. |
| Identity remains separately available | Yes; it is emitted alongside QAOS. |
| Abundance as soft evidence | No; read-level nodes have no bundle-abundance model and abundance does not guide paths. |
| Paired-end information | Mate IDs are retained as metadata; mate pairing does not constrain path search. |
| SHM-tolerant consensus | No explicit SHM model; consensus is quality-posterior aggregation. |
| V/D/J/C-aware path resolution | No; TRUST4 annotation is downstream and does not resolve graph paths. |
| Original TRUST4 implementation intact | Yes; original extractor/annotator binaries are reused without source edits. |

Two output diagnostics should not be trusted as currently serialized:

* `consensus_support.tsv` has a selected-base field-format bug at
  `clean_graph.cpp:480`. Of 10,425,308 base rows, 2,424,884 have an empty
  selected-base field and 8,000,424 have only seven fields despite the
  eight-field header. The emitted FASTA consensus is separate; this defect
  invalidates the per-position support table and its selected-base column.
* `graph_edges.bed` uses read-relative intervals but does not transform
  intervals for reverse-oriented endpoints. It is unsuitable for
  orientation-faithful coordinate interpretation. The visualization folder
  also contains DOTs from earlier runs alongside the current selections.

## Interpretation and next step

**Computational assembly result:** the subset workflow completed successfully,
used 64 OpenMP threads, emitted 51,777 paths, and produced internally
count-consistent FASTA and annotation files.

**Biological reconstruction result:** there are zero complete IGH V+J+CDR3
annotations, zero eligible unique IGH clonotypes, and zero reference matches
under the corrected nucleotide key. The extra contigs and higher annotation
assignment totals do not demonstrate improved recall. The primary recall
measure is sensitivity, which is zero for this subset comparison.

**Recommendation:** do not scale this implementation to full FZ-116 yet.
Debug the existing subset's 0 IGH V/J co-assignment and low read-node path
coverage, and correct/validate the support-table serialization before any
further run. No thresholds or assembly parameters were changed for this
evaluation; no new assembly was launched.
