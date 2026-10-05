# Graph-TRUST4-QAOS design

## Scope and integration boundary

`Graph-TRUST4-QAOS` is an isolated experimental assembler. It does not
modify `algorithms/TRUST4` or the frozen `Graph-TRUST4-pure` implementation.
The intended production boundary is:

1. TRUST4 performs read parsing, paired-read handling, k-mer filtering, and
   candidate extraction.
2. This executable consumes the resulting paired candidate FASTQs.
3. Existing `algorithms/TRUST4/annotator` annotates
   `assembled_contigs.fa` and writes AIRR output.

The current executable accepts FASTQ pairs directly, so callers must ensure
that the files represent the candidate-read interface. It does not silently
reimplement TRUST4's rough V/J candidate filter.

## Graph model

Identical sequences are compressed into abundance-bearing bundles. K-mer
postings propose bounded candidate pairs. Each accepted directed edge stores
overlap length, ordinary identity, quality-aware overlap score, mismatch
quality counts, candidate support, and paired-end support. Complete edge and
node provenance is emitted in TSV files.

An edge is retained when **either** identity meets `min_identity` **or** QAOS
meets `min_qaos`; it is rejected only when both fail. Thus QAOS is evidence,
not the sole hard filter, and plausible alternatives remain available for path
enumeration.

## QAOS

For aligned overlap positions, the per-base agreement probability is

`P(match) = 1 - 10^(-Q_left/10) - 10^(-Q_right/10) +
10^(-Q_left/10-Q_right/10)` for matching bases.

For mismatches, the corresponding probability is the product of the two
error probabilities. The edge QAOS is the geometric mean of these
probabilities, with a small numerical floor. High-quality mismatches are
reported separately and are not hidden by the aggregate score.

## Paired-end evidence

The loader maps each mate sequence to its abundance bundle and records
observed bundle pairs. An edge receives paired support when its endpoints were
observed as a mate pair; path statistics sum this support. This is a
conservative compatibility signal: it does not infer insert size or rescue
missing mates, and reverse-complement/orientation constraints remain encoded
by the overlap edge.

## Paths and consensus

The assembler enumerates bounded simple paths per connected component without
greedily selecting one outgoing edge. Each path becomes one consensus contig.
Consensus overlap conflicts are resolved using the higher base quality. Path
and contig provenance retain branch status, edge scores, bundle abundance, and
paired support.

V/D/J/C and CDR3 constraints are intentionally deferred to the existing
TRUST4 annotator. This prevents premature biological filtering and keeps
annotation conventions comparable with the published TRUST4 benchmark.

## Resource controls

`-k`, `-m`, `-t`, `-M`, candidate/posting caps, batch size, and path bounds are
explicit CLI parameters. The run writes requested/OpenMP thread counts,
runtime, peak RSS, and graph/path statistics. No full FZ-116 run is launched
by the build or test targets.
