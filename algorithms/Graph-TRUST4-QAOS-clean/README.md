# Graph-TRUST4-QAOS-clean

This is an independent, minimal experimental read-overlap assembler. It does
not modify or depend on the Graph-TRUST4-QAOS v2/v3 module. TRUST4 is used
only for candidate extraction and downstream annotation.

## Workflow

```text
raw paired FASTQ
  -> algorithms/TRUST4/fastq-extractor
  -> one graph node per candidate read (no read bundling)
  -> bounded k-mer candidate generation
  -> explicit suffix/prefix overlap, identity, and per-base QAOS
  -> directed bidirected overlap graph retaining all accepted alternatives
  -> bounded deterministic path traversal
  -> base-quality probability aggregation across all path reads
  -> existing TRUST4 annotator
```

Candidate extraction is performed by the unmodified
`algorithms/TRUST4/fastq-extractor` (`FastqExtractor.cpp`, including
`IsGoodCandidate` and `ProcessReads_Thread`). R1/R2 record identity and mate
IDs are preserved in node and path provenance. Pairing is metadata only; it
does not constrain path selection.

## QAOS formula

For an observed base `x` with Phred quality Q, the probability of true base b
is:

```text
P(b | x,Q) = 1-e              if b=x
             e/3              otherwise
e = 10^(-Q/10)
```

For each aligned overlap position i, agreement probability is
`p_i = sum_b P(b|x_i,Qx) P(b|y_i,Qy)`. The edge QAOS is the geometric mean
`exp(sum(log(p_i))/overlap_length)`. Thus high-quality mismatches are stronger
contradictions than low-quality mismatches. Identity is the unweighted
fraction of matching bases and remains a separate metric. An edge is accepted
only when overlap, identity, and QAOS each meet their configured minimum.

## Graph search

The k-mer index stores oriented read-prefix seeds and is only used to generate
candidate overlaps. Candidate lists are capped per read and postings are
capped per seed. Each candidate then receives explicit overlap evaluation.
Per-read oriented candidate storage is capped while seeds are consumed; the
output distinguishes raw k-mer candidate-hit events from unique oriented
pairs actually scored and reports hit events not retained after the cap.
All qualifying overlap-length/orientation edges are retained, including
parallel edges between the same reads when overlap length or strand evidence
differs; there is no best-next-read selection, global ranking, BCR evidence,
path dominance pruning, or abundance filtering.

Traversal enumerates directed simple paths from oriented graph boundary
states. Cycles are detected using the visited oriented-node set. Maximum
nodes/path and paths/component are configurable and capped components are
counted in `graph_statistics.tsv`; path termination reason and read-level
provenance are emitted for each path.

## Consensus

Every read in a path is placed using the overlap offsets. At each consensus
position, each observed base contributes its four-base posterior probability
vector derived from its Phred quality. Vectors are summed across all reads;
the maximum-support nucleotide is selected and its fraction of total support
is reported as confidence in `consensus_support.tsv`. This is a conservative
quality-aware consensus; no somatic-hypermutation model or abundance weighting
is applied.

## Parameters

Defaults: `-k 9 -m 31 --max-overlap 150 -i 0.90 -q 0.90
--max-postings 1000 --max-candidates 50 --max-path-nodes 64
--max-paths 100 -t 64`.

## Output and visualization

The output includes candidate statistics, read nodes, overlap edges, paths,
FASTA consensus, per-position consensus support, BED/BEDPE-style files, and
small Graphviz component snapshots. BED files use read- or path-relative
coordinates, not invented genomic coordinates. `graph_edges.bed` is
BEDPE-style because overlap links connect two independent read sequences and
cannot be represented faithfully as a single-interval BED feature.

Selected `.dot` files can be rendered with Graphviz, for example:

```bash
dot -Tsvg results/Graph-TRUST4-QAOS-clean/FZ-116/visualization/component_0.dot \
  -o component_0.svg
```

The diagnostic visualizations include small representative components; the
largest component is labelled only as a long-path candidate, not asserted to
be a BCR.

## Annotation

`scripts/run_graph_qaos_clean.sh` runs the existing TRUST4 annotator with
`--fasta --needReverseComplement --noImpute --outputFormat 1`. Annotation is
batched into bounded FASTA chunks with identical options to control memory.
The same batching helper is used on the original TRUST4 subset assembly for
the side-by-side annotation comparison.

The first required test is the FZ-116 subset. Full FZ-116 execution is not
part of any build or test target.
