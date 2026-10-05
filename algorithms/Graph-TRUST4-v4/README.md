# Graph-TRUST4-v4

Graph-TRUST4-v4 is an isolated successor to Graph-TRUST4-v3-debug. The v3 and v3-debug directories are unchanged.

## Fixed computational parameters

- k-mer size: 9
- minimum variable overlap: 31 bp
- minimum identity: 0.90
- minimum QAOS: 0.90
- high-quality mismatch threshold: Q20
- maximum candidates per bundle: 500
- maximum k-mer postings: 1000
- batch size: 10,000 bundles
- memory limit: 32 GB RSS target
- OpenMP threads: 64
- branch QAOS delta: 0.01
- maximum alternative paths per branch: 2

## Differences from v3-debug

1. Candidate IDs are ranked before full QAOS scoring by local k-mer support, cheap suffix overlap length, and approximate identity. Only the top 500 are fully scored.
2. Accepted edges retain overlap, identity, QAOS, high-quality mismatches, mean/minimum quality, orientation, and paired-support fields.
3. Paired-end support is explicitly zero/unavailable because the inherited compact loader does not retain mate identifiers. No unsupported mate evidence is fabricated.
4. Path extension does not stop at every branch. Outgoing edges are ranked by QAOS, overlap, high-quality mismatches, identity, paired support, and target abundance. Branches with QAOS delta >= 0.01 resolve to the best edge; ambiguous branches retain at most two alternatives.
5. Each path has local repeated-node and cycle protection, with no uncontrolled path expansion.
6. Contig consensus uses Phred-weighted per-base log likelihoods instead of unweighted sequence concatenation.
7. Every path validates expected length against actual consensus length.

The implementation preserves compact k-mer postings, bounded thread-local candidate/edge buffers, batch processing, variable-length overlaps, QAOS, and RSS progress reporting. No biological annotation or iRepertoire evaluation is run automatically.
