# Graph-TRUST4-v4

This isolated implementation preserves v3-debug memory limits, k=9 compact indexing, variable QAOS overlaps, and 64-thread batched candidate generation. v4 adds cheap candidate ranking, branch-aware bounded alternatives, lexicographic edge ranking, and Phred-weighted consensus.

Paired-end support is recorded as zero/unavailable because the inherited loader does not retain mate identifiers. No biological evaluation is run automatically.
