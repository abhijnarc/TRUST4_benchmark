# Graph-TRUST4-QAOS v2 subset validation

This report covers only the required `FZ116_test_1.fq` / `FZ116_test_2.fq`
subset. No full FZ-116 run was started. Original TRUST4 source files remain
unchanged.

## A. Test status

The TRUST4 candidate extractor, identity-only graph ablation, and combined
identity+QAOS graph run all exited with status 0 and no signals. Both graph
FASTAs and both annotation files are non-empty and record-complete:

| Run | FASTA records | Annotation records | Path records |
|---|---:|---:|---:|
| Identity-only | 464,610 | 464,610 | 464,610 |
| Identity + QAOS | 498,351 | 498,351 | 498,351 |

The output logs are in
`logs/Graph-TRUST4-QAOS/v2_qaos.time`,
`logs/Graph-TRUST4-QAOS/v2_identity.time`,
`logs/Graph-TRUST4-QAOS/v2_qaos_annotation.time`, and
`logs/Graph-TRUST4-QAOS/v2_identity_annotation.time`.

## B. Computational performance

| Measurement | Identity-only | Identity + QAOS |
|---|---:|---:|
| Requested/OpenMP threads | 64 / 64 | 64 / 64 |
| Graph elapsed | 1,734.02 s | 1,736.74 s |
| Graph user CPU | 103,904.54 s | 103,945.87 s |
| Graph CPU utilization | 5,678% | 5,664% |
| Approximate average cores | 56.8 | 56.6 |
| Graph peak RSS | 982 MB | 887 MB |
| Annotation elapsed | 54:19.03 | 52:01.38 |
| Annotation user CPU | 73,708.74 s | 83,588.30 s |
| Annotation CPU utilization | 2,278% | 2,695% |
| Annotation peak RSS | ~36.2 GiB | ~37.3 GiB |

The graph stage used approximately 56–57 cores on average, not all 64
continuously. The graph remained below the 32 GiB assembler budget. The
downstream TRUST4 annotator still exceeded the desired 32 GiB complete
pipeline target.

## C. Candidate extraction

The v2 runner calls the existing TRUST4 `fastq-extractor` executable. Its
implementation is in `algorithms/TRUST4/FastqExtractor.cpp`, using
`IsGoodCandidate`, `ProcessReads_Thread`, and its paired `-1/-2` output.

For this particular subset:

| Statistic | Value |
|---|---:|
| Raw mate reads | 200,000 |
| Candidate mate reads | 200,000 |
| Candidate pairs | 100,000 |
| Candidate fraction | 1.0 |

Thus candidate extraction is now used, but it does not reduce this input.
This is an observed property of the subset, not an assumption that raw reads
are equivalent to candidates. Candidate provenance is recorded in
`candidate_statistics.tsv`.

## D. Graph and QAOS ablation

| Statistic | Identity-only | Identity + QAOS |
|---|---:|---:|
| Graph nodes/bundles | 59,148 | 59,148 |
| Candidate pairs | 233,810,604 | 233,810,604 |
| Accepted edges | 7,871,713 | 6,830,654 |
| Connected components | 1,998 | 2,370 |
| Isolated nodes | 620 | 620 |
| Branch events | 158,492 | 179,408 |
| Paths/contigs | 464,610 | 498,351 |
| N50 | 5,121 | 5,088 |

The combined score accepted 1,041,059 fewer edges than identity-only, a
13.23% reduction. Therefore QAOS now demonstrably affects edge decisions.
The combined run retained 4,859 edges below identity 0.90, showing that the
score is not equivalent to an identity-only gate.

Accepted-edge QAOS statistics in the combined run:

| Statistic | Value |
|---|---:|
| Mean overlap | 83.6795 bp |
| Mean identity | 0.974813 |
| Mean QAOS | 0.867360 |
| QAOS-passing edges (QAOS >= 0.90) | 3,102,387 |
| Low-identity accepted edges | 4,859 |

The edge score is a configurable weighted combination of identity, QAOS,
normalized overlap, and quality support. Identity and QAOS remain separate
columns in `graph_edges.tsv`.

## E. Biological reconstruction

| Measurement | Identity-only | Identity + QAOS |
|---|---:|---:|
| V assignments | 211,789 | 236,595 |
| J assignments | 55,792 | 41,520 |
| V+J assignments | 5,795 | 14,420 |
| CDR3/junction assignments | 31,644 | 28,454 |
| Complete V+J+CDR3 | 5,577 | 14,086 |
| Productive IGH | 416 | 1,166 |
| Unique clonotypes | 103 | 90 |
| Reference matches | 0 | 0 |
| Precision | 0 | 0 |
| Sensitivity | 0 | 0 |
| Pearson abundance correlation | NOT AVAILABLE | NOT AVAILABLE |

The combined score changes the biological yield materially: V+J and complete
V+J+CDR3 counts increase despite fewer accepted edges, while total CDR3 calls
decrease. This is an ablation observation, not a superiority claim.

## F. Same-subset TRUST4 comparison

The original TRUST4 subset produced 707 assembled records in its existing
run; its AIRR output contained 149 annotation records after the downstream
reporting path. The established subset clonotype comparison is retained as
previously computed.

| Metric | TRUST4 | Identity-only | Identity + QAOS |
|---|---:|---:|---:|
| Reconstructed contigs | 707 | 464,610 | 498,351 |
| V assignments | 149 | 211,789 | 236,595 |
| J assignments | 148 | 55,792 | 41,520 |
| CDR3s | 149 | 31,644 | 28,454 |
| Complete V+J+CDR3 | 148 | 5,577 | 14,086 |
| Unique clonotypes | 127 | 103 | 90 |
| Reference matches | 0 | 0 | 0 |
| Precision | 0 | 0 | 0 |
| Sensitivity | 0 | 0 | 0 |
| Abundance correlation | NOT AVAILABLE | NOT AVAILABLE | NOT AVAILABLE |
| Runtime | not captured | 1,734.02 s | 1,736.74 s |
| Peak graph memory | not captured | 982 MB | 887 MB |

No method is declared better or worse.

## G. Main failure point

The v2 result confirms that edge scoring can affect graph topology and
biological yield, but it has not yet implemented the intended BCR-aware path
resolver. Paths are still generated by bounded simple-path traversal. V/D/J/C
and CDR3 information is not propagated to nodes or used to prune/rescore
candidate paths before contig emission.

The principal remaining drop is therefore:

```text
graph paths
  -> annotation records with V/J/CDR3 compatibility
```

The combined run has 498,351 paths but only 14,086 complete V+J+CDR3
records. The downstream annotation is complete by record count; the loss is
biological compatibility, not missing output rows.

## H. Implementation verification

| Requirement | Current status |
|---|---|
| Existing TRUST4 candidate extraction | Integrated through `fastq-extractor`; subset fraction is 1.0 |
| Non-greedy graph | Bounded simple-path enumeration retains alternatives |
| QAOS overlap evaluation | Implemented and now affects edge acceptance |
| Separate identity | Preserved in edge/path output |
| Abundance soft evidence | Stored and accumulated as a path diagnostic, but not yet used for pruning or selection |
| Paired-end information | Binary exact bundle-pair endpoint support only |
| SHM-aware consensus | Not complete; consensus still chooses the higher-quality overlap base |
| V/D/J/C path resolution | Not implemented; annotation remains downstream |
| Original TRUST4 intact | Confirmed |
| Annotation memory <=32 GiB | Not met; annotator peaks around 36–37 GiB |

The v2 implementation is therefore a more informative graph prototype, not
yet the complete BCR-aware/SHM-aware design specified in the experiment.

## I. Recommendation

Do not start full FZ-116. The next step should be targeted debugging on this
same subset:

1. annotate graph nodes or candidate bundles before path resolution;
2. propagate V/J/C/locus evidence through paths;
3. use that evidence for dominance pruning rather than emitting all simple
   paths;
4. make abundance and paired consistency contribute to path ranking without
   deleting alternatives solely by abundance;
5. replace higher-quality-base consensus with explicit quality- and
   abundance-weighted base support;
6. reduce annotation memory through batching or path reduction.

The v2 run is computationally stable and QAOS is now demonstrably active, but
the full-scale readiness condition is not met.
