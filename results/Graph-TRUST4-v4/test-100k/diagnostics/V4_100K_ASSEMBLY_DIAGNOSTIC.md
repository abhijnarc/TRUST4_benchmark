# Graph-TRUST4-v4 100K Assembly Diagnostic

## Dataset
- Raw reads: 200000
- Bundles: 59,148
- Files used: bundles, graph edges/statistics, branch decisions, assembled FASTA, path/assembly debug and run metadata in this directory; raw FASTQ absent.
- Date/time: NOT AVAILABLE in metadata.

## Computational profile
- Candidate pairs: 29573809
- Accepted edges: 6,141,188
- Peak RSS: 666 MB
- Threads: 64
- Wall time, CPU time, CPU utilization, swap: NOT AVAILABLE (no in-directory timing/log fields).

## Graph structure
- Components: 6,903 weak components; isolated nodes: 5,494
- Branches/resolved/ambiguous: 3271559/1003314/2268245
- Cycles: 0
- Accepted edges / candidate pairs: 0.2076563083; edges per bundle: 103.8274836

## Assembly structure
- Contigs: 131,032; min/Q1/median/mean/Q3/max bp: 150 / 237 / 348 / 382.296149 / 469 / 1621
- N50/N90: 451/221 bp; total bases: 50,093,029
- Nodes/path min/median/mean/max: 2/21/48.03204561/122; edges/path NOT AVAILABLE.

## Branch resolution
- Resolution/ambiguity rates: 0.3066776421/0.6933223579
- QAOS_delta n/min/Q1/median/mean/Q3/max: 3271559 / 0 / 0.00120467 / 0.00449699 / 0.007802919134 / 0.0124659 / 0.0993344
- See BRANCH_DECISION_SUMMARY.md and QAOS_DELTA_DISTRIBUTION.tsv for categories and bins.

## Edge quality
- overlap_length: n=6,141,188, min=31, Q1=97, median=117, mean=113.9743543, Q3=135, max=150.
- identity: n=6,141,188, min=0.933333, Q1=0.984615, median=0.989796, mean=0.9892963483, Q3=0.992908, max=1.
- qaos: n=6,141,188, min=0.9, Q1=0.920419, median=0.942705, mean=0.9476431199, Q3=0.969399, max=0.999831.
- high_quality_mismatches: n=6,141,188, min=0, Q1=0, median=0, mean=0.2574790089, Q3=1, max=2.
- paired_support: n=6,141,188, min=0, Q1=0, median=0, mean=0, Q3=0, max=0.

## V3-debug comparison

| metric | v3_debug | v4 | absolute_change | percent_change | interpretation |
| raw_reads | NOT AVAILABLE | 200000 | NOT AVAILABLE | NOT AVAILABLE | NOT AVAILABLE |
| bundles | NOT AVAILABLE | 59148 | NOT AVAILABLE | NOT AVAILABLE | NOT AVAILABLE |
| candidate_pairs | NOT AVAILABLE | 29573809 | NOT AVAILABLE | NOT AVAILABLE | NOT AVAILABLE |
| accepted_edges | 270989 | 6141188 | 5870199 | 2166.21302 | increased by 5870199 |
| components | 46495 | 6903 | -39592 | -85.15324228 | decreased by 39592 |
| isolated_nodes | 45790 | 5494 | -40296 | -88.00174711 | decreased by 40296 |
| cycles | 59 | 0 | -59 | -100 | decreased by 59 |
| peak_rss_mb | 138 | 666 | 528 | 382.6086957 | increased by 528 |
| contigs | 2909 | 131032 | 128123 | 4404.365761 | increased by 128123 |
| min_contig_length | 150 | 150 | 0 | 0 | unchanged at 150 |
| median_contig_length | 179 | 348 | 169 | 94.41340782 | increased by 169 |
| mean_contig_length | 191.56 | 382.296149 | 190.736149 | 99.56992537 | increased by 190.736149 |
| N50 | 191 | 451 | 260 | 136.1256545 | increased by 260 |
| max_contig_length | 394 | 1621 | 1227 | 311.4213198 | increased by 1227 |
| contigs_gt_300 | 10 | 77880 | 77870 | 778700 | increased by 77870 |
| contigs_gt_500 | 0 | 23355 | 23355 | NOT AVAILABLE | increased by 23355 |
| contigs_gt_1000 | 0 | 408 | 408 | NOT AVAILABLE | increased by 408 |
| repeated_node_paths | 0 | 0 | 0 | NOT AVAILABLE | unchanged at 0 |
| length_mismatches | 0 | 0 | 0 | NOT AVAILABLE | unchanged at 0 |

## Potential issues
- Path-edge quality, orientation transitions, overlap validity, and branch membership are absent from supplied path tables; not inferred.
- Timing/swap data and a log are absent. Length/node outliers are statistical candidates, not evidence of biological error.

## Evidence needed before full FZ-116
- Capture repeatable wall/CPU time, peak memory, and swap.
- Validate orientation, overlap, cycles, repeated nodes, and expected lengths with explicit path-edge membership.
- Inspect flagged length/node outliers and their edge evidence.
- Review branch and edge distributions against thresholds, and verify resource headroom/output completeness before scaling.

No biological benchmarking or iRepertoire evaluation was performed; computational statistics alone do not establish biological improvement.
