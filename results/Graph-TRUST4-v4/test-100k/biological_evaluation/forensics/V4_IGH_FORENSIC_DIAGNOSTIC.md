# Graph-TRUST4-v4 IGH Forensic Diagnostic

## 1. Observed problem
The assembly has 131,032 contigs, but AIRR reports 712 IGH V-hit contigs, 0 with J, 0 with C, 1 with a formal junction/CDR3, and 1 productive IGH calls. Primary V-J-C/CDR3 completeness is zero.

## 2. V-hit sequence characteristics
The 712 V-hit contigs have length min/median/mean/N50/max 150/246.0/274.94/273/556 bp; all v4 contigs have median 348.0 bp, mean 382.30, N50 451, max 1621.
TRUST4 CIGAR query coordinates show V-alignment length min/median/mean/max 18/24.0/25.70/56 nt; identity min/median/mean/max 71.43/85.19/84.99/95.65%. CIGAR query lengths reconcile for all 712 records.
V-hit length-bin counts: 150-199=178, 200-249=199, 250-299=146, 300-399=79, 400-499=85, 500-749=25, 750-999=0, 1000-1499=0, >=1500=0.

## 3. V-to-J/CDR3 continuity
The median query sequence after the end of the V alignment is 138 bp (range 0-397); downstream-bin counts: <50=120, 50-99=71, 100-199=392, 200-299=32, 300-499=97, >=500=0. 120 V hits have <50 bp downstream, while 129 have >=200 bp downstream but no formal J call. This demonstrates extension beyond the short V alignment in many records, but does not prove that the extension is IGH V-D-J context.
No IGH J alignment or partial-junction motif was independently established. The same TRUST4 annotator/reference produced no IGH J calls, including on reverse-complement copies; no local BLAST/minimap2/vsearch executable and no exact path membership are available, so a separate IGHJ sequence search is NOT AVAILABLE rather than guessed. CDR3 context categories: FORMAL_CDR3=1 (0.14%); UNCERTAIN=711 (99.86%); partial/no-junction categories are not assigned without a reliable J alignment.

## 4. Graph/path evidence
`PATH_MEMBERSHIP_AVAILABLE=false`. `path_statistics.tsv` stores contig length, node count, validity, and aggregate flags but no node IDs. `graph_edges.tsv` identifies global source/target and edge metrics; `branch_decisions.tsv` records global branch/edge IDs but not the contig path that used them. Therefore no exact IGH path nodes/edges, branch counts, branch-specific QAOS, or IGH-path edge-quality distribution can be assigned. No heuristic path reconstruction was used.
Across all graph paths, v4 has 131032/131032 valid, 0 repeated-node, 0 cycle, and 0 length-mismatch paths. Global branch totals are 3271559 (1003314 resolved; 2268245 ambiguous; 69.33% ambiguous); median QAOS_delta=0.00449699. Branch IDs cannot be tied to the 712 V-hit paths. These global checks do not establish biological correctness.

## 5. Edge-quality evidence
For all accepted edges (not specifically IGH paths): overlap mean/median/min 113.974/117/31 bp; identity 0.989296/0.989796/0.933333; QAOS 0.947643/0.942705/0.900000. IGH-path edge metrics and overlap bins are NOT AVAILABLE because path edge IDs were not saved. Global edge quality cannot certify IGH joins.

## 6. Orientation evidence
All 712 V-hit AIRR records have `rev_comp=F`. Re-annotating only their reverse-complement sequences with the same TRUST4 annotator/reference yielded 0 J calls, 0 C calls, and 1 formal junction; 0 contigs gained J/CDR3 evidence. Thus reverse-complement rescue is NOT SUPPORTED as the main explanation. All 6,141,188 accepted graph edges have orientation +1 because the v4 scorer hard-codes that value; path-level orientation continuity cannot be evaluated without membership.

## 7. Paired-end evidence
Status: NOT_RETAINED. The FASTQ reader discards record names; bundling keys reads by sequence and keeps abundance/quality but no mate IDs. Every saved edge and branch support value is zero. Paired support appears as a tie-breaker in branch ordering but cannot affect decisions when all values are zero. This is a real loss of evidence, not proof of causation. See `PAIRED_END_SUPPORT_AUDIT.md`.

## 8. Comparison with v3
v3-debug 100K: N50=191, max=394, valid=2850/2909, repeated paths=0, cycles=59, length mismatches=0; 217 contig lengths are exact multiples of 150.
v4 100K: N50=451, max=1621, valid=131032/131032, repeated paths=0, cycles=0, length mismatches=0; 1113 exact multiples of 150. This supports that v4's stored path/overlap length checks are structurally cleaner and the multiple-of-150 pattern is much rarer proportionally (7.46% vs 0.85%); structural validity is not biological completeness.
BIOLOGICAL_COMPLETENESS_PROBLEM=YES (observed missing IGH J/C/CDR3 completeness); mechanism is not fully resolved.

Top 100 longest contigs: IGH=0, with V=0, J=100, CDR3=0, V+J=0, V+J+CDR3=0, complete V+J+C+CDR3=0; chain annotations={'TRB': 100}. All 100 are annotated TRB, with J calls TRBJ2-6 and C calls TRBC2 (98) or TRBC1 (2), and no V/CDR3 calls. This is evidence that the longest outputs are T-cell-receptor-associated contigs, not evidence that they are complete IGH assemblies or chimeras.

## 9. Most strongly supported failure mode
- A. Paths too short: SUPPORTED for the V-hit subset (median length 246 bp, N50 273 bp; 120 have <50 bp after the short V alignment), but not a complete explanation because many have longer downstream sequence.
- B. Unrelated sequence merging: UNCLEAR; no concrete sequence/path evidence of a chimera was found, and exact paths are unavailable.
- C. Incorrect branch resolution: UNCLEAR; branch decisions cannot be tied to these contigs.
- D. Premature termination before V-D-J: UNCLEAR; the contigs lack formal J/CDR3 calls, but long downstream sequence exists and independent J sequence evidence is unavailable.
- E. J/CDR3 annotator failure despite sequence presence: UNCLEAR; no independent IGHJ alignment established sequence presence.
- F. Orientation problem: NOT SUPPORTED as the main cause; all 712 annotator flags are forward and reverse-complement re-annotation did not restore J/C or additional junctions. Graph edge orientation is still hard-coded +1.
- G. Paired-end problem: SUPPORTED as a missing source of path evidence (mate IDs discarded; support always zero); causal link is unproven.
- H. Candidate-edge problem: UNCLEAR; all-edge scores are strong but cannot be associated with IGH paths.
- I. Wrong branch/path extraction: UNCLEAR; exact path membership was not emitted.
- Other: SUPPORTED observation that V hits are short, low-identity partial alignments (median 85.19% identity over median 24 nt), and longer contigs are not translating into IGH complete calls.

## 10. Recommended next algorithmic investigation
First preserve exact ordered node IDs, edge IDs, orientations, branch decisions, and abundance for each contig. Then use those exact paths to align V-hit context against the same TRUST4 IMGT reference in both orientations and distinguish missing J sequence from missed annotation. Retain mate identifiers through bundling and measure mate support at candidate edges. Inspect long non-IGH contigs and V-hit downstream sequence for architecture discontinuities only after path membership exists. These are investigation recommendations only; no algorithm or threshold was changed.
