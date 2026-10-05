# Graph-TRUST4 FZ-116 Biological Output Diagnostic

This diagnostic reads completed Graph-TRUST4 and TRUST4-annotator outputs only. It does not rerun or modify assembly.

## Annotation yield

The annotator produced 5679 rows for 5679 FASTA contigs. The status distribution is: C=2854, J=5, J+C=183, J+C+CDR3=200, V=928, V+C=71, V+C+CDR3=3, V+CDR3=114, V+J=2, V+J+C=17, V+J+C+CDR3=549, V+J+CDR3=190, none=563.

Representative records are in `annotation_diagnostic.tsv`: 20 rows each for valid_IGH, missing_V, missing_J, and missing_CDR3. `v_region_status=annotator_V_hit` means the raw annotator supplied a V call or V CIGAR; it is not an inferred biological rescue.

## CDR3 relationships

- V but no CDR3: 1018
- J but no CDR3: 207
- V+J but no CDR3: 19
- CDR3 but no V: 200
- CDR3 but no J: 117
- V+J+CDR3: 739
- V+J+C+CDR3: 549

These combinations show whether CDR3 loss is caused by missing V/J calls or by partial sequences that do not contain a detectable junction.

## Abundance propagation

`abundance_diagnostic.tsv` reconstructs the existing path selection from `graph_edges.tsv`, looks up each node in `bundles.tsv`, and compares the sum of bundle abundances with the current normalized contig count. `number_of_supporting_reads_if_available` is the same sum; no independent read IDs are retained by the assembler.

The first 20 paths were traced. Their summed path abundances have min=2, median=2.0, mean=2.15, max=3.

## Findings

- Most contigs are short partial sequences (the completed FASTA is 150-376 bp), so missing V/J/CDR3 calls are consistent with insufficient germline/junction context rather than a TSV field-loss problem.
- The raw annotation table contains the fields directly; this does not support an annotation conversion/parser loss.
- The current path abundance is the sum of bundle abundances along the path, and the diagnostic compares that value with the normalized count.
- The corrected assembly validation had zero repeated-node paths and zero expected-vs-actual length mismatches; the 25 cycle flags remain explicitly recorded.

=== ROOT CAUSE SUMMARY ===

IGH yield problem: CONFIRMED - low yield is driven by short partial contigs and incomplete V calls; the raw annotator output contains the same missing calls, with no evidence of FASTA/header loss.

CDR3 yield problem: CONFIRMED - most contigs lack a detectable junction, primarily because they are short partial reconstructions; V/J/CDR3 combination counts quantify the affected groups.

Abundance problem: NOT SUPPORTED - existing path counts equal the sum of bundle abundances for the traced contigs; the distribution is compressed because most paths consist of low-abundance bundles, not because counts are silently replaced by one node.

Annotation/parser problem: NOT SUPPORTED - TRUST4 annotator output is present for all FASTA records and missing fields correspond to partial annotation results, not conversion loss.

Assembly problem: NOT SUPPORTED - path validation shows zero repeated-node paths and zero length mismatches; 25 cycle paths are flagged/broken as designed.
