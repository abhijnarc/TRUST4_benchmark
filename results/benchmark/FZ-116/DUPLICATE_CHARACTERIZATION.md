# Duplicate Clonotype Characterization Report

**Date**: 2026-09-18  
**Purpose**: Understand duplicate structure before choosing deduplication strategy  
**Status**: Observations only - no deduplication decisions made

---

## Executive Summary

| Metric | iRepertoire | TRUST4 IGH |
|--------|-------------|-----------|
| Total records | 737,373 | 11,548 |
| Records with CDR3(pep) = '*' | 5,014 (0.7%) | 0 (0.0%) |
| Productive/valid records | 732,359 (99.3%) | 11,548 (100%) |
| Unique (CDR3aa+V+J+C) | 99,541 | 8,451 |
| Records in duplicate groups | 632,818 (86.4% of productive) | 3,097 (26.8% of total) |
| Duplicate groups | 39,572 | 651 |
| Reads in duplicate groups | N/A | 13,212 (29.1% of total) |

---

## iRepertoire Analysis

### A. Completeness

**Total records**: 737,373

**Records with CDR3(pep) = '*'** (stop codon representation):
- Count: 5,014
- Percentage: 0.7%
- Status: Present in data, exact meaning requires documentation review

**Missing fields**:
- V gene NULL: 0
- J gene NULL: 0
- C gene NULL: 0
- CDR3(pep) NULL/NaN: 0

**Interpretation**: All records have complete V, J, C fields. Stop codon records represent a small fraction (0.7%) of total.

---

### B. Productive Records Analysis (CDR3(pep) ≠ '*')

**Count**: 732,359 records (99.3% of total)

**Unique clonotypes** (defined by CDR3aa + V_gene + J_gene + C_gene):
- Count: 99,541 unique keys
- Productivity rate: 7.4 unique keys per record on average

**Duplicate records**:
- Count: 632,818 records participate in duplicate groups
- Percentage of productive: 86.4%
- Interpretation: Most iRepertoire records share identical (CDR3aa + V + J + C) with other records

---

### C. Duplicate Group Size Distribution (Productive Only)

**Number of duplicate groups**: 39,572

**Group sizes**:
- Minimum: 2 records
- Maximum: 4,532 records
- Mean: 16.99 records per group
- Median: 6 records per group

**Distribution** (showing representative sizes):

| Size | Count | Percentage |
|------|-------|-----------|
| 2 | 8,737 | 22.1% |
| 3 | 4,331 | 10.9% |
| 4 | 3,222 | 8.1% |
| 5 | 2,549 | 6.4% |
| 6 | 1,997 | 5.0% |
| 7-10 | 6,121 | 15.4% |
| 11-20 | 5,087 | 12.9% |
| 21-30 | 1,717 | 4.3% |
| 31-100 | 1,333 | 3.4% |
| 101-500 | 285 | 0.7% |
| 500+ | 1 | 0.0% |

**Observation**: Duplicate groups show a right-skewed distribution (most are small 2-6 records, few are very large).

The single group with 4,532 records represents a highly redundant clonotype (CDR3aa+V+J+C identical across 4,532 separate rows).

---

### D. What Fields Differ Within Duplicate Groups?

For each duplicate group, we examined whether records differ by:
- D gene
- CDR3 nucleotide sequence
- V/J allele (V and J with allele suffixes)
- Other fields

**Results** (among 39,572 productive duplicate groups):

| Difference Type | Count |
|-----------------|-------|
| Differ by D only | 3,325 |
| Differ by CDR3nt only | 12,768 |
| Differ by V/J allele only | 7,483 |
| Differ by multiple fields | 5,895 |
| Differ by other | 23,092 |

**Key observations**:
- **23,092 groups (58%)** differ by "other" fields - meaning records have identical CDR3aa+V+J+C but differ in fields beyond D, CDR3nt, and allele
- **12,768 groups** differ only in CDR3 nucleotide (but same amino acid) - this is expected (synonymous codons)
- **7,483 groups** differ only in V/J allele suffix
- **3,325 groups** differ only in D gene assignment
- **5,895 groups** differ by multiple fields simultaneously

**Interpretation**: The most common scenario (58%) is identical CDR3aa+V+J+C with differences elsewhere in the record.

---

### E. Representative Examples

**Example 1: "Other" field difference**
```
CDR3aa: *RARYFDWSGEDHYLFDY
V: IGHV1-18, J: IGHJ4, C: IGHA2

Record 1: D=IGHD3-9, CDR3nt=..., copy=1
Record 2: D=IGHD3-9, CDR3nt=..., copy=2  (identical CDR3nt)
Record 3: D=IGHD3-9, CDR3nt=..., copy=1  (identical CDR3nt)
```

Observation: All three records have same D, same CDR3nt, but are separate rows with different copy counts.

**Example 2: CDR3nt difference (synonymous)**
```
CDR3aa: CAAGFACMAW
V: IGHV3-23, J: IGHJ4, C: IGHA1

Record 1: D=IGHD2-8, CDR3nt=...GCTTGTAT..., count=3
Record 2: D=IGHD2-8, CDR3nt=...GCGTGCAT..., count=1  (1-2 nucleotide differences)
```

Observation: Amino acid identical, nucleotide sequences differ by 1-2 bases (synonymous variation or sequencing error).

---

### F. Copy Field Behavior

**Total copy counts**: 1,121,909 across all records  
**Total copy in productive**: 1,116,476

**If deduplicated by summing copy values**: 1,116,476 (same as current productive sum)

**Interpretation**: Summing copy counts within duplicate groups would preserve total abundance, but decisions about which records to keep or aggregate require explicit choice.

---

## TRUST4 Analysis

### A. Unique Clonotypes

**Total records**: 11,548 (all IGH)

**Unique clonotypes** (CDR3aa + V_gene + J_gene + C_gene):
- Count: 8,451 unique keys
- Duplicate records: 3,097 (26.8%)

**Interpretation**: Unlike iRepertoire, TRUST4 has fewer duplicates - only ~27% of records share the same benchmark key.

---

### B. Duplicate Groups

**Number of groups**: 651

**Group sizes**:
- Minimum: 2 records
- Maximum: 15 records
- Mean: 2.38 records per group
- Median: 2 records per group

**Distribution**:

| Size | Count | Percentage |
|------|-------|-----------|
| 2 | 506 | 77.7% |
| 3 | 105 | 16.1% |
| 4 | 21 | 3.2% |
| 5+ | 19 | 2.9% |

**Observation**: TRUST4 duplicates are predominantly small (77.7% have exactly 2 records).

---

### C. What Fields Differ Within Duplicate Groups?

| Difference Type | Count |
|-----------------|-------|
| Differ by D only | 59 |
| Differ by CDR3nt only | 594 |
| Differ by V/J allele only | 99 |
| Differ by multiple fields | 83 |
| Differ by other | 0 |

**Key observation**: TRUST4 duplicate groups differ mainly by:
- CDR3 nucleotide (594 groups / 91.2%)
- V/J allele (99 groups / 15.2%)
- D gene (59 groups / 9.1%)
- Multiple fields (83 groups / 12.8%)

**Notably**: No TRUST4 duplicates differ by "other" fields (unlike iRepertoire's 58%).

---

### D. Representative Examples

**Example 1: CDR3nt difference**
```
CDR3aa: CAAFSDRNFGDYRYYFEYW
V: IGHV1-58, J: IGHJ4, C: IGHA1

Record 1: D=IGHD4-17, CDR3nt=...TCGGATCGC..., count=2
Record 2: D=IGHD4-17, CDR3nt=...TCGGATCGC..., count=1  (very similar, 1-2 bases differ)
```

Observation: Same CDR3aa, same D, very similar CDR3nt, but separate read counts.

**Example 2: CDR3nt difference (different V assignment)**
```
CDR3aa: CAAGFACMAW
V: IGHV3-23, J: IGHJ4, C: IGHA1 (record 1)
V: IGHV3-53, J: IGHJ4, C: IGHA1 (record 2)

Both records have different V genes, indicating ambiguous V gene assignment.
```

---

### E. Read Count Analysis

**Total reads**: 45,406  
**Reads in duplicate groups**: 13,212 (29.1%)

**Interpretation**: About 3 in 10 reads belong to clonotypes that have multiple entries in the report. Deduplication decision will affect read count distribution.

---

## Implications for Benchmarking

### Deduplication Decisions Required

**iRepertoire challenge**: 86.4% of productive records are in duplicate groups
- Most duplicates share identical CDR3aa+V+J+C with different "other" fields
- Copy counts would need aggregation strategy (sum, max, mean?)
- Stop codon records (0.7%) need explicit handling

**TRUST4 challenge**: 26.8% of records are in duplicate groups
- Most duplicates differ by CDR3 nucleotide (somatic variants?)
- Some differ by V/J allele (ambiguous gene calls?)
- Read counts in duplicates represent 29.1% of total

### Primary Benchmark Key Definition

The benchmark key is defined as:
```
(CDR3aa_norm, V_gene_norm, J_gene_norm, C_gene_norm)
```

This key treats:
- **Identical** regardless of: D gene, CDR3nt, allele suffixes, other fields
- **Different** if: CDR3aa changes, or V/J/C gene identity changes (after normalization)

**This definition is used uniformly for both datasets.**

---

## Observations (Not Interpretations)

### What We Observe About iRepertoire

1. **Stop codon records**: 5,014 records (0.7%) have CDR3(pep) = '*'
   - These exist in the data
   - What they represent (non-productive rearrangements, sequencing errors, other) is undocumented

2. **Duplicate groups**: 39,572 groups with 632,818 records (86.4% of productive)
   - Most groups are small (median = 6 records)
   - One group has 4,532 records
   - Many differ only in fields beyond CDR3aa+V+J+C

3. **Copy field**: Summing within duplicate groups preserves total abundance

### What We Observe About TRUST4

1. **Duplicate groups**: 651 groups with 3,097 records (26.8% of total)
   - Most groups are very small (77.7% have exactly 2 records)
   - Predominantly differ by CDR3nt (91.2%) or V/J allele (15.2%)

2. **Read distribution**: 29.1% of reads are in duplicate groups
   - Aggregation will change abundance distribution

3. **No stop codons**: All TRUST4 records have valid CDR3aa

---

## Technically Defensible Deduplication Approaches

These approaches preserve scientific rigor without making premature interpretive claims:

### Approach 1: Aggregate by Summing Abundance

**iRepertoire**: Sum `copy` values for records with identical CDR3aa+V+J+C  
**TRUST4**: Sum `count` values for records with identical CDR3aa+V+J+C

**Rationale**:
- Unambiguous: every record contributes to exactly one clonotype
- Preserves total abundance
- Simple to implement
- Defensible as clonotype-level summary

**Consequence**:
- iRep: 732,359 → 99,541 unique clonotypes
- TRUST4: 11,548 → 8,451 unique clonotypes
- Loses information about nucleotide/allele variants

### Approach 2: Keep High-Count Variant Only

**iRepertoire**: For each key, keep record with highest `copy`  
**TRUST4**: For each key, keep record with highest `count`

**Rationale**:
- Unambiguous: one record per clonotype
- Assumes most-abundant variant is most representative
- Simple to implement

**Consequence**:
- Discards lower-abundance variants
- Risk of losing biologically meaningful somatic variants
- Same final size as Approach 1

### Approach 3: Filter Then Aggregate

**iRepertoire Step 1**: Remove records with CDR3(pep) = '*'  
**iRepertoire Step 2**: Aggregate remaining by CDR3aa+V+J+C

**TRUST4**: Aggregate as in Approach 1

**Rationale**:
- Removes stop-codon records explicitly
- Cleaner iRepertoire dataset for benchmarking
- Aggregate preserves abundance for valid sequences

**Consequence**:
- iRep: 737,373 → 732,359 (remove 5,014 stop codons) → 99,541 unique
- TRUST4: 11,548 → 8,451 unique
- Depends on whether stop codons should be in benchmark

### Approach 4: Keep Variants, Benchmark Separately

**Primary benchmark**: Use aggregated counts (Approach 1)  
**Secondary analysis**: Document how many clonotypes had variants, variant distribution

**Rationale**:
- Primary benchmark is unambiguous (aggregated clonotypes)
- Variant information is available for post-hoc analysis
- Allows investigating whether variants are important

**Consequence**:
- More work to implement
- Richer understanding of duplicate structure
- More defensible interpretation

---

## Recommendation Summary

**Before calculating metrics, choose:**

1. **iRepertoire treatment of stop codons**: Include or exclude from benchmark?
2. **Aggregation strategy**: Sum, max, or other?
3. **Reporting of variants**: Documented or implicit?

**No single approach is "correct" without domain context or publication precedent.**

The current data structure is clear: the duplication and aggregation decision should be explicit and documented, not silent.

---

## Data Files Available for Inspection

✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/trust4_normalized.tsv` (11,548 records)
✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/irep_normalized.tsv` (737,373 records)
✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/duplicate_characterization.tsv` (summary table)
✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/duplicate_characterization_detailed.tsv` (full metrics)

All normalized files preserve original data - no deduplication has been applied.

---

**Next step**: Choose deduplication strategy based on scientific/publication context, then update benchmarking script to implement the choice explicitly.
