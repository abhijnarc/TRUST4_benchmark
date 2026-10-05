#!/usr/bin/env python3
"""
Detailed characterization of duplicate clonotypes in iRepertoire and TRUST4.

Purpose: Understand data structure before deciding on deduplication strategy.

Does NOT:
- Make interpretation claims
- Calculate benchmark metrics
- Deduplicate or modify data
- Make decisions

Only describes observations about duplicate structure.
"""

import sys
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import gzip

def load_files():
    """Load normalized datasets"""
    print("Loading normalized datasets...")

    trust4_file = '/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/trust4_normalized.tsv'
    irep_file = '/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/irep_normalized.tsv'

    df_t4 = pd.read_csv(trust4_file, sep='\t')
    df_ir = pd.read_csv(irep_file, sep='\t')

    print(f"  TRUST4: {len(df_t4)} records")
    print(f"  iRepertoire: {len(df_ir)} records")

    return df_t4, df_ir

def analyze_irep(df_ir):
    """Detailed analysis of iRepertoire data"""

    print("\n" + "="*80)
    print("iRepertoire ANALYSIS")
    print("="*80)

    results = {}

    # A. Total records
    n_total = len(df_ir)
    print(f"\nA. TOTAL RECORDS: {n_total}")
    results['n_total_irep'] = n_total

    # B. Missing/empty records
    print(f"\nB. COMPLETENESS:")

    n_stop_codon = (df_ir['cdr3aa_norm'] == '*').sum()
    pct_stop = 100 * n_stop_codon / n_total
    print(f"   CDR3(pep) == '*' (stop codon): {n_stop_codon} ({pct_stop:.1f}%)")
    results['n_stop_codon'] = n_stop_codon

    n_missing_cdr3 = df_ir['cdr3aa_norm'].isna().sum()
    print(f"   CDR3(pep) == NULL/NaN: {n_missing_cdr3}")
    results['n_missing_cdr3'] = n_missing_cdr3

    n_missing_v = df_ir['v_gene_norm'].isna().sum()
    print(f"   V gene == NULL: {n_missing_v}")
    results['n_missing_v'] = n_missing_v

    n_missing_j = df_ir['j_gene_norm'].isna().sum()
    print(f"   J gene == NULL: {n_missing_j}")
    results['n_missing_j'] = n_missing_j

    n_missing_c = df_ir['c_gene_norm'].isna().sum()
    print(f"   C gene == NULL: {n_missing_c}")
    results['n_missing_c'] = n_missing_c

    # C. Productive (non-stop) analysis
    print(f"\nC. AMONG NON-STOP RECORDS:")

    df_prod = df_ir[df_ir['cdr3aa_norm'] != '*'].copy()
    n_productive = len(df_prod)
    print(f"   Total productive records: {n_productive}")
    results['n_productive_irep'] = n_productive

    # Unique keys in productive
    n_unique_keys = df_prod.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']).ngroups
    print(f"   Unique (CDR3aa + V + J + C) keys: {n_unique_keys}")
    results['n_unique_keys_prod'] = n_unique_keys

    n_duplicate_groups = len(df_prod) - n_unique_keys
    print(f"   Records in duplicate groups: {n_duplicate_groups}")
    results['n_dup_records_prod'] = n_duplicate_groups

    # Distribution of duplicate group sizes
    print(f"\n   DUPLICATE GROUP SIZE DISTRIBUTION (productive only):")
    dup_sizes = df_prod.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']).size()
    dup_only = dup_sizes[dup_sizes > 1]

    print(f"   - Number of duplicate groups: {len(dup_only)}")
    results['n_dup_groups_prod'] = len(dup_only)

    if len(dup_only) > 0:
        print(f"   - Min records per group: {dup_only.min()}")
        print(f"   - Max records per group: {dup_only.max()}")
        print(f"   - Mean records per group: {dup_only.mean():.2f}")
        print(f"   - Median records per group: {dup_only.median():.1f}")

        # Distribution
        print(f"\n   DISTRIBUTION:")
        for size in sorted(dup_only.unique()):
            count = (dup_only == size).sum()
            pct = 100 * count / len(dup_only)
            print(f"     {size} records per group: {count} groups ({pct:.1f}%)")

        results['dup_size_min'] = int(dup_only.min())
        results['dup_size_max'] = int(dup_only.max())
        results['dup_size_mean'] = float(dup_only.mean())

    # What fields differ in duplicates?
    print(f"\n   FIELD DIFFERENCES IN DUPLICATE GROUPS (productive):")

    differ_by_d = 0
    differ_by_nt = 0
    differ_by_allele = 0
    differ_by_multi = 0
    differ_by_other = 0

    for key, group in df_prod.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']):
        if len(group) == 1:
            continue

        # Check what differs
        d_values = group['d_gene_norm'].nunique()
        nt_values = group['cdr3nt_norm'].nunique()
        allele_v = group['V'].nunique()  # Original with allele
        allele_j = group['J'].nunique()  # Original with allele

        diff_count = 0
        diff_fields = []

        if d_values > 1:
            differ_by_d += 1
            diff_count += 1
            diff_fields.append('D')

        if nt_values > 1:
            differ_by_nt += 1
            diff_count += 1
            diff_fields.append('CDR3nt')

        if allele_v > 1 or allele_j > 1:
            differ_by_allele += 1
            diff_count += 1
            diff_fields.append('allele')

        if diff_count > 1:
            differ_by_multi += 1
        elif diff_count == 0:
            differ_by_other += 1

    print(f"     Differ by D only: {differ_by_d}")
    print(f"     Differ by CDR3nt only: {differ_by_nt}")
    print(f"     Differ by V/J allele only: {differ_by_allele}")
    print(f"     Differ by multiple fields: {differ_by_multi}")
    print(f"     Differ by other: {differ_by_other}")

    results['irep_differ_d'] = differ_by_d
    results['irep_differ_nt'] = differ_by_nt
    results['irep_differ_allele'] = differ_by_allele
    results['irep_differ_multi'] = differ_by_multi

    # D. Representative examples
    print(f"\n   REPRESENTATIVE EXAMPLES (productive duplicates):")
    shown = 0
    for key, group in df_prod.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']):
        if len(group) > 1 and shown < 3:
            print(f"\n     Example {shown + 1}:")
            print(f"       Key: CDR3aa={key[0]}, V={key[1]}, J={key[2]}, C={key[3]}")
            print(f"       Records: {len(group)}")
            for idx, row in group.iterrows():
                d_val = row['d_gene_norm'] if pd.notna(row['d_gene_norm']) else 'NA'
                nt = row['cdr3nt_norm'][:20] + '...' if len(str(row['cdr3nt_norm'])) > 20 else row['cdr3nt_norm']
                print(f"         - D={d_val}, CDR3nt={nt}, copy={row['copy']}")
            shown += 1

    # E. Copy field consideration
    print(f"\nE. COPY FIELD ANALYSIS:")
    print(f"   Total copy counts across all records: {df_ir['copy'].sum()}")
    print(f"   Total copy counts in productive: {df_prod['copy'].sum()}")

    # If duplicates exist, would summing change totals?
    copy_by_key = df_prod.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm'])['copy'].sum()
    copy_summed = copy_by_key.sum()
    print(f"   If deduplicated by summing: total={copy_summed}")
    print(f"   (Should equal productive total: {copy_summed == df_prod['copy'].sum()})")

    return results

def analyze_trust4(df_t4):
    """Detailed analysis of TRUST4 data"""

    print("\n" + "="*80)
    print("TRUST4 ANALYSIS")
    print("="*80)

    results = {}

    # A. Total and unique keys
    print(f"\nA. UNIQUE CLONOTYPES:")
    n_total = len(df_t4)
    print(f"   Total records: {n_total}")
    results['n_total_trust4'] = n_total

    n_unique = df_t4.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']).ngroups
    print(f"   Unique (CDR3aa + V + J + C): {n_unique}")
    results['n_unique_trust4'] = n_unique

    # B. Duplicate groups
    print(f"\nB. DUPLICATE GROUPS:")
    n_in_duplicates = n_total - n_unique
    print(f"   Records in duplicate groups: {n_in_duplicates}")
    results['n_dup_records_trust4'] = n_in_duplicates

    # Distribution
    dup_sizes = df_t4.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']).size()
    dup_only = dup_sizes[dup_sizes > 1]

    print(f"   Number of duplicate groups: {len(dup_only)}")
    results['n_dup_groups_trust4'] = len(dup_only)

    if len(dup_only) > 0:
        print(f"   Min records per group: {dup_only.min()}")
        print(f"   Max records per group: {dup_only.max()}")
        print(f"   Mean records per group: {dup_only.mean():.2f}")
        print(f"   Median records per group: {dup_only.median():.1f}")

        print(f"\n   DISTRIBUTION:")
        for size in sorted(dup_only.unique()):
            count = (dup_only == size).sum()
            pct = 100 * count / len(dup_only)
            print(f"     {size} records per group: {count} groups ({pct:.1f}%)")

    # C. What differs in duplicates?
    print(f"\nC. FIELD DIFFERENCES IN DUPLICATE GROUPS:")

    differ_by_d = 0
    differ_by_nt = 0
    differ_by_allele = 0
    differ_by_multi = 0
    differ_by_other = 0

    for key, group in df_t4.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']):
        if len(group) == 1:
            continue

        d_values = group['d_gene_norm'].nunique()
        nt_values = group['cdr3nt_norm'].nunique()
        allele_v = group['V'].nunique()  # Original with allele
        allele_j = group['J'].nunique()  # Original with allele

        diff_count = 0
        diff_fields = []

        if d_values > 1:
            differ_by_d += 1
            diff_count += 1
            diff_fields.append('D')

        if nt_values > 1:
            differ_by_nt += 1
            diff_count += 1
            diff_fields.append('CDR3nt')

        if allele_v > 1 or allele_j > 1:
            differ_by_allele += 1
            diff_count += 1
            diff_fields.append('allele')

        if diff_count > 1:
            differ_by_multi += 1
        elif diff_count == 0:
            differ_by_other += 1

    print(f"   Differ by D only: {differ_by_d}")
    print(f"   Differ by CDR3nt only: {differ_by_nt}")
    print(f"   Differ by V/J allele only: {differ_by_allele}")
    print(f"   Differ by multiple fields: {differ_by_multi}")
    print(f"   Differ by other: {differ_by_other}")

    results['trust4_differ_d'] = differ_by_d
    results['trust4_differ_nt'] = differ_by_nt
    results['trust4_differ_allele'] = differ_by_allele
    results['trust4_differ_multi'] = differ_by_multi

    # D. Representative examples
    print(f"\nD. REPRESENTATIVE EXAMPLES (duplicates):")
    shown = 0
    for key, group in df_t4.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']):
        if len(group) > 1 and shown < 3:
            print(f"\n   Example {shown + 1}:")
            print(f"     Key: CDR3aa={key[0]}, V={key[1]}, J={key[2]}, C={key[3]}")
            print(f"     Records: {len(group)}")
            for idx, row in group.iterrows():
                d_val = row['d_gene_norm'] if pd.notna(row['d_gene_norm']) else 'NA'
                nt = row['cdr3nt_norm'][:20] + '...' if len(str(row['cdr3nt_norm'])) > 20 else row['cdr3nt_norm']
                print(f"       - D={d_val}, CDR3nt={nt}, count={row['count']}")
            shown += 1

    # E. Read counts in duplicates
    print(f"\nE. READ COUNT ANALYSIS:")
    total_reads = df_t4['count'].sum()
    print(f"   Total reads across all records: {total_reads}")

    dup_records = df_t4.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']).size()
    dup_mask = dup_records > 1

    reads_in_dups = 0
    for key, group in df_t4.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']):
        if len(group) > 1:
            reads_in_dups += group['count'].sum()

    pct_reads_in_dups = 100 * reads_in_dups / total_reads if total_reads > 0 else 0
    print(f"   Reads in duplicate groups: {reads_in_dups} ({pct_reads_in_dups:.1f}%)")
    results['trust4_reads_in_dups'] = reads_in_dups
    results['trust4_pct_reads_in_dups'] = pct_reads_in_dups

    return results

def create_summary_table(irep_results, trust4_results):
    """Create summary table"""

    summary = pd.DataFrame([
        {
            'Dataset': 'iRepertoire',
            'Total Records': irep_results['n_total_irep'],
            'Stop Codon Records': irep_results['n_stop_codon'],
            'Productive Records': irep_results['n_productive_irep'],
            'Unique Keys': irep_results['n_unique_keys_prod'],
            'Duplicate Groups': irep_results['n_dup_groups_prod'],
            'Records in Duplicates': irep_results['n_dup_records_prod']
        },
        {
            'Dataset': 'TRUST4 IGH',
            'Total Records': trust4_results['n_total_trust4'],
            'Stop Codon Records': 0,
            'Productive Records': trust4_results['n_total_trust4'],
            'Unique Keys': trust4_results['n_unique_trust4'],
            'Duplicate Groups': trust4_results['n_dup_groups_trust4'],
            'Records in Duplicates': trust4_results['n_dup_records_trust4']
        }
    ])

    return summary

def main():
    print("="*80)
    print("DUPLICATE CHARACTERIZATION (No Deduplication)")
    print("="*80)

    # Load
    df_t4, df_ir = load_files()

    # Analyze
    irep_results = analyze_irep(df_ir)
    trust4_results = analyze_trust4(df_t4)

    # Summary table
    summary = create_summary_table(irep_results, trust4_results)

    print("\n" + "="*80)
    print("SUMMARY TABLE")
    print("="*80)
    print(summary.to_string(index=False))

    # Save summary table
    output_dir = '/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116'
    summary.to_csv(f'{output_dir}/duplicate_characterization.tsv', sep='\t', index=False)
    print(f"\nSummary saved to: {output_dir}/duplicate_characterization.tsv")

    # Save full results dict
    results_combined = {**irep_results, **trust4_results}
    results_df = pd.DataFrame([results_combined])
    results_df.to_csv(f'{output_dir}/duplicate_characterization_detailed.tsv', sep='\t', index=False)

    return irep_results, trust4_results, summary

if __name__ == '__main__':
    irep_res, t4_res, summary = main()
