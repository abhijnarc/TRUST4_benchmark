#!/usr/bin/env python3
"""
Rigorous benchmarking of TRUST4 FZ-116 against iRepertoire BCR-seq reference.

Methodology:
- Restrict to IGH (heavy chain) only
- Normalize gene names with documented rules
- Define TP/FP/FN based on CDR3aa + V + J + C agreement
- Check for duplicates before aggregation
- Separate diagnostic metrics from primary benchmark
- Abundance analysis as concordance study (not equivalence)
"""

import sys
import gzip
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.stats import spearmanr
import os

def load_trust4_report(filepath):
    """Load TRUST4 report.tsv"""
    print(f"Loading TRUST4 report from {filepath}...")
    # Note: Header line starts with # but is not a comment line
    # So we need to read without comment parameter
    df = pd.read_csv(filepath, sep='\t')
    # Fix header: remove # prefix if present
    df.columns = [col.lstrip('#') for col in df.columns]
    print(f"  Loaded {len(df)} records")
    print(f"  Columns: {list(df.columns)}")
    return df

def load_irep_reference(filepath):
    """Load iRepertoire gzipped CSV reference"""
    print(f"Loading iRepertoire reference from {filepath}...")
    df = pd.read_csv(filepath, sep=',', compression='gzip')
    print(f"  Loaded {len(df)} records")
    print(f"  Columns: {list(df.columns)}")
    return df

def extract_gene_name(gene_str):
    """
    Extract gene name without allele.
    E.g. IGHV3-23*05 -> IGHV3-23
         hIGHV3-23*05 -> IGHV3-23
    """
    if pd.isna(gene_str) or gene_str == '.':
        return None
    # Handle multiple genes separated by |
    if '|' in str(gene_str):
        genes = str(gene_str).split('|')
        return genes[0]  # Take first (highest ranked)
    gene = str(gene_str)
    # Remove 'h' prefix if present (iRepertoire convention)
    if gene.startswith('h'):
        gene = gene[1:]
    # Remove allele suffix (*01, *02, etc)
    if '*' in gene:
        gene = gene.split('*')[0]
    return gene

def normalize_c_gene(c_str):
    """
    Normalize C gene name.
    Remove only the 'h' prefix (iRepertoire convention).
    Preserve the actual C gene/isotype (e.g. IGHA1, IGHA2, IGHG1, etc).
    E.g. hIGHA1*01 -> IGHA1
         IGHA1 -> IGHA1
    """
    if pd.isna(c_str) or c_str == '.':
        return None
    c = str(c_str)
    # Remove 'h' prefix
    if c.startswith('h'):
        c = c[1:]
    # Remove allele suffix if present
    if '*' in c:
        c = c.split('*')[0]
    return c

def normalize_cdr3_aa(cdr3_str):
    """
    Normalize CDR3 amino acid sequence conservatively.
    Convert to uppercase, but preserve biological information.
    """
    if pd.isna(cdr3_str):
        return None
    return str(cdr3_str).upper()

def normalize_cdr3_nt(cdr3_str):
    """
    Normalize CDR3 nucleotide sequence conservatively.
    Convert to uppercase.
    """
    if pd.isna(cdr3_str):
        return None
    return str(cdr3_str).upper()

def filter_to_igh(df):
    """
    Filter to IGH (heavy chain) only.
    IGH genes have V genes starting with 'IGH'.
    """
    print(f"\nFiltering to IGH (heavy chain) only...")
    print(f"  Before filter: {len(df)} records")

    # Extract first V gene (highest ranked)
    df['v_first'] = df['V'].apply(lambda x: str(x).split('|')[0] if pd.notna(x) else None)
    igh_mask = df['v_first'].str.startswith('IGH', na=False)
    df_igh = df[igh_mask].copy()

    print(f"  After filter: {len(df_igh)} IGH records")

    # Show chain type distribution before filtering
    print(f"\n  Chain type distribution in original data:")
    chain_counts = df['v_first'].apply(
        lambda x: str(x)[:3] if pd.notna(x) else 'NA'
    ).value_counts()
    for chain, count in chain_counts.items():
        print(f"    {chain}: {count}")

    return df_igh.drop('v_first', axis=1)

def create_normalized_trust4(df_trust4):
    """
    Create normalized TRUST4 dataframe with additional columns.
    """
    print(f"\nNormalizing TRUST4 data...")

    df = df_trust4.copy()

    # Normalize CDR3 sequences
    df['cdr3aa_norm'] = df['CDR3aa'].apply(normalize_cdr3_aa)
    df['cdr3nt_norm'] = df['CDR3nt'].apply(normalize_cdr3_nt)

    # Normalize genes - preserve originals and create normalized versions
    df['v_gene_norm'] = df['V'].apply(extract_gene_name)
    df['j_gene_norm'] = df['J'].apply(extract_gene_name)
    df['c_gene_norm'] = df['C'].apply(normalize_c_gene)

    # Extract and normalize D gene
    df['d_gene_first'] = df['D'].apply(
        lambda x: (str(x).split('|')[0] if pd.notna(x) else None) if x != '.' else None
    )
    df['d_gene_norm'] = df['d_gene_first'].apply(extract_gene_name)

    print(f"  Added normalized columns: cdr3aa_norm, cdr3nt_norm, v_gene_norm, j_gene_norm, c_gene_norm, d_gene_norm")

    return df

def create_normalized_irep(df_irep):
    """
    Create normalized iRepertoire dataframe with additional columns.
    """
    print(f"\nNormalizing iRepertoire data...")

    df = df_irep.copy()

    # Rename columns for consistency (iRep uses different naming)
    df.rename(columns={'CDR3(pep)': 'CDR3aa', 'CDR3(nuc)': 'CDR3nt'}, inplace=True)

    # Normalize CDR3 sequences
    df['cdr3aa_norm'] = df['CDR3aa'].apply(normalize_cdr3_aa)
    df['cdr3nt_norm'] = df['CDR3nt'].apply(normalize_cdr3_nt)

    # Normalize genes
    df['v_gene_norm'] = df['V'].apply(extract_gene_name)
    df['j_gene_norm'] = df['J'].apply(extract_gene_name)
    df['c_gene_norm'] = df['C'].apply(normalize_c_gene)
    df['d_gene_norm'] = df['D'].apply(extract_gene_name)

    print(f"  Added normalized columns: cdr3aa_norm, cdr3nt_norm, v_gene_norm, j_gene_norm, c_gene_norm, d_gene_norm")

    return df

def check_duplicates(df, name, groupby_cols):
    """
    Check for duplicate clonotypes.
    """
    print(f"\nChecking for duplicates in {name}...")

    duplicate_mask = df.duplicated(subset=groupby_cols, keep=False)
    n_duplicates = duplicate_mask.sum()

    if n_duplicates > 0:
        print(f"  WARNING: {n_duplicates} records with duplicate {groupby_cols}")
        duplicates = df[duplicate_mask].sort_values(by=groupby_cols)
        print(f"\n  Duplicate groups (showing first 20):")

        grouped = duplicates.groupby(groupby_cols).size().reset_index(name='count')
        for idx, row in grouped.head(20).iterrows():
            print(f"    {dict(row)}")

        return duplicates
    else:
        print(f"  OK: No duplicates found for {groupby_cols}")
        return None

def build_matching_key(row):
    """
    Build the primary matching key for benchmarking.
    Key = CDR3aa + V + J + C (all normalized)
    """
    return (
        row['cdr3aa_norm'],
        row['v_gene_norm'],
        row['j_gene_norm'],
        row['c_gene_norm']
    )

def create_match_index(df, dataset_name):
    """
    Create index for fast matching.
    """
    print(f"\nCreating match index for {dataset_name}...")
    index = defaultdict(list)

    for idx, row in df.iterrows():
        key = build_matching_key(row)
        index[key].append(idx)

    print(f"  Created index with {len(index)} unique keys")

    # Check for keys with multiple entries
    multi_entry_keys = [k for k, v in index.items() if len(v) > 1]
    if multi_entry_keys:
        print(f"  WARNING: {len(multi_entry_keys)} keys with multiple entries")
        print(f"    (showing first 5)")
        for key in multi_entry_keys[:5]:
            print(f"      {key}: {len(index[key])} entries")

    return index

def perform_matching(df_trust4, df_irep, index_irep):
    """
    Match TRUST4 clonotypes to iRepertoire.
    """
    print(f"\nPerforming matching...")

    matches = []

    for idx, row in df_trust4.iterrows():
        key = build_matching_key(row)

        if key in index_irep:
            # Found a match
            irep_indices = index_irep[key]
            if len(irep_indices) == 1:
                irep_idx = irep_indices[0]
                irep_row = df_irep.iloc[irep_idx]

                matches.append({
                    'trust4_idx': idx,
                    'irep_idx': irep_idx,
                    'cdr3aa': row['cdr3aa_norm'],
                    'cdr3nt_trust4': row['cdr3nt_norm'],
                    'cdr3nt_irep': irep_row['cdr3nt_norm'],
                    'cdr3nt_match': row['cdr3nt_norm'] == irep_row['cdr3nt_norm'],
                    'v_gene': row['v_gene_norm'],
                    'd_gene_trust4': row['d_gene_norm'],
                    'd_gene_irep': irep_row['d_gene_norm'],
                    'd_gene_match': row['d_gene_norm'] == irep_row['d_gene_norm'],
                    'j_gene': row['j_gene_norm'],
                    'c_gene': row['c_gene_norm'],
                    'trust4_count': row['count'],
                    'trust4_frequency': row['frequency'],
                    'irep_copy': irep_row['copy'],
                    'match_status': 'TP',
                    'mismatch_reason': 'NONE'
                })
            else:
                # Multiple iRep entries with same key - record as ambiguous
                for irep_idx in irep_indices:
                    irep_row = df_irep.iloc[irep_idx]
                    matches.append({
                        'trust4_idx': idx,
                        'irep_idx': irep_idx,
                        'cdr3aa': row['cdr3aa_norm'],
                        'cdr3nt_trust4': row['cdr3nt_norm'],
                        'cdr3nt_irep': irep_row['cdr3nt_norm'],
                        'cdr3nt_match': row['cdr3nt_norm'] == irep_row['cdr3nt_norm'],
                        'v_gene': row['v_gene_norm'],
                        'd_gene_trust4': row['d_gene_norm'],
                        'd_gene_irep': irep_row['d_gene_norm'],
                        'd_gene_match': row['d_gene_norm'] == irep_row['d_gene_norm'],
                        'j_gene': row['j_gene_norm'],
                        'c_gene': row['c_gene_norm'],
                        'trust4_count': row['trust4_count'],
                        'trust4_frequency': row['trust4_frequency'],
                        'irep_copy': irep_row['copy'],
                        'match_status': 'AMBIGUOUS_MULTI_IREP',
                        'mismatch_reason': f'{len(irep_indices)}_irep_entries'
                    })

    df_matches = pd.DataFrame(matches)
    print(f"  Found {len(df_matches)} matches")

    return df_matches

def calculate_metrics(df_trust4_igh, df_irep, df_matches):
    """
    Calculate benchmarking metrics.
    """
    print(f"\nCalculating metrics...")

    n_trust4 = len(df_trust4_igh)
    n_irep = len(df_irep)
    n_unique_irep_keys = len(df_irep.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']).size())

    # Count match types
    tp = len(df_matches[df_matches['match_status'] == 'TP'])
    ambiguous = len(df_matches[df_matches['match_status'] == 'AMBIGUOUS_MULTI_IREP'])

    fp = n_trust4 - tp - ambiguous
    fn = n_irep - tp - ambiguous

    # Calculate metrics
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0.0

    metrics = {
        'n_trust4_igh': n_trust4,
        'n_irep_total': n_irep,
        'n_irep_unique_keys': n_unique_irep_keys,
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'ambiguous': ambiguous,
        'precision': precision,
        'sensitivity': sensitivity,
        'recall': sensitivity,
        'f1': f1
    }

    return metrics

def calculate_abundance_correlation(df_matches):
    """
    Calculate Spearman correlation for abundance measures.
    Only for TP matches (not ambiguous).
    """
    print(f"\nCalculating abundance correlation (TP matches only)...")

    df_tp = df_matches[df_matches['match_status'] == 'TP'].copy()

    if len(df_tp) < 2:
        print(f"  Insufficient TP matches for correlation ({len(df_tp)})")
        return None

    # Calculate correlations
    corr_count_vs_copy, pval_count_vs_copy = spearmanr(df_tp['trust4_count'], df_tp['irep_copy'])
    corr_freq_vs_copy, pval_freq_vs_copy = spearmanr(df_tp['trust4_frequency'], df_tp['irep_copy'])

    correlations = {
        'n_tp_pairs': len(df_tp),
        'trust4_count_vs_irep_copy_spearman': corr_count_vs_copy,
        'trust4_count_vs_irep_copy_pvalue': pval_count_vs_copy,
        'trust4_frequency_vs_irep_copy_spearman': corr_freq_vs_copy,
        'trust4_frequency_vs_irep_copy_pvalue': pval_freq_vs_copy
    }

    print(f"  Count vs Copy Spearman r={corr_count_vs_copy:.4f}, p={pval_count_vs_copy:.4e}")
    print(f"  Frequency vs Copy Spearman r={corr_freq_vs_copy:.4f}, p={pval_freq_vs_copy:.4e}")

    return correlations

def main():
    # Define paths
    base_dir = '/data1/wetlab/TRUST4_benchmark'
    trust4_file = f'{base_dir}/results/TRUST4/FZ-116/TRUST_FZ-116_report.tsv'
    irep_file = f'{base_dir}/reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz'

    output_dir = f'{base_dir}/results/benchmark/FZ-116'
    os.makedirs(output_dir, exist_ok=True)

    print("="*80)
    print("TRUST4 vs iRepertoire Benchmarking")
    print("="*80)

    # Load data
    df_trust4 = load_trust4_report(trust4_file)
    df_irep = load_irep_reference(irep_file)

    # Filter TRUST4 to IGH
    df_trust4_igh = filter_to_igh(df_trust4)

    # Normalize
    df_trust4_norm = create_normalized_trust4(df_trust4_igh)
    df_irep_norm = create_normalized_irep(df_irep)

    # Check for duplicates before proceeding
    print("\n" + "="*80)
    print("DUPLICATE ANALYSIS")
    print("="*80)

    trust4_dupes = check_duplicates(
        df_trust4_norm,
        'TRUST4',
        ['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']
    )

    irep_dupes = check_duplicates(
        df_irep_norm,
        'iRepertoire',
        ['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']
    )

    if trust4_dupes is not None or irep_dupes is not None:
        print("\n" + "!"*80)
        print("DUPLICATES DETECTED - STOPPING BEFORE METRIC CALCULATION")
        print("!"*80)
        print("\nAction required: Review duplicate structure and define aggregation strategy")

        # Save normalized data for inspection
        df_trust4_norm.to_csv(f'{output_dir}/trust4_normalized.tsv', sep='\t', index=False)
        df_irep_norm.to_csv(f'{output_dir}/irep_normalized.tsv', sep='\t', index=False)

        if trust4_dupes is not None:
            trust4_dupes.to_csv(f'{output_dir}/trust4_duplicates.tsv', sep='\t', index=False)
        if irep_dupes is not None:
            irep_dupes.to_csv(f'{output_dir}/irep_duplicates.tsv', sep='\t', index=False)

        print(f"\nNormalized data saved for inspection:")
        print(f"  {output_dir}/trust4_normalized.tsv")
        print(f"  {output_dir}/irep_normalized.tsv")
        if trust4_dupes is not None:
            print(f"  {output_dir}/trust4_duplicates.tsv ({len(trust4_dupes)} records)")
        if irep_dupes is not None:
            print(f"  {output_dir}/irep_duplicates.tsv ({len(irep_dupes)} records)")

        return 1

    # If no duplicates, proceed with matching
    print("\n" + "="*80)
    print("MATCHING AND METRICS")
    print("="*80)

    # Build match index
    index_irep = create_match_index(df_irep_norm, 'iRepertoire')

    # Perform matching
    df_matches = perform_matching(df_trust4_norm, df_irep_norm, index_irep)

    # Calculate metrics
    metrics = calculate_metrics(df_trust4_igh, df_irep_norm, df_matches)

    # Abundance correlation
    corr = calculate_abundance_correlation(df_matches)

    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)

    df_trust4_norm.to_csv(f'{output_dir}/trust4_normalized.tsv', sep='\t', index=False)
    df_irep_norm.to_csv(f'{output_dir}/irep_normalized.tsv', sep='\t', index=False)
    df_matches.to_csv(f'{output_dir}/matches.tsv', sep='\t', index=False)

    # Save metrics
    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(f'{output_dir}/metrics.tsv', sep='\t', index=False)

    # Save correlation
    if corr is not None:
        corr_df = pd.DataFrame([corr])
        corr_df.to_csv(f'{output_dir}/abundance_correlation.tsv', sep='\t', index=False)

    print(f"\nResults saved to {output_dir}/")
    print(f"  - trust4_normalized.tsv")
    print(f"  - irep_normalized.tsv")
    print(f"  - matches.tsv")
    print(f"  - metrics.tsv")
    if corr is not None:
        print(f"  - abundance_correlation.tsv")

    # Print summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"\nTRUST4 IGH records:        {metrics['n_trust4_igh']}")
    print(f"iRepertoire records:       {metrics['n_irep_total']}")
    print(f"iRep unique clonotypes:    {metrics['n_irep_unique_keys']}")
    print(f"\nTrue Positives:            {metrics['tp']}")
    print(f"False Positives:           {metrics['fp']}")
    print(f"False Negatives:           {metrics['fn']}")
    print(f"Ambiguous matches:         {metrics['ambiguous']}")
    print(f"\nPrecision:                 {metrics['precision']:.4f}")
    print(f"Sensitivity/Recall:        {metrics['sensitivity']:.4f}")
    print(f"F1-score:                  {metrics['f1']:.4f}")

    if corr is not None:
        print(f"\nAbundance correlation (TP only):")
        print(f"  Count vs Copy Spearman:  {corr['trust4_count_vs_irep_copy_spearman']:.4f}")
        print(f"  Frequency vs Copy Spearman: {corr['trust4_frequency_vs_irep_copy_spearman']:.4f}")

    return 0

if __name__ == '__main__':
    sys.exit(main())
