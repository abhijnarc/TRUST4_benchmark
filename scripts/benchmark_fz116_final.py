#!/usr/bin/env python3
"""
TRUST4 vs iRepertoire Benchmarking - Final Implementation

Methodology:
- Exclude iRep CDR3='*' from primary benchmark (documented separately)
- Aggregate duplicates by CDR3aa+V+J+C (sum counts/copy)
- Match on CDR3aa+V+J+C only (D evaluated separately)
- Calculate precision/sensitivity/F1
- Diagnostic metrics for gene agreement and variants
- Abundance concordance analysis

Does NOT:
- Modify source files
- Make silent discards
- Claim biological causation
- Run other samples
"""

import sys
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import os
from pathlib import Path

def load_normalized_data():
    """Load the normalized datasets"""
    print("Loading normalized datasets...")

    trust4_file = '/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/trust4_normalized.tsv'
    irep_file = '/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/irep_normalized.tsv'

    df_t4 = pd.read_csv(trust4_file, sep='\t')
    df_ir = pd.read_csv(irep_file, sep='\t')

    print(f"  TRUST4: {len(df_t4)} records")
    print(f"  iRepertoire: {len(df_ir)} records")

    return df_t4, df_ir

def report_data_counts(df_t4, df_ir):
    """Report data counts before processing"""
    print("\n" + "="*80)
    print("DATA COUNTS BEFORE PROCESSING")
    print("="*80)

    n_irep_total = len(df_ir)
    n_irep_stop = (df_ir['cdr3aa_norm'] == '*').sum()
    n_irep_productive = n_irep_total - n_irep_stop

    print(f"\n1. Raw iRepertoire records: {n_irep_total}")
    print(f"2. Records with CDR3(pep) == '*': {n_irep_stop} ({100*n_irep_stop/n_irep_total:.1f}%)")
    print(f"3. Productive iRepertoire records: {n_irep_productive}")

    n_trust4_total = len(df_t4)
    print(f"\n5. Raw TRUST4 IGH records: {n_trust4_total}")

    # Unique before aggregation
    n_irep_unique_before = df_ir[df_ir['cdr3aa_norm'] != '*'].groupby(
        ['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']
    ).ngroups
    n_trust4_unique_before = df_t4.groupby(
        ['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']
    ).ngroups

    print(f"4. Unique productive iRepertoire clonotypes (before agg): {n_irep_unique_before}")
    print(f"6. Unique TRUST4 clonotypes (before agg): {n_trust4_unique_before}")

    return {
        'n_irep_total': n_irep_total,
        'n_irep_stop': n_irep_stop,
        'n_irep_productive': n_irep_productive,
        'n_irep_unique_before': n_irep_unique_before,
        'n_trust4_total': n_trust4_total,
        'n_trust4_unique_before': n_trust4_unique_before
    }

def aggregate_clonotypes(df_t4, df_ir):
    """Aggregate duplicate clonotypes by CDR3aa+V+J+C"""
    print("\n" + "="*80)
    print("AGGREGATING CLONOTYPES")
    print("="*80)

    # Filter iRepertoire to productive only
    print(f"\nFiltering iRepertoire to exclude CDR3(pep) == '*'...")
    df_ir_prod = df_ir[df_ir['cdr3aa_norm'] != '*'].copy()
    print(f"  Removed: {len(df_ir) - len(df_ir_prod)} stop-codon records")
    print(f"  Remaining: {len(df_ir_prod)} productive records")

    # Aggregate TRUST4
    print(f"\nAggregating TRUST4 by CDR3aa+V+J+C...")
    t4_agg = df_t4.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']).agg({
        'CDR3nt': 'first',  # Keep first for reference
        'V': 'first',  # Keep original for traceability
        'D': 'first',
        'J': 'first',
        'C': 'first',
        'd_gene_norm': 'first',
        'count': 'sum',
        'frequency': 'sum',  # Will recalculate
        'cdr3nt_norm': lambda x: '|'.join(x.drop_duplicates()),  # Preserve variants
    }).reset_index()

    # Recalculate frequency
    total_counts = t4_agg['count'].sum()
    t4_agg['frequency'] = t4_agg['count'] / total_counts if total_counts > 0 else 0

    print(f"  Before aggregation: {len(df_t4)} records")
    print(f"  After aggregation: {len(t4_agg)} clonotypes")

    # Aggregate iRepertoire
    print(f"\nAggregating iRepertoire by CDR3aa+V+J+C...")
    ir_agg = df_ir_prod.groupby(['cdr3aa_norm', 'v_gene_norm', 'j_gene_norm', 'c_gene_norm']).agg({
        'CDR3nt': 'first',
        'V': 'first',
        'D': 'first',
        'J': 'first',
        'C': 'first',
        'd_gene_norm': 'first',
        'copy': 'sum',
        'cdr3nt_norm': lambda x: '|'.join(x.drop_duplicates()),  # Preserve variants
    }).reset_index()

    print(f"  Before aggregation: {len(df_ir_prod)} records")
    print(f"  After aggregation: {len(ir_agg)} clonotypes")

    return t4_agg, ir_agg, df_ir_prod

def match_clonotypes(t4_agg, ir_agg):
    """Match TRUST4 and iRepertoire clonotypes"""
    print("\n" + "="*80)
    print("MATCHING CLONOTYPES")
    print("="*80)

    # Create match index
    ir_index = {}
    for idx, row in ir_agg.iterrows():
        key = (row['cdr3aa_norm'], row['v_gene_norm'], row['j_gene_norm'], row['c_gene_norm'])
        ir_index[key] = idx

    print(f"\nMatching TRUST4 to iRepertoire...")
    matches = []

    for idx, t4_row in t4_agg.iterrows():
        key = (t4_row['cdr3aa_norm'], t4_row['v_gene_norm'], t4_row['j_gene_norm'], t4_row['c_gene_norm'])

        if key in ir_index:
            ir_idx = ir_index[key]
            ir_row = ir_agg.iloc[ir_idx]

            # Check agreement on D and CDR3nt
            d_match = t4_row['d_gene_norm'] == ir_row['d_gene_norm'] if pd.notna(t4_row['d_gene_norm']) and pd.notna(ir_row['d_gene_norm']) else None
            cdr3nt_match = t4_row['cdr3nt_norm'] == ir_row['cdr3nt_norm']

            matches.append({
                'cdr3aa': t4_row['cdr3aa_norm'],
                'trust4_v': t4_row['v_gene_norm'],
                'irep_v': ir_row['v_gene_norm'],
                'trust4_j': t4_row['j_gene_norm'],
                'irep_j': ir_row['j_gene_norm'],
                'trust4_c': t4_row['c_gene_norm'],
                'irep_c': ir_row['c_gene_norm'],
                'trust4_d': t4_row['d_gene_norm'],
                'irep_d': ir_row['d_gene_norm'],
                'd_match': d_match,
                'trust4_cdr3nt': t4_row['cdr3nt_norm'][:50],  # First 50 chars
                'irep_cdr3nt': ir_row['cdr3nt_norm'][:50],
                'cdr3nt_match': cdr3nt_match,
                'trust4_count': t4_row['count'],
                'trust4_frequency': t4_row['frequency'],
                'irep_copy': ir_row['copy'],
                'trust4_variants': t4_row['cdr3nt_norm'].count('|') + 1,
                'irep_variants': ir_row['cdr3nt_norm'].count('|') + 1,
            })

    df_matches = pd.DataFrame(matches)

    print(f"\n  Matches found: {len(df_matches)}")
    print(f"  TRUST4 clonotypes unmatched: {len(t4_agg) - len(df_matches)}")
    print(f"  iRepertoire clonotypes unmatched: {len(ir_agg) - len(df_matches)}")

    return df_matches

def calculate_metrics(t4_agg, ir_agg, df_matches):
    """Calculate primary and diagnostic metrics"""
    print("\n" + "="*80)
    print("CALCULATING METRICS")
    print("="*80)

    tp = len(df_matches)
    fp = len(t4_agg) - tp
    fn = len(ir_agg) - tp

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0

    print(f"\nPRIMARY METRICS:")
    print(f"  True Positives (TP): {tp}")
    print(f"  False Positives (FP): {fp}")
    print(f"  False Negatives (FN): {fn}")
    print(f"  Precision: {precision:.4f}")
    print(f"  Sensitivity (Recall): {sensitivity:.4f}")
    print(f"  F1-score: {f1:.4f}")

    # Diagnostic metrics
    print(f"\nDIAGNOSTIC METRICS (among matched clonotypes):")

    if len(df_matches) > 0:
        # V agreement
        v_agree = (df_matches['trust4_v'] == df_matches['irep_v']).sum()
        j_agree = (df_matches['trust4_j'] == df_matches['irep_j']).sum()
        c_agree = (df_matches['trust4_c'] == df_matches['irep_c']).sum()
        d_agree = df_matches['d_match'].sum()
        cdr3nt_agree = df_matches['cdr3nt_match'].sum()

        print(f"  V gene agreement: {v_agree}/{len(df_matches)} ({100*v_agree/len(df_matches):.1f}%)")
        print(f"  J gene agreement: {j_agree}/{len(df_matches)} ({100*j_agree/len(df_matches):.1f}%)")
        print(f"  C gene agreement: {c_agree}/{len(df_matches)} ({100*c_agree/len(df_matches):.1f}%)")
        print(f"  D gene agreement: {d_agree}/{len(df_matches)} ({100*d_agree/len(df_matches):.1f}%)")
        print(f"  CDR3nt exact match: {cdr3nt_agree}/{len(df_matches)} ({100*cdr3nt_agree/len(df_matches):.1f}%)")

        # Variant information
        t4_multi_variant = (df_matches['trust4_variants'] > 1).sum()
        ir_multi_variant = (df_matches['irep_variants'] > 1).sum()
        print(f"\n  TRUST4 clonotypes with multiple CDR3nt variants: {t4_multi_variant}")
        print(f"  iRepertoire clonotypes with multiple records: {ir_multi_variant}")

    metrics = {
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'precision': precision,
        'sensitivity': sensitivity,
        'f1': f1,
        'v_agreement': v_agree if len(df_matches) > 0 else 0,
        'j_agreement': j_agree if len(df_matches) > 0 else 0,
        'c_agreement': c_agree if len(df_matches) > 0 else 0,
        'd_agreement': d_agree if len(df_matches) > 0 else 0,
        'cdr3nt_agreement': cdr3nt_agree if len(df_matches) > 0 else 0,
    }

    return metrics

def calculate_abundance_correlation(df_matches):
    """Calculate Spearman correlation for abundance"""
    print(f"\nABUNDANCE CONCORDANCE (matched clonotypes):")

    if len(df_matches) < 2:
        print(f"  Insufficient matched pairs for correlation")
        return None

    # Remove NaN values
    df_corr = df_matches[['trust4_count', 'trust4_frequency', 'irep_copy']].dropna()

    if len(df_corr) < 2:
        print(f"  Insufficient data for correlation")
        return None

    corr_count, pval_count = spearmanr(df_corr['trust4_count'], df_corr['irep_copy'])
    corr_freq, pval_freq = spearmanr(df_corr['trust4_frequency'], df_corr['irep_copy'])

    print(f"  Spearman r (count vs copy): {corr_count:.4f} (p={pval_count:.4e})")
    print(f"  Spearman r (frequency vs copy): {corr_freq:.4f} (p={pval_freq:.4e})")
    print(f"  Note: This is CONCORDANCE analysis, not equivalence of molecule counts")

    return {
        'n_pairs': len(df_corr),
        'count_vs_copy_r': corr_count,
        'count_vs_copy_p': pval_count,
        'frequency_vs_copy_r': corr_freq,
        'frequency_vs_copy_p': pval_freq,
    }

def save_results(t4_agg, ir_agg, df_matches, metrics, corr, data_counts):
    """Save all output files"""
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)

    output_dir = '/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116'
    os.makedirs(output_dir, exist_ok=True)

    # Prepare output files
    t4_agg_out = t4_agg.rename(columns={
        'cdr3aa_norm': 'CDR3aa',
        'v_gene_norm': 'V_gene',
        'j_gene_norm': 'J_gene',
        'c_gene_norm': 'C_gene',
        'd_gene_norm': 'D_gene',
        'count': 'read_count',
        'cdr3nt_norm': 'CDR3nt_variants',
    })

    ir_agg_out = ir_agg.rename(columns={
        'cdr3aa_norm': 'CDR3aa',
        'v_gene_norm': 'V_gene',
        'j_gene_norm': 'J_gene',
        'c_gene_norm': 'C_gene',
        'd_gene_norm': 'D_gene',
        'copy': 'copy_count',
        'cdr3nt_norm': 'CDR3nt_variants',
    })

    # Save clonotype tables
    t4_agg_out.to_csv(f'{output_dir}/trust4_clonotypes.tsv', sep='\t', index=False)
    ir_agg_out.to_csv(f'{output_dir}/irep_clonotypes.tsv', sep='\t', index=False)

    # Save matches
    df_matches.to_csv(f'{output_dir}/matches.tsv', sep='\t', index=False)

    # Save metrics
    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(f'{output_dir}/metrics.tsv', sep='\t', index=False)

    # Save correlation
    if corr is not None:
        corr_df = pd.DataFrame([corr])
        corr_df.to_csv(f'{output_dir}/abundance_correlation.tsv', sep='\t', index=False)

    # Save data counts
    counts_df = pd.DataFrame([data_counts])
    counts_df.to_csv(f'{output_dir}/data_counts.tsv', sep='\t', index=False)

    print(f"\nResults saved to {output_dir}/:")
    print(f"  - trust4_clonotypes.tsv ({len(t4_agg_out)} clonotypes)")
    print(f"  - irep_clonotypes.tsv ({len(ir_agg_out)} clonotypes)")
    print(f"  - matches.tsv ({len(df_matches)} matches)")
    print(f"  - metrics.tsv")
    print(f"  - abundance_correlation.tsv")
    print(f"  - data_counts.tsv")

    return output_dir

def print_summary(data_counts, metrics):
    """Print final summary"""
    print("\n" + "="*80)
    print("FINAL SUMMARY - FZ-116 BENCHMARK")
    print("="*80)

    print(f"\nDATA COUNTS:")
    print(f"  1. Raw iRepertoire records: {data_counts['n_irep_total']}")
    print(f"  2. Excluded (CDR3 = '*'): {data_counts['n_irep_stop']}")
    print(f"  3. Productive iRepertoire: {data_counts['n_irep_productive']}")
    print(f"  4. iRepertoire clonotypes (aggregated): {data_counts['n_irep_unique_before']} → see irep_clonotypes.tsv")
    print(f"  5. Raw TRUST4 IGH records: {data_counts['n_trust4_total']}")
    print(f"  6. TRUST4 clonotypes (aggregated): {data_counts['n_trust4_unique_before']} → see trust4_clonotypes.tsv")
    print(f"  7. Primary benchmark matches: {metrics['tp']}")
    print(f"  8. Unmatched TRUST4 clonotypes: {metrics['fp']}")
    print(f"  9. Unmatched iRepertoire clonotypes: {metrics['fn']}")

    print(f"\nBENCHMARK METRICS:")
    print(f"  Precision: {metrics['precision']:.4f}")
    print(f"  Sensitivity (Recall): {metrics['sensitivity']:.4f}")
    print(f"  F1-score: {metrics['f1']:.4f}")

    print(f"\nVGene agreement (among {metrics['tp']} matches): {metrics['v_agreement']}/{metrics['tp']} ({100*metrics['v_agreement']/metrics['tp'] if metrics['tp'] > 0 else 0:.1f}%)")
    print(f"J-gene agreement: {metrics['j_agreement']}/{metrics['tp']} ({100*metrics['j_agreement']/metrics['tp'] if metrics['tp'] > 0 else 0:.1f}%)")
    print(f"C-gene agreement: {metrics['c_agreement']}/{metrics['tp']} ({100*metrics['c_agreement']/metrics['tp'] if metrics['tp'] > 0 else 0:.1f}%)")
    print(f"D-gene agreement: {metrics['d_agreement']}/{metrics['tp']} ({100*metrics['d_agreement']/metrics['tp'] if metrics['tp'] > 0 else 0:.1f}%)")

def main():
    print("="*80)
    print("TRUST4 vs iRepertoire Benchmarking - FZ-116 FINAL")
    print("="*80)

    # Load data
    df_t4, df_ir = load_normalized_data()

    # Report counts
    data_counts = report_data_counts(df_t4, df_ir)

    # Aggregate
    t4_agg, ir_agg, df_ir_prod = aggregate_clonotypes(df_t4, df_ir)

    # Match
    df_matches = match_clonotypes(t4_agg, ir_agg)

    # Calculate metrics
    metrics = calculate_metrics(t4_agg, ir_agg, df_matches)

    # Abundance correlation
    corr = calculate_abundance_correlation(df_matches)

    # Save
    output_dir = save_results(t4_agg, ir_agg, df_matches, metrics, corr, data_counts)

    # Print summary
    print_summary(data_counts, metrics)

    print("\n" + "="*80)
    print("BENCHMARK COMPLETE")
    print("="*80)
    print(f"\nAll results in: {output_dir}/")

    return 0

if __name__ == '__main__':
    sys.exit(main())
