#!/usr/bin/env python3
"""
Clinical Survival Validation Module for DiCE-Duan
=================================================
This module analyzes the top-ranked DiCE disease genes to determine
their clinical prognostic relevance. It performs:
1. Kaplan-Meier survival analysis (High vs. Low expression)
2. Log-rank statistical testing
3. Cox Proportional Hazards regression

It expects:
- Expression Matrix (Gene x Patient)
- Clinical Metadata (patient_id, days_to_death, event)
- DiCE Ranked Genes list

Outputs:
- Kaplan-Meier survival curve plots with confidence intervals.
- A summary CSV containing Hazard Ratios and P-values.
"""

import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
except ImportError:
    print("CRITICAL ERROR: The 'lifelines' library is required for survival analysis.")
    print("Please install it using: pip install lifelines")
    sys.exit(1)

logging.basicConfig(level=logging.INFO, format='INFO: %(message)s')
logger = logging.getLogger(__name__)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Clinical Survival Validation Module")
    parser.add_argument('--expression', required=True, help="Expression matrix file (CSV or TSV, Gene x Patient)")
    parser.add_argument('--clinical', required=True, help="Clinical data file (CSV or TSV with patient_id, days_to_death, event)")
    parser.add_argument('--dice-genes', required=True, help="Ranked DiCE genes file (e.g., dice_genes_top500.txt)")
    parser.add_argument('--top-n', type=int, default=20, help="Number of top genes to analyze (default: 20)")
    parser.add_argument('--out-dir', default='data/results/survival', help="Output directory for plots and results")
    return parser.parse_args()


def load_dataset(expr_path, clin_path, sep=None):
    """
    Robustly loads expression and clinical data, inferring separators.
    """
    logger.info(f"Loading expression data from {expr_path}")
    # Auto-detect separator if not provided
    if str(expr_path).endswith('.tsv') or str(expr_path).endswith('.txt'):
        expr_sep = '\t'
    else:
        expr_sep = ','

    expr_df = pd.read_csv(expr_path, sep=expr_sep, index_col=0)
    
    logger.info(f"Loading clinical data from {clin_path}")
    if str(clin_path).endswith('.tsv') or str(clin_path).endswith('.txt'):
        clin_sep = '\t'
    else:
        clin_sep = ','
    
    clin_df = pd.read_csv(clin_path, sep=clin_sep)
    
    # Standardize column names for clinical Dataframe
    clin_df.columns = clin_df.columns.str.lower().str.strip()

    required_cols = ['patient_id', 'days_to_death', 'event']
    for col in required_cols:
        if col not in clin_df.columns:
            logger.error(f"Clinical data is missing required column: '{col}'")
            logger.error(f"Available columns: {clin_df.columns.tolist()}")
            sys.exit(1)
            
    clin_df.set_index('patient_id', inplace=True)
    
    return expr_df, clin_df


def validate_and_merge(expr_df, clin_df, target_genes):
    """
    Validates patient IDs, filters for target genes, and merges expression with clinical data.
    """
    # Transpose expression matrix so patients are rows and genes are columns
    expr_t = expr_df.T
    
    # Strip any whitespace from patient IDs in both matrices
    expr_t.index = expr_t.index.astype(str).str.strip()
    clin_df.index = clin_df.index.astype(str).str.strip()

    # Highly common in TCGA data: Expression patient IDs (e.g. TCGA-A8-A08G-01A-...)
    # are longer than clinical IDs (e.g. TCGA-A8-A08G). Let's attempt to match cleanly.
    common_patients = set(expr_t.index).intersection(set(clin_df.index))
    
    if len(common_patients) == 0:
        logger.warning("No exact patient ID overlaps found between expression and clinical data.")
        logger.info("Attempting TCGA-style barcode matching (first 12 characters)...")
        expr_t.index = expr_t.index.str[:12]
        clin_df.index = clin_df.index.str[:12]
        
        # Remove duplicates caused by barcode truncation (e.g. multiple samples per patient)
        expr_t = expr_t[~expr_t.index.duplicated(keep='first')]
        clin_df = clin_df[~clin_df.index.duplicated(keep='first')]
        
        common_patients = set(expr_t.index).intersection(set(clin_df.index))
        
    if len(common_patients) == 0:
        logger.error("Failed to match patient IDs between expression and clinical files.")
        logger.error(f"Expression IDs sample: {list(expr_t.index[:5])}")
        logger.error(f"Clinical IDs sample: {list(clin_df.index[:5])}")
        sys.exit(1)
        
    logger.info(f"Successfully matched {len(common_patients)} patients.")

    # Filter for target genes (only keeping those present in the expression matrix)
    valid_genes = [g for g in target_genes if g in expr_t.columns]
    missing = len(target_genes) - len(valid_genes)
    if missing > 0:
        logger.warning(f"{missing} target genes not found in expression matrix.")
        
    expr_t = expr_t[valid_genes]
    
    # Merge using an inner join on patient ID
    merged_df = expr_t.join(clin_df[['days_to_death', 'event']], how='inner')
    
    # Drop rows with missing survival data
    merged_df.dropna(subset=['days_to_death', 'event'], inplace=True)
    
    # Ensure numeric types
    merged_df['days_to_death'] = pd.to_numeric(merged_df['days_to_death'], errors='coerce')
    merged_df['event'] = pd.to_numeric(merged_df['event'], errors='coerce')
    merged_df.dropna(subset=['days_to_death', 'event'], inplace=True)
    
    logger.info(f"Final usable cohort size: {len(merged_df)} patients.")
    
    return merged_df, valid_genes


def plot_kaplan_meier(kmf_high, kmf_low, gene_name, p_value, hr, out_dir):
    """
    Generates and saves a Kaplan-Meier survival plot with confidence intervals.
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot curves with confidence bands
    kmf_high.plot_survival_function(ax=ax, color='#DC2626', label='High Expression', ci_alpha=0.2)
    kmf_low.plot_survival_function(ax=ax, color='#2563EB', label='Low Expression', ci_alpha=0.2)

    # Styling
    ax.set_title(f'Survival Analysis: {gene_name}\nLog-Rank p-value: {p_value:.2e} | HR: {hr:.2f}', fontsize=14, loc='left')
    ax.set_xlabel('Time (Days)', fontsize=12)
    ax.set_ylabel('Survival Probability', fontsize=12)
    ax.grid(alpha=0.3)
    
    # Improve spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    out_path = out_dir / f"KM_plot_{gene_name}.png"
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def perform_survival_analysis(merged_df, genes, out_dir):
    """
    Executes Kaplan-Meier estimation, Log-Rank testing, and Cox regression.
    """
    results = []
    km_dir = out_dir / 'kaplan_meier_plots'
    km_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Performing survival analysis for {len(genes)} genes...")

    for gene in genes:
        try:
            # 1. Stratify patients (High vs Low expression based on median)
            median_expr = merged_df[gene].median()
            
            # Create boolean masks
            high_mask = merged_df[gene] >= median_expr
            low_mask = ~high_mask
            
            time = merged_df['days_to_death']
            event = merged_df['event']
            
            # --- Kaplan-Meier & Log-Rank Test ---
            # Compute Log-Rank Test to get the p-value
            results_logrank = logrank_test(
                time[high_mask], time[low_mask],
                event_observed_A=event[high_mask], event_observed_B=event[low_mask]
            )
            p_value = results_logrank.p_value

            # Fit Kaplan-Meier for plotting
            kmf_high = KaplanMeierFitter()
            kmf_high.fit(time[high_mask], event[high_mask], label='High Expression')
            
            kmf_low = KaplanMeierFitter()
            kmf_low.fit(time[low_mask], event[low_mask], label='Low Expression')

            # --- Cox Proportional Hazards Regression ---
            # To get Hazard Ratio (HR) and its confidence intervals
            cox_df = merged_df[[gene, 'days_to_death', 'event']].copy()
            cph = CoxPHFitter()
            cph.fit(cox_df, duration_col='days_to_death', event_col='event')
            
            # Extract metrics
            hr = cph.hazard_ratios_.iloc[0]
            hr_ci_lower = np.exp(cph.confidence_intervals_.iloc[0, 0])
            hr_ci_upper = np.exp(cph.confidence_intervals_.iloc[0, 1])
            cox_p_value = cph.summary['p'].iloc[0]

            # Store result
            results.append({
                'gene': gene,
                'hazard_ratio': round(hr, 4),
                'hr_ci_lower': round(hr_ci_lower, 4),
                'hr_ci_upper': round(hr_ci_upper, 4),
                'logrank_p_value': p_value,
                'cox_p_value': cox_p_value
            })

            # 2. Plotting
            plot_kaplan_meier(kmf_high, kmf_low, gene, p_value, hr, km_dir)
            
        except Exception as e:
            logger.error(f"Failed to analyze gene {gene}: {e}")

    # Create summary dataframe
    results_df = pd.DataFrame(results)
    
    if not results_df.empty:
        # Sort by statistical significance
        results_df.sort_values(by='logrank_p_value', inplace=True)
        out_csv = out_dir / 'survival_results.csv'
        results_df.to_csv(out_csv, index=False)
        logger.info(f"Saved numerical survival results to {out_csv}")
    else:
        logger.warning("No significant genes were successfully analyzed.")
        
    return results_df


def main():
    args = parse_arguments()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load ranked genes
    logger.info(f"Loading top {args.top_n} genes from {args.dice_genes}")
    try:
        dice_df = pd.read_csv(args.dice_genes, sep='\t')
        target_genes = dice_df['gene'].head(args.top_n).tolist()
    except Exception as e:
        logger.error(f"Failed to read DiCE genes file: {e}")
        sys.exit(1)

    # 2. Load Expression and Clinical Datasets
    expr_df, clin_df = load_dataset(args.expression, args.clinical)

    # 3. Validate, Filter, and Merge
    merged_df, valid_genes = validate_and_merge(expr_df, clin_df, target_genes)

    if len(valid_genes) == 0:
        logger.error("No valid genes to analyze. Exiting.")
        sys.exit(1)

    # 4. Perform Survival Analysis & Plotting
    results_df = perform_survival_analysis(merged_df, valid_genes, out_dir)
    
    logger.info("\n=== Analysis Summary ===")
    if not results_df.empty:
        sig_count = (results_df['logrank_p_value'] < 0.05).sum()
        logger.info(f"Total genes analyzed:     {len(results_df)}")
        logger.info(f"Significant genes (<0.05): {sig_count}")
        print("\nTop 5 Prognostic Genes:")
        print(results_df.head(5).to_string(index=False))
    logger.info("========================")
    

if __name__ == "__main__":
    main()
