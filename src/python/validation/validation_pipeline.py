#!/usr/bin/env python3
"""
DiCE-Duan Validation & Benchmarking Pipeline
=============================================
Publication-quality validation of disease gene predictions.

Usage:
    python validation_pipeline.py --all
    python validation_pipeline.py --analysis 1 2 3
    python validation_pipeline.py --all --skip-pubmed
"""

import argparse
import os
import sys
import warnings
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import hypergeom, mannwhitneyu, ttest_ind

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ─── Style Configuration ────────────────────────────────────────────────────
STYLE = {
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.figsize': (8, 6),
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
}
plt.rcParams.update(STYLE)
sns.set_theme(style='whitegrid', palette='deep')
PALETTE_MAIN = ['#2563EB', '#DC2626', '#16A34A', '#F59E0B', '#7C3AED']


# ═══════════════════════════════════════════════════════════════════════════
#  DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════

def find_project_root() -> Path:
    """Find the DICE project root."""
    # Walk up the directory tree until we find the 'data' folder
    current_dir = Path(__file__).resolve().parent
    while current_dir.name != '' and current_dir.parent != current_dir:
        if (current_dir / 'data' / 'results').exists() or (current_dir / 'data' / 'raw').exists():
            return current_dir
        current_dir = current_dir.parent
    return Path(__file__).resolve().parent


def load_dice_genes(root):
    """Load full ranked DiCE genes table (17,466 LCC universe)."""
    full_path = root / 'data' / 'results' / 'dice_genes.txt'
    top_path = root / 'data' / 'results' / 'dice_genes_top500.txt'
    path = full_path if full_path.exists() else top_path
    df = pd.read_csv(path, sep='\t')
    log.info(f"Loaded {len(df)} DiCE genes from {path.name}")
    return df


def load_centrality(root, condition):
    """Load centrality file (normal or tumor)."""
    path = root / 'data' / 'results' / f'centrality_{condition}.txt'
    df = pd.read_csv(path, sep='\t')
    log.info(f"Loaded {len(df)} genes from centrality_{condition}.txt")
    return df


def load_knockout(root):
    """Load knockout analysis results."""
    path = root / 'data' / 'results' / 'knockout_analysis.txt'
    df = pd.read_csv(path, sep='\t')
    log.info(f"Loaded {len(df)} knockout results")
    return df


def load_expression(root):
    """Load normal and tumor expression matrices."""
    normal_path = root / 'data' / 'raw' / 'expression' / 'expression_normal.txt'
    tumor_path = root / 'data' / 'raw' / 'expression' / 'expression_tumor.txt'

    normal = pd.read_csv(normal_path, sep='\t', index_col=0)
    tumor = pd.read_csv(tumor_path, sep='\t', index_col=0)
    log.info(f"Loaded expression: normal={normal.shape}, tumor={tumor.shape}")
    return normal, tumor


def load_known_genes(root):
    """Load known disease gene lists."""
    known_dir = root / 'data' / 'known_genes'
    result = {}
    for name, fname in [('COSMIC', 'cosmic_census.txt'), ('DisGeNET', 'disgenet_prostate.txt')]:
        path = known_dir / fname
        if path.exists():
            genes = set(l.strip() for l in path.read_text().splitlines() if l.strip() and not l.startswith('#'))
            result[name] = genes
            log.info(f"Loaded {len(genes)} known genes from {fname}")
        else:
            log.warning(f"Known gene file not found: {path}")
    return result


def load_network_edges(root):
    """Load merged network edges properly handling the 4-column TSV header."""
    path = root / 'data' / 'processed' / 'networks' / 'merged_network.txt'
    if not path.exists():
        return None
    edges = pd.read_csv(path, sep='\t', usecols=[0, 1], header=0, names=['gene1', 'gene2'])
    log.info(f"Loaded {len(edges)} network edges")
    return edges


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 1: PATHWAY ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════════

def analysis_pathway_enrichment(dice_df, root, out_fig, out_tbl):
    """Enrichment analysis of top DiCE genes split by Upregulated vs Downregulated."""
    log.info("=" * 60)
    log.info("Analysis 1: Pathway Enrichment (Split Up/Down)")
    log.info("=" * 60)

    try:
        import gseapy
    except ImportError:
        log.error("gseapy not installed. Skipping pathway enrichment.")
        return None

    # Load DEA log2FC data to split genes
    dea_table_path = root / 'data' / 'results' / 'dea_results_table.txt'
    if not dea_table_path.exists():
        log.warning(f"DEA table not found at {dea_table_path}. Cannot perform split enrichment.")
        return None
        
    dea_df = pd.read_csv(dea_table_path, sep='\t')
    log2fc_map = dict(zip(dea_df['gene'], dea_df['log2FC']))
    
    # 1. Split top 150 DiCE genes into Up and Down
    top_genes = dice_df['gene'].head(150).tolist()
    
    up_genes = [g for g in top_genes if log2fc_map.get(g, 0) > 0]
    down_genes = [g for g in top_genes if log2fc_map.get(g, 0) < 0]
    
    log.info(f"Split {len(top_genes)} top genes into {len(up_genes)} Upregulated and {len(down_genes)} Downregulated.")

    gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2021']
    
    # Process both lists
    def run_enrichr(genes, label):
        if len(genes) < 5:
            log.warning(f"Too few {label} genes for enrichr.")
            return None
        all_results = []
        for lib in gene_sets:
            try:
                enr = gseapy.enrichr(gene_list=genes, gene_sets=lib, organism='human', outdir=None, cutoff=0.05)
                res = enr.results
                res['Database'] = lib.split('_')[0]
                res['Regulation'] = label
                all_results.append(res)
                log.info(f"  {label} {lib}: {len(res[res['Adjusted P-value'] < 0.05])} significant terms")
            except Exception as e:
                log.warning(f"  Enrichr failed for {label} {lib}: {e}")
                
        if not all_results: return None
        return pd.concat(all_results, ignore_index=True)

    up_results = run_enrichr(up_genes, 'Upregulated')
    down_results = run_enrichr(down_genes, 'Downregulated')
    
    if up_results is None and down_results is None:
        return None
        
    combined = pd.concat([df for df in [up_results, down_results] if df is not None], ignore_index=True)
    sig = combined[combined['Adjusted P-value'] < 0.05].sort_values(['Regulation', 'Adjusted P-value'])
    sig.to_csv(out_tbl / '01_pathway_enrichment_split.csv', index=False)
    log.info(f"  Saved {len(sig)} significant terms to table")

    # ── Divergent Bar plot ──
    try:
        plot_df = []
        if up_results is not None:
            up_sig = up_results[up_results['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value').head(10)
            plot_df.append(up_sig)
        if down_results is not None:
            down_sig = down_results[down_results['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value').head(10)
            plot_df.append(down_sig)
            
        if not plot_df:
            return combined
            
        plot_df = pd.concat(plot_df, ignore_index=True)
        # Flip values for downregulated for divergent plot
        plot_df['neg_log10_padj'] = -np.log10(plot_df['Adjusted P-value'].clip(lower=1e-50))
        plot_df.loc[plot_df['Regulation'] == 'Downregulated', 'neg_log10_padj'] *= -1
        
        # Sort so Up are on top, Down on bottom
        plot_df = plot_df.sort_values('neg_log10_padj')
        plot_df['Term_short'] = plot_df['Term'].str[:50]
        
        fig, ax = plt.subplots(figsize=(12, 8))
        colors = ['#EF4444' if r == 'Upregulated' else '#3B82F6' for r in plot_df['Regulation']]
        
        bars = ax.barh(range(len(plot_df)), plot_df['neg_log10_padj'], color=colors)
        
        # Add labels dynamically inside/outside bars
        for i, (val, term) in enumerate(zip(plot_df['neg_log10_padj'], plot_df['Term_short'])):
            if val > 0:
                ax.text(0.1, i, term, va='center', ha='left', fontsize=9, color='white', fontweight='bold')
            else:
                ax.text(-0.1, i, term, va='center', ha='right', fontsize=9, color='white', fontweight='bold')
                
        ax.set_yticks([]) # Hide traditional y-ticks
        
        # Format axes
        ax.axvline(0, color='black', linewidth=1)
        ax.set_xlabel('Significance ($-\\log_{10}$(Adjusted p-value))')
        # Labels for axes
        ticks = ax.get_xticks()
        ax.set_xticklabels([str(abs(int(t))) for t in ticks])
        
        ax.set_title('Pathway Enrichment (Upregulated DiCE vs Downregulated DiCE)')
        
        from matplotlib.patches import Patch
        legend_elems = [Patch(color='#EF4444', label='Upregulated (Log2FC > 0)'), 
                        Patch(color='#3B82F6', label='Downregulated (Log2FC < 0)')]
        ax.legend(handles=legend_elems, loc='upper left')
        
        plt.tight_layout()
        fig.savefig(out_fig / '01_pathway_enrichment_divergent_barplot.png')
        plt.close(fig)
        log.info("  Saved divergent barplot")
    except Exception as e:
        log.warning(f"  Failed generating divergent plot: {e}")

    return combined


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 2: HYPERGEOMETRIC SIGNIFICANCE TEST
# ═══════════════════════════════════════════════════════════════════════════

def analysis_hypergeometric(dice_df, known_genes, universe_size, out_fig, out_tbl):
    """Hypergeometric test for overlap with known disease genes."""
    log.info("=" * 60)
    log.info("Analysis 2: Hypergeometric Significance Test")
    log.info("=" * 60)

    predicted = set(dice_df['gene'].head(100).tolist())
    results = []

    for db_name, db_genes in known_genes.items():
        overlap = predicted & db_genes
        k = len(overlap)  # successes in sample
        M = universe_size  # population size
        n = len(db_genes & set(dice_df['gene']))  # total successes in population
        N = len(predicted)  # sample size

        # P(X >= k) — survival function
        pval = hypergeom.sf(k - 1, M, n, N)
        expected = N * n / M

        fold_enrichment = (k / expected) if expected > 0 else float('inf')

        results.append({
            'Database': db_name,
            'Known_genes_in_universe': n,
            'Predicted_top100': N,
            'Overlap': k,
            'Expected_overlap': round(expected, 2),
            'Fold_enrichment': round(fold_enrichment, 2),
            'P_value': pval,
            'Significant': pval < 0.05,
            'Overlapping_genes': ', '.join(sorted(overlap))
        })
        log.info(f"  {db_name}: overlap={k}, expected={expected:.2f}, fold={fold_enrichment:.2f}, p={pval:.2e}")

    res_df = pd.DataFrame(results)
    res_df.to_csv(out_tbl / '02_hypergeometric_results.csv', index=False)

    # ── Plot ──
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for i, (_, row) in enumerate(res_df.iterrows()):
        ax = axes[i]
        vals = [row['Expected_overlap'], row['Overlap']]
        colors_bar = ['#94A3B8', PALETTE_MAIN[0] if row['Significant'] else '#EF4444']
        bars = ax.bar(['Expected', 'Observed'], vals, color=colors_bar, width=0.5, edgecolor='black', linewidth=0.8)
        ax.set_title(f"{row['Database']}\n(p = {row['P_value']:.2e})", fontweight='bold')
        ax.set_ylabel('Gene Overlap Count')

        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                    f'{val:.1f}' if isinstance(val, float) else str(val),
                    ha='center', va='bottom', fontweight='bold')

        ax.set_ylim(0, max(vals) * 1.3 + 1)

    fig.suptitle('Overlap with Known Disease Genes (Top 100 DiCE Genes)', fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(out_fig / '02_hypergeometric_overlap.png')
    plt.close(fig)
    log.info("  Saved overlap plot")
    return res_df


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 3: PRECISION@K
# ═══════════════════════════════════════════════════════════════════════════

def analysis_precision_at_k(dice_df, known_genes, out_fig, out_tbl):
    """Precision@K for K = 10, 20, 50, 100."""
    log.info("=" * 60)
    log.info("Analysis 3: Precision@K Evaluation")
    log.info("=" * 60)

    K_values = [10, 20, 50, 100, 200, 500]
    all_known = set()
    for genes in known_genes.values():
        all_known |= genes

    records = []
    per_db = {db: [] for db in known_genes}
    per_db['Combined'] = []

    for k in K_values:
        top_k = set(dice_df['gene'].head(k))
        # Combined
        prec = len(top_k & all_known) / k
        per_db['Combined'].append(prec)
        row = {'K': k, 'Combined_precision': round(prec, 4)}

        for db_name, db_genes in known_genes.items():
            p = len(top_k & db_genes) / k
            per_db[db_name].append(p)
            row[f'{db_name}_precision'] = round(p, 4)

        records.append(row)
        log.info(f"  K={k}: combined precision = {prec:.4f}")

    res_df = pd.DataFrame(records)
    res_df.to_csv(out_tbl / '03_precision_at_k.csv', index=False)

    # ── Plot ──
    fig, ax = plt.subplots(figsize=(8, 5.5))
    markers = ['o', 's', 'D']
    colors_pk = [PALETTE_MAIN[0], PALETTE_MAIN[1], PALETTE_MAIN[2]]

    for i, (name, vals) in enumerate(per_db.items()):
        ax.plot(K_values, vals, marker=markers[i % len(markers)], color=colors_pk[i % len(colors_pk)],
                label=name, linewidth=2, markersize=7)

    ax.set_xlabel('K (Number of Top Genes)')
    ax.set_ylabel('Precision@K')
    ax.set_title('Precision@K — DiCE-Duan Gene Predictions')
    ax.legend(frameon=True)
    ax.set_xticks(K_values)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_fig / '03_precision_at_k.png')
    plt.close(fig)
    log.info("  Saved Precision@K plot")
    return res_df


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 4: NETWORK TOPOLOGY VALIDATION (Scale vs Magnitude)
# ═══════════════════════════════════════════════════════════════════════════

def analysis_network_topology(dice_df, out_fig, out_tbl):
    """Plot Characteristic Scale vs Magnitude (Aim 1)."""
    log.info("=" * 60)
    log.info("Analysis 4: Scale vs Magnitude Topology")
    log.info("=" * 60)

    # We use dice_df directly because it now contains all structural information
    if 'characteristic_scale' not in dice_df.columns:
        log.error("Missing Aim 1 columns in dice_genes.txt")
        return None

    top_genes = set(dice_df['gene'].head(100))
    dice_df['is_candidate'] = dice_df['gene'].isin(top_genes)

    # ── Scatter plot ──
    fig, ax = plt.subplots(figsize=(8, 7))
    bg = dice_df[~dice_df['is_candidate']]
    fg = dice_df[dice_df['is_candidate']]

    # X-axis is Scale (Log scale), Y-axis is Magnitude
    ax.scatter(bg['characteristic_scale'], bg['ensemble_score'], alpha=0.15, s=8, c='#94A3B8', label='Background')
    ax.scatter(fg['characteristic_scale'], fg['ensemble_score'], alpha=0.8, s=25, c='#DC2626',
               edgecolors='black', linewidths=0.3, label='DiCE Top 100', zorder=5)

    ax.set_xscale('log')
    ax.set_xlabel('Characteristic Scale (Diffusion Time $t^*$)')
    ax.set_ylabel('Magnitude of Structural Shift ($\\Delta$)')
    ax.set_title('Structural Repositioning: Scale vs Magnitude')
    ax.legend(frameon=True, loc='upper right')
    
    plt.tight_layout()
    fig.savefig(out_fig / '04_scale_magnitude_scatter.png')
    plt.close(fig)
    log.info("  Saved Scale vs Magnitude scatter plot")

    return dice_df

# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 5: DIFFERENTIAL NETWORK REWIRING (FDR Volcano Plot)
# ═══════════════════════════════════════════════════════════════════════════

def analysis_rewiring(dice_df, out_fig, out_tbl):
    """Volcano plot of Structural Rewiring (Magnitude vs FDR)."""
    log.info("=" * 60)
    log.info("Analysis 5: Differential Rewiring Volcano Plot")
    log.info("=" * 60)

    if 'q_value_fdr' not in dice_df.columns:
        log.error("Missing FDR column. Skipping Volcano plot.")
        return None

    # Calculate -log10(q-value) for volcano Y-axis
    dice_df['neg_log10_q'] = -np.log10(dice_df['q_value_fdr'].clip(lower=1e-300))
    
    # We use signed magnitude for X-axis (Magnitude * Direction)
    dice_df['signed_magnitude'] = dice_df['ensemble_score'] * dice_df['direction']

    top_genes = set(dice_df['gene'].head(100))
    dice_df['is_candidate'] = dice_df['gene'].isin(top_genes)

    # ── Volcano plot ──
    fig, ax = plt.subplots(figsize=(9, 7))
    bg = dice_df[~dice_df['is_candidate']]
    fg = dice_df[dice_df['is_candidate']]

    ax.scatter(bg['signed_magnitude'], bg['neg_log10_q'], alpha=0.2, s=8, c='#94A3B8', label='Background')
    ax.scatter(fg['signed_magnitude'], fg['neg_log10_q'], alpha=0.8, s=30, c='#DC2626',
               edgecolors='black', linewidths=0.3, label='DiCE Top 100', zorder=5)

    # Label top rewired genes
    top_rewired = fg.nlargest(10, 'neg_log10_q')
    for _, row in top_rewired.iterrows():
        ax.annotate(row['gene'], (row['signed_magnitude'], row['neg_log10_q']),
                    fontsize=8, fontweight='bold', color='#DC2626',
                    xytext=(5, 5), textcoords='offset points')

    # Draw FDR < 0.05 line
    ax.axhline(-np.log10(0.05), color='grey', linestyle='--', alpha=0.7, label='FDR q = 0.05')
    ax.axvline(0, color='grey', linestyle='-', alpha=0.3)
    
    ax.set_xlabel('Signed Structural Shift (Localizing ← vs → Bridging)')
    ax.set_ylabel('$-\\log_{10}$(FDR q-value)')
    ax.set_title('Network Rewiring Volcano Plot (Aim 1)')
    ax.legend(frameon=True, loc='upper right')
    
    plt.tight_layout()
    fig.savefig(out_fig / '05_rewiring_volcano.png')
    plt.close(fig)
    log.info("  Saved FDR volcano plot")

    return dice_df


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 6: LITERATURE VALIDATION
# ═══════════════════════════════════════════════════════════════════════════

def analysis_literature(dice_df, out_fig, out_tbl, skip_pubmed=False):
    """Query PubMed for disease-related publications."""
    log.info("=" * 60)
    log.info("Analysis 6: Literature Validation")
    log.info("=" * 60)

    top_genes = dice_df['gene'].head(20).tolist()
    records = []

    if skip_pubmed:
        log.info("  PubMed queries skipped (--skip-pubmed flag)")
        # Generate mock data for demonstration
        np.random.seed(42)
        for gene in top_genes:
            count = int(np.random.exponential(50)) + 1
            records.append({'Gene': gene, 'Rank': top_genes.index(gene) + 1,
                            'PubMed_count': count, 'Source': 'mock'})
    else:
        try:
            from Bio import Entrez
            import time
            Entrez.email = 'dice.validation@example.com'

            for i, gene in enumerate(top_genes):
                query = f'"{gene}"[Title/Abstract] AND ("prostate cancer"[Title/Abstract] OR "prostate neoplasm"[MeSH])'
                try:
                    handle = Entrez.esearch(db='pubmed', term=query, retmax=0)
                    result = Entrez.read(handle)
                    handle.close()
                    count = int(result['Count'])
                    records.append({'Gene': gene, 'Rank': i + 1,
                                    'PubMed_count': count, 'Source': 'pubmed'})
                    log.info(f"  {gene}: {count} publications")
                    time.sleep(0.35)  # NCBI rate limit
                except Exception as e:
                    log.warning(f"  PubMed query failed for {gene}: {e}")
                    records.append({'Gene': gene, 'Rank': i + 1,
                                    'PubMed_count': 0, 'Source': 'error'})
        except ImportError:
            log.warning("  Biopython not installed. Using mock data.")
            np.random.seed(42)
            for i, gene in enumerate(top_genes):
                count = int(np.random.exponential(50)) + 1
                records.append({'Gene': gene, 'Rank': i + 1,
                                'PubMed_count': count, 'Source': 'mock'})

    res_df = pd.DataFrame(records)
    res_df.to_csv(out_tbl / '06_literature_counts.csv', index=False)

    # ── Bar chart ──
    fig, ax = plt.subplots(figsize=(10, 6))
    colors_lit = [PALETTE_MAIN[0] if c > 10 else '#94A3B8' for c in res_df['PubMed_count']]
    bars = ax.barh(range(len(res_df)), res_df['PubMed_count'], color=colors_lit, edgecolor='black', linewidth=0.3)
    ax.set_yticks(range(len(res_df)))
    ax.set_yticklabels(res_df['Gene'], fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel('Number of Publications (Prostate Cancer)')
    ax.set_title('Literature Validation — Top 20 DiCE Genes')

    source = res_df['Source'].iloc[0] if len(res_df) > 0 else 'unknown'
    if source == 'mock':
        ax.text(0.98, 0.02, 'Mock data (PubMed not queried)', transform=ax.transAxes,
                fontsize=8, ha='right', va='bottom', style='italic', color='grey')

    plt.tight_layout()
    fig.savefig(out_fig / '06_literature_validation.png')
    plt.close(fig)
    log.info("  Saved literature plot")
    return res_df


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 7: VIRTUAL KNOCKOUT SIMULATION (Moved to Aim 2)
# ═══════════════════════════════════════════════════════════════════════════

def analysis_knockout(root):
    """Placeholder for Aim 2 Edge Ablation."""
    log.info("=" * 60)
    log.info("Analysis 7: Virtual Knockout / Edge Ablation")
    log.info("=" * 60)
    log.info("  NOTE: Node-level virtual knockout is deprecated in this version.")
    log.info("  Exact edge-level ablation will be performed via the Siamese GNN in Aim 2.")
    return None


# ═══════════════════════════════════════════════════════════════════════════
#  BENCHMARKING: DiCE vs DE vs Degree Centrality
# ═══════════════════════════════════════════════════════════════════════════

def analysis_benchmarking(dice_df, network_edges, expr_normal, expr_tumor,
                          known_genes, out_fig, out_tbl):
    """Compare DiCE-Flow (Aim 1) against baseline methods."""
    log.info("=" * 60)
    log.info("Benchmarking: Aim 1 vs DE Ranking vs Degree Centrality")
    log.info("=" * 60)

    all_known = set()
    for genes in known_genes.values():
        all_known |= genes

    methods = {}

    # 1. DiCE Aim 1 ranking
    methods['DiCE Aim 1 (HKS)'] = dice_df['gene'].tolist()

    # 2. Differential Expression ranking
    common_genes = list(set(expr_normal.index) & set(expr_tumor.index))
    if common_genes:
        log.info("  Computing differential expression ranking...")
        de_results = []
        for gene in common_genes:
            n_vals = expr_normal.loc[gene].values.astype(float)
            t_vals = expr_tumor.loc[gene].values.astype(float)
            n_vals = n_vals[n_vals > 0]
            t_vals = t_vals[t_vals > 0]
            if len(n_vals) >= 3 and len(t_vals) >= 3:
                stat, pval = ttest_ind(t_vals, n_vals, equal_var=False)
                log2fc = np.mean(t_vals) - np.mean(n_vals) 
                de_results.append({'gene': gene, 'abs_log2FC': abs(log2fc)})

        de_df = pd.DataFrame(de_results).sort_values('abs_log2FC', ascending=False)
        methods['DE Ranking'] = de_df['gene'].tolist()
        log.info(f"  DE ranking: {len(de_df)} genes")

    # 3. Raw Degree Centrality baseline (calculated directly from edges)
    if network_edges is not None:
        log.info("  Computing degree centrality baseline...")
        # Count occurrences in source and target columns
        degree_counts = pd.concat([network_edges['gene1'], network_edges['gene2']]).value_counts()
        methods['Degree Centrality'] = degree_counts.index.tolist()
        log.info(f"  Degree ranking: {len(methods['Degree Centrality'])} genes")

    # Compute Precision@K for all methods
    K_values = [10, 20, 50, 100, 200, 500]
    benchmark_records = []

    for method_name, gene_list in methods.items():
        for k in K_values:
            top_k = set(gene_list[:k])
            prec = len(top_k & all_known) / k
            benchmark_records.append({'Method': method_name, 'K': k, 'Precision': round(prec, 4)})

    bench_df = pd.DataFrame(benchmark_records)
    bench_df.to_csv(out_tbl / '08_benchmark_comparison.csv', index=False)

    # ── Precision@K Comparison Plot ──
    fig, ax = plt.subplots(figsize=(9, 6))
    method_colors = {'DiCE Aim 1 (HKS)': PALETTE_MAIN[0], 'DE Ranking': PALETTE_MAIN[1],
                     'Degree Centrality': PALETTE_MAIN[3]}
    method_markers = {'DiCE Aim 1 (HKS)': 'o', 'DE Ranking': 's', 'Degree Centrality': 'D'}

    for method_name in methods:
        subset = bench_df[bench_df['Method'] == method_name]
        ax.plot(subset['K'], subset['Precision'],
                marker=method_markers.get(method_name, '^'),
                color=method_colors.get(method_name, '#666'),
                label=method_name, linewidth=2.5, markersize=8)

    ax.set_xlabel('K (Number of Top Genes)')
    ax.set_ylabel('Precision@K (Known Disease Gene Overlap)')
    ax.set_title('Benchmarking: DiCE Aim 1 vs Baselines')
    ax.legend(frameon=True, loc='best')
    ax.set_xticks(K_values)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_fig / '08_benchmark_comparison.png')
    plt.close(fig)
    log.info("  Saved benchmark plot")

    return bench_df


# ═══════════════════════════════════════════════════════════════════════════
#  SUMMARY REPORT
# ═══════════════════════════════════════════════════════════════════════════

def write_summary(out_dir, results):
    """Write a text summary of all validation results."""
    summary_path = out_dir / 'validation_summary.txt'
    lines = [
        "=" * 70,
        "DiCE-Duan Validation Pipeline — Summary Report",
        "=" * 70,
        ""
    ]

    if 'hypergeometric' in results and results['hypergeometric'] is not None:
        lines.append("HYPERGEOMETRIC SIGNIFICANCE TEST")
        lines.append("-" * 40)
        for _, row in results['hypergeometric'].iterrows():
            lines.append(f"  {row['Database']}: overlap={row['Overlap']}, "
                         f"expected={row['Expected_overlap']}, "
                         f"fold={row['Fold_enrichment']}, p={row['P_value']:.2e}")
        lines.append("")

    if 'precision' in results and results['precision'] is not None:
        lines.append("PRECISION@K")
        lines.append("-" * 40)
        for _, row in results['precision'].iterrows():
            lines.append(f"  K={int(row['K'])}: combined={row['Combined_precision']:.4f}")
        lines.append("")

    if 'benchmark' in results and results['benchmark'] is not None:
        lines.append("BENCHMARKING")
        lines.append("-" * 40)
        for _, row in results['benchmark'].iterrows():
            lines.append(f"  {row['Method']} @ K={int(row['K'])}: precision={row['Precision']:.4f}")
        lines.append("")

    lines.extend([
        "=" * 70,
        "Output files generated in:",
        f"  Figures: {out_dir / 'figures'}",
        f"  Tables:  {out_dir / 'tables'}",
        "=" * 70
    ])

    summary_path.write_text('\n'.join(lines))
    log.info(f"Summary written to {summary_path}")


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='DiCE-Flow Validation Pipeline')
    parser.add_argument('--all', action='store_true', help='Run all analyses')
    parser.add_argument('--analysis', nargs='+', type=int, choices=[1,2,3,4,5,6,7,8,9,10])
    parser.add_argument('--output-dir', type=str, default=None)
    parser.add_argument('--skip-pubmed', action='store_true')
    args = parser.parse_args()

    if not args.all and not args.analysis:
        args.all = True

    analyses = set(range(1, 11)) if args.all else set(args.analysis)

    root = find_project_root()
    out_dir = Path(args.output_dir) if args.output_dir else root / 'data' / 'results' / 'validation'
    out_fig = out_dir / 'figures'
    out_tbl = out_dir / 'tables'
    out_fig.mkdir(parents=True, exist_ok=True)
    out_tbl.mkdir(parents=True, exist_ok=True)

    log.info("DiCE-Flow Validation Pipeline")
    dice_df = load_dice_genes(root)
    universe_size = 17466 # Approx size of Largest Connected Component

    results = {}

    if 1 in analyses:
        results['enrichment'] = analysis_pathway_enrichment(dice_df, root, out_fig, out_tbl)

    if 2 in analyses:
        known = load_known_genes(root)
        if known:
            results['hypergeometric'] = analysis_hypergeometric(dice_df, known, universe_size, out_fig, out_tbl)

    if 3 in analyses:
        known = load_known_genes(root)
        if known:
            results['precision'] = analysis_precision_at_k(dice_df, known, out_fig, out_tbl)

    if 4 in analyses:
        results['topology'] = analysis_network_topology(dice_df, out_fig, out_tbl)

    if 5 in analyses:
        results['rewiring'] = analysis_rewiring(dice_df, out_fig, out_tbl)

    if 6 in analyses:
        results['literature'] = analysis_literature(dice_df, out_fig, out_tbl, skip_pubmed=args.skip_pubmed)

    if 7 in analyses:
        results['knockout'] = analysis_knockout(root)

    if 8 in analyses:
        known = load_known_genes(root)
        expr_n, expr_t = load_expression(root)
        edges = load_network_edges(root)
        if known:
            results['benchmark'] = analysis_benchmarking(dice_df, edges, expr_n, expr_t, known, out_fig, out_tbl)

    if 9 in analyses:
        try:
            from fitness_analysis import run_fitness_gene_analysis
            results['fitness'] = run_fitness_gene_analysis(dice_df, root, out_fig, out_tbl, universe_size=universe_size)
        except ImportError:
            log.warning("Skipping Analysis 9 (fitness_analysis.py not found).")

    if 10 in analyses:
        # Same as before, looks for secondary dataset
        pass

    write_summary(out_dir, results)

if __name__ == '__main__':
    main()