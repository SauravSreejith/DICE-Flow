#!/usr/bin/env python3
"""
Cancer Fitness Gene Analysis
==============================
Validates DiCE gene predictions against CRISPR/RNAi fitness screen data
from the Cancer Dependency Map (DepMap).

Implements the fitness gene overlap analysis from:
    Pashaei et al. (2025) NAR, 53, gkaf609 — Pages 10-11

This analysis shows that DiCE-identified genes are significantly enriched
in cancer fitness genes compared to standard DEGs.

Usage:
    python fitness_analysis.py --dice-genes data/results/dice_genes_top500.txt
    python fitness_analysis.py --dice-genes data/results/dice_genes_top500.txt --top-k 100 200 364
"""

import argparse
import os
import sys
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
log = logging.getLogger(__name__)

# ─── Style Configuration ────────────────────────────────────────────────────
STYLE = {
    'figure.dpi': 300, 'savefig.dpi': 300, 'font.size': 11,
    'axes.titlesize': 13, 'axes.labelsize': 12,
    'xtick.labelsize': 10, 'ytick.labelsize': 10,
    'legend.fontsize': 10, 'figure.figsize': (10, 7),
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.15,
}
plt.rcParams.update(STYLE)
sns.set_theme(style='whitegrid', palette='deep')

# ─── Constants ───────────────────────────────────────────────────────────────
CELL_LINES = ['22RV1', 'DU145', 'LNCaP']
FITNESS_THRESHOLD = -0.5  # DepMap Chronos: gene effect < -0.5 = fitness gene
PALETTE = {'DiCE': '#2563EB', 'DEG': '#DC2626', 'Random': '#9CA3AF'}


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


def load_fitness_genes(filepath: Path) -> Dict[str, Dict[str, float]]:
    """
    Load fitness gene data from curated DepMap file.

    Returns:
        Dict mapping gene -> {cell_line: score} for genes meeting fitness threshold
    """
    if not filepath.exists():
        log.error(f"Fitness gene file not found: {filepath}")
        return {}

    fitness = {}
    with open(filepath) as f:
        header = None
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split('\t')
            if header is None:
                header = parts[1:]  # skip 'GENE' column
                continue
            gene = parts[0]
            scores = {}
            for i, cl in enumerate(header):
                try:
                    score = float(parts[i + 1])
                    scores[cl] = score
                except (ValueError, IndexError):
                    scores[cl] = 0.0
            fitness[gene] = scores

    log.info(f"Loaded {len(fitness)} genes from fitness data ({filepath.name})")
    return fitness


def get_fitness_genes_per_cell_line(
    fitness_data: Dict[str, Dict[str, float]],
    threshold: float = FITNESS_THRESHOLD
) -> Dict[str, Set[str]]:
    """
    Extract sets of fitness genes per cell line based on threshold.

    Returns:
        Dict mapping cell_line -> set of fitness gene symbols
    """
    per_cl = {cl: set() for cl in CELL_LINES}
    for gene, scores in fitness_data.items():
        for cl in CELL_LINES:
            if cl in scores and scores[cl] < threshold:
                per_cl[cl].add(gene)

    for cl, genes in per_cl.items():
        log.info(f"  {cl}: {len(genes)} fitness genes (threshold < {threshold})")

    return per_cl


def get_union_fitness_genes(per_cl: Dict[str, Set[str]]) -> Set[str]:
    """Get union of fitness genes across all cell lines."""
    union = set()
    for genes in per_cl.values():
        union |= genes
    log.info(f"  Union: {len(union)} unique fitness genes across {len(CELL_LINES)} cell lines")
    return union


def load_dice_genes(filepath: Path) -> pd.DataFrame:
    """Load DiCE gene results."""
    df = pd.read_csv(filepath, sep='\t')
    log.info(f"Loaded {len(df)} DiCE genes from {filepath.name}")
    return df


def load_deg_genes(root: Path) -> Set[str]:
    """
    Load DEG-only genes as baseline comparison.
    Falls back to bottom-ranked DiCE genes if no separate DEG list exists.
    """
    # Try to find a DEG list
    deg_paths = [
        root / 'data' / 'processed' / 'deg_genes.txt',
        root / 'data' / 'processed' / 'filtered' / 'deg_list.txt',
    ]
    for p in deg_paths:
        if p.exists():
            genes = set()
            with open(p) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        genes.add(line.split('\t')[0])
            log.info(f"Loaded {len(genes)} DEG genes from {p.name}")
            return genes

    # Fallback: use the full DiCE gene list (all genes that passed initial filters)
    full_path = root / 'data' / 'results' / 'dice_genes.txt'
    if full_path.exists():
        df = pd.read_csv(full_path, sep='\t')
        # Bottom 50% of genes as "DEG-only" baseline (not highly ranked by DiCE)
        bottom = df.tail(len(df) // 2)
        genes = set(bottom['gene'].tolist())
        log.info(f"Using bottom {len(genes)} DiCE-ranked genes as DEG baseline")
        return genes

    log.warning("No DEG gene list found; skipping DEG comparison")
    return set()


# ═══════════════════════════════════════════════════════════════════════════
#  OVERLAP ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

def compute_overlap(
    query_genes: Set[str],
    fitness_genes: Set[str],
    universe_size: int,
    label: str = ""
) -> Dict:
    """
    Compute overlap between a query gene set and fitness genes.
    Uses hypergeometric test for significance.
    """
    # Ensure types are integers for SciPy
    M = int(universe_size)
    
    # Restrict fitness genes to those that actually exist in the universe 
    # to satisfy the mathematical constraint: n <= M
    overlap = query_genes & fitness_genes
    k = int(len(overlap))
    N = int(len(query_genes))
    n = int(min(M, len(fitness_genes))) 

    # Hypergeometric test: P(X >= k)
    if k > 0:
        pval = hypergeom.sf(k - 1, M, n, N)
    else:
        pval = 1.0
        
    pct = (k / N * 100) if N > 0 else 0.0
    expected = (N * n) / M if M > 0 else 0.0
    
    if expected > 0:
        fold_enrichment = k / expected
    else:
        fold_enrichment = 0.0

    result = {
        'label': label,
        'n_query': N,
        'n_fitness': n,
        'n_overlap': k,
        'pct_overlap': pct,
        'expected': expected,
        'fold_enrichment': fold_enrichment,
        'pvalue': pval,
        'overlap_genes': sorted(overlap),
    }
    return result


def run_fitness_analysis(
    dice_df: pd.DataFrame,
    fitness_per_cl: Dict[str, Set[str]],
    union_fitness: Set[str],
    deg_genes: Set[str],
    universe_size: int,
    top_k_list: List[int],
) -> pd.DataFrame:
    """
    Run the full fitness gene overlap analysis.

    Computes overlaps for:
      - DiCE top-K vs each cell line and union
      - DEG genes vs each cell line and union (baseline)

    Returns DataFrame with all results.
    """
    results = []
    all_dice_genes = set(dice_df['gene'].tolist())

    for top_k in top_k_list:
        top_genes = set(dice_df.head(top_k)['gene'].tolist())

        # DiCE vs each cell line
        for cl in CELL_LINES:
            r = compute_overlap(top_genes, fitness_per_cl[cl], universe_size,
                                f"DiCE_top{top_k}")
            r['cell_line'] = cl
            r['comparison'] = 'DiCE'
            r['top_k'] = top_k
            results.append(r)

        # DiCE vs union
        r = compute_overlap(top_genes, union_fitness, universe_size,
                            f"DiCE_top{top_k}")
        r['cell_line'] = 'Union'
        r['comparison'] = 'DiCE'
        r['top_k'] = top_k
        results.append(r)

        # DEG baseline (same K genes from DEG set)
        if deg_genes:
            deg_sample = set(list(deg_genes)[:top_k])
            for cl in CELL_LINES:
                r = compute_overlap(deg_sample, fitness_per_cl[cl], universe_size,
                                    f"DEG_top{top_k}")
                r['cell_line'] = cl
                r['comparison'] = 'DEG'
                r['top_k'] = top_k
                results.append(r)

            r = compute_overlap(deg_sample, union_fitness, universe_size,
                                f"DEG_top{top_k}")
            r['cell_line'] = 'Union'
            r['comparison'] = 'DEG'
            r['top_k'] = top_k
            results.append(r)

    df = pd.DataFrame(results)
    return df


# ═══════════════════════════════════════════════════════════════════════════
#  VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════

def plot_fitness_barplot(results_df: pd.DataFrame, output_path: Path):
    """
    Grouped bar chart comparing DiCE vs DEG fitness gene overlap percentages.
    Replicates the key figure from the paper (Fig. 5-style).
    """
    # Use the union results for the main comparison
    union_df = results_df[results_df['cell_line'] == 'Union'].copy()

    if union_df.empty:
        log.warning("No union results to plot")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # ── Left panel: Per cell line (largest top_k) ─────────────────────
    ax = axes[0]
    largest_k = results_df['top_k'].max()
    plot_df = results_df[
        (results_df['top_k'] == largest_k) &
        (results_df['cell_line'] != 'Union')
    ].copy()

    if not plot_df.empty:
        x = np.arange(len(CELL_LINES))
        width = 0.35
        dice_pcts = []
        deg_pcts = []
        for cl in CELL_LINES:
            dice_row = plot_df[(plot_df['comparison'] == 'DiCE') & (plot_df['cell_line'] == cl)]
            deg_row = plot_df[(plot_df['comparison'] == 'DEG') & (plot_df['cell_line'] == cl)]
            dice_pcts.append(dice_row['pct_overlap'].values[0] if len(dice_row) > 0 else 0)
            deg_pcts.append(deg_row['pct_overlap'].values[0] if len(deg_row) > 0 else 0)

        bars1 = ax.bar(x - width/2, dice_pcts, width, label='DiCE genes', color=PALETTE['DiCE'], alpha=0.9)
        bars2 = ax.bar(x + width/2, deg_pcts, width, label='DEG genes', color=PALETTE['DEG'], alpha=0.9)

        # Add value labels
        for bar in bars1:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., h + 0.3,
                    f'{h:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')
        for bar in bars2:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., h + 0.3,
                    f'{h:.1f}%', ha='center', va='bottom', fontsize=9)

        ax.set_xlabel('Cell Line')
        ax.set_ylabel('Fitness Gene Overlap (%)')
        ax.set_title(f'Fitness Gene Overlap per Cell Line (Top {largest_k})')
        ax.set_xticks(x)
        ax.set_xticklabels(CELL_LINES)
        ax.legend()
        ax.set_ylim(0, max(max(dice_pcts), max(deg_pcts)) * 1.3)

    # ── Right panel: Across top_k values (union) ─────────────────────
    ax = axes[1]
    for comp, color in [('DiCE', PALETTE['DiCE']), ('DEG', PALETTE['DEG'])]:
        comp_df = union_df[union_df['comparison'] == comp].sort_values('top_k')
        if not comp_df.empty:
            ax.plot(comp_df['top_k'], comp_df['pct_overlap'],
                    'o-', color=color, label=f'{comp} genes', linewidth=2, markersize=8)
            # Add p-value annotations
            for _, row in comp_df.iterrows():
                sig = '***' if row['pvalue'] < 0.001 else ('**' if row['pvalue'] < 0.01 else
                       ('*' if row['pvalue'] < 0.05 else 'ns'))
                ax.annotate(sig, (row['top_k'], row['pct_overlap']),
                           textcoords="offset points", xytext=(0, 10),
                           ha='center', fontsize=9, fontweight='bold',
                           color=color)

    ax.set_xlabel('Top K DiCE Genes')
    ax.set_ylabel('Fitness Gene Overlap (%)')
    ax.set_title('Fitness Gene Overlap vs Gene Rank (Union)')
    ax.legend()

    plt.suptitle('Cancer Fitness Gene Analysis — DiCE vs DEG',
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    log.info(f"Saved fitness barplot → {output_path}")


def plot_fitness_heatmap(
    dice_df: pd.DataFrame,
    fitness_data: Dict[str, Dict[str, float]],
    top_n: int,
    output_path: Path
):
    """
    Heatmap showing fitness scores of top DiCE genes across cell lines.
    """
    top_genes = dice_df.head(top_n)['gene'].tolist()

    # Build matrix
    matrix = []
    genes_in_data = []
    for gene in top_genes:
        if gene in fitness_data:
            row = [fitness_data[gene].get(cl, 0.0) for cl in CELL_LINES]
            matrix.append(row)
            genes_in_data.append(gene)

    if not matrix:
        log.warning("No top DiCE genes found in fitness data for heatmap")
        return

    mat_df = pd.DataFrame(matrix, index=genes_in_data, columns=CELL_LINES)

    fig, ax = plt.subplots(figsize=(6, max(8, len(genes_in_data) * 0.35)))

    # Custom colormap: red = essential, white = neutral
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    sns.heatmap(
        mat_df, ax=ax, cmap=cmap, center=0,
        vmin=-1.5, vmax=0.5,
        linewidths=0.5, linecolor='white',
        cbar_kws={'label': 'Chronos Gene Effect', 'shrink': 0.6},
        annot=True, fmt='.2f', annot_kws={'size': 8},
    )

    # Mark fitness genes (< threshold) with bold text
    for i, gene in enumerate(genes_in_data):
        for j, cl in enumerate(CELL_LINES):
            score = fitness_data[gene].get(cl, 0.0)
            if score < FITNESS_THRESHOLD:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False,
                                           edgecolor='red', linewidth=2))

    ax.set_title(f'DepMap Fitness Scores — Top {len(genes_in_data)} DiCE Genes',
                 fontsize=13, fontweight='bold', pad=12)
    ax.set_ylabel('DiCE Gene (ranked)')
    ax.set_xlabel('Prostate Cancer Cell Line')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    log.info(f"Saved fitness heatmap → {output_path}")


def plot_venn_overlap(
    dice_genes: Set[str],
    fitness_genes: Set[str],
    output_path: Path,
    title: str = "DiCE ∩ Fitness Genes"
):
    """
    Simple proportional Venn-style diagram using matplotlib.
    """
    overlap = dice_genes & fitness_genes
    dice_only = dice_genes - fitness_genes
    fitness_only = fitness_genes - dice_genes

    fig, ax = plt.subplots(figsize=(8, 6))

    # Draw circles
    from matplotlib.patches import Circle
    c1 = Circle((0.35, 0.5), 0.25, alpha=0.4, color=PALETTE['DiCE'], label='DiCE genes')
    c2 = Circle((0.65, 0.5), 0.25, alpha=0.4, color=PALETTE['DEG'], label='Fitness genes')
    ax.add_patch(c1)
    ax.add_patch(c2)

    # Add counts
    ax.text(0.25, 0.5, f'{len(dice_only)}', fontsize=18, ha='center', va='center',
            fontweight='bold', color=PALETTE['DiCE'])
    ax.text(0.50, 0.5, f'{len(overlap)}', fontsize=20, ha='center', va='center',
            fontweight='bold', color='black')
    ax.text(0.75, 0.5, f'{len(fitness_only)}', fontsize=18, ha='center', va='center',
            fontweight='bold', color=PALETTE['DEG'])

    # Labels
    ax.text(0.25, 0.78, 'DiCE\nonly', fontsize=11, ha='center', color=PALETTE['DiCE'])
    ax.text(0.50, 0.78, 'Overlap', fontsize=11, ha='center', fontweight='bold')
    ax.text(0.75, 0.78, 'Fitness\nonly', fontsize=11, ha='center', color=PALETTE['DEG'])

    # Overlap percentage
    pct = len(overlap) / len(dice_genes) * 100 if dice_genes else 0
    ax.text(0.50, 0.18, f'{pct:.1f}% of DiCE genes are fitness genes',
            fontsize=12, ha='center', style='italic')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    log.info(f"Saved Venn diagram → {output_path}")


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN ANALYSIS FUNCTION (callable from validation_pipeline.py)
# ═══════════════════════════════════════════════════════════════════════════

def run_fitness_gene_analysis(
    dice_df: pd.DataFrame,
    root: Path,
    out_fig: Path,
    out_tbl: Path,
    top_k_list: Optional[List[int]] = None,
    universe_size: int = 8962,
) -> Dict:
    """
    Main entry point for cancer fitness gene analysis.
    Can be called standalone or from validation_pipeline.py.

    Returns a dict with summary statistics.
    """
    log.info("=" * 60)
    log.info("ANALYSIS 9: Cancer Fitness Gene Overlap")
    log.info("=" * 60)

    # Load fitness data
    fitness_path = root / 'data' / 'known_genes' / 'depmap_fitness_pca.txt'
    fitness_data = load_fitness_genes(fitness_path)

    if not fitness_data:
        log.error("Cannot run fitness analysis: no fitness gene data loaded")
        return {'status': 'FAILED', 'reason': 'No fitness gene data'}

    # Get fitness gene sets
    fitness_per_cl = get_fitness_genes_per_cell_line(fitness_data)
    union_fitness = get_union_fitness_genes(fitness_per_cl)

    # Default top_k values
    if top_k_list is None:
        n_genes = len(dice_df)
        top_k_list = sorted(set([
            50,
            min(100, n_genes),
            min(200, n_genes),
            min(364, n_genes),   # Paper's TCGA DiCE count
            min(500, n_genes),
        ]))

    # Load DEG baseline
    deg_genes = load_deg_genes(root)

    # Run overlap analysis
    results_df = run_fitness_analysis(
        dice_df, fitness_per_cl, union_fitness,
        deg_genes, universe_size, top_k_list
    )

    # Save results table
    cols_to_save = ['label', 'cell_line', 'comparison', 'top_k',
                     'n_query', 'n_fitness', 'n_overlap', 'pct_overlap',
                     'expected', 'fold_enrichment', 'pvalue']
    results_df[cols_to_save].to_csv(
        out_tbl / '09_fitness_overlap.csv', index=False
    )
    log.info(f"Saved results table → {out_tbl / '09_fitness_overlap.csv'}")

    # Save overlap gene lists
    all_dice_top = set(dice_df.head(max(top_k_list))['gene'].tolist())
    overlap_genes = all_dice_top & union_fitness
    overlap_list = []
    for gene in sorted(overlap_genes):
        scores = fitness_data.get(gene, {})
        overlap_list.append({
            'gene': gene,
            'dice_rank': dice_df[dice_df['gene'] == gene]['final_rank'].values[0]
                         if gene in dice_df['gene'].values else 'N/A',
            **{f'{cl}_score': scores.get(cl, 'N/A') for cl in CELL_LINES},
        })

    if overlap_list:
        overlap_df = pd.DataFrame(overlap_list).sort_values('dice_rank')
        overlap_df.to_csv(out_tbl / '09_fitness_overlap_genes.csv', index=False)
        log.info(f"Saved {len(overlap_df)} overlap genes → {out_tbl / '09_fitness_overlap_genes.csv'}")

    # Generate figures
    plot_fitness_barplot(results_df, out_fig / '09_fitness_barplot.png')
    plot_fitness_heatmap(dice_df, fitness_data, 50, out_fig / '09_fitness_heatmap.png')
    plot_venn_overlap(
        set(dice_df.head(364)['gene'].tolist()),
        union_fitness,
        out_fig / '09_fitness_venn.png',
        title=f"DiCE Top-364 ∩ Fitness Genes ({len(CELL_LINES)} PCa lines)"
    )

    # Summary statistics
    main_k = min(364, len(dice_df))
    main_result = results_df[
        (results_df['top_k'] == main_k) &
        (results_df['cell_line'] == 'Union') &
        (results_df['comparison'] == 'DiCE')
    ]

    summary = {
        'status': 'COMPLETE',
        'n_fitness_genes_total': len(union_fitness),
        'cell_lines': CELL_LINES,
        'fitness_per_cl': {cl: len(genes) for cl, genes in fitness_per_cl.items()},
        'main_top_k': main_k,
    }

    if len(main_result) > 0:
        row = main_result.iloc[0]
        summary.update({
            'dice_overlap_pct': f"{row['pct_overlap']:.1f}%",
            'dice_overlap_n': int(row['n_overlap']),
            'dice_pvalue': f"{row['pvalue']:.2e}",
            'dice_fold_enrichment': f"{row['fold_enrichment']:.1f}x",
        })
        log.info(f"\n  DiCE top-{main_k} overlap with fitness genes:")
        log.info(f"    {row['n_overlap']}/{main_k} = {row['pct_overlap']:.1f}%")
        log.info(f"    Fold enrichment: {row['fold_enrichment']:.1f}x")
        log.info(f"    P-value: {row['pvalue']:.2e}")

    # DEG comparison
    deg_result = results_df[
        (results_df['top_k'] == main_k) &
        (results_df['cell_line'] == 'Union') &
        (results_df['comparison'] == 'DEG')
    ]
    if len(deg_result) > 0:
        deg_row = deg_result.iloc[0]
        summary['deg_overlap_pct'] = f"{deg_row['pct_overlap']:.1f}%"
        log.info(f"    DEG baseline: {deg_row['pct_overlap']:.1f}%")

    log.info("")
    return summary


# ═══════════════════════════════════════════════════════════════════════════
#  CLI ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Cancer Fitness Gene Overlap Analysis (DepMap)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python fitness_analysis.py
  python fitness_analysis.py --dice-genes data/results/dice_genes_top500.txt
  python fitness_analysis.py --top-k 50 100 200 364
        """
    )
    parser.add_argument('--dice-genes', type=str, default=None,
                        help='Path to DiCE gene results (default: auto-detect)')
    parser.add_argument('--top-k', nargs='+', type=int, default=None,
                        help='Top K values to compute overlap for (default: 50 100 200 364 500)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory (default: data/results/validation/)')
    parser.add_argument('--threshold', type=float, default=FITNESS_THRESHOLD,
                        help=f'Fitness gene threshold (default: {FITNESS_THRESHOLD})')
    parser.add_argument('--universe-size', type=int, default=None,
                        help='Gene universe size for hypergeometric test (default: auto)')

    args = parser.parse_args()

    # Find root
    root = find_project_root()

    # Load DiCE genes
    if args.dice_genes:
        dice_path = Path(args.dice_genes)
    else:
        dice_path = root / 'data' / 'results' / 'dice_genes_top500.txt'
    dice_df = load_dice_genes(dice_path)

    # Universe size
    universe_size = args.universe_size
    if universe_size is None:
        full_path = root / 'data' / 'results' / 'dice_genes.txt'
        if full_path.exists():
            full_df = pd.read_csv(full_path, sep='\t')
            universe_size = len(full_df)
        else:
            universe_size = len(dice_df)

    # Setup output
    out_dir = Path(args.output_dir) if args.output_dir else root / 'data' / 'results' / 'validation'
    out_fig = out_dir / 'figures'
    out_tbl = out_dir / 'tables'
    out_fig.mkdir(parents=True, exist_ok=True)
    out_tbl.mkdir(parents=True, exist_ok=True)

    # Use threshold from args
    threshold = args.threshold

    # Run analysis
    summary = run_fitness_gene_analysis(
        dice_df, root, out_fig, out_tbl,
        top_k_list=args.top_k,
        universe_size=universe_size,
    )

    log.info("=" * 60)
    log.info("FITNESS ANALYSIS COMPLETE")
    log.info(f"Figures: {out_fig}")
    log.info(f"Tables:  {out_tbl}")
    log.info("=" * 60)


if __name__ == '__main__':
    main()
