#!/usr/bin/env python3
"""
Post-Analysis Utilities (Updated for Aim 1)
Enrichment analysis and visualization for DiCE results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DiCEVisualizer:
    """Visualization tools for DiCE-Flow Aim 1 results"""
    
    def __init__(self, dice_file: str, output_dir: str = 'data/results/figures'):
        self.dice_df = pd.read_csv(dice_file, sep='\t')
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        sns.set_theme(style="whitegrid", palette="husl")
    
    def plot_top_genes(self, n: int = 20):
        """Bar plot of top N DiCE genes"""
        fig, ax = plt.subplots(figsize=(10, 8))
        top_n = self.dice_df.head(n).iloc[::-1]  # Reverse for bottom-to-top plotting
        
        bars = ax.barh(range(len(top_n)), top_n['ensemble_score'], color='teal')
        ax.set_yticks(range(len(top_n)))
        ax.set_yticklabels(top_n['gene'])
        
        # Color bars by direction
        for i, bar in enumerate(bars):
            if top_n.iloc[i]['direction'] < 0:
                bar.set_color('#DC2626') # Red for localizing
            else:
                bar.set_color('#2563EB') # Blue for bridging
                
        ax.set_xlabel('Magnitude of Structural Shift ($\Delta$)', fontsize=12)
        ax.set_title(f'Top {n} Repositioned Genes', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        output = self.output_dir / f'top_{n}_genes_aim1.png'
        plt.savefig(output, dpi=300, bbox_inches='tight')
        logger.info(f"Saved plot: {output}")
        plt.close()
    
    def plot_rank_distribution(self):
        """Distribution of structural shift scores"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        ax1.hist(self.dice_df['ensemble_score'], bins=50, color='purple', alpha=0.7)
        ax1.set_xlabel('Structural Shift Magnitude', fontsize=12)
        ax1.set_ylabel('Frequency', fontsize=12)
        ax1.set_title('Distribution of Shifts', fontsize=14, fontweight='bold')
        ax1.set_yscale('log')
        
        ax2.plot(range(1, len(self.dice_df) + 1), self.dice_df['ensemble_score'].values, color='purple')
        ax2.set_xlabel('Rank', fontsize=12)
        ax2.set_ylabel('Structural Shift Magnitude', fontsize=12)
        ax2.set_title('Magnitude by Rank', fontsize=14, fontweight='bold')
        ax2.set_yscale('log')
        
        plt.tight_layout()
        output = self.output_dir / 'score_distribution_aim1.png'
        plt.savefig(output, dpi=300, bbox_inches='tight')
        logger.info(f"Saved plot: {output}")
        plt.close()
    
    def generate_all_plots(self):
        logger.info("Generating Aim 1 visualizations...")
        self.plot_top_genes(20)
        self.plot_rank_distribution()
        logger.info(f"All plots saved to {self.output_dir}/")

class EnrichmentAnalyzer:
    def __init__(self, dice_file: str):
        self.dice_df = pd.read_csv(dice_file, sep='\t')
    
    def export_for_enrichr(self, output_file: str, top_n: int = 500):
        top_genes = self.dice_df.head(top_n)['gene'].tolist()
        with open(output_file, 'w') as f:
            for gene in top_genes:
                f.write(f"{gene}\n")
        logger.info(f"Exported {len(top_genes)} genes for Enrichr: {output_file}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description='DiCE post-analysis utilities')
    parser.add_argument('--dice-file', required=True, help='DiCE genes file')
    parser.add_argument('--output-dir', default='data/results/figures')
    parser.add_argument('--visualize', action='store_true')
    parser.add_argument('--export-enrichr', type=str)
    args = parser.parse_args()
    
    if args.visualize:
        viz = DiCEVisualizer(args.dice_file, args.output_dir)
        viz.generate_all_plots()
    
    if args.export_enrichr:
        enrichment = EnrichmentAnalyzer(args.dice_file)
        enrichment.export_for_enrichr(args.export_enrichr)

if __name__ == '__main__':
    main()