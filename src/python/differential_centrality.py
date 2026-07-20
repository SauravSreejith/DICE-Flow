#!/usr/bin/env python3
"""
Differential Centrality Analysis
Implements DiCE Phase V: Ensemble Ranking based on Absolute Differences
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
import argparse

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DifferentialCentralityAnalyzer:
    def __init__(self, normal_centrality: str, tumor_centrality: str):
        self.normal_file = Path(normal_centrality)
        self.tumor_file = Path(tumor_centrality)
    
    def analyze(self) -> pd.DataFrame:
        logger.info("Loading centrality data...")
        try:
            normal_df = pd.read_csv(self.normal_file, sep='\t')
            tumor_df = pd.read_csv(self.tumor_file, sep='\t')
        except Exception as e:
            logger.error(f"Error reading files: {e}")
            return pd.DataFrame()
        
        # Merge
        merged = pd.merge(normal_df, tumor_df, on='gene', suffixes=('_norm', '_tumor'))
        
        # 1. Calculate Absolute Differences (DiCE Metric)
        # We care about magnitude of change, whether up or down
        merged['diff_betweenness'] = (merged['betweenness_tumor'] - merged['betweenness_norm']).abs()
        merged['diff_eigenvector'] = (merged['eigenvector_tumor'] - merged['eigenvector_norm']).abs()
        
        # 2. Ranking (0 to 1)
        # Rank the differences. pct=True gives a score 0.0 to 1.0
        merged['rank_bet'] = merged['diff_betweenness'].rank(pct=True)
        merged['rank_eig'] = merged['diff_eigenvector'].rank(pct=True)
        
        # 3. Ensemble Score (Product of Ranks)
        # This penalizes genes that only change in one metric but not the other
        merged['ensemble_score'] = merged['rank_bet'] * merged['rank_eig']
        
        # 4. Filter Noise (Phase V Filtering)
        # 4a. Remove extreme high-eigenvector nodes (Top 12%)
        # The paper explicitly notes these represent dense, non-differential 
        # noise clusters (housekeeping core) and must be removed.
        eig_threshold_norm = merged['eigenvector_norm'].quantile(0.88)  # Bottom 88%
        eig_threshold_tumor = merged['eigenvector_tumor'].quantile(0.88)
        
        mask_not_hub = (merged['eigenvector_norm'] <= eig_threshold_norm) & \
                       (merged['eigenvector_tumor'] <= eig_threshold_tumor)

        # 4b. Remove extreme low-betweenness nodes (Peripheral Noise)
        # Exclude genes that have very low centrality in BOTH conditions
        # (They are peripheral in both networks, changes are likely noise)
        avg_bet_norm = merged['betweenness_norm'].mean()
        avg_bet_tumor = merged['betweenness_tumor'].mean()
        
        # Keep if gene is above average in AT LEAST one condition
        mask_not_periph = (merged['betweenness_norm'] > avg_bet_norm) | \
                          (merged['betweenness_tumor'] > avg_bet_tumor)
               
        # Apply combined strict filtering mask
        if getattr(self, 'disable_noise', False):
            filtered_df = merged.copy()
            logger.info("NOISE FILTER DISABLED: Retaining ultra-hubs like EP300 and peripheral nodes.")
        else:
            filtered_df = merged[mask_not_hub & mask_not_periph].copy()
            
            logger.info(f"Filtering: Started with {len(merged)} genes.")
            logger.info(f"  - Removed {len(merged) - mask_not_hub.sum()} high-eigenvector noise genes")
            logger.info(f"  - Removed {len(merged) - mask_not_periph.sum()} low-betweenness peripheral genes")
            logger.info(f"  - Final retained pool: {len(filtered_df)} genes")
        
        # 5. Sort by Ensemble Score (Descending)
        filtered_df = filtered_df.sort_values('ensemble_score', ascending=False)
        
        # Add rank column (1 = top gene)
        filtered_df['final_rank'] = range(1, len(filtered_df) + 1)

        # Add delta columns for visualization (keep sign here to show up/down)
        filtered_df['delta_betweenness'] = filtered_df['betweenness_tumor'] - filtered_df['betweenness_norm']
        filtered_df['delta_eigenvector'] = filtered_df['eigenvector_tumor'] - filtered_df['eigenvector_norm']

        logger.info(f"Ranking complete. Top gene: {filtered_df.iloc[0]['gene']} (Score: {filtered_df.iloc[0]['ensemble_score']:.4f})")
        
        return filtered_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--normal', required=True)
    parser.add_argument('--tumor', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--top-n', type=int, default=500)
    parser.add_argument('--disable-noise', action='store_true', help='Skip Phase 5 12 hub removal')
    
    args = parser.parse_args()
    
    analyzer = DifferentialCentralityAnalyzer(args.normal, args.tumor)
    if getattr(args, 'disable_noise', False):
        analyzer.disable_noise = True
    else:
        analyzer.disable_noise = False
        
    dice_df = analyzer.analyze()
    
    if dice_df.empty:
        logger.error("Analysis produced no results.")
        return

    # Save full
    dice_df.to_csv(args.output, index=False, sep='\t')
    
    # Save top N
    top_output = args.output.replace('.txt', f'_top{args.top_n}.txt')
    dice_df.head(args.top_n).to_csv(top_output, index=False, sep='\t')
    logger.info(f"Saved Top {args.top_n} to {top_output}")

if __name__ == '__main__':
    main()