#!/usr/bin/env python3
"""
DiCE-Flow Aim 1 Orchestrator (Vectorized & GPU Accelerated)
Executes the Multi-Scale Diffusion Descriptor and Null Model Calibration.
"""

import pandas as pd
import numpy as np
import scipy.sparse as sp
import networkx as nx
from pathlib import Path
import logging
import argparse
import time

from diffusion_descriptor import compute_laplacians, hks_diag_cheb_torch, get_scale_bank, get_lazy_walk_profile
from null_models import stratified_pvals, get_differential_profile

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger("Aim1_Pipeline")

def load_and_intersect(expr_norm_path, expr_tum_path, network_path):
    """Loads expression and network, using fast Pandas vectorization."""
    logger.info("Loading expression matrices...")
    expr_n = pd.read_csv(expr_norm_path, sep="\t", index_col=0)
    expr_t = pd.read_csv(expr_tum_path, sep="\t", index_col=0)
    
    logger.info("Loading network edges...")
    edges = pd.read_csv(network_path, sep="\t")
    
    # 1. Fast Set Intersection
    genes = set(expr_n.index) & set(expr_t.index) & set(edges['source']) & set(edges['target'])
    genes = sorted(list(genes))
    
    expr_n = expr_n.loc[genes]
    expr_t = expr_t.loc[genes]
    
    # 2. FAST Vectorized Graph Building (No iterrows!)
    logger.info("Extracting network topology...")
    mask = edges['source'].isin(genes) & edges['target'].isin(genes)
    filtered_edges = edges[mask]
    
    G = nx.from_pandas_edgelist(filtered_edges, 'source', 'target')
    G.add_nodes_from(genes) # Ensure isolated nodes aren't lost
            
    # Keep the largest connected component
    lcc = max(nx.connected_components(G), key=len)
    genes_lcc = sorted(list(lcc))
    logger.info(f"Retained {len(genes_lcc)} genes in Largest Connected Component.")
    
    expr_n = expr_n.loc[genes_lcc]
    expr_t = expr_t.loc[genes_lcc]
    
    return expr_n, expr_t, genes_lcc, G.subgraph(genes_lcc).copy()

def build_weighted_adjacency(expr_matrix, G_topology, genes):
    """Builds the condition-specific adjacency matrix using vectorized NumPy."""
    N = len(genes)
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    
    # Fast Correlation
    data = expr_matrix.values
    
    # Ignore the divide-by-zero warnings from zero-variance genes
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        corr = np.corrcoef(data)
        
    # FIX: Safely replace NaNs (from zero-division) with 0.0
    corr = np.nan_to_num(corr, nan=0.0)
    
    np.fill_diagonal(corr, 0)
    affinity = np.abs(corr) 
    
    # FAST Vectorized Edge Mapping (No for-loops!)
    edge_list = list(G_topology.edges())
    if not edge_list:
        return sp.csr_matrix((N, N))
        
    u_nodes, v_nodes = zip(*edge_list)
    u_idx = np.array([gene_to_idx[u] for u in u_nodes])
    v_idx = np.array([gene_to_idx[v] for v in v_nodes])
    
    # Extract weights and symmetrize
    weights = affinity[u_idx, v_idx]
    
    row = np.concatenate([u_idx, v_idx])
    col = np.concatenate([v_idx, u_idx])
    weight_data = np.concatenate([weights, weights])
            
    return sp.csr_matrix((weight_data, (row, col)), shape=(N, N))
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--norm-expr", default="data/raw/expression/expression_normal.txt")
    parser.add_argument("--tum-expr", default="data/raw/expression/expression_tumor.txt")
    parser.add_argument("--network", default="data/processed/networks/merged_network.txt")
    parser.add_argument("--n-perms", type=int, default=50)
    args = parser.parse_args()

    start_time = time.time()

    # 1. Setup
    expr_n, expr_t, genes, G_topo = load_and_intersect(args.norm_expr, args.tum_expr, args.network)
    
    n_samples = min(expr_n.shape[1], expr_t.shape[1])
    expr_n = expr_n.iloc[:, :n_samples]
    expr_t = expr_t.iloc[:, :n_samples]
    
    # 2. Build Networks
    logger.info("Building condition-specific weighted graphs...")
    A_norm = build_weighted_adjacency(expr_n, G_topo, genes)
    A_tumor = build_weighted_adjacency(expr_t, G_topo, genes)
    
    # 3. Laplacians
    L_norm, P_norm, deg_norm = compute_laplacians(A_norm)
    L_tumor, P_tumor, deg_tumor = compute_laplacians(A_tumor)
    
    ts, lam_max = get_scale_bank(L_norm + L_tumor)
    
    # 4. GPU-Accelerated Descriptor 
    logger.info("Computing HKS (GPU Accelerated)...")
    H_norm = hks_diag_cheb_torch(L_norm, lam_max, ts, K=60, S=512)
    H_tumor = hks_diag_cheb_torch(L_tumor, lam_max, ts, K=60, S=512)
    
    # 5. Differential Profile
    _, Delta, kstar, sign = get_differential_profile(H_norm, H_tumor)
    
    # 6. Null Models 
    logger.info(f"Running {args.n_perms} null permutations for FDR control...")
    merged_expr = pd.concat([expr_n, expr_t], axis=1)
    Delta_null = np.zeros((args.n_perms, len(genes)))
    n_cols = expr_n.shape[1]
    
    for perm in range(args.n_perms):
        shuffled_cols = np.random.permutation(merged_expr.columns)
        
        expr_n_null = merged_expr[shuffled_cols[:n_cols]]
        expr_t_null = merged_expr[shuffled_cols[n_cols:]]
        
        A_n_null = build_weighted_adjacency(expr_n_null, G_topo, genes)
        A_t_null = build_weighted_adjacency(expr_t_null, G_topo, genes)
        
        L_n_null, _, _ = compute_laplacians(A_n_null)
        L_t_null, _, _ = compute_laplacians(A_t_null)
        
        # Lower S for speed in nulls
        H_n_null = hks_diag_cheb_torch(L_n_null, lam_max, ts, S=128) 
        H_t_null = hks_diag_cheb_torch(L_t_null, lam_max, ts, S=128)
        
        _, D_null, _, _ = get_differential_profile(H_n_null, H_t_null)
        Delta_null[perm, :] = D_null
        if (perm + 1) % 10 == 0:
            logger.info(f"  Completed permutation {perm + 1}/{args.n_perms}")

    # 7. FDR
    avg_degree = (deg_norm + deg_tumor) / 2
    p_val, q_val = stratified_pvals(Delta, Delta_null, avg_degree, n_bins=10)
    
    # 8. Export
    results = pd.DataFrame({
        'gene': genes,
        'ensemble_score': Delta, 
        'characteristic_scale': ts[kstar],
        'direction': sign,
        'p_value': p_val,
        'q_value_fdr': q_val
    })
    
    results = results.sort_values('ensemble_score', ascending=False)
    results['final_rank'] = range(1, len(results) + 1)
    
    out_dir = Path("data/results")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    results.to_csv(out_dir / "dice_genes.txt", sep="\t", index=False)
    results.head(500).to_csv(out_dir / "dice_genes_top500.txt", sep="\t", index=False)
    
    sig_count = (results['q_value_fdr'] < 0.05).sum()
    logger.info(f"Pipeline Complete in {time.time() - start_time:.1f}s! Identified {sig_count} genes (q < 0.05).")

if __name__ == "__main__":
    main()