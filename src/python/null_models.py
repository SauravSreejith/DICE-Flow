#!/usr/bin/env python3
"""
Statistical Calibration and Null Models
Implements Degree-Stratified FDR control and Per-Scale Centering for Aim 1.
"""

import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import logging

logger = logging.getLogger(__name__)

def stratified_pvals(Delta, Delta_null, degree, min_pool_size=400000):
    """
    Computes degree-stratified robust p-values.
    Uses Median Absolute Deviation (MAD) to prevent random trace-estimator 
    outliers from artificially capping the empirical p-value.
    """
    logger.info("Computing degree-stratified robust p-values...")
    
    n_genes = len(Delta)
    n_perms = Delta_null.shape[0]
    
    max_bins = max(1, int((n_genes * n_perms) / min_pool_size))
    n_bins = min(10, max_bins)
    logger.info(f"Using {n_bins} degree strata for robust null pooling.")
    
    degree_jitter = degree + np.random.uniform(0, 1e-5, size=degree.shape)
    edges = np.quantile(degree_jitter, np.linspace(0, 1, n_bins + 1)[1:-1])
    bins = np.digitize(degree_jitter, edges)
    
    p = np.ones_like(Delta)
    
    for b in np.unique(bins):
        idx = np.where(bins == b)[0]
        pool = Delta_null[:, idx].ravel()
        
        if len(pool) == 0:
            continue
            
        # Robust Statistics: Median and MAD ignore extreme trace estimator noise spikes
        median_null = np.median(pool)
        mad_null = np.median(np.abs(pool - median_null))
        
        # Convert MAD to standard deviation (1.4826 is the constant for normal dists)
        sigma_null = mad_null * 1.4826 + 1e-8
        
        # Calculate Robust Z-scores
        z_scores = (Delta[idx] - median_null) / sigma_null
        
        # Right-tail survival function (p-value)
        p[idx] = stats.norm.sf(z_scores)
        
    # Benjamini-Hochberg FDR correction
    q = multipletests(p, method="fdr_bh")[1]
    
    return p, q

def get_differential_profile(H_norm, H_tumor):
    """
    Computes the differential profile d = H_D - H_C.
    Implements Per-Scale Centering & Z-scoring to neutralize global density drift.
    """
    H_norm_log = np.log(np.clip(np.abs(H_norm), a_min=1e-8, a_max=None))
    H_tumor_log = np.log(np.clip(np.abs(H_tumor), a_min=1e-8, a_max=None))
    
    mean_N = np.mean(H_norm_log, axis=1, keepdims=True)
    std_N = np.std(H_norm_log, axis=1, keepdims=True) + 1e-8
    H_norm_z = (H_norm_log - mean_N) / std_N
    
    mean_T = np.mean(H_tumor_log, axis=1, keepdims=True)
    std_T = np.std(H_tumor_log, axis=1, keepdims=True) + 1e-8
    H_tumor_z = (H_tumor_log - mean_T) / std_T
    
    d = H_tumor_z - H_norm_z
    
    # Magnitude (Max absolute shift across scales - isolates characteristic signal)
    Delta = np.max(np.abs(d), axis=0)
    
    kstar = np.argmax(np.abs(d), axis=0)
    sign = np.sign(d[kstar, np.arange(d.shape[1])])
    
    return d, Delta, kstar, sign