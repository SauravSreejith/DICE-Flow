#!/usr/bin/env python3
"""
Statistical Calibration and Null Models
Implements Degree-Stratified FDR control for Aim 1.
"""

import numpy as np
from statsmodels.stats.multitest import multipletests
import logging

logger = logging.getLogger(__name__)

def stratified_pvals(Delta, Delta_null, degree, n_bins=10):
    """
    Computes degree-stratified empirical p-values and applies Benjamini-Hochberg FDR.
    (Transcribed exactly from the Professor's proposal - Page 15/16)
    """
    logger.info("Computing degree-stratified empirical p-values...")
    
    # Create degree bins using quantiles
    edges = np.quantile(degree, np.linspace(0, 1, n_bins + 1)[1:-1])
    bins = np.digitize(degree, edges)
    
    p = np.empty_like(Delta)
    
    for b in np.unique(bins):
        idx = np.where(bins == b)[0]
        
        # Pool all null Delta values within this specific degree bin
        pool = np.sort(Delta_null[:, idx].ravel())
        
        if len(pool) == 0:
            p[idx] = 1.0
            continue
            
        # Find rank of actual observation in the null distribution
        rank = np.searchsorted(pool, Delta[idx], side="right")
        
        # Right tail p-value calculation
        p[idx] = (len(pool) - rank + 1) / (len(pool) + 1)
        
    # Benjamini-Hochberg FDR correction
    q = multipletests(p, method="fdr_bh")[1]
    
    return p, q

def get_differential_profile(H_norm, H_tumor):
    """
    Computes the differential profile d = H_D - H_C
    Returns Magnitude, Characteristic Scale, and Direction.
    """
    # Centering (log space recommended by proposal to handle global density drift)
    # FIX: Use np.clip to strictly enforce a minimum value, preventing log(negative)
    H_norm_safe = np.clip(H_norm, a_min=1e-12, a_max=None)
    H_tumor_safe = np.clip(H_tumor, a_min=1e-12, a_max=None)
    
    H_norm_log = np.log(H_norm_safe)
    H_tumor_log = np.log(H_tumor_safe)
    
    # d = H_D - H_C (shape: scales x genes)
    d = H_tumor_log - H_norm_log
    
    # Magnitude: Delta_i = ||d_i|| (Euclidean norm across scales)
    Delta = np.linalg.norm(d, axis=0)
    
    # Characteristic scale per gene
    kstar = np.argmax(np.abs(d), axis=0)
    
    # Direction of movement
    sign = np.sign(d[kstar, np.arange(d.shape[1])])
    
    return d, Delta, kstar, sign