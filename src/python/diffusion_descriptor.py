import torch
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as splinalg
import logging
import multiprocessing
import warnings

# Silence the PyTorch sparse tensor warning
warnings.filterwarnings("ignore", message=".*Sparse invariant checks.*")

logger = logging.getLogger(__name__)

# Cache the device so it only prints once!
_BEST_DEVICE = None

def get_best_device():
    global _BEST_DEVICE
    if _BEST_DEVICE is not None:
        return _BEST_DEVICE
        
    if torch.cuda.is_available():
        device_name = torch.cuda.get_device_name(0)
        logger.info(f"Hardware Acceleration Enabled: GPU detected -> {device_name}")
        _BEST_DEVICE = torch.device('cuda')
    
    elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
        logger.info("Hardware Acceleration Enabled: Apple Silicon (MPS) detected")
        _BEST_DEVICE = torch.device('mps')
    
    else:
        cores = multiprocessing.cpu_count()
        torch.set_num_threads(cores)
        logger.info(f"No GPU detected. Falling back to multi-core CPU ({cores} threads).")
        _BEST_DEVICE = torch.device('cpu')
        
    return _BEST_DEVICE

def compute_laplacians(adj_matrix):
    degrees = np.array(adj_matrix.sum(axis=1)).flatten()
    degrees[degrees == 0] = 1e-10 
    
    D = sp.diags(degrees, format="csr")
    D_inv = sp.diags(1.0 / degrees, format="csr")
    
    L_c = D - adj_matrix
    I = sp.identity(adj_matrix.shape[0], format="csr")
    P_t = 0.5 * (I + D_inv @ adj_matrix)
    
    return L_c, P_t, degrees

def hks_diag_cheb_torch(L, lam_max, ts, K=60, S=512, seed=0):
    """
    Batched Generalized Trace Estimator using PyTorch.
    Seamlessly adapts to NVIDIA, AMD (ROCm), Mac (MPS), or Multi-core CPU.
    """
    device = get_best_device()
    N = L.shape[0]
    
    # 1. Compute Chebyshev Coefficients (Fast enough on CPU)
    j = np.arange(K + 1)
    theta = np.pi * (j + 0.5) / (K + 1)
    lam_grid = (lam_max / 2) * (np.cos(theta) + 1)
    
    C = np.empty((len(ts), K + 1))
    for a, t in enumerate(ts):
        f = np.exp(-t * lam_grid)
        C[a] = (2 / (K + 1)) * (f[None, :] * np.cos(np.outer(j, theta))).sum(1)
        C[a, 0] *= 0.5
    C_pt = torch.tensor(C, dtype=torch.float32, device=device)
    
    # 2. Move Laplacian to PyTorch Sparse Tensor
    L_coo = L.tocoo()
    indices = torch.tensor(np.vstack((L_coo.row, L_coo.col)), dtype=torch.long, device=device)
    values = torch.tensor(L_coo.data, dtype=torch.float32, device=device)
    L_pt = torch.sparse_coo_tensor(indices, values, L_coo.shape).coalesce()
    
    # Identity matrix on target device
    I_pt = torch.sparse_coo_tensor(
        torch.arange(N).unsqueeze(0).repeat(2, 1).to(device),
        torch.ones(N, device=device),
        (N, N)
    ).coalesce()
    
    # Normalized Laplacian [-1, 1]
    Lt_pt = (2.0 / lam_max) * L_pt - I_pt
    
    # 3. Batched Hutchinson Trace Estimation
    torch.manual_seed(seed)
    
    # Generate S random vectors at once: shape (N, S)
    z = torch.randint(0, 2, (N, S), device=device, dtype=torch.float32) * 2.0 - 1.0
    
    Tm2 = z
    Tm1 = torch.sparse.mm(Lt_pt, z)
    
    acc = torch.zeros((len(ts), N), device=device)
    
    # Diag expectation across S samples
    diag_k0 = (z * Tm2).mean(dim=1)
    diag_k1 = (z * Tm1).mean(dim=1)
    
    acc += torch.ger(C_pt[:, 0], diag_k0) + torch.ger(C_pt[:, 1], diag_k1)
    
    # Recurrence relation
    for k in range(2, K + 1):
        Tk = 2.0 * torch.sparse.mm(Lt_pt, Tm1) - Tm2
        diag_k = (z * Tk).mean(dim=1)
        acc += torch.ger(C_pt[:, k], diag_k)
        
        Tm2, Tm1 = Tm1, Tk
        
    return acc.cpu().numpy()

def get_scale_bank(L_c):
    """
    Computes the diffusion scale bank. 
    Bypasses exact Fiedler value calculation to avoid 30+ minute sparse LU decompositions.
    """
    # 1. Largest Eigenvalue (Standard Arnoldi iteration - very fast)
    lam_max = splinalg.eigsh(L_c, k=1, which='LM', return_eigenvectors=False)[0]
    
    # 2. Second Smallest (Fiedler value heuristic)
    lam_2 = 1e-4
        
    t_min = 4 * np.log(10) / lam_max
    t_max = 4 * np.log(10) / lam_2
    
    # Create 16 logarithmic steps
    ts = np.logspace(np.log10(t_min), np.log10(t_max), 16)
    return ts, lam_max

def get_lazy_walk_profile(P_t, k_steps=[2, 4, 8, 16]):
    N = P_t.shape[0]
    rwse = np.zeros((len(k_steps), N))
    P_curr = P_t.copy()
    current_k = 1
    
    for i, target_k in enumerate(k_steps):
        while current_k < target_k:
            P_curr = P_curr @ P_t
            current_k += 1
        rwse[i] = P_curr.diagonal()
    return rwse