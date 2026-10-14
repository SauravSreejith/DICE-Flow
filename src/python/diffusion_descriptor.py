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

# Cache the device so it only prints once
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

def hks_diag_cheb_torch(L, lam_max, ts, K=60, S=256, seed=0):
    """
    Batched Generalized Trace Estimator using PyTorch.
    Uses paired Common Random Numbers (via `seed`) across Normal/Tumor pairs
    to cancel stochastic variance at mid-to-large diffusion scales.
    """
    device = get_best_device()
    N = L.shape[0]
    
    # 1. Compute Chebyshev Coefficients on CPU
    j = np.arange(K + 1)
    theta = np.pi * (j + 0.5) / (K + 1)
    lam_grid = (lam_max / 2.0) * (np.cos(theta) + 1.0)
    
    C = np.empty((len(ts), K + 1))
    for a, t in enumerate(ts):
        f = np.exp(-t * lam_grid)
        C[a] = (2.0 / (K + 1)) * (f[None, :] * np.cos(np.outer(j, theta))).sum(1)
        C[a, 0] *= 0.5
    C_pt = torch.tensor(C, dtype=torch.float32, device=device)
    
    # 2. Move Laplacian to PyTorch Sparse Tensor
    L_coo = L.tocoo()
    indices = torch.tensor(np.vstack((L_coo.row, L_coo.col)), dtype=torch.long, device=device)
    values = torch.tensor(L_coo.data, dtype=torch.float32, device=device)
    L_pt = torch.sparse_coo_tensor(indices, values, L_coo.shape).coalesce()
    
    I_pt = torch.sparse_coo_tensor(
        torch.arange(N, device=device).unsqueeze(0).repeat(2, 1),
        torch.ones(N, device=device, dtype=torch.float32),
        (N, N)
    ).coalesce()
    
    # Rescale spectrum to [-1, 1]
    Lt_pt = (2.0 / lam_max) * L_pt - I_pt
    
    # 3. Batched Hutchinson Trace Estimation with explicit seed
    torch.manual_seed(seed)
    if device.type == 'cuda':
        torch.cuda.manual_seed(seed)
        
    z = torch.randint(0, 2, (N, S), device=device, dtype=torch.float32) * 2.0 - 1.0
    
    Tm2 = z
    Tm1 = torch.sparse.mm(Lt_pt, z)
    
    acc = torch.zeros((len(ts), N), device=device, dtype=torch.float32)
    
    diag_k0 = (z * Tm2).mean(dim=1)
    diag_k1 = (z * Tm1).mean(dim=1)
    
    acc += torch.outer(C_pt[:, 0], diag_k0) + torch.outer(C_pt[:, 1], diag_k1)
    
    for k in range(2, K + 1):
        Tk = 2.0 * torch.sparse.mm(Lt_pt, Tm1) - Tm2
        diag_k = (z * Tk).mean(dim=1)
        acc += torch.outer(C_pt[:, k], diag_k)
        Tm2, Tm1 = Tm1, Tk
        
    # Clamp strictly positive before returning (physical heat cannot be <= 0)
    return acc.cpu().numpy()

def get_scale_bank(L_c):
    """
    Computes the 16-step logarithmic diffusion scale bank (Proposal Step 4).
    """
    lam_max = float(splinalg.eigsh(L_c, k=1, which='LM', return_eigenvectors=False)[0])
    lam_2 = 1e-4  # Safe lower spectral bound for connected PPI networks
        
    t_min = 4.0 * np.log(10.0) / lam_max
    t_max = 4.0 * np.log(10.0) / lam_2
    
    ts = np.logspace(np.log10(t_min), np.log10(t_max), 16)
    return ts, lam_max

def get_lazy_walk_profile(P_t, k_steps=(2, 4, 8, 16)):
    """
    Computes Lazy Random-Walk Structural Encoding (RWSE) diagonals.
    """
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