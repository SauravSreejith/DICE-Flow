#!/usr/bin/env python3
"""
DiCE-Flow EA-GAT Pipeline Runner
================================
Loads normal / tumor PyG graphs, runs EA-GAT, computes Embedding Drift (Δh),
exports rankings + edge attention for the FastAPI / React UI, and prints
embedding shapes for verification.

Paths are hardcoded relative to the project root — no CLI argument parsing.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import torch

# Ensure ``src/python`` is importable when executed as a script
_THIS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _THIS_DIR.parents[1]
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from data_builder import DiCEDataBuilder  # noqa: E402
from models import DiCE_EAGAT, count_parameters  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
)
logger = logging.getLogger("pipeline_runner")

# ---------------------------------------------------------------------------
# Hardcoded configuration (project-root relative)
# ---------------------------------------------------------------------------
PROJECT_ROOT = _PROJECT_ROOT
RESULTS_DIR = PROJECT_ROOT / "data" / "results" / "eagat"
NODE_DIM = 2
EDGE_ATTR_DIM = 3
HIDDEN_DIM = 32
EMBEDDING_DIM = 16
# Export attention for edges among the top-K drift genes (UI skeleton graph)
ATTENTION_TOP_K = 200


def select_device() -> torch.device:
    """Return CUDA device when available, otherwise CPU."""
    if torch.cuda.is_available():
        device = torch.device("cuda")
        logger.info("Using device: cuda (%s)", torch.cuda.get_device_name(0))
    else:
        device = torch.device("cpu")
        logger.info("Using device: cpu")
    return device


def _invert_gene_map(gene_map: Dict[str, int]) -> Dict[int, str]:
    return {idx: gene for gene, idx in gene_map.items()}


def compute_embedding_drift(
    h_normal: torch.Tensor,
    h_tumor: torch.Tensor,
) -> torch.Tensor:
    """Per-gene L2 Embedding Drift: Δh_i = ||h_tumor[i] − h_normal[i]||_2."""
    return torch.norm(h_tumor - h_normal, p=2, dim=1)


def build_gene_rankings(
    gene_map: Dict[str, int],
    drift: torch.Tensor,
) -> pd.DataFrame:
    """
    Rank genes by Embedding Drift and derive UI metrics.

    Columns
    -------
    gene, rank, embedding_drift, egat_prob, sks_disruption
    """
    idx_to_gene = _invert_gene_map(gene_map)
    drift_np = drift.detach().cpu().numpy().astype(np.float64)

    # Min-max normalize drift → egat_prob in [0, 1]
    d_min, d_max = float(drift_np.min()), float(drift_np.max())
    if d_max > d_min:
        egat_prob = (drift_np - d_min) / (d_max - d_min)
    else:
        egat_prob = np.zeros_like(drift_np)

    # SKS channel not yet available — placeholder zeros
    # [WAITING ON UPSTREAM] Structural Knockout Signature arrays
    sks = np.zeros_like(drift_np)

    order = np.argsort(-drift_np)  # descending drift
    rows = []
    for rank, idx in enumerate(order, start=1):
        rows.append(
            {
                "gene": idx_to_gene[int(idx)],
                "rank": rank,
                "embedding_drift": float(drift_np[idx]),
                "egat_prob": float(egat_prob[idx]),
                "sks_disruption": float(sks[idx]),
            }
        )
    return pd.DataFrame(rows)


def export_edge_attention(
    model: DiCE_EAGAT,
    data,
    gene_map: Dict[str, int],
    top_genes: list[str],
    out_path: Path,
) -> int:
    """
    Run Layer-2 attention on ``data`` and export α_ij for edges whose
    endpoints are both in ``top_genes``.

    Returns number of attention edges written.
    """
    idx_to_gene = _invert_gene_map(gene_map)
    top_set = set(top_genes)
    top_idx = {gene_map[g] for g in top_genes if g in gene_map}

    model.eval()
    with torch.no_grad():
        _, attn_edge_index, alpha = model.forward_with_attention(
            data.x, data.edge_index, data.edge_attr
        )

    src = attn_edge_index[0].detach().cpu().numpy()
    tgt = attn_edge_index[1].detach().cpu().numpy()
    att = alpha.detach().cpu().numpy()

    rows = []
    for s, t, a in zip(src, tgt, att):
        s, t = int(s), int(t)
        if s == t:
            continue  # skip self-loops added by GATv2
        if s not in top_idx or t not in top_idx:
            continue
        g_s = idx_to_gene.get(s)
        g_t = idx_to_gene.get(t)
        if g_s is None or g_t is None:
            continue
        if g_s not in top_set or g_t not in top_set:
            continue
        rows.append({"source": g_s, "target": g_t, "attention": float(a)})

    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
    logger.info("Wrote %d attention edges → %s", len(rows), out_path)
    return len(rows)


def export_artifacts(
    model: DiCE_EAGAT,
    gene_map: Dict[str, int],
    rankings: pd.DataFrame,
    tumor_data,
    device: torch.device,
) -> None:
    """Persist rankings, model weights, gene map, and tumor-edge attention."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    rankings_path = RESULTS_DIR / "gene_rankings.tsv"
    rankings.to_csv(rankings_path, sep="\t", index=False)
    logger.info("Wrote gene rankings → %s (%d genes)", rankings_path, len(rankings))

    top_path = RESULTS_DIR / "gene_rankings_top500.tsv"
    rankings.head(500).to_csv(top_path, sep="\t", index=False)

    model_path = RESULTS_DIR / "model.pt"
    torch.save(
        {
            "state_dict": model.state_dict(),
            "config": {
                "node_dim": NODE_DIM,
                "edge_attr_dim": EDGE_ATTR_DIM,
                "hidden_dim": HIDDEN_DIM,
                "embedding_dim": EMBEDDING_DIM,
            },
        },
        model_path,
    )
    logger.info("Wrote model checkpoint → %s", model_path)

    gene_map_path = RESULTS_DIR / "gene_map.json"
    with open(gene_map_path, "w", encoding="utf-8") as fh:
        json.dump(gene_map, fh)
    logger.info("Wrote gene map → %s", gene_map_path)

    top_genes = rankings.head(ATTENTION_TOP_K)["gene"].tolist()
    attn_path = RESULTS_DIR / "edge_attention_tumor.tsv"
    export_edge_attention(model, tumor_data, gene_map, top_genes, attn_path)

    meta = {
        "n_genes": len(gene_map),
        "embedding_dim": EMBEDDING_DIM,
        "attention_top_k": ATTENTION_TOP_K,
        "device": str(device),
        "top_gene": rankings.iloc[0]["gene"] if len(rankings) else None,
        "max_drift": float(rankings.iloc[0]["embedding_drift"]) if len(rankings) else 0.0,
    }
    with open(RESULTS_DIR / "meta.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)


def run_pipeline() -> Tuple[torch.Tensor, torch.Tensor, pd.DataFrame]:
    """
    End-to-end EA-GAT forward pass, drift ranking, and artifact export.

    Returns
    -------
    h_normal, h_tumor, rankings
    """
    device = select_device()

    # ------------------------------------------------------------------
    # 1. Data loading
    # ------------------------------------------------------------------
    logger.info("Instantiating DiCEDataBuilder (root=%s)", PROJECT_ROOT)
    try:
        builder = DiCEDataBuilder(data_root=PROJECT_ROOT)
        graphs = builder.build_graphs()
    except FileNotFoundError as exc:
        logger.error("Data loading failed — missing files:\n%s", exc)
        sys.exit(1)
    except ValueError as exc:
        logger.error("Data loading failed — invalid content:\n%s", exc)
        sys.exit(1)
    except Exception as exc:  # noqa: BLE001
        logger.exception("Unexpected error while building graphs: %s", exc)
        sys.exit(1)

    normal_data = graphs["normal"].to(device)
    tumor_data = graphs["tumor"].to(device)
    gene_map = graphs["gene_map"]

    logger.info(
        "Graphs ready — genes=%d | normal edges=%d | tumor edges=%d",
        len(gene_map),
        normal_data.num_edges,
        tumor_data.num_edges,
    )

    # ------------------------------------------------------------------
    # 2. Model instantiation
    # ------------------------------------------------------------------
    model = DiCE_EAGAT(
        node_dim=NODE_DIM,
        edge_attr_dim=EDGE_ATTR_DIM,
        hidden_dim=HIDDEN_DIM,
        embedding_dim=EMBEDDING_DIM,
    ).to(device)
    model.eval()

    logger.info("Model: %s", model)
    logger.info("Trainable parameters: %s", f"{count_parameters(model):,}")

    # ------------------------------------------------------------------
    # 3. Forward pass (embeddings)
    # ------------------------------------------------------------------
    try:
        with torch.no_grad():
            h_normal = model(
                normal_data.x,
                normal_data.edge_index,
                normal_data.edge_attr,
            )
            h_tumor = model(
                tumor_data.x,
                tumor_data.edge_index,
                tumor_data.edge_attr,
            )
    except RuntimeError as exc:
        logger.error("Forward pass failed: %s", exc)
        sys.exit(1)

    # ------------------------------------------------------------------
    # 4. Embedding Drift + export for UI
    # ------------------------------------------------------------------
    drift = compute_embedding_drift(h_normal, h_tumor)
    rankings = build_gene_rankings(gene_map, drift)
    export_artifacts(model, gene_map, rankings, tumor_data, device)

    # ------------------------------------------------------------------
    # 5. Verification output
    # ------------------------------------------------------------------
    print("=" * 60)
    print("EA-GAT forward pass complete")
    print(f"  device          : {device}")
    print(f"  gene vocabulary : {len(gene_map)}")
    print(f"  h_normal shape  : {tuple(h_normal.shape)}")
    print(f"  h_tumor  shape  : {tuple(h_tumor.shape)}")
    print(f"  h_normal dtype  : {h_normal.dtype}")
    print(f"  h_tumor  dtype  : {h_tumor.dtype}")
    print(f"  results dir     : {RESULTS_DIR}")
    print("  top 5 by Δh     :")
    for _, row in rankings.head(5).iterrows():
        print(
            f"    #{int(row['rank']):3d}  {row['gene']:<12}  "
            f"Δh={row['embedding_drift']:.4f}  p={row['egat_prob']:.3f}"
        )
    print("=" * 60)

    # ==================================================================
    # [HANDOFF TO DOWNSTREAM]
    # h_normal and h_tumor are ready for contrastive loss and Embedding
    # Drift (Δh) calculation. Rankings + attention exported under
    # data/results/eagat/ for FastAPI and the React UI.
    # ==================================================================

    # ==================================================================
    # [WAITING ON UPSTREAM]
    # Awaiting Structural Knockout Signature (SKS) numpy arrays to inject
    # into edge_attr Col 2. Until then sks_disruption remains 0.0.
    # ==================================================================

    return h_normal, h_tumor, rankings


if __name__ == "__main__":
    run_pipeline()
