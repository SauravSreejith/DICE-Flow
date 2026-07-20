#!/usr/bin/env python3
"""
DiCE-Flow EA-GAT Architecture
=============================
Edge-Augmented Graph Attention Network using GATv2Conv.

Design rationale
----------------
* **GATv2** (Brody et al., 2022) provides *dynamic* attention, relaxing the
  static-attention limitation of classic GAT and improving performance under
  heterophily common in biological interaction networks.
* **edge_dim** injects edge attributes (Pearson distance weight, directionality,
  and a reserved SKS channel) into the attention score, enabling anisotropic
  message routing that mitigates oversquashing on long biological pathways.

This module deliberately returns *raw* node embeddings only — no classification
head — so a downstream contrastive Embedding Drift (Δh) loss can operate on
continuous latent representations.
"""

from __future__ import annotations

from typing import Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATv2Conv


class DiCE_EAGAT(nn.Module):
    """
    Two-layer Edge-Augmented Graph Attention Network (EA-GAT).

    Architecture
    ------------
    Layer 1
        ``GATv2Conv(node_dim → hidden_dim, heads=4, concat=True, edge_dim)``
        followed by ELU + Dropout(0.4). Output width = ``hidden_dim * 4``.
    Layer 2
        ``GATv2Conv(hidden_dim*4 → embedding_dim, heads=1, concat=False, edge_dim)``
        returns raw embeddings ``[N, embedding_dim]``.

    Parameters
    ----------
    node_dim :
        Dimensionality of input node features (default 2: mean & variance).
    edge_attr_dim :
        Dimensionality of edge attributes (default 3: weight, directed, SKS).
    hidden_dim :
        Per-head hidden size of the first GATv2 layer (default 32).
    embedding_dim :
        Final embedding dimensionality (default 16).
    dropout :
        Dropout probability applied after Layer 1 (default 0.4).
    heads :
        Number of attention heads in Layer 1 (default 4).
    """

    def __init__(
        self,
        node_dim: int = 2,
        edge_attr_dim: int = 3,
        hidden_dim: int = 32,
        embedding_dim: int = 16,
        dropout: float = 0.4,
        heads: int = 4,
    ) -> None:
        super().__init__()

        if node_dim < 1:
            raise ValueError(f"node_dim must be ≥ 1, got {node_dim}")
        if edge_attr_dim < 1:
            raise ValueError(f"edge_attr_dim must be ≥ 1, got {edge_attr_dim}")
        if hidden_dim < 1 or embedding_dim < 1:
            raise ValueError("hidden_dim and embedding_dim must be ≥ 1")
        if not 0.0 <= dropout < 1.0:
            raise ValueError(f"dropout must be in [0, 1), got {dropout}")
        if heads < 1:
            raise ValueError(f"heads must be ≥ 1, got {heads}")

        self.node_dim = node_dim
        self.edge_attr_dim = edge_attr_dim
        self.hidden_dim = hidden_dim
        self.embedding_dim = embedding_dim
        self.dropout_p = dropout
        self.heads = heads

        # Layer 1: multi-head GATv2 with edge attributes; concat → hidden_dim * heads
        self.conv1 = GATv2Conv(
            in_channels=node_dim,
            out_channels=hidden_dim,
            heads=heads,
            edge_dim=edge_attr_dim,
            concat=True,
            dropout=dropout,
        )

        # Layer 2: single-head projection to embedding space
        self.conv2 = GATv2Conv(
            in_channels=hidden_dim * heads,
            out_channels=embedding_dim,
            heads=1,
            edge_dim=edge_attr_dim,
            concat=False,
            dropout=dropout,
        )

        self.dropout = nn.Dropout(p=dropout)

    def forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute node embeddings for a single graph.

        Parameters
        ----------
        x :
            Node feature matrix of shape ``[N, node_dim]``.
        edge_index :
            COO edge index of shape ``[2, E]`` (``torch.long``).
        edge_attr :
            Edge attribute matrix of shape ``[E, edge_attr_dim]``.

        Returns
        -------
        torch.Tensor
            Raw embeddings of shape ``[N, embedding_dim]``. No Softmax / MLP.
        """
        if x.dim() != 2:
            raise ValueError(f"Expected x of rank 2 [N, F], got shape {tuple(x.shape)}")
        if edge_index.dim() != 2 or edge_index.size(0) != 2:
            raise ValueError(
                f"Expected edge_index of shape [2, E], got {tuple(edge_index.shape)}"
            )
        if edge_attr is None:
            raise ValueError(
                "edge_attr is required for EA-GAT (edge_dim is configured). "
                "Pass a tensor of shape [E, edge_attr_dim]."
            )
        if edge_attr.size(0) != edge_index.size(1):
            raise ValueError(
                f"edge_attr rows ({edge_attr.size(0)}) must match "
                f"edge_index columns ({edge_index.size(1)})"
            )
        if edge_attr.size(1) != self.edge_attr_dim:
            raise ValueError(
                f"edge_attr last dim must be {self.edge_attr_dim}, "
                f"got {edge_attr.size(1)}"
            )

        # --- Layer 1: multi-head anisotropic attention ---
        h = self.conv1(x, edge_index, edge_attr=edge_attr)
        h = F.elu(h)
        h = self.dropout(h)

        # --- Layer 2: embedding projection (raw, no activation head) ---
        embeddings = self.conv2(h, edge_index, edge_attr=edge_attr)

        return embeddings

    def forward_with_attention(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Forward pass that also returns Layer-2 attention coefficients (α_ij).

        Returns
        -------
        embeddings : Tensor [N, embedding_dim]
        attn_edge_index : LongTensor [2, E']
            Edge index corresponding to ``alpha`` (may include self-loops
            inserted by GATv2Conv).
        alpha : Tensor [E']
            Per-edge attention weights (mean over heads if multi-head).
        """
        h = self.conv1(x, edge_index, edge_attr=edge_attr)
        h = F.elu(h)
        h = self.dropout(h)

        embeddings, (attn_edge_index, alpha) = self.conv2(
            h,
            edge_index,
            edge_attr=edge_attr,
            return_attention_weights=True,
        )
        # alpha: [E, heads] or [E, 1] → collapse to [E]
        if alpha.dim() > 1:
            alpha = alpha.mean(dim=-1)
        else:
            alpha = alpha.view(-1)

        return embeddings, attn_edge_index, alpha

    def extra_repr(self) -> str:
        return (
            f"node_dim={self.node_dim}, edge_attr_dim={self.edge_attr_dim}, "
            f"hidden_dim={self.hidden_dim}, embedding_dim={self.embedding_dim}, "
            f"heads={self.heads}, dropout={self.dropout_p}"
        )


def count_parameters(model: nn.Module) -> int:
    """Return the number of trainable parameters in ``model``."""
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


if __name__ == "__main__":
    # Smoke test with synthetic tensors (no data files required)
    torch.manual_seed(0)
    n_nodes, n_edges = 32, 64
    model = DiCE_EAGAT(node_dim=2, edge_attr_dim=3, hidden_dim=32, embedding_dim=16)
    x = torch.randn(n_nodes, 2)
    edge_index = torch.randint(0, n_nodes, (2, n_edges), dtype=torch.long)
    edge_attr = torch.rand(n_edges, 3)
    out = model(x, edge_index, edge_attr)
    print(model)
    print(f"parameters : {count_parameters(model):,}")
    print(f"embeddings : {tuple(out.shape)}")  # expected (32, 16)
