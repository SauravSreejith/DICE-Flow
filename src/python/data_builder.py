#!/usr/bin/env python3
"""
DiCE-Flow Data Builder
======================
Parses TCGA expression matrices, Pearson-distance edge weights, and merged
network directionality into PyTorch Geometric ``Data`` objects for the
Edge-Augmented Graph Attention Network (EA-GAT).

Critical: gene-level intersection produces a universal vocabulary so that
normal and tumor embeddings share identical node indices for Embedding Drift
(Δh) comparison downstream.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
)
logger = logging.getLogger(__name__)

# Default paths relative to the project root
_DEFAULT_ROOT = Path(__file__).resolve().parents[2]


class DiCEDataBuilder:
    """
    Build PyG graphs for normal and tumor tissue conditions.

    Parameters
    ----------
    data_root :
        Project root directory containing ``data/raw`` and ``data/processed``.
        Defaults to two levels above this file (repository root).
    expression_normal :
        Override path for normal expression TSV (genes × samples).
    expression_tumor :
        Override path for tumor expression TSV (genes × samples).
    network_normal :
        Override path for normal weighted edges
        ``[source, target, weight, weight_int]``.
    network_tumor :
        Override path for tumor weighted edges.
    merged_network :
        Override path for merged topology with directionality
        ``[source, target, source_db, directed]``.
    """

    def __init__(
        self,
        data_root: Optional[Union[str, Path]] = None,
        expression_normal: Optional[Union[str, Path]] = None,
        expression_tumor: Optional[Union[str, Path]] = None,
        network_normal: Optional[Union[str, Path]] = None,
        network_tumor: Optional[Union[str, Path]] = None,
        merged_network: Optional[Union[str, Path]] = None,
    ) -> None:
        self.data_root = Path(data_root) if data_root is not None else _DEFAULT_ROOT

        self.expression_normal_path = Path(
            expression_normal
            or self.data_root / "data/raw/expression/expression_normal.txt"
        )
        self.expression_tumor_path = Path(
            expression_tumor
            or self.data_root / "data/raw/expression/expression_tumor.txt"
        )
        self.network_normal_path = Path(
            network_normal
            or self.data_root / "data/processed/weighted/network_normal.txt"
        )
        self.network_tumor_path = Path(
            network_tumor
            or self.data_root / "data/processed/weighted/network_tumor.txt"
        )
        self.merged_network_path = Path(
            merged_network
            or self.data_root / "data/processed/networks/merged_network.txt"
        )

        self.gene_to_idx: Dict[str, int] = {}
        self.idx_to_gene: Dict[int, str] = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def build_graphs(self) -> Dict[str, Union[Data, Dict[str, int]]]:
        """
        Construct normal and tumor PyG ``Data`` objects with a shared gene map.

        Returns
        -------
        dict
            ``{"normal": Data, "tumor": Data, "gene_map": gene_to_idx}``

        Raises
        ------
        FileNotFoundError
            If any required input file is missing.
        ValueError
            If the gene intersection is empty or tensors cannot be formed.
        """
        self._validate_paths()

        logger.info("Loading expression matrices...")
        expr_normal = self._load_expression(self.expression_normal_path)
        expr_tumor = self._load_expression(self.expression_tumor_path)

        logger.info("Loading weighted edge lists...")
        edges_normal = self._load_weighted_edges(self.network_normal_path)
        edges_tumor = self._load_weighted_edges(self.network_tumor_path)

        logger.info("Loading merged network directionality...")
        directionality = self._load_directionality(self.merged_network_path)

        # Universal vocabulary: genes in BOTH expression matrices ∩ edge lists
        self.gene_to_idx = self._build_gene_vocabulary(
            expr_normal, expr_tumor, edges_normal, edges_tumor
        )
        self.idx_to_gene = {idx: gene for gene, idx in self.gene_to_idx.items()}
        n_genes = len(self.gene_to_idx)
        logger.info("Universal gene vocabulary size: %d", n_genes)

        if n_genes == 0:
            raise ValueError(
                "Gene intersection is empty. Check that expression matrices "
                "and edge lists share common gene symbols."
            )

        normal_data = self._build_single_graph(
            condition="normal",
            expression=expr_normal,
            edges=edges_normal,
            directionality=directionality,
        )
        tumor_data = self._build_single_graph(
            condition="tumor",
            expression=expr_tumor,
            edges=edges_tumor,
            directionality=directionality,
        )

        logger.info(
            "Built graphs — normal: %d nodes / %d edges | tumor: %d nodes / %d edges",
            normal_data.num_nodes,
            normal_data.num_edges,
            tumor_data.num_nodes,
            tumor_data.num_edges,
        )

        return {
            "normal": normal_data,
            "tumor": tumor_data,
            "gene_map": self.gene_to_idx,
        }

    # ------------------------------------------------------------------
    # Path / I/O helpers
    # ------------------------------------------------------------------

    def _validate_paths(self) -> None:
        """Raise ``FileNotFoundError`` for any missing required input."""
        required = {
            "expression_normal": self.expression_normal_path,
            "expression_tumor": self.expression_tumor_path,
            "network_normal": self.network_normal_path,
            "network_tumor": self.network_tumor_path,
            "merged_network": self.merged_network_path,
        }
        missing = [name for name, path in required.items() if not path.is_file()]
        if missing:
            details = "\n".join(
                f"  - {name}: {required[name]}" for name in missing
            )
            raise FileNotFoundError(
                f"Missing required data file(s):\n{details}\n"
                f"Expected project root: {self.data_root}"
            )

    @staticmethod
    def _load_expression(path: Path) -> pd.DataFrame:
        """
        Load a genes × samples expression matrix (TSV, gene symbols as index).

        The first column is treated as the gene identifier regardless of name
        (often labeled ``sample`` in TCGA exports).
        """
        logger.info("Reading expression matrix: %s", path)
        df = pd.read_csv(path, sep="\t", index_col=0)
        # Coerce numeric; drop all-NaN rows that cannot form features
        df = df.apply(pd.to_numeric, errors="coerce")
        df = df.dropna(how="all")
        # Normalize gene symbols for robust intersection
        df.index = df.index.astype(str).str.strip()
        # Collapse duplicate gene symbols by mean
        if df.index.duplicated().any():
            n_dup = int(df.index.duplicated().sum())
            logger.warning(
                "Collapsing %d duplicate gene symbols in %s via mean", n_dup, path.name
            )
            df = df.groupby(df.index).mean()
        logger.info(
            "  → %d genes × %d samples", df.shape[0], df.shape[1]
        )
        return df

    @staticmethod
    def _load_weighted_edges(path: Path) -> pd.DataFrame:
        """Load weighted network TSV with columns source, target, weight, weight_int."""
        logger.info("Reading weighted edges: %s", path)
        df = pd.read_csv(
            path,
            sep="\t",
            usecols=["source", "target", "weight"],
            dtype={"source": str, "target": str, "weight": np.float64},
        )
        df["source"] = df["source"].str.strip()
        df["target"] = df["target"].str.strip()
        df = df.dropna(subset=["source", "target", "weight"])
        logger.info("  → %d edges", len(df))
        return df

    @staticmethod
    def _load_directionality(path: Path) -> pd.DataFrame:
        """
        Load merged network directionality.

        Returns a DataFrame with unique (source, target) pairs and a float
        ``directed`` column (1.0 / 0.0).
        """
        logger.info("Reading merged network (directionality): %s", path)
        df = pd.read_csv(
            path,
            sep="\t",
            usecols=["source", "target", "directed"],
            dtype={"source": str, "target": str},
        )
        df["source"] = df["source"].str.strip()
        df["target"] = df["target"].str.strip()
        # Boolean or string True/False → float
        df["directed"] = df["directed"].map(
            lambda v: 1.0 if str(v).strip().lower() in {"true", "1", "yes"} else 0.0
        )
        # Prefer directed=True when duplicates exist across source DBs
        df = (
            df.groupby(["source", "target"], as_index=False)["directed"]
            .max()
        )
        logger.info("  → %d unique directed pairs", len(df))
        return df

    # ------------------------------------------------------------------
    # Vocabulary & feature engineering
    # ------------------------------------------------------------------

    @staticmethod
    def _edge_gene_set(edges: pd.DataFrame) -> Set[str]:
        """Return the set of gene symbols appearing as source or target."""
        return set(edges["source"]).union(set(edges["target"]))

    def _build_gene_vocabulary(
        self,
        expr_normal: pd.DataFrame,
        expr_tumor: pd.DataFrame,
        edges_normal: pd.DataFrame,
        edges_tumor: pd.DataFrame,
    ) -> Dict[str, int]:
        """
        Intersection of genes present in both expression matrices and the
        combined edge lists. Sorted alphabetically for deterministic indices.
        """
        expr_genes = set(expr_normal.index).intersection(set(expr_tumor.index))
        edge_genes = self._edge_gene_set(edges_normal).union(
            self._edge_gene_set(edges_tumor)
        )
        common = sorted(expr_genes.intersection(edge_genes))
        logger.info(
            "Vocabulary intersection — expr∩both: %d | edge genes: %d | final: %d",
            len(expr_genes),
            len(edge_genes),
            len(common),
        )
        return {gene: idx for idx, gene in enumerate(common)}

    def _compute_node_features(
        self, expression: pd.DataFrame
    ) -> torch.Tensor:
        """
        Node features ``x`` of shape ``[N, 2]``: row-wise mean and variance.

        Genes missing from the expression matrix (should not occur after
        intersection) receive zero features.
        """
        n = len(self.gene_to_idx)
        features = np.zeros((n, 2), dtype=np.float32)

        # Align expression rows to vocabulary order efficiently
        present = [g for g in self.gene_to_idx if g in expression.index]
        if not present:
            raise ValueError("No vocabulary genes found in expression matrix.")

        sub = expression.loc[present]
        means = sub.mean(axis=1).to_numpy(dtype=np.float32)
        # ddof=0 population variance for stability with small sample counts
        variances = sub.var(axis=1, ddof=0).to_numpy(dtype=np.float32)
        # Replace NaN (e.g. single-sample or all-missing rows)
        means = np.nan_to_num(means, nan=0.0)
        variances = np.nan_to_num(variances, nan=0.0)

        for gene, mean_val, var_val in zip(present, means, variances):
            idx = self.gene_to_idx[gene]
            features[idx, 0] = mean_val
            features[idx, 1] = var_val

        return torch.from_numpy(features)

    def _build_edge_tensors(
        self,
        edges: pd.DataFrame,
        directionality: pd.DataFrame,
    ) -> tuple[torch.LongTensor, torch.Tensor]:
        """
        Filter edges to the gene vocabulary, merge directionality, and return
        ``edge_index`` ``[2, E]`` and ``edge_attr`` ``[E, 3]``.

        edge_attr columns
        -----------------
        0 : weight (Pearson distance 1 − |r|)
        1 : directed (1.0 / 0.0)
        2 : placeholder zeros for Structural Knockout Signatures (SKS)
        """
        gene_set = self.gene_to_idx
        mask = edges["source"].isin(gene_set) & edges["target"].isin(gene_set)
        filtered = edges.loc[mask].copy()
        if filtered.empty:
            raise ValueError(
                "No edges remain after filtering to the gene vocabulary."
            )

        # Merge directionality on (source, target); default undirected = 0.0
        merged = filtered.merge(
            directionality,
            on=["source", "target"],
            how="left",
        )
        merged["directed"] = merged["directed"].fillna(0.0).astype(np.float32)
        merged["weight"] = merged["weight"].astype(np.float32)

        src_idx = merged["source"].map(gene_set).to_numpy(dtype=np.int64)
        tgt_idx = merged["target"].map(gene_set).to_numpy(dtype=np.int64)

        edge_index = torch.from_numpy(
            np.vstack([src_idx, tgt_idx])
        ).contiguous().long()

        # Col 2: zeros — placeholder for Structural Knockout Signatures (SKS)
        sks_placeholder = np.zeros(len(merged), dtype=np.float32)
        edge_attr = torch.from_numpy(
            np.column_stack(
                [
                    merged["weight"].to_numpy(dtype=np.float32),
                    merged["directed"].to_numpy(dtype=np.float32),
                    sks_placeholder,
                ]
            )
        ).contiguous()

        logger.info(
            "  edge_index %s | edge_attr %s | directed fraction %.3f",
            tuple(edge_index.shape),
            tuple(edge_attr.shape),
            float(merged["directed"].mean()),
        )
        return edge_index, edge_attr

    def _build_single_graph(
        self,
        condition: str,
        expression: pd.DataFrame,
        edges: pd.DataFrame,
        directionality: pd.DataFrame,
    ) -> Data:
        """Assemble a single ``torch_geometric.data.Data`` for one condition."""
        logger.info("Building %s graph...", condition)
        x = self._compute_node_features(expression)
        edge_index, edge_attr = self._build_edge_tensors(edges, directionality)

        data = Data(
            x=x,
            edge_index=edge_index,
            edge_attr=edge_attr,
        )
        data.num_nodes = x.size(0)
        data.condition = condition  # metadata for debugging / logging
        return data


def build_graphs(
    data_root: Optional[Union[str, Path]] = None,
) -> Dict[str, Union[Data, Dict[str, int]]]:
    """
    Convenience wrapper around :class:`DiCEDataBuilder`.

    Parameters
    ----------
    data_root :
        Optional project root override.

    Returns
    -------
    dict
        ``{"normal": Data, "tumor": Data, "gene_map": gene_to_idx}``
    """
    builder = DiCEDataBuilder(data_root=data_root)
    return builder.build_graphs()


if __name__ == "__main__":
    result = build_graphs()
    normal: Data = result["normal"]  # type: ignore[assignment]
    tumor: Data = result["tumor"]  # type: ignore[assignment]
    gene_map: Dict[str, int] = result["gene_map"]  # type: ignore[assignment]
    print(f"gene_map size : {len(gene_map)}")
    print(f"normal        : x={tuple(normal.x.shape)} edges={normal.num_edges}")
    print(f"tumor         : x={tuple(tumor.x.shape)} edges={tumor.num_edges}")
    print(f"edge_attr dim : {normal.edge_attr.shape[1]}")
