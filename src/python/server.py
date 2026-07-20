"""
DiCE-Duan FastAPI backend.

Serves classic NetworkX graphs (hub rewiring, path tracer, knockout) and
EA-GAT deep-learning artifacts (Embedding Drift rankings, GAT attention,
edge directionality) for the React UI.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import networkx as nx
import pandas as pd
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware

# --- CONFIGURATION ---
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DiCE-API")

app = FastAPI(title="DiCE-Duan API", version="3.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- GLOBAL STATE ---
DATA_STORE: Dict[str, Any] = {
    "G_normal": None,
    "G_tumor": None,
    "knockout_df": None,
    "top_genes": [],
    "log2fc_map": {},
    "essential_genes": set(),
    # EA-GAT
    "eagat_rankings": None,  # DataFrame
    "eagat_by_gene": {},  # gene -> row dict
    "edge_attention": {},  # (src, tgt) frozenset or ordered pair -> float
    "edge_directed": {},  # (src, tgt) -> bool
    "eagat_meta": {},
}

BASE_DIR = Path(__file__).parent.parent.parent
PROCESSED_DIR = BASE_DIR / "data/processed/weighted"
NETWORKS_DIR = BASE_DIR / "data/processed/networks"
RESULTS_DIR = BASE_DIR / "data/results"
EAGAT_DIR = RESULTS_DIR / "eagat"
TABLES_DIR = RESULTS_DIR / "validation/tables"


def load_graph(file_path: Path) -> nx.Graph:
    if not file_path.exists():
        logger.error("Graph file not found: %s", file_path)
        return nx.Graph()
    df = pd.read_csv(file_path, sep="\t")
    G = nx.from_pandas_edgelist(df, "source", "target", edge_attr="weight")
    return G


def load_directionality(file_path: Path) -> Dict[Tuple[str, str], bool]:
    """Map (source, target) → directed flag from merged_network.txt."""
    directed_map: Dict[Tuple[str, str], bool] = {}
    if not file_path.exists():
        logger.warning("Merged network not found (no directionality): %s", file_path)
        return directed_map

    logger.info("Loading edge directionality from %s ...", file_path.name)
    # Stream in chunks to avoid holding the full 1.9M-row frame twice
    usecols = ["source", "target", "directed"]
    for chunk in pd.read_csv(file_path, sep="\t", usecols=usecols, chunksize=200_000):
        chunk["source"] = chunk["source"].astype(str).str.strip()
        chunk["target"] = chunk["target"].astype(str).str.strip()
        chunk["directed"] = chunk["directed"].map(
            lambda v: str(v).strip().lower() in {"true", "1", "yes"}
        )
        for src, tgt, d in chunk[["source", "target", "directed"]].itertuples(
            index=False, name=None
        ):
            key = (src, tgt)
            # Prefer True if any source DB marks the edge directed
            directed_map[key] = directed_map.get(key, False) or bool(d)
    logger.info("  → %d directed-pair records", len(directed_map))
    return directed_map


def load_edge_attention(file_path: Path) -> Dict[Tuple[str, str], float]:
    """Map (source, target) → GAT attention α_ij (and undirected reverse)."""
    attn: Dict[Tuple[str, str], float] = {}
    if not file_path.exists():
        logger.warning("EA-GAT attention file not found: %s", file_path)
        return attn

    df = pd.read_csv(file_path, sep="\t")
    for src, tgt, a in df[["source", "target", "attention"]].itertuples(
        index=False, name=None
    ):
        src, tgt = str(src), str(tgt)
        val = float(a)
        attn[(src, tgt)] = val
        # NetworkX graphs are undirected — store reverse for lookup convenience
        attn.setdefault((tgt, src), val)
    logger.info("Loaded %d EA-GAT attention edges", len(df))
    return attn


def load_eagat_rankings(file_path: Path) -> Optional[pd.DataFrame]:
    if not file_path.exists():
        logger.warning("EA-GAT rankings not found: %s", file_path)
        return None
    df = pd.read_csv(file_path, sep="\t")
    required = {"gene", "rank", "embedding_drift", "egat_prob", "sks_disruption"}
    missing = required - set(df.columns)
    if missing:
        logger.error("EA-GAT rankings missing columns: %s", missing)
        return None
    logger.info("Loaded EA-GAT rankings for %d genes", len(df))
    return df


@app.on_event("startup")
async def startup_event():
    logger.info("Loading DiCE-Duan Data...")

    # 1. Load Graphs
    DATA_STORE["G_normal"] = load_graph(PROCESSED_DIR / "network_normal.txt")
    DATA_STORE["G_tumor"] = load_graph(PROCESSED_DIR / "network_tumor.txt")

    # 2. Load Knockout Results
    ko_path = RESULTS_DIR / "knockout_analysis.txt"
    if ko_path.exists():
        DATA_STORE["knockout_df"] = pd.read_csv(ko_path, sep="\t")

    # 3. Load Top Genes (prefer EA-GAT rankings when available)
    eagat_top = EAGAT_DIR / "gene_rankings_top500.tsv"
    eagat_all = EAGAT_DIR / "gene_rankings.tsv"
    rankings = load_eagat_rankings(
        eagat_top if eagat_top.exists() else eagat_all
    )
    if rankings is not None:
        DATA_STORE["eagat_rankings"] = rankings
        DATA_STORE["eagat_by_gene"] = {
            str(row["gene"]).upper(): {
                "gene": str(row["gene"]),
                "rank": int(row["rank"]),
                "embeddingDrift": float(row["embedding_drift"]),
                "egatProb": float(row["egat_prob"]),
                "sksDisruption": float(row["sks_disruption"]),
            }
            for _, row in rankings.iterrows()
        }
        DATA_STORE["top_genes"] = rankings["gene"].astype(str).tolist()
        logger.info("Using EA-GAT top genes for network skeleton")
    else:
        dice_path = RESULTS_DIR / "dice_genes_top500.txt"
        if dice_path.exists():
            df_top = pd.read_csv(dice_path, sep="\t")
            DATA_STORE["top_genes"] = df_top["gene"].tolist()
            logger.info("Fallback: classic DiCE top genes")

    # 4. DEA annotations
    dea_path = RESULTS_DIR / "dea_results_table.txt"
    if dea_path.exists():
        dea_df = pd.read_csv(dea_path, sep="\t")
        DATA_STORE["log2fc_map"] = dict(zip(dea_df["gene"], dea_df["log2FC"]))

    # 5. Essentiality
    essential_path = TABLES_DIR / "09_fitness_overlap_genes.csv"
    if essential_path.exists():
        ess_df = pd.read_csv(essential_path)
        DATA_STORE["essential_genes"] = set(ess_df["gene"].tolist())

    # 6. EA-GAT edge attention + directionality
    DATA_STORE["edge_attention"] = load_edge_attention(
        EAGAT_DIR / "edge_attention_tumor.tsv"
    )
    DATA_STORE["edge_directed"] = load_directionality(
        NETWORKS_DIR / "merged_network.txt"
    )

    meta_path = EAGAT_DIR / "meta.json"
    if meta_path.exists():
        with open(meta_path, encoding="utf-8") as fh:
            DATA_STORE["eagat_meta"] = json.load(fh)

    logger.info("Server Ready.")


# --- HELPERS ---

def _edge_meta(u: str, v: str, weight: float) -> Dict[str, Any]:
    """Attach weight, GAT attention, and directed flag for one edge."""
    attn_map = DATA_STORE["edge_attention"]
    dir_map = DATA_STORE["edge_directed"]

    attention = attn_map.get((u, v))
    if attention is None:
        attention = attn_map.get((v, u))
    # Fallback: use physical weight so the UI still scales edges
    if attention is None:
        attention = float(weight)

    directed = bool(dir_map.get((u, v), False) or dir_map.get((v, u), False))

    return {
        "weight": round(float(weight), 4),
        "attention": round(float(attention), 4),
        "directed": directed,
    }


def build_subgraph_json(G: nx.Graph, nodes: list) -> Dict[str, Any]:
    """Creates a Cytoscape/D3-compatible JSON for a subset of nodes."""
    subgraph = G.subgraph(nodes)

    log2fc_map = DATA_STORE["log2fc_map"]
    essential_genes = DATA_STORE["essential_genes"]
    eagat_by_gene = DATA_STORE["eagat_by_gene"]

    response_nodes = []
    for n in subgraph.nodes():
        rank = nodes.index(n) + 1 if n in nodes else 999
        log2fc = float(log2fc_map.get(n, 0.0))
        is_essential = n in essential_genes
        eagat = eagat_by_gene.get(str(n).upper(), {})

        response_nodes.append(
            {
                "id": n,
                "label": n,
                "rank": rank,
                "log2fc": log2fc,
                "is_essential": is_essential,
                "egatProb": eagat.get("egatProb"),
                "embeddingDrift": eagat.get("embeddingDrift"),
                "sksDisruption": eagat.get("sksDisruption"),
            }
        )

    response_links = []
    for u, v, d in subgraph.edges(data=True):
        meta = _edge_meta(u, v, d.get("weight", 1.0))
        response_links.append({"source": u, "target": v, **meta})

    return {"nodes": response_nodes, "links": response_links}


def get_path_details(G, path):
    log2fc_map = DATA_STORE["log2fc_map"]
    essential_genes = DATA_STORE["essential_genes"]

    nodes = [
        {
            "id": n,
            "label": n,
            "log2fc": float(log2fc_map.get(n, 0.0)),
            "is_essential": n in essential_genes,
        }
        for n in path
    ]
    links = []
    total_cost = 0.0

    for i in range(len(path) - 1):
        u, v = path[i], path[i + 1]
        if G.has_edge(u, v):
            weight = G[u][v].get("weight", 1.0)
            total_cost += weight
            meta = _edge_meta(u, v, weight)
            links.append({"source": u, "target": v, **meta})

    return {
        "found": True,
        "total_cost": round(total_cost, 4),
        "nodes": nodes,
        "links": links,
    }


# --- ENDPOINTS ---

@app.get("/api/health")
async def health():
    return {
        "status": "ok",
        "eagat_loaded": DATA_STORE["eagat_rankings"] is not None,
        "n_top_genes": len(DATA_STORE["top_genes"]),
        "n_attention_edges": len(DATA_STORE["edge_attention"]),
        "meta": DATA_STORE.get("eagat_meta", {}),
    }


@app.get("/api/eagat/leaderboard")
async def eagat_leaderboard(limit: int = Query(50, ge=1, le=500)):
    """
    Top genes by Embedding Drift (Δh) with EA-GAT probability and SKS fields.
    """
    df = DATA_STORE["eagat_rankings"]
    if df is None or df.empty:
        raise HTTPException(
            status_code=404,
            detail=(
                "EA-GAT rankings not found. Run: "
                "python src/python/pipeline_runner.py"
            ),
        )

    subset = df.head(limit)
    records = []
    for _, row in subset.iterrows():
        records.append(
            {
                "rank": int(row["rank"]),
                "gene": str(row["gene"]),
                "egatProb": round(float(row["egat_prob"]), 4),
                "embeddingDrift": round(float(row["embedding_drift"]), 4),
                "sksDisruption": round(float(row["sks_disruption"]), 4),
                # aliases kept for older UI code paths
                "ensembleScore": round(float(row["egat_prob"]), 4),
            }
        )
    return {
        "source": "eagat",
        "count": len(records),
        "genes": records,
        "meta": DATA_STORE.get("eagat_meta", {}),
    }


@app.get("/api/eagat/gene/{gene}")
async def eagat_gene(gene: str):
    """Per-gene EA-GAT metrics for the network explorer sidebar."""
    key = gene.upper().strip()
    # try exact / case-insensitive
    info = DATA_STORE["eagat_by_gene"].get(key)
    if info is None:
        # rankings may preserve original case
        for g, row in DATA_STORE["eagat_by_gene"].items():
            if g.upper() == key:
                info = row
                break
    if info is None:
        raise HTTPException(status_code=404, detail=f"Gene {gene} not in EA-GAT map")
    return info


@app.get("/api/network/initial")
async def get_initial_network(limit: int = 150):
    """
    Returns the skeleton network: Top N genes (EA-GAT ranks when available).
    Edges include weight, attention (α_ij), and directed flags.
    """
    top_genes = DATA_STORE["top_genes"]
    if not top_genes:
        return {"nodes": [], "links": [], "source": "empty"}

    subset_genes = top_genes[:limit]
    G = DATA_STORE["G_tumor"]
    payload = build_subgraph_json(G, subset_genes)
    payload["source"] = "eagat" if DATA_STORE["eagat_rankings"] is not None else "dice"
    return payload


@app.get("/api/knockout")
async def get_knockout():
    """Returns top virtual knockout candidates."""
    df = DATA_STORE["knockout_df"]
    if df is None or df.empty:
        return []

    df = df.copy()
    df["vitality_score"] = pd.to_numeric(df["vitality_score"], errors="coerce")
    sorted_df = df.sort_values(by="vitality_score", ascending=False).head(50)

    records = sorted_df.to_dict(orient="records")
    essential_genes = DATA_STORE["essential_genes"]

    for row in records:
        row["is_essential"] = row.get("gene") in essential_genes

    return records


@app.get("/api/path")
async def get_path(source: str, target: str):
    """Returns full graph data (nodes/edges) for paths in BOTH networks."""
    src = source.upper().strip()
    tgt = target.upper().strip()

    G_norm = DATA_STORE["G_normal"]
    G_tumor = DATA_STORE["G_tumor"]

    res = {
        "normal": {"found": False, "nodes": [], "links": [], "total_cost": 0},
        "tumor": {"found": False, "nodes": [], "links": [], "total_cost": 0},
        "stats": {"cost_diff": 0, "status": "Unknown"},
    }

    try:
        if src in G_norm and tgt in G_norm:
            p = nx.shortest_path(G_norm, src, tgt, weight="weight")
            res["normal"] = get_path_details(G_norm, p)
    except nx.NetworkXNoPath:
        pass

    try:
        if src in G_tumor and tgt in G_tumor:
            p = nx.shortest_path(G_tumor, src, tgt, weight="weight")
            res["tumor"] = get_path_details(G_tumor, p)
    except nx.NetworkXNoPath:
        pass

    if res["normal"]["found"] and res["tumor"]["found"]:
        diff = res["normal"]["total_cost"] - res["tumor"]["total_cost"]
        res["stats"]["cost_diff"] = round(diff, 4)
        if diff > 0.1:
            res["stats"]["status"] = "Shortcut Created (Gain of Function)"
        elif diff < -0.1:
            res["stats"]["status"] = "Path Lengthened (Loss of Efficiency)"
        else:
            res["stats"]["status"] = "Stable"

    return res


# Weight-delta threshold for continuous rewiring labels.
# Topology N/T is identical; status is driven by Δw on Pearson distance (1 − |r|).
REWIRING_TAU = 0.3


def _classify_edge_rewiring(
    gene: str,
    neighbor: str,
    G_n: nx.Graph,
    G_t: nx.Graph,
    tau: float = REWIRING_TAU,
) -> Dict[str, Any]:
    """
    Classify a hub–neighbor edge by continuous weight change.

    Status rules (weight = Pearson distance 1 − |r|):
      * gained   : w_tumor − w_normal > tau   (or tumor-only edge)
      * lost     : w_normal − w_tumor > tau   (or normal-only edge)
      * maintained : |w_tumor − w_normal| ≤ tau
    """
    in_n = G_n.has_edge(gene, neighbor)
    in_t = G_t.has_edge(gene, neighbor)

    w_n = float(G_n[gene][neighbor]["weight"]) if in_n else None
    w_t = float(G_t[gene][neighbor]["weight"]) if in_t else None

    # Topology-only cases (rare with current weighted nets)
    if in_n and not in_t:
        return {
            "status": "lost",
            "weight": w_n,
            "weight_normal": w_n,
            "weight_tumor": None,
            "delta_weight": None,
        }
    if in_t and not in_n:
        return {
            "status": "gained",
            "weight": w_t,
            "weight_normal": None,
            "weight_tumor": w_t,
            "delta_weight": None,
        }

    # Both present — continuous Δw classification
    delta = w_t - w_n  # tumor − normal
    if delta > tau:
        status = "gained"  # upregulated distance in tumor
    elif delta < -tau:
        status = "lost"  # downregulated (distance decreased) / weaker tumor shift
    else:
        status = "maintained"

    return {
        "status": status,
        "weight": (w_n + w_t) / 2.0,
        "weight_normal": w_n,
        "weight_tumor": w_t,
        "delta_weight": delta,
    }


@app.get("/api/rewiring")
async def get_rewiring(
    gene: str,
    tau: float = Query(REWIRING_TAU, ge=0.0, le=2.0),
    limit: int = Query(30, ge=1, le=200),
):
    """
    Ego network for a gene with weight-delta rewiring labels.

    Because Normal/Tumor share the same topology, status is defined on
    continuous edge weights (Pearson distance 1 − |r|):

      gained     : w_tumor − w_normal > tau
      lost       : w_normal − w_tumor > tau
      maintained : |Δw| ≤ tau

    Default ``tau=0.3``. Edges carry attention + directed metadata for the UI.
    """
    gene = gene.upper().strip()
    G_n = DATA_STORE["G_normal"]
    G_t = DATA_STORE["G_tumor"]

    if gene not in G_t and gene not in G_n:
        return {"nodes": [], "links": [], "tau": tau, "stats": {}}

    neighbors_n = set(G_n.neighbors(gene)) if gene in G_n else set()
    neighbors_t = set(G_t.neighbors(gene)) if gene in G_t else set()
    all_neighbors = neighbors_n.union(neighbors_t)

    log2fc_map = DATA_STORE["log2fc_map"]
    essential_genes = DATA_STORE["essential_genes"]
    eagat_by_gene = DATA_STORE["eagat_by_gene"]

    hub_eagat = eagat_by_gene.get(gene, {})
    nodes = [
        {
            "id": gene,
            "type": "hub",
            "log2fc": float(log2fc_map.get(gene, 0.0)),
            "is_essential": gene in essential_genes,
            "egatProb": hub_eagat.get("egatProb"),
            "embeddingDrift": hub_eagat.get("embeddingDrift"),
            "sksDisruption": hub_eagat.get("sksDisruption"),
            "rank": hub_eagat.get("rank"),
        }
    ]

    # Score every neighbor, then keep top-|Δw| for the ego view
    scored: List[Dict[str, Any]] = []
    for n in all_neighbors:
        cls = _classify_edge_rewiring(gene, n, G_n, G_t, tau=tau)
        abs_delta = (
            abs(cls["delta_weight"])
            if cls["delta_weight"] is not None
            else 1.0  # topology-only changes rank high
        )
        scored.append({"neighbor": n, "abs_delta": abs_delta, **cls})

    scored.sort(key=lambda r: r["abs_delta"], reverse=True)
    top = scored[:limit]

    links = []
    counts = {"gained": 0, "lost": 0, "maintained": 0}

    for row in top:
        n = row["neighbor"]
        status = row["status"]
        counts[status] = counts.get(status, 0) + 1

        n_eagat = eagat_by_gene.get(str(n).upper(), {})
        nodes.append(
            {
                "id": n,
                "type": "neighbor",
                "status": status,
                "log2fc": float(log2fc_map.get(n, 0.0)),
                "is_essential": n in essential_genes,
                "egatProb": n_eagat.get("egatProb"),
                "embeddingDrift": n_eagat.get("embeddingDrift"),
            }
        )
        meta = _edge_meta(gene, n, row["weight"] if row["weight"] is not None else 0.5)
        links.append(
            {
                "source": gene,
                "target": n,
                "status": status,
                "weight_normal": (
                    round(row["weight_normal"], 4)
                    if row["weight_normal"] is not None
                    else None
                ),
                "weight_tumor": (
                    round(row["weight_tumor"], 4)
                    if row["weight_tumor"] is not None
                    else None
                ),
                "delta_weight": (
                    round(row["delta_weight"], 4)
                    if row["delta_weight"] is not None
                    else None
                ),
                **meta,
            }
        )

    return {
        "nodes": nodes,
        "links": links,
        "tau": tau,
        "stats": {
            "total": len(links),
            "gained": counts["gained"],
            "lost": counts["lost"],
            "maintained": counts["maintained"],
            "n_neighbors_all": len(all_neighbors),
        },
    }


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
