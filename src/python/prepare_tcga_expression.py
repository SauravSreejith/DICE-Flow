#!/usr/bin/env python3
"""
Prepare full unfiltered TCGA-PRAD expression matrices for DiCE-Flow.

Source archive (already in the repo):
  data/raw/expression/TCGA-PRAD.gz
  → UCSC Xena ``TCGA.PRAD.sampleMap/HiSeqV2`` (gene × sample log2(norm_count+1))

This script does **not** download anything. It splits the archive by TCGA sample
type code (barcode suffix after the third hyphen):

  * ``-11``  Solid Tissue Normal  → expression_normal.txt
  * ``-01``  Primary Solid Tumor  → expression_tumor.txt

All ~20,530 genes are retained (including TP53, EGFR, MYC, …).

Usage (from project root, venv active)::

    python src/python/prepare_tcga_expression.py

    # custom paths / overwrite confirmation
    python src/python/prepare_tcga_expression.py \\
        --input data/raw/expression/TCGA-PRAD.gz \\
        --out-dir data/raw/expression \\
        --force

After regenerating matrices you must re-run network weighting + EA-GAT::

    # optional: full interactome without DEA gene filter
    ./run_pipeline.sh --global --skip-build

    # or only re-weight + EA-GAT if merged_network already exists
    python src/python/weight_network.py \\
        --network data/processed/networks/merged_network.txt \\
        --expression data/raw/expression/expression_normal.txt \\
        --output data/processed/weighted/network_normal.txt --scale 1000
    python src/python/weight_network.py \\
        --network data/processed/networks/merged_network.txt \\
        --expression data/raw/expression/expression_tumor.txt \\
        --output data/processed/weighted/network_tumor.txt --scale 1000
    python src/python/pipeline_runner.py
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger("prepare_tcga")

_ROOT = Path(__file__).resolve().parents[2]
_DEFAULT_INPUT = _ROOT / "data/raw/expression/TCGA-PRAD.gz"
_DEFAULT_OUT = _ROOT / "data/raw/expression"


def _sample_type_code(barcode: str) -> str:
    """
    Return the 2-digit TCGA sample-type code from a barcode.

    Examples
    --------
    TCGA-EJ-7125-11  →  '11'  (solid tissue normal)
    TCGA-XJ-A83F-01A →  '01'  (primary solid tumor)
    """
    parts = str(barcode).split("-")
    if len(parts) < 4:
        return ""
    # 4th field is SS (sample type) + optional vial letter, e.g. 01A → 01
    return parts[3][:2]


def load_xena_hiseqv2(path: Path) -> pd.DataFrame:
    """Load Xena HiSeqV2 matrix (genes × samples) from .gz or plain TSV."""
    logger.info("Reading %s (full matrix — this may take ~30s)…", path)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:  # type: ignore[arg-type]
        df = pd.read_csv(fh, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    logger.info("  → %d genes × %d samples", df.shape[0], df.shape[1])
    return df


def split_normal_tumor(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split columns by TCGA sample-type codes 11 (normal) and 01 (tumor)."""
    codes = {c: _sample_type_code(c) for c in df.columns}
    normal_cols = [c for c, code in codes.items() if code == "11"]
    tumor_cols = [c for c, code in codes.items() if code == "01"]

    # Fallback: if vial-stripped codes are empty, keep any column containing -11- / ending
    if not normal_cols:
        normal_cols = [c for c in df.columns if "-11" in c]
    if not tumor_cols:
        tumor_cols = [c for c in df.columns if "-01" in c]

    logger.info("Normal samples (-11): %d", len(normal_cols))
    logger.info("Tumor  samples (-01): %d", len(tumor_cols))

    if not normal_cols or not tumor_cols:
        raise ValueError(
            "Could not find both normal (-11) and tumor (-01) TCGA barcodes. "
            f"Example columns: {list(df.columns[:5])}"
        )

    normal = df[normal_cols].copy()
    tumor = df[tumor_cols].copy()
    # Drop genes that are all-NaN in a split (keep full gene set otherwise)
    normal = normal.dropna(how="all")
    tumor = tumor.dropna(how="all")
    # Align to shared genes so EA-GAT vocabulary is consistent
    shared = normal.index.intersection(tumor.index)
    normal = normal.loc[shared]
    tumor = tumor.loc[shared]
    logger.info("Shared genes retained: %d", len(shared))
    for g in ("TP53", "EGFR", "MYC", "KRAS"):
        logger.info("  %s present: %s", g, g in shared)
    return normal, tumor


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Split TCGA-PRAD HiSeqV2 into full unfiltered normal/tumor TSVs."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=_DEFAULT_INPUT,
        help="Path to TCGA-PRAD.gz (Xena HiSeqV2)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=_DEFAULT_OUT,
        help="Output directory for expression_*.txt",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing expression_normal/tumor.txt",
    )
    args = parser.parse_args(argv)

    if not args.input.is_file():
        logger.error("Input not found: %s", args.input)
        logger.error(
            "Download from Xena: hub=https://tcga.xenahubs.net "
            "dataset=TCGA.PRAD.sampleMap/HiSeqV2  → save as TCGA-PRAD.gz"
        )
        return 1

    out_n = args.out_dir / "expression_normal.txt"
    out_t = args.out_dir / "expression_tumor.txt"
    if (out_n.exists() or out_t.exists()) and not args.force:
        logger.error(
            "Output exists. Re-run with --force to overwrite:\n  %s\n  %s",
            out_n,
            out_t,
        )
        return 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = load_xena_hiseqv2(args.input)
    normal, tumor = split_normal_tumor(df)

    normal.to_csv(out_n, sep="\t")
    tumor.to_csv(out_t, sep="\t")
    logger.info("Wrote %s  (%d genes × %d samples)", out_n, *normal.shape)
    logger.info("Wrote %s  (%d genes × %d samples)", out_t, *tumor.shape)
    logger.info(
        "Next: re-weight networks and re-run EA-GAT "
        "(see docstring in this script)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
