#!/bin/bash

# ============================================================================
# DiCE-Flow Complete Pipeline (Aim 1)
# Multi-Scale Diffusion Descriptor, Chebyshev Trace Estimation, and FDR Control
# ============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

# ============================================================================
# Configuration & Setup
# ============================================================================

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

DATA_DIR="$PROJECT_ROOT/data"
RAW_DIR="$DATA_DIR/raw"
PROCESSED_DIR="$DATA_DIR/processed"
RESULTS_DIR="$DATA_DIR/results"

KEGG_DIR="$RAW_DIR/kegg"
BIOGRID_FILE="$RAW_DIR/biogrid/BIOGRID-ALL-5.0.262.tab3.txt"
NORMAL_EXPR="$RAW_DIR/expression/expression_normal.txt"
TUMOR_EXPR="$RAW_DIR/expression/expression_tumor.txt"

# ============================================================================
# Pre-Flight Check: Verify all datasets exist
# ============================================================================
log_info "Running Pre-Flight Dataset Checks..."

MISSING_DATA=0

# 1. Check KEGG Directory (must exist and not be empty)
if [ ! -d "$KEGG_DIR" ] || [ -z "$(ls -A "$KEGG_DIR" 2>/dev/null | grep -i "\.xml$" || true)" ]; then
    log_error "KEGG data missing: Directory $KEGG_DIR contains no XML files."
    MISSING_DATA=1
fi

# 2. Check BioGRID File
if [ ! -f "$BIOGRID_FILE" ]; then
    log_error "BioGRID data missing: File $BIOGRID_FILE does not exist."
    MISSING_DATA=1
fi

# 3. Check Expression Files
if [ ! -f "$NORMAL_EXPR" ]; then
    log_error "Normal Expression data missing: $NORMAL_EXPR"
    MISSING_DATA=1
fi
if [ ! -f "$TUMOR_EXPR" ]; then
    log_error "Tumor Expression data missing: $TUMOR_EXPR"
    MISSING_DATA=1
fi

if [ $MISSING_DATA -eq 1 ]; then
    echo -e "${RED}Pre-flight checks failed. Please download the required datasets and fix the paths.${NC}"
    exit 1
else
    log_info "All datasets verified! Proceeding with pipeline..."
fi

PYTHON="python3"

# Create necessary directories
mkdir -p "$PROCESSED_DIR/networks" "$PROCESSED_DIR/weighted" "$RESULTS_DIR/validation"

# Check expression data
if [ ! -f "$NORMAL_EXPR" ] || [ ! -f "$TUMOR_EXPR" ]; then
    log_error "Expression files not found in $RAW_DIR/expression/"
fi

# ============================================================================
# Phase 1: Differential Expression Analysis (For Downstream Metadata)
# ============================================================================
log_info "=== Phase 1: Differential Expression Analysis (DEA) ==="
# We run this just to get the Log2FC annotations for the web-viewer and plots.
# We DO NOT use this to filter the network anymore.

DEGS_FILE="$PROCESSED_DIR/dea_filtered_genes.txt"
DEA_RESULTS="$RESULTS_DIR/dea_results_table.txt"

$PYTHON src/python/dea_filter.py \
    --normal "$NORMAL_EXPR" \
    --tumor "$TUMOR_EXPR" \
    --output-genes "$DEGS_FILE" \
    --output-table "$DEA_RESULTS" \
    --p-val 0.05 \
    --lfc 1.0

# ============================================================================
# Phase 2: Network Scaffold Construction (Unweighted)
# ============================================================================
log_info "=== Phase 2: Network Scaffold Construction ==="

KEGG_OUTPUT="$PROCESSED_DIR/networks/kegg_interactions.txt"
BIOGRID_OUTPUT="$PROCESSED_DIR/networks/biogrid_interactions.txt"
MERGED_NETWORK="$PROCESSED_DIR/networks/merged_network.txt"

log_info "Parsing KEGG pathways..."
$PYTHON src/python/parse_kegg.py --kegg-dir "$KEGG_DIR" --output "$KEGG_OUTPUT"

log_info "Parsing BioGRID interactions..."
$PYTHON src/python/parse_biogrid.py --input "$BIOGRID_FILE" --output "$BIOGRID_OUTPUT" --organism 9606

log_info "Merging KEGG and BioGRID networks..."
$PYTHON src/python/merge_networks.py \
    --kegg "$KEGG_OUTPUT" \
    --biogrid "$BIOGRID_OUTPUT" \
    --output "$MERGED_NETWORK"

# ============================================================================
# Phase 3: Aim 1 - Multi-Scale Diffusion Descriptor & Null Models
# ============================================================================
log_info "=== Phase 3: Aim 1 - Multi-Scale Diffusion & Calibration ==="
# This script loads expression, calculates correlation weights on-the-fly, 
# computes Laplacians, runs HKS trace estimation on GPU, and applies FDR null models.

$PYTHON src/python/aim1_pipeline.py \
    --norm-expr "$NORMAL_EXPR" \
    --tum-expr "$TUMOR_EXPR" \
    --network "$MERGED_NETWORK" \
    --n-perms 50

# ============================================================================
# Phase 4: Extended Validation & Benchmarking
# ============================================================================
log_info "=== Phase 4: Validation & Benchmarking ==="

# 1. Fitness Gene Analysis
if [ -f "src/python/validation/fitness_analysis.py" ]; then
    log_info "Running Cancer Fitness Gene Overlap Analysis..."
    $PYTHON src/python/validation/fitness_analysis.py \
        --dice-genes "$RESULTS_DIR/dice_genes_top500.txt" \
        --output-dir "$RESULTS_DIR/validation"
fi

# 2. Dual-Dataset Validation (if alternative dataset exists)
if [ -f "src/python/validation/dual_dataset_analysis.py" ]; then
    if [ -f "$DATA_DIR/results_gse21032/dice_genes_top500.txt" ]; then
        log_info "Running Dual-Dataset Cross-Validation..."
        $PYTHON src/python/validation/dual_dataset_analysis.py \
            --dataset-a "$RESULTS_DIR/dice_genes_top500.txt" \
            --dataset-b "$DATA_DIR/results_gse21032/dice_genes_top500.txt" \
            --output-dir "$RESULTS_DIR/validation"
    fi
fi

# ============================================================================
# Summary
# ============================================================================
log_info "=== Pipeline Complete ==="
echo ""
echo "Results saved in: $RESULTS_DIR/"
echo "Output files:"
echo "  - dice_genes.txt                  Full candidate list (HKS Magnitude, Scale, FDR)"
echo "  - dice_genes_top500.txt           Top 500 significant candidates"
echo ""

# Display top 10 genes
echo "Top 10 Disease Genes (Sorted by Magnitude of Structural Shift):"
echo "================================================================"
$PYTHON -c "
import pandas as pd; 
try:
    df = pd.read_csv('$RESULTS_DIR/dice_genes_top500.txt', sep='\t')
    print(df[['gene', 'ensemble_score', 'characteristic_scale', 'q_value_fdr']].head(10).to_string(index=False))
except Exception:
    print('Failed to read results file.')
"
echo ""