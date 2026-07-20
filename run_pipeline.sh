#!/bin/bash

# ============================================================================
# DiCE-Duan Complete Pipeline
# Integrates KEGG/BioGRID networks with gene expression for differential
# centrality analysis using the Duan 2025 SSSP algorithm
# ============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Log functions
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# ============================================================================
# Configuration
# ============================================================================

# Project root
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

# Directories
DATA_DIR="$PROJECT_ROOT/data"
RAW_DIR="$DATA_DIR/raw"
PROCESSED_DIR="$DATA_DIR/processed"
RESULTS_DIR="$DATA_DIR/results"

KEGG_DIR="$RAW_DIR/kegg"
BIOGRID_FILE="$RAW_DIR/biogrid/BIOGRID-ALL-5.0.259.tab2.txt"
EXPR_DIR="$RAW_DIR/expression"

# Python and C++ binaries
PYTHON="python3"
CPP_BIN="$PROJECT_ROOT/build/dice_analyzer"

# Create directories
mkdir -p "$PROCESSED_DIR/networks" "$PROCESSED_DIR/weighted" "$RESULTS_DIR"

# ============================================================================
# Parse command line arguments
# ============================================================================

MODE="full"  # full, test, or custom
USE_DUMMY=false
SKIP_BUILD=false
GLOBAL_MODE=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --test)
            MODE="test"
            USE_DUMMY=true
            shift
            ;;
        --global)
            GLOBAL_MODE=1
            shift
            ;;
        --skip-build)
            SKIP_BUILD=true
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --test          Run with dummy data"
            echo "  --global        Run on full interactome (skip filtering)"
            echo "  --skip-build    Skip C++ compilation"
            echo "  --help          Show this help message"
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [ "$GLOBAL_MODE" -eq 1 ]; then
    log_warn "Running in GLOBAL mode. Bypassing DEA and IG filters to compute the entire human interactome."
fi

# ============================================================================
# Step 0: Build C++ Engine
# ============================================================================

if [ "$SKIP_BUILD" = false ]; then
    log_info "Building C++ Duan engine..."
    
    mkdir -p build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j$(nproc)
    cd "$PROJECT_ROOT"
    
    if [ ! -f "$CPP_BIN" ]; then
        log_error "Failed to build C++ engine"
        exit 1
    fi
    
    log_info "C++ engine built successfully"
else
    log_warn "Skipping C++ build"
fi

# ============================================================================
# Step 0.5: Generate dummy data if in test mode
# ============================================================================

if [ "$USE_DUMMY" = true ]; then
    log_info "Generating dummy data for testing..."
    
    $PYTHON src/python/generate_dummy_data.py \
        --output-dir "$DATA_DIR/dummy" \
        --n-genes 100 \
        --n-samples 30
    
    # Update paths to use dummy data
    KEGG_DIR="$DATA_DIR/dummy/kegg"
    BIOGRID_FILE="$DATA_DIR/dummy/biogrid/BIOGRID-dummy.tab2.txt"
    EXPR_DIR="$DATA_DIR/dummy/expression"
fi

# ============================================================================
# Phase 1 & 2: Differential Expression and Information Gain Filtering
# ============================================================================

log_info "=== Phase 1: Differential Expression Analysis (DEA) ==="

# Check for expression files upfront since we need them to filter
NORMAL_EXPR="$EXPR_DIR/expression_normal.txt"
TUMOR_EXPR="$EXPR_DIR/expression_tumor.txt"

if [ ! -f "$NORMAL_EXPR" ] || [ ! -f "$TUMOR_EXPR" ]; then
    log_error "Expression files not found in $EXPR_DIR"
    exit 1
fi

# Step 1.2: Perform DEA
log_info "Performing differential expression analysis (T-Test, FDR < 0.05, |LFC| > 1)..."
DEGS_FILE="$PROCESSED_DIR/dea_filtered_genes.txt"
DEA_RESULTS="$RESULTS_DIR/dea_results_table.txt"

$PYTHON src/python/dea_filter.py \
    --normal "$NORMAL_EXPR" \
    --tumor "$TUMOR_EXPR" \
    --output-genes "$DEGS_FILE" \
    --output-table "$DEA_RESULTS" \
    --p-val 0.05 \
    --lfc 1.0

log_info "=== Phase 2: Information Gain (IG) Filtering ==="

# Filter the DEGs via Information Gain (Mutual Information)
IG_CANDIDATES="$PROCESSED_DIR/ig_candidates.txt"

log_info "Filtering DEGs through Information Gain threshold..."
$PYTHON src/python/pre_filter.py \
    --normal "$NORMAL_EXPR" \
    --tumor "$TUMOR_EXPR" \
    --output "$IG_CANDIDATES"

# Since pre_filter.py takes the raw expression matrix, it calculates IG on all common genes.
# We must take the intersection of DEGs and IG_CANDIDATES to match the paper exactly.
FINAL_CANDIDATES="$PROCESSED_DIR/final_candidates.txt"
comm -12 <(sort "$DEGS_FILE") <(sort "$IG_CANDIDATES") > "$FINAL_CANDIDATES"

log_info "Final filtered candidate pool: $(wc -l < "$FINAL_CANDIDATES") genes"

# ============================================================================
# Phase 3: Data Engineering - Network Construction (FILTERED)
# ============================================================================

log_info "=== Phase 3: Network Construction (Sub-Network) ==="

# Step 3.1: Parse KEGG
log_info "Parsing KEGG pathways..."
KEGG_OUTPUT="$PROCESSED_DIR/networks/kegg_interactions.txt"

$PYTHON src/python/parse_kegg.py \
    --kegg-dir "$KEGG_DIR" \
    --output "$KEGG_OUTPUT"

# Step 3.2: Parse BioGRID
log_info "Parsing BioGRID interactions..."
BIOGRID_OUTPUT="$PROCESSED_DIR/networks/biogrid_interactions.txt"

$PYTHON src/python/parse_biogrid.py \
    --input "$BIOGRID_FILE" \
    --output "$BIOGRID_OUTPUT" \
    --organism 9606

# Step 3.3: Merge networks
if [ "$GLOBAL_MODE" -eq 1 ]; then
    log_info "Merging KEGG and BioGRID networks (GLOBAL: No candidate filter)..."
    MERGED_NETWORK="$PROCESSED_DIR/networks/merged_network.txt"

    $PYTHON src/python/merge_networks.py \
        --kegg "$KEGG_OUTPUT" \
        --biogrid "$BIOGRID_OUTPUT" \
        --output "$MERGED_NETWORK"
else
    log_info "Merging KEGG and BioGRID networks (Filtered to candidates)..."
    MERGED_NETWORK="$PROCESSED_DIR/networks/merged_network.txt"

    $PYTHON src/python/merge_networks.py \
        --kegg "$KEGG_OUTPUT" \
        --biogrid "$BIOGRID_OUTPUT" \
        --candidates "$FINAL_CANDIDATES" \
        --output "$MERGED_NETWORK"
fi

if [ ! -f "$MERGED_NETWORK" ]; then
    log_error "Network merging failed"
    exit 1
fi

# ============================================================================
# Phase 4: Network Weighting with Expression Data
# ============================================================================

log_info "=== Phase 4: Network Weighting ==="

# Step 4.1: Weight network with normal expression
log_info "Weighting network with normal expression..."
NORMAL_WEIGHTED="$PROCESSED_DIR/weighted/network_normal.txt"

$PYTHON src/python/weight_network.py \
    --network "$MERGED_NETWORK" \
    --expression "$NORMAL_EXPR" \
    --output "$NORMAL_WEIGHTED" \
    --scale 1000

# Step 2.2: Weight network with tumor expression
log_info "Weighting network with tumor expression..."
TUMOR_WEIGHTED="$PROCESSED_DIR/weighted/network_tumor.txt"

$PYTHON src/python/weight_network.py \
    --network "$MERGED_NETWORK" \
    --expression "$TUMOR_EXPR" \
    --output "$TUMOR_WEIGHTED" \
    --scale 1000

# ============================================================================
# Phase 5: Centrality Computation (C++ Engine)
# ============================================================================

log_info "=== Phase 5: Computing Network Centralities ==="

NORMAL_INT="${NORMAL_WEIGHTED/.txt/_int.txt}"
TUMOR_INT="${TUMOR_WEIGHTED/.txt/_int.txt}"

if [ ! -f "$CPP_BIN" ]; then
    log_error "C++ engine not found. Please build first."
    exit 1
fi

log_info "Running Duan 2025 algorithm on normal and tumor networks..."

$CPP_BIN \
    "$NORMAL_INT" \
    "$TUMOR_INT" \
    "$RESULTS_DIR"

if [ $? -ne 0 ]; then
    log_error "C++ analysis failed"
    exit 1
fi

# ============================================================================
# Phase 6: Differential Centrality Analysis
# ============================================================================

log_info "=== Phase 6: Differential Centrality Analysis ==="

NORMAL_CENTRALITY="$RESULTS_DIR/centrality_normal.txt"
TUMOR_CENTRALITY="$RESULTS_DIR/centrality_tumor.txt"

if [ ! -f "$NORMAL_CENTRALITY" ] || [ ! -f "$TUMOR_CENTRALITY" ]; then
    log_error "Centrality files not found"
    exit 1
fi

log_info "Computing differential centrality..."
DICE_GENES="$RESULTS_DIR/dice_genes.txt"

if [ "$GLOBAL_MODE" -eq 1 ]; then
    log_info "GLOBAL MODE: Disabling Phase 5 noise filter to retain hyper-hubs like EP300..."
    $PYTHON src/python/differential_centrality.py \
        --normal "$NORMAL_CENTRALITY" \
        --tumor "$TUMOR_CENTRALITY" \
        --output "$DICE_GENES" \
        --top-n 500 \
        --disable-noise
else
    $PYTHON src/python/differential_centrality.py \
        --normal "$NORMAL_CENTRALITY" \
        --tumor "$TUMOR_CENTRALITY" \
        --output "$DICE_GENES" \
        --top-n 500
fi

if [ ! -f "$DICE_GENES" ]; then
    log_error "Differential centrality analysis failed"
    exit 1
fi

# ============================================================================
# Phase 7: Virtual Knockout Analysis (Optional)
# ============================================================================

log_info "=== Phase 7: Virtual Knockout Analysis ==="

# Create gene list from top DiCE genes
TOP_GENES="$RESULTS_DIR/genes_for_knockout.txt"
tail -n +2 "$RESULTS_DIR/dice_genes_top500.txt" | cut -f1 | head -50 > "$TOP_GENES"

log_info "Running virtual knockout analysis on top 50 DiCE genes..."

$CPP_BIN \
    "$NORMAL_INT" \
    "$TUMOR_INT" \
    "$RESULTS_DIR" \
    "knockout" \
    "$TOP_GENES"

# ============================================================================
# Phase 8: Extended Validation (Fitness Genes + Dual Dataset)
# ============================================================================

log_info "=== Phase 8: Extended Validation ==="

# Run cancer fitness gene analysis
if [ -f "fitness_analysis.py" ]; then
    log_info "Running cancer fitness gene overlap analysis..."
    $PYTHON fitness_analysis.py \
        --dice-genes "$RESULTS_DIR/dice_genes_top500.txt" \
        --output-dir "$RESULTS_DIR/validation" \
        2>&1 | while read line; do log_info "  $line"; done
    log_info "Fitness analysis complete."
else
    log_info "fitness_analysis.py not found, skipping."
fi

# Run dual-dataset analysis if second dataset exists
if [ -f "dual_dataset_analysis.py" ]; then
    ALT_RESULTS=""
    for ALT_DIR in "$DATA_DIR/results_gse21032" "$DATA_DIR/results_alt" "$DATA_DIR/results_metastasis"; do
        if [ -f "$ALT_DIR/dice_genes_top500.txt" ]; then
            ALT_RESULTS="$ALT_DIR/dice_genes_top500.txt"
            break
        fi
    done

    if [ -n "$ALT_RESULTS" ]; then
        log_info "Running dual-dataset cross-validation..."
        $PYTHON dual_dataset_analysis.py \
            --dataset-a "$RESULTS_DIR/dice_genes_top500.txt" \
            --dataset-b "$ALT_RESULTS" \
            --output-dir "$RESULTS_DIR/validation" \
            2>&1 | while read line; do log_info "  $line"; done
        log_info "Dual-dataset analysis complete."
    else
        log_info "No second dataset found. To enable dual-dataset analysis:"
        log_info "  Run the pipeline on GSE21032 data and place results in data/results_gse21032/"
    fi
fi

# ============================================================================
# Summary
# ============================================================================

log_info "=== Pipeline Complete ==="
echo ""
echo "Results saved in: $RESULTS_DIR/"
echo ""
echo "Output files:"
echo "  - dice_genes.txt                  Full DiCE gene list with ranks"
echo "  - dice_genes_top500.txt           Top 500 DiCE genes"
echo "  - centrality_normal.txt           Normal network centralities"
echo "  - centrality_tumor.txt            Tumor network centralities"
if [ -f "$RESULTS_DIR/knockout_analysis.txt" ]; then
echo "  - knockout_analysis.txt           Virtual knockout vitality scores"
fi
echo ""

# Display top 10 DiCE genes
echo "Top 10 DiCE genes (Score 0-1 implies Rank Product):"
echo "==================="

# Python one-liner to safely print the correct columns by name
$PYTHON -c "
import pandas as pd; 
df = pd.read_csv('$RESULTS_DIR/dice_genes_top500.txt', sep='\t'); 
print(df[['gene', 'ensemble_score']].head(10).to_string(index=False))
"
echo ""
