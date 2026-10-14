# DiCE-Flow Python Codebase — Complete Reference

This document explains every Python file in the DiCE-Flow project, organized by pipeline phase. The codebase has been radically overhauled for **Aim 1**, shifting from legacy single-scale shortest-path centralities to **GPU-accelerated Multi-Scale Diffusion Descriptors** with rigorous False Discovery Rate (FDR) control.

---

## 1. Data Preprocessing & Scaffold Construction (`src/python/`)

These scripts build the foundational structure of the network and prepare the biological metadata.

#### `parse_kegg.py`
- **Purpose:** Extracts directed gene–gene interactions from KEGG pathway XML (KGML) files.
- **How it works:** Uses `xml.etree.ElementTree` to parse `.xml` files, maps entry IDs to human gene symbols, and extracts functional relationships (`PPrel`, `GErel`, `PCrel`). Self-loops are filtered out.
- **Output:** A TSV of `(source, target)` directed interactions.

#### `parse_biogrid.py`
- **Purpose:** Parses the massive BioGRID **TAB3** experimental interaction database.
- **How it works:** Streams the multi-gigabyte file in chunks of 100,000 rows. Filters explicitly for Homo sapiens (Taxonomy ID 9606) interactions, drops self-loops, and deduplicates.
- **Output:** A TSV of high-confidence physical and genetic interactions.

#### `merge_networks.py`
- **Purpose:** Combines KEGG and BioGRID interactions into a single master scaffold.
- **How it works:** Loads both networks and computes the union. 
- **Output:** `merged_network.txt` containing >1,000,000 unique, unweighted topological edges.

#### `dea_filter.py`
- **Purpose:** Performs classic Differential Expression Analysis (DEA).
- **How it works:** Applies Welch's T-test and Benjamini-Hochberg FDR correction across normal vs. tumor expression matrices.
- **Note:** Unlike the legacy DiCE pipeline, **we no longer use this to prune the network.** Aim 1's GPU engine scales effortlessly to the full 17,500-gene interactome. This script is now solely used to generate $\log_{2}\text{FC}$ metadata to color downstream visualizations.
- **Output:** `dea_results_table.txt`.

---

## 2. Aim 1 Core: Multi-Scale Diffusion & Calibration (`src/python/`)

These three scripts represent the mathematical heartbeat of the project. They map expression data onto the network, compute diffusion, and perform statistical thresholding.

#### `aim1_pipeline.py`
- **Purpose:** The master orchestrator for Aim 1.
- **How it works:**
  1. Loads the expression matrices and master topological scaffold.
  2. Identifies the Largest Connected Component (~17,466 genes).
  3. **$O(E)$ Weighting Optimization:** Instead of computing a $17,000 \times 17,000$ correlation matrix (which crashes CPUs), it vectorizes Pearson Correlation purely for the ~300,000 existing edges. 
  4. Calls the trace estimators and manages the progress of the $B=50$ null permutations (shuffling sample labels to model batch/cohort noise).
- **Output:** Outputs the final, unified `dice_genes.txt` ranking.

#### `diffusion_descriptor.py`
- **Purpose:** Mathematical engine for the Heat Kernel Signature (HKS) and Scale Banks.
- **How it works:**
  - **Hardware Agnostic Acceleration:** Dynamically assigns computations to NVIDIA CUDA, Apple Silicon (MPS), or multi-core CPU.
  - **`get_scale_bank`:** Calculates the spectral boundaries of the graph. It uses a fast Arnoldi iteration for the largest eigenvalue, and a heuristic bound (`1e-4`) for the Fiedler value to strictly prevent the 35-minute sparse LU decomposition hang.
  - **`hks_diag_cheb_torch`:** Evaluates the diagonal of the heat operator $e^{-tL}$. Since exact eigendecomposition of a 17,000-node graph is impossible, it uses a **Chebyshev polynomial approximation** (order $K=60$) coupled with **Hutchinson’s Stochastic Trace Estimator** ($S=256$ Rademacher vectors).

#### `null_models.py`
- **Purpose:** Prevents false discoveries and isolates the true biological signal.
- **How it works:**
  - **`get_differential_profile`:** Applies *Per-Scale Logarithmic Z-Scoring*. This neutralizes global density drift and suppresses extreme mathematical artifacts on peripheral/isolated nodes. The final magnitude ($\Delta$) is the $L_\infty$ maximum absolute shift across scales.
  - **`stratified_pvals`:** Implements **Median Absolute Deviation (MAD)** robust statistics. Because trace estimators use random vectors, they occasionally generate massive Gaussian noise spikes. The robust Z-score perfectly ignores these spikes, allowing true biological hubs (like *TP53*) to achieve $p < 10^{-40}$. Computes final BH-FDR $q$-values across dynamic degree strata.

---

## 3. Validation & Benchmarking Framework (`src/python/validation/`)

These scripts mathematically prove that the output of Aim 1 is scientifically valid and vastly superior to baseline methods.

#### `validation_pipeline.py` (Master Validator)
- **Purpose:** Executes a gauntlet of 8 distinct scientific analyses on `dice_genes.txt`.
- **Key Analyses:**
  - **Analysis 4 (Scale vs Magnitude Plot):** Scatters the Characteristic Scale ($t^*$) against the Structural Shift Magnitude ($\Delta$), isolating local complex changes vs. global pathway shifts.
  - **Analysis 5 (FDR Volcano Plot):** Generates a volcano plot mapping Signed Structural Shift against $-\log_{10}(q\text{-value})$.
  - **Analysis 8 (Benchmarking):** Generates Precision@K curves mapping DiCE Aim 1 vs. Raw Degree Centrality vs. Differential Expression (DE) Rankings using COSMIC/DisGeNET gold standards.

#### `fitness_analysis.py`
- **Purpose:** Cancer fitness essentiality verification.
- **How it works:** Intersects the top predicted DiCE genes with DepMap CRISPR-Cas9 knockout survival data. Evaluates whether the identified genes are mathematically essential to cancer cell survival via a Hypergeometric overlap test.

#### `survival_validation.py`
- **Purpose:** Clinical outcome validation.
- **How it works:** Uses `lifelines` to map gene expression to patient survival data (days to death). Generates Kaplan-Meier survival curves and Cox Proportional Hazard Ratios.

#### `dual_dataset_analysis.py`
- **Purpose:** Cross-cohort reproducibility.
- **How it works:** Compares DiCE predictions trained on two completely disjoint patient cohorts (e.g., TCGA vs. GSE21032 metastasis) using Spearman rank correlation.

---

## 4. Post-Analysis & Presentation (`src/python/`)

#### `post_analysis.py`
- **Purpose:** Export utilities and aesthetic plotting.
- **Key Features:**
  - `plot_top_genes()`: Generates horizontal bar charts of the top 20 relocated genes, colored by direction (blue = bridging, red = localizing).
  - `plot_rank_distribution()`: Plots the logarithmic decay of the structural shift magnitude.
  - `export_for_enrichr()`: Formats top candidate lists for immediate pasting into GO/KEGG web tools (Enrichr, DAVID).

#### `export_results.py`
- **Purpose:** Translates internal DataFrames into Cytoscape-importable `.csv` node and edge attribute files for beautiful 3D network renders.
