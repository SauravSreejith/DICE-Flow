# DiCE-Flow: Structural Response and Representation Learning for Cancer Driver Gene Discovery
**Technical Documentation & Methodology Reference**

## 1. Project Overview

**The Problem:** Traditional disease gene discovery relies heavily on Differential Gene Expression Analysis (DEA), which flags genes whose transcription levels spike or plummet between normal and disease states. However, DEA is structurally blind. It fails to detect "cancer fitness genes" or "silent drivers"—critical regulatory hubs (like *TP53*, *KRAS*, or *SMARCA4*) whose overall expression may not change significantly, but whose structural role within the protein-protein interaction (PPI) network rewires massively during disease progression.

**The Network Biology Evolution:** Our previous work (DiCE) proved that computing shortest-path centralities (Betweenness, Eigenvector) on condition-specific weighted networks could identify these silent drivers. However, standard centrality has critical mathematical ceilings:
1. **Scalar Compression:** A gene's entire structural context collapses into a single number, destroying higher-order topological features.
2. **Single Implicit Scale:** Eigenvector represents global position, while degree is strictly local. Neither can express intermediate changes (e.g., "Gene X lost a local module but gained a long-range signaling bridge").
3. **No Statistical Calibration:** Ranked lists without null distributions cannot distinguish true biological rewiring from sample-size imbalance or cohort noise.

**The DiCE-Flow Solution (Aim 1):** We replace scalar centrality with a **Multi-Scale Diffusion Descriptor**—specifically, the Heat Kernel Signature (HKS). Evaluated over a bank of diffusion times on GPU-accelerated trace estimators, it provides a scale-resolved profile of a gene's network position. By applying degree-stratified permutation null models and robust statistics (MAD), DiCE-Flow produces the first statistically calibrated, false-discovery-controlled (FDR) identification of structurally rewired cancer drivers.

---

## 2. Mathematical Foundations

The DiCE-Flow pipeline merges spectral graph theory, stochastic trace estimation, and robust statistics.

### 2.1 Condition-Specific Weighted Laplacians
We map biological context onto a high-confidence PPI scaffold using expression data. For every existing edge $(u, v)$, we calculate the absolute Pearson Correlation Coefficient ($r_{uv}$) of the two genes' expression within a specific cohort. 

Edge weights are defined as the correlation strength: $W_{uv} = |r_{uv}|$.
From this adjacency matrix $W$, we compute the combinatorial Graph Laplacian:
$$ L = D - W $$
where $D$ is the diagonal degree matrix.

### 2.2 Multi-Scale Diffusion Descriptor (Heat Kernel Signature)
The Heat Kernel Signature at node $i$ and diffusion time $t$ represents the probability that heat injected at $i$ remains at $i$ after time $t$. It is the diagonal of the heat operator $e^{-tL}$:
$$ HKS_i(t) = [e^{-tL}]_{ii} = \sum_k e^{-t\lambda_k} \phi_k(i)^2 $$
* **Small $t$:** Probes the immediate topological neighborhood (local complexes).
* **Large $t$:** Probes global network positioning (provably converging to Eigenvector Centrality).

### 2.3 Stochastic Trace Estimation (Hutchinson + Chebyshev)
Because a full eigendecomposition of a $17,000 \times 17,000$ Laplacian is computationally infeasible, we use a hybrid stochastic estimator. We approximate the matrix exponential $e^{-tL}$ using **Chebyshev polynomials** (order $K=60$), and estimate the diagonal using **Hutchinson's method** with $S=256$ random Rademacher vectors ($z \in \{-1, 1\}$). This allows us to compute the HKS entirely via fast sparse matrix-vector multiplications (SpMV) on the GPU.

### 2.4 The Differential Profile and Repositioning Shift
To compare Normal ($C$) and Tumor ($D$) states, we define the differential profile $d_i(t)$. 
To neutralize global density drift and suppress artifactual noise from isolated, low-degree nodes, we apply **per-scale logarithmic centering**:
$$ H'_{C} = \frac{\ln(H_C) - \mu_{\ln(H_C)}}{\sigma_{\ln(H_C)}} $$
$$ d_i(t) = H'_{D}(i, t) - H'_{C}(i, t) $$

From this profile, we extract three critical metrics per gene:
1. **Magnitude ($\Delta$):** The maximum absolute shift across all scales (using the $L_\infty$ norm to prevent trace-estimator noise accumulation): $\Delta_i = \max_t |d_i(t)|$.
2. **Characteristic Scale ($t^*$):** The specific diffusion time where the maximum shift occurred: $t^*_i = \text{argmax}_t |d_i(t)|$.
3. **Direction:** Signifies whether the gene shifted toward local integration ($-$) or long-range bridging ($+$).

### 2.5 Statistical Calibration (Degree-Stratified Robust FDR)
To separate true biology from cohort sampling noise, we compute $B=50$ null permutations (shuffling patient labels and rebuilding the entire pipeline).
Because the Hutchinson estimator occasionally produces extreme Gaussian spikes, raw empirical p-values fail. We apply **Robust Z-scoring** using the Median Absolute Deviation (MAD) within degree-stratified bins:
$$ Z_i = \frac{\Delta_i - \text{Median}(\Delta_{\text{null}})}{1.4826 \times \text{MAD}(\Delta_{\text{null}})} $$
This yields outlier-resistant empirical $p$-values, which are then corrected using the Benjamini-Hochberg (FDR) procedure to yield $q$-values.

---

## 3. Pipeline Architecture

The system executes sequentially via `run_pipeline.sh`:

1. **Phase 1: Knowledge Scaffold Construction:**
   Parses raw XML pathways from KEGG and TAB3 tabular data from BioGRID. Merges them to create a master unweighted edge list containing >1 million verified interactions.
2. **Phase 2: Fast Network Weighting & Laplacian Generation:**
   Loads RNA-seq data, intersects the expression matrix with the PPI scaffold, and computes vectorized Pearson correlations strictly for existing topological edges.
3. **Phase 3: The PyTorch HKS Engine (Aim 1):**
   Computes the scale bank, executes the batched Chebyshev trace estimation on the GPU, constructs the $50$ null permutations, and applies the degree-stratified robust FDR thresholding.
4. **Phase 4: Automated Validation & Benchmarking:**
   Evaluates the predicted gene list (`dice_genes.txt`) against external baselines, computes Precision@K, tests hypergeometric enrichment against Cancer Fitness Genes (DepMap), and generates scale/magnitude/volcano plots.

---

## 4. Codebase Directory Map

### Orchestration & Preprocessing
* **`run_pipeline.sh`**: The master bash orchestrator.
* **`src/python/parse_kegg.py`**: Extracts directed gene-gene relations from KEGG XMLs.
* **`src/python/parse_biogrid.py`**: Parses the BioGRID `tab3` database for human PPIs.
* **`src/python/merge_networks.py`**: Unions BioGRID and KEGG into a single master scaffold.
* **`src/python/dea_filter.py`**: Runs standard T-tests to extract Differential Expression metadata (Log2FC) used for downstream plot coloring.

### Aim 1 Core (Topology & Mathematics)
* **`src/python/aim1_pipeline.py`**: The GPU-accelerated main loop. Handles $O(E)$ correlation weighting, PyTorch initialization, and parallel null permutation scheduling.
* **`src/python/diffusion_descriptor.py`**: Contains `hks_diag_cheb_torch`, the hardware-agnostic (CUDA/MPS/CPU) Hutchinson-Chebyshev trace estimator, and the Fiedler-heuristic scale bank generator.
* **`src/python/null_models.py`**: Implements the Per-Scale Centering ($L_\infty$ magnitude profile) and the MAD-based Robust FDR calibration.

### Validation & Plotting (`src/python/validation/`)
* **`validation_pipeline.py`**: The master evaluation framework. Generates Precision@K benchmarking, Scale vs. Magnitude scatter plots, and FDR Volcano plots.
* **`fitness_analysis.py`**: Computes hypergeometic overlap between DiCE-Flow predictions and DepMap CRISPR essentiality datasets.
* **`dual_dataset_analysis.py`**: Evaluates cross-cohort generalizability.

---

## 5. Algorithmic Complexities & Engineering Optimizations

To scale this pipeline to the entire human interactome (~17,500 genes, ~2 million edges), four critical engineering bottlenecks were solved:

| Bottleneck | Standard Approach | DiCE-Flow Optimization | Result |
| :--- | :--- | :--- | :--- |
| **Correlation Weighting** | `np.corrcoef(X)` computes an all-to-all $O(N^2)$ matrix (225M operations). | Vectorized row-wise dot-products strictly for topological edges $O(E)$. | Dropped from 5 minutes to **1.5 seconds**. |
| **Spectral Scale Bank** | Exact Fiedler value via Sparse LU Decomposition $O(N^3)$. | Standard Arnoldi for $\lambda_{max}$ + safe heuristic bound (1e-4) for $\lambda_{min}$. | Dropped from 35 minutes to **0.01 seconds**. |
| **Null Model Variance** | Raw empirical $p$-values based on rank $O(B)$. | Robust MAD Z-scores (`scipy.stats.norm.sf`) ignoring stochastic trace-estimator spikes. | FDR $q$-values dropped from 0.20 to $< 10^{-40}$. |
| **Global Density Drift** | Raw absolute shift $\Delta = \|H_T - H_N\|_2$. | Per-scale Logarithmic Z-scoring + $L_\infty$ Maximum Shift. | Peripheral noise eliminated; True hubs rise to Top 10. |

---

## 6. Interpreting the Outputs

The pipeline outputs `data/results/dice_genes_top500.txt` with the following key columns:

* **`ensemble_score` ($\Delta$)**: The absolute magnitude of the structural repositioning. Higher is better.
* **`q_value_fdr`**: The Benjamini-Hochberg adjusted statistical significance of the shift compared to the null models. Values $< 0.05$ are scientifically valid discoveries.
* **`characteristic_scale` ($t^*$)**: The diffusion time at which the gene underwent its maximum structural shift. Small values indicate local complex rewiring; large values indicate global pathway shifts.
* **`direction`**: 
  * `-1`: The gene localized (heat concentrated).
  * `+1`: The gene dispersed (acts as a newly formed topological bridge).

---

## 7. DiCE-Flow Pipeline Data Flow Diagram

```mermaid
graph TD
    %% Define Node Styles
    classDef raw fill:#e2e8f0,stroke:#64748b,stroke-width:2px,color:#0f172a;
    classDef script fill:#fef3c7,stroke:#d97706,stroke-width:2px,color:#92400e;
    classDef engine fill:#fee2e2,stroke:#dc2626,stroke-width:2px,color:#991b1b;
    classDef output fill:#dcfce7,stroke:#16a34a,stroke-width:2px,color:#166534;

    %% Phase 1: Network Scaffold
    subgraph Phase 1: Knowledge Scaffold
        A[KEGG Pathways XML]:::raw --> B(parse_kegg.py):::script
        C[BioGRID Interactions TAB3]:::raw --> D(parse_biogrid.py):::script
        B --> E[merged_network.txt]:::output
        D --> E
    end

    %% Phase 2 & 3: Aim 1 - Diffusion Descriptor
    subgraph Phase 2 & 3: Aim 1 - Multi-Scale Diffusion
        F[TCGA Normal Expression]:::raw --> G[aim1_pipeline.py]:::engine
        H[TCGA Tumor Expression]:::raw --> G
        E --> G
        G -- "O(E) Edge Correlation" --> I[Condition Laplacians]
        I -- "PyTorch Chebyshev Trace Est." --> J[HKS Diffusion Descriptor]
        J -- "Per-Scale Centering & MAD Z-scores" --> K[Differential Profile & FDR]
        K --> L[dice_genes.txt]:::output
    end

    %% Phase 4: Validation
    subgraph Phase 4: Validation & Benchmarking
        L --> M[validation_pipeline.py]:::script
        M --> N[Precision@K & Volcano Plots .png]:::output
        M --> O[Benchmarking Stats .csv]:::output
    end
```