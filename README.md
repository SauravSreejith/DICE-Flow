# DiCE-Flow: Structural Response and Representation Learning for Cancer Driver Gene Discovery

Traditional disease gene discovery relies heavily on Differential Gene Expression (DEA), which flags genes whose transcription levels spike or plummet. However, DEA is structurally blind—it fails to detect "silent" cancer driver genes (like *TP53*, *KRAS*, or *SMARCA4*) whose overall expression remains stable, but whose structural role within the protein-protein interaction (PPI) network rewires massively to support the tumor.

**DiCE-Flow** is a highly scalable, GPU-accelerated computational pipeline that identifies these silent drivers by measuring their **Multi-Scale Structural Repositioning** between Normal and Disease states. 

It replaces standard, single-scale graph centralities with a **Heat Kernel Signature (HKS) Diffusion Descriptor**, evaluated via PyTorch Chebyshev trace estimation, and guarded by rigorous Degree-Stratified Robust False Discovery Rate (FDR) control.

---

## 🌟 Key Features

- **Multi-Scale Resolution:** Probes both immediate local complexes (small $t$) and global pathway positioning (large $t$).
- **Hardware-Agnostic Acceleration:** Automatically scales across NVIDIA CUDA, Apple Silicon (MPS), or multi-core CPUs.
- **Massively Scalable:** $O(E)$ memory-optimized Pearson correlation and batched Hutchinson trace estimators process 17,500+ nodes and >1M edges in minutes.
- **Rigorous Statistical Control:** Degree-stratified robust FDR (using Median Absolute Deviation) cleanly separates true biological rewiring from stochastic trace-estimator noise and cohort batch effects.
- **Automated Validation:** Built-in benchmarking against CRISPR-Cas9 essentiality screens (DepMap), Literature (PubMed), and baseline methods.

---

## ⚙️ Installation

### 1. Requirements
* **Python:** 3.9+
* **Hardware:** 16GB+ System RAM. (NVIDIA GPU or Apple M-Series chip highly recommended for 10x speedup).

### 2. Setup Environment
```bash
git clone https://github.com/your-repo/dice-flow.git
cd dice-flow

# Create a virtual environment
python3 -m venv env
source env/bin/activate

# Install dependencies (Make sure to install the correct PyTorch version for your hardware)
pip install -r requirements.txt
```
*(Note: If you have an NVIDIA GPU, ensure you install `torch` with CUDA support via the official PyTorch website).*

---

## 📥 Data Acquisition

Place the required datasets into the `data/raw/` directories before running the pipeline.

### 1. KEGG Pathways (`data/raw/kegg/`)
Download Human (`hsa`) KGML (XML) pathway maps from the [KEGG Pathway Database](https://www.genome.jp/kegg/pathway.html).

### 2. BioGRID Interactions (`data/raw/biogrid/`)
Download the latest **TAB3** format file from [BioGRID](https://downloads.thebiogrid.org/BioGRID/Release-Archive/).
* *Target File:* `BIOGRID-ALL-5.0.262.tab3.txt`

### 3. Expression Matrices (`data/raw/expression/`)
You need two Tab-Separated Value (TSV) matrices: `expression_normal.txt` and `expression_tumor.txt`.
* **Rows:** Gene Symbols (e.g., `TP53`, `MYC`).
* **Columns:** Patient Samples (e.g., TCGA Barcodes).
* **Values:** Log2 Normalized Expression (e.g., $\log_2(\text{TPM} + 1)$).
* *Recommended Source:* [UCSC Xena Browser (TCGA Pan-Cancer)](https://xenabrowser.net/datapages/).

---

## 🚀 Usage

### 1. Execute the Main Pipeline (Aim 1)
The master bash script will automatically parse the databases, build the PPI scaffold, compute $O(E)$ correlations, fire up the GPU Trace Estimator, run the null permutations, and generate the final discoveries.

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```

### 2. Generate Validation & Benchmarking Reports
Once the pipeline finishes, run the automated validation suite to mathematically prove the results and generate publication-ready plots.

```bash
python src/python/validation/validation_pipeline.py --all --skip-pubmed
```
*(Figures will be saved to `data/results/validation/figures/`)*.

---

## 📊 Interpreting the Results

The primary output is `data/results/dice_genes_top500.txt`. 

| Metric | Description |
| :--- | :--- |
| **`ensemble_score`** ($\Delta$) | The absolute magnitude of the structural repositioning (using the $L_\infty$ maximum absolute shift across scales). **Higher is better.** |
| **`q_value_fdr`** | The Benjamini-Hochberg adjusted statistical significance (Robust MAD). **Values < 0.05 are significant discoveries.** |
| **`characteristic_scale`** ($t^*$) | The specific diffusion time where the gene shifted the most. Small values = local complex rewiring. Large values = global pathway shifts. |
| **`direction`** | `-1` = Gene localized (heat concentrated). `<br>` `+1` = Gene dispersed (acts as a newly formed topological bridge). |

**A successful run looks like this:**
Genes like *TP53*, *KRAS*, *SMARCA4*, and *EGFR* will rise to the Top 10 with $q < 10^{-20}$. The validation pipeline will show a **Fitness Gene Overlap > 18%**, massively outperforming the Differential Expression baseline.

---

## 📈 Performance

Benchmarks calculated on a 17,466-node, 2,000,000-edge human interactome with $S=256$ trace estimation vectors and $B=50$ null permutations:

* **NVIDIA RTX 4050 / RTX 3060:** ~3.0 - 4.0 Minutes.
* **Apple Silicon (M2/M3 Pro):** ~4.5 - 6.0 Minutes.
* **Multi-Core CPU (16 Cores):** ~15 - 20 Minutes.

---

## 📖 Citation & References

If you use DiCE-Flow in your research, please refer to the foundational methodologies:

1. **DiCE (Base Premise):** Pashaei E. et al. (2025). DiCE: Differential Centrality-Ensemble Analysis. *Nucleic Acids Research*.
2. **Heat Kernel Signatures:** Sun J. et al. (2009). A Concise and Provably Informative Multi-Scale Signature Based on Heat Diffusion. *Computer Graphics Forum*.
3. **Trace Estimation:** Hutchinson, M. F. (1990). A stochastic estimator of the trace of the influence matrix. *Communications in Statistics-Simulation and Computation*.

## 📄 License
MIT License - See LICENSE file for details.