// DiCE-Duan Analyzer - Mock Data
// Realistic bioinformatics data for cancer network analysis

// Top 10 Rewired Genes with ensemble scores
export const topRewiredGenes = [
  {
    rank: 1,
    gene: "ZRANB1",
    ensembleScore: 0.94,
    egatProb: 0.94,
    embeddingDrift: 0.8,
    sksDisruption: 0.68,
    deltaBetweenness: 0.21,
    deltaEigenvector: 0.11,
    normalBetweenness: 0.03,
    normalEigenvector: 0.05,
    tumorBetweenness: 0.24,
    tumorEigenvector: 0.16,
    category: "top",
    description: "Zinc finger RANBP2-type containing 1 - Shows significant rewiring in tumor networks"
  },
  {
    rank: 2,
    gene: "TP53",
    ensembleScore: 0.91,
    egatProb: 0.91,
    embeddingDrift: 0.77,
    sksDisruption: 0.66,
    deltaBetweenness: 0.18,
    deltaEigenvector: 0.09,
    normalBetweenness: 0.12,
    normalEigenvector: 0.14,
    tumorBetweenness: 0.30,
    tumorEigenvector: 0.23,
    category: "top",
    description: "Tumor protein p53 - Master regulator frequently mutated in cancer"
  },
  {
    rank: 3,
    gene: "EGFR",
    ensembleScore: 0.88,
    egatProb: 0.88,
    embeddingDrift: 0.74,
    sksDisruption: 0.64,
    deltaBetweenness: 0.16,
    deltaEigenvector: 0.08,
    normalBetweenness: 0.08,
    normalEigenvector: 0.11,
    tumorBetweenness: 0.24,
    tumorEigenvector: 0.19,
    category: "top",
    description: "Epidermal growth factor receptor - Key oncogene in multiple cancer types"
  },
  {
    rank: 4,
    gene: "MYC",
    ensembleScore: 0.85,
    egatProb: 0.85,
    embeddingDrift: 0.71,
    sksDisruption: 0.62,
    deltaBetweenness: 0.15,
    deltaEigenvector: 0.07,
    normalBetweenness: 0.15,
    normalEigenvector: 0.18,
    tumorBetweenness: 0.30,
    tumorEigenvector: 0.25,
    category: "top",
    description: "MYC proto-oncogene - Transcription factor driving cell proliferation"
  },
  {
    rank: 5,
    gene: "KRAS",
    ensembleScore: 0.82,
    egatProb: 0.82,
    embeddingDrift: 0.69,
    sksDisruption: 0.59,
    deltaBetweenness: 0.14,
    deltaEigenvector: 0.06,
    normalBetweenness: 0.09,
    normalEigenvector: 0.10,
    tumorBetweenness: 0.23,
    tumorEigenvector: 0.16,
    category: "top",
    description: "KRAS proto-oncogene - Frequently mutated in pancreatic and colorectal cancers"
  },
  {
    rank: 6,
    gene: "PTEN",
    ensembleScore: 0.79,
    egatProb: 0.79,
    embeddingDrift: 0.66,
    sksDisruption: 0.57,
    deltaBetweenness: 0.13,
    deltaEigenvector: 0.06,
    normalBetweenness: 0.11,
    normalEigenvector: 0.13,
    tumorBetweenness: 0.24,
    tumorEigenvector: 0.19,
    category: "candidate",
    description: "Phosphatase and tensin homolog - Tumor suppressor regulating PI3K pathway"
  },
  {
    rank: 7,
    gene: "BRCA1",
    ensembleScore: 0.76,
    egatProb: 0.76,
    embeddingDrift: 0.63,
    sksDisruption: 0.55,
    deltaBetweenness: 0.12,
    deltaEigenvector: 0.05,
    normalBetweenness: 0.07,
    normalEigenvector: 0.09,
    tumorBetweenness: 0.19,
    tumorEigenvector: 0.14,
    category: "candidate",
    description: "Breast cancer 1 - DNA repair protein associated with hereditary cancer"
  },
  {
    rank: 8,
    gene: "PIK3CA",
    ensembleScore: 0.74,
    egatProb: 0.74,
    embeddingDrift: 0.62,
    sksDisruption: 0.53,
    deltaBetweenness: 0.11,
    deltaEigenvector: 0.05,
    normalBetweenness: 0.06,
    normalEigenvector: 0.08,
    tumorBetweenness: 0.17,
    tumorEigenvector: 0.13,
    category: "candidate",
    description: "Phosphatidylinositol-4,5-bisphosphate 3-kinase - Key signaling enzyme"
  },
  {
    rank: 9,
    gene: "BRAF",
    ensembleScore: 0.71,
    egatProb: 0.71,
    embeddingDrift: 0.59,
    sksDisruption: 0.51,
    deltaBetweenness: 0.10,
    deltaEigenvector: 0.04,
    normalBetweenness: 0.05,
    normalEigenvector: 0.07,
    tumorBetweenness: 0.15,
    tumorEigenvector: 0.11,
    category: "candidate",
    description: "B-Raf proto-oncogene - Serine/threonine kinase in MAPK pathway"
  },
  {
    rank: 10,
    gene: "CDH1",
    ensembleScore: 0.68,
    egatProb: 0.68,
    embeddingDrift: 0.56,
    sksDisruption: 0.49,
    deltaBetweenness: 0.09,
    deltaEigenvector: 0.04,
    normalBetweenness: 0.08,
    normalEigenvector: 0.10,
    tumorBetweenness: 0.17,
    tumorEigenvector: 0.14,
    category: "candidate",
    description: "Cadherin 1 - Cell adhesion molecule, loss associated with metastasis"
  }
];

// Network nodes (genes)
export const networkNodes = [
  // Top rewired genes (larger nodes)
  { id: "ZRANB1", label: "ZRANB1", group: "top", size: 24 },
  { id: "TP53", label: "TP53", group: "top", size: 28 },
  { id: "EGFR", label: "EGFR", group: "top", size: 26 },
  { id: "MYC", label: "MYC", group: "top", size: 25 },
  { id: "KRAS", label: "KRAS", group: "top", size: 22 },
  { id: "PTEN", label: "PTEN", group: "candidate", size: 20 },
  { id: "BRCA1", label: "BRCA1", group: "candidate", size: 20 },
  { id: "PIK3CA", label: "PIK3CA", group: "candidate", size: 18 },
  { id: "BRAF", label: "BRAF", group: "candidate", size: 18 },
  { id: "CDH1", label: "CDH1", group: "candidate", size: 16 },
  
  // Additional network genes
  { id: "AKT1", label: "AKT1", group: "normal", size: 14 },
  { id: "MTOR", label: "MTOR", group: "normal", size: 14 },
  { id: "RB1", label: "RB1", group: "normal", size: 14 },
  { id: "ATM", label: "ATM", group: "normal", size: 12 },
  { id: "CHEK2", label: "CHEK2", group: "normal", size: 12 },
  { id: "MLH1", label: "MLH1", group: "normal", size: 12 },
  { id: "MSH2", label: "MSH2", group: "normal", size: 12 },
  { id: "ERBB2", label: "ERBB2", group: "normal", size: 14 },
  { id: "FGFR1", label: "FGFR1", group: "normal", size: 12 },
  { id: "VEGFA", label: "VEGFA", group: "normal", size: 14 },
  { id: "STAT3", label: "STAT3", group: "normal", size: 13 },
  { id: "JAK2", label: "JAK2", group: "normal", size: 12 },
  { id: "SRC", label: "SRC", group: "normal", size: 13 },
  { id: "CTNNB1", label: "CTNNB1", group: "normal", size: 14 },
  { id: "SMAD4", label: "SMAD4", group: "normal", size: 12 },
  { id: "TGFBR1", label: "TGFBR1", group: "normal", size: 11 },
  { id: "NOTCH1", label: "NOTCH1", group: "normal", size: 13 },
  { id: "WNT1", label: "WNT1", group: "normal", size: 11 },
  { id: "HIF1A", label: "HIF1A", group: "normal", size: 12 },
  { id: "FOXO3", label: "FOXO3", group: "normal", size: 11 },
  { id: "NFKB1", label: "NFKB1", group: "normal", size: 13 },
  { id: "RELA", label: "RELA", group: "normal", size: 12 },
  { id: "APC", label: "APC", group: "normal", size: 14 },
  { id: "CDKN2A", label: "CDKN2A", group: "normal", size: 13 },
  { id: "MDM2", label: "MDM2", group: "normal", size: 12 },
  { id: "CDK4", label: "CDK4", group: "normal", size: 11 },
  { id: "E2F1", label: "E2F1", group: "normal", size: 11 },
  { id: "RBL2", label: "RBL2", group: "normal", size: 10 },
  { id: "BCL2", label: "BCL2", group: "normal", size: 12 },
  { id: "CASP9", label: "CASP9", group: "normal", size: 10 },
  { id: "PARP1", label: "PARP1", group: "normal", size: 11 },
  { id: "ATR", label: "ATR", group: "normal", size: 11 },
  { id: "BRCA2", label: "BRCA2", group: "normal", size: 13 },
  { id: "RAD51", label: "RAD51", group: "normal", size: 11 },
  { id: "MSH6", label: "MSH6", group: "normal", size: 11 },
  { id: "PMS2", label: "PMS2", group: "normal", size: 10 },
  { id: "EPCAM", label: "EPCAM", group: "normal", size: 10 },
  { id: "STK11", label: "STK11", group: "normal", size: 11 },
  // Missing nodes added here:
  { id: "BAX", label: "BAX", group: "normal", size: 12 },
  { id: "MAPK1", label: "MAPK1", group: "normal", size: 13 }
];

// Network edges (interactions)
export const networkEdges = [
  // TP53 hub connections
  { source: "TP53", target: "MDM2", weight: 0.9, attention: 0.88, directed: true },
  { source: "TP53", target: "CDKN2A", weight: 0.85, attention: 0.83, directed: true },
  { source: "TP53", target: "BAX", weight: 0.8, attention: 0.78, directed: true },
  { source: "TP53", target: "ATM", weight: 0.75, attention: 0.73, directed: true },
  { source: "TP53", target: "CHEK2", weight: 0.7, attention: 0.68, directed: true },
  
  // EGFR signaling
  { source: "EGFR", target: "PIK3CA", weight: 0.88, attention: 0.86, directed: true },
  { source: "EGFR", target: "KRAS", weight: 0.82, attention: 0.8, directed: true },
  { source: "EGFR", target: "STAT3", weight: 0.78, attention: 0.76, directed: true },
  { source: "EGFR", target: "SRC", weight: 0.75, attention: 0.73, directed: true },
  { source: "EGFR", target: "ERBB2", weight: 0.85, attention: 0.83, directed: true },
  
  // PI3K/AKT pathway
  { source: "PIK3CA", target: "AKT1", weight: 0.92, attention: 0.89, directed: true },
  { source: "PIK3CA", target: "PTEN", weight: 0.88, attention: 0.86, directed: true },
  { source: "AKT1", target: "MTOR", weight: 0.86, attention: 0.84, directed: true },
  { source: "AKT1", target: "FOXO3", weight: 0.72, attention: 0.7, directed: true },
  { source: "PTEN", target: "MTOR", weight: 0.68, attention: 0.67, directed: true },
  
  // KRAS/MAPK pathway
  { source: "KRAS", target: "BRAF", weight: 0.9, attention: 0.88, directed: true },
  { source: "BRAF", target: "MAPK1", weight: 0.85, attention: 0.83, directed: true },
  { source: "KRAS", target: "PIK3CA", weight: 0.76, attention: 0.74, directed: true },
  
  // MYC network
  { source: "MYC", target: "CDK4", weight: 0.82, attention: 0.8, directed: true },
  { source: "MYC", target: "E2F1", weight: 0.78, attention: 0.76, directed: true },
  { source: "MYC", target: "BCL2", weight: 0.7, attention: 0.68, directed: true },
  
  // Cell cycle
  { source: "CDK4", target: "RB1", weight: 0.88, attention: 0.86, directed: true },
  { source: "RB1", target: "E2F1", weight: 0.85, attention: 0.83, directed: true },
  { source: "CDKN2A", target: "CDK4", weight: 0.9, attention: 0.88, directed: true },
  
  // DNA repair
  { source: "BRCA1", target: "BRCA2", weight: 0.87, attention: 0.85, directed: true },
  { source: "BRCA1", target: "RAD51", weight: 0.82, attention: 0.8, directed: true },
  { source: "BRCA1", target: "ATM", weight: 0.78, attention: 0.76, directed: true },
  { source: "BRCA2", target: "RAD51", weight: 0.85, attention: 0.83, directed: true },
  
  // TGF-beta
  { source: "TGFBR1", target: "SMAD4", weight: 0.8, attention: 0.78, directed: true },
  { source: "SMAD4", target: "CTNNB1", weight: 0.68, attention: 0.67, directed: true },
  
  // Wnt signaling
  { source: "WNT1", target: "CTNNB1", weight: 0.82, attention: 0.8, directed: true },
  { source: "CTNNB1", target: "APC", weight: 0.88, attention: 0.86, directed: true },
  { source: "APC", target: "CDH1", weight: 0.72, attention: 0.7, directed: true },
  
  // Notch
  { source: "NOTCH1", target: "MYC", weight: 0.7, attention: 0.68, directed: true },
  
  // Hypoxia
  { source: "HIF1A", target: "VEGFA", weight: 0.86, attention: 0.84, directed: true },
  { source: "HIF1A", target: "MTOR", weight: 0.74, attention: 0.72, directed: true },
  
  // NF-kB
  { source: "NFKB1", target: "RELA", weight: 0.88, attention: 0.86, directed: true },
  { source: "NFKB1", target: "STAT3", weight: 0.72, attention: 0.7, directed: true },
  
  // Apoptosis
  { source: "BCL2", target: "CASP9", weight: 0.78, attention: 0.76, directed: true },
  { source: "PARP1", target: "CASP9", weight: 0.7, attention: 0.68, directed: true },
  
  // MMR
  { source: "MLH1", target: "MSH2", weight: 0.9, attention: 0.88, directed: true },
  { source: "MLH1", target: "MSH6", weight: 0.85, attention: 0.83, directed: true },
  { source: "MSH2", target: "PMS2", weight: 0.82, attention: 0.8, directed: true },
  
  // Additional connections for density
  { source: "ZRANB1", target: "TP53", weight: 0.65, attention: 0.64, directed: true },
  { source: "ZRANB1", target: "MYC", weight: 0.62, attention: 0.61, directed: true },
  { source: "ZRANB1", target: "EGFR", weight: 0.58, attention: 0.57, directed: true },
  { source: "JAK2", target: "STAT3", weight: 0.84, attention: 0.82, directed: true },
  { source: "FGFR1", target: "PIK3CA", weight: 0.76, attention: 0.74, directed: true },
  { source: "STK11", target: "MTOR", weight: 0.8, attention: 0.78, directed: true },
  { source: "EPCAM", target: "CDH1", weight: 0.74, attention: 0.72, directed: true },
  { source: "RBL2", target: "E2F1", weight: 0.68, attention: 0.67, directed: true },
  { source: "ATR", target: "ATM", weight: 0.82, attention: 0.8, directed: true },
  { source: "ATR", target: "CHEK2", weight: 0.76, attention: 0.74, directed: true }
];

// Stats for the dashboard
export const dashboardStats = {
  nodesAnalyzed: 2180,
  edgesProcessed: 12400,
  topCandidate: "ZRANB1",
  tissueTypes: ["Normal", "Tumor"],
  analysisDate: "2024-01-15"
};

// Scatter plot data (centrality comparison)
export const scatterData = [
  // Top rewired genes (gold)
  { gene: "ZRANB1", deltaBetweenness: 0.21, deltaEigenvector: 0.11, category: "top", ensembleScore: 0.94 },
  { gene: "TP53", deltaBetweenness: 0.18, deltaEigenvector: 0.09, category: "top", ensembleScore: 0.91 },
  { gene: "EGFR", deltaBetweenness: 0.16, deltaEigenvector: 0.08, category: "top", ensembleScore: 0.88 },
  { gene: "MYC", deltaBetweenness: 0.15, deltaEigenvector: 0.07, category: "top", ensembleScore: 0.85 },
  { gene: "KRAS", deltaBetweenness: 0.14, deltaEigenvector: 0.06, category: "top", ensembleScore: 0.82 },
  
  // Candidates
  { gene: "PTEN", deltaBetweenness: 0.13, deltaEigenvector: 0.06, category: "candidate", ensembleScore: 0.79 },
  { gene: "BRCA1", deltaBetweenness: 0.12, deltaEigenvector: 0.05, category: "candidate", ensembleScore: 0.76 },
  { gene: "PIK3CA", deltaBetweenness: 0.11, deltaEigenvector: 0.05, category: "candidate", ensembleScore: 0.74 },
  { gene: "BRAF", deltaBetweenness: 0.10, deltaEigenvector: 0.04, category: "candidate", ensembleScore: 0.71 },
  { gene: "CDH1", deltaBetweenness: 0.09, deltaEigenvector: 0.04, category: "candidate", ensembleScore: 0.68 },
  
  // Other genes (filler for visualization)
  { gene: "AKT1", deltaBetweenness: 0.08, deltaEigenvector: 0.035, category: "other", ensembleScore: 0.62 },
  { gene: "MTOR", deltaBetweenness: 0.075, deltaEigenvector: 0.032, category: "other", ensembleScore: 0.58 },
  { gene: "RB1", deltaBetweenness: 0.07, deltaEigenvector: 0.03, category: "other", ensembleScore: 0.55 },
  { gene: "ATM", deltaBetweenness: 0.065, deltaEigenvector: 0.028, category: "other", ensembleScore: 0.52 },
  { gene: "CHEK2", deltaBetweenness: 0.06, deltaEigenvector: 0.025, category: "other", ensembleScore: 0.48 },
  { gene: "MLH1", deltaBetweenness: 0.055, deltaEigenvector: 0.023, category: "other", ensembleScore: 0.45 },
  { gene: "MSH2", deltaBetweenness: 0.05, deltaEigenvector: 0.02, category: "other", ensembleScore: 0.42 },
  { gene: "ERBB2", deltaBetweenness: 0.048, deltaEigenvector: 0.019, category: "other", ensembleScore: 0.40 },
  { gene: "FGFR1", deltaBetweenness: 0.045, deltaEigenvector: 0.018, category: "other", ensembleScore: 0.38 },
  { gene: "VEGFA", deltaBetweenness: 0.042, deltaEigenvector: 0.017, category: "other", ensembleScore: 0.36 },
  { gene: "STAT3", deltaBetweenness: 0.04, deltaEigenvector: 0.016, category: "other", ensembleScore: 0.34 },
  { gene: "JAK2", deltaBetweenness: 0.038, deltaEigenvector: 0.015, category: "other", ensembleScore: 0.32 },
  { gene: "SRC", deltaBetweenness: 0.035, deltaEigenvector: 0.014, category: "other", ensembleScore: 0.30 },
  { gene: "CTNNB1", deltaBetweenness: 0.032, deltaEigenvector: 0.013, category: "other", ensembleScore: 0.28 },
  { gene: "SMAD4", deltaBetweenness: 0.03, deltaEigenvector: 0.012, category: "other", ensembleScore: 0.26 },
  { gene: "TGFBR1", deltaBetweenness: 0.028, deltaEigenvector: 0.011, category: "other", ensembleScore: 0.24 },
  { gene: "NOTCH1", deltaBetweenness: 0.025, deltaEigenvector: 0.01, category: "other", ensembleScore: 0.22 },
  { gene: "WNT1", deltaBetweenness: 0.022, deltaEigenvector: 0.009, category: "other", ensembleScore: 0.20 },
  { gene: "HIF1A", deltaBetweenness: 0.02, deltaEigenvector: 0.008, category: "other", ensembleScore: 0.18 },
  { gene: "FOXO3", deltaBetweenness: 0.018, deltaEigenvector: 0.007, category: "other", ensembleScore: 0.16 },
  { gene: "NFKB1", deltaBetweenness: 0.015, deltaEigenvector: 0.006, category: "other", ensembleScore: 0.14 },
  { gene: "RELA", deltaBetweenness: 0.012, deltaEigenvector: 0.005, category: "other", ensembleScore: 0.12 },
  { gene: "APC", deltaBetweenness: 0.01, deltaEigenvector: 0.004, category: "other", ensembleScore: 0.10 },
  { gene: "CDKN2A", deltaBetweenness: 0.008, deltaEigenvector: 0.003, category: "other", ensembleScore: 0.08 },
  { gene: "MDM2", deltaBetweenness: 0.005, deltaEigenvector: 0.002, category: "other", ensembleScore: 0.06 }
];

// Methodology steps
export const methodologySteps = [
  {
    id: 1,
    title: "Build Networks",
    description: "Construct Normal and Tumor co-expression networks from BioGRID/KEGG with Pearson-distance edge weights.",
    icon: "network"
  },
  {
    id: 2,
    title: "EA-GAT Embedding",
    description: "Edge-augmented GATv2 learns anisotropic attention over heterophilic biological graphs.",
    icon: "calculator"
  },
  {
    id: 3,
    title: "Embedding Drift (Δh)",
    description: "Contrastive self-supervised loss ranks genes by normal→tumor latent displacement.",
    icon: "git-merge"
  },
  {
    id: 4,
    title: "Rank & Export",
    description: "Generate ranked gene lists with EA-GAT probability, Δh, and SKS disruption scores.",
    icon: "download"
  }
];

// Helper function to get gene details
export const getGeneDetails = (geneName) => {
  return topRewiredGenes.find(g => g.gene === geneName) || null;
};

// Helper function to get node color based on group
export const getNodeColor = (group, tissue = "normal") => {
  const colors = {
    top: tissue === "normal" ? "#00F0FF" : "#FF2D8D",
    candidate: tissue === "normal" ? "#00F0FF" : "#FF2D8D",
    normal: tissue === "normal" ? "#00F0FF" : "#FF2D8D",
    selected: "#FFD166"
  };
  return colors[group] || colors.normal;
};
