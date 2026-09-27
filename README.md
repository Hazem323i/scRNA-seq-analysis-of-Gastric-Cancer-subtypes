# scRNA-seq-analysis-of-Gastric-Cancer-subtypes
Gastric Cancer scRNA-seq Analysis (MSI Subtype)

Single-cell RNA-seq analysis of a gastric adenocarcinoma sample belonging to the MSI (microsatellite instability-high) molecular subtype, characterizing the tumor microenvironment's cell-type composition using a standard Seurat workflow in R. This is part of a group project (6 contributors); this repository covers the per-sample QC → clustering → annotation pipeline for one sample, with integration of additional MSI-subtype samples planned as a next step.

Data
Sample: GSM5573480 (10x Genomics single-cell RNA-seq, gastric cancer, MSI subtype)
Format: standard 10x Genomics output (barcodes / features / matrix)
(Add the parent GEO series accession and study citation here once confirmed.)
Pipeline

All analysis is in R using Seurat.

Load data — Read10X() → CreateSeuratObject(min.cells = 3, min.features = 200)
QC — mitochondrial content via PercentageFeatureSet(pattern = "^MT-"); cells filtered to nCount_RNA between 200–5500 and percent_mt < 10%
Doublet removal — scDblFinder run on a SingleCellExperiment conversion of the object; doublet scores/classes transferred back to the Seurat object, singlets retained
Normalization — NormalizeData() (LogNormalize, scale factor 10,000)
Feature selection & scaling — FindVariableFeatures() (vst, top 2,000 genes) → ScaleData() on variable features
Dimensionality reduction — RunPCA() on variable features; PCs inspected via DimHeatmap() and ElbowPlot()
Clustering — FindNeighbors() / FindClusters() (dims 1:17, resolution 0.5), UMAP via RunUMAP()
Cell-type annotation — three complementary approaches, cross-checked against each other:
Canonical markers — DotPlot() across T/NK, B/Plasma, Myeloid/DC, and Epithelial/Stromal marker panels → cell_type
Data-driven markers — FindAllMarkers() per cluster, top hits validated against known cell-type markers → cell_type_DE
Reference-based annotation — SingleR + celldex (in progress)
Key Findings
14 transcriptionally distinct clusters identified across the sample
Cross-checking canonical-marker calls against each cluster's top DE genes corrected several initial labels — most notably uncovering a mast cell population and a regulatory T cell (Treg) population that canonical markers alone had mislabeled
One cluster driven entirely by cell-cycle genes (no lineage markers) was further resolved via targeted marker panels (PTPRC, EPCAM, MKI67, then CD3D, CD8A, FOXP3, NKG7) as proliferating cytotoxic CD8+ T / NK cells — consistent with the strong tumor-reactive T-cell expansion expected in MSI-high disease
A rich immune infiltrate was identified alongside tumor/normal epithelial, endothelial, and fibroblast compartments
Repository Contents
File	Description
analysis.R (rename as appropriate)	Full analysis script
all_markers.csv	Full FindAllMarkers() output, all clusters
top5.csv	Top 5 marker genes per cluster
Dependencies
R ≥ 4.5
Seurat
scDblFinder
SingleR, celldex
tidyverse (dplyr, ggplot2)
Next Steps
Complete SingleR/celldex reference-based annotation
Integrate additional MSI-subtype samples (multi-sample integration, e.g. Seurat anchor-based integration or Harmony)
Copy-number variation inference (inferCNV / CopyKAT) on epithelial clusters to distinguish malignant from normal gastric epithelium
How to Run
Update the file path in Read10X() to point to your local 10x data folder
Run the script top to bottom in R
Outputs: marker tables (all_markers.csv, top5.csv) and diagnostic/result plots (QC violins, PCA, elbow plot, UMAP, DotPlot, DoHeatmap)
Contributors

Group project — 6 contributors 
Hazem Sayed
Bassem Elsayed 
Neama ALaa-eldin
Retag Mohamed 
Iman Manar 
Yossef Mostafa
