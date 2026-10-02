# =============================================================================
# Single-Cell RNA-Seq Analysis of Gastric Cancer Subtypes — MSI Subtype
# Samples: sample15 (Stage II) and sample33 (Stage III)
# Author:  Hazem Sayed Abdelmoaty
# Part of a 6-member team project testing a published SPP1+/C1QC+ macrophage -
# CD8+ exhausted T cell immunosuppressive barrier model (Ma et al. 2025,
# Front. Immunol.) across the GS, CIN, and MSI molecular subtypes of
# gastric cancer.
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(Seurat)
library(scDblFinder)
library(tidyverse)
library(ggplot2)
# install.packages("harmony")
library(harmony)
library(AUCell)


# =============================================================================
# SAMPLE 15 (MSI, Cancer Stage II)
# =============================================================================

# ---- Load data and compute QC metrics ---------------------------------------
sc_data <- Read10X("C:/Users/DELL/Downloads/10X_samples/GSM5573480_sample15.csv/")
sample15 <- CreateSeuratObject(project = "Project1", counts = sc_data, min.cells = 3, min.features = 200)
sample15[["percent_mt"]] <- PercentageFeatureSet(sample15, pattern = "^MT-")

VlnPlot(sample15, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)
FeatureScatter(sample15, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# ---- QC filtering -------------------------------------------------------
sample15 <- subset(sample15,
                   subset = nCount_RNA > 200
                   & nCount_RNA < 5500
                   & percent_mt < 10)

VlnPlot(sample15, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)  # post-filter check

# ---- Doublet detection (scDblFinder) -----------------------------------------
scp15 <- as.SingleCellExperiment(sample15)
set.seed(223)
scp15 <- scDblFinder(scp15)

sample15$doublet_score <- colData(scp15)$scDblFinder.score
sample15$doublet_class <- colData(scp15)$scDblFinder.class
sample15 <- subset(sample15, subset = doublet_class == "singlet")

# ---- Normalization, variable features, scaling -------------------------------
sample15 <- NormalizeData(sample15, normalization.method = "LogNormalize", scale.factor = 10000)
sample15 <- ScaleData(sample15)

sample15 <- FindVariableFeatures(sample15, selection.method = "vst", nfeatures = 2000)
sample15 <- ScaleData(sample15, features = VariableFeatures(sample15))

# ---- PCA and dimensionality selection -----------------------------------------
sample15 <- RunPCA(sample15, features = VariableFeatures(sample15))
DimPlot(sample15, reduction = "pca")
DimHeatmap(sample15, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(sample15)

# ---- Clustering and UMAP -----------------------------------------------------
sample15 <- FindNeighbors(sample15, dims = 1:17)
sample15 <- FindClusters(sample15, resolution = 0.5)
sample15 <- RunUMAP(sample15, dims = 1:17)
DimPlot(sample15, reduction = "umap", label = TRUE)

# ---- Cell-type annotation: canonical lineage markers --------------------------
lineage_markers <- c(
  "CD3D", "CD4", "CD8A", "NCAM1", "NKG7",       # T / NK
  "MS4A1", "IGHA1",                             # B / Plasma
  "CD14", "CD68", "FCGR3B", "LILRA4",           # Myeloid / Neutrophil / DC
  "EPCAM", "MUC5AC", "COL1A1", "PECAM1"         # Epithelial / Stromal
)

DotPlot(sample15, features = lineage_markers) + RotatedAxis()

sample15$cell_type <- dplyr::recode(
  as.character(sample15$seurat_clusters),
  "0"  = "CD8+ T Cells",
  "1"  = "CD8+ T / NK Cells",
  "2"  = "T Cells",
  "3"  = "NK Cells",
  "4"  = "Monocytes / Macrophages",
  "5"  = "B Cells",
  "6"  = "Myeloid Cells",
  "7"  = "CD4+ T Cells",
  "8"  = "Endothelial Cells",
  "9"  = "Plasma Cells",
  "10" = "Mucosal Epithelial Cells",
  "11" = "Cytotoxic T / NK Cells",
  "12" = "Epithelial Cells",
  "13" = "Fibroblasts"
)
DimPlot(object = sample15, reduction = "umap", group.by = "cell_type", repel = TRUE, label = TRUE) + ggtitle("sample15")

# ---- Cell-type annotation: validated against per-cluster DE markers -----------
markers_sample15 <- FindAllMarkers(sample15, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
top_5_mark <- markers_sample15 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  arrange(cluster, desc(avg_log2FC))

DoHeatmap(sample15, features = unique(top_5_mark$gene), size = 3)

write.csv(top_5_mark, "top5.csv")
write.csv(markers_sample15, "all.markers.csv")

sample15$cell_type_DE <- dplyr::recode(
  as.character(sample15$seurat_clusters),
  "0"  = "CD8+ T Cells",              # CD8A/B, GZMH, GZMB, NKG7 — matches
  "1"  = "CD8+ T / NK Cells",         # CD8A/B, GZMK, GZMM — mostly CD8 T signal, see note
  "2"  = "CD4+ T Cells",              # CD4, CD28, ICOS, CTLA4, no FOXP3 — confirmed
  "3"  = "NK Cells",                  # GNLY, KLRC1, KIR2DL4, GZMA — matches, some γδ T signal
  "4"  = "Monocytes / Macrophages",   # CD14, C1QA/B/C, CD68, APOE — matches
  "5"  = "B Cells",                   # CD19, MS4A1, CD79A/B — matches
  "6"  = "Mast Cells",                # CPA3, TPSAB1, KIT, MS4A2 — strongly confirmed
  "7"  = "Regulatory T Cells (Treg)", # FOXP3, IL2RA, CTLA4, TNFRSF18 — strongly confirmed
  "8"  = "Endothelial Cells",         # CDH5, VWF, CLDN5, PLVAP — matches
  "9"  = "Plasma Cells",              # IGHG/IGL/IGK genes, MZB1, JCHAIN — matches
  "10" = "Mucosal Epithelial Cells",  # DUOX2, GKN1, GKN2, MUC4 — matches (gastric-specific)
  "11" = "Proliferating Cells",       # MKI67, TOP2A, CDK1 — no lineage genes present, see note
  "12" = "Epithelial Cells",          # EPCAM, TFF3, AGR2 — matches
  "13" = "Fibroblasts"                # PDGFRA, LUM, POSTN, PDPN — matches
)
DimPlot(sample15, reduction = "umap", group.by = "cell_type_DE", repel = TRUE, label = TRUE) + ggtitle("Sample15_DE")


# =============================================================================
# SAMPLE 33 (MSI, cancer Stage III)
# =============================================================================

# ---- Load data, QC, and filter -----------------------------------------------
sc_data2 <- Read10X("C:/Users/DELL/Downloads/10X_samples/GSM5573498_sample33.csv/")
sample33 <- CreateSeuratObject(project = "Project2", counts = sc_data2, min.cells = 3, min.features = 200)
sample33[["percent_mt"]] <- PercentageFeatureSet(sample33, pattern = "^MT-")

VlnPlot(sample33, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)
FeatureScatter(sample33, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

sample33 <- subset(sample33, subset = nCount_RNA > 200
                   & nCount_RNA < 4500
                   & percent_mt < 10)

# ---- Doublet detection (scDblFinder) -----------------------------------------
scp33 <- as.SingleCellExperiment(sample33)
set.seed(123)
scp33 <- scDblFinder(scp33)

sample33$doublet_score <- colData(scp33)$scDblFinder.score
sample33$doublet_class <- colData(scp33)$scDblFinder.class
sample33 <- subset(sample33, subset = sample33$doublet_class == "singlet")

# ---- Normalization, variable features, scaling -------------------------------
sample33 <- NormalizeData(sample33, normalization.method = "LogNormalize", scale.factor = 10000)
sample33 <- ScaleData(sample33)

sample33 <- FindVariableFeatures(sample33, selection.method = "vst", nfeatures = 2000)
sample33 <- ScaleData(sample33, features = VariableFeatures(sample33))


# =============================================================================
# INTEGRATION (Harmony) — sample15 + sample33
# =============================================================================

combined <- merge(sample15, sample33, add.cell.ids = c("s15", "s33"))

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(combined))
DimPlot(combined, reduction = "pca", group.by = "orig.ident")   # pre-integration batch check

combined <- RunHarmony(combined, group.by.vars = "orig.ident")
ElbowPlot(combined)

# ---- Clustering and UMAP on the Harmony-corrected embedding -------------------
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:17)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:17)

DimPlot(combined, reduction = "umap", group.by = "orig.ident")  # check sample mixing
DimPlot(combined, reduction = "umap", label = TRUE)

# ---- Join layers before differential expression -------------------------
combined <- JoinLayers(combined)

markers_combined <- FindAllMarkers(combined, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(markers_combined, "C:/Users/DELL/Downloads/Combined_markers.csv")

DotPlot(combined, features = lineage_markers) + RotatedAxis()

top_5_mark_combined <- markers_combined %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  arrange(cluster, desc(avg_log2FC))

# ---- Cell-type annotation: validated against per-cluster DE markers -----------
combined$cell_type_DE <- dplyr::recode(
  as.character(combined$seurat_clusters),
  "0"  = "CD8+ Exhausted / Resident T Cells", # CD8A, CD8B, GZMH, GZMB, LAG3, ITGAE, CXCR6
  "1"  = "CD8+ Effector Memory T Cells",      # GZMK, CD8A, CD8B, IL7R, CXCR3
  "2"  = "Regulatory T Cells",                # FOXP3, IL2RA, CTLA4, TNFRSF4/18, TIGIT, CD4
  "3"  = "CD4+ Naive / Memory T Cells",       # CCR7, IL7R, KLF2, RORA, KLRB1, CCR6
  "4"  = "B Cells",                           # CD19, MS4A1, CD79A/B, BANK1
  "5"  = "Gastric Epithelial Cells",          # CLDN18, CLDN4, ELF3, KLF5
  "6"  = "NK / γδ T Cells",                   # GNLY, KLRC1, KIR2DL4, SH2D1B, TRDC
  "7"  = "Monocytes / Macrophages",           # CD14, CD68, C1QA/B/C, APOE, CSF1R
  "8"  = "CD16+ Cytotoxic NK Cells",          # FCGR3A, KLRF1, KLRD1, GNLY, PRF1
  "9"  = "Mast Cells",                        # CPA3, TPSAB1, TPSB2, KIT, HDC, MS4A2
  "10" = "Plasma Cells",                      # IGHG genes, JCHAIN, MZB1, TNFRSF17
  "11" = "Endothelial Cells",                 # CDH5, VWF, CLDN5, PLVAP, CD34
  "12" = "Mucous Epithelial Cells",           # TFF3, AGR2, KRT18, FCGBP, DEFB1
  "13" = "Proliferating Cells",               # MKI67, TOP2A, CDK1, RRM2, UBE2C
  "14" = "Fibroblasts",                       # PDGFRA, LUM, POSTN, PDPN, DCN
  "15" = "Neutrophils"                        # S100A8, S100A9, CSF3R, FCGR3B, CXCL8
)
DimPlot(combined, reduction = "umap", group.by = "cell_type_DE", repel = TRUE, label = TRUE) + ggtitle("combined_samples_DE")


# =============================================================================
# MACROPHAGE / MONOCYTE SUB-CLUSTERING
# =============================================================================

mac_cells <- subset(combined, subset = cell_type_DE == "Monocytes / Macrophages")
DefaultAssay(mac_cells) <- "RNA"
mac_cells <- FindVariableFeatures(mac_cells, nfeatures = 2000)
mac_cells <- ScaleData(mac_cells)
mac_cells <- RunPCA(mac_cells, npcs = 20)
DimPlot(mac_cells, reduction = "pca", group.by = "orig.ident")
ElbowPlot(mac_cells)

mac_cells <- RunHarmony(mac_cells, group.by.vars = "orig.ident")
mac_cells <- FindNeighbors(mac_cells, reduction = "harmony", dims = 1:6)
mac_cells <- FindClusters(mac_cells, resolution = 0.8)
mac_cells <- RunUMAP(mac_cells, reduction = "harmony", dims = 1:6)

# ---- Marker checks: SPP1/C1Q macrophages, FCN1 monocytes, dendritic cells ----
VlnPlot(mac_cells, features = c("SPP1", "C1QA", "C1QB", "C1QC"), group.by = "seurat_clusters")
DimPlot(mac_cells, reduction = "umap")

VlnPlot(mac_cells, features = c("FCN1", "LYZ", "CD14", "VCAN"), group.by = "seurat_clusters")         # classical/inflammatory monocytes
VlnPlot(mac_cells, features = c("FCER1A", "CD1C", "CLEC9A", "LAMP3"), group.by = "seurat_clusters")   # dendritic cells

# ---- Fine-grained macrophage/monocyte subtype labels --------------------------
mac_cells$mac_subtype <- dplyr::recode(as.character(mac_cells$seurat_clusters),
                                       "0" = "Mono_FCN1_CD14lo",     # FCN1+, VCAN+, CD14-low — monocyte-lineage, atypical CD14 status
                                       "1" = "Macro_C1Q",            # C1QA/B/C-high, CD14+, SPP1-negative
                                       "2" = "Macro_C1Q",            # C1QA/B/C-high, CD14+, SPP1-negative
                                       "3" = "Macro_C1Q",            # C1QA/B/C-high, CD14+, SPP1-negative
                                       "4" = "Macro_C1Q_Quiescent",  # C1Q-high but low nCount/nFeature, no doublet/mito signal — distinct lower-activity state
                                       "5" = "Mono_FCN1"             # FCN1+, VCAN+, CD14+ — matches paper's Mono_FCN1 (responder-associated)
)


# ---- Merge fine labels back into the main integrated object -------------------
combined$cell_type_fine <- as.character(combined$cell_type_DE)
combined$cell_type_fine[colnames(mac_cells)] <- mac_cells$mac_subtype
table(combined$cell_type_fine)   # confirm the merge landed correctly

DimPlot(mac_cells, reduction = "umap", group.by = "mac_subtype", label = TRUE)
DimPlot(combined, reduction = "umap", group.by = "cell_type_fine")


# =============================================================================
# CD8+ T CELL SUB-CLUSTERING
# =============================================================================

cd8_cells <- subset(combined, subset = cell_type_DE %in% c("CD8+ Exhausted / Resident T Cells", "CD8+ Effector Memory T Cells"))
DefaultAssay(cd8_cells) <- "RNA"
cd8_cells <- FindVariableFeatures(cd8_cells, nfeatures = 2000)
cd8_cells <- ScaleData(cd8_cells)
cd8_cells <- RunPCA(cd8_cells, npcs = 20)
DimPlot(cd8_cells, group.by = "orig.ident", reduction = "pca")
ElbowPlot(object = cd8_cells)

cd8_cells <- RunHarmony(cd8_cells, group.by.vars = "orig.ident")
cd8_cells <- FindNeighbors(cd8_cells, reduction = "harmony", dims = 1:9)
cd8_cells <- FindClusters(cd8_cells, resolution = 0.8)
cd8_cells <- RunUMAP(cd8_cells, reduction = "harmony", dims = 1:9)

# ---- Marker checks: naive / effector / exhaustion panel -----------------------
VlnPlot(cd8_cells, features = c("CCR7", "IL7R", "GZMK", "GZMB", "LAG3", "PDCD1", "HAVCR2", "TOX"), group.by = "seurat_clusters")
VlnPlot(cd8_cells, features = c("TCF7", "SELL"), group.by = "seurat_clusters")   # distinguishes true naive from stem-like/progenitor-exhausted, which can share some markers with both ends

# ---- QC robustness checks on sub-clusters -------------------------------------
VlnPlot(cd8_cells, features = c("nCount_RNA", "nFeature_RNA"), group.by = "seurat_clusters")
VlnPlot(cd8_cells, features = c("percent_mt", "doublet_score"), group.by = "seurat_clusters")
table(cd8_cells$orig.ident[cd8_cells$seurat_clusters == "3"])

# ---- Fine-grained CD8+ T cell subtype labels -----------------------------------
cd8_cells$cd8_subtype <- dplyr::recode(as.character(cd8_cells$seurat_clusters),
                                       "0" = "CD8_Tpex_stemlike",   # progenitor-exhausted, not memory — tumor = chronic antigen
                                       "1" = "CD8_Teff_GZMK",
                                       "2" = "CD8_Tnaive",
                                       "3" = "CD8_Teff_GZMK",       # low-depth but real — confirmed via mito%/doublet score
                                       "4" = "CD8_Teff_GZMK",
                                       "5" = "CD8_Teff_GZMK",
                                       "6" = "CD8_Tex_terminal",    # terminally exhausted — HAVCR2+/LAG3-high/PDCD1+
                                       "7" = "CD8_Teff_GZMK"
)

DimPlot(cd8_cells, reduction = "umap", group.by = "cd8_subtype", label = TRUE)

# ---- Merge fine labels back into the main integrated object -------------------
combined$cell_type_fine[colnames(cd8_cells)] <- cd8_cells$cd8_subtype


# =============================================================================
# AUCELL FUNCTIONAL SCORING (Exhaustion / Cytotoxicity / M2-macrophage)
# =============================================================================

gene_sets <- list(
  Exhaustion    = c("LAG3", "PDCD1", "HAVCR2", "CTLA4", "TOX", "TIGIT"),
  Cytotoxicity  = c("GZMB", "GZMA", "PRF1", "NKG7", "GNLY"),
  M2_macrophage = c("SPP1", "C1QA", "C1QB", "C1QC", "APOE", "MRC1")
)

expr_matrix <- GetAssayData(combined, assay = "RNA", layer = "data")
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)

combined <- AddMetaData(combined, metadata = as.data.frame(t(getAUC(cells_AUC))))

VlnPlot(combined, features = c("Exhaustion", "Cytotoxicity", "M2_macrophage"), group.by = "cell_type_fine") & Seurat::RotatedAxis()

# ---- CD8 substates only — Exhaustion & Cytotoxicity ---------------------------
cd8_labels <- c("CD8_Tnaive", "CD8_Teff_GZMK", "CD8_Tpex_stemlike", "CD8_Tex_terminal")
VlnPlot(subset(combined, subset = cell_type_fine %in% cd8_labels),
        features = c("Exhaustion", "Cytotoxicity"),
        group.by = "cell_type_fine") & Seurat::RotatedAxis()

# ---- Macrophage/monocyte substates only — M2_macrophage, then individual genes
mac_labels <- c("Mono_FCN1_CD14lo", "Macro_C1Q", "Macro_C1Q_Quiescent", "Mono_FCN1")
VlnPlot(subset(combined, subset = cell_type_fine %in% mac_labels),
        features = "M2_macrophage",
        group.by = "cell_type_fine") & Seurat::RotatedAxis()

VlnPlot(subset(combined, subset = cell_type_fine %in% mac_labels),
        features = c("APOE", "MRC1", "SPP1"),
        group.by = "cell_type_fine") & Seurat::RotatedAxis()