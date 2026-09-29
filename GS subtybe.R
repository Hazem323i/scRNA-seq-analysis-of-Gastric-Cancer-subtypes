library(Seurat)
sc_data=Read10X("C:/Users/DELL/Downloads/10X_samples/GSM5573493_sample28.csv/")
sample28=CreateSeuratObject(project = "Project1",counts = sc_data,min.cells = 3,min.features = 200)
sample28[["percent_mt"]]=PercentageFeatureSet(sample28,pattern = "^MT-")
VlnPlot(sample28,features = c("nCount_RNA","nFeature_RNA","percent_mt"),ncol=3)
FeatureScatter(sample28,feature2="nCount_RNA",feature1 = "nFeature_RNA")

sample28=subset(sample28,
                subset = nCount_RNA > 200 
                & nCount_RNA< 5000
                & percent_mt< 10) 


VlnPlot(sample28,features = c("nCount_RNA","nFeature_RNA","percent_mt"),ncol=3)

library(scDblFinder)
scp28=as.SingleCellExperiment(sample28)
set.seed(298)
scp28= scDblFinder(scp28)

sample28$doublet_score=colData(scp28)$scDblFinder.score
sample28$doublet_class=colData(scp28)$scDblFinder.class
sample28=subset(sample28,subset = doublet_class =="singlet")
sample28=NormalizeData(sample28,normalization.method = "LogNormalize",scale.factor = 10000)
sample28=ScaleData(sample28)

sample28=FindVariableFeatures(sample28,selection.method = "vst",nfeatures = 2000)
sample28=ScaleData(sample28,features = VariableFeatures(sample28))
sample28=RunPCA(sample28,features = VariableFeatures(sample28))
DimPlot(sample28,reduction = "pca")

DimHeatmap(sample28,dims = 1:10,cells = 500,balanced = TRUE)


ElbowPlot(sample28)
sample28 = FindNeighbors(sample28, dims = 1:15)
sample28 = FindClusters(sample28, resolution = 0.5)
sample28 = RunUMAP(sample28, dims = 1:15)
DimPlot(sample28,reduction = "umap",label = TRUE)
lineage_markers <- c(
  "CD3D", "CD4", "CD8A", "NCAM1", "NKG7",       # T / NK
  "MS4A1", "IGHA1",                             # B / Plasma
  "CD14", "CD68", "FCGR3B", "LILRA4",           # Myeloid / Neutrophil / DC
  "EPCAM", "MUC5AC", "COL1A1", "PECAM1"         # Epithelial / Stromal
)
library(tidyverse)

DotPlot(sample28, features = lineage_markers) + RotatedAxis()

DimPlot(object = sample28,reduction = "umap",group.by = "cell_type",repel = T,label = T ) + ggtitle("sample28")
library(ggplot2)
markers_sample28=FindAllMarkers(sample28,min.pct = 0.25,only.pos = T,logfc.threshold = 0.25)
top_5_mark=markers_sample28 %>% group_by(cluster)%>% slice_max(order_by = avg_log2FC,n=5
)%>% arrange(cluster,desc(avg_log2FC))
DoHeatmap(sample28,features = unique(top_5_mark$gene),size = 3)
getwd()
write.csv(top_5_mark, "top5.csv")
write.csv(markers_sample28,"all.markers28.csv")
# Recode clusters to cell types for sample28
sample28$cell_type_DE = dplyr::recode(
  as.character(sample28$seurat_clusters),
  "0"  = "CD8+ Cytotoxic / Resident T Cells",
  "1"  = "Plasma Cells",
  "2"  = "CD4+ T Cells",
  "3"  = "B Cells",
  "4"  = "Epithelial Cells",
  "5"  = "Fibroblasts",
  "6"  = "NK Cells",
  "7"  = "Monocytes / Macrophages",
  "8"  = "Mast Cells",
  "9"  = "Endothelial Cells"
)

# Visualize the annotated UMAP
DimPlot(sample28, reduction = "umap", group.by = "cell_type_DE", repel = TRUE, label = TRUE) + 
  ggtitle("Sample28_Annotated")# Recode clusters to cell types for sample28
# Recode clusters to cell types for sample28
sample28$cell_type_DE = dplyr::recode(
  as.character(sample28$seurat_clusters),
  "0"  = "CD8+ Cytotoxic / Resident T Cells",
  "1"  = "Plasma Cells",
  "2"  = "CD4+ T Cells",
  "3"  = "B Cells",
  "4"  = "Epithelial Cells",
  "5"  = "Fibroblasts",
  "6"  = "NK Cells",
  "7"  = "Monocytes / Macrophages",
  "8"  = "Mast Cells",
  "9"  = "Endothelial Cells"
)

# Visualize the annotated UMAP
DimPlot(sample28, reduction = "umap", group.by = "cell_type_DE", repel = TRUE, label = TRUE) + 
  ggtitle("Sample28_Annotated")


#sample40
sc_data2=Read10X("C:/Users/DELL/Downloads/10X_samples/GSM5573505_sample40.csv/")
sample40=CreateSeuratObject(project = "Project2",counts = sc_data2,min.cells = 3,min.features = 200)
sample40[["percent_mt"]]=PercentageFeatureSet(sample40,pattern = "^MT-")
VlnPlot(sample40,features = c("nCount_RNA","nFeature_RNA","percent_mt"),ncol=3)
FeatureScatter(sample40,feature1="nCount_RNA",feature2 = "nFeature_RNA") 
sample40=subset(sample40,subset = nCount_RNA > 200 
                & nCount_RNA< 5000
                & percent_mt< 10)
scp40=as.SingleCellExperiment(sample40)
set.seed(123)
scp40=scDblFinder(scp40)
sample40$doublet_score=colData(scp40)$scDblFinder.score
sample40$doublet_class=colData(scp40)$scDblFinder.class
sample40=subset(sample40, subset=sample40$doublet_class == "singlet")
sample40=NormalizeData(sample40,normalization.method = "LogNormalize",scale.factor = 10000)
sample40=ScaleData(sample40)

sample40=FindVariableFeatures(sample40,selection.method = "vst",nfeatures = 2000)
sample40=ScaleData(sample40,features = VariableFeatures(sample40))
#integration
combined <- merge(sample28, sample40, add.cell.ids = c("s28", "s40"))


combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(combined))
DimPlot(combined, reduction = "pca", group.by = "orig.ident")

library('harmony')
combined <- RunHarmony(combined, group.by.vars = "orig.ident")
ElbowPlot(combined)
# 5. Cluster and UMAP on the harmony reduction
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:17)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:17)

DimPlot(combined, reduction = "umap", group.by = "orig.ident")  # check mixing
DimPlot(combined, reduction = "umap", label = TRUE)     


combined <- JoinLayers(combined)

markers_combined <- FindAllMarkers(combined, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(markers_combined,"C:/Users/DELL/Downloads/Combined_markersR.csv")


DotPlot(combined, features = lineage_markers) + RotatedAxis()
top_5_mark_combined=markers_combined %>% group_by(cluster)%>% slice_max(order_by = avg_log2FC,n=5
)%>% arrange(cluster,desc(avg_log2FC))
DoHeatmap(sample28,features = unique(top_5_mark$gene),size = 3)
#####
table(combined$orig.ident[combined$seurat_clusters == "15"])
VlnPlot(combined, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), idents = "15")
VlnPlot(combined, features = c("S100A8", "S100A9", "CSF3R", "FCGR3B"), idents = "15", split.by = "orig.ident")
####


# Recode clusters to cell types for the combined dataset
combined$cell_type_DE = dplyr::recode(
  as.character(combined$seurat_clusters),
  "0"  = "CD8+ T Cells / Resident Memory T Cells",
  "1"  = "Regulatory T Cells (Tregs)",
  "2"  = "Plasma Cells",
  "3"  = "Monocytes",
  "4"  = "B Cells",
  "5"  = "NK / γδ T Cells",
  "6"  = "Endothelial Cells",
  "7"  = "Pericytes / Stromal Cells",
  "8"  = "Epithelial Cells",
  "9"  = "Fibroblasts",
  "10" = "Mast Cells",
  "11" = "Neutrophils"
)

# Visualize the annotated UMAP
DimPlot(combined, reduction = "umap", group.by = "cell_type_DE", repel = TRUE, label = TRUE) + 
  ggtitle("Combined_Samples_GS_subtype")
#subset Macrophage
mac_cells <- subset(combined, subset = cell_type_DE == "Monocytes")
DefaultAssay(mac_cells) <- "RNA"
mac_cells <- FindVariableFeatures(mac_cells, nfeatures = 2000)
mac_cells <- ScaleData(mac_cells)
mac_cells <- RunPCA(mac_cells)
DimPlot(mac_cells,reduction = "pca",group.by = "orig.ident")
ElbowPlot(mac_cells)
mac_cells <- RunHarmony(mac_cells, group.by.vars = "orig.ident")
mac_cells <- FindNeighbors(mac_cells, reduction = "harmony", dims = 1:10)
mac_cells <- FindClusters(mac_cells, resolution = 0.8)
mac_cells <- RunUMAP(mac_cells, reduction = "harmony", dims = 1:10)
VlnPlot(mac_cells, features = c("SPP1", "C1QA", "C1QB", "C1QC"), group.by = "seurat_clusters")
DimPlot(mac_cells,reduction = "umap")

VlnPlot(mac_cells, features = c("FCN1", "LYZ", "CD14", "VCAN"), group.by = "seurat_clusters")   # classical/inflammatory monocytes
VlnPlot(mac_cells, features = c("FCER1A", "CD1C", "CLEC9A", "LAMP3"), group.by = "seurat_clusters")  # dendritic cells
VlnPlot(mac_cells, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), group.by = "seurat_clusters")  # rule out cluster 5 being low-quality, given how small/distinct it is
table(mac_cells$orig.ident[mac_cells$seurat_clusters == "4"])   # is it skewed to one sample?
VlnPlot(mac_cells, features = "doublet_score", group.by = "seurat_clusters")  # any residual signal?

mac_cells$mac_subtype <- dplyr::recode(
  as.character(mac_cells$seurat_clusters),
  "0" = "Mono_FCN1",             # Classical monocytes (FCN1+, VCAN+, CD14+)
  "1" = "SPP1_Mac_or_Mono",      # SPP1-expressing high-intensity subset / intermediate state
  "2" = "Macro_C1Q",             # C1QA/B/C-high tumor-associated macrophages
  "3" = "Dendritic_Cells"        # Dendritic cells (FCER1A+, LAMP3+)
)

FeaturePlot(mac_cells, features = c("SPP1", "C1QC"), reduction = "umap", label = TRUE)
# Check if SPP1 and C1QC co-express in the same cells
FeaturePlot(mac_cells, features = c("SPP1", "C1QC"), blend = TRUE, reduction = "umap")
# Sanity check the assignment
table(mac_cells$mac_subtype, mac_cells$seurat_clusters)

# Merge fine labels back into the main integrated object
combined$cell_type_fine <- as.character(combined$cell_type_DE)
combined$cell_type_fine[colnames(mac_cells)] <- mac_cells$mac_subtype
# Confirm the merge landed correctly
table(combined$cell_type_fine)
DimPlot(mac_cells,reduction = "umap",group.by = "mac_subtype",label = T)
DimPlot(combined,reduction = "umap",group.by = "cell_type_fine")


cd8_cells <- subset(combined, subset = cell_type_DE == "CD8+ T Cells / Resident Memory T Cells")
DefaultAssay(cd8_cells) <- "RNA"
cd8_cells <- FindVariableFeatures(cd8_cells, nfeatures = 2000)
cd8_cells <- ScaleData(cd8_cells)
cd8_cells <- RunPCA(cd8_cells, npcs = 20)
DimPlot(cd8_cells,group.by = "orig.ident",reduction = "pca")
ElbowPlot(object = cd8_cells)

cd8_cells <- RunHarmony(cd8_cells, group.by.vars = "orig.ident")
cd8_cells <- FindNeighbors(cd8_cells, reduction = "harmony", dims = 1:10)
cd8_cells <- FindClusters(cd8_cells, resolution = 0.8)
cd8_cells <- RunUMAP(cd8_cells, reduction = "harmony", dims = 1:10)
VlnPlot(cd8_cells, features = c("CCR7", "IL7R", "GZMK", "GZMB", "LAG3", "PDCD1", "HAVCR2", "TOX"), group.by = "seurat_clusters")
table(cd8_cells$seurat_clusters)   # cluster sizes — is cluster 0 tiny, making its violin unreliable?
VlnPlot(cd8_cells, features = c("TCF7", "SELL"), group.by = "seurat_clusters")   # distinguishes true naive from stem-like/progenitor-exhausted (Tpex), which can share some markers with both ends



cd8_cells$cd8_subtype <- dplyr::recode(
  as.character(cd8_cells$seurat_clusters),
  "0" = "CD8_Effector_Exhausted",
  "1" = "CD8_Advanced_Dysfunctional_Tim3", # Higher HAVCR2
  "2" = "CD8_Effector_Exhausted",
  "3" = "CD8_Effector_Exhausted",
  "4" = "CD8_TCF7_Progenitor_Exhausted",  # Stem-like Tpex
  "5" = "CD8_TOX_Terminally_Exhausted"    # Deep exhaustion
)

# Visualize the refined CD8 UMAP
DimPlot(cd8_cells, reduction = "umap", group.by = "cd8_subtype", label = TRUE, repel = TRUE) + 
  ggtitle("Refined CD8+ T-Cell States")

combined$cell_type_fine[colnames(cd8_cells)] <- cd8_cells$cd8_subtype
VlnPlot(combined, features = c("Exhaustion", "Cytotoxicity", "M2_macrophage"), group.by = "cell_type_fine") & Seurat::RotatedAxis()

# 1. Define gene sets
library(AUCell)
gene_sets <- list(
  Exhaustion    = c("LAG3", "PDCD1", "HAVCR2", "CTLA4", "TOX", "TIGIT"),
  Cytotoxicity  = c("GZMB", "GZMA", "PRF1", "NKG7", "GNLY"),
  M2_macrophage = c("SPP1", "C1QA", "C1QB", "C1QC", "APOE", "MRC1")
)

# 2. Build rankings and compute AUC scores
expr_matrix <- GetAssayData(combined, assay = "RNA", layer = "data")
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)

# 3. Attach scores to the object
combined <- AddMetaData(combined, metadata = as.data.frame(t(getAUC(cells_AUC))))

# CD8 T cell substates only — Exhaustion & Cytotoxicity
cd8_labels <- c(
  "CD8_Effector_Exhausted",
  "CD8_Advanced_Dysfunctional_Tim3",
  "CD8_TCF7_Progenitor_Exhausted",
  "CD8_TOX_Terminally_Exhausted"
)
VlnPlot(subset(combined, subset = cell_type_fine %in% cd8_labels),
        features = c("Exhaustion", "Cytotoxicity"),
        group.by = "cell_type_fine") & Seurat::RotatedAxis()

# Macrophage/monocyte substates 
mac_labels <- c(
  "Mono_FCN1",
  "SPP1_Mac_or_Mono",
  "Macro_C1Q",
  "Dendritic_Cells"
)


VlnPlot(subset(combined, subset = cell_type_fine %in% mac_labels),
        features = c("APOE", "MRC1", "SPP1"),
        group.by = "cell_type_fine") & Seurat::RotatedAxis()

##########