# scRNA-seq Analysis of Gastric Cancer Subtypes

Team project (6 members) investigating whether a published immunosuppressive stromal barrier — SPP1⁺/C1QC⁺ macrophages signaling to exhausted CD8⁺ T cells (Ma et al., 2025, *Frontiers in Immunology*) — is present across three gastric cancer (GC) molecular subtypes: Genomically Stable (GS), Chromosomal Instability (CIN), and Microsatellite Instability (MSI), excluding EBV⁺ tumors.

## My Contribution

While this is a team project, **the analysis in this repository — the complete pipeline from QC through functional scoring, for both the MSI and GS subtypes — was performed entirely by me.** MSI was my originally assigned subtype; I additionally carried out a full, independent re-analysis of the GS subtype. Together these cover 4 of the project's 6 patient samples. The CIN subtype was analyzed separately by other team members; their results are referenced below for the team's cross-subtype comparison but were not produced or independently reproduced in this repository.

| Subtype | Samples | Analysis by |
|---|---|---|
| MSI | sample15 (Stage II), sample33 (Stage III) | **Hazem Sayed Abdelmoaty** (this repo) |
| GS | sample28 (Stage III), sample40 (Stage III) | **Hazem Sayed Abdelmoaty** (this repo) |
| CIN | sample34 (Stage III), sample39 (Stage III) | Other team members |

## Research Question

Is the SPP1⁺/C1QC⁺ macrophage–CD8⁺ Tex immunosuppressive barrier present across the GS, CIN, and MSI molecular subtypes of gastric cancer, or is it subtype-selective — as the source paper's own stratification would predict?

## Pipeline (R / Seurat v5)

Applied independently to each of the four samples I analyzed:

1. QC filtering (`nCount_RNA`, `nFeature_RNA`, mitochondrial %) and doublet removal (scDblFinder)
2. Normalization, variable feature selection, and scaling
3. Sample-pair integration via Harmony; joint clustering (Louvain) and UMAP
4. Cell-type annotation via canonical marker panels, cross-validated against per-cluster differential expression (`FindAllMarkers`)
5. Targeted sub-clustering of the macrophage/monocyte and CD8⁺ T cell compartments to resolve fine-grained states
6. AUCell functional scoring (exhaustion, cytotoxicity, M2-macrophage gene signatures)
7. CellChat ligand–receptor analysis, testing macrophage → CD8⁺ Tex signaling directly

## Key Findings

**MSI (my analysis):** No SPP1⁺ macrophage population detected across any check performed. The closest macrophage subset (C1Q-high, APOE-high) showed an atypical, partial M2 signature with inconsistent MRC1 detection. A distinct, AUCell-confirmed terminally exhausted CD8⁺ T cell population (HAVCR2⁺, LAG3-high, PDCD1⁺) was identified alongside naive, GZMK⁺ effector, and TCF7⁺/SELL⁺ progenitor-exhausted (Tpex) states.

**GS (my analysis):** A distinct SPP1⁺ macrophage population was identified, spatially separated from a C1QC⁺ population with no double-positive overlap — the clean two-population architecture the source paper describes. The CD8⁺ compartment resolved into four exhaustion-associated states (Dysfunctional_Tim3, Effector_Exhausted, TCF7_Progenitor_Exhausted, TOX_Terminally_Exhausted).

**CIN (teammates' analysis, referenced for comparison):** A reproducible SPP1⁺/TREM2⁺ lipid-associated macrophage population and a distinct exhausted/tumor-reactive CD8⁺ cluster were reported, consistent across both patient samples.

**Synthesis:** GS and CIN both show the SPP1⁺ macrophage component of the barrier; MSI does not — consistent with the source paper's own characterization of MSI⁺ tumors as exhaustion-dominant rather than barrier-dominant.

## Tools

R, RStudio, Seurat (v5), scDblFinder, Harmony, AUCell, CellChat, tidyverse (dplyr, ggplot2), Git/GitHub

## Repository Contents

- `gastric_cancer_MSI_analysis.R` — full MSI subtype analysis script
- GS subtype analysis script
- Marker gene output tables (`top5.csv`, `all.markers.csv`, and GS equivalents)

## Team

Neama Alaa eldin Ali, Hazem Sayed Abdelmoaty, Retag Mohamed Amer, Youssef Mostafa Mohamed, Bassem El-Sayed Ahmed, Iman Manar El-Islam

## Reference

Ma G, Liu X, Jiang Q, et al. (2025) Identification of a stromal immunosuppressive barrier orchestrated by SPP1⁺/C1QC⁺ macrophages and CD8⁺ exhausted T cells driving gastric cancer immunotherapy resistance. *Front. Immunol.* 16:1618591. doi: 10.3389/fimmu.2025.1618591