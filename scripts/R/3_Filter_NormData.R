# ─────────────────────────────────────────────────────────────────────────────
# 1. Load Required Libraries
# ─────────────────────────────────────────────────────────────────────────────
set.seed(453)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DESeq2)
library(SingleR)
library(gridExtra)
library(rtracklayer)
library(infercnv)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Load Data
# ─────────────────────────────────────────────────────────────────────────────
# Load a list of Seurat objects
conteo <- readRDS('seurat_object.rds')

# ─────────────────────────────────────────────────────────────────────────────
# 3. Inspect Data Structure
# ─────────────────────────────────────────────────────────────────────────────
typeof(conteo)        # Check the object type
length(conteo)        # Number of Seurat objects (samples)
typeof(conteo[[1]])   # Type of the first element (should be Seurat)

# ─────────────────────────────────────────────────────────────────────────────
# 4. Merge All Seurat Objects into One
# ─────────────────────────────────────────────────────────────────────────────
lista_nombres <- names(conteo)
merged <- merge(
  x = conteo[[1]],
  y = conteo[2:length(conteo)],
  add.cell.ids = lista_nombres,
  project = "sc"
)

# ─────────────────────────────────────────────────────────────────────────────
# 5. Quality Control (QC): Calculate Mitochondrial and Ribosomal Content
# ─────────────────────────────────────────────────────────────────────────────
merged$percent.MT <- PercentageFeatureSet(merged, pattern="^MT-")           # Mitochondrial genes
merged$percent.Ribosomal <- PercentageFeatureSet(merged, pattern="^RP[LS]") # Ribosomal genes

# ─────────────────────────────────────────────────────────────────────────────
# 6. Visualize QC Metrics
# ─────────────────────────────────────────────────────────────────────────────
options(repr.plot.height=5, repr.plot.width=25)
VlnPlot(
  merged,
  features = c("nFeature_RNA", "nCount_RNA", "percent.MT"),
  ncol = 3,
  group.by = "individualID"
)

# ─────────────────────────────────────────────────────────────────────────────
# 7. Filter Cells Based on QC Metrics
# ─────────────────────────────────────────────────────────────────────────────
merge_filter <- subset(merged, subset =
                         nCount_RNA > 500 &
                         nFeature_RNA > 200 &
                         percent.MT < 30 &
                         percent.Ribosomal < 5)

# ─────────────────────────────────────────────────────────────────────────────
# 8. Normalization and Feature Selection
# ─────────────────────────────────────────────────────────────────────────────
merge_filter <- NormalizeData(merge_filter)
merge_filter <- FindVariableFeatures(merge_filter)
merge_filter <- ScaleData(merge_filter)

# ─────────────────────────────────────────────────────────────────────────────
# 9. Principal Component Analysis (PCA)
# ─────────────────────────────────────────────────────────────────────────────
merge_filter <- RunPCA(merge_filter)
ElbowPlot(merge_filter)

# ─────────────────────────────────────────────────────────────────────────────
# 10. Determine Optimal Number of PCs to Use
# ─────────────────────────────────────────────────────────────────────────────
pct <- merge_filter[["pca"]]@stdev / sum(merge_filter[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
pcs <- min(co1, co2)

# ─────────────────────────────────────────────────────────────────────────────
# 11. Clustering and Dimensionality Reduction
# ─────────────────────────────────────────────────────────────────────────────
merge_filter <- FindNeighbors(object = merge_filter, dims = 1:pcs)
merge_filter <- FindClusters(object = merge_filter)
merge_filter <- RunUMAP(object = merge_filter, dims = 1:pcs)
merge_filter <- RunTSNE(object = merge_filter, dims = 1:pcs)

# ─────────────────────────────────────────────────────────────────────────────
# 12. Visualize Clusters
# ─────────────────────────────────────────────────────────────────────────────
options(repr.plot.width = 8, repr.plot.height = 8)
DimPlot(merge_filter, label = TRUE) +
  ggtitle("UMAP of 50,142 Brain Cells") +
  labs(subtitle = "Unintegrated Samples")

# ─────────────────────────────────────────────────────────────────────────────
# 13. CCA-Based Data Integration
# ─────────────────────────────────────────────────────────────────────────────
integration <- IntegrateLayers(
  object = merge_filter,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

# Merge raw RNA layers
integration[["RNA"]] <- JoinLayers(integration[["RNA"]])

# ─────────────────────────────────────────────────────────────────────────────
# 14. Post-Integration Clustering
# ─────────────────────────────────────────────────────────────────────────────
integration <- FindNeighbors(integration, reduction = "integrated.cca", dims = 1:18)
integration <- FindClusters(integration, resolution = 0.5)

# ─────────────────────────────────────────────────────────────────────────────
# 15. Visualize Integrated Clustering
# ─────────────────────────────────────────────────────────────────────────────
DimPlot(integration, label = TRUE) +
  ggtitle("UMAP of 50,142 Brain Cells") +
  labs(subtitle = "CCA Integration")

# ─────────────────────────────────────────────────────────────────────────────
# 16. Define Marker Genes for Cell Types (Edit Manually)
# ─────────────────────────────────────────────────────────────────────────────
markers <- list()
markers[['']] <- c('')
markers[['']] <- c('')
markers[['']] <- c('')
markers[['']] <- c('')

# ─────────────────────────────────────────────────────────────────────────────
# 17. Visualize Marker Gene Expression
# ─────────────────────────────────────────────────────────────────────────────
options(repr.plot.height=8, repr.plot.width=20)
DotPlot(integration, features = markers)

# ─────────────────────────────────────────────────────────────────────────────
# 18. Find Differentially Expressed Genes (Markers)
# ─────────────────────────────────────────────────────────────────────────────
markers <- FindAllMarkers(integration, only.pos = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 19. Save Processed Objects
# ─────────────────────────────────────────────────────────────────────────────
saveRDS(merged, "/object_merged.rds")
saveRDS(integration, "/seurat_integ.rds")

# ─────────────────────────────────────────────────────────────────────────────
# 20. Prepare Object for Cell Type Annotation with SingleR
# ─────────────────────────────────────────────────────────────────────────────
neurocounts <- as.SingleCellExperiment(DietSeurat(integration))
saveRDS(neurocounts, file = "counts.rds")
