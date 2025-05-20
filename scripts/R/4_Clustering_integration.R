# ─────────────────────────────────────────────────────────────────────────────
# 1. Load Required Libraries
# ─────────────────────────────────────────────────────────────────────────────
set.seed(453)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(SingleR)
library(SingleCellExperiment)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Load Preprocessed Seurat Object
# ─────────────────────────────────────────────────────────────────────────────
seurat_object <- readRDS("seurat_object_filtered.rds")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Principal Component Analysis (PCA)
# ─────────────────────────────────────────────────────────────────────────────
seurat_object <- RunPCA(seurat_object)
ElbowPlot(seurat_object)

# ─────────────────────────────────────────────────────────────────────────────
# 4. Determine Optimal Number of PCs to Use
# ─────────────────────────────────────────────────────────────────────────────
pct <- seurat_object[["pca"]]@stdev / sum(seurat_object[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
pcs <- min(co1, co2)

# ─────────────────────────────────────────────────────────────────────────────
# 5. Clustering and Dimensionality Reduction
# ─────────────────────────────────────────────────────────────────────────────
seurat_object <- FindNeighbors(seurat_object, dims = 1:pcs)
seurat_object <- FindClusters(seurat_object)
seurat_object <- RunUMAP(seurat_object, dims = 1:pcs)
seurat_object <- RunTSNE(seurat_object, dims = 1:pcs)

# ─────────────────────────────────────────────────────────────────────────────
# 6. Visualize Clusters
# ─────────────────────────────────────────────────────────────────────────────
options(repr.plot.width = 8, repr.plot.height = 8)
DimPlot(seurat_object, label = TRUE) +
  ggtitle("UMAP of Filtered Cells") +
  labs(subtitle = "Unintegrated Samples")

# ─────────────────────────────────────────────────────────────────────────────
# 7. CCA-Based Data Integration
# ─────────────────────────────────────────────────────────────────────────────
integration <- IntegrateLayers(
  object = seurat_object,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)
integration[["RNA"]] <- JoinLayers(integration[["RNA"]])

# ─────────────────────────────────────────────────────────────────────────────
# 8. Post-Integration Clustering
# ─────────────────────────────────────────────────────────────────────────────
integration <- FindNeighbors(integration, reduction = "integrated.cca", dims = 1:18)
integration <- FindClusters(integration, resolution = 0.5)

# ─────────────────────────────────────────────────────────────────────────────
# 9. Visualize Integrated Clustering
# ─────────────────────────────────────────────────────────────────────────────
DimPlot(integration, label = TRUE) +
  ggtitle("UMAP of Integrated Cells") +
  labs(subtitle = "CCA Integration")

# ─────────────────────────────────────────────────────────────────────────────
# 10. Save Results
# ─────────────────────────────────────────────────────────────────────────────
saveRDS(seurat_object, "seurat_object_clustered.rds")
saveRDS(integration, "seurat_object_integrated.rds")

