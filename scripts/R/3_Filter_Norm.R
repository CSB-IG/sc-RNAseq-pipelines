# ─────────────────────────────────────────────────────────────────────────────
# 1. Load Required Libraries
# ─────────────────────────────────────────────────────────────────────────────
set.seed(453)
library(Seurat)
library(ggplot2)
library(tidyverse)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Load Data
# ─────────────────────────────────────────────────────────────────────────────
seurat_object <- readRDS('seurat_object.rds')

# ─────────────────────────────────────────────────────────────────────────────
# 3. Quality Control (QC): Calculate Mitochondrial and Ribosomal Content
# ─────────────────────────────────────────────────────────────────────────────
seurat_object$percent.MT <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
seurat_object$percent.Ribosomal <- PercentageFeatureSet(seurat_object, pattern = "^RP[LS]")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Visualize QC Metrics
# ─────────────────────────────────────────────────────────────────────────────
options(repr.plot.height = 5, repr.plot.width = 25)
VlnPlot(
  seurat_object,
  features = c("nFeature_RNA", "nCount_RNA", "percent.MT"),
  ncol = 3,
  group.by = "individualID"
)

# ─────────────────────────────────────────────────────────────────────────────
# 5. Filter Cells Based on QC Metrics
# ─────────────────────────────────────────────────────────────────────────────
seurat_object_filtered <- subset(seurat_object, subset =
                                   nCount_RNA > 500 &
                                   nFeature_RNA > 200 &
                                   percent.MT < 30 &
                                   percent.Ribosomal < 5)

# ─────────────────────────────────────────────────────────────────────────────
# 6. Normalization and Feature Selection
# ─────────────────────────────────────────────────────────────────────────────
seurat_object_filtered <- NormalizeData(seurat_object_filtered)
seurat_object_filtered <- FindVariableFeatures(seurat_object_filtered)
seurat_object_filtered <- ScaleData(seurat_object_filtered)

# ─────────────────────────────────────────────────────────────────────────────
# 7. Save Filtered & Normalized Object
# ─────────────────────────────────────────────────────────────────────────────
saveRDS(seurat_object_filtered, file = "seurat_object_filtered.rds")
