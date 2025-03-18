library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(data.table)
library(stringr)
library(cowplot)
library(ggpubr)

# List filtered .h5 files
data_files <- list.files("../clean_adata/", pattern = "filtered.h5$")

data_files

# Function to load data
load_it <- function(file) {
  samp <- strsplit(file, "_")[[1]][1]  # Extract patient/sample identifier
  dx <- strsplit(file, "_")[[1]][2]  # Extract diagnosis code
  seurat_obj <- Read10X_h5(paste0("../clean_adata/", file)) %>%
    CreateSeuratObject()  # Load the dataset and create a Seurat object
  seurat_obj$Patient <- samp  # Assign patient ID
  seurat_obj$DX <- dx  # Assign diagnosis code
  seurat_obj$Sample <- paste0(samp, "_", dx)  # Create a unique sample identifier
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = paste0(samp, "_", dx))  # Rename cells for uniqueness
  return(seurat_obj)
}

# Load data from all files
datasets <- lapply(data_files, load_it)

# Quality control function
qc <- function(seurat_obj) {
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200)  # Filter out low-quality cells
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)  # Normalize expression data
  return(seurat_obj)
}

# Apply QC to all datasets
datasets <- lapply(datasets, qc)

# Merge all datasets into a single Seurat object
seurat_combined <- merge(datasets[[1]], y = datasets[-1])

# Identify mitochondrial gene percentage
seurat_combined["percent.mt"] <- PercentageFeatureSet(seurat_combined, pattern = "^MT-")

# Filter out cells with high mitochondrial content
seurat_combined <- subset(seurat_combined, subset = percent.mt < 25)

# Identify and remove doublets
seurat_combined <- SCTransform(seurat_combined)  # Normalize using SCTransform
seurat_combined <- RunPCA(seurat_combined)  # Perform PCA for dimensionality reduction
seurat_combined <- RunUMAP(seurat_combined, dims = 1:10)  # Perform UMAP for visualization

# Save final Seurat object
saveRDS(seurat_combined, file = "~/doublets_removed.rds")

# Save filtered count matrix
writeMM(seurat_combined@assays$RNA@counts, file = "~/filtered_counts.mtx")
