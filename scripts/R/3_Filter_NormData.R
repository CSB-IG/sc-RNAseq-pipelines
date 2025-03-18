if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(Seurat, dplyr, ggplot2, patchwork, Matrix, data.table, stringr, cowplot, ggpubr)

# Set the working directory
setwd("~/Prueba_Get_Data/")

# Download and extract data
system("wget https://cf.10xgenomics.com/samples/cell-exp/1.1.0/aml035_pre_transplant/aml035_pre_transplant_raw_gene_bc_matrices.tar.gz")
system("tar -xvzf aml035_pre_transplant_raw_gene_bc_matrices.tar.gz")

# Load extracted count matrix
gene_matrix <- readMM("matrices_mex/hg19/matrix.mtx")
gene_names <- fread("matrices_mex/hg19/genes.tsv", header = FALSE)
barcode_names <- fread("matrices_mex/hg19/barcodes.tsv", header = FALSE)

# Ensure unique row names
gene_names$V2 <- make.unique(gene_names$V2)

# Assign row and column names
rownames(gene_matrix) <- gene_names$V2
colnames(gene_matrix) <- barcode_names$V1

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = gene_matrix)

# Normalize and process data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Save processed Seurat object
saveRDS(seurat_obj, file = "~/Prueba_Get_Data/processed_seurat_obj.rds")

# Save filtered count matrix
writeMM(seurat_obj@assays$RNA@counts, file = "~/Prueba_Get_Data/filtered_matrix.mtx")
write.table(rownames(seurat_obj), file = "~/Prueba_Get_Data/filtered_genes.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(seurat_obj), file = "~/Prueba_Get_Data/filtered_barcodes.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)
