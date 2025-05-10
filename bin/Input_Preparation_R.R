library(Seurat)
library(SeuratDisk)
library(data.table)

# Take positional arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Error: Two positional arguments required: <seurat_obj_dir> <output_dir>")
}

seurat_obj_dir <- args[1]  # The directory of the .rds file
output_dir <- args[2]  # The output directory

# Load the Seurat object
seurat_obj <- readRDS(seurat_obj_dir)

# apply log-normalization 
if (sum(GetAssayData(seurat_obj, layer = "data", assay = "RNA")) == 0) {
  message("\nLog-normalization has not been run. Running it now...")
  seurat_obj <- NormalizeData(seurat_obj)
} else {
  message("\nLog-normalization has already been run.")
}

# select top 2000 variance genes
if (length(VariableFeatures(seurat_obj)) == 0) {
  # Run FindVariableFeatures if it hasn't been run yet
  message("Performing high variance genes selection ...\n")
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
} else {
  message("High variance genes already selected.\n")
}

# Attempt to convert the Seurat object to h5ad and save it to the output directory
h5seurat_dir <- file.path(output_dir, "data.h5seurat")
h5ad_dir <- file.path(output_dir, "data.h5ad")

conversion_result <- tryCatch({
  SaveH5Seurat(seurat_obj, filename = h5seurat_dir)  # Convert to .h5seurat
  Convert(h5seurat_dir, dest = "h5ad")  # Convert to .h5ad
  file.remove(h5seurat_dir)
  list(success = TRUE)
}, error = function(e) {
  cat("Error during h5ad conversion: ", e$message, "\n")
  list(success = FALSE)
})

# If the conversion failed, remove the .h5ad and .h5seurat files
if (!conversion_result$success) {
  cat("Conversion to h5ad failed. Removing .h5ad and .h5seurat files...\n")
  if (file.exists(h5ad_dir)) {
    file.remove(h5ad_dir)
  }
  if (file.exists(h5seurat_dir)) {
    file.remove(h5seurat_dir)
  }
  
  # Extract metadata and gene expression matrix as CSVs instead
  cat("Extracting sparse matrix and metadata instead...\n")
  log_normalized_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  # Save gene names (rows)
  write.table(rownames(log_normalized_data), file = file.path(output_dir, "genes.tsv"), 
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

  # Save cell names (columns)
  write.table(colnames(log_normalized_data), file = file.path(output_dir, "barcodes.tsv"), 
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  # save gene expression sparse matrix
  Matrix::writeMM(log_normalized_data, file = file.path(output_dir, "expression.mtx"))


  # Extract the metadata
  meta_path <- file.path(output_dir, "metadata.csv")
  meta_data <- seurat_obj@meta.data
  write.csv(meta_data, file = meta_path, row.names = TRUE)

  # Free memory 
  rm(meta_data, log_normalized_data)
  gc()
}

# Free memory for Seurat object and other large variables
rm(seurat_obj, h5seurat_dir, h5ad_dir)
gc()

cat("Data Preparation completed.\n")