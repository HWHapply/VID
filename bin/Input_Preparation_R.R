library(Seurat)
library(SeuratDisk)

# Take positional arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Error: Two positional arguments required: <seurat_obj_dir> <output_dir>")
}

seurat_obj_dir <- args[1]  # The directory of the .rds file
output_dir <- args[2]  # The output directory

# Load the Seurat object
seurat_obj <- readRDS(seurat_obj_dir)

# normalize and select top 2000 variance genes
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

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
  cat("Extracting gene expression and metadata as CSVs instead...\n")
  
  # Extract the metadata
  meta_path <- file.path(output_dir, "metadata.csv")
  meta_data <- seurat_obj@meta.data
  write.csv(meta_data, file = meta_path, row.names = TRUE)
  
  # Free memory after metadata extraction
  rm(meta_data)
  gc()
  
  # Extract the log-normalized gene expression matrix
  output_file <- file.path(output_dir, "dmatrix.csv")
  num_genes <- 2000
  top_genes <- VariableFeatures(seurat_obj)[1:num_genes]
  log_normalized_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[top_genes, ]
  
  # Free memory after extracting log-normalized data
  rm(top_genes)
  gc()
  
  # Handle large matrices by saving in chunks if necessary
  tryCatch({
    dmatrix <- as.data.frame(as.matrix(log_normalized_data))
    write.csv(dmatrix, file = output_file, row.names = TRUE)
    # Free memory after saving
    rm(dmatrix)
    gc()
  }, error = function(e) {
    cat("Error during CSV export: ", e$message, "\n")
    # Save in chunks
    chunk_size <- 1000
    chunk <- log_normalized_data[, 1:chunk_size]
    chunk_df <- as.data.frame(as.matrix(chunk))
    fwrite(chunk_df, file = output_file, row.names = TRUE)
    
    # Free memory after first chunk
    rm(chunk, chunk_df)
    gc()
    
    for (i in seq(chunk_size + 1, ncol(log_normalized_data), by = chunk_size)) {
      chunk <- log_normalized_data[, i:min(i + chunk_size - 1, ncol(log_normalized_data))]
      chunk_df <- as.data.frame(as.matrix(chunk))
      fwrite(chunk_df, file = output_file, row.names = TRUE, append = TRUE, col.names = FALSE)
      
      # Free memory after each chunk
      rm(chunk, chunk_df)
      gc()
    }
  })
  
  # Free memory after processing log-normalized data
  rm(log_normalized_data)
  gc()
}

# Free memory for Seurat object and other large variables
rm(seurat_obj, h5seurat_dir, h5ad_dir)
gc()

cat("Data Preparation completed.\n")