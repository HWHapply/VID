# save the predicted infection status to seurat object
library(Seurat)

# Take positional arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Error: Two positional arguments required: <rawdata_dir> <output_dir>")
}

rawdata_dir <- args[1]  # Directory of the Seurat .rds file
output_dir <- args[2]   # Output directory to check for files

cat("Saving result to seurat object...")
# Detect files in the output directory
files <- list.files(output_dir, full.names = TRUE)

# Check if only 'metadata.csv' file is present
metadata_file <- grep("metadata\\.csv$", files, value = TRUE)

if (length(metadata_file) == 1) {
  # Read the metadata CSV
  cat("Reading 'metadata.csv'...\n")
  metadata <- read.csv(metadata_file, row.names = 1)
  
  # Load the Seurat object from rawdata_dir
  cat("Loading Seurat object from rawdata_dir...\n")
  seurat_obj <- readRDS(rawdata_dir)
  
  # Replace the metadata in the Seurat object
  cat("Replacing metadata in the Seurat object...\n")
  seurat_obj@meta.data <- metadata
  
  # Save the updated Seurat object as 'data.rds'
  saveRDS(seurat_obj, file = file.path(output_dir, "data.rds"))
  cat("'data.rds' with updated metadata saved successfully in the output directory.\n")
  
  # Remove variables to free up memory
  rm(seurat_obj, metadata, metadata_file, files)
  gc()
  
} else {
  stop("Error: Unexpected number or type of files in the output directory. Expecting only 'metadata.csv'.")
}

cat("Result saved to seurat object!\n")