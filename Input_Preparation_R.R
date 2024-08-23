"
Prepare the input for VID
Two options for VID input:
  A. h5ad file (.h5ad)
  B. gene expression table (.csv) + metadata table (.csv)
"

# Install packages: Seurat, SeuratDisk (optional)
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("satijalab/seurat", ref = "v4.4.0")
# remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(SeuratDisk)

# Load seurat object
seurat_obj_dir <- './demo/data/demo.rds' # seurat object directory
seurat_obj <- readRDS(seurat_obj_dir) # load seurat object 

# Please choose the method A or method B according to your seurat version to generate the input for vid

# A. If your Seurat object version is not higher than V5, convert to h5ad file with the code below
#   seurat object (.rds) -> h5seurat (.h5seurat) -> h5ad (.h5ad)
h5seurat_dir <- '' # h5seurat file directory
h5ad_dir <- '' # h5ad file directory (VID input)
SaveH5Seurat(seurat_obj, filename = h5seurat_dir) # convert and save seurat object to h5seurat file
Convert(h5ad_dir, dest = "h5ad") # convert h5seurat file to h5ad file

# B. For seurat V5 or above, Extract the gene expression and meta data as two csv tables from seurat object
#   seurat object (.rds) -> gene expression matrix (.csv) + metadata(.csv)

# Extract the metadata
meta_path <- './demo/data/metadata.csv' # metadata table directory
meta_data <- seurat_obj@meta.data # extract meta-data from seurat object
write.csv(meta_data, file = meta_path, row.names = TRUE) # save

# Extract the log normalized gene expression matrix
output_file <- './demo/data/dmatrix.csv' # gene expression table output directory
# First, identify the top high variance genes
num_genes <- 2000 # set the number of top genes
top_genes <- VariableFeatures(seurat_obj)[1:num_genes] # get high variance genes
log_normalized_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[top_genes, ] # filter the gene expression

# save the gene expression table directly if sample size is not too big, otherwise, the table will be saved with chunk
result <- tryCatch({
  # Code that might throw an error
  dmatrix <- as.data.frame(as.matrix(log_normalized_data))
  write.csv(dmatrix, file = output_file, row.names = TRUE)
}, warning = function(w) {
  # Handle warnings
  cat("A warning occurred: ", w$message, "\n")
  list(success = FALSE, data = NULL)
}, error = function(e) {
  # Handle errors
  cat("An error occurred: ", e$message, "\n")
  list(success = FALSE, data = NULL)
  # Code to execute regardless of error
  # Define the chunk size
  chunk_size <- 1000
  
  # First chunk: write with header
  chunk <- log_normalized_data[, 1:chunk_size]
  chunk_df <- as.data.frame(as.matrix(chunk))
  fwrite(chunk_df, file = output_file, row.names = TRUE)
  
  # Subsequent chunks: append without header
  for (i in seq(chunk_size + 1, ncol(log_normalized_data), by = chunk_size)) {
    chunk <- log_normalized_data[, i:min(i + chunk_size - 1, ncol(log_normalized_data))]
    chunk_df <- as.data.frame(as.matrix(chunk))
    fwrite(chunk_df, file = output_file, row.names = TRUE, append = TRUE, col.names = FALSE)
  }
}, finally = {
  cat("Execution completed.")
})

# If both methods failed, please consider downgrade the seurat to V4.