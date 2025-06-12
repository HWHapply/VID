# get DEG on ground truth
library(Seurat)
library(dplyr)

# take the input 
args <- commandArgs(trailingOnly = TRUE)

input_dir <- args[1]  # the directory of the input
random_state <- as.integer(args[2]) # the random seed of DEG

# Load expression matrix and metadata
expr <- read.csv(file.path(input_dir, "dmatrix.csv"), row.names = 1)
meta <- read.csv(file.path(input_dir, "metadata.csv"), row.names = 1)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = as(as.matrix(expr), "dgCMatrix"), meta.data = meta) 
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], layer = "data", new.data = as(as.matrix(expr), "dgCMatrix")) 
rm(expr, meta)
gc()
# Set identities based on 'label' column (make sure it's a factor or character)
seurat_obj$label <- as.factor(seurat_obj$label)

# Set the identity class to 'label' column
Idents(seurat_obj) <- "label"

# Run differential expression analysis between label = 1 and label = 0
deg <- FindMarkers(
  object = seurat_obj,
  slot = "data",
  ident.1 = "1",
  ident.2 = "0",
  test.use = "wilcox",       
  logfc.threshold = 0.25,     
  min.pct = 0.01,           
  only.pos = FALSE,          
  verbose = TRUE,
  random.seed = random_state
)

# Filter for statistically significant DEGs (adjusted p-value < 0.05)
deg_sig <- deg[deg$p_val_adj < 0.05, ]

# Save the significant gene names to a text file (one gene per line)
write.csv(
  deg_sig,
  file = file.path(input_dir, "DEGs.csv"),
)
