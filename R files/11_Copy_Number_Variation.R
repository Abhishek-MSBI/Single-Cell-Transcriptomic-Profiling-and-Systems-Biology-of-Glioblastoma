# ---------------------------------------------------------------------------
# Script 11: Copy Number Variation (inferCNV)
# ---------------------------------------------------------------------------

source("00_setup_packages.R")
# BiocManager::install("infercnv")
library(Seurat)
library(infercnv)

# 1. Load the FULL annotated object (We need the healthy cells as a reference)
obj <- readRDS("output/05_annotated.rds")

# 2. Prepare the Inputs
message("Extracting raw counts and annotations...")
# inferCNV needs raw, un-normalized counts
counts_matrix <- GetAssayData(obj, assay = "RNA", layer = "counts")

# Create a simple 2-column dataframe for annotations (Cell ID and Cell Type)
annotations <- data.frame(Cell_Type = obj$cell_type)
rownames(annotations) <- colnames(obj)

# 3. Create the inferCNV Object
message("Building inferCNV object...")
# NOTE: You must provide a path to a standard GRCh38 gene ordering file
gene_ordering_file_path <- "data/gencode_v21_gen_pos.txt" # Update this path to your downloaded file

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file = annotations,
  delim = "\t",
  gene_order_file = gene_ordering_file_path,
  # Define the healthy cells to serve as the normal baseline
  ref_group_names = c("Astrocyte", "Microglia_Macrophage", "T_cells")
)

# 4. Run the CNV Analysis
# Warning: This is computationally heavy and may take an hour or more depending on your RAM.
message("Running inferCNV... (This will take a while)")
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1, # 0.1 for 10x Genomics data
  out_dir = "output/11_inferCNV_results", 
  cluster_by_groups = TRUE, 
  denoise = TRUE,
  HMM = FALSE, # Set to TRUE if you want strict mathematical bounds, but it vastly increases run time
  output_format = "png"
)

message("inferCNV Complete. Heatmaps are saved in the output/11_inferCNV_results directory.")