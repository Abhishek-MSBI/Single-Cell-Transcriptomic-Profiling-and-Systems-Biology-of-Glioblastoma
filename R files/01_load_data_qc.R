# ---------------------------------------------------------------------------
# Script 01: Data Loading and Quality Control
# ---------------------------------------------------------------------------

# Load required packages
source("00_setup_packages.R")
library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(1234)

# 1. Define paths and load data
data_dir <- "C:/Users/abhir/Desktop/ScRNAseq"
h5_file <- file.path(data_dir, "Parent_SC3v3_Human_Glioblastoma_filtered_feature_bc_matrix.h5")

if (!file.exists(h5_file)) stop("File not found! Check the file path.")

counts <- Read10X_h5(h5_file)

# 2. Create Seurat Object and calculate Mitochondrial percentage
obj <- CreateSeuratObject(counts = counts, project = "GBM_10x", min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# ---------------------------------------------------------------------------
# 3. Generate Publication-Quality QC Plots
# ---------------------------------------------------------------------------

# Violin Plots
p1 <- VlnPlot(obj, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, 
              pt.size = 0) & 
  theme_classic(base_size = 14) & 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1), # Angles sample names for readability
    axis.ticks.x = element_blank(),
    legend.position = "none"        
  )

ggsave("output/qc_violin.png", p1, width = 12, height = 5, dpi = 600)

# Scatter Plots
p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

p3 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.5) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# Combine scatters with patchwork
scatter_combined <- p2 + p3 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave("output/qc_scatter.png", scatter_combined, width = 12, height = 5, dpi = 600)

# ---------------------------------------------------------------------------
# 4. Filter Cells and Save Object
# ---------------------------------------------------------------------------

# Apply QC thresholds
obj <- subset(obj, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt <= 20)

# Save the filtered object for downstream analysis
saveRDS(obj, "output/01_qc_seurat.rds")