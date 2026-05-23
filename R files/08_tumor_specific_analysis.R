# ---------------------------------------------------------------------------
# Script 08: Tumor-Specific Analysis and Sub-Clustering
# ---------------------------------------------------------------------------

source("00_setup_packages.R")
library(Seurat)
library(ggplot2)
library(patchwork)

# Load the fully annotated master object
obj <- readRDS("output/05_annotated.rds")

# ---------------------------------------------------------------------------
# 1. Isolate the Tumor Compartment
# ---------------------------------------------------------------------------
message("Subsetting Tumor cells...")
tumor <- subset(obj, idents = "Tumor")

# ---------------------------------------------------------------------------
# 2. Cell Cycle Scoring
# ---------------------------------------------------------------------------
# Seurat has built-in lists of cell cycle markers (cc.genes.updated.2019)
message("Calculating Cell Cycle Scores...")
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

tumor <- CellCycleScoring(tumor, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Plot A: Cell Cycle Distribution
p1 <- DimPlot(tumor, reduction = "umap", group.by = "Phase", pt.size = 1) +
  theme_classic(base_size = 14) +
  labs(title = "Tumor Cells: Pre-Reclustering Cell Cycle Phase") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("output/08a_tumor_cellcycle_pre_reclustering.png", p1, width = 8, height = 6, dpi = 600)

# ---------------------------------------------------------------------------
# 3. Re-Processing the Tumor Object (Critical Step)
# ---------------------------------------------------------------------------
# Because we removed all other cell types, we must recalculate the variance
message("Re-calculating variance and PCA for Tumor isolation...")
tumor <- FindVariableFeatures(tumor, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
tumor <- ScaleData(tumor, verbose = FALSE)
tumor <- RunPCA(tumor, npcs = 30, verbose = FALSE)

# Re-cluster to find tumor sub-clones
tumor <- FindNeighbors(tumor, dims = 1:20, verbose = FALSE)
# Using a slightly lower resolution since this is already a subset
tumor <- FindClusters(tumor, resolution = 0.3, verbose = FALSE) 
tumor <- RunUMAP(tumor, dims = 1:20, verbose = FALSE)

# ---------------------------------------------------------------------------
# 4. Visualize the Intra-Tumoral Heterogeneity
# ---------------------------------------------------------------------------
# Plot B: The new UMAP showing distinct Tumor Sub-Clusters
p2 <- DimPlot(tumor, reduction = "umap", group.by = "seurat_clusters", 
              label = TRUE, label.size = 5, repel = TRUE, pt.size = 1) +
  theme_classic(base_size = 14) +
  labs(title = "UMAP: Tumor Sub-Clones", subtitle = "Resolution = 0.3") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# Plot C: Cell Cycle Overlay on the NEW Sub-Clusters
p3 <- DimPlot(tumor, reduction = "umap", group.by = "Phase", pt.size = 1) +
  theme_classic(base_size = 14) +
  labs(title = "Cell Cycle Distribution Across Sub-Clones") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# Combine sub-cluster and cell cycle plots side-by-side
combined_tumor_plot <- p2 + p3
ggsave("output/08b_tumor_subclones_and_cellcycle.png", combined_tumor_plot, width = 14, height = 6, dpi = 600)

# ---------------------------------------------------------------------------
# 5. Export Data
# ---------------------------------------------------------------------------
message("Exporting metadata and saving isolated object...")

# Export metadata including the new cell cycle scores and sub-cluster IDs
write.csv(as.data.frame(tumor@meta.data), "output/08_tumor_metadata_annotated.csv")

# Save the newly processed, tumor-only Seurat object
saveRDS(tumor, "output/08_tumor_only_reprocessed.rds")
message("Tumor-specific analysis complete.")
