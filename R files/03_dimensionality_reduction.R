# ---------------------------------------------------------------------------
# Script 03: Dimensionality Reduction and Clustering
# ---------------------------------------------------------------------------

source("00_setup_packages.R")
library(Seurat)
library(ggplot2)
library(patchwork)

# Load the SCTransformed object
obj <- readRDS("output/02_sct_rds.rds")

# 1. Correct Workflow Order: PCA -> Neighbors -> Clusters -> UMAP/tSNE
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30)
obj <- RunTSNE(obj, dims = 1:30)

# ---------------------------------------------------------------------------
# 2. Publication-Quality Plots
# ---------------------------------------------------------------------------

# PCA Plot
# Explicitly grouping by clusters so you see how they separate linearly
# PCA Plot with direct numbering
p1 <- DimPlot(obj, reduction = "pca", group.by = "seurat_clusters", 
              label = TRUE, label.size = 5, repel = TRUE, pt.size = 0.5) +
  theme_classic(base_size = 14) +
  labs(title = "PCA: Principal Component Analysis") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right", 
    legend.text = element_text(size = 12)
  )

# UMAP Plot
# label = TRUE and repel = TRUE puts the cluster numbers directly on the "islands"
p2 <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", 
              label = TRUE, label.size = 5, repel = TRUE, pt.size = 0.5) +
  theme_classic(base_size = 14) +
  labs(title = "UMAP: Cluster Visualization") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right", # Restoring your legend!
    legend.text = element_text(size = 12)
  )

# tSNE Plot
p3 <- DimPlot(obj, reduction = "tsne", group.by = "seurat_clusters", 
              label = TRUE, label.size = 5, repel = TRUE, pt.size = 0.5) +
  theme_classic(base_size = 14) +
  labs(title = "t-SNE: Cluster Visualization") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 12)
  )

# Bonus: Combine UMAP and tSNE for a side-by-side figure (often requested by journals)
combined_plot <- p2 + p3 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

# 3. Save Outputs (600 dpi for crisp text and graphics)
ggsave("output/pca_plot.png", p1, width = 8, height = 6, dpi = 600)
ggsave("output/umap_plot.png", p2, width = 8, height = 6, dpi = 600)
ggsave("output/tsne_plot.png", p3, width = 8, height = 6, dpi = 600)
ggsave("output/combined_umap_tsne.png", combined_plot, width = 14, height = 6, dpi = 600)

saveRDS(obj, "output/03_dimred.rds")
