# ---------------------------------------------------------------------------
# Script 04: Clustering and Marker Gene Identification
# ---------------------------------------------------------------------------

source("00_setup_packages.R")
library(Seurat)
library(ggplot2)
library(dplyr) # Required for easy extraction of top markers

# Load the dimensionally reduced object
obj <- readRDS("output/03_dimred.rds")

# 1. Clustering
# (Note: FindNeighbors was already run in step 03, but is fine to rerun here 
# if you ever want to experiment with different dimensions)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)

# ---------------------------------------------------------------------------
# 2. Publication-Quality UMAP
# ---------------------------------------------------------------------------
p1 <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, label.size = 5, repel = TRUE, pt.size = 0.5) +
  theme_classic(base_size = 14) +
  labs(
    title = "UMAP: Cell Clustering",
    subtitle = "Resolution = 0.4"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
    legend.position = "right"
  )

ggsave("output/04_clusters_umap.png", p1, width = 9, height = 7, dpi = 600)

# ---------------------------------------------------------------------------
# 3. Find Marker Genes
# ---------------------------------------------------------------------------

# Increase the allowable size for parallel processing to 8 GB
options(future.globals.maxSize = 8000 * 1024^2)

message("Calculating cluster markers. This may take a moment...")
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)

# Save the full data table
write.csv(markers, "output/cluster_markers.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Visualize Top Markers (Crucial for Publication)
# ---------------------------------------------------------------------------
# Extract the top 5 genes per cluster based on fold change
top5_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Create a DotPlot of the top genes
p2 <- DotPlot(obj, features = unique(top5_markers$gene)) +
  theme_classic(base_size = 12) +
  labs(title = "Top Marker Genes per Cluster", x = "Marker Gene", y = "Cluster") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), # Genes should be italicized in print
    legend.position = "right"
  )

ggsave("output/04_cluster_dotplot.png", p2, width = 14, height = 6, dpi = 600)

# Save the clustered object
saveRDS(obj, "output/04_clustered.rds")
