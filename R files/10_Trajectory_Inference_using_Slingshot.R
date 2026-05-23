# ---------------------------------------------------------------------------
# Script 10: Trajectory Inference using Slingshot
# ---------------------------------------------------------------------------

source("00_setup_packages.R")
BiocManager::install(c("slingshot", "SingleCellExperiment", "RColorBrewer"))
library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(ggplot2)

# 1. Load the Reprocessed Tumor-Only Object
tumor <- readRDS("output/08_tumor_only_reprocessed.rds")

# 2. Convert Seurat object to SingleCellExperiment (required for Slingshot)
sce <- as.SingleCellExperiment(tumor)

# 3. Run Slingshot
message("Calculating developmental trajectories...")
# We use the UMAP embeddings and the sub-clusters we generated in Script 08
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')

# ---------------------------------------------------------------------------
# 4. Publication-Quality Trajectory Plot (Final Legend Fix)
# ---------------------------------------------------------------------------
message("Generating Trajectory Plot...")

# Widened the canvas to 10 inches to guarantee room for the text
png("output/10_tumor_pseudotime_trajectory.png", width = 10, height = 7, units = "in", res = 600)

# Increased the right margin (the 4th number) from 8 to 12
par(mar = c(5, 5, 4, 12), xpd = TRUE) 

# Set up the color palette
cluster_ids <- sort(unique(tumor$seurat_clusters))
colors <- colorRampPalette(brewer.pal(11, 'Spectral'))(length(cluster_ids))
plotcol <- colors[as.numeric(tumor$seurat_clusters)]

# Plot the UMAP coordinates
plot(reducedDims(sce)$UMAP, col = plotcol, pch = 16, cex = 0.5, 
     xlab = "UMAP 1", ylab = "UMAP 2", 
     main = "Tumor Trajectory Inference (Slingshot)",
     cex.main = 1.5, cex.lab = 1.2, bty = "l")

# Overlay the Principal Curves (The Trajectory Lines)
lines(SlingshotDataSet(sce), lwd = 3, col = 'black')

# Reduced the inset slightly (from -0.3 to -0.25) so it anchors safely inside the new wider margin
legend("topright", inset = c(-0.25, 0), 
       legend = paste("Sub-clone", cluster_ids),
       col = colors, pch = 16, pt.cex = 1.5, 
       title = "Tumor Clusters", bty = "n", cex = 1.1)

dev.off()

# 5. Extract Pseudotime Values to Metadata
# This allows you to see how "far along" the trajectory each cell is
tumor$pseudotime <- slingPseudotime(sce)[,1] # Takes the primary lineage
write.csv(as.data.frame(tumor@meta.data), "output/10_tumor_metadata_with_pseudotime.csv")

saveRDS(tumor, "output/10_tumor_slingshot.rds")
message("Trajectory Inference Complete.")
