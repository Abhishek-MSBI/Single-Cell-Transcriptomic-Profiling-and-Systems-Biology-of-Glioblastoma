# ---------------------------------------------------------------------------
# Script 06: Complete Differential Expression Analysis
# ---------------------------------------------------------------------------

source("00_setup_packages.R")
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

# Load the annotated object
obj <- readRDS("output/05_annotated.rds")

# ---------------------------------------------------------------------------
# 1. Run Differential Expression & Export CSV
# ---------------------------------------------------------------------------
message("Running DGE: Tumor vs Astrocyte...")
dge_results <- FindMarkers(obj, ident.1 = "Tumor", ident.2 = "Astrocyte", 
                           logfc.threshold = 0.25, min.pct = 0.1)

dge_results$gene <- rownames(dge_results)

# ESSENTIAL: Save the raw results for your manuscript's Supplementary Tables
write.csv(dge_results, "output/06_Tumor_vs_Astrocyte_DGE.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Extract Top Genes Automatically
# ---------------------------------------------------------------------------
dge_results$Significance <- "Not Significant"
dge_results$Significance[dge_results$avg_log2FC > 1 & dge_results$p_val_adj < 0.05] <- "Upregulated in Tumor"
dge_results$Significance[dge_results$avg_log2FC < -1 & dge_results$p_val_adj < 0.05] <- "Downregulated in Tumor"

# Get the top 4 highly upregulated genes in the Tumor for specific plotting
top4_tumor_genes <- dge_results %>% 
  filter(Significance == "Upregulated in Tumor") %>% 
  slice_max(order_by = avg_log2FC, n = 4) %>% 
  pull(gene)

# Get top 10 up/down for the Volcano labels
top_volcano_genes <- dge_results %>% 
  filter(Significance != "Not Significant") %>% 
  group_by(Significance) %>% 
  slice_max(order_by = abs(avg_log2FC), n = 10)

# ---------------------------------------------------------------------------
# 3. Plot A: The Volcano Plot (Global View)
# ---------------------------------------------------------------------------
p_volcano <- ggplot(dge_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated in Tumor" = "#D62728", 
                                "Downregulated in Tumor" = "#1F77B4", 
                                "Not Significant" = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(data = top_volcano_genes, aes(label = gene), 
                  size = 4, color = "black", fontface = "italic", max.overlaps = 20) +
  theme_classic(base_size = 14) +
  labs(title = "Volcano Plot: Tumor vs. Astrocytes", x = "Log2 Fold Change", y = "-Log10(Adj P-Value)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom", legend.title = element_blank())

ggsave("output/06a_volcano_plot.png", p_volcano, width = 10, height = 8, dpi = 600)

# ---------------------------------------------------------------------------
# 4. Plot B: Feature Plots (Spatial View)
# ---------------------------------------------------------------------------
# This shows exactly where the top 4 tumor-driving genes are expressed on the UMAP
p_feature <- FeaturePlot(obj, features = top4_tumor_genes, pt.size = 0.5, ncol = 2) &
  theme(plot.title = element_text(face = "bold.italic")) # Italicizes gene names

ggsave("output/06b_top_genes_featureplot.png", p_feature, width = 10, height = 8, dpi = 600)

# ---------------------------------------------------------------------------
# 5. Plot C: Violin Plots (Distribution View)
# ---------------------------------------------------------------------------
# We subset the object to ONLY plot Tumor and Astrocytes for a clean, direct comparison
obj_subset <- subset(obj, idents = c("Tumor", "Astrocyte"))

p_violin <- VlnPlot(obj_subset, features = top4_tumor_genes, pt.size = 0, ncol = 2) &
  theme_classic(base_size = 14) &
  theme(
    plot.title = element_text(face = "bold.italic", hjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none"
  )
ggsave("output/06c_top_genes_violin.png", p_violin, width = 10, height = 8, dpi = 600)
# ---------------------------------------------------------------------------
# 6. Plot D: The Single-Cell Heatmap
# ---------------------------------------------------------------------------
# Extract the top 10 upregulated and top 10 downregulated genes
top20_heatmap_genes <- dge_results %>% 
  filter(Significance != "Not Significant") %>% 
  group_by(Significance) %>% 
  slice_max(order_by = abs(avg_log2FC), n = 10) %>% 
  pull(gene)

# Scale the data for the subset (required for DoHeatmap to look correct)
obj_subset <- ScaleData(obj_subset, features = top20_heatmap_genes, verbose = FALSE)

# Generate the heatmap
p_heatmap <- DoHeatmap(obj_subset, features = top20_heatmap_genes, size = 4, angle = 90) +
  theme(
    axis.text.y = element_text(size = 10, face = "italic"), # Italicize gene names
    legend.position = "right"
  ) +
  scale_fill_gradientn(colors = c("#1F77B4", "white", "#D62728")) # Blue to White to Red

ggsave("output/06d_top_genes_heatmap.png", p_heatmap, width = 12, height = 8, dpi = 600)

# ---------------------------------------------------------------------------
# 7. Export Directional Gene Lists for Pathway Analysis
# ---------------------------------------------------------------------------
# Isolate just the upregulated genes
upregulated_only <- dge_results %>% filter(Significance == "Upregulated in Tumor")
write.csv(upregulated_only, "output/06_Tumor_Upregulated_Only.csv", row.names = FALSE)

# Isolate just the downregulated genes
downregulated_only <- dge_results %>% filter(Significance == "Downregulated in Tumor")
write.csv(downregulated_only, "output/06_Tumor_Downregulated_Only.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 8. Save the Final Seurat Object
# ---------------------------------------------------------------------------
saveRDS(obj, "output/06_dge_completed.rds")
message("DGE Analysis Complete. Object saved to output/06_dge_completed.rds")
