library(Seurat)
library(ggplot2)

# Load the QC-filtered object
obj <- readRDS("output/01_qc_seurat.rds")

# Run SCTransform and PCA
obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)

# 1. Create the base plot with better colors and point sizing
p1 <- VariableFeaturePlot(obj, 
                          pt.size = 1.2, 
                          cols = c("grey80", "#D62728")) + # Soft grey for background genes, bold red for variable ones
  theme_classic(base_size = 14) +
  labs(
    title = "Highly Variable Genes (SCTransform)",
    subtitle = "Genes with high residual variance drive biological heterogeneity.",
    x = "Mean Expression",
    y = "Residual Variance"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    legend.position = "none" # The colors are self-explanatory; legend is redundant
  )

# 2. Extract top 10 genes to label
top10_genes <- head(VariableFeatures(obj), 10)

# 3. Add clean, repelled labels
p2 <- LabelPoints(plot = p1, 
                  points = top10_genes, 
                  repel = TRUE, 
                  size = 4.5,            # Slightly larger text
                  color = "black", 
                  fontface = "bold",     # Makes the gene names pop
                  max.overlaps = Inf)    # Forces ggrepel to label everything

# Save the enhanced plot and object
ggsave("output/variable_features.png", p2, width = 8, height = 6, dpi = 600)
saveRDS(obj, "output/02_sct_rds.rds")
