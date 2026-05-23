# ---------------------------------------------------------------------------
# Script 05: Cell Type Annotation
# ---------------------------------------------------------------------------

source("00_setup_packages.R")
library(Seurat)
library(ggplot2)

# Load the clustered object
obj <- readRDS("output/04_clustered.rds")

# ---------------------------------------------------------------------------
# 1. Visualize Canonical Markers
# ---------------------------------------------------------------------------
# Plotting known marker genes to see which numeric cluster expresses them
canonical_markers <- c("PTPRC", "CD3D", "NKG7", "MS4A1", "LYZ", "PECAM1", "PDGFRA", "AQP4", "GFAP")

p1 <- DotPlot(obj, features = canonical_markers) + 
  theme_classic(base_size = 14) +
  labs(
    title = "Expression of Canonical Markers", 
    x = "Marker Gene", 
    y = "Cluster"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), # Genes italicized for print
    legend.position = "right"
  )

ggsave("output/05_canonical_markers_dotplot.png", p1, width = 12, height = 6, dpi = 600)

# ---------------------------------------------------------------------------
# 2. Assign Biological Cell Types (Corrected for 11 Clusters)
# ---------------------------------------------------------------------------

# Providing exactly 11 names corresponding to Clusters 0 through 10
new.ids <- c(
  "Tumor",                # Cluster 0: Low canonical expression
  "Astrocyte",            # Cluster 1: GFAP+, AQP4+
  "Tumor",                # Cluster 2: Low canonical expression
  "Astrocyte",            # Cluster 3: Strong GFAP+, AQP4+
  "Astrocyte",            # Cluster 4: Strong GFAP+, AQP4+
  "Microglia_Macrophage", # Cluster 5: PTPRC+, LYZ+
  "Microglia_Macrophage", # Cluster 6: PTPRC+, LYZ+
  "T_cells",              # Cluster 7: PTPRC+, CD3D+
  "OPC",                  # Cluster 8: PDGFRA+
  "OPC",                  # Cluster 9: PDGFRA+
  "Tumor"                 # Cluster 10: Low canonical expression
)

# Safety check (This will now pass smoothly!)
if(length(new.ids) != length(levels(obj))) {
  warning(paste("Mismatch! You provided", length(new.ids), "names, but there are", length(levels(obj)), "clusters."))
}

# Map the new names to the current numerical identities
names(new.ids) <- levels(obj)
obj <- RenameIdents(obj, new.ids)
obj$cell_type <- Idents(obj)

# ---------------------------------------------------------------------------
# 3. Generate the Final Annotated UMAP
# ---------------------------------------------------------------------------
p2 <- DimPlot(obj, reduction = "umap", group.by = "cell_type", 
              label = TRUE, label.size = 4.5, repel = TRUE, pt.size = 0.5) +
  theme_classic(base_size = 14) +
  labs(title = "UMAP: Annotated Cell Types") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "right",
    legend.text = element_text(size = 12)
  )

ggsave("output/05_annotated_umap.png", p2, width = 10, height = 7, dpi = 600)

# ---------------------------------------------------------------------------
# Script 05b: Exporting Cell Type Counts and Proportions
# ---------------------------------------------------------------------------

library(Seurat)
library(dplyr)

# Load the annotated object (if not already loaded in your session)
obj <- readRDS("output/05_annotated.rds")

# 1. Extract the raw counts for each cell type
cell_counts <- as.data.frame(table(obj$cell_type))
colnames(cell_counts) <- c("Cell_Type", "Count")

# 2. Calculate the percentage of total cells for each type
total_cells <- sum(cell_counts$Count)
cell_counts$Percentage <- round((cell_counts$Count / total_cells) * 100, 2)

# 3. Sort the table from most abundant to least abundant
cell_counts <- cell_counts %>% arrange(desc(Count))

# 4. Save to a CSV file
write.csv(cell_counts, "output/05_cell_type_summary.csv", row.names = FALSE)

# Print to console so you can see it immediately
print(cell_counts)

# Save the fully annotated object
saveRDS(obj, "output/05_annotated.rds")
