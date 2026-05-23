# ---------------------------------------------------------------------------
# Script 12: Transcription Factor Activity (DoRothEA & Viper)
# ---------------------------------------------------------------------------

source("00_setup_packages.R")
# BiocManager::install(c("dorothea", "viper"))
library(Seurat)
library(dorothea)
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Load the FULL annotated object
obj <- readRDS("output/05_annotated.rds")

# We only want to compare Tumor vs Astrocytes for this specific analysis
obj_subset <- subset(obj, idents = c("Tumor", "Astrocyte"))

# 2. Load the Regulon Database
message("Loading Human TF database...")
# We use only highly confident interactions (A, B, and C levels)
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
regulon <- dorothea_regulon_human %>% filter(confidence %in% c("A", "B", "C"))

# 3. Run Viper Analysis
message("Calculating Transcription Factor activities...")
# This calculates the activity of every TF in every single cell
obj_subset <- run_viper(obj_subset, regulon, 
                        options = list(method = "scale", minsize = 4, 
                                       eset.filter = FALSE, cores = 1, verbose = FALSE))

# 4. Differential TF Activity (Tumor vs Astrocyte)
message("Identifying top tumor-driving Transcription Factors...")
DefaultAssay(obj_subset) <- "dorothea" # Switch Seurat's focus to the new TF data
obj_subset <- ScaleData(obj_subset)

tf_markers <- FindMarkers(obj_subset, ident.1 = "Tumor", ident.2 = "Astrocyte", 
                          only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
tf_markers$TF <- rownames(tf_markers)

# Save the raw results
write.csv(tf_markers, "output/12_Tumor_TF_Activity.csv", row.names = FALSE)

# 5. Publication-Quality Visualization
message("Generating TF Heatmap...")
# Extract the top 20 most active Master Regulators in the Tumor
top20_tfs <- head(tf_markers$TF, 20)

p1 <- DoHeatmap(obj_subset, features = top20_tfs, size = 4, angle = 90) +
  scale_fill_gradient2(low = "#1F77B4", mid = "white", high = "#D62728", midpoint = 0) +
  theme(axis.text.y = element_text(size = 10, face = "bold")) +
  labs(title = "Master Regulators: Top Active TFs in Tumor Cells")

ggsave("output/12_top_TFs_heatmap.png", p1, width = 12, height = 8, dpi = 600)

saveRDS(obj_subset, "output/12_dorothea_scored.rds")
message("Gene Regulatory Network Analysis Complete.")