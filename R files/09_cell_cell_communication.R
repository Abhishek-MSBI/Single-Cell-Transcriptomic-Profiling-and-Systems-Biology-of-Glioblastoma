# ---------------------------------------------------------------------------
# Script 09: Cell-Cell Communication Analysis (CellChat)
# ---------------------------------------------------------------------------

source("00_setup_packages.R")
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)

# 1. Load the FULL annotated object (we need the microenvironment too!)
obj <- readRDS("output/05_annotated.rds")

# ---------------------------------------------------------------------------
# 2. Initialize CellChat (Seurat v5 Bypass Method)
# ---------------------------------------------------------------------------
message("Extracting data and initializing CellChat object...")

# Extract the normalized data matrix using the modern 'layer' syntax
data.input <- GetAssayData(obj, assay = "SCT", layer = "data")

# Extract the metadata 
meta <- obj@meta.data

# Create the CellChat object directly from the matrix and metadata
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")

# Set the reference database to Human
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# ---------------------------------------------------------------------------
# 3. Preprocessing and Inferring the Network
# ---------------------------------------------------------------------------
message("Identifying overexpressed genes and inferring communication...")
# Subset data to only include the genes present in the CellChat database to save memory
cellchat <- subsetData(cellchat)

# Identify overexpressed ligands and receptors in each cell group
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Calculate the probability of communication between cell types
# raw.use = TRUE uses the raw counts, which is standard practice for CellChat
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

# Filter out interactions where there are fewer than 10 cells in a group to prevent false positives
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the pathway-level communication (grouping individual genes into overall signaling pathways)
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the global aggregated network
cellchat <- aggregateNet(cellchat)

# ---------------------------------------------------------------------------
# 4. Publication-Quality Visualizations
# ---------------------------------------------------------------------------
message("Generating visual outputs...")

# Plot A: Global Number of Interactions Circle Plot
png("output/09a_global_interactions_circle.png", width = 8, height = 8, units = "in", res = 600)
# Increased the top margin (the 3rd number) from 2 to 6
par(mar = c(2, 2, 6, 2), xpd = TRUE) 
netVisual_circle(cellchat@net$count, weight.scale = TRUE, label.edge= FALSE, 
                 title.name = "Number of Interactions")
dev.off()

# Plot B: Global Interaction Strength Circle Plot
png("output/09b_global_interaction_strength.png", width = 8, height = 8, units = "in", res = 600)
# Increased the top margin here as well
par(mar = c(2, 2, 6, 2), xpd = TRUE)
netVisual_circle(cellchat@net$weight, weight.scale = TRUE, label.edge= FALSE, 
                 title.name = "Interaction Strength (Weight)")
dev.off()

# Plot C: Bubble Plot of Sender-Receiver Dynamics (Tumor to Macrophages)
# Shows the exact Ligand-Receptor pairs driving communication from Tumor -> Immune
p1 <- netVisual_bubble(cellchat, sources.use = "Tumor", 
                       targets.use = c("Microglia_Macrophage", "T_cells"), 
                       remove.isolate = FALSE) +
  theme_classic(base_size = 12) +
  labs(title = "Signaling: Tumor to Immune Microenvironment") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave("output/09c_Tumor_to_Immune_BubblePlot.png", p1, width = 8, height = 10, dpi = 600)

# ---------------------------------------------------------------------------
# 5. Export Data and Save
# ---------------------------------------------------------------------------
# Save the inferred communication network to a CSV
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "output/09_Inferred_Communication_Network.csv", row.names = FALSE)

# Save the CellChat object for future specialized plotting
saveRDS(cellchat, "output/09_cellchat_object.rds")
message("Cell-Cell Communication Analysis Complete.")
