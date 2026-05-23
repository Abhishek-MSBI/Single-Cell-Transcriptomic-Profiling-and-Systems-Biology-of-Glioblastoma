# ---------------------------------------------------------------------------
# Script 07: Pathway Enrichment Analysis (GO & KEGG)
# ---------------------------------------------------------------------------

source("00_setup_packages.R")

# Ensure these packages are installed:
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# ---------------------------------------------------------------------------
# 1. Load Data & Convert Gene IDs
# ---------------------------------------------------------------------------
# We load the clean, pre-filtered upregulated genes from Step 06
deg_up <- read.csv("output/06_Tumor_Upregulated_Only.csv")

# clusterProfiler requires ENTREZ IDs instead of Gene Symbols. 
# bitr (Biological Id TRanslator) converts our gene names to these database IDs.
gene_map <- bitr(deg_up$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Create a named vector of Fold Changes (required later for the network plot)
fold_changes <- deg_up$avg_log2FC
names(fold_changes) <- deg_up$gene

# ---------------------------------------------------------------------------
# 2. Gene Ontology (GO) Enrichment
# ---------------------------------------------------------------------------
message("Running GO Enrichment (Biological Processes)...")
ego <- enrichGO(gene          = gene_map$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",          # BP = Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)          # CRITICAL: Converts Entrez IDs back to Gene Symbols in the CSV!

write.csv(as.data.frame(ego), "output/07_GO_Enrichment_Results.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. KEGG Pathway Enrichment
# ---------------------------------------------------------------------------
message("Running KEGG Enrichment...")
ekegg <- enrichKEGG(gene         = gene_map$ENTREZID,
                    organism     = 'hsa',      # hsa = Homo sapiens
                    pvalueCutoff = 0.05)

# Translate KEGG IDs back to gene symbols for human readability
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(as.data.frame(ekegg), "output/07_KEGG_Enrichment_Results.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Publication-Quality Visualizations
# ---------------------------------------------------------------------------

# Plot A: Enhanced DotPlot for GO terms
p1 <- dotplot(ego, showCategory = 10, title = "GO: Biological Processes Enriched in Tumor") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 11)
  )
ggsave("output/07a_GO_dotplot.png", p1, width = 10, height = 7, dpi = 600)

# Plot B: Category Netplot (Cnetplot)
# This connects the top 5 enriched pathways to the specific genes driving them.
p2 <- cnetplot(ego, foldChange = fold_changes, showCategory = 5) +
  ggtitle("Gene-Pathway Network (Top 5 GO Terms)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))

ggsave("output/07b_GO_cnetplot.png", p2, width = 12, height = 10, dpi = 600)

# Plot C: KEGG Dotplot
p3 <- dotplot(ekegg, showCategory = 10, title = "KEGG Pathways Enriched in Tumor") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
ggsave("output/07c_KEGG_dotplot.png", p3, width = 10, height = 7, dpi = 600)

message("Pathway Enrichment Complete.")
