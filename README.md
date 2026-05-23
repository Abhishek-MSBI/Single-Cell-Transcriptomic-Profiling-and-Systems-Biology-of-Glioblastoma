# Single-Cell Transcriptomic Profiling and Systems Biology of Glioblastoma

**Author:** Abhishek S R 
**Domain:** Bioinformatics, Cancer Genomics, Single-Cell RNA-seq  

## Project Overview

This repository contains a modular, end-to-end bioinformatics pipeline for single-cell RNA sequencing (scRNA-seq) analysis of human Glioblastoma Multiforme (GBM). The dataset consists of 5,604 high-quality cells profiled on the 10x Genomics Chromium platform (3' v3 chemistry) from a 57-year-old male donor. Sequencing was performed on an Illumina NovaSeq 6000 at a mean depth of 44,736 reads per cell, with a median of 3,094 detected genes and 9,434 UMIs per cell, providing deep coverage suitable for advanced systems biology. 

The primary goals are to characterize intra-tumoral heterogeneity, map the tumor microenvironment (TME), and identify tumor-intrinsic oncogenic programs, pathways, and master regulators.

## Pipeline Architecture & Results Discussion

The pipeline is implemented in R and organized into sequential scripts under `R/`, leveraging Seurat v5 as the core framework alongside specialized systems biology packages. The generated visual outputs are stored in `results/plots/`.

### Phase 1: Pre-processing & Core Analytics

**1. Quality Control (`01_load_data_qc.R`)**
* **Why it was done:** To remove empty droplets, stressed, or dying cells that could skew the analysis.
* **Packages used:** `Seurat v5`.
* **Results & Interpretation:** Cells were filtered based on total RNA counts, detected features, and mitochondrial read percentages. This enriched the dataset for robust transcriptional profiles, yielding 5,604 high-quality cells. 
* **Plots Generated:** `qc_scatter.png`, `qc_violin.png`.

**2. Normalization & Feature Selection (`02_normalization_feature_selection.R`)**
* **Why it was done:** To regress out technical variation (e.g., sequencing depth) and identify highly variable genes that drive true biological signals.
* **Packages used:** `Seurat v5` (SCTransform).
* **Results & Interpretation:** SCTransform successfully normalized the counts and defined the subset of variable genes used for downstream clustering. 
* **Plots Generated:** `variable_features.png`.

**3. Dimensionality Reduction (`03_dimensionality_reduction.R`)**
* **Why it was done:** To reduce the high-dimensional gene expression data into a manageable format and visually group transcriptionally similar cells.
* **Packages used:** `Seurat v5`.
* **Results & Interpretation:** PCA was computed on the variable features, followed by UMAP and t-SNE embeddings to map the cellular landscape in a nonlinear 2D space.
* **Plots Generated:** `pca_plot.png`, `tsne_plot.png`, `umap_plot.png`, `combined_umap_tsne.png`.

**4. Clustering & Annotation (`04_clustering.R`, `05_annotation.R`)**
* **Why it was done:** To map the tumor microenvironment (TME) and assign discrete biological identities to the cell clusters using canonical markers.
* **Packages used:** `Seurat v5`.
* **Results & Interpretation:** The analysis revealed a classic GBM ecosystem dominated by an astroglial-like compartment within an immune-rich environment. The final composition consisted of **Astrocytes** (1,921 cells, 41.08%), **Tumor cells** (1,885 cells, 40.31%), **Microglia/Macrophages** (476 cells, 10.18%), **OPCs** (252 cells, 5.39%), and **T cells** (142 cells, 3.04%).
* **Plots Generated:** `04_cluster_dotplot.png`, `05_annotated_umap.png`, `05_canonical_markers_dotplot.png`.

### Phase 2: Differential Expression & Functional Enrichment

**5. Differential Expression (`06_differential_expression.R`)**
* **Why it was done:** To delineate tumor-intrinsic oncogenic programs from normal glial biology by comparing the malignant compartment directly against reference astrocytes.
* **Packages used:** `Seurat v5`.
* **Results & Interpretation:** The DE analysis successfully isolated tumor-upregulated genes (capturing candidate oncogenic effectors, proliferation, and stemness markers) and tumor-downregulated genes (highlighting astrocytic homeostatic programs lost during malignant transformation). 
* **Plots Generated:** `06a_volcano_plot.png`, `06b_top_genes_featureplot.png`, `06c_top_genes_violin.png`, `06d_top_genes_heatmap.png`.

**6. Pathway Enrichment (`07_pathway_enrichment.R`)**
* **Why it was done:** To translate the differentially expressed gene signatures into functional biological networks and therapeutic vulnerabilities.
* **Packages used:** `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`.
* **Results & Interpretation:** Gene Ontology (GO) enrichment of the upregulated genes showed strong aggregation in cell cycle regulation, DNA replication, and mitotic spindle organization, confirming a highly proliferative GBM phenotype. KEGG enrichment further linked these signatures to canonical oncogenic signaling circuits.
* **Plots Generated:** `07a_GO_dotplot.png`, `07b_GO_cnetplot.png`, `07c_KEGG_dotplot.png`.

### Phase 3: Advanced Systems Biology & Tumor Heterogeneity

**7. Sub-clonal Architecture (`08_tumor_specific_analysis.R`)**
* **Why it was done:** To dissect intra-tumoral heterogeneity and investigate functional differences between distinct tumor cell states.
* **Packages used:** `Seurat v5`.
* **Results & Interpretation:** By performing cell cycle scoring (S, G2/M phases) and reclustering exclusively on the tumor subset, we identified distinct proliferative versus quiescent subpopulations. This fine-grained stratification links specific sub-clones to unique functional niches.
* **Plots Generated:** `08a_tumor_cellcycle_pre_reclustering.png`, `08b_tumor_subclones_and_cellcycle.png`.

**8. Cell-Cell Communication (`09_cell_cell_communication.R`)**
* **Why it was done:** To model intercellular cross-talk within the GBM microenvironment, uncovering how malignant cells interact with stromal and immune components.
* **Packages used:** `CellChat`.
* **Results & Interpretation:** The inferred communication network highlighted key ligand-receptor interactions, demonstrating how tumor cells, astrocytes, and macrophages communicate via specific growth factors, cytokines, and immune regulatory axes. 
* **Plots Generated:** `09a_global_interactions_circle.png`, `09b_global_interaction_strength.png`, `09c_Tumor_to_Immune_BubblePlot.png`.

**9. Pseudotime & Trajectory Inference (`10_pseudotime.R`)**
* **Why it was done:** To reconstruct putative evolutionary or phenotypic trajectories, mapping how tumor cells transition between different states.
* **Packages used:** `slingshot`, `SingleCellExperiment`.
* **Results & Interpretation:** Pseudotime analysis distinguished early-like from late-like tumor states. The resulting lineage trajectories suggest potential evolutionary paths, which may correspond to dedifferentiation or therapy-induced adaptations.
* **Plots Generated:** `10_tumor_pseudotime_trajectory.png`.

**10. Master Regulators (`11_Transcription_factor_activity.R`)**
* **Why it was done:** To identify the core transcription factors (TFs) driving the observed malignant phenotypes.
* **Packages used:** `DoRothEA`, `viper`.
* **Results & Interpretation:** By inferring TF activity from the expression profiles, we pinpointed candidate master regulators with strongly perturbed activity. These TFs likely orchestrate the GBM-specific oncogenic programs, providing mechanistic hypotheses for targeted drug interventions.
* **Plots Generated:** `11_top_TFs_heatmap.png`.

## Installation and Usage

1. Clone the repository:
   ```bash
   git clone [https://github.com/](https://github.com/)<your-username>/SingleCell_GBM_Pipeline.git
   cd SingleCell_GBM_Pipeline
