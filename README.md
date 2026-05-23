# Single-Cell Transcriptomic Profiling and Systems Biology of Glioblastoma

**Author:** Abhishek S R  
**Organization:** Sequensolutions  
**Domain:** Bioinformatics, Cancer Genomics, Single-Cell RNA-seq  

## Project Overview

This repository contains a modular, end-to-end bioinformatics pipeline for single-cell RNA sequencing (scRNA-seq) analysis of human Glioblastoma Multiforme (GBM). The dataset consists of 5,604 cells profiled on the 10x Genomics Chromium platform from a 57-year-old male donor and sequenced on an Illumina NovaSeq 6000 at a mean depth of 44,736 reads per cell, with a median of 3,094 detected genes and 9,434 UMIs per cell [file:11]. The project spans from raw 10x count matrices through advanced systems biology, cell–cell communication, pseudotime, and transcription factor activity inference.

The primary goals are to:

- Characterize intra-tumoral heterogeneity within the GBM tumor compartment. [file:5]  
- Map the tumor microenvironment (TME), including astrocytes, microglia/macrophages, OPCs, and T cells. [file:11]  
- Identify tumor-intrinsic oncogenic programs, pathways, and master regulators. [file:8][file:6][file:2]

The dataset and derived materials are distributed under the [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license, and 10x Genomics should be cited following their guidelines: <https://www.10xgenomics.com/support/software/cell-ranger/latest/miscellaneous/cr-citations>.

## Input Data

- **Platform:** 10x Genomics Chromium Single Cell 3’ v3
- **Parent library:** Whole-transcriptome gene expression library used as input for downstream target enrichment
- **Sequencer:** Illumina NovaSeq 6000
- **Configuration:** Paired-end 28 x 91, single-indexed (16 bp barcode, 10 bp UMI, 91 bp transcript)
- **Cell Ranger:** Used for demultiplexing, alignment, and gene–cell matrix generation; web summary and output files are available.

## Pipeline Architecture

The pipeline is implemented in R and organized into 12 sequential scripts under `R/`, using Seurat v5 as the core framework.

### Phase 1: Pre-processing & Core Analytics

1. **00_setup_packages.R**  
   - Installs and loads required CRAN/Bioconductor packages.  
   - Configures Seurat v5 and helper options for reproducibility (e.g., seed, parallelism).

2. **01_load_data_qc.R**  
   - Reads 10x matrices using `Read10X` and constructs a Seurat object.  
   - Performs QC filtering using thresholds on:
     - Total RNA counts per cell,
     - Number of detected features,
     - Percentage of mitochondrial reads.  
   - Generates basic QC plots (UMI vs features, percent MT, etc.).

3. **02_normalization_feature_selection.R**  
   - Applies **SCTransform** to normalize counts, regress out sequencing depth and other technical covariates, and identify highly variable genes.  
   - Stores the SCTransformed assay as the default for downstream steps.

4. **03_dimensionality_reduction.R**  
   - Runs PCA on variable genes and selects principal components based on standard diagnostics (e.g., elbow plots).  
   - Computes UMAP embeddings for nonlinear visualization of the cellular landscape.

5. **04_clustering.R**  
   - Builds a shared nearest neighbor (SNN) graph in the SCTransformed space.  
   - Performs graph-based clustering (Louvain algorithm) across resolutions, identifying discrete transcriptional clusters.

6. **05_annotation.R**  
   - Annotates clusters into biologically meaningful cell types using canonical markers (e.g. astrocytic, immune, oligodendroglial, and tumor markers).  
   - Summary of inferred cell-type composition (counts and percentages) is stored in `05_cell_type_summary.csv`. [file:11]

### Phase 2: Differential Expression & Functional Enrichment

7. **06_differential_expression.R**  
   - Focuses on the tumor compartment versus astrocytes as a reference glial population.  
   - Performs differential expression using Seurat’s DE framework to derive:
     - **Full tumor vs astrocyte DE table** (`06_Tumor_vs_Astrocyte_DGE.csv`) with log fold changes, adjusted p-values, and expression metrics. [file:8]  
     - **Tumor-upregulated genes** (`06_Tumor_Upregulated_Only.csv`), representing candidate oncogenic drivers and GBM markers. [file:9]  
     - **Tumor-downregulated genes** (`06_Tumor_Downregulated_Only.csv`), indicating glial programs that are lost or suppressed in tumor cells. [file:10]

8. **07_pathway_enrichment.R**  
   - Uses **clusterProfiler** with **org.Hs.eg.db** to perform over-representation and/or GSEA-based enrichment.  
   - Enrichment is run on tumor-upregulated signatures against:
     - **Gene Ontology (GO)** – biological processes and molecular functions (`07_GO_Enrichment_Results.csv`). [file:7]  
     - **KEGG pathways** – signaling and metabolic circuits (`07_KEGG_Enrichment_Results.csv`). [file:6]  
   - Results are visualized with dotplots, barplots, and network-style plots (via enrichplot), highlighting hallmark GBM-associated pathways and microenvironment interactions.

### Phase 3: Advanced Systems Biology & Tumor Heterogeneity

9. **08_tumor_specific_analysis.R**  
   - Subsets the Seurat object to tumor cells only.  
   - Performs cell cycle scoring (S-phase, G2/M) and annotates cycling vs non-cycling populations. [file:5]  
   - Re-clusters tumor cells at higher resolution to reveal sub-clonal architecture and intra-tumoral heterogeneity, storing metadata in `08_tumor_metadata_annotated.csv`. [file:5]  

10. **09_cell_cell_communication.R**  
    - Converts Seurat objects into CellChat inputs, handling Seurat v5 layers and assays to avoid matrix extraction issues.  
    - Infers ligand–receptor-mediated signaling programs between:
      - Tumor sub-clones,
      - Astrocytes,
      - Microglia/macrophages,
      - OPCs and T cells. [file:11]  
    - Exports an inferred communication network in `09_Inferred_Communication_Network.csv`, including source and target cell types and pathway-level signaling strength. [file:4]

11. **10_pseudotime.R**  
    - Constructs a **SingleCellExperiment** from the tumor subset, embedding it in the previously computed low-dimensional space.  
    - Uses **slingshot** to infer lineage trajectories across tumor sub-clones, representing putative evolutionary or phenotypic paths.  
    - Stores pseudotime values alongside metadata in `10_tumor_metadata_with_pseudotime.csv`. [file:3]

12. **11_Transcription_factor_activity.R**  
    - Integrates **DoRothEA** regulons with the **viper** algorithm to compute transcription factor (TF) activity scores from single-cell expression profiles.  
    - Tailors the workflow to work around Seurat v5’s layered architecture, ensuring expression matrices are correctly passed into viper.  
    - Outputs TF activity scores for tumor cells, including candidate “master regulators”, in `11_Tumor_TF_Activity_Scores.csv`. [file:2]

## Technologies & Dependencies

- **Language:** R (≥ 4.2 recommended)
- **Core single-cell framework:** Seurat v5
- **Systems biology and signaling:**
  - CellChat
  - DoRothEA
  - viper
- **Trajectory inference:**
  - slingshot
  - SingleCellExperiment
- **Functional enrichment:**
  - clusterProfiler
  - org.Hs.eg.db
- **Visualization:**
  - ggplot2
  - patchwork
  - enrichplot

Additional packages include standard data handling and plotting libraries as configured in `00_setup_packages.R`.

## Installation and Usage

1. Clone the repository:

   ```bash
   git clone https://github.com/<your-username>/SingleCell_GBM_Pipeline.git
   cd SingleCell_GBM_Pipeline
   ```

2. Open R or RStudio in this directory.

3. Run the setup script:

   ```r
   source("R/00_setup_packages.R")
   ```

4. Place your 10x Genomics count matrices under `data/raw/` and adjust paths in `01_load_data_qc.R` if necessary.

5. Run the pipeline stepwise, or create a wrapper script to run all stages:

   ```r
   source("R/01_load_data_qc.R")
   source("R/02_normalization_feature_selection.R")
   source("R/03_dimensionality_reduction.R")
   source("R/04_clustering.R")
   source("R/05_annotation.R")
   source("R/06_differential_expression.R")
   source("R/07_pathway_enrichment.R")
   source("R/08_tumor_specific_analysis.R")
   source("R/09_cell_cell_communication.R")
   source("R/10_pseudotime.R")
   source("R/11_Transcription_factor_activity.R")
   ```

Intermediate objects, metadata tables, and figures are written to `data/metadata/` and `results/` as shown in the repository tree.

## Key Technical Highlights

- Seamless compatibility with **Seurat v5** layer structure for legacy systems biology tools such as CellChat and viper.  
- Modular scripts allow targeted re‑runs (e.g., only the tumor-specific analysis or enrichment) without recomputing earlier steps.  
- Outputs include cluster markers, cell-type proportions, tumor vs astrocyte DE, enriched pathways, tumor sub-cluster annotations, communication networks, pseudotime, and TF activity scores. [file:11][file:8][file:6][file:4][file:3][file:2]

## Licensing and Attribution

This work uses a 10x Genomics GBM dataset distributed under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) license. Please follow 10x Genomics’ citation recommendations for any publications built on this code or dataset.

Unless otherwise noted, the code in this repository is released under the MIT License (see `LICENSE`).
