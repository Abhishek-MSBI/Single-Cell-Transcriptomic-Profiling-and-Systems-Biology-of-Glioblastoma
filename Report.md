# Single-Cell Transcriptomic Profiling and Systems Biology of Glioblastoma

**Author:** Abhishek S R (Sequensolutions)  

## 1. Dataset and Experimental Design

Glioblastoma Multiforme (GBM) cells from a 57‑year‑old male donor were profiled using the 10x Genomics Chromium Single Cell 3' v3 chemistry. Cell Ranger preprocessing produced a filtered gene–cell matrix with 5,604 high-quality cells, a median of 3,094 detected genes, and 9,434 UMIs per cell, indicating deep coverage at the single-cell level suitable for downstream systems biology analysis. [file:11] Sequencing was performed on an Illumina NovaSeq 6000 with approximately 44,736 reads per cell, providing sufficient depth to detect both housekeeping and lower-abundance regulatory transcripts.

## 2. Quality Control and Pre-processing

Standard QC metrics were computed, including total RNA counts, number of detected features, and mitochondrial read percentages per cell, and low-quality cells were removed based on empirically defined cutoffs. This filtering step enriched for cells with robust transcriptional profiles while removing empty droplets and stressed or dying cells. SCTransform normalization was then applied to regress out technical variation and to derive a set of highly variable genes for dimensionality reduction and clustering.

## 3. Cell Type Composition and Cluster Annotation

Graph-based clustering on the SCTransformed space followed by canonical marker-based annotation revealed major cell compartments within the tumor and its microenvironment. The cell type summary indicates that the sample is dominated by astrocytes (3,921 cells, 41.08%) and tumor cells (3,885 cells, 40.31%), with additional contributions from microglia/macrophages (476 cells, 10.18%), oligodendrocyte progenitor cells (252 cells, 5.39%), and T cells (142 cells, 3.04%). [file:11] This distribution reflects a typical GBM ecosystem characterized by a large malignant astroglial-like compartment embedded in an immune-rich microenvironment.

## 4. Tumor vs Astrocyte Differential Expression

To delineate tumor-intrinsic programs from normal glial biology, differential expression analysis was performed between the tumor compartment and astrocyte population. The global tumor vs astrocyte contrast is summarized in `06_Tumor_vs_Astrocyte_DGE.csv`, while tumor-specific upregulated and downregulated genes are provided in `06_Tumor_Upregulated_Only.csv` and `06_Tumor_Downregulated_Only.csv`. [file:8][file:9][file:10] Tumor-upregulated genes capture candidate oncogenic effectors, proliferative and stemness-associated genes, and pathways associated with invasion and treatment resistance, whereas tumor-downregulated genes highlight astrocytic homeostatic programs lost in malignant transformation.

## 5. Functional Enrichment of Tumor Programs

Tumor-upregulated genes were subjected to Gene Ontology and KEGG pathway enrichment using clusterProfiler. GO results in `07_GO_Enrichment_Results.csv` demonstrate that tumor genes aggregate into coherent biological processes and molecular functions such as cell cycle regulation, DNA replication, and mitotic spindle organization, consistent with a highly proliferative GBM phenotype. [file:7] KEGG enrichment in `07_KEGG_Enrichment_Results.csv` further reveals tumor association with canonical oncogenic pathways and signaling circuits that underlie GBM biology and therapeutic vulnerability. [file:6]

## 6. Tumor Heterogeneity and Sub-clonal Architecture

The tumor compartment was isolated and reanalyzed to dissect intra-tumoral heterogeneity. Cell cycle scoring for S and G2/M phases, followed by reclustering, identified proliferative versus quiescent tumor subpopulations and finer-grained sub-clusters with distinct transcriptional states, as captured in `08_tumor_metadata_annotated.csv`. [file:5] This stratification provides a basis for linking sub-clones to specific pathways, microenvironmental niches, or inferred evolutionary trajectories.

## 7. Cell–Cell Communication in the Tumor Microenvironment

CellChat-based ligand–receptor analysis was used to model intercellular communication between tumor sub-clones and non-malignant compartments, including astrocytes, microglia/macrophages, OPCs, and T cells. The resulting communication network, exported in `09_Inferred_Communication_Network.csv`, summarizes signaling edges between sender and receiver cell types and the activity of specific signaling pathways. [file:4] These results offer a systems-level view of how malignant and stromal cells may cross-talk through growth factors, cytokines, and immune regulatory axes within the GBM microenvironment.

## 8. Trajectory Inference and Pseudotime

To reconstruct putative evolutionary or phenotypic trajectories within the tumor, the re-clustered tumor subset was embedded into a SingleCellExperiment object and analyzed using slingshot. Pseudotime values were inferred along lineages spanning tumor sub-clusters, and integrated with cluster and cell cycle metadata in `10_tumor_metadata_with_pseudotime.csv`. [file:3] These pseudotime trajectories help distinguish early-like and late-like tumor states and suggest potential transition paths that may correspond to dedifferentiation, lineage switching, or therapy-induced adaptations.

## 9. Transcription Factor Activity and Master Regulators

Finally, transcription factor activity was inferred using DoRothEA regulons and the viper algorithm, applied to tumor expression data extracted from the Seurat v5 object. TF activity scores summarized in `11_Tumor_TF_Activity_Scores.csv` highlight candidate master regulators with strongly perturbed activity, which are likely to orchestrate the GBM-specific oncogenic programs uncovered in earlier analyses. [file:2] These regulators provide mechanistic hypotheses for targeted perturbation, drug repurposing, or combinatorial therapy design.

## 10. Reproducibility and Extensibility

The analysis is fully scripted in R, with each processing stage encapsulated in a dedicated script under `R/`, producing standardized intermediate objects and tabular outputs. By design, the pipeline can be rerun end-to-end or partially, allowing users to swap in different QC thresholds, clustering resolutions, or enrichment databases while preserving the overall workflow. The modular architecture and reliance on widely used packages (Seurat v5, CellChat, slingshot, clusterProfiler, DoRothEA, viper) facilitate adaptation to other GBM samples and to broader single-cell oncology projects.
