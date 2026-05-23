packages <- c(
  "Seurat", "SeuratObject", "ggplot2", "dplyr", "patchwork",
  "Matrix", "future", "scales", "RColorBrewer", "SingleR",
  "celldex", "clusterProfiler", "org.Hs.eg.db", "msigdbr",
  "fgsea", "infercnv", "CellChat", "monocle3", "slingshot",
  "MatrixGenerics", "scater", "scran", "DoubletFinder"
)

installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

invisible(lapply(packages, library, character.only = TRUE))

cran_pkgs <- c(
  "Seurat", "SeuratObject", "ggplot2", "dplyr", "patchwork",
  "Matrix", "future", "scales", "RColorBrewer",
  "clusterProfiler", "msigdbr", "fgsea"
)

bioc_pkgs <- c("org.Hs.eg.db")
