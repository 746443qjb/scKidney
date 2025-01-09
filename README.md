# scKidney
scKidney is an R package designed for kidney single-cell RNA-seq data analysis. It provides tools for data quality control (QC), principal component selection (PCA), doublet detection, cell annotation, and data format conversion (Seurat ↔ H5AD).
## Installation and use it
Run the following command to install the scKidney package:
```R
devtools::install_github("746443qjb/scKidney")
library(scKidney)
# Perform data quality control
SCE <- KidneyQC(seurat = your_seurat_object, species = "human", Plot = TRUE)
# Select optimal principal components
pcs <- KidneyPC(seurat = your_seurat_object, cum = 90, var = 5)
# Detect and filter doublets
seurat_filtered <- KidneyDoublet(seurat = your_seurat_object, PC = pcs, rate = 8, select = "high")
# Annotate cells
seurat_annotated <- KidneyCell(seurat = your_seurat_object, species = "human", method = "Aucell", class = 1, plot = "heatmap", addcelltype = TRUE)
# Convert Seurat to H5AD(Note: It is recommended to run joinlayers on Seurat v5 objects before using this function to avoid any issues during conversion)
library(reticulate)
anndata <- import("anndata")
np <- import("numpy")
kidneyH5(seurat_object, "sce.h5ad")
#Convert H5AD to Seurat
sce=scSeurat(h5ad_path)




