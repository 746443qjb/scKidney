run the following command to install scKidney:
devtools::install_github("746443qjb/scKidney")

QC:SCE <- KidneyQC(seurat = your_seurat_object, species = "human", Plot = TRUE)

Find the best PCA: pcs <- KidneyPC(seurat = your_seurat_object, cum = 90, var = 5)

Find double cells:seurat_filtered <- KidneyDoublet(seurat = your_seurat_object, PC = pcs, rate = 8, select = "high")

cell annotation: seurat_annotated <- KidneyCell(seurat = your_seurat_object, species = "human", method = "Aucell", class = 1, plot = "heatmap", addcelltype = TRUE)
