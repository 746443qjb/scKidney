#' KidneyCell
#'
#' This function annotates cell types in a Seurat object based on marker gene scoring, currently supporting human samples only.
#'
#' @param seurat A Seurat object that has been processed with FindNeighbors and FindClusters.
#' @param species Either 'human' or 'mouse', currently only 'human' is supported.
#' @param method The scoring method to use, either 'Aucell' or 'Vision', default is 'Aucell'.
#' @param class Either 1 or 2, specifying the level of cell type identification. 1 for preliminary, 2 for detailed annotation.
#' @param plot The type of plot to visualize the results, either 'heatmap' or 'bubble'.
#' @param addcelltype Logical value, whether to add a 'celltype' metadata to clusters, default is TRUE.
#' @param colors A vector of colors to use for plotting, default is c("#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#D65076").
#' @return A Seurat object with cell type annotations.
#' @export
KidneyCell <- function(seurat, species = "human", method = "Aucell", class = 1, plot = "heatmap", addcelltype = TRUE, colors = c("#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#D65076")) {
  # 加载 marker 列表
  markers_file <- system.file("extdata", "markers_list.rds", package = "YourPackageName")
  markers_list <- readRDS(markers_file)

  # 选择 markers
  if (class == 1) {
    markers <- markers_list$class_1
  } else if (class == 2) {
    markers <- markers_list$class_2
  }

  # 评分计算部分
  if (method == "Aucell") {
    library(AUCell)
    gene_sets <- lapply(markers, function(marker) as.character(marker))
    cells_rankings <- AUCell_buildRankings(seurat@assays$RNA@data, plotStats = FALSE)
    auc <- AUCell_calcAUC(gene_sets, cells_rankings)
    seurat <- AddMetaData(seurat, as.data.frame(t(auc@assays$AUC)), col.name = "celltype_scores")
  } else if (method == "Vision") {
    library(VISION)
    gene_sets <- lapply(markers, function(marker) Vision::GeneSet(marker, name = names(marker)))
    vision_obj <- Vision(seurat, gene_sets)
    vision_obj <- analyze(vision_obj)
    seurat <- AddMetaData(seurat, vision_obj@MetaData$SignatureScores, col.name = "celltype_scores")
  }

  if (addcelltype) {
    max_scores <- apply(seurat@meta.data$celltype_scores, 1, which.max)
    seurat$celltype <- names(markers)[max_scores]
  }

  if (plot == "heatmap") {
    p <- DoHeatmap(seurat, features = unlist(markers), group.by = "celltype") + scale_fill_manual(values = colors)
  } else if (plot == "bubble") {
    p <- DotPlot(seurat, features = unlist(markers), group.by = "celltype") + scale_color_manual(values = colors)
  }
  print(p)

  return(seurat)
}

# 示例用法：
# seurat_annotated <- KidneyCell(seurat = your_seurat_object, species = "human", method = "Aucell", class = 1, plot = "heatmap", addcelltype = TRUE)

