#' KidneyH5
#'
#' This function converts a Seurat object to an H5AD file format.
#'
#' @param seurat_obj A Seurat object to be converted.
#' @param output_path Path where the H5AD file will be saved.
#' @details This function requires the installation of the reticulate package and the Python anndata library.
#' It is recommended to run `joinlayers` on Seurat v5 objects before using this function.
#' @export
kidneyH5 <-  function(seurat_obj, output_path) {
  # 获取 counts 矩阵
  counts_matrix <- GetAssayData(seurat_obj, assay = 'RNA', slot = 'counts')

  # 确保 counts_matrix 是矩阵类型
  if (!inherits(counts_matrix, "matrix") && !inherits(counts_matrix, "dgCMatrix")) {
    counts_matrix <- as.matrix(counts_matrix)
  }

  # 转置 counts 矩阵使其为细胞数 x 基因数
  counts_matrix <- t(counts_matrix)

  # 确保 counts_matrix 是 dgCMatrix 格式的稀疏矩阵
  counts_matrix <- as(counts_matrix, "dgCMatrix")

  # 获取 metadata
  meta_data <- seurat_obj@meta.data
  meta_data$barcode <- rownames(meta_data)

  # 创建 AnnData 对象
  adata <- anndata$AnnData(X = counts_matrix, obs = pd$DataFrame(meta_data))

  # 设置基因名称
  gene_names <- colnames(counts_matrix)
  adata$var_names <- np$array(gene_names)

  # 添加 PCA 降维信息
  if ("pca" %in% names(seurat_obj@reductions)) {
    pca_embeddings <- Embeddings(seurat_obj, reduction = "pca")
    pca_df <- pd$DataFrame(data = pca_embeddings, index = adata$obs_names)
    adata$obsm["X_pca"] <- np$array(pca_embeddings)
  }

  # 添加 Harmony 降维信息
  if ("harmony" %in% names(seurat_obj@reductions)) {
    harmony_embeddings <- Embeddings(seurat_obj, reduction = "harmony")
    harmony_df <- pd$DataFrame(data = harmony_embeddings, index = adata$obs_names)
    adata$obsm["X_harmony"] <- np$array(harmony_embeddings)
  }

  # 保存为 H5AD 文件
  adata$write_h5ad(output_path)
}

# 示例用法
# sce <- kidneyH5(sce,"sce.h5ad")
