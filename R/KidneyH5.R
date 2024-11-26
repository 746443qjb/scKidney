#' KidneyH5
#'
#' This function converts a Seurat object to an H5AD file format.
#'
#' @param seurat_obj A Seurat object to be converted.
#' @param output_path Path where the H5AD file will be saved.
#' @details This function requires the installation of the reticulate package and the Python anndata library.
#' It is recommended to run `joinlayers` on Seurat v5 objects before using this function.
#' @export
kidneyH5 <- function(seurat_obj, output_path) {
  # 获取counts矩阵并转置
  counts_matrix <- t(GetAssayData(seurat_obj, assay = 'RNA', slot = 'counts'))

  # 获取metadata
  meta_data <- seurat_obj@meta.data
  meta_data$barcode <- rownames(meta_data)  # 添加条形码作为索引

  # 确保 metadata 的行名和 counts 矩阵的行名一致
  meta_data <- meta_data[rownames(counts_matrix), , drop = FALSE]

  # 检查并确保没有 NA 行
  if (any(is.na(rownames(meta_data)))) {
    stop("Metadata contains NA values after reordering. Please check the Seurat object.")
  }

  # 转换 counts 矩阵为稀疏矩阵格式
  counts_csr <- as(counts_matrix, "CsparseMatrix")

  # 创建 AnnData 对象
  adata <- anndata$AnnData(X = counts_csr, obs = meta_data)

  # 添加 var（基因名称）
  gene_names <- data.frame('gene' = colnames(counts_matrix))
  adata$var <- data.frame(index = colnames(counts_matrix))
  adata$var_names <- np$array(colnames(counts_matrix))

  # 遍历所有降维结果并添加到 AnnData 对象中
  for (reduction_name in names(seurat_obj@reductions)) {
    embeddings <- Embeddings(seurat_obj, reduction = reduction_name)
    if (nrow(embeddings) == nrow(counts_matrix)) {
      embeddings <- embeddings[rownames(counts_matrix), , drop = FALSE]  # 确保降维结果行与细胞匹配
      adata$obsm[[paste0("X_", tolower(reduction_name))]] <- np$array(embeddings)
    } else {
      warning(paste("Reduction", reduction_name, "has a different number of cells than the counts matrix. Skipping."))
    }
  }

  # 保存为 H5AD 文件
  adata$write_h5ad(output_path)
}
# 示例用法
# sce <- kidneyH5(sce,"sce.h5ad")
