#' scSeurat
#'
#' This function converts an AnnData object stored in an H5AD file into a Seurat object.
#' It supports loading raw counts (if available) or normalized counts and includes metadata
#' and dimensionality reduction embeddings.
#'
#' @param h5ad_path A string specifying the file path to the H5AD file.
#' @return A Seurat object with the following components:
#'   \item{counts}{A sparse matrix containing raw or normalized counts.}
#'   \item{meta.data}{A data frame containing cell metadata.}
#'   \item{dimensional reductions}{Dimensionality reduction embeddings (e.g., PCA, UMAP, t-SNE).}
#' @details The function first attempts to load raw counts from `adata$raw$X`. If raw counts are not available,
#' it falls back to normalized counts from `adata$X`. The function also extracts cell metadata from `adata$obs`,
#' gene names from `adata$var` (or `adata$raw$var` if raw counts are used), and dimensional reduction embeddings from `adata$obsm`.
#'
#' The gene and cell names are added to the counts matrix to ensure compatibility with Seurat.
#' Dimensional reduction embeddings are added to the Seurat object under their respective names.
#'
#' @export
scSeurat <- function(h5ad_path) {
  # 加载必要的库
  library(reticulate)
  library(Seurat)
  library(Matrix)

  # 加载 AnnData 对象
  anndata <- import("anndata", convert = FALSE)
  adata <- anndata$read_h5ad(h5ad_path)

  # 检查并优先使用 adata$raw$X
  use_raw <- FALSE
  var_names <- NULL  # 基因名称
  if (!py_is_null_xptr(adata$raw)) {  # 检查 adata$raw 是否为空
    tryCatch({
      if (!is.null(adata$raw$X)) {  # 检查 raw$X 是否存在
        use_raw <- TRUE
        # 提取 raw 的基因名称
        var_names <- as.character(reticulate::py_to_r(adata$raw$var$index$to_list()))
      }
    }, error = function(e) {
      use_raw <- FALSE  # 如果访问 adata$raw$X 报错，回退到 adata$X
    })
  }

  if (use_raw) {
    message("Using raw counts matrix from adata$raw$X")
    counts <- adata$raw$X
  } else {
    message("Using normalized counts matrix from adata$X")
    counts <- adata$X
    # 如果未使用 raw，则提取标准化计数矩阵的基因名称
    var_names <- as.character(reticulate::py_to_r(adata$var$index$to_list()))
  }

  # 转换计数矩阵为 R 对象
  counts <- reticulate::py_to_r(counts)

  # 检查 counts 是否已经是 R 的稀疏矩阵格式
  if (inherits(counts, "dgRMatrix") || inherits(counts, "dgCMatrix")) {
    # 如果 counts 是稀疏矩阵，直接使用
    message("Counts matrix is already a sparse matrix")
  } else {
    # 如果是 SciPy 稀疏矩阵，转换为 dgCMatrix
    counts <- as(Matrix::sparseMatrix(
      i = counts@i + 1,
      p = counts@p,
      x = counts@x,
      dims = c(counts@shape[[1]], counts@shape[[2]])
    ), "dgCMatrix")
  }

  # 转置计数矩阵（确保行是基因，列是细胞）
  counts <- t(counts)

  # 提取细胞名称
  obs_names <- as.character(reticulate::py_to_r(adata$obs$index$to_list()))

  # 确保基因名称和细胞名称与矩阵的维度匹配
  if (length(var_names) != nrow(counts)) {
    stop("Number of gene names does not match the number of rows in counts")
  }
  if (length(obs_names) != ncol(counts)) {
    stop("Number of cell names does not match the number of columns in counts")
  }

  # 为计数矩阵设置行名和列名
  rownames(counts) <- var_names
  colnames(counts) <- obs_names

  # 提取元数据（adata$obs）
  meta_data <- reticulate::py_to_r(adata$obs)
  rownames(meta_data) <- obs_names  # 确保元数据行名与细胞名称一致

  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = meta_data,
    assay = "RNA"
  )

  # 提取降维数据（adata$obsm）
  obsm_dict <- reticulate::py_to_r(dict(adata$obsm))  # 替换 as_dict 为 dict
  for (key in names(obsm_dict)) {
    # 提取降维数据
    embedding <- obsm_dict[[key]]

    # 设置降维数据的行名和列名
    rownames(embedding) <- colnames(seurat_obj)  # 细胞名称
    colnames(embedding) <- paste0(key, "_", 1:ncol(embedding))  # 特征名称

    # 修复降维 key 名称（去掉 'X_' 前缀并添加下划线）
    fixed_key <- gsub("^X_", "", key)
    fixed_key <- paste0(fixed_key, "_")

    # 创建降维对象并添加到 Seurat 对象
    seurat_obj[[fixed_key]] <- CreateDimReducObject(
      embeddings = embedding,
      key = fixed_key,
      assay = "RNA"
    )
  }

  # 返回 Seurat 对象
  return(seurat_obj)
}
# 示例用法：
# seurat <- scSeurat("sce.h5ad")


