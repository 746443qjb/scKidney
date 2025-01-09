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

  library(reticulate)
  library(Seurat)
  library(Matrix)


  anndata <- import("anndata", convert = FALSE)
  adata <- anndata$read_h5ad(h5ad_path)


  use_raw <- FALSE
  var_names <- NULL
  if (!py_is_null_xptr(adata$raw)) {  # 检查 adata$raw 是否为空
    tryCatch({
      if (!is.null(adata$raw$X)) {  # 检查 raw$X 是否存在
        use_raw <- TRUE

        var_names <- as.character(reticulate::py_to_r(adata$raw$var$index$to_list()))
      }
    }, error = function(e) {
      use_raw <- FALSE
    })
  }

  if (use_raw) {
    message("Using raw counts matrix from adata$raw$X")
    counts <- adata$raw$X
  } else {
    message("Using normalized counts matrix from adata$X")
    counts <- adata$X

    var_names <- as.character(reticulate::py_to_r(adata$var$index$to_list()))
  }


  counts <- reticulate::py_to_r(counts)


  if (inherits(counts, "dgRMatrix") || inherits(counts, "dgCMatrix")) {

    message("Counts matrix is already a sparse matrix")
  } else {

    counts <- as(Matrix::sparseMatrix(
      i = counts@i + 1,
      p = counts@p,
      x = counts@x,
      dims = c(counts@shape[[1]], counts@shape[[2]])
    ), "dgCMatrix")
  }


  counts <- t(counts)


  obs_names <- as.character(reticulate::py_to_r(adata$obs$index$to_list()))


  if (length(var_names) != nrow(counts)) {
    stop("Number of gene names does not match the number of rows in counts")
  }
  if (length(obs_names) != ncol(counts)) {
    stop("Number of cell names does not match the number of columns in counts")
  }


  rownames(counts) <- var_names
  colnames(counts) <- obs_names


  meta_data <- reticulate::py_to_r(adata$obs)
  rownames(meta_data) <- obs_names  # 确保元数据行名与细胞名称一致


  seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = meta_data,
    assay = "RNA"
  )


  obsm_dict <- reticulate::py_to_r(dict(adata$obsm))
  for (key in names(obsm_dict)) {

    embedding <- obsm_dict[[key]]


    rownames(embedding) <- colnames(seurat_obj)
    colnames(embedding) <- paste0(key, "_", 1:ncol(embedding))


    fixed_key <- gsub("^X_", "", key)
    fixed_key <- paste0(fixed_key, "_")


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


