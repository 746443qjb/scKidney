#' KidneyH5
#'
#' This function converts a Seurat object to an H5AD file format.
#'
#' @param seurat_obj A Seurat object to be converted.
#' @param output_path A character string specifying the path where the H5AD file will be saved.
#' @details This function requires the installation of the reticulate package and the Python anndata library.
#' It is recommended to run `joinlayers()` on Seurat v5 objects before using this function to ensure proper compatibility.
#' This function supports dynamic addition of dimensionality reduction embeddings stored in the Seurat object.
#' @return No return value. The H5AD file will be saved at the specified location.
#' @export
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(KidneyH5)
#'
#' # Example usage
#' sce <- CreateSeuratObject(counts = matrix(runif(10000), nrow = 100, ncol = 100))
#' sce <- SCTransform(sce)
#' kidneyH5(sce, "sce.h5ad")
#' }
#'
kidneyH5 <- function(seurat_obj, output_path) {
  # Load required Python libraries
  anndata <- reticulate::import("anndata", delay_load = TRUE)
  np <- reticulate::import("numpy", delay_load = TRUE)
  pd <- reticulate::import("pandas", delay_load = TRUE)

  counts_matrix <- Seurat::GetAssayData(seurat_obj, assay = 'RNA', layer = 'counts')

  if (!inherits(counts_matrix, "matrix") && !inherits(counts_matrix, "dgCMatrix")) {
    counts_matrix <- as.matrix(counts_matrix)
  }

  counts_matrix <- Matrix::t(counts_matrix)

  counts_matrix <- as(counts_matrix, "dgCMatrix")

  meta_data <- seurat_obj@meta.data
  meta_data$barcode <- rownames(meta_data)
  adata <- anndata$AnnData(X = counts_matrix, obs = pd$DataFrame(meta_data))
  gene_names <- colnames(counts_matrix)
  adata$var_names <- np$array(gene_names)

  if (length(seurat_obj@reductions) > 0) {
    for (reduction_name in names(seurat_obj@reductions)) {
      embeddings <- Seurat::Embeddings(seurat_obj, reduction = reduction_name)
      obsm_key <- paste0("X_", reduction_name)
      adata$obsm[obsm_key] <- np$array(embeddings)
    }
  }
  adata$write_h5ad(output_path)
}
# 示例用法：
#kidneyH5(sce, "sce.h5ad")


