#' KidneyPC
#'
#' This function determines the optimal number of Principal Components (PCs) to use in downstream analysis of a Seurat object.
#'
#' @param seurat A Seurat object that has been normalized, scaled, and processed with PCA or other dimensionality reduction methods.
#' @param reduction The dimensionality reduction method to use (e.g., "pca", "harmony"). Must exist in the Seurat object.
#' @param cum The cumulative percentage threshold to use for selecting PCs.
#' @param var The minimum variation percentage for selecting PCs.
#' @return The number of PCs that meet the given criteria.
#' @export
KidneyPC <- function(seurat, reduction = "pca", cum = 90, var = 5) {

  # 参数检查
  if (missing(seurat)) {
    stop("You must provide a Seurat object.")
  }
  if (!reduction %in% names(seurat@reductions)) {
    stop("The specified reduction method does not exist in the Seurat object.")
  }

  # 计算每个 PC 的标准差和累积百分比
  pct <- seurat[[reduction]]@stdev / sum(seurat[[reduction]]@stdev) * 100
  cumu <- cumsum(pct)

  # 找到累积百分比超过 cum% 且对应 PC 的变异度小于 var% 的第一个 PC
  co1 <- which(cumu > cum & pct < var)[1]

  # 找到前后两个 PC 的变异度差异超过 0.1% 的位置
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1

  # 确定最终的最佳 PC 数量
  pcs <- min(co1, co2, na.rm = TRUE)

  # 绘图部分
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

  p <- ggplot(plot_df, aes(x = cumu, y = pct, label = rank, color = rank > pcs)) +
    geom_text(check_overlap = TRUE, size = 3.5, fontface = "bold") +
    geom_vline(xintercept = cum, color = "grey", linetype = "dashed") +
    geom_hline(yintercept = min(pct[pct > var], na.rm = TRUE), color = "grey", linetype = "dashed") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "skyblue")) +
    labs(
      title = "Optimal Principal Components Selection",
      x = "Cumulative Percentage",
      y = "Percentage of Variation"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  print(p) # 确保绘图在函数内部被正确显示

  return(pcs)
}

# 示例用法：
# pcs <- KidneyPC(seurat = your_seurat_object, reduction = "pca", cum = 90, var = 5)


