#' Create Matrix Heatmap.
#'
#' Create a heatmap for the provided matrix. This is essentially just a wrapper for
#' \code{pheatmap}.
#'
#' @param mat the matrix for plotting
#' @param title a string for the plot title
#' @param cluster a boolean indicating whether to perform hierarchical clustering
#' @param names a boolean indicating whether to display column/row names
#' @return a heatmap
#' @export

makeHeatmap <- function(mat, title, cluster = T, names = F) {
  heatmap <- pheatmap::pheatmap(mat,
                                cluster_rows = cluster,
                                cluster_cols = cluster,
                                main = title,
                                grDevices::colorRampPalette(c("navy", "white", "red"))(100),
                                treeheight_row = 0, treeheight_col = 0,
                                show_rownames = names, show_colnames = names,
                                breaks = seq(-1, 1, length.out = 100),
                                fontsize_row = 5, fontsize_col = 5)
  return(heatmap)
}
