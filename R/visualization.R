#' @import Seurat
#' @import RColorBrewer
#' @export

get_umap <- function(h5, dims, res){
  h5 <- FindNeighbors(h5, dims = 1:dims)
  h5 <- FindClusters(h5, resolution = res)
  h5 <- RunUMAP(h5, dims = 1:dims)
  DimPlot(h5, reduction = "umap")
}

