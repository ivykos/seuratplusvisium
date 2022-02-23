#' @import Seurat
#' @import RColorBrewer
#' @import ggplot2
#' @import viridis
#' @export

get_umap <- function(h5, dims, res){
  h5 <- FindNeighbors(h5, dims = 1:dims)
  h5 <- FindClusters(h5, resolution = res)
  h5 <- RunUMAP(h5, dims = 1:dims)
  DimPlot(h5, reduction = "umap")
}

get_visium <- function(){
  
  #Table maniupluation
  #Etc. etc. 
  
  ggplot(pos.ordered, aes(x=pos.ordered$V3, y=pos.ordered$V4, color=as.factor(table.a$V2)))
  + geom_point()
}

get_expression <- function(){
  tmp <- as.matrix(GetAssayData(object = caf_a, slot = "counts"))
  
  
}

