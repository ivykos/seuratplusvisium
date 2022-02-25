#' @import hdf5r
#' @import Seurat
#' @import ggplot2
#' @export

read_h5 <- function(h5, proj, a){
  h5_file <- Read10X_h5(h5, use.names = FALSE) 
  h5.seurat <- CreateSeuratObject(counts = h5_file, project = proj, assay = a)
  h5.seurat <- NormalizeData(h5.seurat)
  h5.seurat <- FindVariableFeatures(h5.seurat, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(h5.seurat)
  h5.seurat <- ScaleData(h5.seurat, features = all.genes)
  h5.seurat <- RunPCA(h5.seurat, features = VariableFeatures(object = h5.seurat))
 
}

#' @export
get_umap <- function(seuratObj, dimensions, res){
  seuratObj <- FindNeighbors(seuratObj, dims = 1:dimensions)
  seuratObj <- FindClusters(seuratObj, resolution = res)
  
  seuratObj <- RunUMAP(seuratObj, dims = 1:dimensions)
  DimPlot(seuratObj, reduction = "umap")
}

#' @export
transfer_clusters <- function(seuratObj, proj, tissue_csv){
  write.table(seuratObj@active.ident, file="tmp.tsv", quote=FALSE, sep="\t", col.names = FALSE)
  idents <- read.delim("tmp.tsv", header = FALSE)
  file.remove("tmp.tsv")
  write.table(seuratObj@meta.data$orig.ident, file="tmp.tsv", quote=FALSE, sep="\t", col.names = FALSE)
  orig <- read.delim("tmp.tsv", header = FALSE)
  file.remove("tmp.tsv")
  
  positions <- read.csv(tissue_csv, header = F)
  positions <- positions[positions$V2 == 1,]
  
  idents["Origin"] = orig$V2
  table <- idents[idents$Origin == proj,]
  cluster.ordered <- table[order(table$V1),]
  pos.ordered <-positions[order(positions$V1),]
  
  ggplot(pos.ordered, aes(x=pos.ordered$V3, y=pos.ordered$V4, color=as.factor(table$V1))) + 
    geom_point() + theme_linedraw() +ggtitle(as.character(proj))
}
