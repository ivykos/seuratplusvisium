#' @import hdf5r
#' @import Seurat
#' @import ggplot2
#' @export

read_h5 <- function(infile, proj, a){
  h5 <- Read10X_h5(infile, use.names = TRUE) 
  h5.seurat <- CreateSeuratObject(counts = h5, project = proj, assay = a)
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
transfer_clusters <- function(seuratObj, pro, tissue_csv){
  write.table(sueratObj@active.ident, file="tmp.tsv", quote=FALSE, sep="\t", col.names = FALSE)
  idents <- read.delim("tmp.tsv", header = FALSE)
  file.remove("tmp.tsv")
  write.table(combined2@meta.data$orig.ident, file="tmp.tsv", quote=FALSE, sep="\t", col.names = FALSE)
  orig <- read.delim("tmp.tsv", header = FALSE)
  file.remove("tmp.tsv")
  
  positions <- read.csv(tissue_csv, header = F)
  positions <- positions[positions$V2 == 1,]
  
  idents["Origin"] = orig$V2
  table <- idents[idents$Origin == "CAF_A",]
  cluster.ordered <- table[order(table$V1),]
  pos.ordered <-positions[order(positions$V1),]
  
  ggplot(pos.ordered, aes(x=pos.ordered$V3, y=pos.ordered$V4, color=as.factor(table$V2))) + 
    geom_point() + theme_linedraw() +ggtitle(as.character(pro))
}
