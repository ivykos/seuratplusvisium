#' @import hdf5r
#' @import Seurat
#' @import ggplot2
#' @import RColorBrewer
#' @import viridis



#' Reads an h5 file from Space Ranger output. Generated a Seurat object for data normalization, feature selection, scaling, and PCA.
#' 
#' @param h5 h5 file
#' @param proj Project name for Seurat object
#' @param a Assay
#' @example SeuratObj <- read_h5("example.h5","my_data","RNA")
#' @export
read_h5 <- function(h5, proj, a){
  h5_file <- Read10X_h5(h5, use.names = TRUE) 
  h5.seurat <- CreateSeuratObject(counts = h5_file, project = proj, assay = a)
  h5.seurat <- NormalizeData(h5.seurat)
  h5.seurat <- FindVariableFeatures(h5.seurat, selection.method = "vst", nfeatures = 2000)
  h5.seurat <- ScaleData(h5.seurat)
  h5.seurat <- RunPCA(h5.seurat, features = VariableFeatures(object = h5.seurat))
  
}

#' After reading an h5 file and performing pre-processing, this function will perform UMAP dimensional reduction using Seurat. 
#' Output from this function will also include a visualization of the Visium slide with each dot colored to correspond to the 
#' assigned cluster.
#'
#' @param sueratObj Seurat object
#' @param proj Project name for Seurat object
#' @param dimensions Dimensionality of the dataset. Use Seurat's JackStraw or ElbowPlot function to determine this
#' @param res Resolution
#' @param tissue_csv tissue_positions_list.csv file from Space Ranger output
#' @example transfer_clusters(SeuratObj,"my_data",13,0.5,"tissue_positions_list.csv")
#' @export
transfer_clusters <- function(seuratObj, proj, dimensions, res, tissue_csv){
  seuratObj <- FindNeighbors(seuratObj, dims = 1:dimensions)
  seuratObj <- FindClusters(seuratObj, resolution = res)
  
  seuratObj <- RunUMAP(seuratObj, dims = 1:dimensions)
  
  
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
  
  
  ggplot(pos.ordered, aes(x=pos.ordered$V3, y=pos.ordered$V4, color=as.factor(table$V2))) + 
    geom_point() + theme_linedraw() +ggtitle(as.character(proj))
}

#' @export
get_expression <- function(seurat_obj, feature){
  tmp <- as.matrix(GetAssayData(object = seurat_obj, slot = "counts"))
  feature_count <- tmp[feature,]
  rotated <- as.data.frame(t(feature_count))
  table.c["Expr"] <- rotated$feature_count
  
}

