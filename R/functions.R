#' @import hdf5r
#' @import Seurat
#' @import ggplot2
#' @import RColorBrewer
#' @import viridis
#' @import dplyr

#' Takes h5 file and creates a Seurat object, normalizes, does variable feature selection, scales the data, and performs PCA.
#' 
#' @param h5 h5 file from space ranger outputt.
#' @param proj Project name
#' @param a Assay
#' @return Seurat object
#' @examples
#' read_h5("example.h5", "my_project", "RNA")
#' @export
read_h5 <- function(h5, proj, a){
  h5_file <- Read10X_h5(h5, use.names = TRUE) 
  h5.seurat <- CreateSeuratObject(counts = h5_file, project = proj, assay = a)
  h5.seurat <- NormalizeData(h5.seurat)
  h5.seurat <- FindVariableFeatures(h5.seurat, selection.method = "vst", nfeatures = 2000)
  h5.seurat <- ScaleData(h5.seurat)
  h5.seurat <- RunPCA(h5.seurat, features = VariableFeatures(object = h5.seurat))
  
}

#' Takes Seurat object and does dimensional reduction using UMAP and visualizes the cluster assignments on the Visium tissue sample.
#' 
#' @param seuratObj Seurat object.
#' @param proj Project name
#' @param dimensions Dimensionality of the data set. Use Seurat's ElbowPlot or JackStraw function for determining this. 
#' @param res Resolution
#' @param tissue_csv "tissue_positions.csv" output from space ranger pipeline
#' @return Visualization of cluster labels on tissue sample
#' @examples
#' transfer_clusters(Obj,"my_project",13,0.5,"tissue_positions.csv")
#' @export
transfer_clusters <- function(seuratObj, proj, dimensions, res, tissue_csv){
  seuratObj <- FindNeighbors(seuratObj, dims = 1:dimensions)
  seuratObj <- FindClusters(seuratObj, resolution = res)
  
  seuratObj <- RunUMAP(seuratObj, dims = 1:dimensions)
  
  #Retrieve the cluster labels and project ID from the Seurat object and put them in tables
  write.table(seuratObj@active.ident, file="tmp.tsv", quote=FALSE, sep="\t", col.names = FALSE)
  idents <- read.delim("tmp.tsv", header = FALSE)
  file.remove("tmp.tsv")
  write.table(seuratObj@meta.data$orig.ident, file="tmp.tsv", quote=FALSE, sep="\t", col.names = FALSE)
  orig <- read.delim("tmp.tsv", header = FALSE)
  file.remove("tmp.tsv")
  
  #Get the coordinates for plotting the tissue image
  positions <- read.csv(tissue_csv, header = F)
  
  #Make sure only the spots within the tissue are included for this analysis
  positions <- positions[positions$V2 == 1,]
  
  #Put everything in one table for simplicity
  idents["Origin"] = orig$V2
  
  #"table" is the table we will use from now on
  table <- idents[idents$Origin == proj,]
  
  cluster.ordered <- table[order(table$V1),]
  pos.ordered <-positions[order(positions$V1),]
  
  #Plot it
  ggplot(pos.ordered, aes(x=pos.ordered$V3, y=pos.ordered$V4, color=as.factor(table$V2))) + 
    geom_point() + theme_linedraw() +ggtitle(as.character(proj))
}
#' Visualizes spatial expression of a gene
#' 
#' @param seuratObj Seurat object.
#' @param proj Project name
#' @param feature Gene or feature you want to visualize
#' @param tissue_csv "tissue_positions.csv" output from space ranger pipeline
#' @return Visualization of gene expression on tissue sample
#' @examples
#' get_expression(Obj,"my_project","Ripk1","tissue_positions.csv")
#' @export
get_expression <- function(seuratObj, proj, feature, tissue_csv){
  
  #Get tissue positions again 
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
  
  #Extract the counts matrix from the Seurat object
  tmp <- as.matrix(GetAssayData(object = seuratObj, slot = "counts"))
  feature_count <- tmp[feature,]
  rotated <- as.data.frame(t(feature_count))
  rotated <- as.data.frame(t(feature_count))
  table["Expr"] <- rotated$V1
  
  ggplot(pos.ordered, aes(pos.ordered$V3, pos.ordered$V4)) + 
    geom_point(aes(color = table$Expr), size = 2) + 
    scale_color_viridis(option = "inferno") + theme_bw()
}

