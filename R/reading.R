#' @import hdf5r
#' @import Seurat
#' @export
read_h5 <- function(infile){
  h5 <- Read10X_h5(infile, use.names = TRUE) 
  h5.seurat <- CreateSeuratObject(counts = h5, project ="A")
  h5.seurat <- NormalizeData(h5.seurat)
  h5.seurat <- FindVariableFeatures(h5.seurat, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(h5.seurat)
  h5.seurat <- ScaleData(h5.seurat, features = all.genes)
  h5.seurat <- RunPCA(h5.seurat, features = VariableFeatures(object = h5.seurat))
  ElbowPlot(h5.seurat)
}

#' @export
read_tissue_positions <- function(infile){
  positions <- read.csv(infile, header = F)
  positions <- positions[positions$V2 == 1,]
}

