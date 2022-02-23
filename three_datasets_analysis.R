library(Seurat)
library(hdf5r)
library(scales)
library(RColorBrewer)
library(patchwork)
library(ggplot2)
library(paletteer)
library(viridis)


CAF.A <- Read10X_h5("CAF_A.h5", use.names = TRUE)
CAF.B <- Read10X_h5("CAF_B.h5", use.names = TRUE)
CAF.C <- Read10X_h5("CAF_C.h5", use.names = TRUE)

caf_a <- CreateSeuratObject(CAF.A, project="APOE2", assay="RNA")
caf_b <- CreateSeuratObject(CAF.B, project="APOE3", assay="RNA")
caf_c <- CreateSeuratObject(CAF.C, project="APOE4", assay="RNA")

combined1 <- merge(caf_a, y = caf_b, project = "combined")
combined2 <- merge(combined1, y = caf_c, project = "combined2") #The one

#Clustering of the big combined dataset
combined2 <- NormalizeData(combined2)
combined2 <- FindVariableFeatures(combined2, selection.method = "vst", nfeatures = 2000)

combined2 <- ScaleData(combined2)
combined2 <- RunPCA(combined2, features = VariableFeatures(object = combined2))
combined2 <- FindNeighbors(combined2, dims = 1:13)
combined2 <- FindClusters(combined2, resolution = 0.5)

combined2 <- RunUMAP(combined2, dims = 1:13)
DimPlot(combined2, reduction = "umap")

#data.markers <- FindAllMarkers(combined2, test.use = "t")
caf_a <- NormalizeData(caf_a)
caf_a <- FindVariableFeatures(caf_a, selection.method = "vst", nfeatures = 2000)

caf_a <- ScaleData(caf_a)
caf_a <- RunPCA(caf_a, features = VariableFeatures(object = caf_a))
caf_a <- FindNeighbors(caf_a, dims = 1:13)
caf_a <- FindClusters(caf_a, resolution = 0.5)

caf_a <- RunUMAP(caf_a, dims = 1:13)


caf_b <- NormalizeData(caf_b)
caf_b <- FindVariableFeatures(caf_b, selection.method = "vst", nfeatures = 2000)

caf_b <- ScaleData(caf_b)
caf_b <- RunPCA(caf_b, features = VariableFeatures(object = caf_b))
caf_b <- FindNeighbors(caf_b, dims = 1:13)
caf_b <- FindClusters(caf_b, resolution = 0.5)

caf_b <- RunUMAP(caf_b, dims = 1:13)

caf_c <- NormalizeData(caf_c)
caf_c <- FindVariableFeatures(caf_c, selection.method = "vst", nfeatures = 2000)

caf_c <- ScaleData(caf_c)
caf_c <- RunPCA(caf_c, features = VariableFeatures(object = caf_c))
caf_c <- FindNeighbors(caf_c, dims = 1:14)
caf_c <- FindClusters(caf_c, resolution = 0.5)

caf_c <- RunUMAP(caf_c, dims = 1:14)

identities <- levels(combined2@active.ident)
my_color_palette <- hue_pal()(length(identities))
palette(my_color_palette)

#Spatial Data
positions.caf.a <- read.csv("CAF_A_positions.csv", header = F)
positions.caf.b <- read.csv("CAF_B_positions.csv", header = F)
positions.caf.c <- read.csv("CAF_C_positions.csv", header = F)

#Only spots inside the tissue
positions.caf.a <- positions.caf.a[positions.caf.a$V2 == 1,]
positions.caf.b <- positions.caf.b[positions.caf.b$V2 == 1,]
positions.caf.c <- positions.caf.c[positions.caf.c$V2 == 1,]


#CAF_A
plot(positions.caf.a$V3, positions.caf.a$V4)

write.table(combined2@active.ident, file="Cell_Clusters.tsv", quote=FALSE, sep="\t", col.names = FALSE)
idents <- read.delim("Cell_Clusters.tsv", header = FALSE)

write.table(combined2@meta.data$orig.ident, file="Cell_origin.tsv", quote=FALSE, sep="\t", col.names = FALSE)
orig <- read.delim("Cell_origin.tsv", header = FALSE)

idents["Origin"] = orig$V2
table.a <- idents[idents$Origin == "CAF_A",]
cluster.ordered <- table.a[order(table.a$V1),]
pos.ordered <-positions.caf.a[order(positions.caf.a$V1),]

plot(pos.ordered$V3, pos.ordered$V4, col=cluster.ordered$V2)
legend("topright",unique(cluster.ordered$V2),col=1:length(cluster.ordered$V2),pch=16, legend = )

#New plotting method
ggplot(pos.ordered, aes(x=pos.ordered$V3, y=pos.ordered$V4, color=as.factor(table.a$V2))) + geom_point() + theme_linedraw() +ggtitle("APOE2")

#
write.table(as.matrix(GetAssayData(object = caf_a, slot = "counts")), 
            'counts_a.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)



ggplot(pos.ordered, aes(pos.ordered$V3, pos.ordered$V4)) + 
  geom_point(aes(color = table.a$Expr), size = 2) + 
  scale_color_viridis(option = "inferno") + theme_bw()+ ggtitle("APOE2 Gfap Expression")

#CAF_B
plot(positions.caf.b$V3, positions.caf.b$V4)
table.b <- idents[idents$Origin == "CAF_B",]
cluster.ordered <- table.b[order(table.b$V1),]
pos.ordered <-positions.caf.b[order(positions.caf.b$V1),]
plot(pos.ordered$V3, pos.ordered$V4, col=cluster.ordered$V2, pch=16)
legend("topright",unique(cluster.ordered$V2),col=1:length(cluster.ordered$V2),pch=16, legend = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13))
#
write.table(as.matrix(GetAssayData(object = caf_b, slot = "counts")), 
            'counts_b.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

counts_b <- read.csv("counts_b.csv")
mbp_b <- counts_b["Gfap",]
rotated <- as.data.frame(t(mbp_b))
table.b["Expr"] <- rotated$Gfap

ggplot(pos.ordered, aes(pos.ordered$V3, pos.ordered$V4)) + 
  geom_point(aes(color = table.b$Expr), size = 2) +ggtitle("APOE3") + theme_linedraw()

#CAF_C
plot(positions.caf.c$V3, positions.caf.c$V4)
table.c <- idents[idents$Origin == "CAF_C",]
cluster.ordered <- table.c[order(table.c$V1),]
pos.ordered <-positions.caf.c[order(positions.caf.c$V1),]
plot(pos.ordered$V3, pos.ordered$V4, col=cluster.ordered$V2, pch=16)
legend("topright",unique(cluster.ordered$V2),col=1:length(cluster.ordered$V2),pch=16, legend = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13))

#
write.table(as.matrix(GetAssayData(object = caf_c, slot = "counts")), 
            'counts_c.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

counts_c <- read.csv("counts_c.csv")
mbp_c <- counts_c["Gfap",]
rotated <- as.data.frame(t(mbp_c))
table.c["Expr"] <- rotated$Gfap

ggplot(pos.ordered, aes(pos.ordered$V3, pos.ordered$V4)) + 
  geom_point(aes(color = table.c$Expr), size = 2) +ggtitle("APOE4") + scale_fill_gradientn()


# Spatial expression

counts_a <- read.csv("counts_a.csv")
counts_b <- read.csv("counts_b.csv")
counts_c <- read.csv("counts_c.csv")

table.b <- idents[idents$Origin == "CAF_B",]
table.c <- idents[idents$Origin == "CAF_C",]

#CAF_A
pos.ordered <-positions.caf.a[order(positions.caf.a$V1),]
mbp_a <- counts_a["Gad1",]
rotated <- as.data.frame(t(mbp_a))
table.a["Expr"] <- rotated$Gad1
mbp.a <- ggplot(pos.ordered, aes(pos.ordered$V3, pos.ordered$V4)) + 
  geom_point(aes(color = table.a$Expr), size = 2) + 
  scale_color_viridis(option = "inferno") + theme_bw()+ ggtitle("APOE2 Gad1")

#CAF_B
pos.ordered <-positions.caf.b[order(positions.caf.b$V1),]
mbp_b <- counts_b["Gad1",]
rotated <- as.data.frame(t(mbp_b))
table.b["Expr"] <- rotated$Gad1

mbp.b <- ggplot(pos.ordered, aes(pos.ordered$V3, pos.ordered$V4)) + 
  geom_point(aes(color = table.b$Expr), size = 2) + 
  scale_color_viridis(option = "inferno") + theme_bw()+ ggtitle("APOE3 Gad1")

#CAF_C
pos.ordered <-positions.caf.c[order(positions.caf.c$V1),]
mbp_c <- counts_c["Gad1",]
rotated <- as.data.frame(t(mbp_c))
table.c["Expr"] <- rotated$Gad1
mbp.c <- ggplot(pos.ordered, aes(pos.ordered$V3, pos.ordered$V4)) + 
  geom_point(aes(color = table.c$Expr), size = 2) + 
  scale_color_viridis(option = "inferno") + theme_bw()+ ggtitle("APOE4 Gad1")

(mbp.a) / (mbp.b) / (mbp.c)