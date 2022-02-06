#########################################
##### AUTHOR: Yaiza Arnáiz Alcácer ######
#########################################

#########PRACTICE 3######################

## IMPORT ALL THE LIBRARIES

library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(knitr)
library(Seurat)

## CREATE SEURAT OBJECT 
# create two to have one with "default parameters" (I will put the minimum number to 
# be less restrictive as possible and have a better contrast) and another with 
# given parameters

# it is need it to read the date as a sparse matrix. 
data = read.table("SCOplanaria.txt", row.names=1)
matrix = as.matrix(data)
sparse_matrix = Matrix(data=matrix, sparse=TRUE)
# using sparse matrixes result in significant memory and speed savings when many cells are zeros
# that is why they are preferred by Seurat
sco_default = CreateSeuratObject(counts=sparse_matrix, project="SCO", min.cells=1, min.features=1)
sco = CreateSeuratObject(counts=sparse_matrix, project="SCO", min.cells=3, min.features=200)
# keep genes expressed in at least 3 cells and cells with at least 200 features 

##CONTROL, NORMALIZATION, FEATURE SELECTION, SCALING

mito_genes <- grep(pattern = "^MT-", x = rownames(x = sco[["RNA"]]), 
                   value = TRUE)
length(mito_genes)
# I did not continue with the code of the lecture because there is no
# mitocondrial genes
#Further exploration
sco_gene = dim(sco)[1]
sco_cell = dim(sco)[2]
sco_gene_per_cell=mean(sco@meta.data$nFeature_RNA)
sco_read_per_cell=mean(sco@meta.data$nCount_RNA)

qc = data.frame(
  sco = c(sco_gene,sco_cell,sco_gene_per_cell,sco_read_per_cell), 
  row.names = c("Number of genes", "Number of cells", 
                "Mean of genes per cell", "Mean of reads per cell"))
pdf("data_sco.pdf")
grid.table(qc) 
dev.off()

mito_genes <- grep(pattern = "^MT-", x = rownames(x = sco_default[["RNA"]]), 
                   value = TRUE)
length(mito_genes)
# I did not continue with the code of the lecture because there is no
# mitocondrial genes
#Further exploration
sco_gene = dim(sco_default)[1]
sco_cell = dim(sco_default)[2]
sco_gene_per_cell=mean(sco_default@meta.data$nFeature_RNA)
sco_read_per_cell=mean(sco_default@meta.data$nCount_RNA)
qc = data.frame(
  sco_default = c(sco_gene,sco_cell,sco_gene_per_cell,sco_read_per_cell), 
  row.names = c("Number of genes", "Number of cells", 
                "Mean of genes per cell", "Mean of reads per cell"))
pdf("data_sco_dafault.pdf")
grid.table(qc) 
dev.off() 

## VIOLIN PLOTS

pdf("violin_plot_comparison.pdf")
v1 = VlnPlot(object = sco, features= c("nFeature_RNA", "nCount_RNA"), ncol = 3) 
v2 = VlnPlot(object = sco_default, features= c("nFeature_RNA", "nCount_RNA"), ncol = 3) 
grid.arrange(v1, v2, ncol=2)
dev.off()


## FEATURE SCATTER PLOT 

pdf("featurescatter_comparison.pdf")
f1 = FeatureScatter(sco, feature1="nFeature_RNA", feature2="nCount_RNA") 
f2 = FeatureScatter(sco_default, feature1="nFeature_RNA", feature2="nCount_RNA") 
grid.arrange(f1, f2, ncol=2)
dev.off()


##QUALITY CONTROL 

sco_default = subset(sco_default, subset=nFeature_RNA>200&nFeature_RNA<2500)
sco = subset(sco, subset=nFeature_RNA>200&nFeature_RNA<2500)
#Just keep the cells with 200 - 2500 number of genes. 


# LogNormalize with scale factor of 1000
sco = NormalizeData(sco, normalization.method="LogNormalize", scale.factor=10000)
sco_default = NormalizeData(sco_default, normalization.method="LogNormalize", scale.factor=10000)
sco = AddMetaData(sco, colSums(sco[["RNA"]]@data), col.name="Normalised_RNA")
sco_default = AddMetaData(sco_default, colSums(sco[["RNA"]]@data), col.name="Normalised_RNA")

## VIOLIN PLOTS AFTER NORMALIZATION

pdf("violin_plot_norm_comparison.pdf")
vp1 = VlnPlot(sco, features="Normalised_RNA") 
vp2 =  VlnPlot(sco_default, features="Normalised_RNA") 
grid.arrange(vp1, vp2, ncol=2)
dev.off()

sco_default = ScaleData(sco_default, features=rownames(sco_default))
sco = ScaleData(sco, features=rownames(sco))

## FIND VARIABLE FEATURES 

sco_default = FindVariableFeatures(sco_default, selection.method = "vst", nfeatures=3000)
sco = FindVariableFeatures(sco, selection.method = "vst", nfeatures=3000)

pdf("variable_features_comparison.pdf")
plot1 = LabelPoints(plot=VariableFeaturePlot(sco), 
                   points=head(VariableFeatures(sco),5), repel=TRUE)

plot2 = LabelPoints(plot=VariableFeaturePlot(sco_default), 
                    points=head(VariableFeatures(sco_default),5), repel=TRUE)

grid.arrange(vp1, vp2, ncol=2)
dev.off()

## PCA 

sco = RunPCA(sco, features=VariableFeatures(object=sco))

sco_default = RunPCA(sco_default, features=VariableFeatures(object=sco_default))


pdf("Dimheatmap.pdf")
DimHeatmap(sco, dims=1:5, cells=500, balanced = TRUE) # we should use the first 5 principal components
dev.off()


VizDimLoadings(sco, dims = 1:5, reduction = "pca")


DimPlot(sco, reduction = "pca")


#Before do the t-sne clustering, since the Seurat cluster cells are based on their 
#PCA scores, we should think in how many prinicipal components include. 
#To know how much we can implement a JackStraw procedure:

sco <- JackStraw(sco, num.replicate = 100)
sco <- ScoreJackStraw(sco, dims = 1:20)
JackStrawPlot(sco, dims = 1:20)

### CLUSTERING OF THE FIRST 5 PRINCIPAL COMPONENTS

sco = FindNeighbors(sco, dims=1:5)
sco = FindClusters(sco, resolution=0.6)

sco_default = FindNeighbors(sco_default, dims=1:5)
sco_default = FindClusters(sco_default, resolution=0.6)


#Finally to visualize the data we can apply UMAP and TSNE.

sco = RunUMAP(sco, dims=1:5)
sco = RunTSNE(sco, dims=1:5)

pdf("umap_tsne_sco.pdf")
plot1 = DimPlot(sco, reduction="umap", label=TRUE, pt.size=3)
plot2 = DimPlot(sco, reduction="tsne", label=TRUE, pt.size=3) 
grid.arrange(plot1, plot2, nrow=2)
dev.off()


sco_default = RunUMAP(sco_default, dims=1:5)
sco_default = RunTSNE(sco_default, dims=1:5)

pdf("umap_tsne_sco_default.pdf")
plot1 = DimPlot(sco_default, reduction="umap", label=TRUE, pt.size=3)
plot2 = DimPlot(sco_default, reduction="tsne", label=TRUE, pt.size=3) 
grid.arrange(plot1, plot2, nrow=2)
dev.off()


## FIND MARKERS

markers = FindAllMarkers(sco, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)

head(markers, n=5)


top5 = markers %>% group_by(cluster) %>% top_n (n=5, wt= avg_log2FC) 

top5_table = top5 %>% group_by(cluster) %>% 
  summarise(gene=paste(gene,collapse=', '))


### HEATMAP

pdf("heatmap_top5.pdf")
DoHeatmap(sco, features=top5$gene) 
dev.off()

## VIOLIN PLOTS 

pdf("Top5_violin_plots.pdf")
v1 <- VlnPlot(sco, features = c("dd-Smed-v6-61-0","dd-Smed-v6-2178-0"))
v2 <- VlnPlot(sco, features = c("dd-Smed-v6-298-0","dd-Smed-v6-1410-0"))
v3 <- VlnPlot(sco, features = c("dd-Smed-v6-702-0","dd-Smed-v6-2548-0"))
v4 <- VlnPlot(sco, features = c("dd-Smed-v6-9977-0","dd-Smed-v6-48-0"))
v5 <- VlnPlot(sco, features = c("dd-Smed-v6-175-0","dd-Smed-v6-1161-1"))

## FEATURES PLOTS
pdf("Top5_featureplots.pdf")
f1 <- FeaturePlot(sco, features = c("dd-Smed-v6-61-0","dd-Smed-v6-2178-0"))
f2 <- FeaturePlot(sco, features = c("dd-Smed-v6-298-0","dd-Smed-v6-1410-0"))
f3 <- FeaturePlot(sco, features = c("dd-Smed-v6-702-0","dd-Smed-v6-2548-0"))
f4 <- FeaturePlot(sco, features = c("dd-Smed-v6-9977-0","dd-Smed-v6-48-0"))
f5 <- FeaturePlot(sco, features = c("dd-Smed-v6-175-0","dd-Smed-v6-1161-1"))
grid.arrange(f1, f2, f3, f4, f5, ncol=2, nrow=2)
dev.off()

## TSNE with the biomarkers of the table


identities <- c("Neural progenitors", "Early epidermal progenitors", 
                "Muscle body", "Epidermis", "Early epidermal progenitors", "Early epidermal 
progenitors", "Phagocytes", "GABA neurons", "Epidermis","Pigment","Late 
epidermal progenitors","Parenchymal cells", "Late epidermal progenitors")
names(identities) <- levels(sco)
sco <- RenameIdents(sco, identities)
pdf("tsne_table.pdf")
DimPlot(sco, reduction = "tsne", label = TRUE, pt.size = 3) + 
  NoLegend()
dev.off()


## VIOLIN PLOT NEURAL PROGENITORS

pdf("violin_neural.pdf")
VlnPlot(sco, features="dd-Smed-v6-1999-0") 
dev.off()

## FEATURE PLOT

pdf("feature_plot.pdf")
FeaturePlot(sco, features="dd-Smed-v6-1999-0") + NoLegend()
dev.off()


## TSNE WITH NEURAL PROGENITORS


new_identities <- c("Neural progenitors", "Early epidermal progenitors", 
                    "Muscle body", "Neoblast (stem)", "Early epidermal progenitors", "Neoblast 
(stem)", "Neoblast (stem)", "GABA neurons", "Epidermis","Pigment","Late 
epidermal progenitors","Parenchymal cells", " Neoblast (stem)")
names(new_identities) = levels(sco)
sco <- RenameIdents(sco, new_identities)
pdf("last_tsne.pdf")
DimPlot(sco, reduction = "tsne", label = TRUE, pt.size = 3) + NoLegend()
dev.off()







