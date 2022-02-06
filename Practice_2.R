###########################################
##### AUTHOR: YAIZA ARNAIZ ALCACER ########
###########################################
## IMPORT ALL THE LIBRARIES

library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(knitr)

## READ DATA 
data=read.csv("TableS5.csv", row.names=1)
matrix_01 = as.matrix(data)
print(matrix_01)

## PLOT THE DISTRIBUTIONS
library(tidyr)
df=data.frame(matrix_01)
df=tibble::rownames_to_column(df, "Genes")
df=gather(df, names(df[-1]), key="Types", value="Expression")
pdf("ditributions_types.pdf")
ggplot(df, aes(x=Types,y=Expression)) + geom_boxplot(color="blue")
dev.off()

## PLOT PCA 
PCA01 = princomp(matrix_01,cor=FALSE,scores=TRUE)
variance = PCA01$sdev[1:9]/sum(PCA01$sdev[1:9])
scores = PCA01$scores
loadings = PCA01$loadings
df = data.frame(variance)
df = tibble::rownames_to_column(df, "PCs")
pdf("PCA_biomarkers.pdf")
ggplot(df, aes(x=PCs,y=variance)) + geom_col(color='blue', fill='blue')
dev.off()

## PLOT LOADINGS

df = data.frame(Comp.1=loadings[,1], Comp.2=loadings[,2], Comp.3=loadings[,3], Comp.4=loadings[,4])
df = tibble::rownames_to_column(df, "types")
df$color = rep("black", nrow(df))
df$color[df$types %in% c("E30", "S18")] = "blue"
pdf("loadings-biomarkers.pdf")
plot1 = ggplot(df, aes(x=Comp.1,y=Comp.2,label=types)) + 
  geom_point(color=df$color) + geom_text(nudge_y=0.1,color=df$color)
plot2 = ggplot(df, aes(x=Comp.3,y=Comp.4,label=types)) + 
  geom_point(color=df$color) + geom_text(nudge_y=0.1,color=df$color)
grid.arrange(plot1, plot2, ncol=2)
dev.off()

## HEATMAP

library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(4, "YlOrRd"))(256)
stem_cells = subset(data, data$WOX5>1 & rowSums(data[,-which(names(data)=="WOX5")]>1)<1)
matrix_stem= data.matrix(stem_cells)
library(ComplexHeatmap)
pdf("heatmap_biomarkers.pdf")
Heatmap(matrix_stem, 
        name = "Expression", #title of legend
        column_title = "Cell types", row_title = "Genes", # Text size for row names
        col = colors, show_row_names = FALSE) 
dev.off()

## PLOT PCA 2, ONLY LOADINGS

PCA02 = princomp(matrix_stem,cor=FALSE,scores=TRUE)
loadings = PCA02$loadings
df = data.frame(Comp.1=loadings[,1], Comp.2=loadings[,2], Comp.3=loadings[,3], Comp.4=loadings[,4])
df = tibble::rownames_to_column(df, "Types")
df$color = rep('black', nrow(df))
df$color[df$Types %in% c("WOX5")] = "blue"
df$Types[!df$Types %in% c("WOX5")] = " "
pdf("second_pca_biomarkers.pdf")
plot1=ggplot(df, aes(x=Comp.1,y=Comp.2,label=Types)) + 
  geom_point(color=df$color) + geom_text(nudge_y=0.1,color=df$color)
plot2=ggplot(df, aes(x=Comp.3,y=Comp.4,label=Types)) + 
  geom_point(color=df$color) + geom_text(nudge_y=0.1,color=df$color)
grid.arrange(plot1, plot2, ncol=2)
dev.off()


## PLOT SCORES
df = data.frame(Comp.1=scores[,1], Comp.2=scores[,2], Comp.3=scores[,3], Comp.4=scores[,4])
n = rownames(matrix_stem)
df$color = rep("black", nrow(df))
df$color[rownames(df) %in% n]= "violet"
df$size = rep(1, nrow(df))
df$size[rownames(df) %in% n]=5
pdf("scores.pdf")
plot1 = ggplot(df, aes(x=Comp.1,y=Comp.2)) + geom_point(color=df$color, size=df$size)
plot2 = ggplot(df, aes(x=Comp.3,y=Comp.4)) + geom_point(color=df$color, size=df$size)
grid.arrange(plot1, plot2, ncol=2)
dev.off()





