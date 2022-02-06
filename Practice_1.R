#########################################
##### AUTHOR: Yaiza Arnáiz Alcácer ######
#########################################

#########PRACTICE 1######################

## IMPORT ALL THE LIBRARIES

library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(knitr)

## EXPLORE THE DATA 

data = read.csv("table.csv", row.names=1)
print(data)
print(dim(data))

## EXTRACT THE RELEVANT INFORMATION FOR THE FIRST PCA

samples_1 = data[,1:8]
matrix_01 = as.matrix(samples_1)

## FIRST PCA 

PCA_01 = princomp(matrix_01, cor=FALSE, scores=TRUE)
variance = PCA_01$sdev/sum(PCA_01$sdev)
scores = PCA_01$scores
loadings = PCA_01$loadings


## PLOT PRINCIPAL COMPONENTS

library(ggplot2)
df = data.frame(variance)
df=tibble::rownames_to_column(df, "PCs")
pdf("first_PCA.pdf")
ggplot(df, aes(x=PCs,y=variance)) + geom_col(color='skyblue', fill='skyblue')
+ theme_get() 

dev.off()

## PLOT LOADINGS

df = data.frame(Comp.1=loadings[,1], Comp.2=loadings[,2])
df = tibble::rownames_to_column(df, "lines")
df$position = c(-0.01,0.01,-0.015,0.015,-0.015,0.015,-0.015,0.015)
pdf("first_loadings.pdf")
ggplot(df, aes(x=Comp.1,y=Comp.2,label=lines))+ geom_point(color='skyblue')+
  geom_text(nudge_x=df$position)
dev.off()

## SECOND PCA 

#For every gene in the transcriptome you need to recalculate its expression;
# for instance: gene in transcriptome 25% = 0.25* genei(WT/J0571)+ 0.75*
# gene(shr). 
samples_2 = transform(samples_1,
                      comp25=0.25*samples_1$J0571+0.75*samples_1$shrJ0571,
                      comp50=0.50*samples_1$J0571+0.50*samples_1$shrJ0571,
                      comp75=0.75*samples_1$J0571+0.25*samples_1$shrJ0571)

matrix_02 = as.matrix(samples_2)
PCA_02 = princomp(matrix_02,cor=FALSE,scores=TRUE)
variance = PCA_02$sdev/sum(PCA_02$sdev)
scores = PCA_02$scores
loadings = PCA_02$loadings

## PLOT PRINCIPAL COMPONENTS

df = data.frame(variance)
df=tibble::rownames_to_column(df, "PCs")
ggplot(df, aes(x=PCs,y=variance)) + geom_col(color='magenta', fill='magenta')

## PLOT LOADINGS 

df = data.frame(Comp.1=loadings[,1], Comp.2=loadings[,2])
df = tibble::rownames_to_column(df, "lines")
g = c("shrJ0571", "J0571", "comp25", "comp50", "comp75")
df$genecolor <- rep('skyblue', nrow(df))
df$genecolor[df$lines %in% g]= "magenta"
pdf("second_loadings.pdf")
ggplot(df, aes(x=Comp.1,y=Comp.2,label=lines)) + 
  geom_point(color=df$genecolor) + geom_text(nudge_y=-0.03,color=df$genecolor)
+ theme_get() 
dev.off()

## THIRD PCA
#Add SCR domain 

samples_03 = transform(samples_2, SCR=data$SCRdomain)
matrix_03 = as.matrix(samples_03)
PCA_03 = princomp(matrix_03,cor=FALSE,scores=TRUE)
variance = PCA_03$sdev/sum(PCA_03$sdev)
scores = PCA_03$scores
loadings = PCA_03$loadings

## PLOT PCA

df = data.frame(variance)
df=tibble::rownames_to_column(df, "PCs")
ggplot(df, aes(x=PCs,y=variance)) + geom_col(color='chocolate', fill='chocolate')
+ theme_get() 

## PLOT LOADINGS

g = c("shrJ0571.JKD", "shrJ0571.SCR", "SCR")
df = data.frame(PC1=loadings[,1],PC2=loadings[,2])
df = tibble::rownames_to_column(df, "lines")
df$color <- rep('black', nrow(df))
df$color[df$lines %in% g] = "chocolate"
pdf("third_loadings.pdf")
ggplot(df, aes(x=PC1,y=PC2,label=lines)) + 
  geom_point(color=df$color) + geom_text(nudge_x=0.02,color=df$color)
+ theme_get() 
dev.off()

## PLOT THE MOST IMPORTANT GENES
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(4, "YlOrRd"))(256)
down_genes = labels(head(sort(scores[,"Comp.1"]),20))
up_genes = labels(tail(sort(scores[,"Comp.1"]),20))
up_gene_exp=c()
for (gene in up_genes){up_gene_exp=rbind(up_gene_exp,samples_03[gene,])}
top_20 <- as.matrix(up_gene_exp)
heatmap(top_20, col = colors) 
# the problem is that I can not add a legend with the gradual scale of colours
# with this library, so I tried to do it with ggplot but does 
# not show the dendogram, at the end I used ComplexHeatMap
library(ComplexHeatmap)
pdf("up_genes.pdf")
Heatmap(top_20, 
        name = "Expression", #title of legend
        column_title = "Cell lines", row_title = "Genes",
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        col = colors) 
dev.off()

# we can also plot the less expressed
down_gene_exp=c()
for (gene in down_genes){down_gene_exp=rbind(down_gene_exp,samples_03[gene,])}
bad_20 <- as.matrix(down_gene_exp)
pdf("down_genes.pdf")
Heatmap(bad_20, 
        name = "Expression", #title of legend
        column_title = "Cell lines", row_title = "Genes",
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        col = colors)
dev.off()


