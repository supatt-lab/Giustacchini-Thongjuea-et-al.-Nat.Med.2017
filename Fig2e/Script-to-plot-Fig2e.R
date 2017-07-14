#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig2e.
###############################################################################
library(matrixStats)
library(cluster)
library(RColorBrewer)
library(pheatmap)
library("gplots")

my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Patient=="OX1407")
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
###############################################################################
min_expr<-1
S.m<-as.matrix(my.cells)
S.m[S.m < min_expr]<- 0
S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
#############################################################################
ox1407.sig<-read.table("DE-Genes-BCR-ABL+_Vs_BCR-ABL-_OX1407.txt",header=T)
ox1407.sig<-ox1407.sig[1:200,]
ox1407.sig<-ox1407.sig[order(abs(ox1407.sig$log2fc),decreasing=T),]
ox1407.sig<-ox1407.sig[1:75,]
#############################################################################
my.genes<-as.character(ox1407.sig$gene)

final.genes.m<-subset(S.m,rownames(S.m) %in% my.genes)
#############################################################################
my.anno$BCR_ABL_Color<-"brown"
my.anno$BCR_ABL_Color[my.anno$BCR_ABL=="negative"]<-"blue"
##########################Assign Color BCR-ABL###############################
Color1<-c()
for(i in 1:ncol(final.genes.m)){
	my.cells<-colnames(final.genes.m)[i]
	my.dat<-subset(my.anno,Cell==my.cells)
	my.color<-my.dat$BCR_ABL_Color
	Color1<-append(Color1,my.color)
}
h1.m<-final.genes.m
################save data matrix to a file####################################
write.table(h1.m,file="Fig2e.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)
##########################Plot heatmap########################################
hr <- hclust(as.dist(1-cor(t(h1.m), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(h1.m, method="pearson")), method="complete") 


hr.dendrogram<- as.dendrogram(hr)
hc.dendrogram<- as.dendrogram(hc)

Rowv        = rowMeans(h1.m, na.rm = T)
hr.dendrogram  = reorder(hr.dendrogram, Rowv)

Colr        = colMeans(h1.m, na.rm = T)
hc.dendrogram  = reorder(hc.dendrogram, Colr)

pdf(file="Fig2e.pdf", width=9, height=8.5, onefile=T, bg="transparent")

heatmap.2(h1.m,Rowv=hr.dendrogram,Colv=hc.dendrogram, col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
		scale="none",key=TRUE,symkey=FALSE,keysize = 0.8, density.info="none", trace="none",cexRow=0.75,cexCol=0.8,ColSideColors=Color1,
		margin=c(8, 10))
dev.off()
##############################################################################




