###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig6b.
###############################################################################
library(RColorBrewer)
library(pheatmap)
library(gplots)
###############################################################################
my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,(Stage_2=="diagnosis" | Stage_2=="blast_crisis" | 
					     Stage_2=="blast_crisis_1_month_TKI") & Patient !="OX1931cd19" & Patient !="OX2046cd19")
my.anno<-subset(my.anno,BCR_ABL =="positive")
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
remove.index<-grep("MIR|SNORA|SNORD",rownames(my.cells))
my.cells<-my.cells[-c(remove.index),]
###############################################################################
min_expr<-1
S.m<-as.matrix(my.cells)
S.m[S.m < min_expr]<- 0
S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
###############################################################################
groupA<-subset(my.anno,Stage_2=="diagnosis")
groupA.cells<-as.character(groupA$Cell)
groupA.m<-S.m[,colnames(S.m) %in% groupA.cells]

groupB<-subset(my.anno,Stage_2=="blast_crisis")
groupB.cells<-as.character(groupB$Cell)

groupB.m<-S.m[,colnames(S.m) %in% groupB.cells]
###############################################################################
X.m<-cbind(groupA.m,groupB.m)
###############################################################################
dat1.sig<-read.table("CP-CML+_Vs_BC-CML+_without_PreBC.Siggenes.txt",header=T)
dat1.sig<-subset(dat1.sig,abs(log2fc) >=2)
dat1.sig<-dat1.sig[1:40,]
#################
my.genes<-as.character(dat1.sig$gene)
###############################################################################
###############################################################################
X.m<-subset(X.m,rownames(X.m) %in% my.genes)
###############################################################################
################save data matrix to a file#####################################
write.table(X.m,file="Fig6b.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)
###############################################################################
my.anno.x<-as.character(my.anno$Stage_2)

annotation_col<-data.frame(Cell=my.anno.x)
rownames(annotation_col)<-my.anno$Cell
###############################################################################
ann_colors = list(
		Cell = c(blast_crisis = "purple",diagnosis="brown")
)
##################Heatmap######################################################
pdf(file="Fig6b.pdf", width=10, height=8, onefile=F, bg="transparent")

pheatmap(X.m,annotation_col=annotation_col,annotation_colors=ann_colors
		,show_colnames=F,show_rownames=T, fontsize_row=10,key=TRUE,symkey=FALSE,keysize = 0.8,
		col=colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),cexRow=0.75,cexCol=0.75)

dev.off()