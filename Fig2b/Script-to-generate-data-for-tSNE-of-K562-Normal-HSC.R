#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate data for the tSNE analysis.
##############################################################
library(matrixStats)
library(limma)

my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,(Patient=="K562" & BCR_ABL=="primers") | Stage_2=="normal_hsc")
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
###############################################################################
###############################################################################
min_expr<-1
S.m<-as.matrix(my.cells)
S.m[S.m < min_expr]<- 0
S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
###############################################################################
###############################################################################
my_bcr_abl<-as.character(my.anno$Stage_2)
my.design <- model.matrix(~my_bcr_abl)
rownames(my.design)<-colnames(S.m)

batch.x<-as.character(my.anno$process_date)
batch.x[which(batch.x=="batch1" | batch.x=="batch2")]<-"A"##There is no batch effect between batch1 and batch2.
batch.x[which(batch.x=="batch3")]<-"B"
batch.x[which(batch.x=="batch5")]<-"C"
batch.x[which(batch.x=="batch6")]<-"D"

S.m.no.batch = removeBatchEffect(S.m,batch=batch.x,design=my.design)
###############################################################################
load(file="K562.highVariable.Genes.rdata")## read in variable genes found in K562 cells
load(file="Normal_HSCs.highVariable.Genes.rdata")##read in variable genes found in Normal HSCs
my.genes<-c(K562.highVariable.Genes,Normal_HSCs.highVariable.Genes)
my.genes<-unique(my.genes)
###############################################################################	
S.m2<-subset(S.m.no.batch,rownames(S.m.no.batch) %in% my.genes)
###############################################################################
write.table(S.m2,file="NORMAL-HSCs-and-K562.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
###############################################################################