#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate data for the tSNE analysis of Fig3c.
###############################################################################
library(limma)
###############################################################################
my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Stage_2=="normal_hsc"| Stage_2=="diagnosis")
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
###############################################################################
min_expr<-1
S.m<-as.matrix(my.cells)
S.m[S.m < min_expr]<- 0
S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
#################################With Batch effect#############################
my_bcr_abl<-as.character(my.anno$BCR_ABL)
my_bcr_abl[my_bcr_abl=="negative_low_gapdh"]<-"negative"

my_bcr_abl<-as.factor(my_bcr_abl)
my.design<-model.matrix(~my_bcr_abl)
rownames(my.design)<-colnames(S.m)

batch.x<-as.character(my.anno$process_date)
batch.x[which(batch.x=="batch1" | batch.x=="batch2")]<-"A"
batch.x[which(batch.x=="batch5")]<-"B"
batch.x[which(batch.x=="batch6")]<-"C"

S.m.no.batch = removeBatchEffect(S.m,batch=batch.x,design=my.design)
###############################################################################
load(file="Diagnosis.Pos.HSC.highVariable.Genes.rdata")
load(file="Diagnosis.Neg.HSC.highVariable.Genes.rdata")
load(file="Normal_HSCs.highVariable.Genes.rdata")
###############################################################################
my.genes<-c(Normal_HSCs.highVariable.Genes,Diagnosis.Pos.HSC.highVariable.Genes,Diagnosis.Neg.HSC.highVariable.Genes)
my.genes<-unique(my.genes)
###############################################################################	
S.m2<-subset(S.m.no.batch,rownames(S.m.no.batch) %in% my.genes)
###############################################################################
write.table(S.m2,file="Normal_HSCs_and_Diagnosis_Cells.for.TSNE.txt",append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
###############################################################################
my.anno$BCR_ABL_Color<-1
my.anno$BCR_ABL_Color[my.anno$BCR_ABL=="positive"]<-2
my.anno$BCR_ABL_Color[my.anno$BCR_ABL=="negative"]<-3
my.anno$BCR_ABL_Color[my.anno$BCR_ABL=="negative_low_gapdh"]<-4
##########################Assign Color BCR-ABL#################################
Color1<-c()
for(i in 1:ncol(S.m)){
	
	my.cells<-colnames(S.m)[i]
	
	my.dat<-subset(my.anno,Cell==my.cells)
	
	my.color<-my.dat$BCR_ABL_Color
	
	Color1<-append(Color1,my.color)
	
}
write.table(Color1,file="bcr_abl_status_Normal_HSCs_and_Diagnosis_Cells.for.tSNE.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=FALSE)