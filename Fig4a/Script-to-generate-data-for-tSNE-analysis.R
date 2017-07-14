###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate data for the tSNE analysis.#################
library(limma)
###############################################################################
my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Stage_2=="diagnosis" & BCR_ABL=="positive")
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
###############################################################################
min_expr<-1
S.m<-as.matrix(my.cells)
S.m[S.m < min_expr]<- 0
S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
#################################remove batches effect#########################
my_responder<-as.character(my.anno$Responder)
my_responder[is.na(my_responder)==TRUE]<-"unknown"

my_responder<-as.factor(my_responder)
my.design<-model.matrix(~my_responder)
rownames(my.design)<-colnames(S.m)

batch.x<-as.character(my.anno$process_date)
batch.x[which(batch.x=="batch1" | batch.x=="batch2")]<-"A"
batch.x[which(batch.x=="batch5")]<-"B"
batch.x[which(batch.x=="batch6")]<-"C"

S.m.no.batch = removeBatchEffect(S.m,batch=batch.x,design=my.design)
##############################################################################
load(file="Diagnosis.Pos.HSC.highVariable.Genes.rdata")
my.genes<-Diagnosis.Pos.HSC.highVariable.Genes
###############################################################################	
S.m2<-subset(S.m.no.batch,rownames(S.m.no.batch) %in% my.genes)
###############################################################################
write.table(S.m2,file="Diagnosis_BCR-ABL+_5011_genes.for.TSNE.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
###############################################################################