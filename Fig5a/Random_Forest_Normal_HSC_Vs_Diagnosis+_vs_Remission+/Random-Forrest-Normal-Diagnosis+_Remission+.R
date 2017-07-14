###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used for the Random forests analysis.
###############################################################################
library(matrixStats)
library(cluster)
library(RColorBrewer)
library(pheatmap)
library(Rtsne)
library(limma)
library(randomForest)

my.anno<-read.delim("../../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Stage_2=="normal_hsc"| (Stage_2=="diagnosis" & BCR_ABL=="positive")
											 | (Stage_2=="remission" & BCR_ABL=="positive"))
###############################################################################
load(file="../../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
###############################################################################
min_expr<-1
S.m<-as.matrix(my.cells)
S.m[S.m < min_expr]<- 0
S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
#################################With Batch effect#############################
my_bcr_abl<-as.character(my.anno$Stage_2)

my_bcr_abl<-as.factor(my_bcr_abl)
my.design<-model.matrix(~my_bcr_abl)
rownames(my.design)<-colnames(S.m)

batch.x<-as.character(my.anno$process_date)
batch.x[which(batch.x=="batch1" | batch.x=="batch2")]<-"A"
batch.x[which(batch.x=="batch3")]<-"B"
batch.x[which(batch.x=="batch5")]<-"C"
batch.x[which(batch.x=="batch6")]<-"D"

S.m.no.batch = removeBatchEffect(S.m,batch=batch.x,design=my.design)
###########################################################################################
###########################################################################################
########The following genes are the top important genes for each round of the RF analysis
########in 3 ways comparison:
########1) Normal Vs Diagnosis+
########2) Diagnosis+ Vs Remission+
########3) Normal Vs Remission+
########Genes are selected based on the cut off criteria :"mean(my.importance.f$MeanDecreaseGini)+sd(my.importance.f$MeanDecreaseGini)"
###########################################################################################
load(file="RF.genes.rdata")
load(file="RF.genes.diagnosis.remission.rdata")
load(file="RF.genes.normal.remission.rdata")
###########All Genes#######################################################################
my.all.genes<-c(my.genes,my.genes.diagnosis.remission,my.genes.normal.remission)
my.genes<-unique(my.all.genes)
###############################
S.m.no.batch<-S.m.no.batch[rownames(S.m.no.batch) %in% my.genes,]
#########Run Random forrest################################################################
group<-my_bcr_abl

rf500_normal_diagnosis_remission <- randomForest(x=t(S.m.no.batch),y=as.factor(group),ntree=500,importance=TRUE) 
rf1000_normal_diagnosis_remission <- randomForest(x=t(S.m.no.batch),y=as.factor(group),ntree=1000,importance=TRUE) 
rf2000_normal_diagnosis_remission <- randomForest(x=t(S.m.no.batch),y=as.factor(group),ntree=2000,importance=TRUE)

###############################################################################
save(rf500_normal_diagnosis_remission,file="RandomForrest_Normal_Diagnosis_Remission_500tree.rdata")
save(rf1000_normal_diagnosis_remission,file="RandomForrest_Normal_Diagnosis_Remission_1000tree.rdata")
save(rf2000_normal_diagnosis_remission,file="RandomForrest_Normal_Diagnosis_Remission_2000tree.rdata")
###############################################################################