###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
###############################################################################
library(limma)
###############################################################################
my.anno<-read.delim("../../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Stage_2=="normal_hsc"| (Stage_2=="diagnosis" & BCR_ABL=="positive")
											 | (Stage_2=="remission" & BCR_ABL=="positive"))
my.anno<-subset(my.anno,Stage_1 !="1_month_TKI")
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
###############################################################################
Top=500
my.important<-read.table("Random_Forests_important_genes.txt",header=T)
my.important<-my.important[1:Top,]

my.genes<-as.character(my.important$gene)
###############################################################################	
S.m2<-subset(S.m.no.batch,rownames(S.m.no.batch) %in% my.genes)
###############################################################################
write.table(S.m2,file="NormalDiagnosisRemission_500_No1M_TKI.txt",append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
###############################################################################
