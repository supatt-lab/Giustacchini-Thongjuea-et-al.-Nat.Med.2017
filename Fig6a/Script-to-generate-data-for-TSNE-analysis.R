###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate data for TSNE analysis of Fig6a.
###############################################################################
library(limma)
library(Rtsne)

my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T,stringsAsFactors=F)
my.anno<-subset(my.anno,used4analysis==1)

my.anno<-subset(my.anno,(Stage_2=="normal_hsc" | (Stage_2=="diagnosis" & BCR_ABL=="positive") | BCR_ABL=="primers" |
					Stage_2=="pre_blast_crisis" | Stage_2=="blast_crisis" |
					Patient =="CML1203") & Patient !="OX1931cd19" & Patient !="OX2046cd19")

my.anno<-subset(my.anno,BCR_ABL =="positive" | Stage_2=="normal_hsc" | BCR_ABL=="primers")

my.anno$Stage_2[my.anno$Stage_2=="pre_blast_crisis" & my.anno$Patient=="CML1266"]<-"pre_blast_crisis_CML1266"

###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
###############################################################################
min_expr<-1
S.m<-as.matrix(my.cells)
S.m[S.m < min_expr]<- 0
S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
#################################With Batch effect#############################
my_bcr_abl<-as.character(my.anno$Stage_2)
my.design <- model.matrix(~my_bcr_abl)
rownames(my.design)<-colnames(S.m)

batch.x<-as.character(my.anno$process_date)
batch.x[which(batch.x=="batch1" | batch.x=="batch2")]<-"A"
batch.x[which(batch.x=="batch3")]<-"B"
batch.x[which(batch.x=="batch4")]<-"C"
batch.x[which(batch.x=="batch5")]<-"D"
batch.x[which(batch.x=="batch6")]<-"E"

S.m.no.batch = removeBatchEffect(S.m,batch=batch.x,design=my.design)
###############################################################################
dat1<-read.table("DEGenes-Diagnosis+_Vs_BC+.txt",header=T)
dat2<-read.table("DEGenes-Diagnosis+_Vs_preBC+.txt",header=T)
dat3<-read.table("DEGenes-preBC+_Vs_BC+.txt",header=T)
dat4<-read.table("DEGenes-Diagnosis+_Vs_Normal_HSCs.txt",header=T)
###############################################################################
dat1.sig<-subset(dat1,abs(log2fc) >=2.5 & p.adjust < 0.05)
dat2.sig<-subset(dat2,abs(log2fc) >=2 & p.adjust < 0.05)
dat3.sig<-subset(dat3,abs(log2fc) >=2 & p.adjust < 0.05)
dat4.sig<-subset(dat4,abs(log2fc) >=1.5 & p.adjust < 0.05)
###############################################################################
dat1.genes<-as.character(dat1.sig$gene)
dat2.genes<-as.character(dat2.sig$gene)
dat3.genes<-as.character(dat3.sig$gene)
dat4.genes<-as.character(dat4.sig$gene)
################################################################################
my.genes<-c(dat1.genes,dat2.genes,dat3.genes,dat4.genes)
my.genes<-unique(my.genes)
###############################################################################	
S.m2<-subset(S.m.no.batch,rownames(S.m.no.batch) %in% my.genes)
###############################################################################
write.table(S.m2,file="Normal_HSCs_Diagnosis_K562_PreBC_BC.data.for.TSNE.txt",append=FALSE, 
		    sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
###############################################################################