###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate data for the GSEA analysis 
#(Normal HSCs Vs Remission ClassA at 3 months TKI).
###############################################################################
my.anno<-read.delim("tSNE_NormalDiagnosisRemission_500_No1M_TKI.FINAL_V4_with_Kmeans_Assignment.txt",header=T)
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
###############################################################################
###############################################################################
S.m<-as.matrix(my.cells)
###############################################################################
###############################################################################
groupA.cells<-subset(my.anno,remission_class=="A")
groupA.cells<-subset(groupA.cells,Stage_1=="3_months_TKI")
groupA<-as.character(groupA.cells$Cell)

groupB.cells<-subset(my.anno,Stage_1=="normal_hsc")
groupB<-as.character(groupB.cells$Cell)
##########################Get All Data#########################################
groupA.m<-S.m[,colnames(S.m) %in% groupA]
groupB.m<-S.m[,colnames(S.m) %in% groupB]

my.phenotype.s1<-rep("3_MONTHS_TKI",ncol(groupA.m))
my.phenotype.s2<-rep("NormalHSC",ncol(groupB.m))

my.phenotype<-c(my.phenotype.s1,my.phenotype.s2)

my.data<-cbind(groupA.m,groupB.m)

my.info<-data.frame(NAME=rownames(my.data))
my.info$DESCRIPTION<-"na"

my.final<-cbind(my.info,my.data)
###############################################################################
##############################################################################
h1<-paste(ncol(my.data),"2","1",sep=" ")
h2<-paste("#","3_MONTHS_TKI","NormalHSC",sep=" ")
h3<-paste(c(rep("3_MONTHS_TKI",length(groupA)),rep("NormalHSC",length(groupB)),sep=" "))

cat(h1,file=paste("3_MONTHS_Vs_NormalHSC-GSEA-phenotype",".cls",sep=""),sep="\n")
cat(h2,file=paste("3_MONTHS_Vs_NormalHSC-GSEA-phenotype",".cls",sep=""),sep="\n",append=TRUE)
cat(h3,file=paste("3_MONTHS_Vs_NormalHSC-GSEA-phenotype",".cls",sep=""),append=TRUE)

write.table(my.final,file="3_MONTHS_Vs_NormalHSC-GSEA-format.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
