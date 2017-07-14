###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
###############################################################################
###############################################################################
S.m<-as.matrix(my.cells)
###############################################################################
###############################################################################
groupA.cells<-subset(my.anno,Stage_2=="normal_hsc")
groupA<-as.character(groupA.cells$Cell)

groupB.cells<-subset(my.anno,BCR_ABL=="positive" & Stage_2=="diagnosis")
groupB<-as.character(groupB.cells$Cell)
##########################Get All Data#########################################
groupA.m<-S.m[,colnames(S.m) %in% groupA]
groupB.m<-S.m[,colnames(S.m) %in% groupB]

my.phenotype.s1<-rep("Normal_HSCs",ncol(groupA.m))
my.phenotype.s2<-rep("Diagnosis+",ncol(groupB.m))

my.phenotype<-c(my.phenotype.s1,my.phenotype.s2)

my.data<-cbind(groupA.m,groupB.m)

my.info<-data.frame(NAME=rownames(my.data))
my.info$DESCRIPTION<-"na"

my.final<-cbind(my.info,my.data)
###############################################################################
##############################################################################
h1<-paste(ncol(my.data),"2","1",sep=" ")
h2<-paste("#","Normal_HSCs","Diagnosis+",sep=" ")
h3<-paste(c(rep("Normal_HSCs",length(groupA)),rep("Diagnosis+",length(groupB)),sep=" "))

cat(h1,file=paste("NormalHSCs-Diagnosis+_GSEA-phenotype",".cls",sep=""),sep="\n")
cat(h2,file=paste("NormalHSCs-Diagnosis+_GSEA-phenotype",".cls",sep=""),sep="\n",append=TRUE)
cat(h3,file=paste("NormalHSCs-Diagnosis+_GSEA-phenotype",".cls",sep=""),append=TRUE)

write.table(my.final,file="NormalHSCs-Diagnosis+_GSEA.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
