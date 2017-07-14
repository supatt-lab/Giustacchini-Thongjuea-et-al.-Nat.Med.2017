#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig2g.
#############################################################
library(matrixStats)
library(beeswarm)


my.anno<-read.delim("CML-PROJECT-DATA_SET1-INFORMATION-Freezed-3.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Patient=="OX1407" | Stage=="normal_hsc")
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
###############################################################################
min_expr<-1
S.m<-as.matrix(my.cells)
S.m[S.m < min_expr]<- 0
S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
###############################################################################
###############################################################################
my.genes<-c("CLU","FCER1A","GAS2","MZB1","RGS2","CXCR4")

S.m<-S.m[rownames(S.m) %in% my.genes,]

################save data matrix to a file####################################
write.table(S.m,file="Fig2g.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)
##############################################################################


groupA<-as.character(my.anno$Cell[grep("negative",my.anno$BCR_ABL)])
groupA.m<-S.m[,colnames(S.m) %in% groupA]

groupB<-as.character(my.anno$Cell[grep("positive",my.anno$BCR_ABL)])
groupB.m<-S.m[,colnames(S.m) %in% groupB]

groupC<-as.character(my.anno$Cell[grep("normal",my.anno$BCR_ABL)])
groupC.m<-S.m[,colnames(S.m) %in% groupC]

pdf(file="Fig2g.pdf", width=12, height=6, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)


par(mfrow=c(2,6),mar=c(4,4,2,2),oma = c(2, 2, 1, 1))

for(i in 1:length(my.genes)){
	
	my.gene<-my.genes[i]
	my.groupA.m<-groupA.m[which(rownames(groupA.m) %in% my.gene),]
	my.meanA<-mean(my.groupA.m)
	my.groupA.sd<-sd(my.groupA.m)
	my.groupA.box<-my.groupA.m[my.groupA.m !=0]
	my.groupA.mean2<-mean(my.groupA.box)
	my.groupA.sd2<-sd(my.groupA.box)
	
	my.groupB.m<-groupB.m[which(rownames(groupB.m) %in% my.gene),]
	my.meanB<-mean(my.groupB.m)
	my.groupB.sd<-sd(my.groupB.m)
	
	my.groupB.box<-my.groupB.m[my.groupB.m !=0]
	my.groupB.mean2<-mean(my.groupB.box)
	my.groupB.sd2<-sd(my.groupB.box)
	
	
	my.groupC.m<-groupC.m[which(rownames(groupC.m) %in% my.gene),]
	my.meanC<-mean(my.groupC.m)
	my.groupC.sd<-sd(my.groupC.m)
	
	my.groupC.box<-my.groupC.m[my.groupC.m !=0]
	my.groupC.mean2<-mean(my.groupC.box)
	my.groupC.sd2<-sd(my.groupC.box)
	
	##################Get present/total##############
	A.present<-length(my.groupA.m[my.groupA.m > 0])
	A.total<-length(my.groupA.m)
	
	B.present<-length(my.groupB.m[my.groupB.m > 0])
	B.total<-length(my.groupB.m)
	
	C.present<-length(my.groupC.m[my.groupC.m > 0])
	C.total<-length(my.groupC.m)
	#################################################
	text.A<-paste(A.present,A.total,sep="/")
	text.B<-paste(B.present,B.total,sep="/")
	text.C<-paste(C.present,C.total,sep="/")
	my.text<-c(text.A,text.B,text.C)
	
	my.list<-list()
	my.list[["BCR-ABL-"]]<- my.groupA.m
	my.list[["BCR-ABL+"]]<- my.groupB.m
	my.list[["Normal HSC"]]<- my.groupC.m
	
	beeswarm(my.list, corral="wrap",pch = 19,cex=0.6,col = c("blue","brown","black")
			,ylab="Expression level (Log2 RPKM)",main=my.gene,las=2,ylim=c(-2,13))
	boxplot(my.list, add = T,col="#0000ff22",axes=FALSE)
	text(c(1:3),-1,labels=my.text,cex=0.80)
	points(c(my.meanA,my.meanB,my.meanC), pch = 22, col = "red", lwd = 2)
	
}
dev.off()