###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig5g.
###############################################################################
library(beeswarm)

my.anno<-read.delim("tSNE_NormalDiagnosisRemission_500_No1M_TKI.FINAL_V4_with_Kmeans_Assignment.txt",header=T)
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
###############################################################################
min_expr<-1
S.m<-as.matrix(my.cells)
S.m[S.m < min_expr]<- 0
S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
###############################################################################
groupA.cells<-subset(my.anno,Stage_1=="normal_hsc")
groupA<-as.character(groupA.cells$Cell)
groupA.m<-S.m[,colnames(S.m) %in% groupA]

groupB.cells<-subset(my.anno,Stage_2=="diagnosis" & BCR_ABL=="positive")
groupB<-as.character(groupB.cells$Cell)
groupB.m<-S.m[,colnames(S.m) %in% groupB]

groupC.cells<-subset(my.anno,remission_class=="A")
groupC<-as.character(groupC.cells$Cell)
groupC.m<-S.m[,colnames(S.m) %in% groupC]

#################
dat1<-read.table("DEGenes-Normal_HSCs_Vs_Remission_classA.txt",header=T)
#################
#################
my.genes<-c("GAS2","CTNNB1","SKIL","NFKBIA","SQSTM1","HIF1A","WTAP","CXCR4","FOS","CD53")
#################
S.m<-S.m[rownames(S.m) %in% my.genes,]
my.diff.stat<-dat1[dat1$gene %in% my.genes,]
################save data matrix to a file####################################
write.table(S.m,file="Fig5g.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)

write.table(my.diff.stat,file="Fig5g.statistical.pvalues.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
##############################################################################
pdf(file="Fig5g.pdf", width=12, height=6, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)

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
	
	my.list<-list()
	my.list[[""]]<- my.groupA.m
	my.list[[""]]<- my.groupB.m
	my.list[[""]]<- my.groupC.m
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
	
	beeswarm(my.list, corral="wrap",pch = 19,cex=0.6,col = c("gray35","brown","cyan")
			,ylab="Expression (RPKM, Log2)",main=my.gene,cex.main=1.5,cex.lab=1.2,cex.axis=1.2,las=2,ylim=c(-2,14))
	boxplot(my.list, add = T,col="#0000ff22",axes=FALSE)
	text(c(1:3),-1,labels=my.text,cex=0.7)
	points(c(my.meanA,my.meanB,my.meanC), pch = 22, col = "red", lwd = 2)
	
}
dev.off()