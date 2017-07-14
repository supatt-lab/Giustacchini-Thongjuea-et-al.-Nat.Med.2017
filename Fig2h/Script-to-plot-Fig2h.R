#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig2h.
##########################Required R packages################
library(data.table)
library(beeswarm)
library(reshape2)
library(MASS)
library(pheatmap)
library(Biobase)
library(Hmisc)
##############################################################
source("../Data/SingleCellBiomarkFileFormat.R")
source("../Data/SingleCellClasses.R")
source("../Data/SingleCellDataManipulations.R")
source("../Data/SingleCellGenerics.R")
source("../Data/SingleCellPlots.R")
##############################################################
my.genes<-c("CLU","FCER1A","GAS2","MZB1","RGS2","CXCR4")

LOD=30

my.obj<-new("SingleCellFluidigm",
		biomark_files=c("20150127_OX1407.heatmap.csv","OX1407_plate2.csv","NBM91.csv"),
		chip_names=c("plate1","plate2","normal_hsc"),
		data_dir="Biomark_files/",
		LOD=LOD)

getTableOfCtValues(my.obj)
normalizeCtValuesByHouseKeepingGenes(my.obj,houseKeeping.genes=c("GAPDH","B2M"))

##############################################################
my.expr<-my.obj@org.data
my.expr$delta.ct[is.na(my.expr$raw_ct)==T]<-0
##############################################################
getBCR_ABL_Annotation<-function(obj,bcr_abl="BCR-ABL1 FUSION"){
	obj<-my.expr
	my.dat<-subset(obj,gene==bcr_abl)
	my.dat$BCR_ABL_status<-"positive"
	my.dat$BCR_ABL_status[is.na(my.dat$raw_ct)==T]<-"negative"
	my.dat$BCR_ABL_status[my.dat$chip=="normal_hsc"]<-"normal_hsc"
	my.anno<-my.dat["BCR_ABL_status"]
	rownames(my.anno)<-my.dat$cell.name
	return(my.anno)
}
##############################################################
my.cell.anno<-getBCR_ABL_Annotation(my.expr)

groupA.anno<-subset(my.cell.anno,BCR_ABL_status=="negative")
groupA.cells<-rownames(groupA.anno)

groupB.anno<-subset(my.cell.anno,BCR_ABL_status=="positive")
groupB.cells<-rownames(groupB.anno)

groupC.anno<-subset(my.cell.anno,BCR_ABL_status=="normal_hsc")
groupC.cells<-rownames(groupC.anno)


m.m <- dcast(my.expr, gene~cell.name,value.var="delta.ct",na.rm=T)
x.m<-as.matrix(m.m[,2:ncol(m.m)])
rownames(x.m)<-m.m$gene

################save data matrix to a file####################
write.table(x.m,file="Fig2h.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)
##############################################################
groupA.m<-x.m[,colnames(x.m) %in% groupA.cells]
groupB.m<-x.m[,colnames(x.m) %in% groupB.cells]
groupC.m<-x.m[,colnames(x.m) %in% groupC.cells]
##############################################################
pdf(file="Fig2h.pdf", width=12, height=6, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)

par(mfrow=c(2,6),mar=c(4,4,2,2),oma = c(2, 2, 1, 1))
#################plot each gene###############################
for(i in 1:length(my.genes)){
	
	my.gene<-my.genes[i]
	my.groupA.m<-groupA.m[which(rownames(groupA.m) %in% my.gene),]
	my.groupA.sd<-sd(my.groupA.m)
	my.groupA.exp<-my.groupA.m[my.groupA.m !=0]
	my.groupA.mean2<-mean(my.groupA.exp)
	my.groupA.sd2<-sd(my.groupA.exp)
	
	my.groupB.m<-groupB.m[which(rownames(groupB.m) %in% my.gene),]
	my.groupB.sd<-sd(my.groupB.m)
	my.groupB.exp<-my.groupB.m[my.groupB.m !=0]
	my.groupB.mean2<-mean(my.groupB.exp)
	my.groupB.sd2<-sd(my.groupB.exp)
	
	my.groupC.m<-groupC.m[which(rownames(groupC.m) %in% my.gene),]
	my.groupC.sd<-sd(my.groupC.m)
	my.groupC.exp<-my.groupC.m[my.groupC.m !=0]
	my.groupC.mean2<-mean(my.groupC.exp)
	my.groupC.sd2<-sd(my.groupC.exp)
	
	min.value<-min(my.groupA.m,my.groupB.m,my.groupC.m)
	my.groupA.m[my.groupA.m==0]<-min.value-3.5
	my.groupB.m[my.groupB.m==0]<-min.value-3.5
	my.groupC.m[my.groupC.m==0]<-min.value-3.5
	##################Get present/total##############
	A.present<-length(my.groupA.m[my.groupA.m > min.value-3.5])
	A.total<-length(my.groupA.m)
	
	B.present<-length(my.groupB.m[my.groupB.m > min.value-3.5])
	B.total<-length(my.groupB.m)
	
	C.present<-length(my.groupC.m[my.groupC.m > min.value-3.5])
	C.total<-length(my.groupC.m)
	#################################################
	text.A<-paste(A.present,A.total,sep="/")
	text.B<-paste(B.present,B.total,sep="/")
	text.C<-paste(C.present,C.total,sep="/")
	my.text<-c(text.A,text.B,text.C)
	
	my.medianA<-median(my.groupA.m)
	my.medianB<-median(my.groupB.m)
	my.medianC<-median(my.groupC.m)
	
	my.meanA<-mean(my.groupA.m)
	my.meanB<-mean(my.groupB.m)
	my.meanC<-mean(my.groupC.m)
	
	
	my.list<-list()
	my.list[["BCR-ABL-"]]<- my.groupA.m
	my.list[["BCR-ABL+"]]<- my.groupB.m
	my.list[["Normal HSC"]]<- my.groupC.m
	
	beeswarm(my.list, corral="wrap",pch = 19,cex=0.6,col = c("blue","brown","black"),
			,ylab="Relative expression",main=my.gene,las=2)
	abline(h=min.value-3,lty=2,lwd=1.5,col="grey")
	text(c(1:3),min.value-4,labels=my.text,cex=0.75)
	boxplot(my.list, add = T,col="#0000ff22",axes=FALSE)
	points(c(my.meanA,my.meanB,my.meanC), pch = 22, col = "red", lwd = 2)
	
}
dev.off()
