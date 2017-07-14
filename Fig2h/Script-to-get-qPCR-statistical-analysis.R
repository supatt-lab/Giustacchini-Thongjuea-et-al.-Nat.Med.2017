#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to run differentially expressed gene analysis of qPCR data.
#The script will show only the statistical analysis results of selected genes.
##########################Required R packages###############
library(data.table)
library(beeswarm)
library(reshape2)
library(MASS)
library(pheatmap)
library(Biobase)
library(gplots)
library(RColorBrewer)
###############################################################################
source("../Data/SingleCellBiomarkFileFormat.R")
source("../Data/SingleCellClasses.R")
source("../Data/SingleCellDataManipulations.R")
source("../Data/SingleCellGenerics.R")
source("../Data/SingleCellPlots.R")
###############################################################################
LOD=30

my.obj<-new("SingleCellFluidigm",
		biomark_files=c("20150127_OX1407.heatmap.csv","OX1407_plate2.csv"),
		chip_names=c("plate1","plate2"),
		data_dir="Biomark_files/",
		LOD=LOD)

getTableOfCtValues(my.obj)
normalizeCtValuesByHouseKeepingGenes(my.obj,houseKeeping.genes=c("GAPDH","B2M"))

my.expr<-my.obj@org.data
##############################################################
getBCR_ABL_Annotation<-function(obj,bcr_abl="BCR-ABL1 FUSION"){
	obj<-my.expr
	my.dat<-subset(obj,gene==bcr_abl)
	my.dat$BCR_ABL_status<-"positive"
	my.dat$BCR_ABL_status[is.na(my.dat$raw_ct)==T]<-"negative"
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


m.m <- dcast(my.expr, gene~cell.name,value.var="delta.ct",na.rm=T)
x.m<-as.matrix(m.m[,2:ncol(m.m)])
rownames(x.m)<-m.m$gene

groupA.m<-x.m[,colnames(x.m) %in% groupA.cells]
groupB.m<-x.m[,colnames(x.m) %in% groupB.cells]
##############################################################
groupA.mean<-data.frame(groupA_mean=rowMeans(groupA.m,na.rm=T))
groupA.mean$Gene<-rownames(groupA.mean)
groupB.mean<-data.frame(groupB_mean=rowMeans(groupB.m,na.rm=T))
groupB.mean$Gene<-rownames(groupB.mean)

my.f<-merge(groupA.mean,groupB.mean,all=T)

qPCR.df <- data.frame(gene=gene<-my.f$Gene,qPCR.neg = my.f$groupA_mean,qPCR.pos = my.f$groupB_mean)
qPCR.df$delta_delta_ct<-qPCR.df$qPCR.pos-qPCR.df$qPCR.neg
##############################################################
##############################################################
z.m <- dcast(my.expr, gene~cell.name,value.var="log2ex",na.rm=T)
t.m<-as.matrix(z.m[,2:ncol(z.m)])
rownames(t.m)<-z.m$gene
##############################################################
groupA.m<-t.m[,colnames(t.m) %in% groupA.cells]
groupB.m<-t.m[,colnames(t.m) %in% groupB.cells]

fishersMethod <-function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)

result.f<-data.frame()

for(i in 1:nrow(groupA.m)){
	
	n<-groupA.m[i,]
	m<-groupB.m[i,]
	n.f<-n[n>0]
	m.f<-m[m>0]
	
	x<-length(n)
	y<-length(m)
	
	x.1<-length(n.f)
	y.1<-length(m.f)
	x.2<-x-x.1
	y.2<-y-y.1
	
	m.m <- matrix(c(x.1, x.2, y.1, y.2), ncol = 2)
	fisher.test.p <-fisher.test(m.m)$p.value
	
	wilcox.p<-wilcox.test(n,m,paired = FALSE)$p.value
	
	all.p<-c(fisher.test.p,wilcox.p)
	fisher.p<-fishersMethod (all.p)
	
	my.test.f<-data.frame(nCellA=x,nCellB=y,expCellA=x.1,expCellB=y.1,
			expFractionCellA=x.1/x,expFractionCellB=y.1/y,fisher.test.p=fisher.test.p,
			wilcox=wilcox.p,fisher=fisher.p)
	my.combined<-cbind(qPCR.df[i,],my.test.f)
	result.f<-rbind(result.f,my.combined)
}

my.final.result<-result.f[order(result.f$fisher),]
###############################################################
my.genes<-c("CLU","FCER1A","GAS2","MZB1","RGS2","CXCR4")
genes.statitical.info<-my.final.result[my.final.result$gene %in% my.genes,]

write.table(genes.statitical.info,file="Fig2h.statistical.values.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)