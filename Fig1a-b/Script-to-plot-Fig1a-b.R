#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig1a and Fig1b.##########
#############################################################
library(beeswarm)
##############################################################
my.dat<-read.table("K562.qPCR.data.txt",header=T)
##############################################################
#####LOD40 was used in these two cases to measure and to compare the sensitivity
#####of detection between Smart-seq2 and the tSS2 protocols.## 
LOD<-40
my.dat$log2ex<-LOD-my.dat$ct_value
my.dat$log2ex[my.dat$ct_value==0 | my.dat$ct_value > LOD] <-0
##############################################################
##############################################################
sm2<-subset(my.dat,method=="sm2")
sm2.1<-subset(sm2,status=="no_primers")

my.g1<-subset(sm2.1,gene=="BCR_ABL")
my.g2<-subset(sm2.1,gene=="GAPDH")

my.list1<-list()
my.list1[["BCR-ABL"]]<- as.numeric(my.g1$log2ex)
my.list1[["GAPDH"]]<- as.numeric(my.g2$log2ex)

################
pdf(file="Fig1a.pdf", width=3, height=5, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
beeswarm(my.list1, corral="wrap",pch = 19,cex=0.50,col = c("blue","blue"),
		,ylab="Expression level (relative to LOD)",main="Smart-seq2",las=2,ylim=c(0,30))
abline(h=0,lty=2,lwd=2,col="grey")
boxplot(my.list1, add = T,col="#0000ff22",axes=FALSE)

dev.off()
################
sm2.2<-subset(sm2,status=="with_primers")

my.g1<-subset(sm2.2,gene=="BCR_ABL")
my.g2<-subset(sm2.2,gene=="GAPDH")

my.list2<-list()
my.list2[["BCR-ABL"]]<- as.numeric(my.g1$log2ex)
my.list2[["GAPDH"]]<- as.numeric(my.g2$log2ex)

################
pdf(file="Fig1b.pdf", width=3, height=5, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
beeswarm(my.list2, corral="wrap",pch = 19,cex=0.50,col = c("brown","brown"),
		,ylab="Expression level (relative to LOD)",main="BCR-ABL tSS2",las=2,ylim=c(0,30))
abline(h=0,lty=2,lwd=2,col="grey")
boxplot(my.list2, add = T,col="#0000ff22",axes=FALSE)

dev.off()
################