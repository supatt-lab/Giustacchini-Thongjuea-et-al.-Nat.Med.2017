#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig2a.
#############################################################
my.info<-read.table("HSC_Sampling_Count_Genes_Per_Cell_1_RPKM.txt",header=T)
g01<-grep("_01",rownames(my.info))
g02<-grep("_02",rownames(my.info))
g03<-grep("_03",rownames(my.info))
g04<-grep("_04",rownames(my.info))
g05<-grep("_05",rownames(my.info))
g1<-grep("_1",rownames(my.info))
g2<-grep("_2",rownames(my.info))
g3<-grep("_3",rownames(my.info))
g4<-grep("_4",rownames(my.info))
g5<-grep("_5",rownames(my.info))
g6<-grep("_6",rownames(my.info))

##############################################
x01<-my.info[g01,]
x02<-my.info[g02,]
x03<-my.info[g03,]
x04<-my.info[g04,]
x05<-my.info[g05,]
x1<-my.info[g1,]
x2<-my.info[g2,]
x3<-my.info[g3,]
x4<-my.info[g4,]
x5<-my.info[g5,]
x6<-my.info[g6,]
##############################################
my.list<-list()
my.list[["0.1"]]<- x01
my.list[["0.2"]]<- x02
my.list[["0.3"]]<- x03
my.list[["0.4"]]<- x04
my.list[["0.5"]]<- x05
my.list[["1"]]<- x1
my.list[["2"]]<- x2
my.list[["3"]]<- x3
my.list[["4"]]<- x4
my.list[["5"]]<- x5
my.list[["6"]]<- x6

pdf(file="Fig2a.pdf", width=6, height=5, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
boxplot(my.list,col="#0000ff22",ylim=c(1000,5500),xlab="Million reads / cell",
		ylab="Number of detected genes (RPKM >=1)",
		main="Normal HSCs (12 cells)",cex=0.2)
dev.off()

