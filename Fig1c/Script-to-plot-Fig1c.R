#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig1c.
#############################################################
my.dat<-read.table("Fig1c.data.txt",header=T)

x<-my.dat$meanA
y<-my.dat$meanB

col="blue"

cor.txt<-cor(x,y)
cor.txt<-format(cor.txt,digit=4)

pdf(file="Fig1c.pdf", width=6, height=6, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)

plot(x,y,xlab="Smart-seq2, mean of Log2(RPKM)",ylab="BCR-ABL tSS2, mean of Log2(RPKM)",xlim=c(0,14),ylim=c(0,14),
		pch = 19,
		main="",col = col,cex=0.8)

text(10,0,paste("Pearson correlation=",cor.txt,sep=""))

dev.off()