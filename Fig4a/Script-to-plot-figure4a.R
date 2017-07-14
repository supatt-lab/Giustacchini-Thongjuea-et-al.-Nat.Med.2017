###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig_4a.####################################
###############################################################################
my.dat<-read.delim("tSNE_Diagnosis+_good_vs_poor_perplexity20.Result.txt",header=T)
my.dat<-subset(my.dat,Responder=="good" | Responder=="poor")

my.dat$my.color<-"red"
my.dat$my.color[my.dat$Responder=="good"]<-"blue"

################save data matrix to a file######################################
#write.table(my.dat[c("Cell","Dim1","Dim2","BCR_ABL","Responder")],
#		file="Fig_S11.data.matrix.txt",append=FALSE, sep="\t", 
#		quote=FALSE,row.names=FALSE, col.names=TRUE)


pdf(file="Fig_4a.pdf", width=7, height=7, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
plot(my.dat$Dim1,my.dat$Dim2,
		xlab = "Using 5,011 variable genes",
		ylab = "",
		main="436 Diagnosis (BCR-ABL+) cells",col = my.dat$my.color,cex=1.5,pch=19,yaxt='n',xaxt='n')

legend("bottomright",c("Poor responder","Good responder"),col=c("red","blue"),pch=c(19),cex=1)
###############################################################################
dev.off()
