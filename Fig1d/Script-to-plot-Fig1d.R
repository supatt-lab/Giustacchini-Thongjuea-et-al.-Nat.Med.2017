#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig1d.
#############################################################
#############################################################
my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Patient=="K562")

my.tsne<-read.table("tSNE_pp20_pca20_result.txt",header=T)
my.data<-read.table("CML4TSNE-K562.txt",header=T)

###############################################################################
my.tsne$Cell<-colnames(my.data)

my.dat<-merge(my.tsne,my.anno,by.x="Cell",by.y="Cell",all=T)

my.dat$my.color<-"red"
my.dat$my.color[my.dat$BCR_ABL=="no_primers"]<-"blue"
################save the data with adding tSNE of Dim1 and Dim2################
write.table(my.dat[c("Cell","Dim1","Dim2","my.color")],
		file="Fig1d_data.matrix.txt",append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
###############################################################################
pdf(file="Fig1d.pdf", width=6, height=6, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
plot(my.dat$Dim1,my.dat$Dim2,
		xlab = paste("Number of genes : ",nrow(my.data),sep=""),
		ylab = "",
		main="",col = my.dat$my.color,cex=2,pch=19,yaxt='n',xaxt='n')
points(my.dat$Dim1,my.dat$Dim2,col="black",cex=2)

legend('topright',pch=c(19),cex=1,col=c("red","blue"),c("BCR-ABL tSS2","Smart-seq2"))

dev.off()

