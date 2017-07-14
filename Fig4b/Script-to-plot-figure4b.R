###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig4a.
###############################################################################
my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Stage_2=="diagnosis" & (BCR_ABL=="negative" | BCR_ABL=="negative_low_gapdh"))

my.tsne<-read.table("tSNE.5611.genes.perplexity20.txt",header=T)
my.data<-read.table("Diagnosis_BCR_ABL-_5611_Genes.for_tSNE.txt",header=T)

################################################################################
my.tsne$Cell<-colnames(my.data)

my.dat<-merge(my.tsne,my.anno,by.x="Cell",by.y="Cell",all=T)

my.dat<-subset(my.dat,Responder=="good" | Responder=="poor")

my.dat$my.color<-"red"
my.dat$my.color[my.dat$Responder=="good"]<-"blue"

################save data matrix to a file###############
write.table(my.dat[c("Cell","Dim1","Dim2","BCR_ABL","Responder","Patient")],file="Fig4a.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
#########################################################
pdf(file="Fig4b.pdf", width=7, height=6, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
plot(my.dat$Dim1,my.dat$Dim2,
		xlab = "Using 5,611 variable genes",
		ylab = "",
		main="356 Diagnosis (BCR-ABL-) cells",col = my.dat$my.color,cex=1.5,pch=19,yaxt='n',xaxt='n')

legend("topright",c("Poor responder","Good responder"),col=c("red","blue"),pch=c(19),cex=0.9)
###############################################################################
dev.off()