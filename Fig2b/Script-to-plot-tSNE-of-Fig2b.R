#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig2b.
##############################################################
my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,(Patient=="K562" & BCR_ABL=="primers") | Stage_2=="normal_hsc")

my.tsne<-read.table("tSNE.pp20.pca20.result.txt",header=T)
my.data<-read.table("NORMAL-HSCs-and-K562.matrix.txt",header=T)

###############################################################
my.tsne$Cell<-colnames(my.data)

my.dat<-merge(my.tsne,my.anno,by.x="Cell",by.y="Cell",all=T)

my.dat$my.color<-"blue"
my.dat$my.color[my.dat$Patient=="K562"]<-"purple"

my.dat$my.shape[my.dat$Patient=="K562"]<-16
my.dat$my.shape[my.dat$Patient=="NBM_89"]<-15
my.dat$my.shape[my.dat$Patient=="NBM_91"]<-10
my.dat$my.shape[my.dat$Patient=="NBM_ERCC"]<-10
my.dat$my.shape[my.dat$Patient=="NBM_petter"]<-8
my.dat$my.shape[my.dat$Patient=="NBM_lonza_age_22"]<-5
my.dat$my.shape[my.dat$Patient=="NBM_NC"]<-9

################save the data with adding tSNE of Dim1 and Dim2################
write.table(my.dat[c("Cell","Dim1","Dim2","BCR_ABL","Stage_1","my.color","my.shape")],
		file="Fig2b_data.matrix.txt",append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
###############################################################################
pdf(file="Fig2b.pdf", width=7, height=6, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
plot(my.dat$Dim1,my.dat$Dim2,xaxt="n",yaxt="n",xlab=paste("Number of genes:",nrow(my.data),sep=""),ylab="",
		pch = my.dat$my.shape,
		main="",col = my.dat$my.color,cex=1.2)
legend("bottomright",c("NBM_89","NBM_91","NBM_PT","NBM_Lonza22","NBM_NC"),pch=c(15,10,8,5,9),cex=1,col="blue")

dev.off()