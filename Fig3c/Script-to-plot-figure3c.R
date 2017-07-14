###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig3c.
###############################################################################
my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Stage_2=="normal_hsc"| Stage_2=="diagnosis")

my.tsne<-read.table("tSNE.8589.genes.txt",header=T)
my.data<-read.table("Normal_HSCs_and_Diagnosis_Cells.for.TSNE.txt",header=T)

################################################################################
my.tsne$Cell<-colnames(my.data)

my.dat<-merge(my.tsne,my.anno,by.x="Cell",by.y="Cell",all=T)


my.dat$my.color<-"blue"
my.dat$my.color[my.dat$BCR_ABL=="negative_low_gapdh"]<-"blue"
my.dat$my.color[my.dat$BCR_ABL=="positive"]<-"red"
my.dat$my.color[my.dat$BCR_ABL=="normal"]<-"black"


x1<-which(my.dat$my.color=="blue")
x2<-which(my.dat$my.color=="red")
x3<-which(my.dat$my.color=="black")

################save data matrix to a file####################
write.table(my.dat[c("Cell","Dim1","Dim2","BCR_ABL","my.color")],file="Fig3c.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
################################################################################
pdf(file="Fig3c.pdf", width=8, height=8, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
plot(my.dat$Dim1,my.dat$Dim2,xaxt="n",yaxt="n",xlab="",ylab="",type="n",
		main="",col = my.dat$my.color,cex=1.5)

points(my.dat$Dim1[x1],my.dat$Dim2[x1],pch=23,bg="blue",cex=1.5,col="black")
points(my.dat$Dim1[x2],my.dat$Dim2[x2],pch=24,bg ="red",col="black",cex=1.5)
points(my.dat$Dim1[x3],my.dat$Dim2[x3],pch=21,bg ="gray35",col="black",cex=1.5)

legend("topleft",c("Normal HSC","BCR_ABL-","BCR_ABL+"),pch=19,col=c("gray35","blue","red"))
dev.off()