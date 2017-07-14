###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig5a.
###############################################################################
my.dat<-read.delim("tSNE_NormalDiagnosisRemission_500_No1M_TKI.FINAL_V4-with-Kmean-Assignment.txt",header=T)
########################
my.dat$my.color[my.dat$BCR_ABL=="normal"]<-"black"
my.dat$my.color[my.dat$BCR_ABL=="positive" & my.dat$Stage_2=="diagnosis"]<-adjustcolor( "gray", alpha.f = 0.2)
my.dat$my.color[my.dat$remission_class=="A"]<-adjustcolor( "cyan", alpha.f = 0.95)
my.dat$my.color[my.dat$remission_class=="B"]<-adjustcolor( "blue", alpha.f = 0.95)


my.dat$type[my.dat$BCR_ABL=="normal"]<-19
my.dat$type[my.dat$BCR_ABL=="positive" & my.dat$Stage_2=="diagnosis"]<-19
my.dat$type[my.dat$remission_class=="A"]<-18
my.dat$type[my.dat$remission_class=="B"]<-17


x1<-which(my.dat$BCR_ABL=="positive" & my.dat$Stage_2=="diagnosis")
x2<-which(my.dat$remission_class=="B")
x3<-which(my.dat$remission_class=="A")


################save data matrix to a file####################################
write.table(my.dat[c("Cell","Dim1","Dim2","BCR_ABL","Patient","Stage_1",
						"remission_class","my.color","type")],
		file="Fig5a.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
##############################################################################
pdf(file="Fig5a.pdf", width=7, height=7, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
plot(my.dat$Dim1,my.dat$Dim2,pch =my.dat$type,
		xlab = "",
		ylab = "",
		main="",col = as.character(my.dat$my.color),cex=1.5,yaxt='n',xaxt='n')
points(my.dat$Dim1[x2],my.dat$Dim2[x2],pch=24,bg ="blue",col="black",cex=1.5)
points(my.dat$Dim1[x3],my.dat$Dim2[x3],pch=23,bg ="cyan",col="black",cex=1.5)

legend("bottomright",c("Normal HSC","BCR-ABL+ (diagnosis)","BCR-ABL+ (remission) class A","BCR-ABL+ (remission) class B")
		,col=c("black","gray","cyan","blue"),pch=c(19,19,18,17),cex=0.85)

dev.off()
####################