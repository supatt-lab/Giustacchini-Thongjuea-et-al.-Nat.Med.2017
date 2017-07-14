###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig6c.
###############################################################################
###############################################################################
my.dat<-read.delim("tSNE_Diagnosis-PreBC-BC-Without-CML1600_Result.txt",header=T)

my.dat$my.color[my.dat$Stage_2=="normal_hsc"]<-adjustcolor( "gray35", alpha.f = 0.2)
my.dat$my.color[my.dat$Stage_2=="cell_line"]<-adjustcolor( "sienna", alpha.f = 0.2)
my.dat$my.color[my.dat$Stage_2=="diagnosis" | my.dat$Stage_2=="accelerated_phase"]<-adjustcolor( "brown", alpha.f = 0.2)
my.dat$my.color[my.dat$Stage_2=="pre_blast_crisis" & my.dat$Patient=="OX1931"]<-"orange"
my.dat$my.color[my.dat$Stage_2=="pre_blast_crisis_CML1266"]<-adjustcolor( "darkgreen", alpha.f = 0.1)
my.dat$my.color[my.dat$Stage_2=="blast_crisis" & my.dat$Patient=="CML1266"]<-adjustcolor( "purple", alpha.f = 0.1)
my.dat$my.color[my.dat$Stage_2=="blast_crisis" & my.dat$Patient=="OX2046"]<-"cyan"
my.dat$my.color[my.dat$Stage_2=="blast_crisis_1_month_TKI"]<-adjustcolor( "deeppink", alpha.f = 0.1)


my.dat$type[my.dat$Stage_2=="normal_hsc"]<-19
my.dat$type[my.dat$Stage_2=="cell_line"]<-19
my.dat$type[my.dat$Stage_2=="diagnosis" | my.dat$Stage_2=="accelerated_phase"]<-17
my.dat$type[my.dat$Stage_2=="pre_blast_crisis" & my.dat$Patient=="OX1931"]<-23
my.dat$type[my.dat$Stage_2=="pre_blast_crisis_CML1266"]<-18
my.dat$type[my.dat$Stage_2=="blast_crisis" & my.dat$Patient=="CML1266"]<-15
my.dat$type[my.dat$Stage_2=="blast_crisis" & my.dat$Patient=="OX2046"]<-15
my.dat$type[my.dat$Stage_2=="blast_crisis_1_month_TKI"]<-15

x1<-which(my.dat$Stage_2=="pre_blast_crisis" & my.dat$Patient=="OX1931")
x2<-which(my.dat$my.color=="cyan")

################save data matrix to a file#####################################
write.table(my.dat[c("Cell","Dim1","Dim2","BCR_ABL","Patient",
						"Stage_2","my.color","type")],file="Fig6c.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
###############################################################################
pdf(file="Fig6c.pdf", width=7, height=7, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
plot(my.dat$Dim1,my.dat$Dim2,xaxt="n",yaxt="n",xlab="",ylab="",
		pch =my.dat$type,
		main="",col = my.dat$my.color,cex=1.5)

points(my.dat$Dim1[x1],my.dat$Dim2[x1],pch=23,bg="orange",cex=1.5,col="black")
points(my.dat$Dim1[x2],my.dat$Dim2[x2],pch=22,bg="cyan",col="black",cex=1.5)
#text(my.dat$Dim1[x1],my.dat$Dim2[x1],labels=my.dat$Cell[x1],col="black",cex=0.3)

legend("topright",c("Pre-BC OX1931","BC of OX1931")
				,pch=c(18,15),col=c("orange","cyan"),cex=1)

dev.off()