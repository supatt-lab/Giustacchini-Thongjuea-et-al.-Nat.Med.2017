###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig6a.
###############################################################################
###############################################################################
my.dat<-read.delim("tSNE_Diagnosis-PreBC-BC-Without-CML1600_Result.txt",header=T)


my.dat$my.color[my.dat$Stage_2=="normal_hsc"]<-"gray35"
my.dat$my.color[my.dat$Stage_2=="cell_line"]<-"sienna"
my.dat$my.color[my.dat$Stage_2=="diagnosis" | my.dat$Stage_2=="accelerated_phase"]<-"brown"
my.dat$my.color[my.dat$Stage_2=="blast_crisis" & my.dat$Patient=="CML1266"]<-"purple"
my.dat$my.color[my.dat$Stage_2=="blast_crisis" & my.dat$Patient=="OX2046"]<-"cyan"
my.dat$my.color[my.dat$Stage_2=="blast_crisis_1_month_TKI"]<-"deeppink"


my.dat$type[my.dat$Stage_2=="normal_hsc"]<-19
my.dat$type[my.dat$Stage_2=="cell_line"]<-19
my.dat$type[my.dat$Stage_2=="diagnosis"]<-17
my.dat$type[my.dat$Stage_2=="blast_crisis" & my.dat$Patient=="CML1266"]<-15
my.dat$type[my.dat$Stage_2=="blast_crisis" & my.dat$Patient=="OX2046"]<-15
my.dat$type[my.dat$Stage_2=="blast_crisis_1_month_TKI"]<-15

################save data matrix to a file####################################
write.table(my.dat[c("Cell","Dim1","Dim2","BCR_ABL","Patient",
						"Stage_2","my.color","type")],file="Fig6a.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)


pdf(file="Fig6a.pdf", width=7, height=7, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
plot(my.dat$Dim1,my.dat$Dim2,xaxt="n",yaxt="n",xlab="",ylab="",
		pch =my.dat$type,
		main="",col = my.dat$my.color,cex=1.5)
legend("topright",c( "Normal HSC",
					    "CP-CML",
						"BC-CML (CML1266)",
						"BC-CML (OX1931)",
						"BC-CML (CML1203)",
						"K562"),pch=c(19,17,15,15,15),col=c("gray35","brown","purple","cyan","deeppink","sienna"),cex=0.95)

#legend("topright",c("CP-CML","BC-CML"),pch=c(17,15),cex=1.2)
dev.off()
###################################################