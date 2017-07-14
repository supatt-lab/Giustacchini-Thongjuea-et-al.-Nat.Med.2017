#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig2d.
###############################################################################
library(matrixStats)

my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Patient=="OX1407")

###############################################################################
my.cells<-read.table("Human_HSC_CML_DATA1_rpkmforgenes_COUNT-UNIQUE.txt",header=T)
my.cells<-my.cells[,as.character(my.anno$Cell)]
###############################################################################
###############################################################################
groupA<-as.character(my.anno$Cell[grep("negative",my.anno$BCR_ABL)])
groupA.m<-my.cells[,colnames(my.cells) %in% groupA]

groupB<-as.character(my.anno$Cell[grep("positive",my.anno$BCR_ABL)])
groupB.m<-my.cells[,colnames(my.cells) %in% groupB]

groupA.m$ensemble_negative<-rowSums(groupA.m)
groupB.m$ensemble_positive<-rowSums(groupB.m)


group.combined<-data.frame(ensemble_negative=groupA.m$ensemble_negative,
						   ensemble_positive=groupB.m$ensemble_positive)
				   
rownames(group.combined)<-rownames(groupB.m)
group.combined$sum<-rowSums(group.combined)
group.combined<-subset(group.combined,sum > 0)

################################################################################
################save data matrix to a file####################
write.table(group.combined,file="Fig2d.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)
###############################################################################
my.siggenes<-read.table("DE-Genes-BCR-ABL+_Vs_BCR-ABL-_OX1407.txt",header=T)
n.siggenes<-subset(my.siggenes,fisher < 0.05 & abs(log2fc) >=3)
################################################################################
HouseKeeping<-read.table("OX1407.HouseKeepingGenes.txt",header=T)
################################################################################
my.DE.index<-which(rownames(group.combined) %in% c("CLU","FCER1A","TESPA1","GAS2","GOLGA8A","IFITM1","SOD2","CKLF","SAT1","MZB1",
				"RGS2","CXCR4"))
my.H.index<-which(rownames(group.combined) %in% HouseKeeping$Gene)
################################################################################
my.g.index<-which(rownames(group.combined) %in% c("CLU","FCER1A","TESPA1","GAS2","GOLGA8A","IFITM1","SOD2","CKLF","SAT1","MZB1",
				"RGS2","CXCR4"))

geneScatterplot <- function( x, y,gene.names, xlab, ylab, col ) {
	plot( NULL, xlim=c( -.1, 6.2 ), ylim=c( -1, 6.2 ),
			xaxt="n", yaxt="n", xaxs="i", yaxs="i", asp=1,
			xlab=xlab, ylab=ylab )
	#abline( a=-1, b=1, col = "lightgray", lwd=2 )
	
	#abline( a=1, b=1, col = "lightgray", lwd=2 )
	abline( h=c(0,2,4,6), v=c(0,2,4,6), col = "lightgray", lwd=2 )
	points(
			ifelse( x > 0, log10(x), -.7 ),
			ifelse( y > 0, log10(y), -.7 ),
			pch=19, cex=.6, col = col )
	
	points(
			ifelse( x[my.H.index] > 0, log10(x[my.H.index]), -.7 ),
			ifelse( y[my.H.index] > 0, log10(y[my.H.index]), -.7 ),
			pch=19, cex=.6, col = "orange" )
	
	points(
			ifelse( x[my.DE.index] > 0, log10(x[my.DE.index]), -.7 ),
			ifelse( y[my.DE.index] > 0, log10(y[my.DE.index]), -.7 ),
			pch=19, cex=.7, col = "brown" )
	

	text(
			ifelse( x[my.g.index] > 0, log10(x[my.g.index]), -.7 ),
			ifelse( y[my.g.index] > 0, log10(y[my.g.index]), -.7 ),
			labels=gene.names, cex=.75, col = "red" )
			
	
	abline( a=0, b=1, col = "lightgray", lwd=2 )
	
	axis( 1, c( -.7, 0:6 ),
			c( "0", "1", "10", "100", expression(10^3), expression(10^4),
					expression(10^5), expression(10^6) ) )
	axis( 2, c( -.7, 0:6 ),
			c( "0", "1", "10", "100", expression(10^3), expression(10^4),
					expression(10^5), expression(10^6) ), las=2 )
	axis( 1, -.35, "//", tick=FALSE, line=-.7 )
	axis( 2, -.35, "\\\\", tick=FALSE, line=-.7 )
	cor.txt<-cor(x,y)
	cor.txt<-format(cor.txt,digit=3)
	text(5,0,paste("Pearson correlation=",cor.txt,sep=""))
}

pdf(file="Fig2d.pdf", width=6, height=6, onefile=T, 
		bg="transparent",fonts = NULL,useDingbats=FALSE)
geneScatterplot(group.combined$ensemble_negative,group.combined$ensemble_positive,
		rownames(group.combined[my.g.index,]),
		"Ensemble of single BCR-ABL- cells", "Ensemble of single BCR-ABL+ cells","blue")
dev.off()