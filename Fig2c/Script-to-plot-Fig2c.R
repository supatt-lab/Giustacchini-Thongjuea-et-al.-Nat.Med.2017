#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig2c.
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
my.bulk<-read.table("Human_Bulk100Cells_rpkmforgenes_COUNT-UNIQUE.txt",header=T)
my.ox1407<-my.bulk["OX1407_bulk_HSC"]

my.cells$Gene<-rownames(my.cells)
my.ox1407$Gene<-rownames(my.ox1407)

my.combined<-merge(my.cells,my.ox1407,all=T)
my.combined[is.na(my.combined)==TRUE]<-0
rownames(my.combined)<-my.combined$Gene
my.combined<-my.combined[,-c(1)]
my.combined$ensembl_sc<-rowSums(my.combined[,(1:40)])

my.combined$ensembl_sc_rpkm<-1000000*(my.combined$ensembl_sc/sum(my.combined$ensembl_sc))
detected.genes<-subset(my.combined,ensembl_sc_rpkm >=1)

################################################################################
################save data matrix to a file####################
write.table(detected.genes,file="Fig2c.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)
################################################################################
my.siggenes<-read.table("DE-Genes-BCR-ABL+_Vs_BCR-ABL-_OX1407.txt",header=T)
n.siggenes<-subset(my.siggenes,fisher < 0.05 & abs(log2fc) >=3)
################################################################################
HouseKeeping<-read.table("OX1407.HouseKeepingGenes.txt",header=T)
################################################################################
my.DE.index<-which(rownames(my.combined) %in% c("CLU","FCER1A","TESPA1","GAS2","GOLGA8A","IFITM1","SOD2","CKLF","SAT1","MZB1","RGS2","CXCR4"))
my.H.index<-which(rownames(my.combined) %in% HouseKeeping$Gene)
gene.names<-c("CLU","FCER1A","TESPA1","GAS2","GOLGA8A","IFITM1","SOD2","CKLF","SAT1","MZB1","RGS2","CXCR4")

pdf(file="Fig2c.pdf", width=6, height=6, onefile=T, 
		bg="transparent",fonts = NULL,useDingbats=FALSE)

	x<-my.combined$OX1407_bulk_HSC
	y<-my.combined$ensembl_sc
	xlab="Bulk of 100 cells"
	ylab="Ensemble of single cells"
	col="blue"
	
	plot( NULL, xlim=c( -.1, 6.2 ), ylim=c( -1, 6.2 ),
			xaxt="n", yaxt="n", xaxs="i", yaxs="i", asp=1,
			xlab=xlab, ylab=ylab )
	
	abline( h=c(0,2,4,6), v=c(0,2,4,6), col = "lightgray", lwd=2 )
	
	points(
			ifelse( x > 0, log10(x), -.7 ),
			ifelse( y > 0, log10(y), -.7 ),
			pch=19, cex=.5, col = col )
	
	points(
			ifelse( x[my.H.index] > 0, log10(x[my.H.index]), -.7 ),
			ifelse( y[my.H.index] > 0, log10(y[my.H.index]), -.7 ),
			pch=19, cex=.5, col = "orange" )
	
	points(
			ifelse( x[my.DE.index] > 0, log10(x[my.DE.index]), -.7 ),
			ifelse( y[my.DE.index] > 0, log10(y[my.DE.index]), -.7 ),
			pch=19, cex=.7, col = "brown" )
	
	text(
			ifelse( x[my.DE.index] > 0, log10(x[my.DE.index]), -.7 ),
			ifelse( y[my.DE.index] > 0, log10(y[my.DE.index]), -.7 ),
			labels=gene.names, cex=.75, col = "red" )
	
	x<-log10(x+1)
	y<-log10(y+1)
	
	#x<-x[my.all.index]
	#y<-y[my.all.index]
	mylm<-lm(y~x)
	abline(mylm,col="lightgray")
	
	#newx<-seq(0,100)
	#prd<-predict(mylm,newdata=data.frame(x=newx),interval = c("confidence"))
	#lines(newx,prd[,2],col="red",lty=2)
	#lines(newx,prd[,3],col="red",lty=2)

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

dev.off()
