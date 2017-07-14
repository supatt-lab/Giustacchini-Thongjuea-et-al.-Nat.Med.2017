###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig6f-1.
###############################################################################
###############################################################################
data.a<-read.table("Relative-distance-between-RUNX1-binding-sites-and-DE-genes-of-OX1931_8_cells_vs_the_rest.txt",header=T)
data.a$distance2tss<-data.a$Distance_From_Nearest_TSS

pdf(file="Fig6f.1.pdf", width=6, height=6, onefile=T, bg="transparent",family = "Helvetica",fonts = NULL)

hist(abs(data.a$distance2tss),breaks=1000,xlim=c(0,50000),col="cyan",main="RUNX1 binding sites distribution",
		xlab="Distance relative to the nearest TSS (bp)",ylab="Frequency of differentially expressed genes")

dev.off()
