#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to count number of detected genes per cell from the sampling of reads.
###############################################################################
source("../Data/RNA-Seq-S1.R")
my.cells<-read.table("Sampling_Normal_HSCs_hg19.RPKM-UNIQUE.txt",header=T)
###############################################################################
CountGenesPerCell<-getNumberOfDetectedGenesPerCell(my.cells,1)

write.table(CountGenesPerCell,file="HSC_Sampling_Count_Genes_Per_Cell_1_RPKM.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)
################################################################################

