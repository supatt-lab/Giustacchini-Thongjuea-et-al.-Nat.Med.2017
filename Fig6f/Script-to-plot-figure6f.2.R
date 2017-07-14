###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig6f-2.
###############################################################################
bg<-read.table("Bg.fraction.random2000.txt",header=T)
sig<-read.table("Siggenes.fraction.txt",header=T)

my.pvalue<-wilcox.test(sig$siggenes.fraction,bg$bg.fraction)$p.value
################save data matrix to a file#####################################
write.table(my.pvalue,file="Fig6f-2.wilcoxon.pvalue.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=FALSE)
###############################################################################
pdf(file="Fig6f.2.pdf", width=4, height=5, onefile=T, bg="transparent",family = "Helvetica",fonts = NULL)

boxplot(sig$siggenes.fraction,bg$bg.fraction,col=c("brown","gray"),ylab="Fraction of RUNX1 binding sites/windows")

dev.off()