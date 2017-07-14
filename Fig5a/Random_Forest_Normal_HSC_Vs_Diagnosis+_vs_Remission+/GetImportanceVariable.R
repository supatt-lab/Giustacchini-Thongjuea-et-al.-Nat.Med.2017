###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used for the Random forests analysis.
###############################################################################
###############################################################################
library(randomForest)
load(file="RandomForrest_Normal_Diagnosis_Remission_2000tree.rdata")
my.importance<-importance(rf2000_normal_diagnosis_remission,type=2)
my.importance.f<-data.frame(my.importance)
my.importance.f$gene<-rownames(my.importance.f)
my.importance.f<-my.importance.f[order(my.importance.f$MeanDecreaseGini,decreasing=T),]
my.importance.f$rank<-1:nrow(my.importance.f)

pdf(file="RandomForestsImportanceScore.pdf", width=8, height=8, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)

plot(my.importance.f$rank,my.importance.f$MeanDecreaseGini,pch=19,col="red",
		ylab="Importance Score",xlab="Rank",cex=0.5,main="Importance variable genes")
text(my.importance.f$rank+40,my.importance.f$MeanDecreaseGini,labels=my.importance.f$gene,cex=0.5)

dev.off()

write.table(my.importance.f,file="RF_important_genes.txt",append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
#########################################################
