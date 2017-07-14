###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig4d.
###############################################################################
CML_A<-read.delim("GSEA-SUMMARY/gsea_report_for_good_vs_poor_diagnosis_BCR_ABL-.txt",header=T,stringsAsFactors=F)
CML_B<-read.delim("GSEA-SUMMARY/gsea_report_for_good_vs_poor_diagnosis_BCR_ABL+.txt",header=T,stringsAsFactors=F)
###############################################################################
CML_A<-CML_A[c("NAME","FDR_q.val","Enriched_In")]
colnames(CML_A)<-c("NAME","FDR_CML_A","CML_A_Enriched_In")
CML_A$CML_A_Enriched_In<-paste(CML_A$CML_A_Enriched_In,format(CML_A$FDR_CML_A,digits=2),sep=" ")


CML_B<-CML_B[c("NAME","FDR_q.val","Enriched_In")]
colnames(CML_B)<-c("NAME","FDR_CML_B","CML_B_Enriched_In")
CML_B$CML_B_Enriched_In<-paste(CML_B$CML_B_Enriched_In,format(CML_B$FDR_CML_B,digits=2),sep=" ")
###############################################################################
###############################################################################
m<-merge(CML_A,CML_B,all=T)
###############################################################################
info1<-subset(m,FDR_CML_A < 0.25 | FDR_CML_B < 0.25)
###############################################################################
r.names<-c("HALLMARK_APICAL_JUNCTION","HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_UV_RESPONSE_UP","HALLMARK_UV_RESPONSE_DN",
		"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_COAGULATION","HALLMARK_ESTROGEN_RESPONSE_LATE","HALLMARK_ESTROGEN_RESPONSE_EARLY",
		"HALLMARK_COMPLEMENT","HALLMARK_MYOGENESIS","HALLMARK_BILE_ACID_METABOLISM")


r.index<-which(info1$NAME %in% r.names)
info1<-info1[-c(r.index),]

info.x<-info1[c("FDR_CML_A","FDR_CML_B")]
rownames(info.x)<-info1$NAME
colnames(info.x)<-c("BCR-ABL- (diagnosis)","BCR-ABL+ (diagnosis)")

info.y<-info1[c("CML_A_Enriched_In","CML_B_Enriched_In")]
rownames(info.y)<-info1$NAME
my.info<-as.matrix(info.y)

info2<- -log10(info.x+0.00001)
info2[info2 < 0]<-0
info2[is.na(info2)==T]<-0
#############Clustering##################################
m.m<-as.matrix(info2)
m.m<-as.data.frame(m.m)

m.y<-as.matrix(info.y)
m.y<-as.data.frame(m.y)
#########################################################
my.pathways<-read.table("Pathways.txt",header=T)

my.data.m<-data.frame()
for(i in 1:nrow(my.pathways)){
	pathway<-as.character(my.pathways$Pathway[i])
	my.x<-m.m[rownames(m.m)==pathway,]
	my.data.m<-rbind(my.data.m,my.x)
}


my.data.y<-data.frame()
for(i in 1:nrow(my.pathways)){
	pathway<-as.character(my.pathways$Pathway[i])
	my.y<-m.y[rownames(m.y)==pathway,]
	my.data.y<-rbind(my.data.y,my.y)
}
#########################################################
#########################################################
m.m<-as.matrix(my.data.m)
info.y<-as.matrix(my.data.y)
################save data matrix to a file###############
write.table(m.m,file="Fig4e.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)

write.table(info.y,file="Fig4e.enrichment_status_and_pvalues_data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=FALSE)
#########################################################
library("gplots")
library(RColorBrewer)

pdf(file="Fig4e.pdf", width=8, height=8, onefile=T, bg="transparent",family = "Helvetica",fonts = NULL)

heatmap.2(m.m,Colv=F,Rowv=F, dendrogram="none",cellnote=info.y,col=colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
		scale="none",key=TRUE,symkey=FALSE,keysize = 1, density.info="none", trace="none",cexRow=0.80,cexCol=0.70,
		notecex=0.80,
		sepcolor="white",
		sepwidth=c(0.01,0.1),
		notecol="black",margins = c(10,25),main ="")
		
dev.off()