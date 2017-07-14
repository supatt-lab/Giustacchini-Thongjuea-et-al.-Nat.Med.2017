#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig3e.####################
#############################################################
#######OX1249, OX1302, OX1674, and OX2038 are just the variables for my set up and
#######they DO NOT relate to any patients information here!.
OX1249<-read.delim("gsea_report_for_Normal_HSC_Vs_CML_Ensemble.txt",header=T,stringsAsFactors=F)
OX1302<-read.delim("gsea_report_for_Normal_HSC_Vs_Diagnosis-.txt",header=T,stringsAsFactors=F)
OX1674<-read.delim("gsea_report_for_Normal_HSC_Vs_Diagnosis+.txt",header=T,stringsAsFactors=F)
OX2038<-read.delim("gsea_report_for_Diagnosis-_VS_Diagnosis+.txt",header=T,stringsAsFactors=F)
###############################################################################
OX1249<-OX1249[c("NAME","FDR_q.val","Enriched_In")]
colnames(OX1249)<-c("NAME","FDR_OX1249","OX1249_Enriched_In")
OX1249$OX1249_Enriched_In[OX1249$OX1249_Enriched_In=="CML"]<-"EC"
OX1249$OX1249_Enriched_In[OX1249$OX1249_Enriched_In=="normal_hsc"]<-"EN"
OX1249$OX1249_Enriched_In<-paste(OX1249$OX1249_Enriched_In,format(OX1249$FDR_OX1249,digits=2),sep=" ")


OX1302<-OX1302[c("NAME","FDR_q.val","Enriched_In")]
colnames(OX1302)<-c("NAME","FDR_OX1302","OX1302_Enriched_In")
OX1302$OX1302_Enriched_In[OX1302$OX1302_Enriched_In=="BCR_ABL-"]<-"-"
OX1302$OX1302_Enriched_In[OX1302$OX1302_Enriched_In=="normal_hsc"]<-"N"
OX1302$OX1302_Enriched_In<-paste(OX1302$OX1302_Enriched_In,format(OX1302$FDR_OX1302,digits=2),sep=" ")


OX1674<-OX1674[c("NAME","FDR_q.val","Enriched_In")]
colnames(OX1674)<-c("NAME","FDR_OX1674","OX1674_Enriched_In")
OX1674$OX1674_Enriched_In[OX1674$OX1674_Enriched_In=="BCR_ABL+"]<-"+"
OX1674$OX1674_Enriched_In[OX1674$OX1674_Enriched_In=="normal_hsc"]<-"N"
OX1674$OX1674_Enriched_In<-paste(OX1674$OX1674_Enriched_In,format(OX1674$FDR_OX1674,digits=2),sep=" ")


OX2038<-OX2038[c("NAME","FDR_q.val","Enriched_In")]
colnames(OX2038)<-c("NAME","FDR_OX2038","OX2038_Enriched_In")
OX2038$OX2038_Enriched_In[OX2038$OX2038_Enriched_In=="BCR_ABL-"]<-"-"
OX2038$OX2038_Enriched_In[OX2038$OX2038_Enriched_In=="BCR_ABL+"]<-"+"
OX2038$OX2038_Enriched_In<-paste(OX2038$OX2038_Enriched_In,format(OX2038$FDR_OX2038,digits=2),sep=" ")

###############################################################################
###############################################################################
m<-merge(OX1249,OX1302,all=T)
m2<-merge(m,OX1674,all=T)
m3<-merge(m2,OX2038,all=T)
###############################################################################
info1<-subset(m3,FDR_OX1249 < 0.25 | FDR_OX1302 < 0.25 | FDR_OX1674 < 0.25 | FDR_OX2038 < 0.25)

r.names<-c("HALLMARK_UV_RESPONSE_DN","HALLMARK_UV_RESPONSE_UP","HALLMARK_APICAL_SURFACE","HALLMARK_ESTROGEN_RESPONSE_EARLY",
		"HALLMARK_ESTROGEN_RESPONSE_LATE","HALLMARK_MYOGENESIS","HALLMARK_PANCREAS_BETA_CELLS","HALLMARK_ANDROGEN_RESPONSE",
		"HALLMARK_SPERMATOGENESIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_HEME_METABOLISM")


r.index<-which(info1$NAME %in% r.names)
info1<-info1[-c(r.index),]

info.x<-info1[c("FDR_OX1249","FDR_OX1302","FDR_OX1674","FDR_OX2038")]
rownames(info.x)<-info1$NAME
colnames(info.x)<-c("Normal Vs CML","Normal Vs BCR-ABL-","Normal Vs BCR-ABL+","BCR-ABL- Vs BCR-ABL+")


info.y<-info1[c("OX1249_Enriched_In","OX1302_Enriched_In","OX1674_Enriched_In","OX2038_Enriched_In")]
rownames(info.y)<-info1$NAME
my.info<-as.matrix(info.y)


info2<- -log10(info.x+0.0001)
info2[info2 < 0]<-0
info2[is.na(info2)==T]<-0
#############Clustering##################################
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
write.table(m.m,file="Fig3e.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)

write.table(info.y,file="Fig3e.enrichment_status_and_pvalues_data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=FALSE)
#########################################################
library("gplots")
library(RColorBrewer)

pdf(file="Fig3e-v2-FINAL.pdf", width=10, height=8, onefile=T, bg="transparent",family = "Helvetica",fonts = NULL)

heatmap.2(m.m,Colv=F,Rowv=F,dendrogram="none", cellnote=info.y,col=colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
		scale="none",key=TRUE,symkey=FALSE,keysize = 0.8, density.info="none", trace="none",cexRow=0.9,cexCol=0.9,
		notecex=0.80,
		sepcolor="white",
		sepwidth=c(0.01,0.1),
		notecol="black",margins = c(10,25),main ="")

dev.off()
