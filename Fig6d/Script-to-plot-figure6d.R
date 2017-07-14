###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig6d.
###############################################################################
###############################################################################
library(RColorBrewer)
library(pheatmap)
###############################################################################
my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T,stringsAsFactors=F)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Stage_2=="pre_blast_crisis" | Stage_2=="blast_crisis"  
				| Stage_2=="blast_crisis_1_month_TKI" | Stage_2=="normal_hsc" | Stage_2=="diagnosis")
my.anno<-subset(my.anno,Patient !="OX2046cd19" & Patient !="OX1931cd19")
my.anno<-subset(my.anno,BCR_ABL=="positive" | Stage_2=="normal_hsc")
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
remove.index<-grep("MIR|SNORA|SNORD",rownames(my.cells))
my.cells<-my.cells[-c(remove.index),]
###############################################################################
min_expr<-1
S.m<-as.matrix(my.cells)
S.m[S.m < min_expr]<- 0
S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
###############################################################################
###############################################################################
l.genes<-c("BANK1","EBF1","FAIM3","LTB","BLNK","CD19","CD79A","DNTT",
			"IGJ","IRF8","LEF1","LY86","RAG1",
			"VPREB1","VPREB3","CD33","CEBPA","CSF1R","ID2",
			"IL3RA","IL6R","LY6E","MCL1","MPO","MSR1")
###############################################################################
my.dat1<-subset(S.m,rownames(S.m) %in% l.genes)

S.m<-data.frame()
for(i in 1:length(l.genes)){
	x.gene<-l.genes[i]
	x.dat<-subset(my.dat1,rownames(my.dat1)==x.gene)
	S.m<-rbind(S.m,x.dat)
}
###############################################################################
groupA.anno<-subset(my.anno,Stage_2=="blast_crisis")
groupA1.anno<-subset(groupA.anno,Patient=="OX2046")
groupA1.cells<-as.character(groupA1.anno$Cell)

groupA2.anno<-subset(groupA.anno,Patient=="CML1266")
groupA2.cells<-as.character(groupA2.anno$Cell)


my.cells<-c("OX1931D5","OX1931F7","OX1931D10","OX1931H9_1","OX1931A3","OX1931D10_1","OX1931E5_1","OX1931H7_1")
groupB.anno<-subset(my.anno,Patient=="OX1931")
r.index<-which(groupB.anno$Cell %in% my.cells)
groupB.anno<-groupB.anno[-c(r.index),]
groupB.cells<-as.character(groupB.anno$Cell)


groupC.cells<-c("OX1931D5","OX1931F7","OX1931D10","OX1931H9_1","OX1931A3","OX1931D10_1","OX1931E5_1","OX1931H7_1")


groupD.anno<-subset(my.anno,Stage_2=="normal_hsc")
groupD.cells<-as.character(groupD.anno$Cell)

groupE.anno<-subset(my.anno,Stage_2=="diagnosis")
groupE.cells<-as.character(groupE.anno$Cell)

groupF.anno<-subset(my.anno,Stage_2=="pre_blast_crisis" & Patient=="CML1266")
groupF.cells<-as.character(groupF.anno$Cell)


groupG.anno<-subset(my.anno,Stage_2=="blast_crisis_1_month_TKI")
groupG.cells<-as.character(groupG.anno$Cell)


groupA1.m<-S.m[,colnames(S.m) %in% groupA1.cells]
groupA2.m<-S.m[,colnames(S.m) %in% groupA2.cells]

groupB.m<-S.m[,colnames(S.m) %in% groupB.cells]
groupC.m<-S.m[,colnames(S.m) %in% groupC.cells]
groupD.m<-S.m[,colnames(S.m) %in% groupD.cells]
groupE.m<-S.m[,colnames(S.m) %in% groupE.cells]
groupF.m<-S.m[,colnames(S.m) %in% groupF.cells]
groupG.m<-S.m[,colnames(S.m) %in% groupG.cells]
##########################Get All Data##########################################
groupA1.m<-S.m[,colnames(S.m) %in% groupA1.cells]
groupA1.col.sum<-colSums(groupA1.m)
groupA1.m<-groupA1.m[,order(groupA1.col.sum)]

groupA2.m<-S.m[,colnames(S.m) %in% groupA2.cells]
groupA2.col.sum<-colSums(groupA2.m)
groupA2.m<-groupA2.m[,order(groupA2.col.sum)]

groupB.m<-S.m[,colnames(S.m) %in% groupB.cells]
groupB.col.sum<-colSums(groupB.m)
groupB.m<-groupB.m[,order(groupB.col.sum)]

groupC.m<-S.m[,colnames(S.m) %in% groupC.cells]
groupC.col.sum<-colSums(groupC.m)
groupC.m<-groupC.m[,order(groupC.col.sum)]

groupD.m<-S.m[,colnames(S.m) %in% groupD.cells]
groupD.col.sum<-colSums(groupD.m)
groupD.m<-groupD.m[,order(groupD.col.sum)]

groupE.m<-S.m[,colnames(S.m) %in% groupE.cells]
groupE.col.sum<-colSums(groupE.m)
groupE.m<-groupE.m[,order(groupE.col.sum)]

groupF.m<-S.m[,colnames(S.m) %in% groupF.cells]
groupF.col.sum<-colSums(groupF.m)
groupF.m<-groupF.m[,order(groupF.col.sum)]

groupG.m<-S.m[,colnames(S.m) %in% groupG.cells]
groupG.col.sum<-colSums(groupG.m)
groupG.m<-groupG.m[,order(groupG.col.sum)]

my.data.m<-cbind(groupE.m,groupF.m,groupC.m,groupB.m,groupG.m,groupA2.m,groupA1.m,groupD.m)

##################Heatmap######################################################
pdf(file="Fig6d.pdf", width=7, height=10, onefile=F, bg="transparent")

r.index<-which(my.anno$Cell %in% groupC.cells)
my.anno$Stage_2[r.index]<-"pre_BC_8_cells"

F.index<-which(my.anno$Cell %in% groupF.cells)
my.anno$Stage_2[F.index]<-"pre_BC_CML1266"

p1.index<-which(my.anno$Cell %in% groupA1.cells)
my.anno$Stage_2[p1.index]<-"BC_OX1931"

p2.index<-which(my.anno$Cell %in% groupA2.cells)
my.anno$Stage_2[p2.index]<-"BC_CML1266"

p3.index<-which(my.anno$Cell %in% groupG.cells)
my.anno$Stage_2[p3.index]<-"BC_CML1203"

my.anno.x<-as.character(my.anno$Stage_2)

annotation_col<-data.frame(Cell=as.character(my.anno.x))
rownames(annotation_col)<-my.anno$Cell

################save data matrix to a file#####################################
write.table(my.data.m,file="Fig6d.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)
###############################################################################
ann_colors = list(
		Cell = c(BC_CML1203= "deeppink",BC_OX1931= "cyan",
				BC_CML1266="purple",pre_blast_crisis="orange",
				normal_hsc="black",diagnosis="brown",
				pre_BC_8_cells="yellow",pre_BC_CML1266="darkgreen")
)

library("gplots")
library(RColorBrewer)
pheatmap(t(my.data.m),cluster_cols=F,cluster_rows=F,show_colnames=T,show_rownames=F, fontsize_col=10,key=TRUE,symkey=FALSE,keysize = 0.8,annotation_row=annotation_col,
		annotation_colors=ann_colors,col=colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),cexRow=0.75,cexCol=0.75)

dev.off()


