###############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script is used to generate Fig3d.
###############################################################################
library(limma)
library(pheatmap)
library(RColorBrewer)
###############################################################################
my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Stage_2=="diagnosis" | Stage_2=="normal_hsc")
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
###############################################################################
	S.m<-as.matrix(my.cells)
	min_expr<-1
	S.m[S.m < min_expr]<- 0
	S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
################################################################################
	my_bcr_abl<-as.character(my.anno$BCR_ABL)
	my_bcr_abl[my_bcr_abl=="negative_low_gapdh"]<-"negative"
	
	my_bcr_abl<-as.factor(my_bcr_abl)
	my.design<-model.matrix(~my_bcr_abl)
	rownames(my.design)<-colnames(S.m)
	
	batch.x<-as.character(my.anno$process_date)
	batch.x[which(batch.x=="batch1" | batch.x=="batch2")]<-"A"
	batch.x[which(batch.x=="batch5")]<-"B"
	batch.x[which(batch.x=="batch6")]<-"C"
	
	S.m.no.batch = removeBatchEffect(S.m,batch=batch.x,design=my.design)
	
################################################################################
	S.m<-S.m.no.batch
##########################Get all differentially expressed genes################
	dat1<-read.table("DEGenes-Normal-HSC_Vs_Diagnosis-.txt",header=T)
	dat2<-read.table("DEGenes-Normal-HSC_Vs_Diagnosis+.txt",header=T)
	dat3<-read.table("DEGenes-Diagnosis-_Vs_Diagnosis+.txt",header=T)
###################################################################################
	my.top=100
	
	dat1.sig<-dat1[1:my.top,]
	dat2.sig<-dat2[1:my.top,]
	dat3.sig<-dat3[1:my.top,]
	###############################################################################
	###############################################################################
	dat1.genes<-as.character(dat1.sig$gene)
	dat2.genes<-as.character(dat2.sig$gene)
	dat3.genes<-as.character(dat3.sig$gene)
	###############################################################################
	my.genes<-unique(c(dat1.genes,dat2.genes,dat3.genes))
	S.m2<-subset(S.m,rownames(S.m) %in% my.genes)
	################save data matrix to a file####################
	write.table(S.m2,file="Fig3d.data.matrix.txt",
			append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)
	###############################################################################
	my.anno.x<-as.character(my.anno$BCR_ABL)
	my.anno.x[my.anno.x=="negative_low_gapdh"]<-"negative"
			
	my.anno.p<-as.character(my.anno$Patient)
	my.anno.p[my.anno.p=="NBM_89"]<-"Normal_donors"
	my.anno.p[my.anno.p=="NBM_91"]<-"Normal_donors"
	my.anno.p[my.anno.p=="NBM_lonza_age_22"]<-"Normal_donors"
	my.anno.p[my.anno.p=="NBM_NC"]<-"Normal_donors"
	my.anno.p[my.anno.p=="NBM_petter"]<-"Normal_donors"
	my.anno.p[my.anno.p=="NBM_ERCC"]<-"Normal_donors"
	
	annotation_col<-data.frame(Cell=my.anno.x)
	rownames(annotation_col)<-my.anno$Cell
	annotation_col$Patient<-my.anno.p
	
	###############################################################################
	ann_colors = list(
			Cell = c(normal = "black", negative = "blue",positive="red"),
			Patient = c(Normal_donors="black",
					CML15="darkslategray4",
					CML22="tan",
					CML25="darkolivegreen3",
					CML655="lightcyan",
					CML656="gray",
					CML691="navajowhite4",
					CML940="lightsalmon",
					CML960="brown",
					OX1249="red",
					OX1302="green",
					OX1674="yellow",
				    OX2038="cyan",
					OX2286="slateblue1",
					OX664="blue",
					OX710="peachpuff1",
					OX714="orange",
					OX824="purple",
					OX967="darkgreen")
	)
	##################Heatmap######################################################
	pdf(file="Fig3d.pdf", width=7, height=7,onefile=F, bg="transparent")
	
		pheatmap(cor(S.m2,method="pearson"), clustering_distance_rows ="correlation",
		clustering_distance_cols ="correlation",annotation_col=annotation_col,
		fontsize_row=1,fontsize_col=1,color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
		annotation_colors = ann_colors,show_rownames=F,show_colnames=F, fontsize=6)
	dev.off()