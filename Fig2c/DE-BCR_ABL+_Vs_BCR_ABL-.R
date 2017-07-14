##############################################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
#This R script for differentially expressed gene analysis
###############################################################################
###############################################################################
my.anno<-read.delim("../Data/CML_PROJECT_ALL_CELLs.Freezed-5.anno.txt",header=T)
my.anno<-subset(my.anno,used4analysis==1)
my.anno<-subset(my.anno,Patient=="OX1407")
###############################################################################
load(file="../Data/CML_PROJECT_ALL_CELLs.rdata")
my.cells<-CML_PROJECT_ALL_CELLs[,as.character(my.anno$Cell)]
remove.index<-grep("MIR|SNORA|SNORD",rownames(my.cells))
my.cells<-my.cells[-c(remove.index),]
###############################################################################
###############################################################################
min_expr<-1
S.m<-as.matrix(my.cells)
S.m[S.m < min_expr]<- 0
S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
###############################################################################
###############################################################################
groupA<-as.character(my.anno$Cell[grep("negative",my.anno$BCR_ABL)])
groupA.m<-S.m[,colnames(S.m) %in% groupA]

groupB<-as.character(my.anno$Cell[grep("positive",my.anno$BCR_ABL)])
groupB.m<-S.m[,colnames(S.m) %in% groupB]

###############################################################################
fishersMethod <-function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)

result.f<-data.frame()

for(i in 1:nrow(groupA.m)){
	
	n<-groupA.m[i,]
	m<-groupB.m[i,]
	n.f<-n[n>0]
	m.f<-m[m>0]
	
	x<-length(n)
	y<-length(m)
	
	x.1<-length(n.f)
	y.1<-length(m.f)
	x.2<-x-x.1
	y.2<-y-y.1
	
	m.m <- matrix(c(x.1, x.2, y.1, y.2), ncol = 2)
	fisher.test.p <-fisher.test(m.m)$p.value
	chi.test.p <-chisq.test(m.m)$p.value
	
	wilcox.p<-wilcox.test(n,m,paired = FALSE)$p.value
	
	all.p<-c(fisher.test.p,wilcox.p)
	fisher.p<-fishersMethod (all.p)
	
	my.meanA<-mean(groupA.m[i,])
	my.meanB<-mean(groupB.m[i,])
	my.log2.fc<-(my.meanB+1)-(my.meanA+1)
	my.test.f<-data.frame(gene=rownames(groupA.m)[i],meanA=my.meanA,meanB=my.meanB,
			log2fc=my.log2.fc,nCellA=x,nCellB=y,expCellA=x.1,expCellB=y.1,
			expFractionCellA=x.1/x,expFractionCellB=y.1/y,chisq.test.p=chi.test.p,fisher.test.p=fisher.test.p,wilcox=wilcox.p,fisher=fisher.p)		
	result.f<-rbind(result.f,my.test.f)
}

result.f<-result.f[order(result.f$fisher),]
result.f$p.adjust<-p.adjust(result.f$fisher, method ="BH")
#############################################################################
write.table(result.f,file="DE-Genes-BCR-ABL+_Vs_BCR-ABL-_OX1407.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)