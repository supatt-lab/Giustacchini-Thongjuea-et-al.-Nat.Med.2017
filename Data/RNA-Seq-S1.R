#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################

getNumberOfDetectedCellsPerGene <- function(expression_matrix, min_expr){
	
	FM <- expression_matrix
	FM_genes <- do.call(rbind, apply(FM, 1, 
					function(x) {
						return(data.frame(
										num_detected_cells=sum(unlist(as.list(x)) >= min_expr)
								)
						)
					})
	)
	return(FM_genes)
}

getNumberOfDetectedGenesPerCell <- function(expression_matrix, min_expr){
	
	FM <- expression_matrix
	FM_cells <- do.call(rbind, apply(FM, 2, 
					function(x) {
						return(data.frame(
										num_genes_detected=sum(unlist(as.list(x)) >= min_expr)
								)
						)
					})
	)
	return(FM_cells)
}


getSingleCellOutliers_RNAseq<-function(expression_matrix, min_expr){

	CountGenesPerCell<-getNumberOfDetectedGenesPerCell(expression_matrix,min_expr)
	GenePerCell<-data.frame(SampleID=rownames(CountGenesPerCell),N.genes=CountGenesPerCell$num_genes_detected)
	GenePerCell<-GenePerCell[order(GenePerCell[,2],decreasing=T),]
	Half.sample.names<-GenePerCell[c(1:round(nrow(GenePerCell)/2)),]

	Half.cell.m<-expression_matrix[,colnames(expression_matrix) %in% Half.sample.names$SampleID]

	CountCellPerGene<-getNumberOfDetectedCellsPerGene(Half.cell.m,min_expr)
	common.genes<-subset(CountCellPerGene,num_detected_cells ==ncol(Half.cell.m))
	common.genes.m<-expression_matrix[rownames(expression_matrix) %in% rownames(common.genes),]

	CountGenesPerCell.common.genes<-getNumberOfDetectedGenesPerCell(common.genes.m,min_expr)
	###########modified z-score############################
	n.genes<-CountGenesPerCell.common.genes$num_genes_detected
	CountGenesPerCell.common.genes$z.N.Genes<-0.6745*(n.genes-median(n.genes))/(mad(n.genes))


	S.m<-as.matrix(common.genes.m)
	S.m[S.m < min_expr]<- 0
	S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])

	n.median<-colMedians(S.m)
	CountGenesPerCell.common.genes$z.median<-0.6745*(n.median-median(n.median))/(mad(n.median))

	outliers.idx<-which(CountGenesPerCell.common.genes$z.N.Genes < -3 & CountGenesPerCell.common.genes$z.median < -3)

	outliers.candidate<-rownames(CountGenesPerCell.common.genes)[outliers.idx]
	return(outliers.candidate)
}

PlotSingleCellOutliers_RNAseq<-function(expression_matrix, min_expr){
	
	CountGenesPerCell<-getNumberOfDetectedGenesPerCell(expression_matrix,min_expr)
	GenePerCell<-data.frame(SampleID=rownames(CountGenesPerCell),N.genes=CountGenesPerCell$num_genes_detected)
	GenePerCell<-GenePerCell[order(GenePerCell[,2],decreasing=T),]
	Half.sample.names<-GenePerCell[c(1:round(nrow(GenePerCell)/2)),]
	
	Half.cell.m<-expression_matrix[,colnames(expression_matrix) %in% Half.sample.names$SampleID]
	
	CountCellPerGene<-getNumberOfDetectedCellsPerGene(Half.cell.m,min_expr)
	common.genes<-subset(CountCellPerGene,num_detected_cells ==ncol(Half.cell.m))
	common.genes.m<-expression_matrix[rownames(expression_matrix) %in% rownames(common.genes),]
	
	CountGenesPerCell.common.genes<-getNumberOfDetectedGenesPerCell(common.genes.m,min_expr)
	
	S.m<-as.matrix(common.genes.m)
	S.m[S.m < min_expr]<- 0
	S.m[S.m >= min_expr]<-log2(S.m[S.m >= min_expr])
	
	#S.m2plot<-S.m[,order(colMedians(S.m),decreasing=T)]
	S.m2plot<-S.m
	boxplot(S.m2plot,col="red",ylab="Expression (Log2)",las=2,
			main=paste("Number of genes that are highly expressed in 50% of total cells: ",nrow(common.genes),sep=""),cex.axis=0.8,cex=0.4)
	
}


getGenesInOutliersDetection_RNAseq<-function(expression_matrix, min_expr){
	
	CountGenesPerCell<-getNumberOfDetectedGenesPerCell(expression_matrix,min_expr)
	GenePerCell<-data.frame(SampleID=rownames(CountGenesPerCell),N.genes=CountGenesPerCell$num_genes_detected)
	GenePerCell<-GenePerCell[order(GenePerCell[,2],decreasing=T),]
	Half.sample.names<-GenePerCell[c(1:round(nrow(GenePerCell)/2)),]
	Half.cell.m<-expression_matrix[,colnames(expression_matrix) %in% Half.sample.names$SampleID]
	
	CountCellPerGene<-getNumberOfDetectedCellsPerGene(Half.cell.m,min_expr)
	common.genes<-subset(CountCellPerGene,num_detected_cells ==ncol(Half.cell.m))
	return(rownames(common.genes))
}

get_N_violin_layout<-function (n) {
	as <- round(sqrt(n) + 0.5)
	bs <- round(sqrt(n) + 0.5)
	a23 <- round(sqrt(n * 3/2) + 0.5)
	b23 <- round(sqrt(n * 2/3) + 0.5)
	a34 <- round(sqrt(n * 4/3) + 0.5)
	b34 <- round(sqrt(n * 3/4) + 0.5)
	a35 <- round(sqrt(n * 5/3) + 0.5)
	b35 <- round(sqrt(n * 3/5) + 0.5)
	if (a23 * b23 == n) {
		a <- a23
		b <- b23
	}
	else if (a34 * b34 == n) {
		a <- a34
		b <- b34
	}
	else if (a35 * b35 == n) {
		a <- a35
		b <- b35
	}
	else if (as * bs == n) {
		a <- as
		b <- bs
	}
	else {
		if (a23 * b23 > n) {
			a <- a23
			b <- b23
		}
		else {
			a <- a23 + 1
			b <- b23
		}
	}
	my.layout <- c(a, b)
	return(my.layout)
}

ShowRNAseqViolins_Panel<-function(expression_matrix,cell_annotation,genes,groupColumnName,title){
	
	if(length(genes) > 100){
		stop("Too many genes to plot!, please select 10-100 genes")
	}
	#################Transform data to Log2##########
	S.m<-as.matrix(expression_matrix)
	#S.m[S.m < LOD]<-log2(S.m[S.m < LOD]+0.1)
	#S.m[S.m >= LOD]<-log2(S.m[S.m >=LOD])
	S.m<-log2(S.m+0.1)
	#################Read in cell annotation#########
	cell_names<-colnames(S.m)
	my_cell_anno<-data.frame()
	for(i in 1:length(cell_names)){
		F<-subset(cell_sheet,SampleID==cell_names[i])
		my_cell_anno<-rbind(my_cell_anno,F)
	}
	
	my.genes<-S.m[rownames(S.m) %in% genes,]
	my.group<-my_cell_anno[,colnames(my_cell_anno) %in% groupColumnName]
	group.unique<-unique(my.group)
	
	combined_data<-data.frame()
	
	for(i in 1:nrow(my.genes)){
		
		g.name<-rownames(my.genes)[i]
		g.exp<-subset(my.genes,rownames(my.genes)==g.name)
		combined.f<-data.frame(gene_id=as.vector(g.name),gene_exp=as.numeric(g.exp),sample_group=my.group)
		combined_data<-rbind(combined_data,combined.f)
	}
	
	font_size=0.8
	col=1:length(group.unique)
	y_scale <- c(-5,0,5,10,15,20,25,30,35,40)
	
	plot_layout<-get_N_violin_layout(nrow(my.genes))
	
	bwplot(gene_exp ~ sample_group | gene_id,as.table = TRUE, data=combined_data,
			par.strip.text = list(cex = font_size, lines = 1),
			par.settings = list(panel.background = list(col = "white")),
			panel = function(..., bg) {
				panel.abline(h = y_scale, col = "black")
				panel.superpose(...)
			}, panel.groups = panel.violin, groups = as.vector(sample_group),
			scales=list(x=list(rot=90, cex=1),y = list(alternating = 3, at = y_scale, labels = y_scale)),
			layout = plot_layout,col = col,ylab = "Expression (log2)",main=title)
	
}

ShowRNAseqViolins<-function(expression_matrix,cell_annotation,genes,groupColumnName,output_file){
	
	#################Transform data to Log2##########
	S.m<-as.matrix(expression_matrix)
	S.m[S.m < LOD]<-log2(S.m[S.m < LOD]+0.1)
	S.m[S.m >= LOD]<-log2(S.m[S.m >=LOD])
	
	#################Read in cell annotation#########
	cell_names<-colnames(S.m)
	my_cell_anno<-data.frame()
	for(i in 1:length(cell_names)){
		F<-subset(cell_sheet,SampleID==cell_names[i])
		my_cell_anno<-rbind(my_cell_anno,F)
	}
	
	my.genes<-S.m[rownames(S.m) %in% genes,]
	my.group<-my_cell_anno[,colnames(my_cell_anno) %in% groupColumnName]
	group.unique<-unique(my.group)
	
	pdf(file=output_file, width=12, height=8, onefile=T, bg="transparent")
	par(mfrow=c(2,2))
	
	for(i in 1:nrow(my.genes)){
		
		g.name<-rownames(my.genes)[i]
		g.exp<-subset(my.genes,rownames(my.genes)==g.name)
		
		my.list<-list()
		for(j in 1:length(group.unique)){
			j.index<-which(my.group==group.unique[j])
			my.list[[group.unique[j]]]<- g.exp[,c(j.index)]
		}
		names(my.list)<-group.unique
		
		x<-lapply(my.list, function (x) x[!is.na(x)])
		list.check<-listLen(x)
		index.r<-which(list.check < 2)
		
		if(length(index.r) >0 ){
			x<-x[-index.r]
		}else{
			x<-x
		}
		simple.violinplot(x,ylim=c(-10,20))
		title(main = g.name,ylab="expression (log2)")
	}
	dev.off()
}
