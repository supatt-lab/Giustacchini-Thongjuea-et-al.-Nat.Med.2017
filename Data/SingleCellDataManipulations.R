#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017

exportCtValuesToExcel<-function(obj,outfile){
	stopifnot( is(obj, "SingleCellFluidigm" ))
	org.data<-CtValues(obj)
	write.csv(org.data, file =paste(outfile,".csv",sep=""))
}

showOutliers<-function(obj,outfile){
	stopifnot( is(obj, "SingleCellFluidigm" ))
	org.data<-OutlierCtValues(obj)
	input_chips<-SingleCellChipNames(obj)
	org.data$cell.name<-paste(org.data$cell,org.data$chip,org.data$group,sep="_")
	
	if(nrow(org.data) > 0){
		m.m <- dcast(org.data, gene~cell.name,value.var="ct",na.rm=T)
		x.m<-  as.matrix(m.m[,2:ncol(m.m)])
		rownames(x.m)<-m.m$gene
	}
	write.table(x.m,file=outfile,append=FALSE, sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)	
}

getTableOfCtValues<-function (obj){
	
	stopifnot( is(obj, "SingleCellFluidigm" ))
	data_dir<-SingleCellDir(obj)
	input_files<-biomarkFiles(obj)
	limit_ct<-LOD_PCR(obj)
	input_chips<-SingleCellChipNames(obj)
	
	objName <- deparse(substitute(obj))
	
	org.table<-data.frame()
	outliers.table<-data.frame()
			
			for(i in 1:length(input_files)){
				
				express<-readinCtSingleExp(paste(data_dir,input_files[i],sep=""),limit_ct)
				####################################################
				no.outliers<-removeOutliers(express,limit_ct)
				outliers<-getOutliers(express,limit_ct)
				
				my.dat<-data.frame(gene = rep(rownames(no.outliers), ncol(no.outliers)),
						cell = rep(colnames(no.outliers), each = nrow(no.outliers)),
						raw_ct=as.vector(no.outliers),log2ex=limit_ct-as.vector(no.outliers),chip=input_chips[i])
				my.dat$log2ex[is.na(my.dat$log2ex)==T]<-0
				org.table<-rbind(org.table,my.dat)
				
				if(!is.null(outliers)){
					outlier.dat<-data.frame(gene = rep(rownames(outliers), ncol(outliers)),
						cell = rep(colnames(outliers), each = nrow(outliers)),
						ct = as.vector(outliers),chip=input_chips[i])
						outliers.table<-rbind(outliers.table,outlier.dat)
				}
			}
			n.cell<-quantile(table(org.table$gene),0.75)
			x.check<-which(table(org.table$gene) < n.cell)
			if(length(x.check) > 0){
				print(names(x.check))
				stop("Found incompatible format of gene names (see above names)!!, please use the same gene name in the BioMark files!")
			}else{
				CtValues(obj)<-org.table
				OutlierCtValues(obj)<-outliers.table
				assign(objName,obj,envir=parent.frame())
				invisible(1)
				print("The table of Ct values is created.")
			}
}
getTableOfCtValuesReplaceNA<-function (obj){
	
	stopifnot( is(obj, "SingleCellFluidigm" ))
	data_dir<-SingleCellDir(obj)
	input_files<-biomarkFiles(obj)
	limit_ct<-LOD_PCR(obj)
	input_chips<-SingleCellChipNames(obj)
	
	objName <- deparse(substitute(obj))
	
	org.table<-data.frame()
	outliers.table<-data.frame()
	
	for(i in 1:length(input_files)){
		
		express<-readinCtSingleExp(paste(data_dir,input_files[i],sep=""),limit_ct)
		####################################################
		no.outliers<-removeOutliers(express,limit_ct)
		outliers<-getOutliers(express,limit_ct)
		
		my.dat<-data.frame(gene = rep(rownames(no.outliers), ncol(no.outliers)),
				cell = rep(colnames(no.outliers), each = nrow(no.outliers)),
				ct = as.vector(no.outliers),chip=input_chips[i])
		my.dat$ct[is.na(my.dat$ct)==T]<-limit_ct
		org.table<-rbind(org.table,my.dat)
		
		if(!is.null(outliers)){
			outlier.dat<-data.frame(gene = rep(rownames(outliers), ncol(outliers)),
					cell = rep(colnames(outliers), each = nrow(outliers)),
					ct = as.vector(outliers),chip=input_chips[i])
			outliers.table<-rbind(outliers.table,outlier.dat)
		}
	}
	n.cell<-quantile(table(org.table$gene),0.75)
	x.check<-which(table(org.table$gene) < n.cell)
	if(length(x.check) > 0){
		print(names(x.check))
		stop("Found incompatible format of gene names (see above names)!!, please use the same gene name in the BioMark files!")
	}else{
		CtValues(obj)<-org.table
		OutlierCtValues(obj)<-outliers.table
		assign(objName,obj,envir=parent.frame())
		invisible(1)
		print("The table of Ct values is created.")
	}
}

defineGroupOfCells<-function (obj,groupByChip=NULL,groupByAnnotationFile=NULL){
	
	stopifnot( is(obj, "SingleCellFluidigm" ))
	objName <- deparse(substitute(obj))
	
	if(is.null(groupByChip) ==T & is.null(groupByAnnotationFile)==T){
		stop("Please input the data frame or the file name of cell identification!")
	}
	if(is.null(groupByChip) ==F & is.null(groupByAnnotationFile)==F){
		stop("Please choose one way to define group of cells!")
	}
	
	org.data<-CtValues(obj)
	
	if(is.null(groupByChip)==F){
		if(class(groupByChip) != "data.frame"){
			stop("Require data frame for the input!")
		}
		if(nrow(groupByChip) >0){
			if(colnames(groupByChip)[1]=="chip" & colnames(groupByChip)[2]=="group"){
				org.data<-merge(org.data,groupByChip)
			}else{
				stop("Column names are incorrect!")
			}
		}
	}
	if(is.null(groupByAnnotationFile)==F){
		groupByAnnotation<-read.table(groupByAnnotationFile,header=T)
		if(nrow(groupByAnnotation) > 0){
			if(colnames(groupByAnnotation)[1]=="chip" & colnames(groupByAnnotation)[2]=="cell" 
					& colnames(groupByAnnotation)[3]=="group"){
				org.data<-merge(org.data,groupByAnnotation)
			}else{
				stop("Column names are incorrect!")
			}
		}	
	}
	
	CtValues(obj)<-org.data
	assign(objName,obj,envir=parent.frame())
	invisible(1)
	print("The table of Ct values has been updated with group names!.")
}

normalizeCtValuesByHouseKeepingGenes<-function(obj,houseKeeping.genes){
	stopifnot( is(obj, "SingleCellFluidigm" ))
	objName <- deparse(substitute(obj))
	org.data<-CtValues(obj)
	org.data$cell.name<-paste(org.data$cell,org.data$chip,sep="_")
	cell.names.uniq<-unique(org.data$cell.name)
	
	new.org.data<-list()
	
	if(length(houseKeeping.genes)==1){
		for(i in 1:length(cell.names.uniq)){
			my.exp<-subset(org.data,cell.name==cell.names.uniq[i])	
			hit.index<-which(my.exp$gene %in% houseKeeping.genes)
			if(length(hit.index) > 0){
				my.control<-subset(my.exp,gene==houseKeeping.genes[1])
				my.control.log2ex<-my.control$log2ex
				my.exp$delta.ct<-my.exp$log2ex-my.control.log2ex
				#my.exp$fold<-2^((my.exp$delta.ct)*-1)
				new.org.data<-rbindlist(list(new.org.data,my.exp))
			}else{
				stop("Could not find the house keeping genes in the list!!")
			}
		}
	}else if(length(houseKeeping.genes) > 1){
		for(i in 1:length(cell.names.uniq)){
			my.exp<-subset(org.data,cell.name==cell.names.uniq[i])	
			hit.index<-which(my.exp$gene %in% houseKeeping.genes)
			
			if(length(hit.index) != length(houseKeeping.genes)){
				stop("Please check your house keeping gene names!!")
			}
			if(length(hit.index) > 0){
				my.control<-my.exp[hit.index,]
				my.control.log2ex<-mean(my.control$log2ex)
				my.exp$delta.ct<-my.exp$log2ex-my.control.log2ex
				#my.exp$fold<-2^((my.exp$delta.ct)*-1)
				new.org.data<-rbindlist(list(new.org.data,my.exp))
			}else{
				stop("Could not find the house keeping genes in the list!!")
			}
		}
	}else{
		stop("Please input your house keeping genes!!")
	}
	CtValues(obj)<-as.data.frame(new.org.data)
	assign(objName,obj,envir=parent.frame())
	invisible(1)
	print("Done!!. The table of Ct values has been updated.")
}

normalizeCtValuesByHouseKeepingGenesWithRawCt<-function(obj,houseKeeping.genes){
	stopifnot( is(obj, "SingleCellFluidigm" ))
	objName <- deparse(substitute(obj))
	org.data<-CtValues(obj)
	org.data$cell.name<-paste(org.data$cell,org.data$chip,sep="_")
	cell.names.uniq<-unique(org.data$cell.name)
	
	new.org.data<-list()
	
	if(length(houseKeeping.genes)==1){
		for(i in 1:length(cell.names.uniq)){
			my.exp<-subset(org.data,cell.name==cell.names.uniq[i])	
			hit.index<-which(my.exp$gene %in% houseKeeping.genes)
			if(length(hit.index) > 0){
				my.control<-subset(my.exp,gene==houseKeeping.genes[1])
				my.control.raw_ct<-my.control$raw_ct
				my.exp$relative_ct<-my.exp$raw_ct-my.control.raw_ct
				#my.exp$fold<-2^((my.exp$delta.ct)*-1)
				new.org.data<-rbindlist(list(new.org.data,my.exp))
			}else{
				stop("Could not find the house keeping genes in the list!!")
			}
		}
	}else if(length(houseKeeping.genes) > 1){
		for(i in 1:length(cell.names.uniq)){
			my.exp<-subset(org.data,cell.name==cell.names.uniq[i])	
			hit.index<-which(my.exp$gene %in% houseKeeping.genes)
			
			if(length(hit.index) != length(houseKeeping.genes)){
				stop("Please check your house keeping gene names!!")
			}
			if(length(hit.index) > 0){
				my.control<-my.exp[hit.index,]
				my.control.raw_ct<-mean(my.control$raw_ct)
				my.exp$relative_ct<-my.exp$raw_ct-my.control.raw_ct
				#my.exp$fold<-2^((my.exp$delta.ct)*-1)
				new.org.data<-rbindlist(list(new.org.data,my.exp))
			}else{
				stop("Could not find the house keeping genes in the list!!")
			}
		}
	}else{
		stop("Please input your house keeping genes!!")
	}
	CtValues(obj)<-as.data.frame(new.org.data)
	assign(objName,obj,envir=parent.frame())
	invisible(1)
	print("Done!!. The table of Ct values has been updated.")
}

normalizeCtValuesByMedian<-function(obj){
	stopifnot( is(obj, "SingleCellFluidigm" ))
	objName <- deparse(substitute(obj))
	org.data<-CtValues(obj)
	org.data$cell.name<-paste(org.data$cell,org.data$chip,sep="_")
	cell.names.uniq<-unique(org.data$cell.name)
	
	new.org.data<-list()
	
	for(i in 1:length(cell.names.uniq)){
		
		my.exp<-subset(org.data,cell.name==cell.names.uniq[i])
		my.control.ct<-median(my.exp$raw_ct,na.rm=T)
		my.exp$delta.ct<-my.exp$raw_ct-my.control.ct
		new.org.data<-rbindlist(list(new.org.data,my.exp))
	}
	CtValues(obj)<-as.data.frame(new.org.data)
	assign(objName,obj,envir=parent.frame())
	invisible(1)
	print("Done!!. The table of Ct values has been updated.")
}
normalizeCtValuesByMedianWithLog2ex<-function(obj){
	stopifnot( is(obj, "SingleCellFluidigm" ))
	objName <- deparse(substitute(obj))
	org.data<-CtValues(obj)
	org.data$cell.name<-paste(org.data$cell,org.data$chip,sep="_")
	cell.names.uniq<-unique(org.data$cell.name)
	
	new.org.data<-list()
	
	for(i in 1:length(cell.names.uniq)){
		
		my.exp<-subset(org.data,cell.name==cell.names.uniq[i])
		my.control.ct<-median(my.exp$log2ex,na.rm=T)
		my.exp$delta.ct<-my.exp$log2ex-my.control.ct
		new.org.data<-rbindlist(list(new.org.data,my.exp))
	}
	CtValues(obj)<-as.data.frame(new.org.data)
	assign(objName,obj,envir=parent.frame())
	invisible(1)
	print("Done!!. The table of Ct values has been updated.")
}
getPA_matrix<- function(exp.m,min_ct,max_ct) {
	exp<-exp.m
	m.m<-matrix(0, nrow=nrow(exp.m), ncol=ncol(exp.m))
	for (i in seq(nrow(exp.m))){
		for (j in seq(ncol(exp.m))){
			n<-exp.m[i,j]
			if(is.na(n)==T){
				m.m[i,j]<-0
			}else if(n >= min_ct & n <=max_ct){
				m.m[i,j]<-1
			}
		}
	}
	rownames(m.m)<-rownames(exp)
	colnames(m.m)<-colnames(exp)
	return(m.m)
}

removeOutliers<-function(exp,limit_ct){
	
	exp.m<-as.matrix(exp$org_data)
	exp.m[exp.m >=limit_ct]<-NA
	
	#########check min and max ct values###################
	min_ct<-min(exp.m[,1:ncol(exp.m)], na.rm = TRUE)
	max_ct<-max(exp.m[,1:ncol(exp.m)], na.rm = TRUE)
	
	##########get present/absent matrix####################
	m.status<-getPA_matrix(exp.m,min_ct,max_ct)
	g.sum<-colSums(m.status,na.rm=T)
	###########modified z-score############################
	z.checked<-0.6745*(g.sum-median(g.sum))/(mad(g.sum))
	#no.outliers.idx<-which(z.checked > -3 & g.sum > quantile(g.sum,probs=0.15))
	no.outliers.idx<-which(z.checked > -3)
	
	no.outliers<-exp.m[,no.outliers.idx]
	return(no.outliers)
}

getOutliers<-function(exp,limit_ct){
	
	exp.m<-as.matrix(exp$org_data)
	exp.m[exp.m >=limit_ct]<-NA
	
	#########check min and max ct values###################
	min_ct<-min(exp.m[,1:ncol(exp.m)], na.rm = TRUE)
	max_ct<-max(exp.m[,1:ncol(exp.m)], na.rm = TRUE)
	
	##########get present/absent matrix####################
	m.status<-getPA_matrix(exp.m,min_ct,max_ct)
	g.sum<-colSums(m.status,na.rm=T)
	###########modified z-score############################
	z.checked<-0.6745*(g.sum-median(g.sum))/(mad(g.sum))
	#outliers.idx<-which(z.checked <= -3 & g.sum <= quantile(g.sum,probs=0.15))
	outliers.idx<-which(z.checked <= -3)
	if(length(outliers.idx) > 0){
		outliers<-exp.m[,outliers.idx]
		return(outliers)
	}
}

getGeneNames<-function(org.data){
	my.gene.names<-unique(org.data$gene)
	return(as.vector(my.gene.names))
}