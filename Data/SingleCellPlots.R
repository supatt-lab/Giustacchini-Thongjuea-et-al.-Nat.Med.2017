#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
###############################################################################
generateRawCtBoxplotToPdf<-function(obj,output_file){
	stopifnot( is(obj, "SingleCellFluidigm" ))
	org.data<-CtValues(obj)
	input_chips<-SingleCellChipNames(obj)
	org.data$cell.name<-paste(org.data$cell,org.data$chip,org.data$group,sep="_")
	
	pdf(file=output_file, width=12, height=8, onefile=T, bg="transparent")
	
	for(i in 1:length(input_chips)){
		my.chip<-subset(org.data,chip==input_chips[i])
		m.m <- dcast(my.chip, gene~cell.name,value.var="ct",na.rm=T)
		x.m<-  as.matrix(m.m[,2:ncol(m.m)])
		
		boxplot(x.m,col="lightblue",las=2,cex.axis=0.5,xlab="Cell",ylab="Ct Value",main=paste("chip:",input_chips[i],sep=""))
		
	}
	dev.off()
}

generateRawCtToPdf<-function(obj,output_file){
	stopifnot( is(obj, "SingleCellFluidigm" ))
	org.data<-CtValues(obj)
	gene_names<-as.vector(unique(org.data$gene))
	limit_ct<-LOD_PCR(obj)
	input_chips<-SingleCellChipNames(obj)	
	
	pdf(file=output_file, width=12, height=8, onefile=T, bg="transparent")
	
	par(mfrow=c(2,2))
	colors<-rgb(runif(length(input_chips)),runif(length(input_chips)),runif(length(input_chips)))
	
	for(i in 1:length(gene_names)){
		g.name<-gene_names[i]
		g.exp<-subset(org.data,gene==g.name)
		g.exp.col<-data.frame()
		
		info<-c()
		for(j in 1:length(input_chips)){
			my.g.exp<-subset(g.exp,chip==input_chips[j])
			my.g.exp$color<-colors[j]
			g.exp.col<-rbind(g.exp.col,my.g.exp)
			info<-append(info,nrow(my.g.exp))
		}
		
		g.exp.col$ct[is.na(g.exp.col$ct)==T]<-limit_ct
		plot(g.exp.col$ct,ylab="Ct",xlab="Cell",col=g.exp.col$color,pch=15,main=paste("Ct values of : ",g.name,sep=""),ylim=c(5,limit_ct))
		abline(v=cumsum(info),col="grey",lwd=1,lty=3)
		text(cumsum(info)-10,limit_ct,input_chips,cex=0.8)
	}
	dev.off()
}

generateDeltaCtToPdf<-function(obj,output_file){
	stopifnot( is(obj, "SingleCellFluidigm" ))
	org.data.norm<-CtValues(obj)
	gene_names<-as.vector(unique(org.data.norm$gene))
	limit_ct<-LOD_PCR(obj)
	input_chips<-SingleCellChipNames(obj)	

	pdf(file=output_file, width=12, height=8, onefile=T, bg="transparent")
	
	par(mfrow=c(2,2))
	colors<-rgb(runif(length(input_chips)),runif(length(input_chips)),runif(length(input_chips)))
	
	for(i in 1:length(gene_names)){
		g.name<-gene_names[i]
		g.exp<-subset(org.data.norm,gene==g.name)
		g.exp.col<-data.frame()
		
		info<-c()
		for(j in 1:length(input_chips)){
			my.g.exp<-subset(g.exp,chip==input_chips[j])
			my.g.exp$color<-colors[j]
			g.exp.col<-rbind(g.exp.col,my.g.exp)
			info<-append(info,nrow(my.g.exp))
		}
		
		g.exp.col$delta.ct[is.na(g.exp.col$delta.ct)==T]<- -10
		plot(g.exp.col$delta.ct,ylab="Delta Ct",xlab="Cell",col=g.exp.col$color,pch=15,main=paste("Delta Ct values of : ",g.name,sep="")
		,ylim=c(-10,20))
		abline(v=cumsum(info),col="grey",lwd=1,lty=3)
		text(cumsum(info)-10,max(g.exp.col$delta.ct,na.rm=T),input_chips,cex=0.8)

	}
	dev.off()
}

generateBeeswarmPlotToPdf<-function(obj,bubble.size=0.6,output_file){
	stopifnot( is(obj, "SingleCellFluidigm" ))
	org.data.norm<-CtValues(obj)
	gene_names<-as.vector(unique(org.data.norm$gene))
	limit_ct<-LOD_PCR(obj)
	input_chips<-SingleCellChipNames(obj)
	input_groups<-SingleCellGroupNames(obj)
	
	pdf(file=output_file, width=12, height=8, onefile=T, bg="transparent")
	
	par(mfrow=c(2,2))
	
	for(i in 1:length(gene_names)){
		
		g.name<-gene_names[i]
		g.exp<-subset(org.data.norm,gene==g.name)
		
		groups<-unique(input_groups)
		
		my.list<-list()
		for(j in 1:length(groups)){
			my.g.exp<-subset(g.exp,group==groups[j])
			my.list[[groups[j]]]<- log2(my.g.exp$fold)
		}
		if(sum(g.exp$delta.ct,na.rm=T) > 0){
			beeswarm(my.list,pch = 16,col = 1:length(groups),ylim=c(-20,2),
					las=2,cex=bubble.size,main=g.name,ylab="Relative expression (log2)")
		}else{
			print(paste("Could not plot",g.name,sep=":"))
		}
	}
	dev.off()
}

generateViolinPlotToPdf<-function(obj,output_file){
	
	stopifnot( is(obj, "SingleCellFluidigm" ))
	org.data.norm<-CtValues(obj)
	gene_names<-as.vector(unique(org.data.norm$gene))
	input_groups<-SingleCellGroupNames(obj)
	
	pdf(file=output_file, width=12, height=8, onefile=T, bg="transparent")
	
	par(mfrow=c(2,2))
	
		for(i in 1:length(gene_names)){
			
			g.name<-gene_names[i]
			g.exp<-subset(org.data.norm,gene==g.name)
		
			groups<-unique(input_groups)
		
			my.list<-list()
		
			for(j in 1:length(groups)){
				my.g.exp<-subset(g.exp,group==groups[j])
				my.list[[groups[j]]]<- log2(my.g.exp$fold)
			}
	
			x<-lapply(my.list, function (x) x[!is.na(x)])
			list.check<-listLen(x)
			index.r<-which(list.check < 2)
			if(length(index.r) >0 ){
				x<-x[-index.r]
			}else{
				x<-x
			}
			if(sum(list.check) > 5){
				simple.violinplot(x,ylim=c(-20,2))
				title(main = g.name,ylab="Relative expression (log2)")
			}else{
				print(paste("Couldn't plot:", g.name,sep=""))
			}
			
		}
		dev.off()
}

generateCoExpressedGeneFractionPlotToPdf<-function(obj,target.gene,gene.set,gene.set.label,output_file){
	
	stopifnot( is(obj, "SingleCellFluidigm" ))
	target.gene<-toupper(target.gene)
	gene.set<-toupper(gene.set)
	
	org.data.norm<-CtValues(obj)
	input_groups<-SingleCellGroupNames(obj)
	
	my.exp<-subset(org.data.norm,gene==target.gene)	
	my.exp<-my.exp[which(!is.na(my.exp$delta.ct)),]
	
	cell.names<-unique(my.exp$cell.name)
	my.selected.cells<-org.data.norm[org.data.norm$cell.name %in% cell.names,]
	my.selected.cells<-my.selected.cells[which(!is.na(my.selected.cells$delta.ct)),]
	
	my.selected.genes<-my.selected.cells[my.selected.cells$gene %in% gene.set,]
	
	group<-unique(input_groups)
	
	freq.t<-as.data.frame(table(as.character(my.selected.genes$gene),my.selected.genes$group))
	colnames(freq.t)<-c("gene","group","n.cell")
	
	group.f<-data.frame()
	group.a<-data.frame()
	
	for(i in 1:length(unique(input_groups))){
		s.group<-unique(input_groups)[i]
		all.cell<-subset(org.data.norm,group==s.group)
		all.cell.n<-length(unique(all.cell$cell.name))
		group.c<-data.frame(group=s.group,cell.n=all.cell.n)
		group.a<-rbind(group.a,group.c)
	}
	
	for(i in 1:length(unique(input_groups))){
		s.group<-unique(input_groups)[i]
		all.cell<-subset(my.selected.cells,group==s.group)
		all.cell.n<-length(unique(all.cell$cell.name))
		group.c<-data.frame(group=s.group,cell.n=all.cell.n)
		group.f<-rbind(group.f,group.c)
	}
	
	fraction.t<-data.frame()
	for(j in 1:nrow(group.f)){
		my.group<-group.f[j,]
		my.t<-subset(freq.t,group==my.group$group)
		my.t$fraction<-(my.t$n.cell)/(my.group$cell.n)
		fraction.t<-rbind(fraction.t,my.t)
	}
	
	my.f.list<-list()
	for(k in 1:nrow(group.f)){
		my.group<-group.f$group[k]
		my.fraction<-subset(fraction.t,group==my.group)
		my.fraction<-my.fraction[order(my.fraction$gene),]
		f<-data.frame(fraction=my.fraction$fraction)
		colnames(f)<-my.group
		rownames(f)<-my.fraction$gene
		my.f.list[[my.group]]<-f
	}
	
	my.final.fraction<-t(do.call(cbind,my.f.list))
	
	pdf(file=output_file, width=12, height=8, onefile=T, bg="transparent")
	par(mfrow=c(2,2))
	
	barplot(my.final.fraction,beside=T,ylim=c(0,1),
			border="white",col=1:length(rownames(my.final.fraction)),ylab="Fraction",
			xlab="Gene",main=paste(target.gene,"expressing cells",sep=" "),cex.name=0.75)
	
	slices <- group.f$cell.n
	slices.1 <- group.a$cell.n
	lbls <- group.f$group
	lbls.n<- paste(slices,slices.1,sep="/")
	lbls <- paste(lbls, lbls.n,sep=":") 
	
	pie(slices,labels = lbls, col=1:length(rownames(my.final.fraction)),main=paste("Number of ", target.gene," expressing cells",sep="")) 
	###############plot frequency###############
	
	u.group<-unique(input_groups)
	
	my.group.list<-list()
	for(i in 1:length(group)){
		my.group.cell<-subset(my.selected.cells,group==u.group[i])
		my.cell.names<-unique(my.group.cell$cell.name)
		
		g.freq<-c()
		for(j in 1:length(my.cell.names)){
			each.cell<-subset(my.group.cell,cell.name==my.cell.names[j])
			checked.each.cell<-each.cell[each.cell$gene %in% gene.set,]
			if(nrow(checked.each.cell) > 0){
				g.freq<-append(g.freq,nrow(checked.each.cell))
			}else{
				g.freq<-append(g.freq,0)
			}
		}
		my.group.list[[u.group[i]]]<-g.freq
	}
	
	co.g.n<-matrix(0,nrow=length(gene.set)+1,ncol=length(my.group.list))
	rownames(co.g.n)<-c(0:length(gene.set))
	colnames(co.g.n)<-names(my.group.list)
	
	for(i in 1:length(my.group.list)){
		y<-sapply(my.group.list[i], function(x) table(x))
		co.g.n[rownames(y),colnames(y)]<-y
	}
	co.g.n<-t(co.g.n)
	co.g.n.f<-co.g.n/group.f$cell.n
	barplot(co.g.n.f,beside=T,col=1:length(rownames(my.final.fraction)),ylab="Fraction",ylim=c(0:1),
			xlab="Number of genes co-expressed",main=paste("Gene set :",gene.set.label,sep=""))
	dev.off()
}

generateSingleCellHeatMapToPdf<-function(obj,removed.genes.list,fontsize=1,output_file){
	
	stopifnot( is(obj, "SingleCellFluidigm" ))
	removed.genes.list<-toupper(removed.genes.list)
	my.expr<-CtValues(obj)
	z.mean<-mean(my.expr$delta.ct,na.rm=T)
	z.sd<-sd(my.expr$delta.ct,na.rm=T)
	
	my.expr$z_score<-as.numeric((my.expr$delta.ct-z.mean)/z.sd)*-1
	my.expr$z_score[is.na(my.expr$z_score)==T]<-min(my.expr$z_score,na.rm=T)
	#######remove genes from the analysis####################
	if(is.character(removed.genes.list)){
		if(length(removed.genes.list) >0 ){
		removed.index<-which(my.expr$gene %in% removed.genes.list)
		my.expr<-my.expr[-c(removed.index),]
		}
	}	
	m.m <- dcast(my.expr, gene~cell.name,value.var="z_score",na.rm=T)
	x.m<-as.matrix(m.m[,2:ncol(m.m)])
	rownames(x.m)<-m.m$gene
	
	anno<-data.frame(cell=colnames(x.m))
	
	my.corr.group.v<-c()
	for(i in 1:nrow(anno)){
		gr<-subset(my.expr,cell.name==anno[i,])
		gr.id<-as.character(unique(gr$group))
		my.corr.group.v<-append(my.corr.group.v,gr.id)
	}
	anno.f<-data.frame(group=my.corr.group.v)
	rownames(anno.f)<-anno$cell
	
	pdf(file=output_file, width=12, height=8, onefile=T, bg="transparent")
	
	pheatmap(x.m,annotation = anno.f,fontsize = fontsize)
	
	dev.off()
}

generateSingleCellGeneCorrelationToPdf<-function(obj,removed.genes.list,fontsize,output_file){
	
	stopifnot( is(obj, "SingleCellFluidigm" ))
	removed.genes.list<-toupper(removed.genes.list)
	my.expr<-CtValues(obj)
	
	max.value<-ifelse(max(my.expr$delta.ct,na.rm=T) < 35, 35, max(my.expr$delta.ct,na.rm=T))
	my.expr$call_values<-my.expr$delta.ct
	my.expr$call_values[is.na(my.expr$call_values)==T]<-max.value
	####################remove genes from the list#######################
	if(is.character(removed.genes.list)){
		if(length(removed.genes.list) >0 ){
			removed.index<-which(my.expr$gene %in% removed.genes.list)
			my.expr<-my.expr[-c(removed.index),]
		}
	}	
	m.m <- dcast(my.expr, cell.name~gene,value.var="call_values",na.rm=T)
	x.m<-as.matrix(m.m[,2:ncol(m.m)])
	rownames(x.m)<-m.m$cell.name
	#####################################################################
	non.zero.var <- logical()
	for(i in 1:ncol(x.m)) {
		non.zero.var[i] <- var(x.m[,i]) > 0
	}
	x.m.cleaned <- x.m[,non.zero.var]
	cor.x.m<-cor(x.m.cleaned,method="pearson")
	#####################plot############################################
	pdf(file=output_file, width=12, height=8, onefile=T, bg="transparent")
	pheatmap(cor.x.m,fontsize = fontsize,display_numbers = TRUE)
	dev.off()
}

generateSingleCellCorrelationToPdf<-function(obj,removed.genes.list,fontsize,output_file){
	
	stopifnot( is(obj, "SingleCellFluidigm" ))
	removed.genes.list<-toupper(removed.genes.list)
	my.expr<-CtValues(obj)
	
	max.value<-ifelse(max(my.expr$delta.ct,na.rm=T) < 35, 35, max(my.expr$delta.ct,na.rm=T))
	my.expr$call_values<-my.expr$delta.ct
	my.expr$call_values[is.na(my.expr$call_values)==T]<-max.value
	####################remove genes from the list#######################
	if(is.character(removed.genes.list)){
		if(length(removed.genes.list) >0 ){
			removed.index<-which(my.expr$gene %in% removed.genes.list)
			my.expr<-my.expr[-c(removed.index),]
		}
	}	
	m.m <- dcast(my.expr, cell.name~gene,value.var="call_values",na.rm=T)
	x.m<-as.matrix(m.m[,2:ncol(m.m)])
	rownames(x.m)<-m.m$cell.name
	#####################################################################
	non.zero.var <- logical()
	for(i in 1:ncol(x.m)) {
		non.zero.var[i] <- var(x.m[,i]) > 0
	}
	x.m.cleaned <- x.m[,non.zero.var]
	x.m.cleaned.t<-t(x.m.cleaned)
	
	cor.x.m<-cor(x.m.cleaned.t,method="pearson")
	
	anno<-data.frame(cell=colnames(cor.x.m))
	
	my.corr.group.v<-c()
	for(i in 1:nrow(anno)){
		gr<-subset(my.expr,cell.name==anno[i,])
		gr.id<-as.character(unique(gr$group))
		my.corr.group.v<-append(my.corr.group.v,gr.id)
	}
	anno.f<-data.frame(group=my.corr.group.v)
	rownames(anno.f)<-anno$cell
	#####################plot############################################
	pdf(file=output_file, width=12, height=8, onefile=T, bg="transparent")
	pheatmap(cor.x.m,fontsize = fontsize,annotation = anno.f)
	dev.off()
}

generateSingleCellPCAtoPdf<-function(obj,removed.genes.list,output_file){
	
	stopifnot( is(obj, "SingleCellFluidigm" ))
	input_groups<-SingleCellGroupNames(obj)
	my.expr<-CtValues(obj)
	new.expr<-data.frame()
	
	group.unique<-unique(input_groups)
	
	for(i in 1:length(group.unique)){
		my.g.data<-subset(my.expr,group==group.unique[i])
		my.g.data$group_id<-i
		new.expr<-rbind(new.expr,my.g.data)
	}
	
	z.mean<-mean(new.expr$delta.ct,na.rm=T)
	z.sd<-sd(new.expr$delta.ct,na.rm=T)
	
	new.expr$z_score<-as.numeric((new.expr$delta.ct-z.mean)/z.sd)*-1
	new.expr$z_score[is.na(new.expr$z_score)==T]<-min(new.expr$z_score,na.rm=T)
	
	#######remove genes from the analysis####################
	if(is.character(removed.genes.list)){
		if(length(removed.genes.list) >0 ){
			removed.index<-which(new.expr$gene %in% removed.genes.list)
			new.expr<-new.expr[-c(removed.index),]
		}
	}	
	
	m.m <- dcast(new.expr, gene~cell.name,value.var="z_score",na.rm=T)
	m.g <- dcast(new.expr, gene~cell.name,value.var="group_id",na.rm=T)
	x.m<-as.matrix(m.m[,2:ncol(m.m)])
	g.m<-m.g[,2:ncol(m.g)]
	rownames(x.m)<-m.m$gene
	my_pca<-prcomp(t(x.m))
	
	######create data frame with scores###########
	scores = as.data.frame(my_pca$x)
	PCA.dat<-data.frame(comp1=scores[,c(1)],comp2=scores[,c(2)],comp3=scores[,c(3)],group=as.integer(g.m[1,]))
	#######PCA Plots##############################
	pdf(file=output_file, width=12, height=8, onefile=T, bg="transparent")
	
	#######2D plot################################
	plot(PCA.dat$comp1,PCA.dat$comp2,col=PCA.dat$group,pch=19,main="PCA analysis",xlab="Component 1",ylab="Component 2")
	legend("bottomright", legend =group.unique, col=1:length(group.unique),pch = 19 )
	###########3D plot###########################
	library(scatterplot3d)
	s3d<-scatterplot3d(PCA.dat[,1:3],        
			color=PCA.dat$group, pch=19,      
			type="h", lty.hplot=2,angle= 45,     
			scale.y=0.75, box=FALSE,               
			main="PCA Analysis",
			xlab="Component 1",
			ylab="Component 2",
			zlab="Component 3")
	legend("topright", legend =group.unique, col=1:length(group.unique),pch = 19 )
	
	dev.off()
}