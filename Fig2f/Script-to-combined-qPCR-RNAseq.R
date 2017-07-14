###############################################################################
my.dat<-read.delim("qPCR.Diff.Genes.Final.Result.txt",header=T)
##################Get RNA-seq DE genes#########################
ox1407<-read.table("DE-Genes-BCR-ABL+_Vs_BCR-ABL-_OX1407.txt",header=T)
###############################################################################
ox1407$isInqPCR<-"FALSE"
ox1407$isInqPCR[which(ox1407$gene %in% my.dat$gene)]<-"TRUE"
###############################################################
my.merge<-merge(ox1407,my.dat,all.x=T,sort=FALSE)
my.merge<-my.merge[order(my.merge$fisher),]
#my.merge<-na.omit(my.merge)
my.merge<-subset(my.merge,isInqPCR==TRUE)

write.table(my.merge,file="Combined.qPCR.RNA-seq.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
