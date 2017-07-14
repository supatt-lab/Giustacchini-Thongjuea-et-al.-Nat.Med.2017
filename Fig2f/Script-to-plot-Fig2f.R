#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
#############################################################
#This R script is used to generate Fig2f.
##############################################################
df<-read.table("Combined.qPCR.RNA-seq.txt",header=T)
df<-df[c("gene","delta_delta_ct","log2fc")]
colnames(df)<-c("gene","x","y")
df1<-subset(df,abs(y)>=3.32)

selected.genes<-c("B2M","GAPDH","HSP90AB1","CD34","CTNNB1","HPRT1","CD33","BCR","ABL1","IGF1R","SELP","BCL2","CDK6","ATG3","SAT1")
df2<-df[df$gene %in% selected.genes,]

#############combined two gene sets#####
my.df<-rbind(df1,df2)
################save data matrix to a file####################################
write.table(my.df,file="Fig2f.data.matrix.txt",
		append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
##############################################################################

nosig.index<-which(abs(my.df$x) < 2 & abs(my.df$y) < 2)

pdf(file="Fig2f.pdf", width=6, height=6, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)

mod <- lm(y ~ x, data = my.df)
newx <- seq(min(my.df$x), max(my.df$x), length.out=100)
preds <- predict(mod, newdata = data.frame(x=newx),interval = 'confidence')

# plot
plot(my.df$x,my.df$y, type = 'n',xlab="qPCR (BCR-ABL- Vs BCR-ABL+), Log2fc",
		ylab="RNA-seq (BCR-ABL- Vs BCR-ABL+), Log2fc",xlim=c(-4,6),ylim=c(-4,6),
		main="qPCR and RNA-seq correlation")
points(my.df$x,my.df$y,pch=19,cex=1.2,col="brown")
points(my.df$x[nosig.index],my.df$y[nosig.index],pch=19,cex=1.2,col="grey")
abline(mod,lty=2,lwd=1.5,col="red")
abline(v=0,lty=2,lwd=1.5,col="grey")
abline(h=0,lty=2,lwd=1.5,col="grey")
text(my.df$x+0.2,my.df$y+0.2,labels=my.df$gene,cex=0.6)
text(3,-4,paste("Pearson correlation= ",format(cor(my.df$x,my.df$y),digit=3),sep=""))
dev.off()