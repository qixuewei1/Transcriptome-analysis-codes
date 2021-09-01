#!/usr/bin/env Rscript 
times<-Sys.time()
library("scatterplot3d")
library('getopt');
options(bitmapType='cairo')
#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character",
	'group', 'g', 1 , "character",
	'varfile', 'v', 1 , "character"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
Options: 
	--help		NULL 		get this help
	--infile 	character 	the input file [forced]
	--outfile 	character 	the filename for output graph [forced]
	--group		character	the group file for draw
	\n")
	q(status=1);
}
# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }
# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
#if ( is.null(opt$group) )	{ print_usage(spec) }
library(ggplot2)
library(ggrepel)
library(egg)
library(pheatmap)
library("RColorBrewer")
data<-read.table(head=T,sep="\t",file=opt$infile,row.names=1,check.names = FALSE)
colnames(data)<-sub(".isoforms.results","",colnames(data))
colnames(data)<-sub(".genes.results","",colnames(data))
group<-read.table(opt$group,head=F,sep="\t",comment="^",check.names = FALSE)
colnames(group)<-c("id","gid")
pca<-prcomp(data)
write.table(pca$sdev,file=paste(opt$outfile,"_pca_sdev.csv",sep=""),quote=F)
write.table(pca$rotation,file=paste(opt$outfile,"_pca_rotation.csv",sep=""),quote=F)
rotation<-as.data.frame(pca$rotation)
rotation$id=row.names(rotation)
print(rotation)
print(group)
rotation<-merge(rotation,group,by="id")
print(rotation)
sdev=data.frame(sdev=pca$sdev)
sdev$components=paste("PC",rownames(sdev),sep="")
sdev$percent=sdev$sdev/sum(sdev$sdev)
sdev<-ggplot(sdev)+geom_bar(aes(x=components,y=percent),stat="identity",fill="blue")+theme_bw()
pc12<-ggplot(rotation)+geom_point(aes(x=PC1,y=PC2,group=gid,col=gid),shape=0,size=1)+geom_text_repel(aes(x=PC1,y=PC2,label=id))+theme_bw()
pc13<-ggplot(rotation)+geom_point(aes(x=PC1,y=PC3,group=gid,col=gid),shape=0,size=2)+geom_text_repel(aes(x=PC1,y=PC3,label=id))+theme_bw()
pc23<-ggplot(rotation)+geom_point(aes(x=PC2,y=PC3,group=gid,col=gid),shape=0,size=2)+geom_text_repel(aes(x=PC2,y=PC3,label=id))+theme_bw()
pdf(paste(opt$outfile,"pca.pdf",sep="."));
ggarrange(pc12,pc13,pc23,sdev,widths=c(1,1),heights=c(1,1))
dev.off()
png(paste(opt$outfile,"pca.png",sep="."));
ggarrange(pc12,pc13,pc23,sdev,widths=c(1,1),heights=c(1,1))
dev.off()
col <- colorRampPalette(brewer.pal(9, "Blues"))(100)
pearson_cor <- as.matrix(cor(data, method="pearson"))
pdf(paste(opt$outfile,"pearson_cor.pdf",sep="."))
pheatmap(pearson_cor,col=col)
dev.off()
png(paste(opt$outfile,"pearson_cor.png",sep="."))
pheatmap(pearson_cor,col=col)
dev.off()


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
