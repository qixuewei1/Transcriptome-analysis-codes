#!/usr/bin/env Rscript
times<-Sys.time();
options(scipen = 200)
options(bitmapType='cairo');
library('getopt');
spec = matrix(c(
		'matrix','m',1,'character',
		'group','g',1,'character',
		'outdir','o',1,'character',
		'help', 'h', 0, 'logical'
			), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
		cat(getopt(spec, usage=TRUE));
	cat("
		Usage example: 
			Rscript DESeq2.R -m counts.txt -g group.txt -c condition.txt -o outdir 

		Usage:
		   -m	--matrix	base express matrix file
		   -g	--group		group file with header
		   -o	--outdir	output dir
		   -h   --help	usage
			  \n")
				q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$matrix) ) { print_usage(spec) }
if ( is.null(opt$outdir) ) { print_usage(spec) }
if ( is.null(opt$group) ) { print_usage(spec) }
if(!dir.exists(opt$outdir)){dir.create(opt$outdir)}



scatter_and_volcano_pic = function(file_name,names){
   group_list=c(unlist(strsplit(names,"_vs_")))
   data = read.table(file_name,sep="\t",head=T,row.names=1)
   meana=data[,colnames(data)==paste(group_list[1],"mean",sep="_")]
   meanb=data[,colnames(data)==paste(group_list[2],"mean",sep="_")]
   FDR =data[,colnames(data)=="FDR"]
   log2FC = data[,colnames(data)=="log2FC"]
   data = data.frame(meana=meana,meanb=meanb,FDR=FDR,log2FC=log2FC)
    data_up = intersect(which(data$FDR<0.05),which(data$log2FC>=1))
    data_down = intersect(which(data$FDR<0.05),which(data$log2FC<=(-1)))
    significant <-rep("normal",times=nrow(data))
    significant[data_up]<-"up"
    significant[data_down]<-"down"
    significant <- factor(significant,levels=c("up","down","normal"))
    p <- ggplot(data,aes(x=meana,y=meanb,colour=significant))+geom_point()+labs(xlab=paste0(group_list[1],"_mean_count"),ylab=paste0(group_list[2],"_mean_count"),size=I(0.7))
    p <- p+scale_color_manual(values=c("up"="red","normal"="black","down"="green"))
    p <- p+theme_bw()+theme(panel.background=element_rect(colour="black",size=1,fill="white"),panel.grid=element_blank())
    pdf(paste0(names,"_scatter.pdf"))
	print(p)
    dev.off()
    png(paste0(names,"_scatter.png"))
	print(p)
    dev.off()
    p1 <- ggplot(data,aes(x=log2FC,y=-log10(FDR),colour=significant))+geom_point()+labs(xlab="log2FC",ylab="-log10(FDR)",size=I(0.7))
    p1 <- p1+scale_color_manual(values=c("up"="red","normal"="black","down"="green"))
    xline=c(-log2(2),log2(2))
    p1 <- p1+geom_vline(xintercept=xline,lty=2,size=I(0.2),colour="grey11")
    yline = -log(0.05,10)
    p1 <- p1+geom_hline(yintercept=yline,lty=2,size=I(0.2),colour="grey11")
    p1 <- p1+theme_bw()+theme(panel.background=element_rect(colour="black",size=1,fill="white"),panel.grid=element_blank())
    pdf(paste0(names,"_volcano.pdf"))
	print(p1)
    dev.off()
    png(paste0(names,"_volcano.png"))
	print(p1)
    dev.off()
}





DESeq2_DE_for_countmatrix = function(database,myfactors){
	database = ceiling(database)
    group_list = myfactors$group_id[!duplicated(myfactors$group_id)]
    fristFactor = myfactors[c(myfactors$group_id==group_list[1]),]
    fristData = database[,colnames(database) %in% fristFactor$id];
    len1 = length(fristFactor[,1])
    secondFactor = myfactors[c(len1+1:length(myfactors[,1]))-len1,]
    secondData = database[,colnames(database) %in% secondFactor$id];
	group1mean = apply(fristData,1,mean)
	group2mean = apply(secondData,1,mean)
    names = paste(group_list[1],group_list[2],sep="_vs_")
    data_matrix <- DESeqDataSetFromMatrix(countData=database,colData=myfactors, design=~ group_id)
    express<-assay(data_matrix)
    data_matrix<-DESeq(data_matrix,parallel=T)
    gid<-unique(myfactors$group_id)
    rld <- rlogTransformation(data_matrix)
    write.table(file=paste(names,"expres.norml.matrix",sep="_"),assay(rld),quote=F,sep="\t")
	express_new<-assay(rld)
	res<-results(data_matrix,parallel=T)
	res$group1mean = group1mean
	res$group2mean = group2mean
	res <- res[order(res$pvalue),]
    colnames(res) = c("baseMean","log2FC","lfcSE","stat","P_value","FDR","group1mean","group2mean")
	names(res)[names(res)=="group1mean"]<-paste(group_list[1],"mean",sep="_")
	names(res)[names(res)=="group2mean"]<-paste(group_list[2],"mean",sep="_")
    write.table(data.frame("gene_id"=rownames(res),res),file=paste(names,"results.csv",sep="_"),col.names=T,quote =F,row.names=F,sep = "\t")
	up=res[which(res$log2FC > 1 & res$FDR <0.05),]
	down=res[which(res$log2FC < -1 & res$FDR <0.05),]
	diff <- express_new[rownames(express_new) %in% rownames(up) | rownames(express_new) %in% rownames(down),]
	diff=data.frame(diff)
	write.table(data.frame("gene_id"=rownames(diff),diff),file=paste(names,"diff.results.csv",sep="_"),col.names=T,row.names=F,quote=F,sep="\t")
	#write.table(diff,file=paste(names,"diff.results.csv",sep="_"),quote=F,sep="\t")
	#resdt.dds = write.table(as(res, "data.frame"),keep.rownames = TRUE)
	#setnames(resdt.dds, "rn", "Ensembl_ID")
	#resdt.dds[, Significant := padj < .1]
	#resdt.dds[!is.na(Significant)]
	#write.table(res.dds,file=paste(names,"out.xls",sep="_"),quote=F,sep="\t")
	scatter_and_volcano_pic(file_name=paste(names,"results.csv",sep="_"),names=names)
}



setwd(opt$outdir)
library(DESeq2)
library("RColorBrewer")
library("ggplot2")
database = read.table(na.omit(opt$matrix), sep = "\t", header = T, row.names = 1)
myfactors = read.table(opt$group,header=T,sep="\t")
n_sample = length(myfactors[,1])
##################
    colnames(myfactors) = c("id","group_id");
    group_list = myfactors$group_id[!duplicated(myfactors$group_id)];
    if(length(group_list)==2){
		newdatabase = database[,colnames(database) %in% myfactors$id]
        DESeq2_DE_for_countmatrix(newdatabase,myfactors);
	}else{
	n = 1;
	for(i in group_list){
	    n = n+1;
		if(n < length(group_list) ){
	    for(m in group_list[n:length(group_list)]){
            newfactors = myfactors[myfactors$group_id %in% c(i,m),];
	        newdatabase = database[,colnames(database) %in% newfactors$id];
            DESeq2_DE_for_countmatrix(newdatabase,newfactors);
		}
		}
	}
	}
################
