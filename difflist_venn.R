#!/usr/bin/env Rscript
library(VennDiagram)
library('getopt');
spec = matrix(c(
		'matrix','m',1,'character',
		'outdir','o',1,'character',
		'pre','d',1,'character',
		'help', 'h', 0, 'logical'
		), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("
	Usage example: 
	Rscript DESeq2.R -m counts.txt -g group.txt -c condition.txt -o outdir 
	Usage:
				   -m	--matrix	gene_express_matrix
    			   -o	--outdir	output dir
				   -d   --pre      pre for filename
				   -h   --help	usage\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$matrix) ) { print_usage(spec) }
if ( is.null(opt$outdir) ) { print_usage(spec) }
if(!dir.exists(opt$outdir)){dir.create(opt$outdir)}

setwd(opt$outdir)
datalist = read.table(opt$matrix,header=F,sep="\t",row.names =1)
filelist = row.names(datalist)
listall = list()
for(m in filelist){
    nlist = strsplit(m,split="/")
	fname = nlist[[1]][length(nlist[[1]])]
	sname = strsplit(fname,split="_results")[[1]][1]
	data = read.table(fname,header=F,sep="\t",row.names=1)
	listm = row.names(data)
	listall[[sname]] = listm
}
colours = rainbow(length(filelist))
venn.diagram(x=listall,filename=paste0(opt$pre,"_diff_venn.pdf"),fill=colours,alpha = 0.5,fontfamily = "serif")


