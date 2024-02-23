#! /usr/bin/Rscript

###########################
##! @Author: Dehai Wang
##! @Todo: 
##! @Version: 1.0.0
##! @Dep: R
##! @ChangeLog:

###########################
#sub functions

version<-function(){
	cat("author: wangdh
version: 1.0.0
updated date: 2021-04-23\n\n")
	q(status=1)
}

usage<-function(spec){
	cat("This script is used to infer CNV for scRNA.\n",
getopt(spec,usage=TRUE),
"Options:
	-e, --exp	file of genes (rows) vs. cells (columns) containing the raw counts.
	-c, --cell	cell iformation file each row is one cell indicating the cell type classifications (cell<tab>group).
	-g, --gene	file containing the positions of each gene along each chromosome in the genome (gene<tab>chr<tab>start<tab>end).
	-r, --ref	the classifications of the reference (normal) cells to use for infering cnv (e.g. group1,group2... default NULL).
	-l, --cluster	cluster number of the observation (default 2).
	-t, --thread	thread number (default 1).
	-o, --out	prefix of out files.
	-v, --version	display version information.
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1)
}

##########################
#get options
library(getopt)
spec = matrix(c(
	'exp','e',1,'character',
	'gene','g',1,'character',
	'cell','c',1,'character',
	'ref','r',1,'character',
	'thread','t',1,'integer',
	'cluster','l',1,'integer',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$exp) || is.null(opt$gene) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$cluster)) {opt$cluster=2}
if (is.null(opt$thread)) {opt$thread=1}
if (is.null(opt$out)) {opt$out='./'}

dir <- gsub("[^/]*$", '', opt$out, perl=T)
if(! dir.exists(dir) && dir != ''){
	dir.create(dir, recursive=T)
}
##########################
#define functions

my_fun<-function(){

}

##########################
#main programe
library(gtools)
library(infercnv)
library(dplyr)
exp <- as.matrix(read.table(opt$exp, header=T, sep='\t', row.names=1, check.names=F))
gene_inf <- read.table(opt$gene, header=T, sep='\t', row.names=1, stringsAsFactors=F)
#gene_inf <- gene_inf %>% group_by(gene) %>% top_n(n=-1, wt=chr) %>% group_by(gene) %>% top_n(n=-1, wt=start) %>% group_by(gene) %>% top_n(n=1, wt=end) %>% data.frame(stringsAsFactors=F)
gene_inf$chr <- factor(gene_inf$chr, levels=mixedsort(unique(gene_inf$chr)))
gene_inf <- gene_inf[order(as.numeric(gene_inf$chr), gene_inf$start, gene_inf$end),]
gene_inf$chr <- as.character(gene_inf$chr)

if(is.null(opt$cell)){
	cell_inf <- data.frame(cell_type=rep('all', ncol(exp)))
	rownames(cell_inf) <- colnames(cell_inf)
}else{
	cell_inf1 <- read.table(opt$cell, header=T, sep='\t', stringsAsFactors=F)
	cell_inf1 <- cell_inf1[cell_inf1[,1] %in% colnames(exp), ]
	cell_inf <- data.frame(group=cell_inf1[,2], stringsAsFactors=F)
	rownames(cell_inf) <- cell_inf1[,1]
	print(head(cell_inf))
	exp <- exp[, colnames(exp) %in% rownames(cell_inf)]
}
if(is.null(opt$ref)){
	ref <- NULL
}else{
	ref <- strsplit(opt$ref, ',')[[1]]
}
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=exp, annotations_file=cell_inf, gene_order_file=gene_inf, ref_group_names=ref, chr_exclude = c("chrX", "chrY", "MT"))
infercnv_run <- run(infercnv_obj, out_dir=opt$out, k_obs_groups=opt$cluster, denoise=T, output_format='pdf', num_threads=opt$thread, cluster_references=F, HMM=F, cluster_by_groups=F)
