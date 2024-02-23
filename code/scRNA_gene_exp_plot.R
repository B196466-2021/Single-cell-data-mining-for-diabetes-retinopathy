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
updated date: 2020-12-25\n\n")
	q(status=1)
}

usage<-function(spec){
	cat("This script is used to display gene expression in tSNE and UMAP for scRNA datas.\n",
getopt(spec,usage=TRUE),
"Options:
	-r, --rds	rds file result from scRNA_cluster.R.
	-g, --gene	gene (e.g. gene1,gene2,...) or list file of genes or file with Cluster and Genes columns (e.g. cluster1<Tab>gene1//gene2//...).
	-o, --out	prefix of out files default (./).
	-v, --version	display version information.
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1)
}

##########################
#get options
library(getopt)
spec = matrix(c(
	'rds','r',1,'character',
	'gene','g',1,'character',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$rds) || is.null(opt$gene) || !is.null(opt$help)) {usage(spec)}

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
library(Seurat)
library(ggplot2)

pbmc <- readRDS(opt$rds)
if(! 'seurat_clusters' %in% names(pbmc@meta.data)){
	stop('The scRNA data was not yet clustered')
}

if(file_test('-f', opt$gene)){
	gf <- read.table(opt$gene, header=T, sep='\t', quote="", stringsAsFactors=F)
	if(all(c('Cluster', 'Genes') %in% names(gf))){
		gene.list <- lapply(unique(gf$Cluster), function(x){gene=unique(gf$Genes[gf$Cluster==x]); gene2=unique(unlist(strsplit(gene, split="//|,|;|\\|", perl=T))); return(gene2)})
		names(gene.list) <- as.character(unique(gf$Cluster))
	}else if(ncol(gf) == 1){
		gene.list <- list(readLines(opt$gene))
		names(gene.list) <- ''
	}
}else{
	gene.list <- list(unique(unlist(strsplit(opt$gene, split="//|,|;|\\|", perl=T))))
	names(gene.list) <- ''
}
if("tsne" %in% names(pbmc@reductions)){
	for(i in 1:length(gene.list)){
		ncl <- ceiling(sqrt(length(gene.list[[i]])))
		nrw <- ceiling(length(gene.list[[i]])/ncl)
		p <- FeaturePlot(pbmc, features=gene.list[[i]], ncol=ncl, reduction='tsne')
		ggsave(paste(opt$out, names(gene.list)[i], '_marker_gene_exp_tSNE.pdf', sep=''), plot=p, width=ncl*4, height=nrw*4, units='in')
		ggsave(paste(opt$out, names(gene.list)[i], '_marker_gene_exp_tSNE.png', sep=''), plot=p, width=ncl*4, height=nrw*4, units='in', dpi=600)
	}
}

if("umap" %in% names(pbmc@reductions)){
	for(i in 1:length(gene.list)){
		ncl <- ceiling(sqrt(length(gene.list[[i]])))
		nrw <- ceiling(length(gene.list[[i]])/ncl)
		p <- FeaturePlot(pbmc, features=gene.list[[i]], ncol=ncl, reduction='umap')
		ggsave(paste(opt$out, names(gene.list)[i], '_marker_gene_exp_UMAP.pdf', sep=''), plot=p, width=ncl*4, height=nrw*4, units='in')
		ggsave(paste(opt$out, names(gene.list)[i], '_marker_gene_exp_UMAP.png', sep=''), plot=p, width=ncl*4, height=nrw*4, units='in', dpi=600)
	}
}



