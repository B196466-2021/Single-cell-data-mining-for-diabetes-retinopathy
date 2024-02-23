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
updated date: 2019-08-13\n\n")
	q(status=1)
}

usage<-function(spec){
	cat("This script is used to\n",
getopt(spec,usage=TRUE),
"Options:
	-r, --rds	clustered Seurat rds file result from scRNA_cluster.R.
	-g, --gene	gmt file of genes for each cell type (e.g. cell1<tab>NA<tab>gene1<tab>gene2<tab>...) or gene names (e.g. gene1,gene2,...).
	-t, --thread	thread number (default 5).
	-o, --out	directory of out files (default ./).
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
	'pval','p',1,'double',
	'thread','t',1,'integer',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$rds) || is.null(opt$gene) || !is.null(opt$help)) {usage(spec)}
if (is.null(opt$thread)) {opt$thread=5}
if (is.null(opt$out)) {opt$out='./'}

dir <- paste(opt$out, '/', sep='')
if(! dir.exists(dir) && dir != ''){
	dir.create(dir, recursive=T)
}
##########################
#define functions

map_plot <- function(gene, cell, object, diru, dirt){
	pct <- DimPlot(object, reduction='tsne', label=TRUE, label.size=6, group.by='seurat_clusters', combine=F)[[1]] + scale_color_igv() + theme(plot.background=element_rect(fill="white", colour=NA))
	pcu <- DimPlot(object, reduction='umap', label=TRUE, label.size=6, group.by='seurat_clusters', combine=F)[[1]] + scale_color_igv() + theme(plot.background=element_rect(fill="white",  colour=NA))
	mt <- paste(cell, ': ', gene, sep='')
	pt <- FeaturePlot(object, features=gene, reduction='tsne', combine=F)[[1]] + labs(title=mt) + theme(plot.background=element_rect(fill="white", colour=NA))
	pu <- FeaturePlot(object, features=gene, reduction='umap', combine=F)[[1]] + labs(title=mt) + theme(plot.background=element_rect(fill="white", colour=NA))
	p.tsne <- plot_grid(pct, pt, align='hv')
	p.umap <- plot_grid(pcu, pu, align='hv')
	ggsave(paste(diru, '/', gsub(' |\\/', '_', cell), '_', gene, '_UMAP.png', sep=''), plot=p.umap, height=6, width=13, dpi=100, units='in', device="png")
	ggsave(paste(dirt, '/', gsub(' |\\/', '_', cell), '_', gene, '_tSNE.png', sep=''), plot=p.tsne, height=6, width=13, dpi=100, units='in', device="png")
}

##########################
#main programe
library(SeuratObject)
library(Seurat)
library(cowplot)
library(ggplot2)
library(ggsci)
library(parallel)

sub2 <- readRDS(opt$rds)
DefaultAssay(sub2) <- 'RNA'
if(file_test("-f", opt$gene)){
	gels <- readLines(opt$gene)
	gels <- strsplit(gels, '\t')
	for(i in 1:length(gels)){
		names(gels)[i] <- gels[[i]][1]
		gels[[gels[[i]][1]]] <- unique(gels[[gels[[i]][1]]][-1])
	}
	diru <- paste(dir, '/umap', sep='')
	dirt <- paste(dir, '/tsne', sep='')
	if(! dir.exists(diru)){
		dir.create(diru, recursive=T)
	}
	if(! dir.exists(dirt)){
		dir.create(dirt, recursive=T)
	}
	pct <- DimPlot(sub2, reduction='tsne', label=TRUE, label.size=6, group.by='seurat_clusters', combine=F)[[1]] + scale_color_igv() + theme(plot.background=element_rect(fill="white", colour=NA))
	pcu <- DimPlot(sub2, reduction='umap', label=TRUE, label.size=6, group.by='seurat_clusters', combine=F)[[1]] + scale_color_igv() + theme(plot.background=element_rect(fill="white", colour=NA))
	for(n in names(gels)){
		genes <- gels[[n]][gels[[n]] %in% rownames(sub2@assays$RNA@data)]
		if(length(genes) > 0){
			mc <- getOption("mc.cores", opt$thread)
			mclapply(genes, map_plot, cell=n, object=sub2, diru=diru, dirt=dirt,  mc.cores = mc)
		}
	}
}else{
	gels <- strsplit(opt$gene, ',')[[1]]
	for(gene in gels){
		pct <- FeaturePlot(sub2, features=gene, reduction='tsne', combine=F)[[1]] + theme(plot.background=element_rect(fill="white", colour=NA), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
		pcu <- FeaturePlot(sub2, features=gene, reduction='umap', combine=F)[[1]] + theme(plot.background=element_rect(fill="white", colour=NA), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
		ggsave(paste(dir, '/', gene, '_UMAP.pdf', sep=''), plot=pcu, width=5, height=5)
		ggsave(paste(dir, '/', gene, '_UMAP.png', sep=''), plot=pcu, width=5, height=5, units='in', dpi=600)
		ggsave(paste(dir, '/', gene, '_tSNE.pdf', sep=''), plot=pct, width=5, height=5)
		ggsave(paste(dir, '/', gene, '_tSNE.png', sep=''), plot=pct, width=5, height=5, units='in', dpi=600)
		p <- VlnPlot(sub2, features=gene, pt.size=0, combine=F)[[1]] + xlab(NULL) + theme(plot.background=element_rect(fill="white", colour=NA))
		ggsave(paste(dir, '/', gene, '_violin.pdf', sep=''), plot=p, width=5, height=5)
		ggsave(paste(dir, '/', gene, '_violin.png', sep=''), plot=p, width=5, height=5, units='in', dpi=600)
	}
}


