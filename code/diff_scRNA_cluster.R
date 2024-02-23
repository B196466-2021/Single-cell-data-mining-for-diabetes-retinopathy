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
	cat("This script is used to cluster cells of scRNA.\n",
getopt(spec,usage=TRUE),
"Options:
	-f, --file	all rds file(s file1,file2,...) resulted from scRNA_QC.R.
	-t, --thread	thread number (default 1).
	-s, --tissue	tissue type of the samples (http://biocc.hrbmu.edu.cn/CellMarker/index.jsp#).
	-o, --out	out directory (default ./).
	-v, --version	display version information.
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1)
}

##########################
#get options
library(getopt)
spec = matrix(c(
	'file','f',1,'character',
	'tissue','s',1,'character',
	'thread','t',1,'integer',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || is.null(opt$tissue) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$out)) {opt$out='./'}
if (is.null(opt$thread)) {opt$thread=1}

dir <- opt$out
if(! dir.exists(dir)){
	dir.create(dir, recursive=T)
}
##########################
#define functions

my_fun<-function(x){
	markers <- FindConservedMarkers(sub, ident.1 = x, grouping.var = "group", verbose = FALSE, test.use='MAST', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	navlfc <- grep('_avg_logFC', names(markers))
	logFC.all <- data.frame(markers[, navlfc])
	avg_logFC <- apply(logFC.all, 1, mean)
	npvalue <- grep('_p_val$', names(markers))
	pvalue.all <- data.frame(markers[, npvalue])
	max_pval <- apply(pvalue.all, 1, max)
	markers <- data.frame(gene=rownames(markers), cluster=x, avg_logFC=avg_logFC, max_pval=max_pval, stringsAsFactors=F)
	if(nrow(markers[max_pval<=0.05, ]) > 3){
		markers <- markers[max_pval<=0.05,]
	}else if(nrow(markers[max_pval<=0.1, ]) > 3){
		markers <- markers[max_pval<=0.1,]
	}
	if(nrow(markers)>0){
		return(markers)
	}
}

##########################
#main programe
library(Seurat)
library(cowplot)
library(ggplot2)
library(parallel)
library(MAST)
library(dplyr)
#library(monocle)
library(foreach)
library(doParallel)
library(pheatmap)
library(gplots)
source('/home/dev/DEV-wangdh/seurat_functions.R')

set.list<-list()
fs <- strsplit(opt$file, ',')
for(f in fs){
	fr <- readRDS(f)
	if(is.list(fr)){
		for(n in names(fr)){
			set.list[[n]] <- fr[[n]]
		}
	}else{
		set.list[[unique(fr$group)]] <- fr
	}
}

#comp <- read.table(opt$comp, header=T, sep='\t', stringsAsFactors=F)

#for(i in nrow(comp)){
#	c_dir <- paste(opt$out, '/', comp$Test[i], '_vs_', comp$Control[i], sep='')
#	if(!dir.exists(c_dir)){
#		dir.create(c_dir, recursive=T)
#	}
#	sub.list <- set.list[names(set.list) %in% c(comp$Test[i], comp$Control[i])]

sub.list <- set.list
c_dir <- paste(opt$out, '/cluster_cell', sep='')
if(!dir.exists(c_dir)){
	dir.create(c_dir, recursive=T)
}

	vfp.list <- list()
	for(g in names(set.list)){
		sub.list[[g]] <-  FindVariableFeatures(sub.list[[g]], selection.method = "vst", nfeatures = 2000)
		top10 <- head(VariableFeatures(sub.list[[g]]), 10)
		p <- VariableFeaturePlot(sub.list[[g]])
		p <- LabelPoints(plot = p, points = top10, repel = TRUE)
		p <- p + ggtitle(label = g)
		vfp.list[g] <- list(p)
	}
	ngp <- length(set.list)
	ncg <- ceiling(sqrt(ngp))
	nrg <- ceiling(ngp/ncg)
	wg <- ncg*4
	hg <- nrg*4
	p <- CombinePlots(plots = vfp.list, legend='none', ncol=ncg)
	ggsave(paste(c_dir, '/feature_select.pdf', sep=''), plot=p, width=wg, height=hg)
	ggsave(paste(c_dir, '/feature_select.png', sep=''), plot=p, height=hg, width=wg, dpi=600, units='in')
	sub <- FindIntegrationAnchors(object.list = sub.list, dims = 1:20)
	sub <- IntegrateData(anchorset = sub, dims = 1:20)
	DefaultAssay(sub) <- "integrated"
	sub <- ScaleData(sub, verbose = FALSE)
	##PCA
	sub <- RunPCA(sub, verbose = FALSE)
	p <- VizDimLoadings(sub, dims = 1:2, reduction = "pca")
	ggsave(paste(c_dir, '/pc_top_gene.pdf', sep=''), plot=p)
	ggsave(paste(c_dir, '/pc_top_gene.png', sep=''), plot=p, height=7, width=7, dpi=600, units='in')
	p <- DimPlot(sub, reduction = "pca", group.by='group')
	ggsave(paste(c_dir, '/PCA.pdf', sep=''), plot=p)
	ggsave(paste(c_dir, '/PCA.png', sep=''), plot=p, height=7, width=7, dpi=600, units='in')
	##PC select
	y <- sub$pca@stdev
	x <- seq(1, 50)
	d1n <- diff(y)/diff(x)
	right.edge <- which.min(d1n)
	left.edge <- which.max(d1n[seq_len(right.edge)])
	fit <- smooth.spline(x, y)
	d1 <- predict(fit, deriv = 1)$y
	d2 <- predict(fit, deriv = 2)$y
	curvature <- d2/(d1^2)^1.5
	npc <- x[which.min(curvature)]
	p <- ElbowPlot(sub, ndims =50)
	p <- p+geom_vline(xintercept = npc, colour = "red", linetype=2)
	ggsave(paste(c_dir, '/PC_select.pdf', sep=''), p, width=9, device="pdf")
	ggsave(paste(c_dir, '/PC_select.png', sep=''), plot=p, height=7, width=7, dpi=600, units='in')
	##PC heatmap
	if(npc>3){
		ncs <- ceiling(sqrt(npc))
		nrs <- ceiling(npc/ncs)
		w <- ncs*3
		h <- nrs*3
	}else{
		ncs <- npc
		h <- 3
		w <- ncs*3
	}
	p <- DimHeatmap(sub, dims = 1:npc, balanced = TRUE, ncol=ncs, fast=F, cells = 500)
	ggsave(paste(c_dir, '/PC_heatmap.pdf', sep=''), p, height=h, width=w, device="pdf")
	ggsave(paste(c_dir, '/PC_heatmap.png', sep=''), plot=p, height=h, width=w, dpi=600, units='in')
	##cluster
	sub <- FindNeighbors(sub, reduction = "pca", dims = 1:npc)
	res <- 0.5
	sub <- FindClusters(sub, resolution = 0.5)
	while(length(levels(sub@active.ident)) > 15){
		res <- res - 0.1
		if(res < 0.1){
			break
		}
		sub <- FindClusters(sub, res =res)
	}
	##TSNE
	sub <- RunTSNE(sub, reduction = "pca", dims = 1:npc)
	p <- DimPlot(sub, reduction = "tsne", group.by = "group")
	p <- p + labs(x='tSNE1', y='tSNE2')
	ggsave(paste(c_dir, '/tSNE_group.pdf', sep=''), plot=p, height=6, width=7)
	ggsave(paste(c_dir, '/tSNE_group.png', sep=''), plot=p, height=6, width=7, dpi=600, units='in')
	p <- DimPlot(sub, reduction = "tsne", label = TRUE, split.by ='group', ncol=ncg)
	p <- p + theme(legend.position='none') + labs(x='tSNE1', y='tSNE2')
	ggsave(paste(c_dir, '/tSNE_cluster.pdf', sep=''), plot=p, height=hg, width=wg)
	ggsave(paste(c_dir, '/tSNE_cluster.png', sep=''), plot=p, height=hg, width=wg, dpi=600, units='in')
	##UMAP
	sub <- RunUMAP(sub, reduction = "pca", dims = 1:npc)
	p <- DimPlot(sub, reduction = "umap", group.by = "group")
	p <- p + labs(x='UMAP1', y='UMAP2')
	ggsave(paste(c_dir, '/UMAP_group.pdf', sep=''), plot=p, height=6, width=7)
	ggsave(paste(c_dir, '/UMAP_group.png', sep=''), plot=p, height=6, width=7, dpi=600, units='in')
	p <- DimPlot(sub, reduction = "umap", label = TRUE, split.by ='group', ncol=ncg)
	p <- p + theme(legend.position='none') + labs(x='UMAP1', y='UMAP2')
	ggsave(paste(c_dir, '/UMAP_cluster.pdf', sep=''), plot=p, height=hg, width=wg)
	ggsave(paste(c_dir, '/UMAP_cluster.png', sep=''), plot=p, height=hg, width=wg, dpi=600, units='in')
	#Identify conserved cell type markers
	DefaultAssay(sub) <- 'RNA'
	sub <- ScaleData(sub, verbose = FALSE)
	mc <- getOption("mc.cores", opt$thread)
	maker.list <- mclapply(as.numeric(levels(sub@active.ident)), my_fun, mc.cores = mc)
	markers <- do.call(rbind, maker.list)
	top20.mk <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
	#navlfc <- grep('_avg_logFC', names(markers))
	marker_gene <- markers[, c('cluster', 'gene', 'avg_logFC')]
	##Assigning cell type identity to clusters
	write.table(markers, file=paste(c_dir, '/maker_gene.txt', sep=''), quote=F, row.names=F, sep='\t')
	if(any(grepl('tumor|cancer', names(set.list), ignore.case=T))){
		sct <- '-c'
	}else{
		sct <- ''
	}
	comd <- paste('Rscript /home/dev/DEV-cuiy/cellTypeEnrich.r -f=', c_dir, '/maker_gene.txt -t=', opt$tissue, ' -o=maker_anno ', sct, sep='')
	xlsfl <- 'maker_anno.xlsx'
	print(comd)
	system(comd)
	makeh <- read.table('maker_anno.txt', header=T, sep='\t')[, c('Cell_Name', 'Cluster', 'p_value')]
	clt <- makeh %>% group_by(Cluster) %>% top_n(n = 1, wt = -p_value)
	cell_type <- as.character(clt$Cell_Name)
	names(cell_type) <- clt$Cluster
	if(any(! levels(sub) %in% names(cell_type))){
		nct <- as.character(levels(sub)[! levels(sub) %in% names(cell_type)])
		cnct <- paste('Cluster_', nct, sep='')
		names(cnct) <- nct
		cell_type <- c(cell_type, cnct)
	}
	sub2 <- RenameIdents(sub, cell_type)
	p <- DimPlot(sub2, reduction = "tsne", label = TRUE, split.by ='group', ncol=ncg)
	p <- p + labs(x='tSNE1', y='tSNE2')
	ggsave(paste(c_dir, '/tSNE_cell_type.pdf', sep=''), plot=p, height=hg, width=wg+2)
	ggsave(paste(c_dir, '/tSNE_cell_type.png', sep=''), plot=p, height=hg, width=wg+2, dpi=600, units='in')
	p <- DimPlot(sub2, reduction = "umap", label = TRUE, split.by ='group', ncol=ncg)
	p <- p + labs(x='UMAP1', y='UMAP2')
	ggsave(paste(c_dir, '/UMAP_cell_type.pdf', sep=''), plot=p, height=hg, width=wg+2)
	ggsave(paste(c_dir, '/UMAP_cell_type.png', sep=''), plot=p, height=hg, width=wg+2, dpi=600, units='in')
	sub@misc<-list(marker.gene=marker_gene)
	saveRDS(sub, file=paste(c_dir, '/cell_cluster.rds', sep=''))
