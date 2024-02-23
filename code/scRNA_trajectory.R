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
	-f, --file	rds file resulted from scRNA_typing.R.
	-c, --cell	cell names will be used (cell1,cell2,..) default all cell type in rds file.
	-g, --group	Regroup cells into a different identity class by State, sample, group, cell_type or seurat_clusters (e.g. State,group,cell_type. default State).
	-t, --thread	thread number (default 10).
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
	'cell','c',1,'character',
	'group','g',1,'character',
	'thread','t',1,'integer',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$out)) {opt$out='./'}
if (is.null(opt$thread)) {opt$thread=10}
if (is.null(opt$group)) {opt$group='State'}
dir <- opt$out
if(! dir.exists(dir)){
	dir.create(dir, recursive=T)
}
##########################
#define functions

tjtrf <- function(cl, thread){
	library(monocle)
	library(Seurat)
	library(ggplot2)
	library(cowplot)
	subcds <- cds[, rownames(ctd)[ctd$cell %in% cl]]
	subcds <- detectGenes(subcds, min_expr = 0.1)
	Biobase::fData(subcds)$use_for_ordering <- Biobase::fData(subcds)$num_cells_expressed > 0.1 * ncol(subcds)
	expressed_genes <- row.names(subset(Biobase::fData(subcds), num_cells_expressed >= ncol(subcds)/10))
	subcds <- reduceDimension(subcds[expressed_genes, ], max_components = 2, norm_method = 'log', num_dim = 3, reduction_method = 'tSNE', verbose = F, perplexity =10)
	save.image()
	subcds <- clusterCells(subcds, verbose = F)
	clustering_DEG_genes <- differentialGeneTest(subcds[expressed_genes, ], fullModelFormulaStr = '~Cluster', cores=thread)
	ordering_genes <- row.names(clustering_DEG_genes)[clustering_DEG_genes$pval <0.05 ]
	subcds <- setOrderingFilter(subcds, ordering_genes)
	if(ncol(subcds)<500){
		subcds <- reduceDimension(subcds, max_components = 2, method = 'DDRTree', auto_param_selection=F)
	}else{
		subcds <- reduceDimension(subcds, max_components = 2, method = 'DDRTree')
	}
	subcds <- orderCells(subcds)
	od1 <- data.frame(Sample=subcds@phenoData@data$orig.ident, subcds@phenoData@data[, c('Pseudotime', 'State', 'group')])
	od2 <- as.data.frame(t(subcds@reducedDimS))
	names(od2) <- c('Component_1', 'Component_2')
	od3 <- merge(od1, od2, by=0)
	odm <- merge(od3, ctd, by.x=1, by.y=0)
	names(odm)[1] <- 'Cell_Id'
	rownames(odm) <-odm$Cell_Id
	write.table(odm, file=paste(opt$out, '/', gsub(' ', '_', as.character(cl)), '_trajectories.txt', sep=''), quote=F, row.names=F, sep='\t')
	p1 <- plot_cell_trajectory(subcds, color_by = "State")
	p2 <- plot_cell_trajectory(subcds, color_by = "Pseudotime")
	if('cell_type' %in% names(subcds@phenoData@data) && length(unique(subcds@phenoData@data$cell_type)) >1){
		p1$layers[[2]]$mapping <- aes(colour=State, shape=cell_type)
		p2$layers[[2]]$mapping <- aes(colour=Pseudotime, shape=cell_type)
	}else if('group' %in% names(subcds@phenoData@data) && length(unique(subcds@phenoData@data$group)) >1){
		p1$layers[[2]]$mapping <- aes(colour=State, shape=group)
		p2$layers[[2]]$mapping <- aes(colour=Pseudotime, shape=group)
	}
	p1 <- p1+theme(axis.text=element_text(size=15, color='black'), axis.title=element_text(size=20), axis.line=element_blank(), panel.border=element_rect(color='black', fill=NA))
	p2 <- p2+theme(axis.text=element_text(size=15, color='black'), axis.title=element_text(size=20), axis.line=element_blank(), panel.border=element_rect(color='black', fill=NA))
	p <- plot_grid(p1, p2, ncol=2, rel_widths=c(1, 1), align="h", axis="tb")
	ggsave(paste(opt$out, '/', gsub(' ', '_', as.character(cl)), '_trajectories.pdf', sep=''), plot=p, width=13)
	ggsave(paste(opt$out, '/', gsub(' ', '_', as.character(cl)), '_trajectories.png', sep=''), plot=p, width=13, dpi=600)
}

#main programe
library(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)
library(monocle)
library(Hmisc)
#library(foreach)
#library(doParallel)


otherCDS <- readRDS(opt$file)
if(!is.null(opt$cell)){
	cells <- strsplit(opt$cell, ',')[[1]]
	otherCDS <- subset(otherCDS, idents=cells)
}
count <- otherCDS@assays$RNA@counts
pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
f.data <- data.frame(gene_short_name = row.names(count), row.names = row.names(count))
fd <- new("AnnotatedDataFrame", data = f.data)
if(all(count %% 1 == 0)) {
	expressionFamily <- negbinomial.size()
} else if(any(count < 0)){
	expressionFamily <- uninormal()
} else {
	expressionFamily <- tobit()
}
cds <- newCellDataSet(count, phenoData = pd, featureData = fd, expressionFamily=expressionFamily, lowerDetectionLimit=0.1)
cds <- estimateSizeFactors(cds)
ctd <- data.frame(otherCDS@active.ident, stringsAsFactors=F)
names(ctd) <- 'cell'

cds <- detectGenes(cds, min_expr = 0.1)
Biobase::fData(cds)$use_for_ordering <- Biobase::fData(cds)$num_cells_expressed > 0.1 * ncol(cds)
expressed_genes <- rownames(subset(Biobase::fData(cds), num_cells_expressed >= ncol(cds)/10))
cds <- reduceDimension(cds[expressed_genes, ], max_components=2, norm_method='log', num_dim=3, reduction_method='tSNE', verbose=F, perplexity=10)
save.image()
cds <- clusterCells(cds, verbose = F)
clustering_DEG_genes <- differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr='~Cluster', cores=opt$thread)
ordering_genes <- rownames(clustering_DEG_genes)[clustering_DEG_genes$pval < 0.05]
cds <- setOrderingFilter(cds, ordering_genes)
if(ncol(cds)<500){
	cds <- reduceDimension(cds, max_components=2, method='DDRTree', auto_param_selection=F)
}else{
	cds <- reduceDimension(cds, max_components=2, method='DDRTree')
}
cds <- orderCells(cds)
od1 <- data.frame(Sample=cds@phenoData@data$orig.ident, cds@phenoData@data[, c('Pseudotime', 'State', 'group', 'sample', 'seurat_clusters', 'cell_type')], stringsAsFactors=F)
od1$cell_type <- as.character(od1$cell_type)
od1$seurat_clusters <- as.character(as.numeric(od1$seurat_clusters))
od2 <- as.data.frame(t(cds@reducedDimS), stringsAsFactors=F)
names(od2) <- c('Component_1', 'Component_2')
od3 <- merge(od1, od2, by=0)
odm <- merge(od3, ctd, by.x=1, by.y=0)
names(odm)[1] <- 'Cell_Id'
rownames(odm) <- odm$Cell_Id
write.table(odm, file=paste(opt$out, 'trajectories.txt', sep=''), quote=F, row.names=F, sep='\t')

groups <-strsplit(opt$group, ',')[[1]]
for(i in c('Pseudotime', groups)){
	if(!class(cds@phenoData@data[, i]) %in% c('numeric', 'integer')){
		nr <- ceiling(nchar(paste((unique(cds@phenoData@data[, i])), collapse='     '))/60)
		p <- plot_cell_trajectory(cds, color_by=i) + guides(color= guide_legend(title=(capitalize(gsub('_', ' ', i))), nrow=nr))
	}else{
		p <- plot_cell_trajectory(cds, color_by=i)
	}
	p <- p + theme(axis.text=element_text(size=15, color='black'), axis.title=element_text(size=20), axis.line=element_blank(), panel.border=element_rect(color='black', fill=NA))
	ggsave(paste(opt$out, gsub(' ', '_', i), '_trajectories.pdf', sep=''), plot=p, width=6.5, height=7, units='in')
	ggsave(paste(opt$out, gsub(' ', '_', i), '_trajectories.png', sep=''), plot=p, width=6.5, height=7, units='in', dpi=600)
}

diff_test_res <- differentialGeneTest(cds[rownames(cds) %in% rownames(clustering_DEG_genes[!is.na(clustering_DEG_genes$gene_short_name),]),], fullModelFormulaStr = "~sm.ns(Pseudotime)", cores=opt$thread)
sig_gene_names <- rownames(diff_test_res %>% top_n(n=-50, wt=qval))
p <- plot_pseudotime_heatmap(cds[sig_gene_names,], num_clusters = 2,  show_rownames = T, return_heatmap=T)
ggsave(paste(opt$out, 'top_50_gene_pseudotime_heatmap.png', sep=''), plot=p, width=7, height=7, units='in', dpi=600)
ggsave(paste(opt$out, 'top_50_gene_pseudotime_heatmap.pdf', sep=''), plot=p, width=7, height=7, units='in')
save.image()
#BEAM_res <- BEAM(cds, branch_point=3)
#p <- plot_genes_branched_heatmap(cds[rownames(BEAM_res %>% top_n(-50, wt=qval)),], branch_point=3, num_clusters=3, show_rownames=T, return_heatmap=T)
#ggsave(paste(opt$out, 'top_50_gene_branched_heatmap.png', sep=''), plot=p$ph_res, width=7, height=7, units='in', dpi=600)
#ggsave(paste(opt$out, 'top_50_gene_branched_heatmap.pdf', sep=''), plot=p$ph_res, width=7, height=7, units='in')
