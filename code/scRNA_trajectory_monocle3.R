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
	-t, --thread	thread number (default 1).
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
	'thread','t',1,'integer',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$out)) {opt$out='./'}
if (is.null(opt$thread)) {opt$thread=1}

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
	subcds <- cds[, rownames(ctd)[ctd$cell==cl]]
	subcds <- detectGenes(subcds, min_expr = 0.1)
	fData(subcds)$use_for_ordering <- fData(subcds)$num_cells_expressed > 0.05 * ncol(subcds)
	subcds <- reduceDimension(subcds, max_components = 2, norm_method = 'log', num_dim = 3, reduction_method = 'tSNE', verbose = F, perplexity =10)
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
	od1 <- data.frame(Sample=subcds@phenoData@data$orig.ident, subcds@phenoData@data[, c('Pseudotime', 'State')])
	od2 <- as.data.frame(t(subcds@reducedDimS))
	names(od2) <- c('Component_1', 'Component_2')
	odm <- merge(od1, od2, by=0)
	names(odm)[1] <- 'Cell_Id'
	write.table(odm, file=paste(opt$out, '/', gsub(' ', '_', as.character(cl)), '_trajectories.txt', sep=''), quote=F, row.names=F, sep='\t')
	p1 <- plot_cell_trajectory(subcds, color_by = "State")
	p2 <- plot_cell_trajectory(subcds, color_by = "Pseudotime")
	p <- CombinePlots(list(p1, p2))
	ggsave(paste(opt$out, '/', gsub(' ', '_', as.character(cl)), '_trajectories.pdf', sep=''), plot=p, width=13)
	ggsave(paste(opt$out, '/', gsub(' ', '_', as.character(cl)), '_trajectories.png', sep=''), plot=p, width=13, dpi=600)
}

#main programe
#library(Seurat)
#library(cowplot)
#library(ggplot2)
#library(parallel)
library(monocle3)
#library(foreach)
#library(doParallel)

rds <- readRDS(opt$file)
cell_metadata <- data.frame(rds@meta.data, cell.type=rds@active.ident)
cds <- new_cell_data_set(rds@assays$RNA@counts, cell_metadata=cell_metadata)
cds <- preprocess_cds(cds, num_dim=50)
cds <- align_cds(cds, alignment_group='group')
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
p <- plot_cells(cds, color_cells_by = "cell.type", label_groups_by_cluster=F, label_leaves=FALSE, label_branch_points=FALSE, group_label_size=4)

ggsave(paste(opt$out, '/trajectories.pdf', sep=''), plot=p)
ggsave(paste(opt$out, '/trajectories.png', sep=''), plot=p, dpi=600)

