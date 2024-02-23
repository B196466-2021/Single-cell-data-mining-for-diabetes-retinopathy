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
	-f, --file	rds file resulted from scRNA_cell_typing.R.
	-c, --comp	compare information file.
	-n, --cell	names of cell type will be analysed (cell1,cell2,..)
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
	'comp','c',1,'character',
	'cell','n',1,'character',
	'thread','t',1,'integer',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || is.null(opt$comp) || !is.null(opt$help)) {usage(spec)}

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
	markers$cluster <- x
	markers$gene <- rownames(markers)
	markers <- filter(markers, max_pval<=0.05)
	return(markers)
}

dexp <- function(x){
	if(length(x) < 2){
		odir <- paste(opt$out, '/', gsub(' ', '_', x), '/', sep='')
		if(! dir.exists(paste(odir, 'sep', sep=''))){
			dir.create(paste(odir, 'sep', sep=''), recursive=T)
		}
	}else{
		odir <- paste(opt$out, '/', sep='')
		if(! dir.exists(paste(odir, 'sep', sep=''))){
			dir.create(paste(odir, 'sep', sep=''), recursive=T)
		}	
	}
	sub_id <- rownames(sub2@meta.data)[sub2@meta.data$group %in% c(comp$Test[i], comp$Control[i])]
	sub_cells <- subset(sub2, idents = x, cells=sub_id)
	sub_data <- as.data.frame(sub_cells$RNA@data)
	sub_cells@meta.data$group <- factor(sub_cells@meta.data$group, levels=c(comp$Test[i], comp$Control[i]))
	sub_data_out <- data.frame(Gene=rownames(sub_data), sub_data)
	write.table(sub_data_out, file=paste(odir, comp$Test[i], '_vs_', comp$Control[i], '_expression.txt', sep=''), quote=F, row.names=F, sep='\t')
	sub_pheno <- data.frame(Sample=rownames(sub_cells@meta.data), Group=sub_cells@meta.data$group)
	write.table(sub_pheno, file=paste(odir,  comp$Test[i], '_vs_', comp$Control[i], '_pheno.txt', sep=''), quote=F, row.names=F, sep='\t')
	sub_cells <- ScaleData(sub_cells, verbose = FALSE)
	avg_group <- log2(AverageExpression(sub_cells, verbose = FALSE, add.ident='group' )$RNA + 1)
	names(avg_group) <- gsub(paste(x, '_', sep=''), '', names(avg_group), fixed=T)
	avg_group$log2FC <- avg_group[, comp$Test[i]] - avg_group[, comp$Control[i]]
	diff.gene <- FindMarkers(sub_cells, ident.1 = comp$Test[i], ident.2 = comp$Control[i], verbose = FALSE, group.by='group', test.use="MAST", logfc.threshold =0, min.pct=0.25)
	mdg <- merge(diff.gene, avg_group, by=0)
	all_diff <- data.frame(Gene=mdg[,1], log2FC=mdg$log2FC, p_value=mdg$p_val, q_value=mdg$p_val_adj, mdg[,c(comp$Test[i], comp$Control[i])], check.names=F)
	write.table(all_diff, file=paste(odir, comp$Test[i], '_vs_', comp$Control[i], '.txt', sep=''), quote=F, row.names=F, sep='\t')
	up_diff <- filter(all_diff, log2FC>=log2(comp$fc[i]) & p_value<=comp$p_value[i] & q_value<=comp$q_value[i])
	if(nrow(up_diff)>0){
		write.table(up_diff, file=paste(odir, 'sep/up_', comp$Test[i], '_vs_', comp$Control[i], '.txt', sep=''), quote=F, row.names=F, sep='\t')
	}else{
		print(paste('Warring: No up gene for ', x, sep=''))
	}
	down_diff <- filter(all_diff, log2FC<=(-log2(comp$fc[i])) & p_value<=comp$p_value[i] & q_value<=comp$q_value[i])
	if(nrow(down_diff)>0){
		write.table(down_diff, file=paste(odir, 'sep/down_', comp$Test[i], '_vs_', comp$Control[i], '.txt', sep=''), quote=F, row.names=F, sep='\t')
	}else{
		print(paste('Warring: No down gene for ', x, sep=''))
	}

	p <- DoHeatmap(sub_cells, features=as.character(c(up_diff$Gene, down_diff$Gene)), group.by='group', label=F)
	levels(p$data$Identity) <- c(comp$Test[i], comp$Control[i])
	p$labels$colour <- 'Group'
	if(length(c(up_diff$Gene, down_diff$Gene))>80){
		p <- p + theme(axis.text.y=element_blank())
	}
	ggsave(paste(odir, 'heatmap_', comp$Test[i], '_vs_', comp$Control[i], '.pdf', sep=''), p, width=8, height=10)
	ggsave(paste(odir, 'heatmap_', comp$Test[i], '_vs_', comp$Control[i], '.png', sep=''), p, width=8, height=10, dpi=600)
}

##########################
#main programe
library(Seurat)
library(cowplot)
library(ggplot2)
library(parallel)
library(MAST)
library(dplyr)
library(foreach)
library(doParallel)
library(pheatmap)
library(gplots)

sub2 <- readRDS(opt$file)

comp <- read.table(opt$comp, header=T, sep='\t', stringsAsFactors=F)

for(i in 1:nrow(comp)){
	###differential expression
	if('group' %in% names(sub2@meta.data) && comp$Test[i] %in% unique(sub2@meta.data$group) && comp$Control[i] %in% unique(sub2@meta.data$group)){
		if(is.null(opt$cell)){
			cell.type <- levels(sub2@active.ident)
			warning('Each cell type will be analysed.')
			mc <- getOption("mc.cores", opt$thread)
			mclapply(cell.type, dexp, mc.cores = mc)
		}else{
			cell.type <- unique(strsplit(opt$cell, ',')[[1]])
			dexp(cell.type)
		}
		#mc <- getOption("mc.cores", opt$thread)
		#mclapply(cell.type, dexp, mc.cores = mc)
	}
}
