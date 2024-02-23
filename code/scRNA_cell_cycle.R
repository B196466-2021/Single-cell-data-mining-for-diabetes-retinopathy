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
	-f, --file	rds file resulted from scRNA_cell_typing.R.
	-s, --species	species (default human).
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
	'species','s',1,'character',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$out)) {opt$out='./'}
if (is.null(opt$species)) {opt$species='human'}

dir <- opt$out
if(! dir.exists(dir)){
	dir.create(dir, recursive=T)
}
##########################
#define functions

##########################
#main programe
library(Seurat)
library(cowplot)
library(ggplot2)
library(parallel)
library(dplyr)
library(monocle)
library(foreach)
library(doParallel)
library(pheatmap)
library(gplots)

sub <- readRDS(opt$file)
##Assign Cell-Cycle Scores

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sub2 <- CellCycleScoring(sub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ccgs <- rownames(sub2@assays$RNA@scale.data)[toupper(rownames(sub2@assays$RNA@scale.data)) %in% c(s.genes, g2m.genes)]
ncs <- ceiling(sqrt(length(ccgs)))
nrs <- ceiling(length(ccgs)/ncs)
if(length(ccgs)>3){
	w <- ncs*3
	h <- nrs*3
}else{
	h <- 3
	w <- length(ccgs)*3
}
p <- RidgePlot(sub2, features = ccgs, ncol = ncs)
ggsave(paste(opt$out, '/cell_cycle_ridge_plot.pdf', sep=''), p, device="pdf", width=w, height=h)
ggsave(paste(opt$out, '/cell_cycle_ridge_plot.png', sep=''), p, device="png", width=w, height=h)

#Cell-Cycle heatmap
otop10 <- sub2@misc$otop10
sub_cc <- ScaleData(sub2,otop10$gene)
sub_cc <- ScaleData(sub2,otop10$gene,scale.max=-min(sub_cc@assays$RNA@scale.data))
scd <- sub_cc@assays$RNA@scale.data
anno_col <- sub2@meta.data[, c('S.Score', 'G2M.Score', 'Phase')]
anno_col$Cell.Type <- sub@active.ident
p <- ggplot(anno_col, aes(Cell.Type, fill = factor(Phase)))
p <- p + geom_bar(position = position_dodge2(preserve = "single"))
p <- p + theme_classic() + labs(x = 'Cell Type', y = 'Cell Count')
p <- p + theme(axis.text.x=element_text(angle=30, hjust=1))
p <- p + guides(fill = guide_legend(title ='Phase'))
ggsave(paste(opt$out, 'cell_cycle_summary.pdf', sep=''), plot=p, width=8, height=7, units = "in")
ggsave(paste(opt$out, 'cell_cycle_summary.png', sep=''), plot=p, width=8, height=7, dpi=600, units = "in")
pda <- scd[otop10$gene, rownames(anno_col)[order(anno_col$Cell.Type)]]
pheatmap(pda, annotation_col = anno_col, show_colnames=F, cluster_cols=F, cluster_rows=F, fontsize_row=nrow(pda)/50, color = colorpanel(128, 'magenta2', 'black', 'yellow'), filename=paste(opt$out, '/top10_hetmap_cell_cycle.pdf', sep=''), width=(nrow(pda)/20)-1, height=nrow(pda)/20, fontsize=5)
pheatmap(pda, annotation_col = anno_col, show_colnames=F, cluster_cols=F, cluster_rows=F, fontsize_row=nrow(pda)/50, color = colorpanel(128, 'magenta2', 'black', 'yellow'), filename=paste(opt$out, '/top10_hetmap_cell_cycle.png', sep=''), width=(nrow(pda)/20)-1, height=nrow(pda)/20, fontsize=5)
cc_out <- data.frame(Cell_id=rownames(anno_col), anno_col)
write.table(cc_out, file=paste(opt$out, 'cell_cycle.txt', sep=''), quote=F, row.names=F, sep='\t')

saveRDS(sub2, file=paste(opt$out, 'cell_cycle.rds', sep=''))
