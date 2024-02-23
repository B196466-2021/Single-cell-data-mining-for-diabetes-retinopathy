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
	cat("This script is used to cell-typing cluster of scRNA.\n",
getopt(spec,usage=TRUE),
"Options:
	-f, --file	rds file resulted from scRNA cluster.
	-x, --xls	xlsx file of cell type enrichment results.
	-t, --thread	thread number (default 1).
	-o, --out	out directory (default cell_type/).
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
	'xls','x',1,'character',
	'thread','t',1,'integer',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || is.null(opt$xls) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$out)) {opt$out='cell_type/'}
if (is.null(opt$thread)) {opt$thread=1}

dir <- opt$out
if(! dir.exists(dir)){
	dir.create(dir, recursive=T)
}
##########################
#define functions
fdmk <- function(x, pbmc){
	markers <- FindMarkers(pbmc, ident.1 = x, verbose = FALSE, test.use='MAST', only.pos = T, min.pct = 0.25, logfc.threshold = 0.1)
	markers$cell_type <- x
	markers$gene <- rownames(markers)
	markers <- markers[order(markers$avg_log2FC, decreasing=T),]
	return(markers)
}

##########################
#main programe
library(Seurat)
library(cowplot)
library(ggplot2)
library(parallel)
library(dplyr)
library(foreach)
library(pheatmap)
library(gplots)
library(openxlsx)
library(stringr)
library(scales)
library(ggsci)

makeh <- read.xlsx(opt$xls, sheet = 1)
makeh$rn <- as.numeric(rownames(makeh))
sub <- readRDS(opt$file)

#cell typing of cluster
clt <- makeh %>% group_by(Cluster) %>% top_n(n = 1, wt = -rn)
oct <- clt[,c('Cluster', 'Cell_Name', 'Genes')]
write.table(oct, file=paste(opt$out, '/cell_marker_selected.txt', sep=''), row.names=F, quote=F, sep='\t')
scg <- apply(oct, 1, function(x){data.frame(clustere=x[1], cell_type=x[2], gene=strsplit(x[3], '//|,|;', perl=T)[[1]], stringsAsFactors = F)})
dd <- do.call(rbind, scg)
scg <- unique(dd[,c('cell_type',  'gene')])
cell_type <- as.character(clt$Cell_Name)
names(cell_type) <- clt$Cluster
if(any(! levels(sub) %in% names(cell_type))){
	nct <- as.character(levels(sub)[! levels(sub) %in% names(cell_type)])
	cnct <- paste('Cluster_', nct, sep='')
	names(cnct) <- nct
	cell_type <- c(cell_type, cnct)
}

sub2 <- RenameIdents(sub, cell_type)
sub2@meta.data$cell_type <- cell_type[as.character(sub2@meta.data$seurat_clusters)]
p <- DimPlot(sub2, reduction='tsne', label = TRUE, label.size=3)
p2 <- DimPlot(sub2, reduction='umap', label = TRUE, label.size=3)
p <- p + labs(x='tSNE1', y='tSNE2') + scale_color_igv()
p <- p + theme(axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
p2 <- p2 + labs(x='UMAP1', y='UMAP2') + scale_color_igv()
p2 <- p2 + theme(axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
w <- 6 + max(nchar(cell_type))/10
ggsave(paste(opt$out, '/tSNE_cell_type.pdf', sep=''), plot=p, height=6, width=w, units='in')
ggsave(paste(opt$out, '/tSNE_cell_type.png', sep=''), plot=p, height=6, width=w, dpi=600, units='in')

ggsave(paste(opt$out, '/UMAP_cell_type.pdf', sep=''), plot=p2, height=6, width=w, units='in')
ggsave(paste(opt$out, '/UMAP_cell_type.png', sep=''), plot=p2, height=6, width=w, dpi=600, units='in')

if('group' %in% names(sub2@meta.data) && length(unique(sub2@meta.data$group)) > 1){
	w <- 6*length(unique(sub2@meta.data$group)) + 0.5
	p <- DimPlot(sub2, reduction = "tsne", label=T, split.by='group', label.size=3) + scale_color_igv()
	p <- p + theme(plot.title=element_blank(), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
	ggsave(paste(opt$out, '/TSNE_split_group.pdf', sep=''), p, height=6, width=w, units='in')
	ggsave(paste(opt$out, '/TSNE_split_group.png', sep=''), p, height=6, width=w, dpi=600, units='in')
	p <- DimPlot(sub2, reduction = "umap", label=T, split.by='group', label.size=3) + scale_color_igv()
	p <- p + theme(plot.title=element_blank(), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
	ggsave(paste(opt$out, '/UMAP_split_group.pdf', sep=''), p, height=6, width=w, units='in')
	ggsave(paste(opt$out, '/UMAP_split_group.png', sep=''), p, height=6, width=w, dpi=600, units='in')
}
#new marker gene
mc <- getOption("mc.cores", opt$thread)
maker.list <- mclapply(levels(sub2@active.ident), fdmk, pbmc=sub2, mc.cores = mc)
maker_gene <- do.call(rbind, maker.list)
write.table(maker_gene, file=paste(opt$out, '/new_marker_gene_of_each_cell_type.txt', sep=''), quote=F, row.names=F, sep='\t')
#marker gene heatmap
top10 <- maker_gene %>% group_by(cell_type) %>% top_n(n = 10, wt = avg_log2FC)
top10$cell_type <- factor(top10$cell_type, levels=levels(sub2@active.ident))
otop10 <- top10[order(as.numeric(top10$cell_type)), ]
xlsz <- 500/length(otop10$gene)
default.colors <- c(hue_pal()(length(levels(sub2@active.ident))))
names(default.colors) <- levels(sub2@active.ident)
p <- DoHeatmap(sub2, features = as.character(otop10$gene), label=F)
p <- p + scale_colour_manual(values =default.colors, name='Cell Type')
p <- p + theme(axis.text.y=element_text(size=xlsz, color='black'))
ncl <- length(levels(sub2@active.ident))
if(ncl <=7){
	w <- 7
    h <- 7
}else{
	w <- ncl
    h <- ncl
}
ggsave(paste(opt$out, '/top10_new_marker_gene_heatmap.pdf', sep=''), plot=p, height=h, width=w, units='in')
ggsave(paste(opt$out, '/top10_new_marker_gene_heatmap.png', sep=''), plot=p, height=h, width=w, dpi=300, units='in')

#expression of marker gene
odir <- paste(opt$out, '/marker_gene_exp', sep='')
if(!dir.exists(odir)){
	dir.create(odir, recursive=T)
}
for(cl in unique(scg$cell_type)){
	genes <- unique(scg[scg$cell_type==cl, 'gene'])
	genes <- genes[!is.na(genes)]
	if(length(genes) >0 ){
		ncs <- ceiling(sqrt(length(genes)))
		nrs <- ceiling(length(genes)/ncs)
		w <- ncs*3.5
		h <- nrs*3
		p <- FeaturePlot(sub2, features = as.character(genes), ncol=ncs, reduction='umap')
		ggsave(paste(odir, '/', gsub(' ', '_', cl), '_marker_gene_exp_UMAP.pdf', sep=''), plot=p, width=w, height=h, device="pdf")
		ggsave(paste(odir, '/', gsub(' ', '_', cl), '_marker_gene_exp_UMAP.png', sep=''), plot=p, height=h, width=w, dpi=600, units='in', device="png")
		p2 <- FeaturePlot(sub2, features = as.character(genes), ncol=ncs, reduction='tsne')
		ggsave(paste(odir, '/', gsub(' ', '_', cl), '_marker_gene_exp_tSNE.pdf', sep=''), plot=p2, width=w, height=h, device="pdf")
		ggsave(paste(odir, '/', gsub(' ', '_', cl), '_marker_gene_exp_tSNE.png', sep=''), plot=p2, height=h, width=w, dpi=600, units='in', device="png")
	}
}
cell.stat <- data.frame(table(sub2@active.ident))
names(cell.stat) <- c('Cell_type', 'Count')
write.table(cell.stat, file=paste(opt$out, '/cell_stat.txt', sep=''), quote=F, sep='\t', row.names=F)
sub2@misc$otop10 <- otop10
#marker genes violin plot
gene <- unique(dd$gene)
gene <- gene[!is.na(gene)]
if(length(gene) >0){
	exp <- as.matrix(sub2@assays$RNA@data[gene,])
	pd <- data.frame(gene=rep(rownames(exp), ncol(exp)), cluster=rep(sub2@meta.data[colnames(exp),'seurat_clusters'], each=nrow(exp)), 
					 cell=rep(sub2@meta.data[colnames(exp),'cell_type'], each=nrow(exp)), y=as.numeric(exp))
	h <- ifelse(length(gene) > 15, length(gene)/3, 5)
	##by cluster
	pd$gene <- factor(pd$gene, levels=gene)
	p <- ggplot(pd, aes(x=cluster, y=y, fill=cluster, color=cluster)) + geom_violin(show.legend=F, scale='width') + scale_color_igv() + scale_fill_igv()
	p <- p + facet_grid(rows = vars(gene), switch='y', scales = "free") + labs(x='Cluster', y=NULL) + theme_bw()
	p <- p+ theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text.y.left = element_text(angle = 0, size=15, hjust=1), 
				  panel.spacing.y=unit(0, 'mm'), panel.grid=element_blank(), strip.background.y=element_rect(color='white', fill='white'), 
				  axis.text=element_text(size=15, color='black'), axis.title=element_text(size=20))
	w <- ifelse(length(unique(pd$cluster)) > 7, length(unique(pd$cluster))/2 + 1.5, 5)
	ggsave(paste(odir, '/marker_gene_cluster_violin_plot.png', sep=''), plot=p, width=w, height=h, units='in', dpi=600)
	ggsave(paste(odir, '/marker_gene_cluster_violin_plot.pdf', sep=''), plot=p, width=w, height=h, units='in')
}
##by cell type
dc <- scg[order(scg$cell_type),]
gene <- unique(dc$gene)
gene <- gene[!is.na(gene)]
if(length(gene) >0){
	pd$gene <- factor(pd$gene, levels=gene)
	p <- ggplot(pd, aes(x=cell, y=y, fill=cell, color=cell)) + geom_violin(show.legend=F, scale='width') + scale_color_igv() + scale_fill_igv()
	p <- p + facet_grid(rows = vars(gene), switch='y', scales = "free") + labs(x=NULL, y=NULL) + theme_bw()
	p <- p+ theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text.y.left=element_text(angle=0, size=15, hjust=1), 
				  panel.spacing.y=unit(0, 'mm'), panel.grid=element_blank(), strip.background.y=element_rect(color='white', fill='white'), 
				  axis.text=element_text(size=15, color='black', angle=45, hjust=1), axis.title=element_text(size=20))
	w <- ifelse(length(unique(pd$cell)) > 7, length(unique(pd$cell))/2 + 1.5, 5)
	ggsave(paste(odir, '/marker_gene_cell_type_violin_plot.png', sep=''), plot=p, width=w, height=h+1, units='in', dpi=600)
	ggsave(paste(odir, '/marker_gene_cell_type_violin_plot.pdf', sep=''), plot=p, width=w, height=h+1, units='in')
}
#marker genes average expreesion dot plot
scgo <- scg[order(scg$cell_type),]
scgo <- scgo[!is.na(scgo$gene), ]
if(nrow(scgo) >0){
	p <- DotPlot(sub2, features = unique(scgo$gene), group.by='cell_type') + xlab('Genes') + RotatedAxis()
	p <- p + guides(color = guide_colourbar(title='Average Expression', title.position='top', direction='horizontal', title.hjust=0.5, barwidth=unit(2, 'in')), 
					size=guide_legend(title = "Percent Expressed", title.position='top', direction='horizontal', title.hjust=0.5)) 
	p <- p + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='top', legend.justification="center", 
				   axis.line=element_blank(), panel.border=element_rect(color='black', size=1), axis.title.y=element_blank())
	if(length(scgo$gene) < 16){
		w <- 8
	}else{
		w <- length(scgo$gene)*0.4
	}
	if(length(unique(sub2@active.ident)) < 6){
		h <- 5
	}else{
		h <- 1.5 + length(unique(sub2@active.ident))*0.7
	}

	ggsave(paste(opt$out, '/used_marker_gene_dotplot.png', sep=''), plot=p, width=w, height=h, dpi=600, units='in')
	ggsave(paste(opt$out, '/used_marker_gene_dotplot.pdf', sep=''), plot=p, width=w, height=h, units='in')
}
#cell ratio
if('sample' %in% names(sub2@meta.data)){
	cst <- sub2@meta.data %>% group_by(sample, cell_type) %>% summarise(count=length(sample)) %>% group_by(sample) %>% mutate(total=sum(count), Ratio=count/total) %>% data.frame()
	if('group' %in% names(sub2@meta.data) && length(unique(sub2@meta.data$group)) >1){
		gsd <- unique(sub2@meta.data[,c('group', 'sample')])
		gsd <- gsd[order(gsd$group, gsd$sample),]
		cst$sample <- factor(cst$sample, levels=gsd$sample)
		sclor <- pal_npg()(length(unique(gsd$group)))[as.numeric(factor(gsd$group))]
	}else{
		sclor <- 'black'
	}
	p <- ggplot(cst, aes(x=sample, y=Ratio, fill=cell_type)) + geom_col()
	p <- p + scale_y_continuous(expand=c(0, 0)) + scale_fill_igv() + xlab(NULL) + theme_bw()
	p <- p + theme(axis.title=element_text(size=20), axis.text=element_text(size=12, color='black'), axis.text.x=element_text(angle=90, vjust=0.5, color=sclor),
				   axis.ticks.x=element_blank(), legend.title=element_blank(), legend.text=element_text(size=15))
	w2 <- ifelse(length(unique(cst$sample))>8, length(unique(cst$sample))/4 + max(nchar(cst$cell_type))/10 + 1, max(nchar(cst$cell_type))/10 + 3)
	ggsave(paste(opt$out, '/cell_component_in_sample.png', sep=''), plot=p, width=w2, height=5, dpi=600, units='in')
	ggsave(paste(opt$out, '/cell_component_in_sample.pdf', sep=''), plot=p, width=w2, height=5, units='in')
}

if('group' %in% names(sub2@meta.data) && length(unique(sub2@meta.data$group)) >1){
	cst <- sub2@meta.data %>% group_by(cell_type, group) %>% summarise(count=length(group)) %>% group_by(group) %>% mutate(total=sum(count), ratio=count/total) %>% data.frame()
	write.table(cst, paste(opt$out, '/cell_ratio_by_cell_type.txt', sep=''), sep='\t', quote=F, row.names=F)
	p <- ggplot(cst, aes(x=cell_type, y=ratio, fill=group, color=group)) + geom_col(position = "dodge2")
	p <- p + labs(x=NULL, y='Cell ratio') + scale_fill_npg() + theme_bw()
	p <- p + theme(axis.text=element_text(size=12, color='black'), axis.title=element_text(size=15), axis.text.x=element_text(angle=45, hjust=1), 
				   legend.title=element_blank(), legend.text=element_text(size=12), legend.position='top')
	w <- ifelse(length(unique(cst$cell_type)) > 6, 2 + length(unique(cst$cell_type))/2, 5)
	ggsave(paste(opt$out, '/cell_ratio_by_cell_type.png', sep=''), plot=p, width=w, height=5, dpi=600, units='in')
	ggsave(paste(opt$out, '/cell_ratio_by_cell_type.pdf', sep=''), plot=p, width=w, height=5, units='in')

	#by sample
	if(any(table(unique(sub2@meta.data[, c('group', 'sample')])[, 'group']) > 1)){
		cstat <- sub2@meta.data %>% group_by(cell_type, group, sample) %>% summarise(count=length(sample)) %>% group_by(sample) %>% mutate(total=sum(count), Ratio=count/total) %>% data.frame(stringsAsFactors=F)
		if(length(unique(sub2@meta.data$group)) > 2){
			pd <- cstat %>% group_by(cell_type) %>% filter(length(unique(group)) > 1) %>% group_by(cell_type) %>% summarise(p=kruskal.test(formula=Ratio~group)$p.value) %>% data.frame(stringsAsFactors=F)
			pdo <- cstat %>% group_by(cell_type) %>% filter(length(unique(group)) > 1) %>% group_by(cell_type) %>% mutate(p_value=kruskal.test(formula=Ratio~group)$p.value) %>% data.frame(stringsAsFactors=F)
		}else{
			pd <- cstat %>% group_by(cell_type) %>% filter(length(unique(group)) > 1) %>% group_by(cell_type) %>% summarise(p=wilcox.test(formula=Ratio~group)$p.value) %>% data.frame(stringsAsFactors=F)
			pdo <- cstat %>% group_by(cell_type) %>% filter(length(unique(group)) > 1) %>% group_by(cell_type) %>% mutate(p_value=wilcox.test(formula=Ratio~group)$p.value) %>% data.frame(stringsAsFactors=F)
		}
		write.table(pdo, paste(opt$out, '/cell_ratio_by_sample.txt', sep=''), sep='\t', quote=F, row.names=F)
		pd$labe <- ''
		pd$labe[pd$p < 0.001] <- '***'
		pd$labe[pd$p < 0.01 & pd$p >=0.001] <- '**'
		pd$labe[pd$p < 0.05 & pd$p >=0.01] <- '*'
		labe.caption <- expression(paste('*: ', 0.01 <= p, '< 0.05;  **: ', 0.001 <= p, '< 0.01;  ***: p < 0.001'))
		p <- ggplot(cstat, aes(x=cell_type, y=Ratio, fill=group)) + geom_boxplot(outlier.shape=NA)
		gp.stat <- ggplot_build(p)$data[[1]]
		pd.stat <- gp.stat[, c('x', 'ymin', 'ymax')]
		p <- p + geom_text(data=pd, aes(x=cell_type, y=max(pd.stat$ymax)+(max(pd.stat$ymax)-min(pd.stat$ymin))*0.01, label=labe), inherit.aes=F, size=7)
		p <- p + coord_cartesian(ylim=c(min(pd.stat$ymin), max(pd.stat$ymax) + (max(pd.stat$ymax)-min(pd.stat$ymin))*0.03)) + scale_fill_npg() 
		p <- p + labs(caption=labe.caption) + theme_bw()
		p <- p + theme(axis.title=element_text(size=20), axis.title.x=element_blank(), axis.text=element_text(size=15, color='black'), 
					   axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank(), legend.text=element_text(size=15), legend.position='top')
		ggsave(paste(opt$out, '/cell_ratio_by_sample.png', sep=''), plot=p, width=w, height=5.5, dpi=600, units='in')
		ggsave(paste(opt$out, '/cell_ratio_by_sample.pdf', sep=''), plot=p, width=w, height=5.5, units='in')
	}
}
write.table(data.frame(Cell=rownames(sub2@meta.data), sub2@meta.data, stringsAsFactors=F), paste(opt$out, '/cell_meta_data.txt', sep=''), quote=F, sep='\t', row.names=F)
sub2@assays[['group.order']] <- levels(sub2@meta.data$group)
saveRDS(sub2, file=paste(opt$out, '/cell_typing.rds', sep=''))

