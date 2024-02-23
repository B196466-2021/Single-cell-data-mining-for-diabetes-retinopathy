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
updated date: 2021-10-21\n\n")
	q(status=1)
}

usage<-function(spec){
	cat("This script use AUCell to identify cells with an active ‘gene set’ (i.e. gene signatures) in single-cell RNA-seq data.\n",
getopt(spec,usage=TRUE),
"Options:
	-r, --rds	rds file result from scRNA_cell_typing.R.
	-g, --gset	gmt file of gene set.
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
	'rds','r',1,'character',
	'gset','g',1,'character',
	'thread','t',1,'integer',
	'count','c',1,'integer',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$rds) || is.null(opt$gset) || !is.null(opt$help)) {usage(spec)}

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
library(AUCell)
library(Seurat)
library(GSEABase)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(dplyr)
source('/home/dev/DEV-wangdh/R/priv_auc.assignmnetThreshold_v6.R')

rds <- readRDS(opt$rds)
exprMatrix <- rds@assays$RNA@counts
geneSets <- getGmt(opt$gset)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=opt$thread, plotStats=F)
expm <- try(as.matrix(exprMatrix), silent=T)
if(!'try-error' %in% class(expm)){
	png(paste(opt$out, 'genes_cells_stat.png', sep=''), width=7, height=7, units="in", res=600)
	plotGeneCount(expm)
	dev.off()
	pdf(paste(opt$out, 'genes_cells_stat.pdf', sep=''), width=7, height=7)
	plotGeneCount(expm)
	dev.off()
}
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, nCores=opt$thread, aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE, nCores=10, thrP=0.05)

#plot the AUC histograms (per gene-set)
aucMatrix <- getAUC(cells_AUC)
if(!is.matrix(aucMatrix)) stop("cellsAUC should contain the AUC values.")
rowSumAUC <- rowSums(aucMatrix)
if(any(rowSumAUC==0)){
	warning("Skipping genesets with all AUC 0: ", paste(names(rowSumAUC)[which(rowSumAUC==0)], collapse=", "), immediate. = TRUE)
}
aucMatrix <- aucMatrix[rowSumAUC>0,,drop=FALSE]
dauc <- data.frame(cell=colnames(aucMatrix), t(aucMatrix))
write.table(dauc, paste(opt$out, '/auc_score.txt', sep=''), sep='\t', quote=F, row.names=F)
odir <- paste(opt$out, '/auc_histograms/', sep='')
if(! dir.exists(odir)){
	dir.create(odir, recursive=T)
}
for(n in rownames(aucMatrix)){
	png(paste(odir, gsub('\\(|\\)', '', gsub(' |/', '_', n)), '_histplot.png', sep=''), width=7, height=7, units='in', res=600)
	aucThr <- .auc_assignmnetThreshold_v6(aucRow=aucMatrix[n,, drop=FALSE], thrP=0.01, smallestPopPercent=.25, plotHist=TRUE, densAdjust=2, nBreaks=100)
	dev.off()
	pdf(paste(odir, gsub('\\(|\\)', '', gsub(' |/', '_', n)), '_histplot.pdf', sep=''), width=7, height=7)
	aucThr <- .auc_assignmnetThreshold_v6(aucRow=aucMatrix[n,, drop=FALSE], thrP=0.01, smallestPopPercent=.25, plotHist=TRUE, densAdjust=2, nBreaks=100)
	dev.off()
}

cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
mas <- merge(rds@meta.data, assignmentTable, by.x=0, by.y=1)
omas <- mas[,names(mas) %in% c('Row.names', 'sample', 'group', 'seurat_clusters', 'cell_type', 'geneSet')]
names(omas)[1] <- 'Cell'
write.table(omas, paste(opt$out, 'cell_assignment.txt', sep=''), sep='\t', quote=F, row.names=F)

#plot the AUC umap (per gene-set)
rds@assays$RNA@data <- aucMatrix
odir <- paste(opt$out, '/auc_unap/', sep='')
if(! dir.exists(odir)){
	dir.create(odir, recursive=T)
}
odir2 <- paste(opt$out, '/auc_violin/', sep='')
if(! dir.exists(odir2)){
	dir.create(odir2, recursive=T)
}
odir3 <- paste(opt$out, '/diff_boxplot/', sep='')
if(! dir.exists(odir3)){
	dir.create(odir3, recursive=T)
}
for(n in rownames(aucMatrix)){
	p <- FeaturePlot(rds, features=n, slot = "data", cols=c('darkorchid4', 'red', 'khaki1'))
	ggsave(paste(odir, gsub('\\(|\\)', '', gsub(' |/', '_', n)), '_umap.png', sep=''), plot=p, width=7, height=7, units='in', dpi=600)
	ggsave(paste(odir, gsub('\\(|\\)', '', gsub(' |/', '_', n)), '_umap.pdf', sep=''), plot=p, width=7, height=7, units='in')
	p <- VlnPlot(rds, features=n, slot = "data",  pt.size=0)
	p <- p + geom_boxplot(outlier.shape=NA, fill='white', width=0.5) + scale_fill_igv() + labs(x=NULL, y='AUC score')
	p <- p + theme(legend.position='none', axis.line = element_line(arrow = arrow(angle=20, length = unit(0.3, 'cm'), type = "closed")))
	ggsave(paste(odir2, gsub('\\(|\\)', '', gsub(' |/', '_', n)), '_violin.png', sep=''), plot=p, width=7, height=7, units='in', dpi=600)
	ggsave(paste(odir2, gsub('\\(|\\)', '', gsub(' |/', '_', n)), '_violin.pdf', sep=''), plot=p, width=7, height=7, units='in')
	if(all(c('group', 'cell_type') %in% names(rds@meta.data))){
		pd <- data.frame(aus=t(aucMatrix)[, n], rds@meta.data[colnames(aucMatrix), c('group', 'cell_type')])
		pd$cell_type <- factor(pd$cell_type, levels=unique(pd$cell_type)[order(nchar(unique(as.character(pd$cell_type))))])
		if(length(unique(rds@meta.data$group))==2){
			pv <- pd %>% group_by(cell_type) %>% filter(length(unique(group)) > 1) %>% group_by(cell_type) %>% summarise(p=wilcox.test(formula=aus~group)$p.value) %>% data.frame()
		}else if(length(unique(rds@meta.data$group)) > 2){
			pv <- pd %>% group_by(cell_type) %>% filter(length(unique(group)) > 1) %>% group_by(cell_type) %>% summarise(p=kruskal.test(formula=aus~group)$p.value) %>% data.frame()
		}
		pv$labe <- ''
		pv$labe[pv$p < 0.001] <- '***'
		pv$labe[pv$p < 0.01 & pv$p >=0.001] <- '**'
		pv$labe[pv$p < 0.05 & pv$p >=0.01] <- '*'
		labe.caption <- expression(paste('*: ', 0.01 <= p, '< 0.05;  **: ', 0.001 <= p, '< 0.01;  ***: p < 0.001'))
		write.table(pv, paste(odir3, gsub('\\(|\\)', '', gsub(' |/', '_', n)), '_pvalue.txt', sep=''), sep='\t', quote=F, row.names=F)
		y <- max(pd$aus) + (max(pd$aus) - min(pd$aus))*0.02
		p <- ggplot(pd, aes(x=cell_type, y=aus, fill=group)) + geom_boxplot()
		p <- p + geom_text(data=pv, aes(x=cell_type, y=y, label=labe), inherit.aes=F, size=5)
		p <- p + labs(x= NULL, y='AUC score', caption=labe.caption) + theme_bw()
		p <- p + theme(axis.text=element_text(size=12, color='black'), axis.title=element_text(size=15), legend.text=element_text(size=12), 
					   legend.title=element_blank(), axis.text.x=element_text(angle=45, hjust=1))
		w <- ifelse(length(unique(pd$cell_type)) > 6, length(unique(pd$cell_type))/2 + 2, 5)
		ggsave(paste(odir3, gsub('\\(|\\)', '', gsub(' |/', '_', n)), '_boxplot.png', sep=''), plot=p, width=w, height=5, units='in', dpi=600)
		ggsave(paste(odir3, gsub('\\(|\\)', '', gsub(' |/', '_', n)), '_boxplot.pdf', sep=''), plot=p, width=w, height=5, units='in')
	}
}

#Explore cells/clusters based on the signature score
selectedThresholds <- getThresholdSelected(cells_assignment)
odir <- paste(opt$out, '/auc_classify/', sep='')
if(! dir.exists(odir)){
	dir.create(odir, recursive=T)
}
setr <- data.frame(Gene_Set=names(selectedThresholds), selectedThresholds=selectedThresholds)
write.table(setr, paste(odir, 'selectedThresholds.txt', sep=''), sep='\t', quote=F, row.names=F)
nBreaks <- 5
colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
pd <- data.frame(x=seq(1,5), y=seq(1,5), s=seq(1,5))
pl <- ggplot(pd, aes(x=x, y=y, fill=s)) + geom_point() + theme(legend.text=element_blank(), legend.title=element_text(angle=90, hjust=0.5), legend.key.height=unit(8, 'mm'))
ln <- pl + scale_fill_gradientn(colours=colorPal_Neg , guide=guide_colourbar(ticks=FALSE, title='Not pass the threshold', title.position='right'))
lp <- pl + scale_fill_gradientn(colours=colorPal_Pos, guide=guide_colourbar(ticks=FALSE, title='Pass the threshold', title.position='right'))
pb <- ggplot() + scale_x_continuous(limits=c(0,1), expand=c(0, 0)) + scale_y_continuous(limits=c(0,1), expand=c(0, 0))
pb <- pb + theme(axis.text=element_blank(), axis.ticks=element_blank(), plot.margin=margin(0,0,0,0), panel.background=element_blank())
for(geneSetName in names(selectedThresholds)){
	passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
	if(sum(passThreshold) >0){
		aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
		cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
					   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
		rds@meta.data$aucol <- cellColor[rownames(rds@meta.data)]
		p <- DimPlot(rds, cols=rds@meta.data$aucol, reduction='umap', group.by='aucol') + labs(title=gsub('_', ' ',  geneSetName), x='UMAP 1', y='UMAP 2') + theme(legend.position='none')
		p <- pb + annotation_custom(ggplotGrob(p),xmin=0,xmax=0.9,ymin=0,ymax=1)
		p <- p + annotation_custom(ggplotGrob(as_ggplot(get_legend(lp))),xmin=0.9,xmax=1,ymin=0.6,ymax=0.8)
		p <- p + annotation_custom(ggplotGrob(as_ggplot(get_legend(ln))),xmin=0.9,xmax=1,ymin=0.2,ymax=0.5)
		ggsave(paste(odir, gsub('\\(|\\)', '', gsub(' |/', '_', geneSetName)), '_umap.png', sep=''), plot=p, width=7, height=7, units='in', dpi=600)
		ggsave(paste(odir, gsub('\\(|\\)', '', gsub(' |/', '_', geneSetName)), '_umap.pdf', sep=''), plot=p, width=7, height=7, units='in')
	}
}

#detected cells
if('cell_type' %in% names(rds@meta.data)){
	cnl <- 'cell_type'
}else{
	cnl <- 'seurat_clusters'
}
confMatrix <- t(sapply(cells_assignment, function(x) table(rds@meta.data[x$assignment, cnl])[unique(rds@meta.data[,cnl])]))
if(length(unique(rds@meta.data[,cnl]))>1){
	colnames(confMatrix) <- unique(rds@meta.data[,cnl])
	confMatrix[which(is.na(confMatrix), arr.ind=TRUE)] <- 0
	od <- data.frame(Gene_set=rownames(confMatrix), confMatrix)
	write.table(od, paste(opt$out, 'number_of_detected_cells.txt', sep=''), sep='\t', quote=F, row.names=F)
}else if(length(unique(rds@meta.data[,cnl]))==1){
	od <- data.frame(gsub(paste('.', unique(rds@meta.data[,cnl]), '$', sep=''), '', colnames(confMatrix)), as.numeric(confMatrix))
	names(od) <- c('Gene_set', unique(rds@meta.data[,cnl]))
	write.table(od, paste(opt$out, 'number_of_detected_cells.txt', sep=''), sep='\t', quote=F, row.names=F)
}


