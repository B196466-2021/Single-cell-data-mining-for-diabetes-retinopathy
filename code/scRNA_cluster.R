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
	cat("This script is used for single cell RNA-seq analysis.\n",
getopt(spec,usage=TRUE),
"Options:
	-f, --file      rds file(s file1,file2,...) resulted from scRNA_QC.R.
	-t, --thread	number of threads (default 1).
	-s, --tissue	tissue in CellMarker file that samples come from (http://biocc.hrbmu.edu.cn/CellMarker/index.jsp#) default NULL.
	-l, --mcl	maximum clusters (default 15).
	-i, --integrate	integrate datasets of all samples or not (default TRUE).
	-c, --ct	cancer cell type while specified. Default: Normal.
	-o, --out	out directory (default cluster/).
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
	'tissue','s',1,'character',
	'integrate','i',2,'logical',
	'ct','c',2,'logical',
	'mcl','l',1,'integer',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$thread)) {opt$thread=1}
if (is.null(opt$out)) {opt$out='cluster/'}
if (is.null(opt$ct)) {opt$ct=F}
if (is.null(opt$mcl)) {opt$mcl=15}
if (is.null(opt$integrate)) {opt$integrate=T}

dir <- gsub("[^/]*$", '', opt$out, perl=T)
if(! dir.exists(dir) && dir != ''){
	dir.create(dir, recursive=T)
}

##########################
#define functions

fdmk <- function(x, pbmc){
	library(dplyr)
	library(Seurat)
	markers <- FindMarkers(pbmc, ident.1 = x, verbose = FALSE, test.use='MAST', only.pos = T, min.pct = 0.25, logfc.threshold = 0.1)
	markers$cluster <- x
	markers$gene <- rownames(markers)
	markers <- markers %>% dplyr::filter(p_val<0.05) %>% top_n(500, wt=avg_log2FC)
	return(markers)
}

depf <- function(x){
	diff_gene <- FindMarkers(pbmc2, ident.1 = x, verbose = FALSE, test.use='MAST', min.pct = 0.25, logfc.threshold = 0)
	od <- data.frame(Gene_Name=rownames(diff_gene), log2FC=diff_gene$avg_log2FC, FC=2**diff_gene$avg_log2FC, p_value=diff_gene$p_val, q_value=diff_gene$p_val_adj, pct.1=diff_gene$pct.1, pct.2=diff_gene$pct.2)
	write.table(od, paste(g_dir, '/', gsub(' ', '_', x), '_vs_other_cell.txt', sep=''), quote=F, row.names=F, sep='\t')
}

##main program
library(DropletUtils)
library(biomaRt)
library(dplyr)
library(scater)
library(ggplot2)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(scran)
library(pheatmap)
library(gplots)
library(parallel)
library(inflection)
library(scales)
library(ggsci)
library(cowplot)
library(harmony)

rds.list <- list()
files <- strsplit(opt$file, ',')[[1]]
for(f in files){
	rds <- readRDS(f)
	if(is.list(rds)){
		for(n in 1:length(rds)){
			if(class(rds[[n]]) == "Seurat"){
				rds.list[names(rds[n])] <- rds[[n]]
			}
		}
	}else if(class(rds) == "Seurat"){
		if(opt$integrate==T && ! 'integrated' %in% names(rds)){
			rds.list <- c(rds.list, SplitObject(rds, split.by='sample'))
		}else{
			rds.list[gsub('.rds$', '', basename(f))] <- rds
		}
	}
}

reduc <- 'pca'
if(length(rds.list) >1){
	if(opt$integrate==T){
		rds.list2 <- try(lapply(X = rds.list, FUN = function(x){x=NormalizeData(x, assay='RNA', verbose = FALSE); x=FindVariableFeatures(x, assay='RNA', verbose = FALSE, nfeatures=5000)}), silent=T)
		if(!'try-error' %in% class(rds.list2)){
			features <- SelectIntegrationFeatures(object.list = rds.list2, assay=rep('RNA', length(rds.list2)))
			rds.list2 <- try(lapply(X = rds.list2, FUN = function(x) {x <- ScaleData(x, features = features, verbose = FALSE); x <- RunPCA(x, features = features, verbose = FALSE, npcs=50)}), silent=T)
		}
		if(!'try-error' %in% class(rds.list2)){
			anchors <- FindIntegrationAnchors(object.list = rds.list2, reduction = "rpca", dims = 1:50)
			pbmc <- try(IntegrateData(anchorset = anchors, dims = 1:50, k.weight=100), silent=TRUE)
			if('try-error' %in% class(pbmc)){
				pbmc <- try(IntegrateData(anchorset = anchors, dims = 1:50, k.weight=50), silent=TRUE)
				if('try-error' %in% class(pbmc)){
					pbmc <- try(IntegrateData(anchorset = anchors, dims = 1:50, k.weight=20), silent=TRUE)
				}
			}
			DefaultAssay(pbmc) <- 'integrated'
			pbmc <- ScaleData(pbmc, verbose=FALSE)
			pbmc <- RunPCA(pbmc, verbose=F)
		}else{
			pbmc <- merge(x=rds.list[[1]], y=rds.list[-1])
			pbmc <- NormalizeData(pbmc, verbose = FALSE)
			pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
			pbmc <- ScaleData(pbmc, features =rownames(pbmc))
			pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F)
			pbmc <- RunHarmony(pbmc, 'sample', plot_convergence=F)
			reduc <- 'harmony'
		}
	}else{
		pbmc <- merge(x=rds.list[[1]], y=rds.list[-1])
		pbmc <- NormalizeData(pbmc, verbose = FALSE)
		pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
		pbmc <- ScaleData(pbmc, features = rownames(pbmc))
		pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F)
	}
}else{
	if(opt$integrate==T && 'integrated' %in% names(rds.list[[1]])){
		pbmc <- rds.list[[1]]
		DefaultAssay(pbmc) <- 'integrated'
	}else{
		pbmc <- NormalizeData(rds.list[[1]], verbose = FALSE)
		pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
	}
	pbmc <- ScaleData(pbmc, features = rownames(pbmc))
	pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F)
}

g_dir <- opt$out

##feature selection
#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#top10 <- head(VariableFeatures(pbmc), 10)
#plot1 <- VariableFeaturePlot(pbmc)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot2 <- plot2 + theme(legend.position='top', legend.justification="center", axis.line=element_blank(), 
#					   panel.border=element_rect(color='black'), plot.margin=margin(t=5, b=5, l=5, r=20, unit='pt'))
#ggsave(paste(g_dir, '/feature_select.pdf', sep=''), plot2,  height=5, width=5, units='in')
#ggsave(paste(g_dir, '/feature_select.png', sep=''), plot=plot2, height=5, width=5, dpi=600, units='in')

##Scaling the data
if('group.order' %in% names(rds)){
	pbmc@meta.data$group <- factor(pbmc@meta.data$group, levels=rds$group.order)
}
##Dimension reduction
###PCA
pl <- VizDimLoadings(pbmc, dims = 1:2, reduction = reduc, combine=F)
pl[[1]] <- pl[[1]] +  scale_x_continuous(n.breaks=4)
pl[[2]] <- pl[[2]] +  scale_x_continuous(n.breaks=4)
p <- plot_grid(plotlist=pl)
ggsave(paste(g_dir, '/',reduc,'_top_gene.pdf', sep=''), p, height=5, width=6, units='in')
ggsave(paste(g_dir, '/',reduc,'_top_gene.png', sep=''), plot=p, height=5, width=6, dpi=600, units='in')
ppca <- T
if('group' %in% names(pbmc@meta.data) && length(unique(pbmc@meta.data$group)) > 1){
	p <- DimPlot(pbmc, reduction =reduc, group.by='group') + scale_color_npg()
	ggsave(paste(g_dir, '/',reduc,'_by_group.pdf', sep=''), p, height=5, width=6,  units='in')
	ggsave(paste(g_dir, '/',reduc,'_by_group.png', sep=''), p, height=5, width=6,  units='in', dpi=600)
	ppca <- F
}
if('sample' %in% names(pbmc@meta.data) && length(unique(pbmc@meta.data$sample)) > 1){
	p <- DimPlot(pbmc, reduction =reduc, group.by='sample') + scale_color_igv()
	ggsave(paste(g_dir, '/',reduc,'_by_sample.pdf', sep=''), p, height=5, width=6,  units='in')
	ggsave(paste(g_dir, '/',reduc,'_by_sample.png', sep=''), p, height=5, width=6,  units='in', dpi=600)
	ppca <- F
}
if(ppca==T){
	pcd <- data.frame(pbmc@reductions$pca@cell.embeddings)[, 1:2]
	pcd$log10_features <- log10(pbmc@meta.data$nFeature_RNA)
	names(pcd)[1:2] <- c('PC1', 'PC2')
	p <- ggplot(pcd, aes(PC1, PC2, col=log10_features))
	p <- p+geom_point(size=0.5, alpha=0.6) + theme_classic()
	p <- p+scale_colour_gradient(low="lightblue",high="red")
	p <- p+guides(color = guide_legend(override.aes = list(size=5)))
	p <- p+theme(axis.text=element_text(size=15, color='black'), axis.title=element_text(size=20), text=element_text(size=15))
	ggsave(paste(g_dir, '/PCA_by_nFeature.pdf', sep=''), p, height=5, width=6,  units='in')
	ggsave(paste(g_dir, '/PCA_by_nFeature.png', sep=''), plot=p, height=5, width=6, dpi=600, units='in')
}

###Select PCs
y <- pbmc[[reduc]]@stdev	
x <- seq(1, length(y))
infl <- data.frame(ede(x, y, index=1))
if(infl$j2 >30 & infl$chi >5){
	npc <- infl$chi
}else{
	npc <- infl$j2
}
p <- ElbowPlot(pbmc, ndims =50, reduction =reduc)
p <- p+geom_vline(xintercept = npc, colour = "red", linetype=2) + xlab(reduc)
p <- p+annotate("text", x = npc, y = max(y)-(max(y)-min(y))*0.3, label = paste(reduc, npc, sep=''), size=5, color='red')
p <- p+theme(axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
ggsave(paste(g_dir, '/',reduc,'_select.pdf', sep=''), p, height=6, width=6, units='in')
ggsave(paste(g_dir, '/',reduc,'_select.png', sep=''), plot=p, height=6, width=6, dpi=600, units='in')
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
p <- DimHeatmap(pbmc, dims = 1:npc, balanced = TRUE, ncol=ncs, fast=F, cells = 500, nfeatures=20, reduction =reduc)
ggsave(paste(g_dir, '/', reduc ,'_heatmap.pdf', sep=''), p, height=h, width=w, device="pdf")
ggsave(paste(g_dir, '/', reduc ,'_heatmap.png', sep=''), plot=p, height=h, width=w, dpi=600, units='in')
##Cluster the cells
pbmc1 <- FindNeighbors(pbmc, dims = 1:npc, reduction =reduc)
res <- ifelse(opt$mcl > 25, 2, 1)
pbmc <- FindClusters(pbmc1, resolution =res)
print(res)
while(length(levels(pbmc@active.ident)) > opt$mcl){
	if(res >= 0.2){
		res <- res - 0.1
	}
	if(res < 0.2 && res >0.02){
		res <- res - 0.02
	}
	if(res < 0.02){
		break
	}
	pbmc <- FindClusters(pbmc1, resolution=res)
}
###TSNE
pbmc <- RunTSNE(pbmc, dims = 1:npc, check_duplicates=F, reduction =reduc)
if('group' %in% names(pbmc@meta.data) && length(unique(pbmc@meta.data$group)) > 1){
	w <- 6 + max(nchar(as.character(pbmc@meta.data$group)))/8
	p <- DimPlot(pbmc, reduction = "tsne", label=F, group.by='group') +  scale_color_npg()
	p <- p + theme(plot.title=element_blank(), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
	ggsave(paste(g_dir, '/TSNE_by_group.pdf', sep=''), p, height=6, width=w, units='in')
	ggsave(paste(g_dir, '/TSNE_by_group.png', sep=''), p, height=6, width=w, dpi=600, units='in')
	p <- DimPlot(pbmc, reduction = "tsne", label=T, split.by='group', group.by='seurat_clusters') + scale_color_igv()
	p <- p + theme(plot.title=element_blank(), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
	w <- 6*length(unique(pbmc@meta.data$group)) + 0.5
	ggsave(paste(g_dir, '/TSNE_split_group.pdf', sep=''), p, height=6, width=w, units='in')
	ggsave(paste(g_dir, '/TSNE_split_group.png', sep=''), p, height=6, width=w, dpi=600, units='in')
}
if('sample' %in% names(pbmc@meta.data) && length(unique(pbmc@meta.data$sample)) > 1){
	w <- 6 + max(nchar(as.character(unique(pbmc@meta.data$sample))))/8
	p <- DimPlot(pbmc, reduction = "tsne", label=F, group.by='sample') + scale_color_igv()
	p <- p + theme(plot.title=element_blank(), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
	ggsave(paste(g_dir, '/TSNE_by_sample.pdf', sep=''), p, height=6, width=w, units='in')
	ggsave(paste(g_dir, '/TSNE_by_sample.png', sep=''), p, height=6, width=w, dpi=600, units='in')
}
p <- DimPlot(pbmc, reduction = "tsne", label=T, group.by='seurat_clusters') + scale_color_igv()
p <- p + theme(plot.title=element_blank(), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
ggsave(paste(g_dir, '/TSNE.pdf', sep=''), p, height=6, width=6.5, units='in')
ggsave(paste(g_dir, '/TSNE.png', sep=''), p, height=6, width=6.5, dpi=600, units='in')

###UMAP
pbmc <- RunUMAP(pbmc, dims = 1:npc, reduction=reduc)
if('group' %in% names(pbmc@meta.data) && length(unique(pbmc@meta.data$group)) > 1){
	w <- 6 + max(nchar(as.character(pbmc@meta.data$group)))/8
	p <- DimPlot(pbmc, reduction = "umap", label=F, group.by='group') + scale_color_npg()
	p <- p + theme(plot.title=element_blank(), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
	ggsave(paste(g_dir, '/UMAP_by_group.pdf', sep=''), p, height=6, width=w, device="pdf", units='in')
	ggsave(paste(g_dir, '/UMAP_by_group.png', sep=''), p, height=6, width=w, device="png", dpi=600, units='in')
	p <- DimPlot(pbmc, reduction = "umap", label=T, split.by='group', group.by='seurat_clusters') + scale_color_igv()
	p <- p + theme(plot.title=element_blank(), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
	w <- 6*length(unique(pbmc@meta.data$group)) + 0.5
	ggsave(paste(g_dir, '/UMAP_split_group.pdf', sep=''), p, height=6, width=w, units='in')
	ggsave(paste(g_dir, '/UMAP_split_group.png', sep=''), p, height=6, width=w, dpi=600, units='in')
}
if('sample' %in% names(pbmc@meta.data) && length(unique(pbmc@meta.data$sample)) > 1){
	w <- 6 + max(nchar(unique(as.character(pbmc@meta.data$sample))))/8
	p <- DimPlot(pbmc, reduction = "umap", label=F, group.by='sample') + scale_color_igv()
	p <- p + theme(plot.title=element_blank(), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
	ggsave(paste(g_dir, '/UMAP_by_sample.pdf', sep=''), p, height=6, width=w, device="pdf", units='in')
	ggsave(paste(g_dir, '/UMAP_by_sample.png', sep=''), p, height=6, width=w, device="png", dpi=600, units='in')
}
p <- DimPlot(pbmc, reduction = "umap", label=T, group.by='seurat_clusters') + scale_color_igv()
p <- p + theme(plot.title=element_blank(), axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
ggsave(paste(g_dir, '/UMAP.pdf', sep=''), p, height=6, width=6.5, device="pdf", units='in')
ggsave(paste(g_dir, '/UMAP.png', sep=''), p, height=6, width=6.5, device="png", dpi=600, units='in')
#saveRDS(pbmc, file = paste(opt$out, "pbmc_tutorial.rds", sep=''))
###ClusterTree
pbmc <- BuildClusterTree(pbmc, reduction=reduc)
pdf(paste(g_dir, '/cluster_tree.pdf', sep=''))
PlotClusterTree(pbmc)
dev.off()

cid <- paste(g_dir, '/cell_identify_data', sep='')
if(!dir.exists(cid)){
	dir.create(cid, recursive=T)
}

exp <- try(t(as.matrix(pbmc@assays$RNA@data)))
if(!"try-error" %in% class(exp)){
	expo <- data.frame(Cell_id=rownames(exp), exp, stringsAsFactors=F, check.names=F)
	write.csv(expo, paste(cid, '/exp.csv', sep=''), quote=F, row.names=F)
}else{
	saveRDS(pbmc@assays$RNA@data, paste(cid, '/exp.rds', sep=''))
}
cooro <- pbmc@reductions$umap@cell.embeddings
colnames(cooro) <- c('Coordinate_1', 'Coordinate_2')
write.csv(cooro, paste(cid, '/Coor_UMAP.csv', sep=''), quote=F, row.names=F)
cooro <- pbmc@reductions$tsne@cell.embeddings
colnames(cooro) <- c('Coordinate_1', 'Coordinate_2')
write.csv(cooro, paste(cid, '/Coor_tSNE.csv', sep=''), quote=F, row.names=F)
clusterl <- c('Cluster_Label', as.character(pbmc@meta.data$seurat_clusters))
clusterl[nchar(as.character(clusterl))==1] <- paste('0', clusterl[nchar(as.character(clusterl))==1], sep='')
writeLines(clusterl, paste(cid, '/Label.csv', sep=''))
##cluster ratio stat
if('group' %in% names(pbmc@meta.data) && length(unique(pbmc@meta.data$group)) > 1){
	ratio.stat <- pbmc@meta.data %>% dplyr::group_by(seurat_clusters, group) %>% summarise(count=length(group)) %>% group_by(group) %>% dplyr::mutate(total=sum(count), ratio=count/total) %>% data.frame(stringsAsFactors=F)
	write.table(ratio.stat, paste(g_dir, '/cell_ratio_by_clusters.txt', sep=''), sep='\t', quote=F, row.names=F)
	p <- ggplot(ratio.stat, aes(x=seurat_clusters, y=ratio, fill=group)) + geom_col(position = "dodge2")
	p <- p + labs(x='Cluster', y='Cell ratio') + scale_fill_npg() + theme_bw()
	p <- p + theme(axis.text=element_text(size=15, color='black'), axis.title=element_text(size=20), 
				   legend.title=element_blank(), legend.text=element_text(size=15), legend.position='top')
	if(length(unique(ratio.stat$seurat_clusters)) > 15){
		w <- 0.4*length(unique(ratio.stat$seurat_clusters))
		}else{
			w <- 6
	}
	ggsave(paste(g_dir, '/cell_ratio_by_clusters.png', sep=''), plot=p, width=w, height=5, units='in', dpi=600)
	ggsave(paste(g_dir, '/cell_ratio_by_clusters.pdf', sep=''), plot=p, width=w, height=5, units='in')
	##by sample
	if(any(table(unique(pbmc@meta.data[, c('group', 'sample')])[, 'group']) > 1)){
		cstat <- pbmc@meta.data %>% group_by(seurat_clusters, group, sample) %>% summarise(count=length(sample)) %>% group_by(sample) %>% mutate(total=sum(count)) %>% group_by(seurat_clusters) %>% mutate(Ratio=count/total) %>% data.frame(stringsAsFactors=F)
		if(any(table(unique(pbmc@meta.data[, c('group', 'sample')])[, 'group']) > 2)){
			pd <- cstat %>% group_by(seurat_clusters) %>% filter(length(unique(group)) > 1) %>% group_by(seurat_clusters) %>% summarise(p=kruskal.test(formula=Ratio~group)$p.value) %>% data.frame(stringsAsFactors=F)
			pdo <- cstat %>% group_by(seurat_clusters) %>% filter(length(unique(group)) > 1) %>% group_by(seurat_clusters) %>% mutate(p_value=kruskal.test(formula=Ratio~group)$p.value) %>% data.frame(stringsAsFactors=F)
		}else{
			pd <- cstat %>% group_by(seurat_clusters) %>% filter(length(unique(group)) > 1) %>% group_by(seurat_clusters) %>% summarise(p=wilcox.test(formula=Ratio~group)$p.value) %>% data.frame(stringsAsFactors=F)
			pdo <- cstat %>% group_by(seurat_clusters) %>% filter(length(unique(group)) > 1) %>% group_by(seurat_clusters) %>% mutate(p_value=wilcox.test(formula=Ratio~group)$p.value) %>% data.frame(stringsAsFactors=F)
		}
		names(pd)[1] <- names(pdo)[1] <- names(cstat)[1] <- 'Cluster'
		write.table(pdo, paste(g_dir, '/cell_ratio_by_sample.txt', sep=''), sep='\t', quote=F, row.names=F)
		pd$labe <- ''
		pd$labe[pd$p < 0.001] <- '***'
		pd$labe[pd$p < 0.01 & pd$p >=0.001] <- '**'
		pd$labe[pd$p < 0.05 & pd$p >=0.01] <- '*'
		labe.caption <- expression(paste('*: ', 0.01 <= p, '< 0.05;  **: ', 0.001 <= p, '< 0.01;  ***: p < 0.001'))
		p <- ggplot(cstat, aes(x=Cluster, y=Ratio, fill=group)) + geom_boxplot(outlier.shape=NA)
		gp.stat <- ggplot_build(p)$data[[1]]
		pd.stat <- gp.stat[, c('x', 'ymin', 'ymax')]
		p <- p + geom_text(data=pd, aes(x=Cluster, y=max(pd.stat$ymax)+(max(pd.stat$ymax)-min(pd.stat$ymin))*0.01, label=labe), inherit.aes=F, size=7)
		p <- p + coord_cartesian(ylim=c(min(pd.stat$ymin), max(pd.stat$ymax))) + scale_fill_npg() 
		p <- p + labs(caption=labe.caption) + theme_bw()
		p <- p + theme(axis.title=element_text(size=20), axis.text=element_text(size=15, color='black'),
					   legend.title=element_blank(), legend.text=element_text(size=15), legend.position='top')
		if(length(unique(cstat$Cluster)) > 15){
			w <- 0.4*length(unique(cstat$Cluster))
		}else{
			w <- 6
		}
		ggsave(paste(g_dir, '/cell_ratio_by_sample.png', sep=''), plot=p, width=w, height=5, dpi=600, units='in')
		ggsave(paste(g_dir, '/cell_ratio_by_sample.pdf', sep=''), plot=p, width=w, height=5, units='in')
	}
}
##cluster biomarkers
#pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,  return.thresh=0.05, test.use='MAST')
if(DefaultAssay(pbmc)!="RNA"){
	DefaultAssay(pbmc) <- "RNA"
	all.genes <- rownames(pbmc)
	pbmc <- ScaleData(pbmc, features = all.genes)
}
saveRDS(pbmc, file = paste(g_dir, '/cell_cluster.rds', sep=''))

mc <- getOption("mc.cores", opt$thread)
maker.list <- mclapply(levels(pbmc@active.ident), fdmk, pbmc=pbmc, mc.cores = mc)
maker_gene <- do.call(rbind, maker.list)
write.table(maker_gene, file=paste(g_dir, '/marker_gene.txt', sep=''), quote=F, row.names=F, sep='\t')
top_maker_gene <- maker_gene %>% dplyr::group_by(cluster) %>% top_n(n = 1, wt = abs(avg_log2FC))
ntg <- nrow(top_maker_gene)
ncs <- ceiling(sqrt(ntg))
nrs <- ceiling(ntg/ncs)
p <- FeaturePlot(pbmc, features =top_maker_gene$gene, ncol=ncs)
if(ntg>2){
	w <- ncs*4
	h <- nrs*3
}else{
	h <- 3
	w <- ntg*3
}
ggsave(paste(g_dir, '/top_marker_gene_exp.pdf', sep=''), p, width=w, height=h, device="pdf")
ggsave(paste(g_dir, '/top_marker_gene_exp.png', sep=''), p, width=w, height=h, device="png", dpi=600, units='in')

top10 <- maker_gene %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
ncl <- length(levels(pbmc@active.ident))
if(ncl <=7){
	w <- 7
	h <- 7
}else{
	w <- ncl
	h <- ncl
}
xlsz <- 500/length(top10$gene)
p <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
p <- p + theme(axis.text.y=element_text(size=xlsz, color='black'))
ggsave(paste(g_dir, '/top10_marker_gene_heatmap.pdf', sep=''), p, width=w, height=h, units='in')
ggsave(paste(g_dir, '/top10_marker_gene_heatmap.png', sep=''), p, width=w, height=h, units='in', dpi=300)
##Assigning cell type identity to clusters
if(!is.null(opt$tissue)){
	if(opt$ct){
		sct <- '-c'
	}else{
		sct <- ''
	}
	comd <- paste('Rscript /home/dev/DEV-cuiy/cellTypeEnrich.r -f=', g_dir, '/marker_gene.txt -t=', opt$tissue, ' -o=marker_anno ', sct, sep='')
	print(comd)
	system(comd)
	makeh <- read.table(paste('marker_anno.txt', sep=''), header=T, sep='\t', stringsAsFactors=F)[, c('Cell_Name', 'Cluster', 'p_value')]
	clt <- makeh %>% dplyr::group_by(Cluster) %>% top_n(n = 1, wt = -p_value)
	cell_type <- as.character(clt$Cell_Name)
	names(cell_type) <- clt$Cluster
	if(any(! levels(pbmc) %in% names(cell_type))){
		nct <- as.character(levels(pbmc)[! levels(sub) %in% names(cell_type)])
		cnct <- paste('Cluster_', nct, sep='')
		names(cnct) <- nct
		cell_type <- c(cell_type, cnct)
	}
	pbmc2 <- RenameIdents(pbmc, cell_type)
	p <- DimPlot(pbmc2, reduction = "tsne", label = TRUE, pt.size = 0.5)
	p <- p + scale_color_igv()
	ggsave(paste(g_dir, '/tSNE_cell_type.pdf', sep=''), p, height=6, width=8, units='in')
	ggsave(paste(g_dir, '/tSNE_cell_type.png', sep=''), plot=p, height=6, width=8, dpi=600, units='in')
	p <- DimPlot(pbmc2, reduction = "umap", label = TRUE, pt.size = 0.5)
	#ggsave(paste(g_dir, '/UMAP_cell_type.pdf', sep=''), p, height=6, width=8, units='in')
	#ggsave(paste(g_dir, '/UMAP_cell_type.png', sep=''), plot=p, height=8, width=9, dpi=600, units='in')
}
pbmc@misc<-list(marker.gene=maker_gene)
if(!file_test('-f', paste(g_dir, '/cell_cluster.rds', sep=''))){
	saveRDS(pbmc, file = paste(g_dir, '/cell_cluster.rds', sep=''))
}


