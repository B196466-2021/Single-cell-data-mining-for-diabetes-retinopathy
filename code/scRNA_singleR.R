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
	-f, --file	rds file result from scRNA_cluster.R.
	-r, --ref	a SummarizedExperiment object: HumanPrimaryCellAtlasData, MouseRNAseqData, NovershternHematopoieticData, ImmGenData, DatabaseImmuneCellExpressionData, MonacoImmuneData, BlueprintEncodeData or a rds file of a SummarizedExperiment object (default HumanPrimaryCellAtlasData).
	-o, --out	prefix of out files (default ./).
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
	'ref','r',1,'character',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$ref)) {opt$ref='HumanPrimaryCellAtlasData'}
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
library(celldex)
library(SingleR)
library(SummarizedExperiment)
library(ggplot2)
library(ggsci)
library(openxlsx)
library(ggplotify)
seud <- readRDS(opt$file)
colData <- DataFrame(seud@meta.data[colnames(seud@assays$RNA@counts),])
sedo <- SummarizedExperiment(assays=list(counts=seud@assays$RNA@counts), colData=colData)
if(file_test('-f', opt$ref)){
	ref <- readRDS(opt$ref)
}else if(opt$ref %in% c('HumanPrimaryCellAtlasData', 'MouseRNAseqData', 'NovershternHematopoieticData', 'ImmGenData', 'DatabaseImmuneCellExpressionData', 'MonacoImmuneData', 'BlueprintEncodeData')){
	ref <- get(opt$ref)()
}
pred <- SingleR(tes=sedo, ref=ref, assay.type.test=1, labels=ref$label.main, clusters=sedo@colData@listData$seurat_clusters)
od <- data.frame(Cluster=rownames(pred), Cell_Name=pred$labels, Genes='NA')
write.xlsx(od, paste(opt$out, 'cluster2celltype.xlsx', sep=''), overwrite=T)

p <- plotScoreHeatmap(pred, show_colnames=T, silent=T)
ggsave(paste(opt$out, "heatmap.png", sep=""), plot=as.ggplot(p), width=10, height=7, units='in', dpi=600, bg='white')
ggsave(paste(opt$out, "heatmap.pdf", sep=""), plot=as.ggplot(p), width=10, height=7, units='in', bg='white')

cell_type <- pred$labels
names(cell_type) <- rownames(pred)
sub2 <- RenameIdents(seud, cell_type)
p <- DimPlot(sub2, reduction='tsne', label = TRUE, label.size=3)
p2 <- DimPlot(sub2, reduction='umap', label = TRUE, label.size=3)
p <- p + labs(x='tSNE1', y='tSNE2') + scale_color_igv()
p <- p + theme(axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
p2 <- p2 + labs(x='UMAP1', y='UMAP2') + scale_color_igv()
p2 <- p2 + theme(axis.text=element_text(size=15, color='black'), axis.title=element_text(size=15))
w <- 6 + max(nchar(cell_type))/10
ggsave(paste(opt$out, 'tSNE_cell_type.pdf', sep=''), plot=p, height=6, width=w, units='in')
ggsave(paste(opt$out, 'tSNE_cell_type.png', sep=''), plot=p, height=6, width=w, dpi=600, units='in')

ggsave(paste(opt$out, 'UMAP_cell_type.pdf', sep=''), plot=p2, height=6, width=w, units='in')
ggsave(paste(opt$out, 'UMAP_cell_type.png', sep=''), plot=p2, height=6, width=w, dpi=600, units='in')

