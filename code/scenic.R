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
updated date: 2021-01-13\n\n")
	q(status=1)
}

usage<-function(spec){
	cat("This script is used to run SCENIC with scRNA-seq data.\n",
getopt(spec,usage=TRUE),
"Options:
	-f, --file	rds file of a Seurat object or matrix file of expression counts (gene in row with rownames and cell in column with colnames).
	-g, --genome	hg19/hg38/mm9/mm10 (default hg38).
	-b, --group	which information to group the cells in Seurat object metadata (sample, group, seurat_clusters, default NULL) or group information file (cell<tab>sample<tab>group).
	-t, --thread	thread number (default 5).
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
	'genome','g',1,'character',
	'group','b',1,'character',
	'thread','t',1,'integer',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || !is.null(opt$help)) {usage(spec)}
if (is.null(opt$genome) || opt$genome=='hg38'){
	org <- 'hgnc'
	dbs <- c('hg38-500bp_up_and_100bp_down_tss.mc9nr.feather', 'hg38-tss-centered-10kb.mc9nr.feather')
}else if(opt$genome=='hg19'){
	org <- 'hgnc'
	dbs <- c('hg19-500bp-upstream-7species.mc9nr.feather', 'hg19-tss-centered-10kb-7species.mc9nr.feather')
}else if(opt$genome=='mm9'){
	org <- 'mgi'
	dbs <- c('mm9-500bp-upstream-7species.mc9nr.feather', 'mm9-tss-centered-10kb-7species.mc9nr.feather')
}else if(opt$genome=='mm10'){
	org <- 'mgi'
	dbs <- c('mm10-500bp_up_and_100bp_down_tss.mc9nr.feather', 'mm10-tss-centered-10kb.mc9nr.feather')
}
names(dbs) <- c('500bp', '10kb')

if(is.null(opt$thread)) {opt$thread=5}
if(is.null(opt$out)) {opt$out='./'}

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
library(SCopeLoomR)
library(SCENIC)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(AUCell)
library(ggsci)
if(grepl('\\.rds$', opt$file, ignore.case=T, perl=T)){
	pbmc <- readRDS(opt$file)
	exprMat <- as.matrix(pbmc@assays$RNA@counts)
	cellInfo <- pbmc@meta.data
	rm(pbmc)
}else{
	exprMat <- as.matrix(read.table(opt$file, header=T, sep='\t', row.names=1, stringsAsFactors=F, check.names=F))
}

if(!is.null(opt$group) && file_test("-f", opt$group)){
	group <- read.table(opt$group, header=T, sep='\t', stringsAsFactors=F, check.names=F)
	group <- group[group[,1] %in% colnames(exprMat),]
	cellInfo <- matrix(group[, 2], ncol=1)
	rownames(cellInfo) <- group[,1]
	colnames(cellInfo)[1] <- colnames(opt$group)[2]
	opt$group <- colnames(cellInfo)[1]
}

scenicOptions <- initializeScenic(org=org, dbDir='/opt/databases/transcription_factor/cistarget', dbs=dbs, datasetTitle='my_SCENIC', nCores=opt$thread)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["500bp"]
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.rds1")
scenicOptions@settings$nCores <- 10
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC")
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
export2loom(scenicOptions, exprMat)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
auc <- getAUC(regulonAUC)
auco <- data.frame(TF=rownames(auc), auc, check.names=F, stringsAsFactors=F)
write.table(auco, paste(opt$out, '/regulon_activity.txt', sep=''), sep='\t', quote=F, row.names=F)
saveRDS(regulonAUC, paste(opt$out, '/regulonAUC.rds', sep=''))
if(!is.null(opt$group)){
	if(length(unique(cellInfo[, opt$group]))>2){
		regulonActivity_byCellType <- sapply(split(rownames(cellInfo), as.character(cellInfo[, opt$group])), function(cells){da=data.frame(getAUC(regulonAUC)[,cells]); rownames(da)=rownames(getAUC(regulonAUC)); rowMeans(da)})
		regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
		write.table(data.frame(names=rownames(regulonActivity_byCellType_Scaled), regulonActivity_byCellType_Scaled), paste(opt$out, '/regulonActivity_byCellType_Scaled.txt', sep=''), sep='\t', quote=F, row.names=F)
		pd <- regulonActivity_byCellType_Scaled[!grepl('_extended', rownames(regulonActivity_byCellType_Scaled)), ]
		if(ncol(pd) > 20){
			show_col_names <- F
		}else{
			show_col_names <- T
		}
		if(all(c('sample', 'group') %in% colnames(cellInfo))){
			ga <- unique(cellInfo[, c('sample', 'group')])
			gao <- ga[order(ga$group),]
			ga <- matrix(gao$group, ncol=1)
			rownames(ga) <- gao$sample
			colnames(ga) <- 'group'
			acl <- pal_igv()(length(unique(gao$group)))
			names(acl) <- unique(gao$group)
			rz <- nrow(pd)*0.15
			p <- Heatmap(pd, heatmap_legend_param=list(title = ""), column_names_rot=45, row_names_gp=gpar(fontsize=rz), cluster_columns=F, column_split=gao$group,
						 top_annotation=HeatmapAnnotation(df=ga, show_legend=F, show_annotation_name=F, col=list(group=acl)), column_gap = unit(0.5, "mm"), 
						 column_names_gp=gpar(fontsize=10), show_column_names=show_col_names)
		}else{
			p <- Heatmap(pd, heatmap_legend_param=list(title = ""), column_names_rot=45, row_names_gp=gpar(fontsize=rz), 
						 column_names_gp=gpar(fontsize=10), show_column_names=show_col_names)
		}
		png(paste(opt$out, '/regulonActivity_heatmap.png', sep=''), width=7, height=7, units = "in", res=600)
		p
		dev.off()
		pdf(paste(opt$out, '/regulonActivity_heatmap.pdf', sep=''), width=7, height=7)
		p
		dev.off()
		#ggsave(paste(opt$out, '/regulonActivity_heatmap.png', sep=''), plot=ggplotify::as.ggplot(p), width=7, height=7, units = "in", dpi=600)
		#ggsave(paste(opt$out, '/regulonActivity_heatmap.pdf', sep=''), plot=ggplotify::as.ggplot(p), width=7, height=7, units = "in")
	}else{
		phod <- data.frame(Sample=rownames(cellInfo), Group=cellInfo[, opt$group], stringsAsFactors=F)
		write.table(phod, paste(opt$out, '/group_inf.txt', sep=''), sep='\t', quote=F, row.names=F)
	}
	#rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "Cell_Type"])
	#write.table(data.frame(names=rownames(rss), rss), paste(opt$out, '/rss.txt', sep=''), sep='\t', quote=F, row.names=F)
}

