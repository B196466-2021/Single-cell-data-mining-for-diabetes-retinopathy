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
updated date: 2021-05-08\n\n")
	q(status=1)
}

usage<-function(spec){
	cat("This script is used to infer and analyze cell-cell communication\n",
getopt(spec,usage=TRUE),
"Options:
	-f, --file	rds file of Seurat object after cell type.
	-s, --species	species: human or mouse (default human).
	-t, --thread	thread number (default 4).
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
	'species','s',1,'character',
	'pval','p',1,'double',
	'thread','t',1,'integer',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$thread)) {opt$thread=4}
if (is.null(opt$species)) {opt$species='human'}
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
library(CellChat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(stringr)
options(stringsAsFactors = FALSE)

sda <- readRDS(opt$file)
cellchat <- createCellChat(object=sda, group.by="cell_type")
if(opt$species=='human'){
	CellChatDB <- CellChatDB.human
	ppiDB <- PPI.human
}else if(opt$species=='mouse'){
	CellChatDB <- CellChatDB.mouse
	ppiDB <- PPI.mouse
}
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers=opt$thread)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, ppiDB)
cellchat <- computeCommunProb(cellchat, trim=0.1)
cellchat <- filterCommunication(cellchat, min.cells = round(min(table(sda@meta.data$cell_type)*0.25)))
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
colnames(cellchat@net$count) <- str_wrap(colnames(cellchat@net$count), width=20)
rownames(cellchat@net$count) <- str_wrap(rownames(cellchat@net$count), width=20)
write.table(data.frame(cells=rownames(cellchat@net$count), cellchat@net$count, check.names=F), paste(opt$out, 'network_count.txt', sep=''), sep='\t', quote=F, row.names=F)
png(paste(opt$out, 'communication_network.png', sep=''), width=7, height=7, units='in', res=600)
p <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name="Number of interactions")
dev.off()
pdf(paste(opt$out, 'communication_network.pdf', sep=''), width=7, height=7)
p <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name="Number of interactions")
dev.off()
#ggsave(paste(opt$out, 'communication_network.png', sep=''), plot=ggdraw(p), width=7, height=7, units='in', dpi=600)
#ggsave(paste(opt$out, 'communication_network.pdf', sep=''), plot=ggdraw(p), width=7, height=7, units='in')

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
p <- netAnalysis_signalingRole_scatter(cellchat)
p <- p + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=12, color='black'), text=element_text(size=12))
ggsave(paste(opt$out, 'communication_statistics.png', sep=''), plot=p, width=6, height=5, units='in', dpi=600)
ggsave(paste(opt$out, 'communication_statistics.pdf', sep=''), plot=p, width=6, height=5, units='in')

write.table(cellchat@LR$LRsig, paste(opt$out, 'communication_LR_annotation.txt', sep=''), sep='\t', quote=F, row.names=F)
df.net <- subsetCommunication(cellchat)
write.table(df.net, paste(opt$out, 'communication_sig_result.txt', sep=''), sep='\t', quote=F, row.names=F)
df.net <- subsetCommunication(cellchat, thresh=NULL)
write.table(df.net, paste(opt$out, 'communication_all_result.txt', sep=''), sep='\t', quote=F, row.names=F)
save.image()
