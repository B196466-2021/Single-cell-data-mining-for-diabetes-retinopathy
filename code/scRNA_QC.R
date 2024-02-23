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
	-d, --dir	directory containing all single cell RNA-seq matrices (default ./) or matrix file with all samples and groups.
	-p, --pheno	group information file.
	-r, --rib	minimum ribosome gene expression persent used to filter cell (default 10)
	-m, --mit	maximum mitochondria gene expression persent used to filter cell (default 10)
	-t, --thread	number of threads (default 1).
	-b, --double	remove doublets or not (default TRUE).
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
	'dir','d',1,'character',
	'pheno','p',1,'character',
	'rib','r',1,'double',
	'mit','m',1,'double',
	'thread','t',1,'integer',
	'double','b',2,'logical',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$pheno) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$rib)) {opt$rib=10}
if (is.null(opt$mit)) {opt$mit=10}
if (is.null(opt$thread)) {opt$thread=1}
if (is.null(opt$out)) {opt$out='./'}
if (is.null(opt$dir)) {opt$dir='./'}
if (is.null(opt$double)) {opt$double=TRUE}
##########################
#define functions
my_fun<-function(s=NULL, ga, dir, ribp=10, mitp=10, odir='./'){
	##loading data
	if(is.null(s) && file_test('-f', dir)){
		print('All input data are in one matrix.')
		readsCount <- read.table(dir, header = T, row.names = 1, check.names=F, na.strings=c("NA", "na", "NULL", "null"))
		sce <- SingleCellExperiment(assays = list(counts = as.matrix(readsCount)))
		s='All'
	}else if(!is.null(s) && file_test('-d', dir)){
		psd <- paste(dir, '/', s, sep='')
		if(file_test('-d', psd)){
			print('Input data are 10X genomics format.')
			red10x <- Read10X(data.dir = psd)
			sce <- SingleCellExperiment(assays = list(counts = red10x))
		}else if(file_test('-f', paste(psd, '.txt', sep=''))){
			print('Input data are matrix format.')
			readsCount <- read.table(paste(psd, '.txt', sep=''), header = T, check.names=F, na.strings=c("NA", "na", "NULL", "null"), sep='\t', stringsAsFactors=F)
			if(class(readsCount[,1])=='character'){
				if(any(duplicated(readsCount[,1]))){
					dam <- aggregate(readsCount[,-1], by=list(gene=readsCount[,1]), mean)
					readsCount <- dam[, -1]
					rownames(readsCount) <- dam$gene
				}else{
					rownames(readsCount) <- readsCount[,1]
					readsCount <- readsCount[,-1]
				}
			}
			readsCount <- as.matrix(readsCount)
			readsCount[is.na(readsCount)] <- 0
			sce <- SingleCellExperiment(assays = list(counts = as.matrix(readsCount)))
		}
	}
	cell.raw.count <- sce@colData@nrows
	bcrank <- barcodeRanks(counts(sce))
	uniq <- !duplicated(bcrank$rank)
	
	pdf(paste(odir, '/', s, '_UMI_count_depth.pdf', sep=''))
	plot(bcrank@listData$rank[uniq], bcrank@listData$total[uniq], log="xy", xlab="Rank", ylab="Total UMI count", cex=0.5, cex.lab=1.2, col='steelblue3')
	abline(h=bcrank@metadata$inflection, col="darkgreen", lty=2)
	abline(h=bcrank@metadata$knee, col="dodgerblue", lty=2)
	legend("left", legend=c("Inflection", "Knee"), bty="n",col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
	dev.off()
	
	if(any(bcrank@listData$total<200)){
		print('Filtering cells with UMI and FDR...')
		set.seed(100)
		e.out <- try(emptyDrops(counts(sce), lower=200), silent=T)
		if(!'try-error' %in% class(e.out)){
			is.cell <- (e.out$FDR <= 0.01)
			pdf(paste(odir, '/', s, '_cells_filter.pdf', sep=''))
			plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"), xlab="Total UMI count", ylab="-Log Probability", cex=0.2)
			abline(v = bcrank@metadata$inflection, col="darkgreen")
			abline(v = bcrank@metadata$knee, col="dodgerblue")
			legend("bottomright", legend=c("Inflection", "Knee"), bty="n", col=c("darkgreen", "dodgerblue"), lty=1, cex=1.2)
			dev.off()
			w2kp <- which(is.cell & e.out$Total >= bcrank@metadata$inflection)
			sce <- sce[, w2kp]
		}
	}
	
	##Mitochondrial or Ribosomal genes
	ribo.file <- '/opt/Genome/GRCh38/Annotation/ribosomal_gene.txt'
	if(file.exists(ribo.file)){
		ribo <- read.table(ribo.file, header=T, sep='\t', stringsAsFactors=F)
		riboh <- rownames(counts(sce))[(rownames(counts(sce)) %in% ribo$Approved.symbol)]
	}else{
		riboh <- c()
	}
	
	is.mito <- rownames(counts(sce))[grepl('^MT-', rownames(counts(sce)), ignore.case=T, perl=T)]
	is.ribo <- rownames(counts(sce))[grepl('^MRPL|^MRPS|^RIMKL|^RNA5|^RNA18S|^RNA28S|^RNA45S|^RPF|^RPL|^RPS|^RRP|^RSL', rownames(counts(sce)), ignore.case=T, perl=T)]
	is.ribo <- unique(c(is.ribo, riboh))
	
	#sce = calculateQCMetrics(sce, feature_controls=list(Mt=is.mito, Ri=is.ribo))
	sce.cell.stat = perCellQCMetrics(sce, subsets=list(Mt=is.mito, Ri=is.ribo))
	sce <- sce[,which(sce.cell.stat$subsets_Mt_percent <= mitp | sce.cell.stat$subsets_Ri_percent >= ribp)]

	pdf(paste(odir, '/', s, '_counts_in_feature.pdf', sep=''))
	par(mfrow=c(2,2), mar=c(5, 4, 1, 1), bty="n")
	hist(log10(sce.cell.stat$total), xlab="log10(Counts)", main="", breaks=20, col="grey80", ylab="Number of cells")
	hist(log10(sce.cell.stat$detected), xlab="log10(# of expressed genes)", main="", breaks=20, col="grey80", ylab="Number of cells")
	hist(sce.cell.stat$subsets_Ri_percent, xlab="Ribosome prop. (%)", ylab="Number of cells", breaks=40, main="", col="grey80")
	hist(sce.cell.stat$subsets_Mt_percent, xlab="Mitochondrial prop. (%)", ylab="Number of cells", breaks=80, main="", col="grey80")
	dev.off()

	pdf(paste(odir, '/', s, '_counts_vs_feature.pdf', sep=''))
	par(mfrow=c(2, 2), mar=c(5, 4, 1, 1), bty="n")
	smoothScatter(log10(sce.cell.stat$total), log10(sce.cell.stat$detected),  xlab="log10(Library sizes)", ylab="log10(# of expressed genes)", nrpoints=500, cex=0.5)
	smoothScatter(log10(sce.cell.stat$total), sce.cell.stat$subsets_Ri_percent, xlab="log10(Library sizes)", ylab="Ribosome prop. (%)", nrpoints=500, cex=0.5)
	abline(h=ribp, lty=1)
	smoothScatter(log10(sce.cell.stat$total), sce.cell.stat$subsets_Mt_percent, xlab="log10(Library sizes)", ylab="Mitochondrial prop. (%)", nrpoints=500, cex=0.5)
	abline(h=mitp, lty=1)
	smoothScatter(sce.cell.stat$subsets_Ri_percent, sce.cell.stat$subsets_Mt_percent, xlab="Ribosome prop. (%)", ylab="Mitochondrial prop. (%)", nrpoints=500, cex=0.5)
	abline(h=mitp,  lty=1)
	abline(v=ribp, lty=1)
	dev.off()
	
	##Summarize gene-level information
	#sce.feature.stat = perFeatureQCMetrics(sce)
	#pdf(paste(opt$out, '/QC/', s, '_UMI_vs_cells.pdf', sep=''), width=20)
	#par(mfrow=c(1,3), mar=c(5,4,1,1))
	#hist(log10(sce.feature.stat$mean+1e-6), col="grey80",  main="", breaks=40, xlab="log10(ave # of UMI + 1e-6)")
	#hist(log10(sce.feature.stat$detected+1), col="grey80", main="", breaks=40, xlab="log10(# of expressed cells + 1)")
	#smoothScatter(log10(sce.feature.stat$mean+1e-6), log10(sce.feature.stat$detected + 1), xlab="log10(ave # of UMI + 1e-6)", ylab="log10(# of expressed cells + 1)")
	#dev.off()
	#od1 = order(sce.feature.stat$mean, decreasing = TRUE)
	#sgd <- data.frame(gene=rownames(sce.feature.stat)[od1[20:1]], exp=sce.feature.stat$mean[od1[20:1]], ncell=sce.feature.stat$detected[od1[20:1]])
	#sgd$gene <- factor(sgd$gene, levels=rownames(rowData(sce))[od1[20:1]])
	#p <- ggplot(sgd, aes(gene, exp, fill=exp))+geom_col(show.legend=F)+coord_flip()
	#p <- p+geom_text(aes(y=exp+(max(exp*0.03)), label=ncell))+theme_bw()
	#p <- p+labs(x='', y='ave # of UMI')
	#p <- p+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
	#ggsave(paste(opt$out, '/QC/', s, '_highly_expressed_genes.pdf', sep=''), p, width=9, device="pdf")
	sce <- sce[!grepl("MALAT1", rownames(sce), ignore.case=T), ]	
	sce <- sce[!grepl("^MT-", rownames(sce), ignore.case=T), ]
	##Normalization
	if(s != 'All'){
		pbmc <- CreateSeuratObject(counts = counts(sce), project = s, min.cells = 3, min.features = 200)
		pbmc@meta.data$sample <- s
		pbmc@meta.data$group <- ga$Group[ga$Sample==s]
		cell.clean.count <- nrow(pbmc@meta.data)
		retl <- list(data=pbmc, cst=data.frame(Sample=s, cell.raw.count=cell.raw.count, cell.clean.count=cell.clean.count))
	}else{
		pbmc <- CreateSeuratObject(counts = counts(sce), project = 'all', min.cells = 3, min.features = 200)
		pbmc@meta.data$sample <- gsub('(.*)_[A-Z]*?$', '\\1', rownames(pbmc@meta.data))
		while(length(unique(pbmc@meta.data$sample)) > 100 && all(grepl('_[A-Z]', pbmc@meta.data$sample))){
			pbmc@meta.data$sample <- gsub('(.*)_[A-Z]*?$', '\\1', pbmc@meta.data$sample)
		}
		cells <- rownames(pbmc@meta.data)[pbmc@meta.data$sample %in% ga$Sample]
		pbmc <- subset(pbmc, cells=cells)
		pbmc@meta.data$group <- ga[pbmc@meta.data$sample, 'Group']
		pbmc@meta.data$orig.ident <- pbmc@meta.data$sample
		cell.clean.count <- nrow(pbmc@meta.data)
		retl <- list(data=pbmc, cst=data.frame(Sample='All', cell.raw.count=cell.raw.count, cell.clean.count=cell.clean.count))
	}

	#Predict doublets
	if(opt$double){
		sweep.res <- paramSweep_v3(pbmc, PCs=1:15, num.cores=round(20/opt$thread))
		sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
		bcmvn <- find.pK(sweep.stats)
		graphics.off()
		pk_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
		DoubletRate = ncol(pbmc)*8*1e-6
		nExp_poi <- round(DoubletRate*ncol(pbmc))
		pbmc <- doubletFinder_v3(pbmc, pN = 0.25, pK = pk_bcmvn, nExp = nExp_poi, PCs = 1:15)
		DF.name = colnames(pbmc@meta.data)[grepl("DF.classification", colnames(pbmc@meta.data))]
		p <- VlnPlot(pbmc, features="nFeature_RNA", group.by=DF.name, pt.size=0) + theme(axis.text=element_text(size=15, color='black'), axis.title=element_blank(), plot.title=element_blank(), legend.position='none')
		ggsave(paste(odir, '/', s, '_doublet_violin.pdf', sep=''), p, height=6, width=6, units='in')
		ggsave(paste(odir, '/', s, '_doublet_violin.png', sep=''), p, height=6, width=6, units='in', dpi=600)
		pbmc <- pbmc[, pbmc@meta.data[, DF.name]=="Singlet"]
	}
	retl$data <- pbmc
	retl$cst$cell.clean.count <- nrow(pbmc@meta.data)
	return(retl)
}

##main program
library(DropletUtils)
library(biomaRt)
library(dplyr)
library(scater)
library(ggplot2)
library(Matrix)
library(Seurat)
library(scran)
library(gplots)
library(parallel)
library(DoubletFinder)
source('/home/dev/DEV-wangdh/R/paramSweep_v3.R')
source('/home/dev/DEV-wangdh/R/doubletFinder_v3.R')

group <- read.table(opt$pheno, header=T, sep='\t', stringsAsFactors=F)
group.order <- unique(group$Group)
group$order <- 1:nrow(group)
osd <- group %>% group_by(Group) %>% summarise(ma=max(order), mi=min(order)) %>% data.frame()
if(nrow(group) >2){
	tf <- c()
	for(i in 2:length(group.order)){
		tf <- c(tf, osd[osd$Group==group.order[i-1], 'ma'] < osd[osd$Group==group.order[i], 'mi'])
	}
	if(!all(tf)){
		group.order <- sort(unique(group$Group))
	}
}else{
	group.order <- sort(group$Group)
}

rownames(group) <- group$Sample
qc_dir <- paste(opt$out, '/QC', sep='')
if(! dir.exists(qc_dir)){
	dir.create(qc_dir, recursive=T)
}
if(file_test('-d', opt$dir)){
	mc <- getOption("mc.cores", opt$thread)
	ret <- mclapply(group$Sample, my_fun, dir=opt$dir, ga=group, ribp=opt$rib, mitp=opt$mit, odir=qc_dir, mc.cores=mc)
}else if(file_test('-f', opt$dir)){
	ret <- list(my_fun(dir=opt$dir, ga=group, ribp=opt$rib, mitp=opt$mit, odir=qc_dir))
}
all.list <- list()
cell.stat <- data.frame()
for(l in ret){
	cell.stat <- rbind(cell.stat, l$cst)
	if(length(unique(l$data@meta.data$orig.ident)) > 1){
		pbmc.list <- SplitObject(l$data, split.by='orig.ident')
		pbmc.list <- lapply(pbmc.list, function(x){x@meta.data$sample=x@meta.data$orig.ident; return(x)})
		names(pbmc.list) <- unlist(lapply(pbmc.list, function(x){return(unique(x@meta.data$sample))}))
	}else{
		pbmc.list <- list(l$data)
		names(pbmc.list) <- unique(l$data@meta.data$orig.ident)
	}
	all.list <- c(all.list, pbmc.list)
}
all.list[['group.order']] <- group.order
cell.stat <- rbind(cell.stat, data.frame(Sample='Total', cell.raw.count=sum(cell.stat$cell.raw.count), cell.clean.count=sum(cell.stat$cell.clean.count)))
write.table(cell.stat, paste(qc_dir, '/cells_statistics.txt', sep=''), sep='\t', quote=F, row.names=F)
saveRDS(all.list, file = paste(opt$out, '/all_seurat_objects.rds', sep=''))
