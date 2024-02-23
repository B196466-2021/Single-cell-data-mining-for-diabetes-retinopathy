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
	-p, --pheno	group information file.
	-t, --thread	number of threads (default 1).
	-s, --tissue	tissue in CellMarker file that samples come from (http://biocc.hrbmu.edu.cn/CellMarker/index.jsp#).
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
	'pheno','p',1,'character',
	'thread','t',1,'integer',
	'tissue','s',1,'character',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$pheno) || is.null(opt$tissue) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$thread)) {opt$thread=1}
if (is.null(opt$out)) {opt$out='./'}

##########################
#define functions
my_fun<-function(s){
	##loading matrix
	
	readsCount <- read.table(paste(s, '.txt', sep=''), header = T, row.names = 1)
	sce <- SingleCellExperiment(assays = list(counts = as.matrix(readsCount)))
	
	bcrank <- barcodeRanks(counts(sce))
	uniq <- !duplicated(bcrank$rank)
	
	pdf(paste(opt$out, '/QC/', s, '_UMI_count_depth.pdf', sep=''))
	plot(bcrank@listData$rank[uniq], bcrank@listData$total[uniq], log="xy", xlab="Rank", ylab="Total UMI count", cex=0.5, cex.lab=1.2, col='steelblue3')
	abline(h=bcrank@metadata$inflection, col="darkgreen", lty=2)
	abline(h=bcrank@metadata$knee, col="dodgerblue", lty=2)
	legend("left", legend=c("Inflection", "Knee"), bty="n",col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
	dev.off()
	
	if(any(bcrank@listData$total<200)){
		print('Filtering cells with UMI and FDR...')
		set.seed(100)
		e.out <- emptyDrops(counts(sce), lower=200)
		is.cell <- (e.out$FDR <= 0.01)
		pdf(paste(opt$out, '/QC/', s, '_cells_filter.pdf', sep=''))
		plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"), xlab="Total UMI count", ylab="-Log Probability", cex=0.2)
		abline(v = bcrank@metadata$inflection, col="darkgreen")
		abline(v = bcrank@metadata$knee, col="dodgerblue")
		legend("bottomright", legend=c("Inflection", "Knee"), bty="n", col=c("darkgreen", "dodgerblue"), lty=1, cex=1.2)
		dev.off()
		w2kp <- which(is.cell & e.out$Total >= bcrank@metadata$inflection)
		sce <- sce[, w2kp]
	}
	
	##Mitochondrial or Ribosomal genes
	ribo.file <- '/opt/Genome/GRCh38/Annotation/ribosomal_gene.txt'
	if(file.exists(ribo.file)){
		ribo <- read.table(ribo.file, header=T, sep='\t', stringsAsFactors=F)
		riboh <- which(row.names(counts(sce)) %in% ribo$Approved.symbol)
	}else{
		riboh <- c()
	}
	
	is.mito <- grep('^MT-', row.names(counts(sce)))
	is.ribo <- grep('^RPL|^MRPL|^MRPS|^RPS', row.names(counts(sce)), ignore.case=T, perl=T)
	is.ribo <- unique(c(is.ribo, riboh))
	
	sce = calculateQCMetrics(sce, feature_controls=list(Mt=is.mito, Ri=is.ribo))
	
	pdf(paste(opt$out, '/QC/', s, '_counts_in_feature.pdf', sep=''))
	par(mfrow=c(2,2), mar=c(5, 4, 1, 1), bty="n")
	hist(log10(sce$total_counts), xlab="log10(Library sizes)", main="", breaks=20, col="grey80", ylab="Number of cells")
	hist(log10(sce$total_features_by_counts), xlab="log10(# of expressed genes)", main="", breaks=20, col="grey80", ylab="Number of cells")
	hist(sce$pct_counts_Ri, xlab="Ribosome prop. (%)", ylab="Number of cells", breaks=40, main="", col="grey80")
	hist(sce$pct_counts_Mt, xlab="Mitochondrial prop. (%)", ylab="Number of cells", breaks=80, main="", col="grey80")
	dev.off()
	
	pdf(paste(opt$out, '/QC/', s, '_counts_vs_feature.pdf', sep=''))
	par(mfrow=c(2, 2), mar=c(5, 4, 1, 1), bty="n")
	smoothScatter(log10(sce$total_counts), log10(sce$total_features_by_counts),  xlab="log10(Library sizes)", ylab="log10(# of expressed genes)", nrpoints=500, cex=0.5)
	smoothScatter(log10(sce$total_counts), sce$pct_counts_Ri, xlab="log10(Library sizes)", ylab="Ribosome prop. (%)", nrpoints=500, cex=0.5)
	abline(h=10, lty=1)
	smoothScatter(log10(sce$total_counts), sce$pct_counts_Mt, xlab="log10(Library sizes)", ylab="Mitochondrial prop. (%)", nrpoints=500, cex=0.5)
	abline(h=5, lty=1)
	smoothScatter(sce$pct_counts_Ri, sce$pct_counts_Mt, xlab="Ribosome prop. (%)", ylab="Mitochondrial prop. (%)", nrpoints=500, cex=0.5)
	abline(h=5,  lty=1)
	abline(v=10, lty=1)
	dev.off()
	
	sce <- sce[,which(sce$pct_counts_Mt <= 10 | sce$pct_counts_Ri >= 10)]
	
	##Summarize gene-level information
	pdf(paste(opt$out, '/QC/', s, '_UMI_vs_cells.pdf', sep=''), width=20)
	par(mfrow=c(1,3), mar=c(5,4,1,1))
	hist(log10(rowData(sce)$mean_counts+1e-6), col="grey80",  main="", breaks=40, xlab="log10(ave # of UMI + 1e-6)")
	hist(log10(rowData(sce)$n_cells_by_counts+1), col="grey80", main="", breaks=40, xlab="log10(# of expressed cells + 1)")
	smoothScatter(log10(rowData(sce)$mean_counts+1e-6), log10(rowData(sce)$n_cells_by_counts + 1), xlab="log10(ave # of UMI + 1e-6)", ylab="log10(# of expressed cells + 1)")
	dev.off()
	od1 = order(rowData(sce)$mean_counts, decreasing = TRUE)
	sgd <- data.frame(gene=rownames(rowData(sce))[od1[20:1]], exp=rowData(sce)$mean_counts[od1[20:1]], ncell=rowData(sce)$n_cells_by_counts[od1[20:1]])
	sgd$gene <- factor(sgd$gene, levels=rownames(rowData(sce))[od1[20:1]])
	p <- ggplot(sgd, aes(gene, exp, fill=exp))+geom_col(show.legend=F)+coord_flip()
	p <- p+geom_text(aes(y=exp+(max(exp*0.03)), label=ncell))+theme_bw()
	p <- p+labs(x='', y='ave # of UMI')
	p <- p+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
	ggsave(paste(opt$out, '/QC/', s, '_highly_expressed_genes.pdf', sep=''), p, width=9, device="pdf")
	
	##Normalization
	pbmc <- CreateSeuratObject(counts = counts(sce), project = s, min.cells = 3, min.features = 200)
	pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
	return(pbmc)
}

fdmk <- function(x){
	markers <- FindMarkers(pbmc, ident.1 = x, verbose = FALSE, test.use='MAST', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	markers$cluster <- x
	markers$gene <- rownames(markers)
	markers <- dplyr::filter(markers, p_val<=0.05)
	return(markers)
}

depf <- function(x){
	diff_gene <- FindMarkers(pbmc2, ident.1 = x, verbose = FALSE, test.use='MAST', min.pct = 0.25, logfc.threshold = 0)
	od <- data.frame(Gene=rownames(diff_gene), log2FC=diff_gene$avg_logFC/log(2), p_value=diff_gene$p_val, q_value=diff_gene$p_val_adj)
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
library(scran)
library(pheatmap)
library(gplots)
library(MAST)
library(monocle)
library(parallel)

group <- read.table(opt$pheno, header=T, sep='\t', stringsAsFactors=F)
gr <- unlist(strsplit(as.character(group$Group[1]), ","))
ga <- data.frame(rep(group$Sample[1], length(gr)), gr, stringsAsFactors=F)
names(ga) <- names(group)
if(nrow(group)>=2){
	for(i in 2:nrow(group)){
		gr <- unlist(strsplit(as.character(group$Group[i]), ","))
        g <- data.frame(rep(group$Sample[i], length(gr)), gr, stringsAsFactors=F)
		names(g) <- names(group)
		ga <- rbind(ga,g)
	}
}
qc_dir <- paste(opt$out, '/QC', sep='')
if(! dir.exists(qc_dir)){
	dir.create(qc_dir, recursive=T)
}

group.list <- list()
for(g in unique(ga$Group)){
	g_dir <- paste(opt$out, '/', g, sep='')
	if(! dir.exists(g_dir)){
		dir.create(g_dir, recursive=T)
	}
	samples <- ga$Sample[ga$Group==g]
	mc <- getOption("mc.cores", opt$thread)
	pbmcs <- mclapply(samples, my_fun, mc.cores = mc)
	pbmc <- merge(pbmcs[[1]], pbmcs[-1])

	##feature selection
	pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
	top10 <- head(VariableFeatures(pbmc), 10)
	plot1 <- VariableFeaturePlot(pbmc)
	plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	ggsave(paste(g_dir, '/feature_selection.pdf', sep=''), plot2, width=9, device="pdf")

	##Scaling the data
	all.genes <- rownames(pbmc)
	pbmc <- ScaleData(pbmc, features = all.genes)
	
	##Dimension reduction
	###PCA
	pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F)
	p <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
	ggsave(paste(g_dir, '/top_dimensional_reduction_genes.pdf', sep=''), p, width=8, device="pdf")
	pcd <- data.frame(pbmc@reductions$pca@cell.embeddings)[, 1:2]
	pcd$log10_features <- log10(pbmc@meta.data$nFeature_RNA)
	names(pcd)[1:2] <- c('PC1', 'PC2')
	p <- ggplot(pcd, aes(PC1, PC2, col=log10_features))
	p <- p+geom_point(size=0.2,alpha=0.6) + theme_classic()
	p <- p+scale_colour_gradient(low="lightblue",high="red")
	p <- p+guides(color = guide_legend(override.aes = list(size=3)))
	ggsave(paste(g_dir, '/PCA.pdf', sep=''), p, width=9, device="pdf")
	
	###Heatmap
	pdf(paste(g_dir, '/heatmap_PC1.pdf', sep=''))
	DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
	dev.off()
	
	pdf(paste(g_dir, '/heatmap_PC2.pdf', sep=''))
	DimHeatmap(pbmc, dims = 2, cells = 500, balanced = TRUE)
	dev.off()
	
	pdf(paste(g_dir, '/heatmap_PC3.pdf', sep=''))
	DimHeatmap(pbmc, dims = 3, cells = 500, balanced = TRUE)
	dev.off()
	
	pdf(paste(g_dir, '/heatmap_PC1_9.pdf', sep=''))
	DimHeatmap(pbmc, dims = 1:9, cells = 500, balanced = TRUE)
	dev.off()
	
	###Select PCs
	pbmc <- JackStraw(pbmc, num.replicate = 50, dims =50)
	pbmc <- ScoreJackStraw(pbmc, dims = 1:50)
	pc_pv <- as.data.frame(pbmc@reductions$pca@jackstraw@overall.p.values)
	npc1 <- max(which(pc_pv$Score < 0.01))
	npc2 <- min(which(cumsum(pbmc@reductions$pca@stdev)/sum(pbmc@reductions$pca@stdev) > 0.7))
	npc <- min(npc1, npc2)
	p <- ElbowPlot(pbmc, ndims =50)
	p <- p+geom_vline(xintercept = npc, colour = "red", linetype=2)
	ggsave(paste(g_dir, '/deviation_vs_PC.pdf', sep=''), p, width=9, device="pdf")
	
	##Cluster the cells
	pbmc <- FindNeighbors(pbmc, dims = 1:npc)
	pbmc <- FindClusters(pbmc, resolution = 0.5)
	###TSNE
	pbmc <- RunTSNE(pbmc, dims = 1:npc)
	p <- DimPlot(pbmc, reduction = "tsne", label=T, group.by='seurat_clusters')
	ggsave(paste(g_dir, '/TSNE.pdf', sep=''), p, width=9, device="pdf")
	
	#saveRDS(pbmc, file = paste(opt$out, "pbmc_tutorial.rds", sep=''))
	###ClusterTree
	pbmc <- BuildClusterTree(pbmc)
	pdf(paste(g_dir, '/cluster_tree.pdf', sep=''))
	PlotClusterTree(pbmc)
	dev.off()
	
	##cluster biomarkers
	#pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,  return.thresh=0.05, test.use='MAST')
	mc <- getOption("mc.cores", opt$thread)
	maker.list <- mclapply(levels(pbmc@active.ident), fdmk, mc.cores = mc)
	maker_gene <- do.call(rbind, maker.list)
	top1_maker_gene <- maker_gene %>% dplyr::group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
	ncs <- ceiling(sqrt(npc))
	nrs <- ceiling(npc/ncs)
	p <- FeaturePlot(pbmc, features =top1_maker_gene$gene, ncol=ncs)
	if(npc>3){
	w <- ncs*4
	h <- nrs*3
	}else{
		h <- 3
		w <- npc*3
	}

	ggsave(paste(g_dir, '/top1_marker_gene_exp.pdf', sep=''), p, width=w, height=h, device="pdf")
	
	top10 <- maker_gene %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
	
	if(npc <=7){
		w <- 7
		h <- 7
	}else{
		w <- npc
		h <- npc
	}
	p <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
	ggsave(paste(g_dir, '/top10_marker_gene_heatmap.pdf', sep=''), p, width=w, height=h, device="pdf")

	##Assigning cell type identity to clusters
	write.table(maker_gene, file=paste(g_dir, '/maker_gene.txt', sep=''), quote=F, row.names=F, sep='\t')
	if(grepl('tumor|cancer', g, ignore.case=T)){
		sct <- '-c'
	}else{
		sct <- ''
	}
	comd <- paste('Rscript /home/dev/DEV-cuiy/cellTypeEnrich.r -f=', g_dir, '/maker_gene.txt -t=', opt$tissue, ' -o=', g, '_marker_anno ', sct, sep='')
	print(comd)
	system(comd)
	makeh <- read.table(paste(g, '_marker_anno.txt', sep=''), header=T, sep='\t')[, c('Cell_Name', 'Cluster', 'p_value')]
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
	ggsave(paste(g_dir, '/TSNE_with_cell_type.pdf', sep=''), p, device="pdf", width=9)

	pbmc2@meta.data$group<-rep(g, nrow(pbmc2@meta.data))
	saveRDS(pbmc2, file = paste(g_dir, '/', g, "_final.rds", sep=''))
	group.list[[g]] <- pbmc2
	##Differential Expression Analysis
	#diff.gene <- FindAllMarkers(pbmc2, min.pct = 0.25, return.thresh=1, logfc.threshold=0, test.use='MAST')
	#for(ct in unique(diff.gene$cluster)){
	#	sudg <- data.frame(Gene=diff.gene$gene, log2FC=diff.gene$avg_logFC/log(2), p_vlaue=diff.gene$p_val, q_vlaue=diff.gene$p_val_adj)[diff.gene$cluster==ct, ]
	#	write.table(sudg, file=paste(g_dir, '/', gsub(' ', '_', ct), '_vs_other_cells.txt', sep=''), quote=F, row.names=F, sep='\t')
	#}
	mc <- getOption("mc.cores", opt$thread)
	mclapply(levels(pbmc2@active.ident), depf, mc.cores = mc)
	##Assign Cell-Cycle Scores
	cc_dir <- paste(g_dir, '/cell_cycle/', sep='')
	if(!dir.exists(cc_dir)){
		dir.create(cc_dir, recursive=T)
	}
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes
	pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
	ccgs <- rownames(pbmc@assays$RNA@scale.data)[toupper(rownames(pbmc@assays$RNA@scale.data)) %in% c(s.genes, g2m.genes)]
	ncs <- ceiling(sqrt(length(ccgs)))
	nrs <- ceiling(length(ccgs)/ncs)
	if(length(ccgs)>3){
		w <- ncs*3
		h <- nrs*3
	}else{
		h <- 3
		w <- length(ccgs)*3
	}
	p <- RidgePlot(pbmc, features = ccgs, ncol = ncs)
	ggsave(paste(g_dir, '/cell_cycle_ridge_plot.pdf', sep=''), p, device="pdf", width=w, height=h)
	
	pdcc <- data.frame(cluster=pbmc[[]]$seurat_clusters, phase=pbmc[[]]$Phase, check.names=F)
	p <- ggplot(pdcc, aes(cluster, fill=phase))
	p <- p+geom_bar(position = "dodge2")
	p <- p+theme_classic()+guides(fill = guide_legend(title = "Cell Cycel Phase"))
	p <- p+labs(x='Clusters', y='Cell Count')
	
	ggsave(paste(cc_dir, '/cell_cycle_summary.pdf', sep=''), p, device="pdf", width=8)
	
	scd <- pbmc@assays$RNA@scale.data
	top10scd <- scd[top10$gene, ]
	anno_col <- pbmc@meta.data[, c('S.Score', 'G2M.Score', 'Phase', 'seurat_clusters')]
	names(anno_col)[4] <- 'Cluster'
	pda <- top10scd[, rownames(anno_col)[order(anno_col$Cluster)]]
	pheatmap(pda, annotation_col = anno_col, show_colnames=F, cluster_cols=F, cluster_rows=F, fontsize_row=nrow(pda)/50, color = colorpanel(128, 'blue3', 'yellow', 'red3'), legend=F, filename=paste(cc_dir, '/top10_hetmap_cell_cycle.pdf', sep=''), width=(nrow(pda)/20)-1, height=nrow(pda)/20)

	marrow <- RunPCA(pbmc, features = c(s.genes, g2m.genes), verbose=F)
	p <- DimPlot(marrow, reduction='pca')
	ggsave(paste(g_dir, '/cell_cycle_PCA.pdf', sep=''), p, device="pdf", width=8)
	
	mo <- merge(pbmc@meta.data, data.frame(pbmc2@active.ident), by=0)
	wo <- mo[c('Row.names', 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'seurat_clusters', 'S.Score', 'G2M.Score', 'Phase', 'pbmc2.active.ident')]
	names(wo)<-c('Cell_Id', 'Sample', 'nCount_RNA', 'nFeature_RNA', 'Clusters', 'S_Score', 'G2M_Score', 'Phase', 'Cell_Type')
	write.table(wo, file=paste(cc_dir, '/cell_sample_informations.txt', sep=''), quote=F, row.names=F, sep='\t')

	##Constructing Single Cell Trajectories
	tjdir <- paste(g_dir, '/trajectory/', sep='')
	if(!dir.exists(tjdir)){
		dir.create(tjdir, recursive=T)
	}
	cds <- as.CellDataSet(pbmc)
	cds <- estimateSizeFactors(cds)
	cds <- estimateDispersions(cds)
	cds <- detectGenes(cds, min_expr = 0.1)
	expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))	
	ctd <- data.frame(pbmc2@active.ident)
	names(ctd) <- 'cell'
	tjtrf <- function(cl){
		subcds <- cds[, rownames(ctd)[ctd$cell==cl]]
		subcds <- detectGenes(subcds, min_expr = 0.1)
		fData(subcds)$use_for_ordering <- fData(subcds)$num_cells_expressed > 0.05 * ncol(subcds)
		subcds <- reduceDimension(subcds, max_components = 2, norm_method = 'log', num_dim = 3, reduction_method = 'tSNE', verbose = F, perplexity =10)
		subcds <- clusterCells(subcds, verbose = F)
		#subcds <- clusterCells(subcds, rho_threshold = 2, delta_threshold = 4, skip_rho_sigma = T, verbose = F)
		clustering_DEG_genes <- differentialGeneTest(subcds[expressed_genes, ], fullModelFormulaStr = '~Cluster', cores=opt$thread)
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
		write.table(odm, file=paste(tjdir, '/', gsub(' ', '_', as.character(cl)), '_single_cell_trajectories.txt', sep=''), quote=F, row.names=F, sep='\t')
		p1 <- plot_cell_trajectory(subcds, color_by = "State")
		p2 <- plot_cell_trajectory(subcds, color_by = "Pseudotime")
		pdf(paste(tjdir, '/', gsub(' ', '_', as.character(cl)), '_single_cell_trajectories.pdf', sep=''), width=13)
		multiplot(p1, p2, cols=2)
		dev.off()
	}
	for(cl in unique(ctd$cell)){
		tjtrf(cl)
	}
	#mc <- getOption("mc.cores", opt$thread)
	#mclapply(unique(ctd$cell), tjtrf, mc.cores = mc)

		#if(length(unique(subcds@phenoData@data$State))==1){
		#	subcds <- orderCells(subcds, root_state=1)
		#}
		
		#BEAM_res <- BEAM(subcds, branch_point = 1, cores = opt$thread)
		#BEAM_res <- BEAM_res[order(BEAM_res$qval),]
		#BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
		#BEAM_res <- BEAM_res %>% top_n(n = 100, wt = -pval)
		#h <- nrow(BEAM_res)/10
		#if(h <7){
		#	h <- 7
		#}
		#if(nrow(BEAM_res)>1){
		#	pdf(paste(g_dir, '/', gsub(' ', '_', as.character(cl)), '_branches_in_single_cell_trajectories.pdf', sep=''), height=h, width=9)
		#	plot_genes_branched_heatmap(subcds[BEAM_res$gene_short_name,], branch_point = 1, show_rownames = T, cores=opt$thread, use_gene_short_name = T)
		#	dev.off()
		#}
		#top10BEAM_res <- BEAM_res %>% top_n(n = 10, wt = -pval)
}

saveRDS(group.list, file = paste(opt$out, '/all_seurat_objects.rds', sep=''))
