########################################
##! @Author: Yi Cui
##! @Todo: Gene Set Cell Type Enrichment Analysis
##! @Version: 1.0.0
##! @Dep: R
##! @ChangLog:
##!   Latest version:
##!   Older versions:
########################################

########################################
## Functions
usage<-function(spec){
	cat("This script is used to\n",
		getopt(spec,usage=TRUE),
"Options:
	-f, --file	input file cell marker file (cluster<tab>gene).
	-d, --dat	gene set file with gmt format.
	-o, --out	prefix of out files (default Anno_Cell).
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
				'dat','d',1,'character',
				'out','o',1,'character',
				'version','v',0,'logical',
				'help','h',0,'logical'
				),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || is.null(opt$dat) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$out)) {opt$out='Anno_Cell'}


inverseList <- function(ilist){
  rId <- unlist(ilist, use.names = FALSE)
  lId <- rep(names(ilist), sapply(ilist, length))
  return(split(lId, rId))
}

#  buildList(funcPar$type, funcPar$cType)
buildList <- function(){
	dblines <- readLines(opt$dat)
	cell2Gene <- list()
	for (i in dblines){
		dbline <- unlist(strsplit(i, split = "\t"))
		cell2Gene[[dbline[1]]] <- dbline[-c(1, 2)]
	}
	return(cell2Gene)
}

### cell enrich
# results <- GOProc(dbList, geneInt, g, funcPar$pValue, funcPar$mCount, funcPar$Lev)
enrichProc <- function(dbList, geneInt, cluster){
  cell2Gene <- dbList
  gene2Cell <- inverseList(cell2Gene)
  intgene2Cell <- gene2Cell[names(gene2Cell) %in% geneInt]
  if (length(intgene2Cell) >= 1){
    lenInt <- length(intgene2Cell)
    intcell2Gene <- inverseList(intgene2Cell)
    ## filter enrichment gene number
    intcell2Gene <- intcell2Gene[sapply(intcell2Gene, length) > 0]
    intgene2Cell <- inverseList(intcell2Gene)
    results <- data.frame(Cell_Name=names(intcell2Gene), Cluster=cluster, stringsAsFactors=F)
    ## Build value lists
    symlist <- sapply(intcell2Gene, function(ids){
                      paste(ids, collapse="//")
                      })
    results$Count = 0
    results$Size = 0
    for (i in c(1:dim(results)[1])){
      tmpCName <- results$Cell_Name[i]
      results[i, ]$Count = length(unlist(intcell2Gene[tmpCName], use.names = F))
      results[i, ]$Size = length(unlist(cell2Gene[tmpCName], use.names = F))
    }
    results$numInt <- lenInt
    results$numTotal <- length(gene2Cell)

    results$p_value <- apply(results[, 3:6], 1, function(li){
                             fisher.test(cbind(matrix(c(li["Count"],
                             li["Size"] - li["Count"],
                             li["numInt"] - li["Count"],
                             li["numTotal"] - li["numInt"] - li["Size"]
                             + li["Count"]), nrow = 2)),
                             alternative="greater")$p.value
                             })
    results$p_value[results$p_value < 1e-299] <- 1e-299
    results$fdr <- p.adjust(results$p_value, "fdr")
    results$Enrichment_Score <- -log10(results$p_value)
    #results$Fold_Enrichment <- (results$Count / results$Size) / (results$numInt / results$numTotal)
    #results$GeneRatio <- results$Count / results$numInt
    results$Genes <- symlist
    # results <- results[results$p_value <= pVal, ]
    # results <- results[results$Count >= minCount, ]
    results <- results[order(results$p_value, decreasing = F), ]
  } else {
    results <- NULL
  }
  return(results)
}

#########################################
### Main
library(stringr)
expTable <- read.table(opt$file, stringsAsFactors = F, header = T,
                       sep = "\t", check.names = F)
expTable <- expTable[, c("gene", "cluster")]
# if (! "gene" %in% colnames(expTable) && ! "cluster" %in% colnames(expTable)){
# }
dbList <- buildList()
results <- NULL
for (i in unique(expTable$cluster)){
  geneInt <- expTable$gene[expTable$cluster == i]
  results <- rbind(results, enrichProc(dbList, geneInt, i))
}
write.table(results, paste(opt$out, ".txt", sep=""), sep = "\t", row.names = FALSE, qmethod = "escape", quote = FALSE)
system(paste("python -c \"from Yilib import table2excel; table2excel(['", opt$out, ".txt'], '", opt$out, ".xlsx')\"", sep = ""))

print(paste("[INFO]", opt$file, "Cell Name Annotation done", sep = " "))

