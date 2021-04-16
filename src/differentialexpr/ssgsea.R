#!/usr/bin/env Rscript

#dependencies
library(GSEABase)
library(GSVA)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(gplots)
library(textshape)
library(optparse)

option_list = list(
  make_option(c("-e", "--expression"), type="character", default=NULL,
              help="signature reference", metavar="character"),
  make_option(c("-g", "--gmt"), type="character", default=NULL,
              help="pathway database gmt file", metavar="character"),
  make_option(c("-c", "--comparisons"), type="character",
              help="column names of phenotypes want to compare"),
  make_option(c("-m", "--meta"), type="character", default=NULL,
              help="meta info", metavar="character"),
  make_option(c("-n", "--top"), type="character", default=NULL,
              help="top n terms for displaying", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("-s", "--order"), type="character", default=NULL,
              help="sample order for multiqc report", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###parameters
geneSets <- getGmt(opt$gmt)
expr.dat <- read.table(opt$expression, sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
samples <- read.table(file = opt$meta, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
top_n <- opt$top
outdir <- opt$outdir
comparisons <- opt$comparisons
col <- opt$order

#temporarily for multiqc report
meta <- samples

#####read in data
ssgsea <- gsva(as.matrix(expr.dat), geneSets,
               method="ssgsea",
              # rnaseq = TRUE,
               ssgsea.norm=TRUE,
               verbose=TRUE)
###min-max normalization(make scale from -1 to 1)
ssgsea.nor <- t(apply(ssgsea,1,function(x) (2*(x-min(x))/(max(x)-min(x)))-1))

###heatmap for every comparison
for(comparison in unlist(strsplit(comparisons,","))){
  ###get the difference between Clinic traits for each term
  annot.col <- data.frame(Clinic = samples[,comparison],row.names = rownames(samples))
  gr.samples <- lapply(unique(annot.col$Clinic), function(x) return(rownames(annot.col)[which(annot.col$Clinic == x)]))
  ###calculate between-group variation for filtering terms with top variance
  SelTerm <- lapply(c(1:dim(ssgsea.nor)[1]), function(t) {
    grandmean <- mean(ssgsea.nor[t,])
    bgV <- lapply(c(1:length(gr.samples)), function(g){
      bg <- ((ssgsea.nor[t,gr.samples[[g]]] - grandmean)^2)*length(gr.samples[[g]])
      return (bg)
    })
    bgsum <- sum(unlist(bgV))
    return (bgsum)
  })
  names(SelTerm) <- rownames(ssgsea.nor)
  gr.diff.sort <- sort(unlist(SelTerm), decreasing = TRUE)
  ssgsea.sub <- ssgsea.nor[names(gr.diff.sort)[1:top_n],]
}

##Preprocess the data for multiqc report
meta <- meta[order(meta[[col]], decreasing = TRUE),]

#re-order table function
order_ta <- function(x) {
  #col_order
  orde_col <- data.frame("samples" = rownames(meta), num = 1:nrow(meta))
  re_col <- data.frame("samples" = colnames(x), re_num = 1:nrow(meta))

  orde_col <- merge(orde_col, re_col, by = 1)
  orde_col <- orde_col[order(orde_col[["num"]], decreasing = FALSE),]
  t <- x[orde_col[["re_num"]]]

  #row_order
  t$ave <- rowSums(t)
  t <- t[order(t$ave, decreasing = TRUE),]
  t <- t[-ncol(t)]
  return(t)
}

##Extract the heatmap matrix
multiqc_table =  as.data.frame(ssgsea.sub)
multiqc_table <- order_ta(multiqc_table)
rownames(multiqc_table)=gsub("KEGG_","",rownames(multiqc_table))
multiqc_table_cluster= as.data.frame(t(cluster_matrix(t(multiqc_table),dim='col')))
write.table(data.frame("Sample"=rownames(multiqc_table_cluster),multiqc_table_cluster), file = paste(opt$outdir,'ssgsea.txt',sep= ""),quote=FALSE,sep="\t",row.names=FALSE)
