#!/usr/bin/env Rscript

#dependencies
library(ggplot2)
library(ggpubr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(optparse)
library(tidyverse)

option_list = list(
  make_option(c("-m", "--deseq2_mat"), type="character", default=NULL,
              help="signature reference", metavar="character"),
  make_option(c("-p", "--pcut"), type="character", default=NULL,
              help="pcut value", metavar="character"),
  make_option(c("-s", "--minsize"), type="character", default=NULL,
              help="size ", metavar="character"),
  make_option(c("-n", "--npermutation"), type="character", default=NULL,
              help="number of permutations", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

####read in data and convert to entrez ID
data <- read.table(opt$deseq2_mat, sep = "\t", header = TRUE, row.names = 1)
data <- na.omit(data)
geneList <- sign(data$log2FoldChange) * (-log10(data$pvalue))
GeneIDSymbol <- toTable(org.Hs.egSYMBOL)
names(geneList) <- GeneIDSymbol[match(data[,1],GeneIDSymbol$symbol),'gene_id']
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]

####function of GO enrichment
GSEAGO <- function(geneList,ont,minGSSize,nPerm,pcut){
  print(ont)
  set.seed(1234)
  gsea.go <- gseGO(geneList = geneList,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = ont,
                   nPerm = nPerm, # number permutations
                   minGSSize = minGSSize,
                   pAdjustMethod = "BH",
                   pvalueCutoff = pcut, # padj cutoff value
                   verbose = FALSE,
                   seed = TRUE)
  go.res <- gsea.go@result %>%
    mutate(group = ifelse(NES > 0,"Up-regulated","Down-regulated")) %>%
    mutate(significance = cut(p.adjust,breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))) %>%
    dplyr::filter(p.adjust < pcut) %>%
    arrange(desc(NES))
    return(go.res)
  res <- NULL
  if(dim(subset(go.res, NES > 0))[1] >= 10){
    res <- rbind(res,go.res[1:10,])
  }
  if(dim(subset(go.res, NES > 0))[1] < 10){
    res <- rbind(res,subset(go.res, NES > 0))
  }
  if(dim(subset(go.res, NES < 0))[1] >= 10){
    total <- dim(go.res)[1]
    res <- rbind(res,go.res[(total-9):total,])
  }
  if(dim(subset(go.res, NES < 0))[1] < 10){
    res <- rbind(res,subset(go.res, NES < 0))
  }
  res$Description <- factor(res$Description,levels = levels(reorder(res$Description,res$NES)))
  res.new <- res[match(levels(res$Description),res$Description),]
  return(res.new)
}

####function of KEGG pathway enrichment
GSEAKEGG <- function(geneList,minGSSize,nPerm,pcut){
  set.seed(1234)
  gsea.kegg <- gseKEGG(geneList = geneList,
                       organism = "hsa",
                       keyType = "kegg",
                       nPerm = nPerm, # number permutations
                       minGSSize = minGSSize,
                       pAdjustMethod = "BH",
                       pvalueCutoff = pcut, # padj cutoff value
                       verbose = FALSE,
                       seed = TRUE)
  kegg.res <- gsea.kegg@result %>%
    mutate(group = ifelse(NES > 0,"Up-regulated","Down-regulated")) %>%
    mutate(significance = cut(p.adjust,breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))) %>%
    dplyr::filter(p.adjust <= pcut)
    #return(kegg.res)
  res <- NULL
  if(dim(subset(kegg.res, NES > 0))[1] >= 10){
    res <- rbind(res,kegg.res[1:10,])
  }
  if(dim(subset(kegg.res, NES > 0))[1] < 10){
    res <- rbind(res,subset(kegg.res, NES > 0))
  }
  if(dim(subset(kegg.res, NES < 0))[1] >= 10){
    total <- dim(kegg.res)[1]
    res <- rbind(res,kegg.res[(total-9):total,])
  }
  if(dim(subset(kegg.res, NES < 0))[1] < 10){
    res <- rbind(res,subset(kegg.res, NES < 0))
  }
  res$Description <- factor(res$Description,levels = levels(reorder(res$Description,res$NES)))
  res.new <- res[match(levels(res$Description),res$Description),]
  return(res.new)
}

####functional enrichment
pcut <- as.numeric(opt$pcut)
minGSSize <- as.numeric(opt$minsize)
nPerm <- as.numeric(opt$npermutation)
tsize <- 4   ###term size in plot
msize <- 5   ###significance label size in plot
###For GO
go.mf <- GSEAGO(geneList,"MF",minGSSize,nPerm,pcut)
go.bp <- GSEAGO(geneList,"BP",minGSSize,nPerm,pcut)
go.cc <- GSEAGO(geneList,"CC",minGSSize,nPerm,pcut)
###For KEGG
kegg <- GSEAKEGG(geneList,minGSSize,nPerm,pcut)
###write enriched term table
write.table(go.mf, paste(opt$outdir,"_GO_MF_terms.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(go.bp, paste(opt$outdir,"_GO_BP_terms.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(go.cc, paste(opt$outdir,"_GO_CC_terms.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(kegg, paste(opt$outdir,"_KEGG_terms.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

#############multiqc file
go.mf.multiqc <- subset(go.mf,  select = c(Description,NES,p.adjust))
colnames(go.mf.multiqc) <- c("Description","MF","Pvalue")
go.mf.multiqc  <- go.mf.multiqc[order(-go.mf.multiqc$Pvalue),]
go.mf.multiqc <- subset(go.mf.multiqc[1:10,],select=c(Description,MF))
write.table(go.mf.multiqc,paste(opt$outdir,"_GO_MF_report.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

go.bp.multiqc <- subset(go.bp,  select = c(Description,NES,p.adjust))
colnames(go.bp.multiqc) <- c("Description", "BP","Pvalue")
go.bp.multiqc  <- go.bp.multiqc[order(-go.bp.multiqc$Pvalue),]
go.bp.multiqc <- subset(go.bp.multiqc[1:10,],select=c(Description,BP))
write.table(go.bp.multiqc,paste(opt$outdir,"_GO_BP_report.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

go.cc.multiqc <- subset(go.cc,  select = c(Description,NES,p.adjust))
colnames(go.cc.multiqc) <- c("Description", "CC", "Pvalue")
go.cc.multiqc <- go.cc.multiqc[order(-go.cc.multiqc$Pvalue),]
go.cc.multiqc <- subset(go.cc.multiqc[1:10,],select=c(Description,CC))
write.table(go.cc.multiqc,paste(opt$outdir,"_GO_CC_report.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

go.kegg.multiqc <- subset( kegg,  select = c(Description,NES,p.adjust))
colnames(go.kegg.multiqc) <- c("Description", "KEGG","Pvalue")
go.kegg.multiqc <- go.kegg.multiqc[order(-go.kegg.multiqc$Pvalue),]
go.kegg.multiqc <- subset(go.kegg.multiqc[1:10,],select=c(Description,KEGG))
write.table(go.kegg.multiqc,paste(opt$outdir,"_GO_KEGG_report.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)