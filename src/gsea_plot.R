#!/usr/bin/env Rscript

#dependencies
library(ggplot2)
library(ggpubr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(fgsea)
library(optparse)
library(tidyverse)

option_list = list(
  make_option(c("-m", "--deseq2_mat"), type="character", default=NULL, 
              help="signature reference", metavar="character"),
  make_option(c("-p", "--pcut"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-s", "--species"), type="character", default=NULL, 
              help="species", metavar="character"),
  make_option(c("-n", "--npermutation"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

####read in data and convert to entrez ID
data <- read.table(opt$deseq2_mat, sep = "\t", header = TRUE, row.names = 1)
# raw <- read.table("~/Documents/rnaseq/data/rsem_DESeq2.txt", sep = "\t", header = TRUE, row.names = 1)
# data <- read.table("~/Documents/rnaseq/data/DESeq2_ConvertID.txt", sep = "\t", header = TRUE, row.names = 1)
data <- na.omit(data)
geneList <- sign(data$log2FoldChange) * (-log10(data$pvalue))
geneList <- sort(geneList, decreasing = TRUE)
if(opt$species == "hg38") {
  GeneIDSymbol <- toTable(org.Hs.egSYMBOL)
  dataset <- org.Hs.eg.db
  sp <- "hsa"
} else {
  GeneIDSymbol <- toTable(org.Mm.egSYMBOL)
  dataset <- org.Mm.eg.db
  sp <- "mmu"
}
names(geneList) <- GeneIDSymbol[match(data[,1],GeneIDSymbol$symbol),'gene_id']
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]

####function of GO enrichment
GSEAGO <- function(geneList,ont,minGSSize,nPerm,pcut){
  print(ont)
  set.seed(1234)
  gsea.go <- gseGO(geneList = geneList, 
                   OrgDb = dataset, 
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
                       organism = sp, 
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

####function of plotting enrichment results
GSEAPlot <- function(res,category){
  nes.pos <- length(which(res$NES > 0))
  nes.neg <- length(which(res$NES < 0))
  min.nes <- -ceiling(abs(min(res$NES[which(res$NES < 0)])))-1
  max.nes <- ceiling(max(res$NES[which(res$NES > 0)]))+1
  if(nes.pos == 0){
    max.nes <- -min.nes
  }
  if(nes.neg == 0){
    min.nes <- -max.nes
  }
  gp <- ggplot(res,aes(x=Description,y=NES,fill=as.factor(group))) + 
    xlab("Enrichment Terms") +
    ylab("Normalized Enrichment Score")+
    ylim(min.nes, max.nes)+
    ggtitle(category)+
    geom_bar(stat='identity') +
    theme_bw()+
    theme(#axis.text.x=element_blank(),
      plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x=element_text(size=12,face = "bold",hjust=1),
      axis.title.x = element_text(size = 12,face = "bold"),
      axis.title.y = element_text(size = 12,face = "bold"),
      legend.position='bottom',
      panel.background = element_rect(colour = "black", size=1))+
    scale_fill_manual(name = "", values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8"))
  
  if(nes.pos == 0 || nes.neg > 0){
    gp <- gp + annotate("text",label=res$Description[which(res$NES < 0)],
                        x=c(1:nes.neg),
                        y=rep(0,nes.neg),size=tsize,hjust=0)   #negative score
  }
  if(nes.neg == 0 || nes.pos > 0){
    gp <- gp + annotate("text",label=res$Description[which(res$NES > 0)],
                        x=c((nes.neg+1):dim(res)[1]),
                        y=rep(0,nes.pos),size=tsize,hjust=1)  #positive score
  }
  gp <- gp + annotate("text",label=res$significance,
                      x=c(1:dim(res)[1]),
                      y=c(rep(min.nes,nes.neg),rep(max.nes,nes.pos)),size=msize,color = "black")+
    coord_flip(expand = TRUE)
  return (gp)
}


####functional enrichment
pcut <- as.numeric(opt$pcut)
minGSSize <- 5
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


###To plot
#pdf(paste(opt$outdir,"diff_gsea_plot.pdf", sep = ""),width = 14, height = 10)
#png("~/Documents/rnaseq/data/diff_gsea_plot.png", res = 320, width = 4700, height = 3000)
png(paste(opt$outdir,"_diff_gsea_plot.png", sep = ""), res = 320, width = 4700, height = 7000)
#g1 <- GSEAPlot(go.mf,"Molecular Function")
#g2 <- GSEAPlot(go.bp,"Biological Process")
#g3 <- GSEAPlot(go.cc,"Cellular Component")
#g4 <- GSEAPlot(kegg,"KEGG")
#ggarrange(g1, g2, g3, g4,ncol = 2, nrow = 2)

if(dim(go.mf)[1] != 0) {
GSEAPlot(go.mf,"Molecular Function")
}
if(dim(go.bp)[1] != 0) {
GSEAPlot(go.bp,"Biological Process")
}

if(dim(go.cc)[1] != 0) {
GSEAPlot(go.cc,"Cellular Component")
}

if(dim(kegg)[1] != 0) {
GSEAPlot(kegg,"KEGG")
}

dev.off()
