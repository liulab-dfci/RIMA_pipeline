#!/usr/bin/env Rscript

#dependencies
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(fgsea))
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))

option_list = list(
  make_option(c("-m", "--deseq2_mat"), type="character", default=NULL,
              help="signature reference", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("-t", "--treatment"), type="character", default=NULL,
              help="meta info", metavar="character"),
  make_option(c("-c", "--control"), type="character", default=NULL,
              help="control", metavar="character"),
  make_option(c("-g", "--condition"), type="character", default=NULL,
              help="condition", metavar="character"),
  make_option(c("-p", "--pcut"), type="character", default=NULL,
              help="p value cutoff", metavar="character"),
  make_option(c("-n", "--npermutation"), type="character", default=NULL,
              help="number of permutation", metavar="character"),
  make_option(c("-s", "--minsize"), type="character", default=NULL,
              help="minimal gsea size", metavar="character"),
  make_option(c("-l", "--hallmark"), type="character", default=NULL,
              help="hallmark gmt file", metavar="character")
); 

###parameters
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
Outdir <- opt$outdir
h_file <- opt$hallmark
Condition <- opt$condition
Treatment <- opt$treatment
Control <- opt$control
pcut <- as.numeric(opt$pcut)
minGSSize <- as.numeric(opt$minsize)
nPerm <- as.numeric(opt$npermutation)
#tsize <- 4   ###term size in plot
#msize <- 5   ###significance label size in plot


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
                   
  return (gsea.go) 
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
 return (gsea.kegg)
}

####function of mSigDB hallmark geneset enrichment
GSEAHALLMARK <- function(geneList,minGSSize,pcut,hallmark){
  set.seed(1234)
  gsea.hallmark <- GSEA(geneList = geneList,
                       TERM2GENE = hallmark, 
                       minGSSize = minGSSize,
                       pAdjustMethod = "BH",
                       pvalueCutoff = pcut, # padj cutoff value
                       verbose = FALSE,
                       seed = TRUE)
 return (gsea.hallmark)
}


format <- function(gseaRes){  
  enrich.res <- gseaRes@result %>%
    mutate(group = ifelse(NES > 0,paste("Enrich in",Treatment) ,paste("Enrich in",Control))) %>%
    mutate(significance = cut(p.adjust,breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))) %>%
    dplyr::filter(p.adjust <= pcut)

  res <- NULL
  if(dim(subset(enrich.res, NES > 0))[1] >= 10){
    res <- rbind(res,enrich.res[1:10,])
  }
  if(dim(subset(enrich.res, NES > 0))[1] < 10){
    res <- rbind(res,subset(enrich.res, NES > 0))
  }
  if(dim(subset(enrich.res, NES < 0))[1] >= 10){
    total <- dim(enrich.res)[1]
    res <- rbind(res,enrich.res[(total-9):total,])
  }
  if(dim(subset(enrich.res, NES < 0))[1] < 10){
    res <- rbind(res,subset(enrich.res, NES < 0))
  }
  res$Description <- factor(res$Description,levels = levels(reorder(res$Description,res$NES)))
  res.new <- res[match(levels(res$Description),res$Description),]
  return(res.new)
}

#new gesa plot
bubble_plot <- function(data, treatment, control, limit) {

  final_up <- data[data$NES > 0,]
  final_down <- data[data$NES < 0,]

  final_up <- final_up[order(final_up$pvalue, decreasing = FALSE),]
  final_down <- final_down[order(final_down$pvalue, decreasing = FALSE),]

  final <- rbind(head(final_up, 15), head(final_down, 15))
  final <- final[order(final$NES, decreasing = TRUE),]
  final$qvalues <- ifelse(final$qvalues <= 0.15, final$qvalues, limit)
  final$Description <- factor(final$Description, levels = rev(unique(final$Description)))

  final$condition <- ifelse(final$NES > 0, paste0("Enriched in ", treatment), paste0("Enriched in ", control))
  final$condition <- factor(final$condition, levels = c(paste0("Enriched in ", control), paste0("Enriched in ", treatment)))
  p <- ggplot(final, aes(x= NES, y= factor(Description), label = qvalues, color=qvalues,size=setSize)) + geom_point() +
      facet_grid(~condition, scales = "free_x") + theme_bw()

  p <- p + scale_color_gradient(low = "red", high = "blue", limits = c(0,limit)) + 
		                theme(axis.title.y = element_blank(),
                                      axis.text = element_text(size = 18),
                                      strip.text = element_text(size = 16),
                                      axis.title.x = element_text(size = 18))
				
  return(p)
}


writeoutput <- function(result,onto){
write.table(format(result), paste(Outdir,Condition, "_",Treatment,"_vs_",Control,"_",onto,"_terms.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(result@result, paste(Outdir,Condition, "_",Treatment,"_vs_",Control,"_",onto,"_allterms.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
p <- bubble_plot(result@result, Treatment, Control, limit = 0.15)
png(paste(Outdir,Condition, "_",Treatment,"_vs_",Control,"_",onto,"_terms.png", sep = ""),width = 3000, height = 1000)
print(p)
dev.off()
}


###For GO
go.mf <- GSEAGO(geneList,"MF",minGSSize,nPerm,pcut)
go.bp <- GSEAGO(geneList,"BP",minGSSize,nPerm,pcut)
go.cc <- GSEAGO(geneList,"CC",minGSSize,nPerm,pcut)

###For KEGG
kegg <- GSEAKEGG(geneList,minGSSize,nPerm,pcut)

###For HALLMARK
hallmark_gene <- read.gmt(h_file)
hallmark <- GSEAHALLMARK(geneList,minGSSize,pcut,hallmark_gene)

###write enriched term table
writeoutput(go.mf,"MF")
writeoutput(go.bp,"BP")
writeoutput(go.cc,"CC")
writeoutput(kegg,"KEGG")
writeoutput(hallmark,"HALLMARK")





#############multiqc file
#go.mf.multiqc <- subset(go.mf,  select = c(Description,NES,p.adjust))
#colnames(go.mf.multiqc) <- c("Description","MF","Pvalue")
#go.mf.multiqc  <- go.mf.multiqc[order(-go.mf.multiqc$Pvalue),]
#go.mf.multiqc <- subset(go.mf.multiqc[1:10,],select=c(Description,MF))
#write.table(go.mf.multiqc,paste(opt$outdir,"_GO_MF_report.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

#go.bp.multiqc <- subset(go.bp,  select = c(Description,NES,p.adjust))
#colnames(go.bp.multiqc) <- c("Description", "BP","Pvalue")
#go.bp.multiqc  <- go.bp.multiqc[order(-go.bp.multiqc$Pvalue),]
#go.bp.multiqc <- subset(go.bp.multiqc[1:10,],select=c(Description,BP))
#write.table(go.bp.multiqc,paste(opt$outdir,"_GO_BP_report.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

#go.cc.multiqc <- subset(go.cc,  select = c(Description,NES,p.adjust))
#colnames(go.cc.multiqc) <- c("Description", "CC", "Pvalue")
#go.cc.multiqc <- go.cc.multiqc[order(-go.cc.multiqc$Pvalue),]
#go.cc.multiqc <- subset(go.cc.multiqc[1:10,],select=c(Description,CC))
#write.table(go.cc.multiqc,paste(opt$outdir,"_GO_CC_report.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

#go.kegg.multiqc <- subset( kegg,  select = c(Description,NES,p.adjust))
#colnames(go.kegg.multiqc) <- c("Description", "KEGG","Pvalue")
#go.kegg.multiqc <- go.kegg.multiqc[order(-go.kegg.multiqc$Pvalue),]
#go.kegg.multiqc <- subset(go.kegg.multiqc[1:10,],select=c(Description,KEGG))
#write.table(go.kegg.multiqc,paste(opt$outdir,"_GO_KEGG_report.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
