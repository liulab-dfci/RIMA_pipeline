library(vegan)
library(pheatmap)

clustering_heatmap <- function(dds, res, pgenes, outpath, treatment, control) {
  
  print("loding the human protein coding genes...")
  pgenes <- read.csv(pgenes)	
  deseq_result <- res
  deseq_result <- as.data.frame(deseq_result)
  
  deseq_result <- deseq_result[rownames(deseq_result) %in% pgenes[[1]],]
  deseq_result <- na.omit(deseq_result)
  
  #deseq_result <-deseq_result[order(as.numeric(deseq_result[,"stat"])),]
  
  #get the normalized count matrix 
  cdata <- as.data.frame(counts(dds, normalized = TRUE))
  cdata <- cdata[rownames(cdata) %in% rownames(deseq_result),]
  
  #cluster the sample based on phenotype 
  coldata <- as.data.frame(colData(dds))
  coldata <- coldata[order(coldata$Condition, decreasing = TRUE),]
  
  cdata <- cdata[rownames(coldata)]
  
  print("Calculating the expression density...")
  div <- diversity(cdata, index = "invsimpson")
  
  cdata <- cbind(cdata, div)
  
  deseq_result <- merge(deseq_result, cdata[c("div")], by = 0, all = FALSE)
  
  #calculate the median of diversity for up-regulated genes and down-regulated genes
  up_div <- median(deseq_result[deseq_result$log2FoldChange > 0 & deseq_result$pvalue <= 0.05,]$div)
  down_div <- median(deseq_result[deseq_result$log2FoldChange < 0 & deseq_result$pvalue <= 0.05,]$div)
  
  #rank the genes
  deseq_result <-deseq_result[order(as.numeric(deseq_result[,"stat"])),]
  
  down_50 <- deseq_result[deseq_result$div > down_div,]$Row.names[1:50]
  up_50 <- tail(deseq_result[deseq_result$div > up_div,], 50)$Row.names
  
  # #cluster the sample based on phenotype 
  # coldata <- as.data.frame(colData(dds))
  # coldata <- coldata[order(coldata$Condition, decreasing = TRUE),]
  
  
  df <- cdata[c(down_50, up_50),]
  #df <- df[rownames(coldata)]
  
  df <- df[-ncol(df)]
  #calculate the z-score across the samples 
  a <- t(scale(t(df)))
  
  pdf(paste0(outpath,'heatmap_',treatment,'_vs_',control,'.pdf'), width = 5, height = 20)
  pheatmap(a, cluster_cols = FALSE, cluster_rows = FALSE, annotation = coldata["Condition"])
  dev.off()
}
#deseq_result <- res
#deseq_result <- as.data.frame(deseq_result)
#
#deseq_result <- deseq_result[rownames(deseq_result) %in% pgenes$symbol,]
#deseq_result <- na.omit(deseq_result)
#
##deseq_result <-deseq_result[order(as.numeric(deseq_result[,"stat"])),]
#
##get the normalized count matrix 
#cdata <- as.data.frame(counts(dds, normalized = TRUE))
#cdata <- cdata[rownames(cdata) %in% rownames(deseq_result),]
#
##cluster the sample based on phenotype 
#coldata <- as.data.frame(colData(dds))
#coldata <- coldata[order(coldata$Condition, decreasing = TRUE),]
#
#cdata <- cdata[rownames(coldata)]
#
#div <- diversity(cdata, index = "invsimpson")
#
#cdata <- cbind(cdata, div)
#
#deseq_result <- merge(deseq_result, cdata[c("div")], by = 0, all = FALSE)
#
##calculate the median of diversity for up-regulated genes and down-regulated genes
#up_div <- median(deseq_result[deseq_result$log2FoldChange > 0 & deseq_result$pvalue <= 0.05,]$div)
#down_div <- median(deseq_result[deseq_result$log2FoldChange < 0 & deseq_result$pvalue <= 0.05,]$div)
#
##rank the genes
#deseq_result <-deseq_result[order(as.numeric(deseq_result[,"stat"])),]
#
#down_50 <- deseq_result[deseq_result$div > down_div,]$Row.names[1:50]
#up_50 <- tail(deseq_result[deseq_result$div > up_div,], 50)$Row.names
#
## #cluster the sample based on phenotype 
## coldata <- as.data.frame(colData(dds))
## coldata <- coldata[order(coldata$Condition, decreasing = TRUE),]
#
#
#df <- cdata[c(down_50, up_50),]
##df <- df[rownames(coldata)]
#
#df <- df[-ncol(df)]
##calculate the z-score across the samples 
#a <- t(scale(t(df)))
#
#pdf("/Users/linyang/Documents/Rplot08.pdf", width = 5, height = 20)
#pheatmap(a, cluster_cols = FALSE, cluster_rows = FALSE, annotation = coldata["Condition"])
#dev.off()
#
