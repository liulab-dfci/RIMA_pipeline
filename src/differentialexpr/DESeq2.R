rm(list=ls())
suppressMessages(library(DESeq2))
suppressMessages(library(optparse))
suppressMessages(library(tximport))
suppressMessages(library(dplyr))
suppressMessages(library(sva))
suppressMessages(library(limma))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="list input files", metavar="character"),
  make_option(c("-b", "--batch"), type="character", default=NULL,
              help="control batch effect or not", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL,
              help="metasheet info", metavar="character"),
  make_option(c("-x", "--tx2gene"), type="character", default=NULL,
              help="transcript annotation", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./",
              help="output file path and prefix ", metavar="character"),
  make_option(c("-g", "--condition"), type="character", default="./",
              help="Condition to do comparison", metavar="character"),
  make_option(c("-r", "--treatment"), type="character", default="./",
              help="Treatment", metavar="character"),
  make_option(c("-c", "--control"), type="character", default="./",
              help="Control", metavar="character")
)
	      

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

batch <- opt$batch 
Condition <- strsplit(opt$condition, "\\,")[[1]]
metadata <- opt$meta
gene <- opt$tx2gene
Treatment <- strsplit(opt$treatment, "\\,")[[1]]
print(Treatment)
Control <- strsplit(opt$control, "\\,")[[1]]
print(Control)
input <- opt$input

print("Reading meta file ...")
meta <- read.table(file = metadata, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1)


print ("Reading tx2gene file ...")
tx2gene <- read.table(file = gene,sep = ",",header = TRUE)


if (is.null(opt$input) || is.null(opt$tx2gene)){
  print_help(opt_parser)
  stop("At least 2 arguments must be supplied ", call.=FALSE)
}


######################################--------------programme-------------##################################
###----If you have transcript quantification files, as produced by Salmon, Sailfish, or kallisto, you would use DESeqDataSetFromTximport.

Transcript <- function(files,samples,tx2gene,Type,batch,condition,tx,con){

  filelist <- strsplit(files, "\\,")[[1]]
  print(filelist)
  print(rownames(samples))
  filelist.samples <- sapply(rownames(samples), function(x) grep(x, filelist, value = TRUE))
  
  filelist.samples <- filelist.samples[lapply(filelist.samples,length)>0]
  print(paste("There are ",length(filelist.samples), " samples to be compared", sep = ""))
  print(filelist.samples)
  
  txi <- tximport(filelist.samples, type="salmon", tx2gene=tx2gene)
  txi$length[txi$length == 0] <- 1
  
  print(head(txi$counts))
  exprsn <- txi$counts
  #exprsn_log <- log2(exprsn + 1)

  
  if(batch != "False"){
    print (paste("Running DESeq2 with batch effects ", batch, " on ", condition, sep=""))  
     
    colnames(samples) <- c("Batch","Condition")
    print(samples)
       
    ddsTxi <- DESeqDataSetFromTximport(txi,
                                       colData = samples,
                                       design = ~ Batch + Condition)
                                       
  }else{
    print(paste("Running DESeq2 on ", condition,sep=""))
    
    colnames(samples) <- c('Condition')
    print(samples)
    
    ddsTxi <- DESeqDataSetFromTximport(txi,
                                       colData = samples,
                                       design = ~ Condition)
                                       
    print ("Generating TPM matrix ...")
  }
  dds <- DESeq(ddsTxi)
  write.table(exprsn,paste(opt$outpath, condition,'_', tx,'_vs_', con,'_estimated_genelevel_count.txt',sep = ""),quote = FALSE,sep = "\t") 
  return (dds)
}

#multi-comparison
n = 1
for (c in Condition) {
  tx <- Treatment[n]
  con <- Control[n]
  if(batch != "False") {
    samples <- meta[meta[,c] %in% c(tx,con),][,c(batch,c)]
  }else{
    samples <- meta[meta[,c] %in% c(tx,con),][c]
  }
  
  tmp_id <- as.character(rownames(samples))
  n = n + 1
  if(length(tmp_id) == 0) {
    next
  }
  
  print ("Star running DESeq2 ...")
  dds <- Transcript(files = input, samples, tx2gene = tx2gene, batch = batch, condition = c, tx = tx, con = con)
  print(class(dds))
  
  print (paste("Comparing ", tx, " VS ", con, sep = ""))
  res <- results(dds, contrast = c("Condition", c(tx,con)))
  
  
  res_final <- as.data.frame(res)
  res_final$Gene_name <- rownames(res_final)
  res_final <- res_final[c(7, 1:6)]
  res_final$`-log10(padj)` <- -log10(res_final$padj)
  write.table(res_final,file = paste(opt$outpath, c, '_', tx, '_vs_', con, '_DESeq2.txt', sep = ""), quote = FALSE, sep = "\t", row.names = F)
}

#print ("Star running DESeq2 ...")
#dds <- Transcript(files = input,samples, tx2gene = tx2gene, Type = Type, batch = batch)
#print(class(dds))
#
#print (paste("Comparing ",opt$treatment , " VS ", opt$control, sep = ""))
#res <- results(dds, contrast = c("Condition",c(opt$treatment,opt$control)))
#
#
#res_final <- as.data.frame(res)
#res_final$Gene_name <- rownames(res_final)
#res_final <- res_final[c(7, 1:6)]
#res_final$`-log10(padj)` <- -log10(res_final$padj)
#write.table(res_final,file = paste(opt$outpath,opt$condition,'_',opt$treatment,'_vs_',opt$control,'_DESeq2.txt',sep = ""), quote = FALSE,sep = "\t", row.names = F)
#sub_res <- subset(res_final, abs(log2FoldChange) > 1)
#sub_res <- sub_res[order(sub_res$`-log10(padj)`, decreasing = TRUE),]
#sub_res <- head(sub_res, 100)
#
#write.table(sub_res,file = paste(opt$outpath,opt$condition,'_',opt$treatment,'_vs_',opt$control,'_DESeq2_sub.txt',sep = ""), quote = FALSE,sep = "\t", row.names = F)
