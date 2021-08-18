rm(list=ls())
suppressMessages(library(DESeq2))
suppressMessages(library(optparse))
suppressMessages(library(tximport))
suppressMessages(library(dplyr))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="list input files", metavar="character"),
  make_option(c("-b", "--batch"), type="character", default=NULL, 
              help="control batch effect or not", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="data source", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="metasheet info", metavar="character"),
  make_option(c("-c", "--comparison"), type="character", default=NULL, 
              help="comparison type", metavar="character"),
  make_option(c("-x", "--tx2gene"), type="character", default=NULL, 
              help="transcript annotation", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./",
              help="output file path and prefix ", metavar="character"),
  make_option(c("-p", "--condition"), type="character", default="./",
              help="Condition to do comparison", metavar="character"),
  make_option(c("-q", "--multiqc"), type="character", default="./",
              help="Path of multiqc folder", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


if (is.null(opt$input) || is.null(opt$tx2gene) || is.null(opt$type)){
  print_help(opt_parser)
  stop("At least three arguments must be supplied ", call.=FALSE)
}


######################################--------------programme-------------##################################
###----If you have transcript quantification files, as produced by Salmon, Sailfish, or kallisto, you would use DESeqDataSetFromTximport.

Transcript <- function(files,samples,tx2gene,Type,batch){
  filelist <- strsplit(files, "\\,")[[1]]
  filelist.samples <- sapply(rownames(samples), function(x) grep(x, filelist, value = TRUE))
  filelist.samples <- filelist.samples[lapply(filelist.samples,length)>0]
  print(filelist.samples)
  print(length(filelist.samples))
  txi <- tximport(filelist.samples, type=Type, tx2gene=tx2gene)
  txi$length[txi$length == 0] <- 1
  
  if(batch != "False"){
    ddsTxi <- DESeqDataSetFromTximport(txi,
                                       colData = samples,
                                       design = ~ Batch + Condition)
  }else{
    ddsTxi <- DESeqDataSetFromTximport(txi,
                                       colData = samples,
                                       design = ~ Condition)
  }
  dds <- DESeq(ddsTxi)
  return (dds)
}


IDConversion <- function(Type, res, tx2gene){
	res.convertID <- cbind.data.frame(gene_name =rownames(res),res)
  return (res.convertID)
}

tx2gene <- read.table(file = opt$tx2gene,sep = ",",header = TRUE)

comp.type <- opt$comparison

comp <- opt$meta
condition <- opt$condition

samples <- read.table(file = comp, sep='\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
if(ncol(samples) == 0) {
  samples <- read.table(file = comp, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
}

if(opt$batch != "False") {
  batch <- opt$batch
#  samples.all <- subset(samples, Batch != "" & Condition != "")
  samples.all <- samples[c(batch, condition)]
  colnames(samples.all) <- c("Batch", "Condition")
  samples.all$Batch <- factor(samples.all$Batch)
} else {
  samples.all <- samples[condition]
  colnames(samples.all) <- c("Condition")
}
#samples.all <- subset(samples, Batch != "" & Condition != "")
print(paste("There are ", dim(samples.all)[1], " samples to be compared", sep = ""))



if(comp.type == "between"){
  samples.all$Condition <- factor(samples.all$Condition)
  ss <- t(as.matrix(combn(levels(samples.all$Condition),2)))
  dds <- Transcript(files = opt$input, samples = samples.all, tx2gene = tx2gene, Type = opt$type, batch = opt$batch)  
  print(class(dds))

  for(i in 1:dim(ss)[1]){
    print (paste("Compare ", ss[i,1], " and ", ss[i,2], sep = ""))
    res <- results(dds, contrast = c("Condition",c(ss[i,1],ss[i,2])))
    res.convertID <- IDConversion(opt$type, res, tx2gene)
    res_final <- as.data.frame(res)
    res_final$Gene_name <- rownames(res_final)
    res_final <- res_final[c(7, 1:6)]
    res_final$`-log10(padj)` <- -log10(res_final$padj)
    sub_res <- subset(res_final, abs(log2FoldChange) > 1)
    sub_res <- sub_res[order(sub_res$`-log10(padj)`, decreasing = TRUE),]
    sub_res <- head(sub_res, 100)
    write.table(sub_res, file = paste(opt$multiqc,"_",ss[i,1],"_VS_",ss[i,2],'_DESeq2_sub.txt',sep = ""), quote = FALSE,sep = "\t", row.names = FALSE)
    write.table(res_final,file = paste(opt$outpath,"_",ss[i,1],"_VS_",ss[i,2],'_DESeq2_raw.txt',sep = ""), quote = FALSE,sep = "\t", row.names = FALSE)
    write.table(res.convertID,file = paste(opt$outpath,"_",ss[i,1],"_VS_",ss[i,2],'_DESeq2_ConvertID.txt',sep = ""), quote = FALSE,sep = "\t")
  } 
}

if(comp.type == "loop"){
  for(gr in unique(as.character(samples.all$Condition))){
    samples.loop <- samples.all %>% mutate(Condition_loop = factor(ifelse(Condition == gr, gr, "others"))) %>% mutate(Condition = Condition_loop)
    print(samples.loop$Condition)
    dds <- Transcript(files = opt$input, samples = samples.loop, tx2gene = tx2gene, Type = opt$type, batch = opt$batch)  
    print (paste("Compare ", gr, " and others", sep = ""))
    res <- results(dds, contrast = c("Condition",c(gr,"others")))
    res.convertID <- IDConversion(opt$type, res, tx2gene)
    res_final <- as.data.frame(res)
    res_final$Gene_name <- rownames(res_final)
    res_final <- res_final[c(7, 1:6)]
    res_final$`-log10(padj)` <- -log10(res_final$padj)
    sub_res <- subset(res_final, abs(log2FoldChange) > 1)
    sub_res <- sub_res[order(sub_res$`-log10(padj)`, decreasing = TRUE),]
    sub_res <- head(sub_res, 100)
    write.table(sub_res,file = paste(opt$multiqc,"_",gr,"_VS_others_DESeq2_sub.txt",sep = ""), quote = FALSE,sep = "\t", row.names = FALSE)
    write.table(res_final,file = paste(opt$outpath,"_",gr,"_VS_others_DESeq2_raw.txt",sep = ""), quote = FALSE,sep = "\t", row.names = FALSE)
    write.table(res.convertID,file = paste(opt$outpath,"_",gr,"_VS_others_DESeq2_ConvertID.txt",sep = ""), quote = FALSE,sep = "\t")
  }
}




