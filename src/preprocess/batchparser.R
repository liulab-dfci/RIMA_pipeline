####----------Script to subset metasheet into separate "Pre" and "Post" category for running TIDE---------####
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="meta information", metavar="character"),
  make_option(c("-e", "--expression"), type="character",
              help="batch tpm expression data"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory for all the samples", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


meta <- read.table(opt$meta,sep = ",", header = TRUE, row.names = 1,stringsAsFactors = FALSE)
tpm.batch <- as.data.frame(t(read.table(opt$expression,sep="\t",header=T,row.names=1,stringsAsFactors = FALSE)))
outdir <- opt$outdir
#cancer <- "Melanoma"


#meta= read.table("/Users/aashna/Desktop/preprocess_RIMA/immune_response/new_metasheet.csv",sep = ",", header = TRUE, row.names = 1,stringsAsFactors = FALSE)
#tpm.batch =  as.data.frame(t(read.table("/Users/aashna/Desktop/preprocess_RIMA/immune_response/tpm_convertID_batch_Entrez.txt",sep="\t",header=T,row.names=1,stringsAsFactors = FALSE)))
#outdir <- "/Users/aashna/Desktop/preprocess_RIMA/immune_response/"
#cancer<- "LUAD"



if (("Pre" %in% meta[["Timing"]]) && ("Post" %in% meta[["Timing"]])) {
  pre.metasheet <- subset(meta, Timing =="Pre")
  post.metasheet <- subset(meta, Timing == "Post")
  input.tide.pre <- as.data.frame(t(tpm.batch[meta$Timing=="Pre", ]))
  input.tide.post <- as.data.frame(t(tpm.batch[meta$Timing=="Post",]))
  #expression= input.tide.pre
  #cancer <- "Melanoma"
  #pre.meta=data.frame("SampleName"=rownames(pre.metasheet),pre.metasheet)
  write.table(data.frame("SampleName"=rownames(pre.metasheet),pre.metasheet),paste(opt$outdir,"pre/pre.metasheet.csv",sep=""), quote=F,sep=",",row.names=FALSE)
  write.table(data.frame("gene_id"=rownames(input.tide.pre),input.tide.pre),paste(opt$outdir,"pre/tideinput.pre.txt",sep=""),quote=F,sep="\t",row.names=FALSE)
  write.table(data.frame("SampleName"=rownames(post.metasheet),post.metasheet),paste(opt$outdir,"post/post.metasheet.csv",sep=""), quote=F,sep=",",row.names=FALSE)
  write.table(data.frame("gene_id"=rownames(input.tide.post),input.tide.post),paste(opt$outdir,"post/tideinput.post.txt",sep=""),quote=F,sep="\t",row.names=FALSE)
  system(("/home/aashna/miniconda3/envs/rna/bin/tidepy  analysis/batchremoval/pre/tideinput.txt  -o outdir -c Melanoma  --force_normalize --pretreat > analysis/tide/pre_tpm_convertID_batch_tide_score.txt "), intern = TRUE,ignore.stderr = TRUE,ignore.stdout=TRUE, wait= TRUE)
}  else {
  metasheet=meta
  tide= as.data.frame(t(tpm.batch))
  write.table(data.frame("SampleName"=rownames(metasheet),metasheet),paste(opt$outdir,"all.metasheet.csv",sep=""), quote=F,sep=",",row.names=FALSE)
  write.table(data.frame("gene_id"=rownames(tide),tide),paste(opt$outdir,"all.tide.txt",sep=""),quote=F,sep="\t",row.names=FALSE)
}

















