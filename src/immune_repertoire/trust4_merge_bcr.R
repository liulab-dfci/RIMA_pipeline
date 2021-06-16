##load package
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(reshape))
suppressMessages(library(ggpubr))
suppressMessages(library(optparse))

## make option list and parse command line
option_list = list(
  make_option(c("-l", "--bcr_light"), type="character", default=NULL,
              help="list input files", metavar="character"),
  make_option(c("-e", "--bcr_heavy"), type="character", default=NULL,
              help="list input files", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL,
              help="metasheet info", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./",
              help="output file path and prefix ", metavar="character"),
  make_option(c("-g", "--condition"), type="character", 
              help="Condition to do comparison", metavar="character"),
  make_option(c("-t", "--treatment"), type="character", 
              help="Treatment", metavar="character"),
  make_option(c("-c", "--control"), type="character", 
              help="Control", metavar="character")
            )
	      

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

##set options
meta <- opt$meta
Condition <- opt$condition
Treatment <- opt$treatment
Control <- opt$control
outdir <- opt$outdir
files_light <- unlist(strsplit(opt$bcr_light,","))
files_heavy <- unlist(strsplit(opt$bcr_heavy,","))
source("src/immune_repertoire/trust4_metric_functions.R")

##meta
meta <- read.csv(file = meta, sep=",",  header = TRUE, row.names=1)
samples <- subset(meta,meta[,Condition] == opt$treatment | meta[,Condition] == opt$control)
print(paste("There are ", length(rownames(samples)), "samples to  merge ...",sep=''))

filelist.samples.light <- sapply(rownames(samples), function(x) grep(x, files_light, value = TRUE))
filelist.samples.heavy <- sapply(rownames(samples), function(x) grep(x, files_heavy, value = TRUE))

###merge samples
merge_process <- function(filelist){
  output <- list()
  
  for( i in names(filelist)){
    print (filelist[i])
    file <- filelist[[i]]
    dat <- read.table(file,header= T,sep ='\t')  
    
    if(dim(dat)[1] != 0 ) {
      output[[i]]<- dat  %>%
        mutate(clinic = as.character(samples[i,Condition]))
    } else { output[[i]] <- dat}
  
    mat <- do.call(rbind.data.frame, output)
  }
  return(mat)
}

bcr.light <- merge_process(filelist.samples.light)
write.table(bcr.light, paste(outdir,Condition,"_",Treatment,"_vs_",Control,"_TRUST4_BCR_light.txt",sep=''),sep = '\t',quote=FALSE,row.names = FALSE)

bcr.heavy <- merge_process(filelist.samples.heavy)
write.table(bcr.heavy, paste(outdir,Condition,"_",Treatment,"_vs_",Control,"_TRUST4_BCR_heavy.txt",sep=''),sep = '\t',quote=FALSE,row.names = FALSE)


