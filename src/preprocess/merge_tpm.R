suppressMessages(library(DESeq2))
suppressMessages(library(optparse))
suppressMessages(library(tximport))
suppressMessages(library(dplyr))
suppressMessages(library(sva))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="list input files", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL,
              help="salmon", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL,
              help="metasheet info", metavar="character"),
  make_option(c("-x", "--tx2gene"), type="character", default=NULL,
              help="transcript annotation", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./",
              help="output file path and prefix ", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

metadata <- opt$meta
gene <- opt$tx2gene
Type <- opt$type
input <- opt$input
    
print("Reading meta file ...")
meta <- read.table(file = metadata, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

print ("Reading tx2gene file ...")
tx2gene <- read.table(file = gene,sep = ",",header = TRUE)

filelist <- strsplit(input, "\\,")[[1]]
print(filelist)
print(rownames(meta))
filelist.samples <- sapply(rownames(meta), function(x) grep(paste0("\\b",x,"\\b"), filelist, value = TRUE))
filelist.samples <- filelist.samples[lapply(filelist.samples,length)>0]
print(paste("There are ",length(filelist.samples), " samples to be merge ...", sep = ""))
print(filelist.samples)

txi <- tximport(filelist.samples, type=Type, tx2gene=tx2gene)
  
print(head(txi$abundance))
exprsn <- txi$abundance
write.table(exprsn,paste(opt$outpath,'tpm.genesymbol.csv',sep=""),quote = FALSE,sep = ",")
