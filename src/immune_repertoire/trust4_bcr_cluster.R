suppressMessages(library(dplyr))
suppressMessages(library(reshape))
suppressMessages(library(optparse))
option_list <- list( 
  make_option(c("-i", "--cdr3"), type="character",
              help="input cdr3 file "),
  make_option(c("-c", "--clinic_col"), type="character",
              help="column number of clinic phenotype traits in meta file[Required]"),
  make_option(c("-m", "--meta"), type="character",
              help="meta info[Required]"),
  make_option(c("-o","--output"),type="character",
              help="Output path [Required]")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);


##set options
meta <- opt$meta
clinic.col <- opt$clinic_col
outdir <- opt$output
file <- opt$cdr3

##set work directory
source("src/immune_repertoire/trust4_metric_functions.R")
meta <- read.csv(file = meta,  header = TRUE,  row.names = 1)


###function of processing cdr3 for each separated sample
cdr3_process <- function(file, meta){
  print(file)
  cdr3 <- read.table(file = file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cdr3 <- subset(cdr3, count > 0) %>% 
    mutate(V = as.character(V), J = as.character(J), C = as.character(C), CDR3aa = as.character(CDR3aa)) %>% 
    mutate(clinic = as.character(meta[sample,clinic.col]))
  print(table(cdr3$clinic))
  cdr3$is_complete <- sapply(cdr3$CDR3aa, function(x) ifelse(x != "partial" && x != "out_of_frame" && !grepl("^_",x) && !grepl("^\\?", x),"Y","N"))
  cdr3.bcr <- subset(cdr3, grepl("^IG",V) | grepl("^IG",J) | grepl("^IG",C))
  ##add lib size and clinic traits
  cdr.bcr.new <- cdr3.bcr %>% mutate(lib.size = sum(count)) 
  cdr3.bcr.heavy <- subset(cdr.bcr.new, grepl("^IGH",V) | grepl("^IGH",J) | grepl("^IGH",C))
  ##save bcr 
  return(cdr3.bcr.heavy)
}

##main function of merging processed cdr3 data
cdr3.bcr.heavy <- cdr3_process(file,meta)
ss <- unlist(strsplit(file,"/"))[3]
print(ss)
sample_bcr_cluster <- BuildBCRlineage(sampleID = ss, Bdata = cdr3.bcr.heavy, start=3, end=10)
save(sample_bcr_cluster,file = paste(outdir, "_TRUST4_BCR_heavy_cluster.Rdata", sep = ""))

