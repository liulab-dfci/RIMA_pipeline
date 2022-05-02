suppressMessages(library(optparse))

## make option list and parse command line
option_list = list(
  make_option(c("-f", "--input"), type="character", default=NULL,
              help="list input files", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL,
              help="metasheet info", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./",
              help="output file path and prefix ", metavar="character"),
  make_option(c("-g", "--condition"), type="character", 
              help="Condition to do comparison", metavar="character")
            )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

##set options
meta <- opt$meta
Condition <- opt$condition
outdir <- opt$outdir
files <- unlist(strsplit(opt$input,","))

meta <- read.table(file = meta, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
samples <- subset(meta, meta[,Condition] != 'NA')
print(paste("There are ", length(rownames(samples)), " samples to  merge ...",sep=''))

filelist.samples<- sapply(rownames(samples), function(x) grep(paste0("\\b",x,"\\b"), files, value = TRUE))
print (filelist.samples)

merge_process <- function(filelist){
  output <- list()
  
  for( i in names(filelist))
   {
    print (filelist[i])
    dat <- read.table(filelist[[i]],header= T,sep ='\t')  
    output[[i]] <- dat
    }
    
    mat <- do.call(rbind.data.frame,output)
    return(mat)
}


msi_score <- merge_process(filelist.samples)
colnames (msi_score) <- c('Total_Number_of_Sites','Number_of_Somatic_Sites','MSI_score')
write.table(msi_score, paste(outdir,Condition,"_msi_score.txt",sep=''),sep = '\t',quote=FALSE,row.names = TRUE)
