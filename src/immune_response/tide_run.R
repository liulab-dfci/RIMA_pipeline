rm(list=ls())

suppressMessages(library(optparse))
suppressMessages(library(org.Hs.eg.db))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="list input files", metavar="character"),
  make_option(c("-d", "--design"), type="character", default=NULL,
              help="design column", metavar="character"),
  make_option(c("-t", "--cancer"), type="character", default="./",
              help="cancer type", metavar="character"),
  make_option(c("-p", "--treated"), type="character", default="./",
              help="pretreated", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./",
              help="output file", metavar="character")            
              )
	      
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


input <- opt$input
cancer_type <- opt$cancer
outdir <- opt$outdir
treat <- opt$treated
design <- opt$design

exprsn <- read.table(file = input, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
tide_input <- exprsn - rowMeans(exprsn)
GeneIDSymbol <- toTable(org.Hs.egSYMBOL)
entrizid <- GeneIDSymbol[match(rownames(tide_input),GeneIDSymbol$symbol),'gene_id']
write.table(cbind(entrizid, tide_input),paste(outdir,design,"_tide_input.txt",sep=''),sep='\t',quote=FALSE,row.names = F)

cancers <-  c('NSCLC','Melanoma')
if (is.element(cancer_type,cancers) ) {
cancer <- cancer_type} else { 
cancer <- 'Other'}

infile <- paste(outdir,design,"_tide_input.txt",sep='')
outfile <- paste(outdir,design,"_tide_output.txt",sep='')


if (treat == 'True' ) {
command = paste ('tidepy -c ',cancer,' --pretreat -o ',outfile,' ', infile , sep= '')
print (command)
} else {
command = paste ('tidepy -c ',cancer,' -o ',outfile,' ', infile , sep= '')
print (command)
}






