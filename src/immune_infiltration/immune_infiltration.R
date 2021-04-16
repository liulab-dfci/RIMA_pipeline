
#############################Script to process immune infiltration #####################
library(rlang)
library(R.utils)
library(dplyr)
library(immunedeconv)
library(optparse)


# make option list and parse command line
library(optparse)
option_list <- list(
  make_option(c("-e", "--expression_dat"), type="character",
              help="Input path of expression file. [Required]"),
  make_option(c("-o", "--output_dat"), type="character",
              help="Ouput path of corresponding tables. [Required]"),
  make_option(c("-p", "--permutation"), type="numeric", default=NULL,
              help="permutation for cibersort", metavar="character"),
  make_option(c("-q", "--quantile_normalization"), type="character", default=NULL,
              help="whether enable quantile normalization for cibersort", metavar="character"),
  make_option(c("-a", "--absolute"), type="character", default=NULL,
              help="absolute mode for cibersort", metavar="character"),
  make_option(c("-m", "--abs_method"), type="character", default=NULL,
              help="absolute method for cibersort", metavar="character"),
  make_option(c("-t", "--cancer_type"), type="character",
              help="cancer type used in timer.  [Required]")
)
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);


# paramenter checking
if(is.null(opts$expression_dat)) stop('Expression file required.')

input <- opts$expression_dat
output <- opts$output_dat
cancertype <- opts$cancer_type
cibersort_perm <- opts$permutation
cibersort_qn <- opts$quantile_normalization
cibersort_abs <- opts$absolute
cibersort_abs_method <- opts$abs_method

###read in data
df <- read.table(input, sep = ",", row.names = 1,header = TRUE, check.names = FALSE)
print(df)

#quantiseq
res_quant = as.data.frame(deconvolute(df, "quantiseq"))
write.table(res_quant,paste0(output,'/quantiseq.txt'),sep='\t',row.names=FALSE, quote = FALSE)
print("run quantiseq")


#xcell
res_xcell = as.data.frame(deconvolute(df, "xcell"))
write.table(res_xcell,paste0(output,'/xcell.txt'),sep='\t',row.names=FALSE, quote = FALSE)
print("run xcell")


#epic
res_epic = as.data.frame(deconvolute(df, "epic"))
write.table(res_epic,paste0(output,'/epic.txt'),sep='\t',row.names=FALSE, quote = FALSE)
print("run epic")

#mcp_counter
res_mcp = as.data.frame(deconvolute(df, "mcp_counter"))
write.table(res_mcp,paste0(output,'/mcp_counter.txt'),sep='\t',row.names=FALSE, quote = FALSE)
print('run mcp_counter')

#cibersort##removing this as it is replaced by cibersort_abs
#print('s')
#results <- CIBERSORT("LM22.txt",tempfile,100)
#write.table(results, file = paste0(output,'/cibersort.txt'),col.names=NA, quote = FALSE, sep="\t")
#cibersort_abs
source('src/immune_infiltration/CIBERSORT.R')
res_ciber <- CIBERSORT("static/cibersort/LM22.txt",input,perm = cibersort_perm, QN = cibersort_qn, absolute = cibersort_abs, abs_method = cibersort_abs_method)
write.table(res_ciber, file = paste0(output,'/cibersort_abs.txt'), sep="\t", quote = FALSE)
print("run cibersort absolution mode")

res_timer = as.data.frame(deconvolute(df, "timer",indications=rep(tolower(cancertype),ncol(df))))
write.table(res_timer,paste0(output,'/timer.txt'),sep='\t', row.names=FALSE, quote = FALSE)
