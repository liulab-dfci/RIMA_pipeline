suppressMessages(library(sva))
suppressMessages(library(limma))
suppressMessages(library(optparse))



# make option list and parse command line
option_list <- list(
  make_option(c("-e", "--expression_dat"), type="character",
              help="Input path of expression file. [Required]"),
  make_option(c("-c", "--covariates"), type="character",
              help="covariates needs to be adjusted for"),
  make_option(c("-m", "--metasheet"), type="character",
              help="covariates needs to be adjusted for"),   
  make_option(c("-o","--output",type="character", help="Output files [Required]"))
)


opt_parser <- OptionParser(option_list=option_list);
opts <- parse_args(opt_parser);

# paramenter checking
if(is.null(opts$expression_dat)) stop('Expression file required.')  ###if not provide batch file output log expression matrix

###functions for inputing and outputing
ssgsvaFormat <- function(dat){
  dat <- cbind(Gene_ID=rownames(dat),dat)
  return(dat)
}

writeDF <- function(dat,path){
  write.table(dat,path,quote = FALSE, sep=',', row.names = FALSE)
}

# load data
  expr.dat <- read.table(opts$expression_dat,sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1,check.names = FALSE)
  expr.dat <- log2(expr.dat + 1)

print('Load data done!')


###filtering out genes with low variance among samples

CVFILTER <- 0
mean_nolym <- apply(expr.dat,1,mean)
var_nolym <- apply(expr.dat,1,var)
cv_nolym <- abs(var_nolym/mean_nolym)
filt_genes <- subset(cv_nolym, cv_nolym > CVFILTER)

## Select those genes that pass variance filtering
exprZero <- expr.dat
expr.dat <- expr.dat[rownames(expr.dat) %in% names(filt_genes),]
exprZero <- subset(exprZero, !(rownames(exprZero) %in% names(filt_genes)))

batch.dat = read.table(opts$metasheet, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1,check.names = FALSE)
overlap_sample = intersect(colnames(expr.dat),rownames(batch.dat))
print('Load meta done!')

# stop at low sample numbers
#if(length(exprZero) <- dim(expr.dat)[2]/2) stop("too few samples")

print('Processing samples:')
print(overlap_sample)


if(opts$covariates == "False"){
  expr.limma <- expr.dat
    } else {
      expr.limma = tryCatch(
      removeBatchEffect(as.matrix(expr.dat[,overlap_sample]),
                          batch.dat$opt$covariates),
      error = function(e){
      print(e)
      })
    }

expr.limma = rbind(expr.limma,exprZero)
writeDF(ssgsvaFormat(expr.limma),opts$output)


