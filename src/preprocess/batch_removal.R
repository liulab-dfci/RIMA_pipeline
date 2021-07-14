suppressMessages(library(sva))
suppressMessages(library(limma))
suppressMessages(library(optparse))



# make option list and parse command line
option_list <- list(
  make_option(c("-e", "--expression_dat"), type="character",
              help="Input path of expression file. [Required]"),
  make_option(c("-c", "--covariates"), type="character",
              help="covariates needs to be adjusted for"),
  make_option(c("-d", "--design"), type="character",
              help="sample to include in the design column"),
  make_option(c("-m", "--metasheet"), type="character",
              help="metasheet"),   
  make_option(c("-b","--output_before",type="character", 
              help="Output files [Required]")),
  make_option(c("-a","--output_after",type="character",
              help="Output files [Required]"))
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

# Get samples
Condition <- opts$design
print(Condition)
print(opts$covariates)
meta <- read.table(file = opts$metasheet, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
samples <- subset(meta, meta[,Condition] != 'NA')
print(samples)

# load data
expr.dat <- read.table(opts$expression_dat,sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1,check.names = FALSE)
expr.dat <- log2(expr.dat + 1)
expr.dat <- expr.dat[,rownames(samples)]
writeDF(ssgsvaFormat(expr.dat),opts$output_before)
print('Load data done!')
print(head(expr.dat))


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


if(opts$covariates == "False"){
        print('No batches used !')
	expr.limma <- expr.dat
    } else {
	print('Running limma for batch removal')
      expr.limma = tryCatch(
      removeBatchEffect(as.matrix(expr.dat),
                          samples$opts$covariates),
      error = function(e){
      print(e)
      })
    }
print(head(expr.limma))
expr.limma = rbind(expr.limma,exprZero)
writeDF(ssgsvaFormat(expr.limma),opts$output_after)


