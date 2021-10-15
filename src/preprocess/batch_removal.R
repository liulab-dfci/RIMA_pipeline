suppressMessages(library(sva))
suppressMessages(library(limma))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggfortify))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggpubr))


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

pca_plot <- function(exprTable, annot,title,Batch) {
  batch_n <- length(unique(as.character(annot[colnames(exprTable),"Batch"])))
  print(paste("there are ", batch_n, " batches in your data"))
  df <- cbind.data.frame(t(exprTable),batch = as.character(annot[colnames(exprTable),"Batch"]))
  pca_plot <- autoplot(prcomp(t(exprTable)), data = df, col = Batch, size = 1, frame = TRUE, frame.type = 'norm')+
    labs(title=title)+
    scale_color_manual(values = brewer.pal(name = "Set1", n = 9)[1:batch_n])+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=12,face = "bold",hjust=1),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"))
  return(pca_plot)
}


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
      expr.limma = rbind(expr.limma,exprZero)
      writeDF(ssgsvaFormat(expr.limma),opts$output_after)
    
  } else {
      samples$Batch <- samples[,opts$covariates] 
      print(samples)
      
      print('Running limma for batch removal')
      expr.limma = tryCatch(
      removeBatchEffect(as.matrix(expr.dat),
                          samples$Batch),
      error = function(e){
      print(e)
      })
      expr.limma = rbind(expr.limma,exprZero)
      writeDF(ssgsvaFormat(expr.limma),opts$output_after)
      
      #png(file = paste(opts$out_before,'.png',sep=''), res = 300, height = 1200, width = 1500)
      #pca_plot(expr.dat, annot = samples, title = "PCA plot Before Batch Removal",Batch=Batch)
      #dev.off()
     
      #png(file = paste(opt$out_after,'.png',sep=''), res = 300, height = 1200, width = 1500)
      #pca_plot(expr.limma, annot = samples, title = "PCA plot Before Batch Removal",Batch=Batch)
      #dev.off()

    }

