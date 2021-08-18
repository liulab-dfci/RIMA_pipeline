library(mMCPcounter)
library(optparse)

# make option list and parse command line
library(optparse)
option_list <- list(  
  make_option(c("-e", "--expression_dat"), type="character", 
              help="Input path of expression file. [Required]"),
  make_option(c("-o", "--output_dat"), type="character",
              help="Ouput path of corresponding tables. [Required]")
)
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);

# paramenter checking
if(is.null(opts$expression_dat)) stop('Expression file required.')

exp <- read.table(opts$expression_dat, sep = "\t", row.names = 1, header = TRUE)
results <- mMCPcounter.estimate(exp, features = c("Gene.Symbol","ENSEMBL.ID","Probes")[1])

write.table(results,paste0(opts$output_dat, '/mMCP.txt', sep=""),sep='\t', row.names=TRUE, quote = FALSE, col.names=NA)
