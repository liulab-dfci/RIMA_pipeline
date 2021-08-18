library(optparse)

option_list <- list( 
  make_option(c("-d", "--deg"), type="character",
              help="input cdr3 file path"),
  make_option(c("-l", "--lisa_gene"), type="character",
              help="input stat file path"),
  make_option(c("-o", "--outdir"), type="character",
              help="column number of clinic phenotype traits in meta file[Required]")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

##test
#deseq_result <- read.table("~/Desktop/comp1_DESeq2_ConvertID.txt", stringsAsFactors = FALSE)
#lisa_gene <- read.table("~/Desktop/lisa_gene.txt", stringsAsFactors = FALSE)
deseq_result <- read.table(opt$deg, stringsAsFactors = FALSE)
lisa_gene <- read.table(opt$lisa_gene, stringsAsFactors = FALSE)
colnames(lisa_gene) <- "gene_name"
merged_gene <- merge(deseq_result,lisa_gene)
merged_gene <- merged_gene[!is.na(merged_gene[,3]),]
uniq_gene <- merged_gene[!duplicated(merged_gene$gene_name),]
down <- as.data.frame(uniq_gene[order(uniq_gene[,3]),][c(1:200),1])
up <- as.data.frame(uniq_gene[order(uniq_gene[,3],decreasing = TRUE),][c(1:200),1])
design <- unlist(strsplit(opt$deg,"/"))[3]
write.table(up,paste0(opt$outdir, "/", design, "/", design, "_upRegGenes/", design, "_upRegGenes.txt"),col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(down,paste0(opt$outdir, "/", design, "/", design, "_downRegGenes/", design, "_downRegGenes.txt"),col.names = FALSE, row.names = FALSE, quote = FALSE)
