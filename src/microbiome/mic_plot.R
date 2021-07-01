library(vegan)
library(ggplot2)
library(dplyr)
library(data.table)
library(optparse)

# make option list and parse command line
option_list <- list(
  make_option(c("-i","--input"), type="character",
              help="Path of merged microbiota abundance data. [Required]"),
  make_option(c("-o","--outdir"),type="character",
              help="Output files [Required]"),
  make_option(c("-c","--clinic_col"),type = "character",
              help = "Phenotypes to be compared"),
  make_option(c("-m","--meta"),type = "character",
              help = "Path of meta data")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

####set parameters
outdir<-opt$outdir
col <- opt$clinic_col
input <- opt$input
meta <- opt$meta

#multiqc heatmap re-order function 
order_ta <- function(x, col) {
  #col_order
  meta <- meta[order(meta[[col]], decreasing = TRUE),]
  
  orde_col <- data.frame("samples" = meta$SampleName, num = 1:nrow(meta), pheo = meta[[col]])
  re_col <- data.frame("samples" = colnames(x), re_num = 1:nrow(meta))
  
  orde_col <- merge(orde_col, re_col, by = 1)
  orde_col <- orde_col[order(orde_col[["num"]], decreasing = FALSE),]
  
  #cluster by the column
  t <- rownames(x)
  for (sa in unique(orde_col$pheo)) {
    tmp_sa <- subset(orde_col, pheo == sa)
    tmp_1 <- x[tmp_sa[["re_num"]]]
    nw_order <- data.frame("ncol" = 1:ncol(tmp_1), "colsum" = colSums(tmp_1))
    nw_order <- nw_order[order(nw_order$colsum, decreasing = TRUE),]
    tmp_1 <- tmp_1[nw_order[["ncol"]]]
    t <- cbind(t, tmp_1)
  }
  t <- t[-1]
  
  #row_order based on expression
  t$ave <- rowSums(t)
  t <- t[order(t$ave, decreasing = TRUE),]
  t <- t[-ncol(t)]
  return(t)
}

####read in merged microbiota abundance data
micro <- fread(input, sep = "\t", header = TRUE)

####read in meta data
meta <- read.table(meta, sep = ",", header = TRUE)
meta <- subset(meta, meta[,col] != 'NA')


#normalization
tem <- subset(micro, taxID != 9606 & taxID != 32630)
tem$relative_ab <- (tem$abundance - min(tem$abundance))/(max(tem$abundance) - min(tem$abundance))
tem$ratio_nor <- tem$relative_ab/sum(tem$relative_ab)
tem <- subset(tem, ratio_nor >= 0.0001)


#convert the result table into lefse
t <- unique(tem$name)

new_matrix <- matrix( nrow = length(unique(tem$taxID)), ncol = nrow(meta), dimnames = list(t, meta$SampleName))
for (sam in colnames(new_matrix)) {
  for (mi in rownames(new_matrix)) {
    n <- subset(tem, sample == sam & name == mi)
    n <- n$relative_ab
    if (length(n) == 0) {
      n = 0
    }
    new_matrix[mi,sam] <- n
  }
}

#diversity
div <- diversity(new_matrix, index = "invsimpson")
new_matrix <- cbind(new_matrix, div)
new_matrix <- new_matrix[order(new_matrix[,"div"], decreasing = TRUE),]

###output selected microbes with abundance
micros.out <- micro %>% dplyr::filter(name %in% rownames(new_matrix)) %>% dplyr::select(name, taxID) %>% distinct()
write.table(micros.out, file = paste(outdir, "selected_microbes_taxonID.txt", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#calculate top 15 microbiome abundance
mi_up <- head(new_matrix,15)
mi_down <- new_matrix[-1:-15,]
other <- colSums(mi_down[,1:nrow(meta)])
mi_up <- mi_up[,1:(ncol(mi_up)-1)]
final_matrix <- rbind(mi_up, other)

mi_up <- as.data.frame(mi_up)
multiqc_heatmap <- order_ta(mi_up, col)
write.table(multiqc_heatmap, paste(outdir, "Microbiome_abundance_ratio.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")

#convert to the table
final_ta <- NULL
for (t in 1:ncol(final_matrix)) {
  #tmp <- final_matrix[,t]
  tmp_ta <- data.frame("sample" = colnames(final_matrix)[t], "abundance" = final_matrix[,t], "species" = names(final_matrix[,t]))
  final_ta <- rbind(final_ta, tmp_ta)
}

#reorder sample 
ord_sa <- NULL
for (sa in unique(final_ta$sample)) {
  tmp <- subset(final_ta, sample == sa)
  sa_tem <- data.frame("sample" = sa, "sum_ab" = colSums(tmp[2]))
  ord_sa <- rbind(ord_sa, sa_tem)
}

final_ta <- merge(final_ta, ord_sa, by = "sample")
final_ta <- final_ta[order(final_ta$sum_ab, decreasing = FALSE),]
final_ta$sample <- factor(final_ta$sample, levels = unique(final_ta$sample))
#ggplot(final_ta, aes(x= abundance, y = sample, fill = species)) + geom_bar(stat = "identity")

#calculate the ratio
final_rat <- NULL
for (rat in unique(final_ta$sample)) {
  tmp <- subset(final_ta, sample == rat)
  tmp$ratio <- tmp$abundance/sum(tmp$abundance)
  if (is.na(sum(tmp$ratio))) {
    tmp$ratio <- 0
  }
  final_rat <- rbind(final_rat, tmp)
}

png(paste(outdir, "microbes_abundance.png", sep = ""), res = 300, width = 1600, height = 1300, pointsize = 12)
ggplot(final_rat, aes(x= ratio, y = sample, fill = species)) + geom_col() +
  theme(
    axis.text.x=element_text(size=6, face = "bold"),
    axis.text.y=element_text(size=6, face = "bold"),
    axis.title.y = element_text(size = 6, face = "bold"),
    axis.title.x = element_text(size = 6, face = "bold"),
    legend.title = element_text(size=6, face = "bold"),
    legend.text = element_text(size=6, face = "bold")
  )
  
dev.off()

#########################
#multiqc table: re-order the sample
re_info <- meta[c("SampleName", col)]

colnames(re_info) <- c("sample", "pheo")
multiqc <- merge(final_rat, re_info, by = "sample")


##order by column 
t <- NULL
for (sa in unique(multiqc$pheo)) {
  tmp_sa <- subset(multiqc, pheo == sa)
  tmp_sa <- tmp_sa[order(tmp_sa$sum_ab, decreasing = TRUE),]
  t <- rbind(t, tmp_sa)
}

multiqc <- t

multiqc <- multiqc[c("species","ratio","sample")]
multiqc$species <- gsub(" ","_", multiqc$species)

#convert into multiqc table matrix
t <- unique(multiqc$species)
multiqc_matrix <- matrix(nrow = nrow(meta), ncol = length(t), dimnames = list(unique(multiqc$sample),t))
for (sam in colnames(multiqc_matrix)) {
  for (mi in rownames(multiqc_matrix)) {
    n <- subset(multiqc, species == sam & sample == mi)
    n <- n$ratio * 100
    if (length(n) == 0) {
      n = 0
    }
    multiqc_matrix[mi,sam] <- n
  }
}
samples <- rownames(multiqc_matrix)
multiqc_matrix <- cbind(samples, multiqc_matrix)

write.table(multiqc_matrix, paste(outdir, "Microbiome_abundance.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")



