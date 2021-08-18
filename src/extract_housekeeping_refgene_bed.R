##read in house keeping genes
load("~/Documents/rnaseq/data/Zuber_Essential.rdata")
write.table(Zuber_Essential, file = "~/Documents/rnaseq/data/Zuber_house_keeping_gene_list.txt",quote = FALSE, sep = "\t", row.names = FALSE)

##read in refgene bed file
ref.bed <- read.table(file = "~/Documents/rnaseq/data/refseqGenes.bed",sep = "\t")
sub.ref.bed <- subset(ref.bed, V4 %in% Zuber_Essential$GeneSymbol)
write.table(sub.ref.bed, file = "~/Documents/rnaseq/data/housekeeping_refseqGenes.bed",sep = "\t", quote = FALSE, row.names = FALSE)
