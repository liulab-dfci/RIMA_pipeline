#!/usr/bin/env Rscript

#dependencies
suppressMessages(library(chimeraviz))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(yaml))

option_list = list(
  make_option(c("-f", "--fusion"), type="character", default=NULL, 
              help="merged fusion prediction file", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="metasheet file", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NULL, 
              help="genome type", metavar="character"),
  make_option(c("-l", "--limit"), type="character", default=NULL, 
              help="showing limited fusion number", metavar="numeric"),
  make_option(c("-p", "--phenotype"), type="character", default=NULL, 
              help="phenotypes needs to be compared", metavar="int"),
  make_option(c("-out", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


###read in fusion data
# fusion <- read.table(file = "~/Documents/rnaseq/data/merged_predictions.abridged_addSample.tsv", header = TRUE)
# meta <- read.table("~/Documents/rnaseq/snakemake_files/metasheet_test.csv", sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
fusion <- read.table(file = opt$fusion, header = TRUE)
meta <- read.table(file = opt$meta, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
genome <- opt$genome
limit <- as.numeric(opt$limit)
outdir <- opt$outdir
phenotypes<- opt$phenotype
phenotypes <- strsplit(phenotypes, "\\,")[[1]]

###remove fusion pairs with low junction reads in each sample
spe <- split(fusion, fusion$sample) 
out <- lapply(1:length(spe),function(x){
  sub <- spe[[x]]
  filter.sub <- subset(sub, JunctionReadCount >= median(sub$JunctionReadCount))
})
filter.fusion <- do.call("rbind", out)

###read in meta information
# phenotypes<-yaml.load_file("./config.yaml")$fusion_clinical_phenotypes
print(phenotypes)
for (phenotype in phenotypes){
  gr <- unique(meta[[phenotype]])
  ###show fusion genes in different treatment experiments
  for(t in gr){
    gr.fusion <- subset(filter.fusion, sample %in% rownames(meta)[which(meta[[phenotype]] == t)]) %>% arrange(desc(JunctionReadCount))
    write.table(gr.fusion, file = paste(outdir,"filter_fusion_",phenotype, "_", t,".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
    # write.table(gr.fusion, file = paste("~/Documents/rnaseq/data/","filter_fusion_",t,".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
    ###import fusion
    dat <- import_starfusion(paste(outdir,"filter_fusion_",phenotype, "_", t,".txt", sep = ""), genome, limit)
    
    png(file = paste(outdir,phenotype, "_" ,t,"_fusion_circle_plot.png",sep = ""), res = 350, width = 3000, height = 3000)
    # dat <- import_starfusion(paste("~/Documents/rnaseq/data/","filter_fusion_",t,".txt", sep = ""), "hg38", limit = 150)
    # png(file = paste("~/Documents/rnaseq/data/",t,"_fusion_circle_plot.png",sep = ""), res = 350, width = 3000, height = 3000)
    plot_circle(dat)
    dev.off()
}
}



