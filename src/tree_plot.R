library(tidyverse)
library(ggtree)
library(optparse)
###make option list and parse command line
option_list <- list(
  make_option(c("-i","--input"), type = "character",
              help = "path of nwk file"),
  make_option(c("-o","--outdir"), type = "character",
              help = "output directory")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

###read file
tree <- read.tree(opt$input)
edge <- as.data.frame(tree$edge)
###estimate max layer number
count_list <- c()
for (i in (length(tree$tip.label)+1):length(tree$edge.length)){
  count=0
  row <- which(edge[,2]==i)
  while(length(row)!=0){
    j <- edge[row,1]
    row <- which(edge[,2]==j)
    count <- count+1
    print(row)
  }
  count_list[(i-(length(tree$tip.label)))] <- count
}
maxLength <- max(count_list)
###plot
png(paste(opt$outdir, "/", "selected_microbes_phylogenetic.png", sep = ""), res = 350, width = 500*(maxLength+4), height = 70*length(tree$tip.label))

ggtree(tree, branch.length = "tree") + geom_nodepoint() +  #geom_text(aes(label=node), hjust=-.3, size = 3)+
  ggtitle("Phylogenetic Tree") +geom_tiplab(linesize=.5) +
  theme_tree2(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"))+
  xlim_tree(maxLength+4)
dev.off()