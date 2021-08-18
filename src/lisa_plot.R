library(optparse)
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)

###setting options
option_list <- list(
  make_option(c("-i", "--inputdir"), type = "character",
              help = "directory of lisa results"),
  make_option(c("-c", "--clinc"), type = "character",
              help = "clinic name of deseq2 comparison results"),
  make_option(c("-n", "--top_n"), type = "integer",
              help = "number of TFs will be labeled"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "output directory"),
  make_option(c("-d", "--deg"), type = "character",
              help = "input pathway of differentially expressed genes file")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
###setting parameters
inputdir <- opt$inputdir
clinc <- opt$clinc
n <- opt$top_n
outdir <- opt$outdir
inputdeg <- opt$deg

###read in lisa result
##TFs for up-regulated targets
TF.up <- read.table(file = paste(opt$inputdir,"/",clinc,"/",clinc,"_upRegGenes","/",clinc,"_upRegGenes.txt_chipseq_cauchy_combine_dedup.csv",sep = ""),sep = ",", header = TRUE, row.names = 1)
TF.up.sig.p <- TF.up[,"pval"]
names(TF.up.sig.p) <- TF.up[,"factor"]
TF.up.sig.p <- na.omit(TF.up.sig.p)
TF.up.sig.p <- cbind.data.frame(TF = names(TF.up.sig.p),TF.up.sig.p,rank = c(rep("top.up",n),rep("others",length(TF.up.sig.p)-n)))
##TFs for down-regulated targets
TF.down <- read.table(file = paste(opt$inputdir,"/",clinc,"/",clinc,"_downRegGenes","/",clinc,"_downRegGenes.txt_chipseq_cauchy_combine_dedup.csv",sep = ""),sep = ",", header = TRUE, row.names = 1)
TF.down.sig.p <- TF.down[,"pval"]
names(TF.down.sig.p) <- TF.down[,"factor"]
TF.down.sig.p <- na.omit(TF.down.sig.p)
TF.down.sig.p <- cbind.data.frame(TF = names(TF.down.sig.p),TF.down.sig.p,rank = c(rep("top.down",n),rep("others",length(TF.down.sig.p)-n)))

####read in differentially expressed genes
deg <- read.table(paste(inputdeg,"/",clinc,"/",clinc,"_DESeq2_ConvertID.txt",sep = ""),sep = "\t", header = TRUE)
deg.data <- deg %>% 
  mutate(col = ifelse(padj >= 0.05, "no-significant", 
                      ifelse(log2FoldChange >= 1, yes = "Up", 
                             no = ifelse(log2FoldChange <= -1, yes = "Down","no-DE")))) # %>% mutate(label = ifelse(abs(logFC) >= 2.5 & logCPM >= median(logCPM) & FDR <= 0.01, genes, "")) 
###combine TC pvalue and foldchange from up and down profile 
DEs <- as.character(subset(deg.data,select = gene_name)[,1])  # col %in% c("Down","Up"), 
TF.fc <- cbind.data.frame(merge(TF.up.sig.p, TF.down.sig.p, by= "TF")) 
TF.fc$logFC <- deg.data[match(TF.fc$TF, deg.data$gene_name),"log2FoldChange"]
TF.fc <- na.omit(TF.fc) %>%
  mutate_at(vars(TF),funs(as.character)) %>%
  dplyr::filter(TF %in% DEs) %>%   ######only remain DEGs
  mutate(Up.logP = -log10(TF.up.sig.p),
         Down.logP = -log10(TF.down.sig.p)) %>%
  mutate_at(vars(rank.x,rank.y),funs(as.character)) %>%
  mutate(type = ifelse(rank.x != "others" & rank.y == "others", rank.x,
                       ifelse(rank.x == "others" & rank.y != "others", rank.y,
                              ifelse(rank.x != "others" & rank.y != "others","top.common","others"))))
##vocano plot
TF.fc.new <- TF.fc # %>% mutate(target = ifelse(TF == target, "target", ifelse(TF %in% cofs,"cofactors","others")))
min_value <- ceiling(min(TF.fc$logFC))
mean_value <- 0
max_value <- ceiling(max(TF.fc$logFC))
png(file = paste(outdir,"/",clinc,"/",clinc,"_lisa_up_down_TF.png",sep = ""),res = 300, width = 2000, height = 2000)
ggplot(TF.fc.new, aes(x = Up.logP, y = Down.logP)) + 
  geom_point(aes(color = logFC), size = 4, na.rm = TRUE) + # add gene points
  theme_bw() + # clean up theme
  theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
        axis.text.x=element_text(size=12,face = "bold",hjust=1),
        axis.text.y=element_text(size=12,face = "bold",hjust=1),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "top")+
  xlab(expression(-log[10]("up-regulated-Pval-LISA"))) + # x-axis label
  ylab(expression(-log[10]("down-regulated-Pval-LISA"))) + # y-axis label
  scale_color_gradientn(name = "DEG_logFC",
                        colours = c("navy","white","firebrick3"), 
                        values = rescale(c(min_value,mean_value,max_value)),
                        guide = "colorbar", limits=c(min_value,max_value),breaks=c(min_value,mean_value,max_value))+
  # guides(fill = guide_colorbar(title='DEG_logFC',barwidth = 0.5, barheight = 3))+
  geom_label_repel(data = subset(TF.fc.new,type != "others"),
                   aes(label = TF, fill = factor(type)), alpha = 0.8,# fill = "white",
                   fontface = 'bold',box.padding = 0.4, color = 'black',#label.size = 0.08,
                   point.padding = 0.5,segment.color = 'black')+
  scale_fill_manual(values = setNames(c("pink", "lightblue","gold"), c("top.up","top.down","top.common")), name = "TFs_Type")

