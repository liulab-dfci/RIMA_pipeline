#!/usr/bin/env Rscript

#dependencies
suppressMessages(library(dplyr))
suppressMessages(library(ggrepel))
suppressMessages(library(ggnewscale))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape))
suppressMessages(library(stringr))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
               help="merged fusion prediction file", metavar="character"),
  make_option(c("-f", "--pradafusion"), type="character",default=NULL,
              help="prada fusion prediction file", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="metasheet file", metavar="character"),
  make_option(c("-e", "--expression"), type="character", default=NULL, 
              help="expression file", metavar="character"),
  make_option(c("-p", "--phenotype"), type="character", default=NULL, 
              help="showing limited fusion number", metavar="int"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-a", "--annot"), type="character", default=NULL, 
              help="annotation file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


###read in input data
fusion <- read.table(file = opt$input, header = TRUE)
pheno <- opt$phenotype
meta <- read.table(file = opt$meta, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
meta <- subset(meta, meta[,pheno] != 'NA')
gcat <- fread(file=opt$annot, header = TRUE, sep = "\t")
pradafusion <- read.table(file = opt$pradafusion, sep="\t")
expr <- read.table(file=opt$expression, header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
expr <- expr[,rownames(meta)]
outdir <- opt$outdir



#-------------------star fusion plot-------------------#
### read in gene category from cancer driver gene list
colnames(gcat) <- gsub(" ","_", colnames(gcat))
gcat.new <- gcat %>% mutate(type = ifelse(Is_Oncogene == "Yes" & Is_Tumor_Suppressor_Gene == "No", "Oncogene",
                                          ifelse(Is_Oncogene == "No" & Is_Tumor_Suppressor_Gene == "Yes", "Tumor_Suppressor",
                                                 ifelse(Is_Oncogene == "Yes" & Is_Tumor_Suppressor_Gene == "Yes", "Both", "Cancer_Driver"))))
oncogene <- data.frame(subset(gcat.new, type == "Oncogene", select = "Hugo_Symbol"))[,1]
tsg <- data.frame(subset(gcat.new, type == "Tumor_Suppressor", select = "Hugo_Symbol"))[,1]
both <- data.frame(subset(gcat.new, type == "Both", select = "Hugo_Symbol"))[,1]
driver <- data.frame(subset(gcat.new, type == "Cancer_Driver", select = "Hugo_Symbol"))[,1]

###remove fusion pairs based on the criteria: 
###0. remove the false positive fusion genes 
fusion <- subset(fusion, LargeAnchorSupport != "NO_LDAS")
###1. with TCR or immunoglobulin genes; 2. mitochondrial genes; 3. same gene or paralog gene; 4. fusions in normals; 5. FFPM>0.1
filter.fusion <- fusion %>% mutate_at(c("LeftGene", "RightGene"), as.character) %>%
  filter(!grepl("^TRAV|^TRBV",LeftGene)) %>% filter(!grepl("^TRAV|^TRBV",RightGene)) %>%
  filter(!grepl("^IGH|^IGL|^IGK",LeftGene)) %>% filter(!grepl("^IGH|^IGL|^IGK",RightGene)) %>%
  filter(!grepl("^MT-",LeftGene)) %>% filter(!grepl("^MT-",RightGene)) %>%
  filter(LeftGene != RightGene) %>%
  filter(FFPM > 0.1)

### add two columns recording gene type and gene fusion gene expression
TypeExpr <- apply(filter.fusion,1, function(x) {
  name <- as.character(x[["FusionName"]])
  ss <- as.character(x[["sample"]])
  fgs <- unlist(strsplit(name,"--"))
  gexpr <- expr[fgs, ss]
  ## fusion genes in essential gene list
  g1_type <- ifelse(fgs[1] %in% oncogene, "Oncogene", ifelse(fgs[1] %in% tsg, "Tumor_Suppressor",ifelse(fgs[1] %in% both, "Both", 
                    ifelse(fgs[1] %in% driver, "Cancer_Driver", "Others"))))
  
  g2_type <- ifelse(fgs[2] %in% oncogene, "Oncogene", ifelse(fgs[2] %in% tsg, "Tumor_Suppressor",ifelse(fgs[2] %in% both, "Both", 
                    ifelse(fgs[2] %in% driver, "Cancer_Driver", "Others"))))
  mapl <- length(grep("Others", c(g1_type,g2_type), value = TRUE))
  
    if(mapl == 2){
    tmp <- "mean"  ## two fusion genes are both not essential gene
    value <- mean(gexpr)
    type <- "Others"
    pos <- "Others"
  }
  if(mapl == 1){
    tmp <- setdiff(c(1,2),match("Others", c(g1_type,g2_type))) ## one fusion gene is essential gene, get the essential gene
    value <- gexpr[tmp]
    type <- ifelse(tmp == 1, g1_type, g2_type)
    if(tmp == 1){
      pos = "left_target"
    } else{
      pos = "right_target"
    }
  }
  if(mapl == 0){
     ## two fusion genes are essential genes
    value <- gexpr
    type <- c(g1_type,g2_type)
    pos <- c("left_target","right_target")
  }
  res <- data.frame(Gene = name, Type = type, Expression = value, Pos = pos, Sample = ss, Phenotype = meta[ss,pheno])
  print(mapl)
  return(res)
})
filter.fusion.new <- data.frame(do.call("rbind", TypeExpr) %>% 
  group_by(Gene, Type, Pos, Sample, Phenotype) %>% summarise(Expression = mean(Expression)) %>%
  mutate(Group = "Fusion"))
filter.fusion.genes <- unique(unlist(lapply(as.character(filter.fusion.new$Gene), function(x) unlist(strsplit(x,"--")))))
merge.fusion.df <- rbind(filter.fusion.new)
merge.fusion.df$Gene <- as.character(merge.fusion.df$Gene)
merge.fusion.df$Type <- as.character(merge.fusion.df$Type)
merge.fusion.df$Pos <- as.character(merge.fusion.df$Pos)

if(pheno == "Timing") {
  merge.fusion.df$Phenotype <- ifelse(merge.fusion.df$Phenotype == "Pre", "APre", 
                                      ifelse(merge.fusion.df$Phenotype == "Post", "Post", "Unknown"))
}

#------------------Prada fusion plot-----------------------------
###read in input data
colnames(pradafusion) <- c("Gene1", "Gene2", "Transcript1_Len", "Transcript2_Len", "Identity", "Align_Len", 
                           "Evalue", "BitScore")
###process the file for plotting
pradafusion<-na.omit(pradafusion)
pradafusion$BitScore <- log(pradafusion$BitScore, 10)

pradafusion<-pradafusion%>%
  unite(FusionName, Gene1, Gene2, sep="--")
highlight_genes <- pradafusion %>% 
  filter((Evalue < 0.05 & Align_Len > 1000) | (Evalue <0.05 & BitScore > 2))


png(paste(outdir, pheno, "_prada_homology.png", sep = ""), width = 1600, height = 1400, res = 300)
pradafusion %>% 
  ggplot(aes(x=BitScore,y=Evalue,size=Align_Len)) + 
  geom_point(alpha=0.3) +
  geom_point(data=highlight_genes, 
             aes(x=BitScore,y=Evalue), 
             color="#e6550d", alpha=0.6)+
  geom_label_repel(data = highlight_genes, aes(label = FusionName), size = 1.5) +
  labs(x ="log10(BitScore)", y = "Evalue") + theme_bw() + 
  theme(axis.text=element_text(angle=0,size=8,face = "bold",hjust=0.5),
      axis.title = element_text(size = 8,face = "bold"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8))
dev.off()

#remove the homogous fusion gene pair 

merge.fusion.df <- merge(merge.fusion.df, highlight_genes, by = 1, all = TRUE)
merge.fusion.df <- subset(merge.fusion.df, is.na(merge.fusion.df$BitScore))
merge.fusion.df <- merge.fusion.df[c(1:7)]
merge.fusion.df <- na.omit(merge.fusion.df)
write.table(merge.fusion.df, paste(outdir,pheno,"_fusion_gene_table.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)



### extract expression of oncogene,tsg and protein kinase genes
png(paste(outdir, pheno,"_fusion_gene_plot.png", sep = ""), width = 1600, height = 1400, res = 300)
#pdf(paste(outdir, pheno,"_fusion_gene_plot.pdf", sep = ""), width = 8, height = 7)
ggplot(merge.fusion.df, aes(x=Phenotype, y=log10(Expression+1))) + 
  geom_violin(trim=TRUE) +
  geom_jitter(aes(color = Type),shape=16#, size = 4#, position=position_jitter(0.15)
  )+
  theme_bw()+
  scale_colour_manual(name = "Type", values = c(Oncogene="#a6bddb",Tumor_Suppressor="#99d8c9", Both="#4DAF4A",Cancer_Driver="grey"))+
  new_scale_color() +
  geom_label_repel(data = subset(merge.fusion.df, Type != "Others"), 
                   aes(label = Gene, fill = Type), #alpha = 0.5,#fill = "white", 
                   fontface = 'bold',label.size = 0.15, family = "serif",size = 2.5
  )+
  scale_fill_manual(values = c(Oncogene="#a6bddb",Tumor_Suppressor="#99d8c9", Both="#4DAF4A",Cancer_Driver="grey"),guide = FALSE)+
  scale_colour_manual(name = "Target", values = c(left_target="#636363",right_target="#e6550d"))+
  scale_x_discrete(name ="", labels=c("APre" = "Pre")) + 
  theme(axis.text.x=element_text(angle=0,size=8,face = "bold",hjust=0.5),
        axis.text.y=element_text(size=8,face = "bold",hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8,face = "bold"),
        legend.text = element_text(size = 8),
	legend.title = element_text(size = 8),
        strip.text = element_text(size = 8, face = "bold"),
        legend.position = "right") 
dev.off()




##fusion box plot
type <- c("Oncogene", "Tumor_Suppressor", "Cancer_Driver", "Others")

ta <- merge.fusion.df

#preprocess the data
processed_ta <- NULL
for (i in unique(ta$Phenotype)) {
  tmp_ta <- ta[ta$Phenotype == i,]
  tmp_ta_others <- tmp_ta[tmp_ta$Type == "Others",]
  tmp_ta_sig <- tmp_ta[tmp_ta$Type != "Others",]

  #remove the duplicates for others pairs, and calculating their average expression
  tmp_ta_others <- tmp_ta_others %>% group_by(Gene) %>% mutate(Expression = mean(Expression))
  tmp_ta_others <- tmp_ta_others[!duplicated(tmp_ta_others$Gene),]
  tmp_ta_others <- tmp_ta_others %>% group_by(Gene) %>% mutate(ID=1:n())

  #group multiple hits for signifcant gene pairs
  tmp_ta_sig <- tmp_ta_sig %>% group_by(Gene) %>% mutate(ID=1:n())

  tmp_ta <- rbind(tmp_ta_others, tmp_ta_sig)

  processed_ta <- rbind(processed_ta, tmp_ta)
}

#re-order the figure followed: Oncogenes, Tumor suppressor, Cancer driver, Others
reorder_ta <- NULL
for (t in type) {
  tmp <- processed_ta[processed_ta$Type ==t,]
  tmp <- tmp[order(tmp$Expression, decreasing = TRUE),]
  reorder_ta <- rbind(reorder_ta, tmp)
}


reorder_ta$Gene <- factor(reorder_ta$Gene, levels = rev(unique(reorder_ta$Gene)))
reorder_ta <- reorder_ta[reorder_ta$Expression != 0,]


####labeling
label <- reorder_ta[reorder_ta$Pos != "Others",]
test <- str_split_fixed(label$Gene, "--", 2)
label$left <- test[,1]
label$right <- test[,2]
label$Pos <- ifelse(label$Pos == "left_target", label$left, ifelse(label$Pos == "right_target", label$right, "both"))
label <- label[colnames(reorder_ta)]
label <- label[!duplicated(label$Gene),]
###


###define the color
color_ta <- data.frame(terms = type, color = c("#e34a33", "#a6bddb", "#fc9272", "#bdbdbd"))
iterms <- as.character(color_ta[color_ta$terms %in% unique(reorder_ta$Type),]$terms)
color <- as.character(color_ta[color_ta$terms %in% iterms,]$color)

png(paste(outdir, "fusion_boxplot.png", sep = ""), width = 2000, height = 1400, res = 300)
p <- ggplot(reorder_ta, aes(x = Gene, y = Expression, fill = Type, group=factor(ID))) +
       geom_bar(stat = "identity",position=position_dodge(), color = "#636363") +
       facet_grid(Phenotype ~ ., scales = "free_y", space = "free_y") + #, scales = "free_y", space = "free_y") +
       #facet_wrap(~Phenotype, ncol = 1, strip.position="right", scales = "free_y") + coord_flip() + theme_bw() +
       geom_text(data = label, aes(label=Pos), position = position_dodge(.9), hjust = -0.1, size = 2) +
       coord_flip() + theme_bw() +
       scale_fill_manual(breaks = iterms, values = color) +
       labs(x = "Gene Pairs") +
        theme(axis.text.x= element_text(size=8),
              axis.text.y = element_text(size=6),
              axis.title = element_text(size = 8,face = "bold"),
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 6),
              strip.text = element_text(size = 8, face = "bold"))
p + expand_limits(y = c(min(reorder_ta$Expression), max(reorder_ta$Expression)*1.1))
dev.off()




