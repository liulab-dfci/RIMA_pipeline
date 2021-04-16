library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)
library(ggseqlogo)
library(optparse)

###options
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input expression file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-m", "--multiqc"), type="character", default=NULL,
              help="multiqc directory", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="file directory", metavar="character"),
  make_option(c("-s", "--meta"), type="character", default=NULL,
              help="metasheet.scv", metavar="character"),
  make_option(c("-c", "--condition"), type="character", default=NULL,
              help="penotype of condition", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###parameters setting
input <- opt$input
outdir <- opt$outdir
condition <- opt$col
meta <- read.delim("/Users/linyang/Dropbox/RIMA/pvcaseq/new/sub_metasheet.csv", sep = ",")
meta <- meta[c("SampleName", condition)]
colnames(meta)[2] <- "condition"

##############################################################################
### Reading matrix
read_err <- tryCatch({
  input_table = read.table("/Users/linyang/Dropbox/RIMA/pvcaseq/new/Merged.filtered.condensed.ranked.addSample.tsv", header=TRUE,sep = '\t')},
  error = function(e) {return(paste("ERROR:",conditionMessage(e),"\n"))})

input_table <- merge(input_table, meta, by = 1, all = FALSE)
###generate figures only when input is not empty, otherwise output empty figures
if(!grepl("ERROR", read_err)){
  input_table$is_novel = ifelse(input_table$Median.Fold.Change>1 & input_table$Median.MT.Score <500, "novel","non_novel")
  COLORS<- brewer.pal(n = 3, name = "Set1")
  #########---------------------write pvacseq epitopes prediction file----------------###########
  #print(input_table)
  write.table(input_table, file = paste(opt$file, "pvacseq_epitope_prediction.txt",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
  
  #########---------------------count the number of epitopes in each patient----------------###########
  pat.epit <- input_table %>% group_by(Sample,is_novel) %>% dplyr::summarise(count = n())
  Pwidth <- 800+80*length(unique(pat.epit[,1]))
  png(file = paste(outdir, "Patient_count_epitopes_plot.png",sep = ""), res = 300, width = Pwidth, height = 1000)
  p1 <- ggplot(pat.epit, aes(x = Sample,y=count,fill = is_novel)) +
    geom_bar(stat = 'identity', position = 'dodge')+
    theme_minimal()+
    scale_fill_manual(values = brewer.pal(n=8, name = "Set1")[2:1])+
    ylab("#Neoepitopes")+xlab("Patient")+
    theme(legend.position = "top",
          plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_text(angle=70, size=8,face = "bold",hjust=1),
          axis.text.y = element_text(size = 8,face = "bold"),
          axis.title.x = element_text(size = 10,face = "bold"),
          axis.title.y = element_text(size = 10,face = "bold"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8))
  print(p1)
  dev.off()
  #########---------------------the ratio of mutated genes----------------###########
  sub <- subset(input_table, is_novel == "novel")
  
  final_sub <- NULL
  for (i in unique(sub$Gene.Name)) {
    tmp <- subset(sub, sub$Gene.Name == i)
    tmp$count <- nrow(tmp)
    final_sub <- rbind(final_sub, tmp)
    
  }
  
  final_sub$ratio <- final_sub$count / nrow(final_sub)
  
  final_ta <- NULL
  for (c in unique(final_sub$condition)) {
    final_tmp <- subset(final_sub, condition == c)
    final_tmp <- final_tmp[!duplicated(final_tmp$Gene.Name),]
    final_ta <- rbind(final_ta, final_tmp)
  }
  
  final_ta <- final_ta[order(final_ta$Corresponding.Fold.Change, decreasing = TRUE),]
  final_ta$Gene.Name <- factor(final_ta$Gene.Name, levels = unique(final_ta$Gene.Name))
  p <- ggplot(final_ta, aes(x = Gene.Name, y = Median.Fold.Change)) + 
    geom_histogram(stat = "identity", aes(fill = ratio))+facet_grid(~condition) + 
    theme_bw() + scale_color_gradient(low = "blue", high = "red") + 
    theme(plot.margin = unit(c(0.5,0.3,0.3,0.65),"cm"),
          plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x= element_text(angle=70,size=8,face = "bold",hjust=1),
          axis.text.y=element_text(size=10,face = "bold",hjust=1),
          axis.title.x = element_text(size = 10, face = "bold", angle = 0.5),
          axis.title.y = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6))
  
  ##########--------------------predicted epitope Affinity scatter plot------------------##############
  png(file = paste(outdir, "epitopes_affinity_plot.png",sep = ""), res = 300, width = 1600, height = 1200)
  text <- subset(input_table, is_novel == "novel")
  p2 <- ggplot(input_table, aes(y = Median.MT.Score, x = Median.WT.Score,colour = condition)) +
    geom_point()+ geom_text(data = text, aes(label=Gene.Name), size = 2, hjust=0, vjust=0) +
    xlab("Normal IC50 binding affinity(nM)")+ylab("Tumor IC50 binding affinity(nM)")+
    theme_minimal()+
    scale_x_log10() +
    scale_y_log10()+
    geom_hline(aes(yintercept=500)) +
    geom_vline(aes(xintercept=500))+
    geom_abline(intercept=0,slope=1 ,linetype='dashed')+ 
    annotate(geom="text", x=250, y=5, label="x=500")+
    annotate(geom="text", x=50, y=450, label="y=500")+
    #scale_color_manual(values = c(novel = "red", non_novel = "blue"))+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_text(size=8,face = "bold",hjust=1),
          axis.text.y = element_text(size = 8,face = "bold"),
          axis.title.x = element_text(size = 10,face = "bold"),
          axis.title.y = element_text(size = 10,face = "bold"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8))
  print(p2)
  dev.off()
  
  ##########--------------------the HLA allele fraction for patients(how many patients have this HLA?)------------------##############
  hla.fraction <- input_table %>% subset(select=c(Sample,HLA.Allele,is_novel)) %>% distinct() %>%
    group_by(HLA.Allele, is_novel) %>% dplyr::summarise(count = n()) %>% 
    group_by(HLA.Allele, is_novel) %>% mutate(fraction = count/sum(count)) %>%
    mutate(HLA_type = substring(HLA.Allele, 5,5))
  hla.fraction$HLA.Allele <- factor(hla.fraction$HLA.Allele,levels = as.character(unique(hla.fraction$HLA.Allele)))
  g1 <- ggplot(hla.fraction, aes(x = HLA.Allele, y = fraction, fill = HLA_type))+
    geom_bar(aes(alpha = is_novel),stat="identity", position=position_dodge())+
    scale_fill_manual(values = COLORS)+
    scale_alpha_manual(values = c('novel' = 0.9,'non_novel' = 0.4), name = "Type")+
    scale_y_reverse()+
    theme_minimal()+
    ylab("Patient Fraction")+xlab("")+
    theme(legend.position = "none",
          plot.margin = unit(c(0.2,0.3,0.3,0.4),"cm"),
          plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_text(angle=70,size=12,face = "bold",hjust=1),
          axis.text.y=element_text(size=12,face = "bold",hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"))
  
  #########-------------the number of epitopes for each HLA(each HLA can be binded with how many epitopes)------------------##############
  hla.epit <- input_table %>% 
    group_by(HLA.Allele,is_novel) %>% 
    dplyr::summarise(count = n()) %>% 
    mutate(HLA_type = substring(HLA.Allele, 5,5))
  hla.epit$HLA.Allele <- factor(hla.epit$HLA.Allele, levels = as.character(unique(hla.fraction$HLA.Allele)))
  g2 <- ggplot(hla.epit, aes(x = HLA.Allele, y = count, fill = HLA_type))+
    geom_bar(aes(alpha = is_novel),stat="identity", position=position_dodge())+
    scale_fill_manual(values = COLORS)+
    scale_alpha_manual(values = c('novel' = 0.9,'non_novel' = 0.4), name = "Type")+
    theme_minimal()+
    ylab("#Neoepitopes")+xlab("")+
    theme(legend.position = "top")+
    theme(plot.margin = unit(c(0.5,0.3,0.3,0.65),"cm"),
          plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x= element_text(angle=70,size=8,face = "bold",hjust=1),
          axis.text.y=element_text(size=10,face = "bold",hjust=1),
          axis.title.x = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6))
  
  ###merge g1 and g2
  Pwidth <- length(unique(hla.fraction[[1]]))*60+1440
  png(file = paste(outdir, "HLA_epitopes_fraction_plot.png",sep = ""), res = 300, width = Pwidth, height = 2000)
  p3 <- ggarrange(g2,g1,ncol=1,nrow=2)
  print(p3)
  dev.off()
  
  ##########--------------------logo sequence plot for top neo-epitopes------------------##############
  ## select top expressed novel epitopes across patients
  novel.epitopes <- input_table %>% subset(., is_novel == "novel") %>% mutate(epitope_len = nchar(as.character(MT.Epitope.Seq))) 
  top.epitopes <- novel.epitopes %>%
    group_by(epitope_len, Gene.Name) %>% summarise(count = n()) %>% arrange(desc(count)) %>% filter(count > 1) %>%
    group_by(epitope_len) %>% slice(1:5)
  if(nrow(top.epitopes) > 0){
    lens <- unique(top.epitopes$epitope_len)
    epitope_list <- list()
    for(len in lens){
      tmp <- NULL
      for(gene in as.character(unlist(top.epitopes[top.epitopes$epitope_len == len, "Gene.Name"]))){
        # len <- 9
        # gene <- "HLA-B"
        sub.epitopes <- novel.epitopes %>% subset(epitope_len == len & Gene.Name == gene)
        tmp <- rbind(tmp,sub.epitopes)
      }
      name <- paste("Neo-epitope_len_", len, sep = "")
      np <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
      p <- ggplot(tmp) + geom_logo(as.character(tmp$MT.Epitope.Seq), method='bits', seq_type='aa',namespace = np) +
        theme_logo()+
        ggtitle(name)+
        facet_wrap(~Gene.Name, ncol = 1)+
        theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=50,r=50,t=20,b=10),size = 8,face = "bold", colour = "black"))
      epitope_list[[name]] <- p
    }
    # generate seqlogo for nwo-epitopes from frequently predicted genes
    Pwidth = 1500*length(lens)
    png(file = paste(outdir, "neo-epitopes_seqlogo_plot.png",sep = ""), res = 300, width = Pwidth, height = 2000)
    ggarrange(plotlist = epitope_list)
    dev.off()
  }else{
    empty_plot <-   ggplot() + theme_void()
    print("empty neo-epitopes_seqlogo_plot.png due to no top predicted novel epitopes")
    png(file = paste(outdir, "neo-epitopes_seqlogo_plot.png",sep = ""), res = 300)
    print(empty_plot)
    dev.off()
  }
}else{  
  empty_plot <-   ggplot() + theme_void()
  
  print("empty Patient_count_epitopes_plot.png due to empty input")
  png(file = paste(outdir, "Patient_count_epitopes_plot.png",sep = ""), res = 300)
  print(empty_plot)
  dev.off()
  
  print("empty epitopes_affinity_plot.png due to empty input")
  png(file = paste(outdir, "epitopes_affinity_plot.png",sep = ""), res = 300)
  print(empty_plot)
  dev.off()
  
  print("empty HLA_epitopes_fraction_plot.png due to empty input")
  png(file = paste(outdir, "HLA_epitopes_fraction_plot.png",sep = ""), res = 300)
  print(empty_plot)
  dev.off()
  
  print("empty neo-epitopes_seqlogo_plot.png due to empty input")
  png(file = paste(outdir, "neo-epitopes_seqlogo_plot.png",sep = ""), res = 300)
  print(empty_plot)
  dev.off()
  
}


#multiqc output

#multiqc table
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

multiqc_table <- input_table
multiqc_table$HGVSc <- as.character(multiqc_table$HGVSc)
multiqc_table$HGVSc <- substrRight(multiqc_table$HGVSc, 3)

multiqc_table <- multiqc_table[c(1,2,3,5,7,8,9,23,21)]
colnames(multiqc_table)[7] <- c("|__________MT.Epitope.Seq__________|")
write.table(multiqc_table, file = paste(opt$multiqc, "pvacseq_multiqc.txt",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)


#figure
png(file = paste(opt$multiqc, "neoantigen_plot.png",sep = ""),res = 300, width = 3600, height = 1000)
#par(oma=c(1,1,1,1))
ggarrange(p2, p, g2,ncol = 3, nrow = 1)
dev.off()

