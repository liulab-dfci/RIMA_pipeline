library(patchwork)
library(plyr)
suppressMessages(library(optparse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

###make option list and parse command line
option_list <- list(
  make_option(c("-t","--tide"), type = "character",
              help = "tide result"),
  make_option(c("-s", "--msi"), type = "character",
              help = "msi result"),
  make_option(c("-m", "--meta"), type = "character",
              help = "metasheet"),
  make_option(c("-o", "--output"), type = "character",
              help = "output folder"),
  make_option(c("-p", "--col"), type = "character",
              help = "Whcih phenotype is used to do comparison")
  
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

#meta <- read.csv("/Users/linyang/Documents/Liulab/analysis/metasheet_latest.csv")
#msi<-read.table("/Users/linyang/Documents/Liulab/msi/msi_score.txt",stringsAsFactors = FALSE)
#reading the data
meta <- read.table(opt$meta, sep = ",", stringsAsFactors = FALSE, header = TRUE)
msi <- read.table(paste(opt$msi, "msi_score.txt", sep = ""), stringsAsFactors = FALSE)
colnames(msi)<-c("SampleName","MSI")
#decide whether the pre condition exist

if ("Timing" %in% colnames(meta)) {
 tide_pre <- read.table(paste0(opt$tide, "pre/", "pre_tpm_convertID_batch_tide_score.txt", sep = ""), sep = "\t", header = TRUE)
 tide_post <- read.table(paste0(opt$tide, "post/", "post_tpm_convertID_batch_tide_score.txt", sep = ""), sep = "\t", header = TRUE)
 tide <- rbind(tide_pre, tide_post)
 colnames(tide)[c(1,3)] <- c("SampleName", "Predicted_Response") 
 tide <- tide[c(-2,-6)]
} else {
  tide <- read.table(paste0(opt$tide, "tpm_convertID_batch_tide_score.txt", sep = ""), sep = "\t", header = TRUE)
}

#Combine the MSI score
tide_table <- merge(tide, msi, by = 1, all = FALSE)
tide_table <- tide_table[c(1,2,3,7,8,9,10,11,12,4,5,6,13)]
tide_table$Predicted_Response <- ifelse(tide_table$Predicted_Response == "True", "R","NR")

#whether the data have Responder info
if ("Timing" %in% colnames(meta)) {
  tide_table <- merge(meta[c("SampleName","Responder")], tide_table, by = 1, all = FALSE)
  colnames(tide_table)[2] <- "Actual_Response"
}

write.table(tide_table, file = paste(opt$output,"tide_table.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

#Biomarker evaluation plot
final_table <- merge(tide_table, meta, by = 1, all = FALSE)
final_table$all <- "z_samples"
final_table <- rename(final_table, c("MSI" = "MSI Score"))
phenotype <- opt$col

bioma <- c("MSI Score", "TIDE", "CD8", "CD274")
new_bioma <- c("Dysfunction", "Exclusion","MDSC","CAF")

biomark_process <- function(bioma, phenotype) {
  if ("Timing"%in%phenotype) {
    #preprocess the timing table
    tmp_table <- final_table[c("SampleName", bioma,"PatName", "Responder","Timing", "all")]
    colnames(tmp_table)[2] <- "biomarker"
    tmp_table$Timing <- ifelse(tmp_table$Timing == "pre" | tmp_table$Timing == "Pre", "Apre", 
                               ifelse(tmp_table$Timing == "post" | tmp_table$Timing == "Post", "Post", "other"))
    compare_list<-list(unique(as.character(tmp_table[,"Timing"])))
    
    p <- ggplot(tmp_table, aes(x=Timing, y=biomarker),color=Timing) +
      #  scale_fill_manual(values = getPalette(colourCount), name = "P
      #  geom_violin(aes(fill=factor(Phenotype)), alpha=0.3, size= 0.6,color="black")+
      geom_boxplot(aes(fill=factor(Timing)), alpha=0.3, size=0.25,outlier.size= -1, width = 0.5) + theme_bw() +
      stat_compare_means(comparisons = compare_list,aes(label = ..p.format..), tip.length = 0.01,size = 2) +
      geom_point(aes(shape = Responder, color = Responder),alpha=0.7) + geom_line(aes(group = PatName, color = Responder), alpha = 0.5) + 
      #  geom_violin(aes(x = all, y=msi), trim=FALSE, size = 0, alpha=0) +
#      geom_boxplot(aes(x = all, y=biomarker, fill=factor(all)), alpha=0.3, size=0.3,outlier.size= -1, width = 0.5)+
      #  geom_violin(aes(x = tem, y=msi), trim=FALSE, size = 0, alpha=0) +
      #  geom_violin(aes(x = all, y=msi, fill=factor(all)), trim=TRUE, size = 0.6, alpha=0.3)+
      #  stat_summary(aes(x = all, y=msi), fun=mean, geom="point", size=1, color = "black") +
      #  geom_violin(aes(x = tem, y=msi), trim=FALSE, size = 0, alpha=0) +
      scale_x_discrete(name ="", labels=c("Apre" = "Pre", "z_samples" = "All Samples")) +
      theme(
        axis.text.x=element_text(size=6, face = "bold"),
        axis.text.y=element_text(size=6, face = "bold"),
        axis.title.y = element_text(size = 6, face = "bold"),
        legend.position="bottom",
        legend.title = element_text(size=8, face = "bold"),
        legend.text = element_text(size=8, face = "bold")
      ) + labs(y = bioma) + guides(fill=FALSE)
  } else {
    tmp_table <- final_table[c("SampleName", bioma,"PatName", phenotype, "all")]
    colnames(tmp_table)[2] <- "biomarker"
    colnames(tmp_table)[4] <- "phenotype"
    #tmp_table$Timing <- ifelse(tmp_table$Timing == "pre" | tmp_table$Timing == "Pre", "Apre", 
    #                           ifelse(tmp_table$Timing == "post" | tmp_table$Timing == "Post", "Post", "other"))
    compare_list<-list(unique(as.character(tmp_table[, "phenotype"])))
    p <- ggplot(tmp_table, aes(x=phenotype, y=biomarker)) +
      #  scale_fill_manual(values = getPalette(colourCount), name = "P
      #  geom_violin(aes(fill=factor(Phenotype)), alpha=0.3, size= 0.6,color="black")+
      geom_boxplot(aes(fill=factor(phenotype)), alpha=0.3, size=0.25,outlier.size= -1, width = 0.5) + theme_bw() +
      stat_compare_means(comparisons = compare_list,aes(label = ..p.format..), tip.length = 0.01,size = 2) +
#      geom_point(aes(shape = Responder, color = Responder),alpha=0.7) + geom_line(aes(group = PatName, color = Responder), alpha = 0.5) + 
      # geom_violin(aes(x = all, y=biomarker), trim=FALSE, size = 0, alpha=0) +
       geom_boxplot(aes(x = all, y=biomarker, fill=factor(all)), alpha=0.3, size=0.3,outlier.size= -1, width = 0.5)+
      #  geom_violin(aes(x = tem, y=msi), trim=FALSE, size = 0, alpha=0) +
      #  geom_violin(aes(x = all, y=msi, fill=factor(all)), trim=TRUE, size = 0.6, alpha=0.3)+
      #  stat_summary(aes(x = all, y=msi), fun=mean, geom="point", size=1, color = "black") +
      #  geom_violin(aes(x = tem, y=msi), trim=FALSE, size = 0, alpha=0) +
      scale_x_discrete(name ="", labels=c("z_samples" = "All Samples")) +
      theme(
        axis.text.x=element_text(size=6, face = "bold"),
        axis.text.y=element_text(size=6, face = "bold"),
        axis.title.y = element_text(size = 6, face = "bold"),
        #legend.position="bottom",
        legend.title = element_text(size=8, face = "bold"),
        legend.text = element_text(size=8, face = "bold")
      ) + labs(y = bioma) + guides(fill=FALSE)
  }
}

plot_list<-lapply(bioma, biomark_process, phenotype)
plot_list2 <- lapply(new_bioma, biomark_process, phenotype)

combined_plot <- c(plot_list, plot_list2)

png(filename = paste0(opt$output,"Biomarker-Evaluation-by-",phenotype,"_mqc.png", sep = ""),
    res = 300 ,width= 3000, height = 1500, pointsize = 12)
ggarrange(plotlist = combined_plot, ncol = 4, nrow = 2, common.legend = TRUE, legend = "bottom")
dev.off()
