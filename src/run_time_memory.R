library(lubridate)
library(dplyr)
library(reshape2)
library(ggplot2)
library(optparse)
library(stringr)
option_list <- list(
  make_option(c("-i","--input"),type="character", 
              help="Input path of benchmarks. [Required]"))
opt_parser <- OptionParser(option_list=option_list);
opts <- parse_args(opt_parser);
input_benchmarks <- unlist(strsplit(opts$input,","))

###the order of modules
modules <- c("star_align","microbiota","optitype",
                      "rseqc_genebody_cvg","rseqc_read_qual","rseqc_read_distrib","rseqc_tin_score","rseqc_insert_size","rseqc_junction_saturation",
                      "quantification","varscan","vep","pvacseq","fusion","msisensor","trust4",
                      "deseq2","lisa","ssgsea","WGCNA","deconv","tide_score","msi_est_score")
# input_dir <- "benchmarks/consumption/"

Collect_consumption <- lapply(input_benchmarks, function(benchmark){
  # benchmark <- "star_align"
  ss.bck <- read.table(file = benchmark, header = TRUE, sep = "\t", check.names = TRUE)
  ss.bck.cvtime <- ss.bck %>% mutate(Time_minutes = as.numeric(hms(h.m.s))/60) %>%
    mutate(Memory_G = as.numeric(max_rss)/1000) %>%
    mutate(Run = benchmark)
  return (ss.bck.cvtime)
})
all.consumptions <- do.call("rbind", Collect_consumption) %>%
  subset(select = c( "Time_minutes","Memory_G","Run" )) %>%
  group_by(Run) %>%
  dplyr::summarise(Median_time_minutes = median(Time_minutes), Median_memory_G = median(Memory_G))
all.consumptions$Run <- str_match(str_subset(unlist(strsplit(all.consumptions$Run,"/")),".benchmark"),"([a-zA-Z_0-9]+)(.benchmark)")[,2] 
all.consumptions.melt <- melt(all.consumptions, id.vars = "Run")
all.consumptions.melt$Run <- factor(all.consumptions.melt$Run, levels = modules)
png("images/consumption/run_time_memory.png",res = 300, width = 2000, height = 2000)
ggplot(all.consumptions.melt,mapping = aes(x = Run, y = value, group = variable)) +
  geom_line(aes(color = variable), size = 1, show.legend = TRUE) +
  geom_point(aes(color = variable), size = 2, show.legend = TRUE) +
  scale_color_brewer(palette = "Set1", guide = 'legend') + 
  # geom_line(data = subset(df1, variable %in% c('norm_ratio')), aes(color = variable), col = 'red', size = 1) +
  # geom_point(data = subset(df1, variable %in% c('norm_ratio')), aes(color = variable), col = 'red', size = 2) +
  # scale_y_continuous(sec.axis = sec_axis(trans = ~ . * 0.3,
  #                                        name = 'Max Memory(G)')) +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=10, face="bold", angle = 30, hjust = 1), 
        axis.title.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10, face="bold", color = 'black'),
        axis.text.y.left  = element_text(size=10, face="bold", color = "black"),
        axis.title.y.left = element_blank(),
        axis.text.y.right = element_text(size=10, face="bold", color = "#377EB8"),
        axis.title.y.right = element_text(size=10, face="bold", color = "#377EB8"),
        axis.title.y = element_text(size=10, face="bold"),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position = "top",
        legend.text = element_text(size = 16),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.width=unit(2,"line"),
        legend.key.height=unit(2,"line"),
        legend.title = element_blank(),
        strip.text.x = element_text(size=10, face="bold", color = "black", angle = 0)) +
  # ylab('Median_Time(min)') + 
  xlab('Run')
dev.off()


# ggplot(run.time,aes(x=Run, y=value, fill=variable))+
#   geom_boxplot(outlier.size = 0.5)+
#   ylab("")+
#   scale_fill_manual(values = c("#377EB8","#E41A1C"))+
#   theme_bw()+
#   theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
#         axis.text.x=element_text(angle=30,size=12,face = "bold",hjust=1),
#         axis.text.y = element_text(size = 12,face = "bold"),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 12,face = "bold"),
#         legend.text = element_text(size = 12),
#         legend.position='bottom')
