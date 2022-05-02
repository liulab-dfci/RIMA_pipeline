#ta <- merge.fusion.df

fusion_bar <- function(ta) {
  type <- c("Both", "Oncogene", "Tumor_Suppressor", "Cancer_Driver", "Others")
  color_ta <- data.frame(terms = type, color = c("#e333cc","#e34a33", "#a6bddb", "#fc9272", "#bdbdbd"))
  
  #ta <- ta[ta$Sample %in% meta[!is.na(meta$Screen_R_vs_NR),][[1]],]
  
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
    if(nrow(tmp_ta_sig) == 0) {
      tmp_ta_sig <- NULL
    } else {
      tmp_ta_sig <- tmp_ta_sig %>% group_by(Gene) %>% mutate(ID=1:n())
    }
    #tmp_ta_sig <- tmp_ta_sig %>% group_by(Gene) %>% mutate(ID=1:n())
    
    tmp_ta <- rbind(tmp_ta_others, tmp_ta_sig)
    
    processed_ta <- rbind(processed_ta, tmp_ta)
  }
  
  #re-order the figure followed: Oncogenes, Tumor suppressor, Cancer driver, Others
  reorder_ta <- NULL
  for (t in type) {
    tmp <- processed_ta[processed_ta$Type ==t,]
    tmp <- tmp[order(tmp$Expression, decreasing = TRUE),]
    if (t == "Others" & nrow(tmp) > 20) {
      tmp_other <- NULL
      for (i in unique(tmp$Phenotype)) {
        tmp_1 <- tmp[tmp$Phenotype == i,]
        tmp_1 <- head(tmp_1, 10)
        tmp_other <- rbind(tmp_other, tmp_1)
      }
      tmp <- tmp_other
    }
    reorder_ta <- rbind(reorder_ta, tmp)
  }
  
  #calculate the number of gene pairs that are "Others" type
  npairs <- nrow(reorder_ta[reorder_ta$Type != "Others",])
  if (npairs > 30) {
    reorder_ta <- reorder_ta[reorder_ta$Type != "Others",]
    
  }
  
  
  reorder_ta$Gene <- factor(reorder_ta$Gene, levels = rev(unique(reorder_ta$Gene)))
  reorder_ta <- reorder_ta[reorder_ta$Expression != 0,]
  
  #seperate the gene pairs if the same gene pair belongs to different categories
  reorder_ta <- reorder_ta %>% group_by(Gene) %>% mutate(nTypes = length(unique(Type)))
  
  #seperate the gene pairs if the same gene pair belongs to different categories
  if(2 %in% reorder_ta$nTypes) {
    tmp_phe_r2 <- reorder_ta[reorder_ta$nTypes != 2,]
    tmp_phe_2 <- reorder_ta[reorder_ta$nTypes == 2,]
    
    tmp_phe_2$Gene <- ifelse(tmp_phe_2$Pos == "left_target", paste0(tmp_phe_2$Gene, "(1)"),  paste0(tmp_phe_2$Gene, "(2)"))
    reorder_ta <- rbind(tmp_phe_r2, tmp_phe_2)
  }
  
  
  ####labeling
  label <- reorder_ta[reorder_ta$Pos != "Others",]
  test <- str_split_fixed(label$Gene, "--", 2)
  label$left <- test[,1]
  label$right <- test[,2]
  label$Pos <- ifelse(label$Pos == "left_target", label$left, ifelse(label$Pos == "right_target", label$right, "both"))
  label <- label[colnames(reorder_ta)]
  
  
  final_label <- NULL
  phe <- unique(label$Phenotype)
  for (p in phe) {
    tmp_phe <- label[label$Phenotype == p ,]
    tmp_phe <- tmp_phe[!duplicated(tmp_phe$Gene),]
    
    final_label <- rbind(final_label, tmp_phe)
    final_label$Pos <- sub("\\(2)", "", final_label$Pos)
  }
  #label <- label[!duplicated(label$Gene),]
  ###
  
  
  iterms <- as.character(color_ta[color_ta$terms %in% unique(reorder_ta$Type),]$terms)
  color <- as.character(color_ta[color_ta$terms %in% iterms,]$color)
  
  
  final_ta <- NULL
  for (t in type) {
    tmp <- reorder_ta[reorder_ta$Type ==t,]
    tmp <- tmp[order(tmp$Expression, decreasing = TRUE),]
    final_ta <- rbind(final_ta, tmp)
  }
  final_ta$Gene <- factor(final_ta$Gene, levels = rev(unique(final_ta$Gene)))
  final_ta <- final_ta[final_ta$Expression != 0,]
  
  mpairs <- final_ta[final_ta$ID > 8,]
  mpairs_list <- NULL
  for (m in unique(mpairs$Gene)) {
    mpairs_tmp <- mpairs[mpairs$Gene == m,]
    
    if(length(unique(mpairs_tmp$Phenotype)) > 1) {
      for (ph in unique(mpairs_tmp$Phenotype)) {
        mpairs_tmp_phe <- mpairs_tmp[mpairs_tmp$Phenotype == ph,]
        mpairs_extract <- data.frame(Gene = m, num = max(mpairs_tmp_phe$ID), phenotype = ph)
        mpairs_list <- rbind(mpairs_list, mpairs_extract)
      }
    } else {
      mpairs_extract <- data.frame(Gene = m, num = max(mpairs_tmp$ID), phenotype = unique(mpairs_tmp$Phenotype))
      mpairs_list <- rbind(mpairs_list, mpairs_extract)
    }
  }
  
  p <- ggplot(final_ta, aes(x = Gene, y = Expression, fill = Type, group=factor(ID))) +
    geom_bar(stat = "identity",position=position_dodge(), color = "#636363") +
    facet_grid(Phenotype ~ ., scales = "free_y", space = "free_y") + #, scales = "free_y", space = "free_y") +
    #facet_wrap(~Phenotype, ncol = 1, strip.position="right", scales = "free_y") + coord_flip() + theme_bw() +
    geom_text(data = final_label, aes(label=Pos), position = position_dodge(.9), hjust = -0.1, size = 2) +
    coord_flip() + theme_bw() +
    scale_fill_manual(breaks = iterms, values = color) +
    labs(x = "Gene Pairs") +
    theme(axis.text.x= element_text(size=8),
          axis.text.y = element_text(size=6),
          axis.title = element_text(size = 8,face = "bold"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          strip.text = element_text(size = 8, face = "bold"))
  p + expand_limits(y = c(min(final_ta$Expression), max(final_ta$Expression)*1.1))
  
  return(p)
}
