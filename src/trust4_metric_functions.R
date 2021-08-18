suppressMessages(library(seqinr))
suppressMessages(library(ape))
suppressMessages(library(msa))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))

##function of BCR clustering 
SeqDist <- function(x, y) {
  if (nchar(x) != nchar(y))
    return(NA)
  x.l = unlist(strsplit(x, ''))
  y.l = unlist(strsplit(y, ''))
  nn1 = length(which(x.l != y.l))
  nn2 = length(which(x.l == '-' & y.l != '-'))
  nn3 = length(which(x.l != '-' & y.l == '-'))
  return(nn1 - nn2 - nn3)
}
SeqDist.AA <- function(x, y) {
  if (nchar(x) != nchar(y))
    return(NA)
  x.l = unlist(strsplit(x, ''))
  y.l = unlist(strsplit(y, ''))
  tmp.vv = which(x.l != '-' & y.l != '-')
  tmp.st.x = min(which(x.l != '-'))
  tmp.ed.x = max(which(x.l != '-'))
  tmp.st.y = min(which(y.l != '-'))
  tmp.ed.y = max(which(y.l != '-'))
  gap.vv = which(x.l == '-' &
                   y.l != '-' | x.l != '-' & y.l == '-')## needs work
  gap.vv = gap.vv[gap.vv >= max(c(tmp.st.x, tmp.st.y)) &
                    gap.vv <= min(c(tmp.ed.x, tmp.ed.y))]
  tmp = sum(diag(BLOSUM62[x.l[tmp.vv], y.l[tmp.vv]]))
  tmp0 = sum(diag(BLOSUM62[x.l[tmp.vv], x.l[tmp.vv]]))
  score = (tmp0 - tmp + length(gap.vv)) / length(tmp.vv)
  return(score)
}
MergeSeqs <- function(seqs, dd) {
  
  tmp.seqs = gsub('-', '', seqs)
  tmp.vv = order(nchar(tmp.seqs))
  seqs = seqs[tmp.vv]
  dd = dd[tmp.vv, ]
  tmp.seqs = tmp.seqs[tmp.vv]
  newseqs = tmp.seqs
  seq.labs = rep(0, length(newseqs))
  nn = length(seq.labs)
  while (T) {
    for (ii in 1:nn) {
      if (length(grep(newseqs[ii], newseqs[-ii])) > 0)
        seq.labs[ii] = 1
    }
    if (sum(seq.labs) == 0)
      break
    vv0 = which(seq.labs == 0)
    newseqs = newseqs[vv0]
    seq.labs = seq.labs[vv0]
    nn = length(newseqs)
    seqs = seqs[vv0]
    dd = dd[vv0, ]
  }
  return(list(SS = seqs, DD = dd))
}
CreateMotifList <- function(mm) {
  tmp.mat = matrix(unlist(strsplit(rep(mm, 8, ''), '')), nrow = 8, byrow =
                     T)
  diag(tmp.mat) = '.'
  mm.list = apply(tmp.mat, 1, paste, collapse = '')
  return(mm.list)
}
MergeMotifs <- function(motif.list) {
  ## Merge motifs by allowing one mismatch
  unique.motifs = c()
  for (mm in motif.list) {
    mm.list = CreateMotifList(mm)
    sign = 0
    for (tmp.mm in mm.list) {
      if (length(grep(tmp.mm, unique.motifs)) > 0) {
        sign = 1
        break
      }
    }
    if (sign == 0)
      unique.motifs = c(unique.motifs, mm)
  }
  return(unique.motifs)
}

BuildBCRlineage <- function(sampleID, Bdata = BCRdata, start=3, end=10) {
  ## Given sample ID, start and end position of complete CDR3,  return all the lineages in the sample
  # sampleID <- "SRR3184301"
  # Bdata = cdr3.bcr.heavy
  # start=3
  # end=10
  tmp.dd.ss = subset(Bdata, sample == sampleID)
  tmp.dd.ss = tmp.dd.ss[!duplicated(tmp.dd.ss[,"CDR3nt"]),]
  if (is.null(dim(tmp.dd.ss)))
    return(NA)
  tmp.comp.vv <- which(tmp.dd.ss[, "is_complete"] == "Y")
  comp.CDR3.ss = data.frame(CDR3aa = tmp.dd.ss[tmp.comp.vv, "CDR3aa"])
  if (length(comp.CDR3.ss) == 0)
    return(NA)
  tmp.tt = table(substr(comp.CDR3.ss$CDR3aa, start, end))
  tmp.tt = sort(tmp.tt, decreasing = T)
  tmp.tt <- tmp.tt[which(nchar(names(tmp.tt))==(end-start+1))]
  tmp.motifs = MergeMotifs(names(tmp.tt))
  count = 0
  BCRlineage = c() ## a list of BCR lineage trees
  kept.motifs = c()
  for (mm in tmp.motifs) {
    mm.list = CreateMotifList(mm)
    tmp.vv.ss = c()
    for (tmp.mm in mm.list) {
      tmp.vv.ss = c(tmp.vv.ss, grep(tmp.mm, tmp.dd.ss$CDR3aa))
    }
    tmp.vv.ss = unique(tmp.vv.ss)
    if (length(tmp.vv.ss) < 2)
      next
    SEQs = unique(as.character(tmp.dd.ss[tmp.vv.ss, "CDR3nt"]))
    #SEQs = SEQs$CDR3nt
    tmp.dd0 = tmp.dd.ss[tmp.vv.ss, ]
    setDF(tmp.dd0)   ###format as dataframe
    rownames(tmp.dd0) = tmp.dd0$CDR3nt   ####same cdr3dna, same cdr3aa, different Ig gene and totaldna
    tmp.dd0 = tmp.dd0[SEQs, ]
    if (length(SEQs) < 3)
      next
    MSAalign = msa(DNAStringSet(SEQs), 'ClustalW')
    seqs = as.character(attributes(MSAalign)$unmasked)
    seqs0 = gsub('-', '', seqs)
    tmp.dd0 = tmp.dd0[match(seqs0, SEQs),]
    tmp = MergeSeqs(seqs, tmp.dd0)
    seqs = tmp$SS
    tmp.dd0 = tmp$DD
    if (is.null(dim(tmp.dd0)))
      next
    nn = nrow(tmp.dd0)
    if (nn <= 3)
      next
    sDist = matrix(0, nn, nn)
    for (ii in 1:nn) {
      for (jj in ii:nn) {
        if (jj == ii)
          next
        tmp.dist = SeqDist(seqs[ii], seqs[jj])
        sDist[ii, jj] = sDist[ii, jj] + tmp.dist
      }
    }
    kept.motifs = c(kept.motifs, mm)
    rownames(tmp.dd0) = NULL
    lineage.obj = list(distMat = sDist,
                       Sequences = seqs,
                       data = tmp.dd0)
    BCRlineage = c(BCRlineage, list(lineage.obj))
    count = count + 1
  }
  names(BCRlineage) = kept.motifs
  return(BCRlineage)
}

##function of computing clonality
getClonality <- function(sampleID, Bdata = BCRdata, start=3, end=10) {
  ## Given sample ID, start and end position of complete CDR3,  return all the lineages in the sample
  # sampleID <- "SRR3184301"
  # Bdata = cdr3.bcr.heavy
  # start=3
  # end=10
  cluster.ID <- c()
  cluster.list <- list()
  tmp.dd.ss = subset(Bdata, sample == sampleID)
  tmp.dd.ss = tmp.dd.ss[!duplicated(tmp.dd.ss[,"CDR3nt"]),]
  if (is.null(dim(tmp.dd.ss)))
    return(NA)
  tmp.comp.vv <- which(tmp.dd.ss[, "is_complete"] == "Y")
  comp.CDR3.ss = data.frame(CDR3aa = tmp.dd.ss[tmp.comp.vv, "CDR3aa"])
  if (length(comp.CDR3.ss) == 0)
    return(NA)
  tmp.tt = table(substr(comp.CDR3.ss$CDR3aa, start, end))
  tmp.tt = sort(tmp.tt, decreasing = T)
  tmp.tt <- tmp.tt[which(nchar(names(tmp.tt))==(end-start+1))]
  tmp.motifs = MergeMotifs(names(tmp.tt))
  count = 0
  BCRlineage = c() ## a list of BCR lineage trees
  kept.motifs = c()
  for (mm in tmp.motifs) {
    mm.list = CreateMotifList(mm)
    
    tmp.vv.ss = c()
    for (tmp.mm in mm.list) {
      tmp.vv.ss = c(tmp.vv.ss, grep(tmp.mm, tmp.dd.ss$CDR3aa))
    }
    tmp.vv.ss = unique(tmp.vv.ss)
    if (length(tmp.vv.ss) < 2)
      next
    SEQs = unique(as.character(tmp.dd.ss[tmp.vv.ss, "CDR3nt"]))
    tmp.dd0 = tmp.dd.ss[tmp.vv.ss, ]
    setDF(tmp.dd0)   ###format as dataframe
    rownames(tmp.dd0) = tmp.dd0$CDR3nt   ####same cdr3dna, same cdr3aa, different Ig gene and totaldna
    tmp.dd0 = tmp.dd0[SEQs, ]
    if (length(SEQs) < 3)
      next
    MSAalign = msa(DNAStringSet(SEQs), 'ClustalW')
    seqs = as.character(attributes(MSAalign)$unmasked)
    seqs0 = gsub('-', '', seqs)
    tmp.dd0 = tmp.dd0[match(seqs0, SEQs),]
    tmp = MergeSeqs(seqs, tmp.dd0)
    seqs = tmp$SS
    tmp.dd0 = tmp$DD
    if (is.null(dim(tmp.dd0)))
      next
    nn = nrow(tmp.dd0)
    if (nn <= 3)
      next
    cluster.ID <- c(cluster.ID,seqs0)
    cluster.freq <- cbind.data.frame(sample = sampleID, CDR3aa = mm, frequency = sum(tmp.dd0$frequency))
    cluster.list[[mm]] <- cluster.freq 
  }
  ###get non-cluster clones
  non.clonal <- subset(Bdata, sample == sampleID & !(CDR3nt %in% cluster.ID),
                       select = c(sample,CDR3aa,frequency))
  ###total cluster clones
  clonal <- do.call("rbind",cluster.list)
  ###combine
  clone.frequency <- rbind(clonal,non.clonal) %>% mutate(normalized_est_clonal_exp = frequency/sum(frequency) )
  ###calculate clonality
  entropy <- -sum(clone.frequency$normalized_est_clonal_exp*log2(clone.frequency$normalized_est_clonal_exp))
  clonality <- 1 - entropy/log2(dim(clone.frequency)[1])
  res <- c(sample = sampleID,clonality = clonality)
  print(sampleID)
  return(res)
}

##function of SHM ratio
getSHMratio <-function(tmp){
  all.dist <- 0
  mm <- 0
  for (kk in tmp){
    no.gap <- sapply(kk$Sequences,function(x) nchar(gsub("-","",x)))
    base.l <- max(no.gap)                
    dist=kk$distMat
    dist[dist > 1] <- 0
    dist.s <- sum(dist)
    all.dist <- all.dist + dist.s
    mm <- mm + base.l
    print(base.l);print(dist.s)
  }
  s.ratio <- all.dist/mm
  return (s.ratio)
}

##function of computing TCR clonality
getClonalityTCR <- function(sampleID, Tdata = TCRdata) {
  ## Given sample ID, start and end position of complete CDR3,  return all the lineages in the sample
  #   sampleID <- "FZ-100Aligned.sortedByCoord.out.bam"
  #   Bdata = BCRdata
  ###get clones and frequency
  clone.frequency <- subset(Tdata, sample == sampleID ,select = c(sample,CDR3aa,frequency)) %>% 
    mutate(normalized_est_clonal_exp = frequency/sum(frequency) )
  ###calculate clonality
  entropy <- -sum(clone.frequency$normalized_est_clonal_exp*log2(clone.frequency$normalized_est_clonal_exp))
  clonality <- 1 - entropy/log2(dim(clone.frequency)[1])
  res <- c(sample = sampleID,clonality = clonality)
  return(res)
}

### BCR cluster & isotype class switch in each sample
get.bcr.cluster.classswitch <- function(bcr_clusters){
  bcr.cluster.isotypes <- NULL
  for(i in 1:length(bcr_clusters)){
    tmp <- bcr_clusters[[i]]
    if(is.null(tmp)==T){
      next
    }
    if(is.na(tmp)==T){
      next
    }
    tmp.cw <- matrix(0, nrow=length(tmp), ncol=12)
    colnames(tmp.cw) <- c('filename', 'motif', 'IGHA1','IGHA2','IGHE','IGHD','IGHM','IGHG1','IGHG2','IGHG3','IGHG4', 'Unidentified')
    tmp.cw[,'filename'] <- names(bcr_clusters)[i]
    tmp.cw[,'motif'] <- names(tmp)
    for(j in 1:length(tmp)){
      # j <- 1
      tmp_cluster <- tmp[[j]]
      tmp.is <- as.character(tmp_cluster$data$C)
      id <- which(tmp.is %in% c('IGHA1','IGHA2','IGHE','IGHD','IGHM','IGHG1','IGHG2','IGHG3','IGHG4'))
      if(length(id) > 0){
        tmp.is[-id] = 'Unidentified'
      }else{
        tmp.is = rep('Unidentified', length(tmp.is))
      }
      isotype.count <- table(tmp.is)
      tmp.cw[j, names(isotype.count)]=isotype.count
    }
    bcr.cluster.isotypes <- rbind(bcr.cluster.isotypes, tmp.cw)
  }
  return(bcr.cluster.isotypes)
}


######### function of computing Jacard Similarity
getCDR3Jacard <- function(meta, cdr3.complete){
  CDR3Jacard <- lapply(rownames(meta), function(x){
    n.bcr <- subset(cdr3.complete,sample == x)
    TBCR <- lapply(rownames(meta), function(y){
      if(x != y){
        t.bcr <- subset(cdr3.complete, sample == y)
        if(dim(t.bcr)[1] == 0){
          nt <- 0
        }
        else{
          nt <- signif(length(intersect(n.bcr$CDR3aa,t.bcr$CDR3aa))/length(na.omit(unique(t.bcr$CDR3aa))),5)
        }
        share <- cbind.data.frame(s1 = x, s2 = y, jacard = nt, 
                                  tIg = paste0(intersect(n.bcr$C,t.bcr$C),collapse = ","),
                                  nIg = paste0(intersect(n.bcr$C,t.bcr$C),collapse = ","),
                                  share = paste0(intersect(n.bcr$CDR3aa,t.bcr$CDR3aa),collapse = ","))
        return (share)
      }
    })
    similarity <- do.call("rbind",TBCR)
    return (similarity)
  })
  share.mat <- do.call("rbind",CDR3Jacard)
  return(share.mat)
}
