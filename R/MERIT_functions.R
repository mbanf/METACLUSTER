


#' install dependencies
#'
#'
#' @keywords 
#' @export
#' @examples
install_and_load_libraries <- function(){
  
  # CRAN
  
  list.of.packages <- c("ggplot2", "reshape2","doParallel", "pheatmap", "igraph", "seqinr", "networkD3", "taRifx", "parmigene", "ranger", "yaml", "plotROC")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  # bioconductor
  source("https://bioconductor.org/biocLite.R")
  
  list.of.packages <- c("Biostrings", "TFBSTools","seqLogo", "PWMEnrich", "BSgenome.Athaliana.TAIR.TAIR9")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) biocLite(new.packages)
  
  
  library(TFBSTools)
  # library(BSgenome.Athaliana.TAIR.TAIR9)
  
  library(networkD3)
  library(reshape2)
  library(foreach)
  library(doParallel)
  library(igraph)
  library(seqinr)
  library(Biostrings)
  library(PWMEnrich)
  library(seqLogo)
  
  library(ggplot2) 
  library(pheatmap)
  library(parmigene)
  library(ranger)
  
  library(taRifx)
  library(yaml)
  
  library(plotROC)
  
  
  
}


#' stabilitySelectionProcedure
#'
#'
#' @keywords 
#' @export
#' @examples
stabilitySelection <- function(x,y,nbootstrap=500,nsteps=5,alpha=0.2,plotme=FALSE){
  
  dimx <- dim(x)
  n <- dimx[1]
  p <- dimx[2]
  halfsize <- as.integer(n/2)
  freq <- matrix(0,nsteps+1,p)
  
  for (i in seq(nbootstrap)) {
    
    # Randomly reweight each variable
    xs <- t(t(x)*runif(p,alpha,1))
    
    # Ramdomly split the sample in two sets
    perm <- sample(dimx[1])
    i1 <- perm[1:halfsize]
    i2 <- perm[(halfsize+1):n]
    
    # run the randomized lasso on each sample and check which variables are selected
    r <- lars(xs[i1,],y[i1], max.steps=nsteps,normalize=FALSE, use.Gram=FALSE) #, eps = .Machine$double.eps)
    freq <- freq + abs(sign(coef.lars(r)))
    r <- lars(xs[i2,],y[i2],max.steps=nsteps,normalize=FALSE, use.Gram=FALSE)
    freq <- freq + abs(sign(coef.lars(r)))
    
  }
  
  # normalize frequence in [0,1]
  freq <- freq/(2*nbootstrap)
  
  if (plotme) {
    matplot(freq,type='l',xlab="LARS iteration",ylab="Frequency")
  }
  
  # the final stability score is the maximum frequency over the steps
  result <- apply(freq,2,max)
}

run_comparative_evaluation <- function(){
  
  n.min_hit_links = 5
  
  #  df.transcriptionFactorAnnotation=l.data$df.transcriptionFactorAnnotation
  #  df.geneGroups=l.data$df.geneGroups
  #  tb.geneGroups=l.data$tb.geneGroups
  #  v.geneGroups=l.data$v.geneGroups
  #  l.geneGroups=l.data$l.geneGroups
  #                                        
  # m.rf_grn = l.res.grn$m.rf_grn
  # m.lr_grn = l.res.grn$m.lr_grn
  # m.clr_grn = l.res.grn$m.clr_grn
  # m.MERIT_grn = l.res.link_annotation$m.grn * l.res.grn_tfbs$m.lead_suppport_w_motif.grn
  # 
  # 
 
  
  df.transcriptionFactorAnnotation = subset(df.transcriptionFactorAnnotation, df.transcriptionFactorAnnotation$with_geneExpression == "yes")
  df.geneGroups = subset(df.geneGroups, df.geneGroups$with_geneExpression == "yes")
  v.tfs = unique(df.transcriptionFactorAnnotation$TF_ID)
  v.genes =  unique(c(v.tfs, rownames(df.geneGroups)))
  

  ####  
  
  df.ATRM <- read.table("A:/junkDNA.ai/MERIT/datasets/evaluation/Regulations_in_ATRM.txt", header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)[,1:2]
  names(df.ATRM) <- c("TF", "Target")
  
  df.AGRIS <- read.table("A:/junkDNA.ai/MERIT/datasets/evaluation/reg_net.txt", header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")
  df.AGRIS$V2 <- toupper(df.AGRIS$V2)
  df.AGRIS$V5 <- toupper(df.AGRIS$V5)
  df.AGRIS <- subset(df.AGRIS, df.AGRIS$V8 == "Confirmed")
  df.AGRIS <- df.AGRIS[,c(2,5)]
  names(df.AGRIS) <- c("TF", "Target")
  #   tb.tfs <- table(df.AGRIS$TF)
  #   tb.tfs <- tb.tfs[tb.tfs <= 1000] # test this
  #   df.AGRIS <- subset(df.AGRIS, df.AGRIS$TF %in% names(tb.tfs))
  # goldstandard - dataframe with two columns, labelled - TF, Target
  
  df.AGRIS <- subset(df.AGRIS, df.AGRIS$TF != df.AGRIS$Target)
  m.AGRIS <- acast(df.AGRIS, TF~Target)
  m.AGRIS[is.na(m.AGRIS)] <- 0
  m.AGRIS[m.AGRIS != 0] <- 1
  class(m.AGRIS) <- "numeric"
  
  
  df.ATRM <- subset(df.ATRM, df.ATRM$TF != df.ATRM$Target)
  m.ATRM <- acast(df.ATRM, TF~Target)
  m.ATRM[is.na(m.ATRM)] <- 0
  m.ATRM[m.ATRM != 0] <- 1
  class(m.ATRM) <- "numeric"
  
  
  df.AGRIS_ATRM <- unique(rbind(df.AGRIS, df.ATRM))
  # removing self regulation
  df.AGRIS_ATRM <- subset(df.AGRIS_ATRM, df.AGRIS_ATRM$TF != df.AGRIS_ATRM$Target)
  m.AGRIS_ATRM <- acast(df.AGRIS_ATRM, TF~Target)
  m.AGRIS_ATRM[is.na(m.AGRIS_ATRM)] <- 0
  m.AGRIS_ATRM[m.AGRIS_ATRM != 0] <- 1
  class(m.AGRIS_ATRM) <- "numeric"
  
  df.Chip500 = read.table("A:/junkDNA.ai/MERIT/datasets/evaluation/Chip500.txt", sep = "\t", stringsAsFactors = F, header = T)
  m.Chip500 <- acast(df.Chip500, TF.Gene.ID~Target.Gene.ID)
  m.Chip500[is.na(m.Chip500)] <- 0
  m.Chip500[m.Chip500 != 0] <- 1
  class(m.Chip500) <- "numeric"  
  
  df.DE = read.table("A:/junkDNA.ai/MERIT/datasets/evaluation/DE.txt", sep = "\t", stringsAsFactors = F, header = T)
  m.DE <- acast(df.DE, TF.Gene.ID~Target.Gene.ID)
  m.DE[is.na(m.DE)] <- 0
  m.DE[m.DE != 0] <- 1
  class(m.DE) <- "numeric"  
  
  ### 
  
  # df.aranet <- read.table("datasets/evaluation/AraNet_V2.txt", header = FALSE, sep = "\t")
  # names(df.aranet) <- c("G1", "G2", "val")
  # m.aranet <- acast(df.aranet, G1~G2, value.var = "val")
  # m.aranet[is.na(m.aranet)] <- 0
  # m.aranet[m.aranet != 0] <- 1
  # class(m.aranet) <- "numeric"   
  
  df.dna = read.table("A:/junkDNA.ai/MERIT/datasets/dna_binding_motifs/df.grn.dapSeq.csv", header = T, stringsAsFactors = F, sep = ";")
  m.dna <- acast(df.dna, TF~Target)
  m.dna[is.na(m.dna)] <- 0
  m.dna[m.dna != 0] <- 1
  class(m.dna) <- "numeric"
  
  # df.Nature2015 = read.table("datasets/evaluation/Nature2015.txt", header = T, stringsAsFactors = F, sep = "\t", fill = T)
  # m.Nature2015 <- acast(df.Nature2015, Transcription.factor~Promoter)
  # m.Nature2015[is.na(m.Nature2015)] <- 0
  # m.Nature2015[m.Nature2015 != 0] <- 1
  # class(m.Nature2015) <- "numeric"
  
  
  df.NCA = read.csv("A:/junkDNA.ai/MERIT/datasets/evaluation/NCA.csv", stringsAsFactors = F, sep = ";")
  df.NCA <- subset(df.NCA, df.NCA$TFGenes %in% v.tfs & df.NCA$TG %in% v.genes)
  m.NCA <- acast(df.NCA, TFGenes~TG)
  m.NCA[is.na(m.NCA)] <- 0
  m.NCA[m.NCA != 0] <- 1
  class(m.NCA) <- "numeric"
  
  
  df.pcc_pre_NCA = read.csv("A:/junkDNA.ai/MERIT/datasets/evaluation/pcc_pre_NCA.csv", stringsAsFactors = F, sep = ",")[,1:3]
  df.pcc_pre_NCA <- subset(df.pcc_pre_NCA, df.pcc_pre_NCA$TFid %in% v.tfs & df.pcc_pre_NCA$TGid %in% v.genes)
  m.pcc_pre_NCA <- acast(df.pcc_pre_NCA, TFid~TGid, value.var = "PCC")
  m.pcc_pre_NCA[is.na(m.pcc_pre_NCA)] <- 0
  # m.pcc_pre_NC[m.pcc_pre_NC != 0] <- 1
  class(m.pcc_pre_NCA) <- "numeric"
  m.pcc_pre_NCA <- abs(m.pcc_pre_NCA)
  
  
  
  df.Vermeissen2014 = read.table("A:/junkDNA.ai/MERIT/datasets/evaluation/Vermeissen2014.txt", header = T, stringsAsFactors = F, sep = "\t")
  df.Vermeissen2014 <- subset(df.Vermeissen2014, df.Vermeissen2014$TF %in% v.tfs & df.Vermeissen2014$Target.gene %in% v.genes)
  m.Vermeissen2014 <- acast(df.Vermeissen2014, TF~Target.gene, value.var = "Rank")
  m.Vermeissen2014[is.na(m.Vermeissen2014)] <- 0
  # m.Vermeissen2014[m.Vermeissen2014 != 0] <- 1
  class(m.Vermeissen2014) <- "numeric"
  
  
  # df.GGM2013 = read.table("datasets/evaluation/GGM2013.txt", header = T, stringsAsFactors = F, sep = "\t")
  # m.GGM2013 <- acast(df.GGM2013, Gene_A_Name~Gene_B_Name)
  # m.GGM2013[is.na(m.GGM2013)] <- 0
  # m.GGM2013[m.GGM2013 != 0] <- 1
  # class(m.GGM2013) <- "numeric"
  
  # df.Carrera2009 = read.table("datasets/evaluation/Carrera2009.txt", header = F, stringsAsFactors = F, sep = "\t")
  # df.Carrera2009 <- subset(df.Carrera2009, df.Carrera2009$V1 %in% v.tfs.on_chip & df.Carrera2009$V2 %in% v.tgs.on_chip)
  # m.Carrera2009 <- acast(df.Carrera2009, V1~V2)
  # m.Carrera2009[is.na(m.Carrera2009)] <- 0
  # m.Carrera2009[m.Carrera2009 != 0] <- 1
  # class(m.Carrera2009) <- "numeric"
  
  df.CNS2014 = read.table("A:/junkDNA.ai/MERIT/datasets/evaluation/CNS2014.txt", header = T, stringsAsFactors = F, sep = "\t")[,1:2]
  df.CNS2014 <- subset(df.CNS2014, df.CNS2014$Transcription.factor %in% v.tfs & df.CNS2014$Target.gene %in% v.genes)
  m.CNS2014 <- acast(df.CNS2014, Transcription.factor~Target.gene)
  m.CNS2014[is.na(m.CNS2014)] <- 0
  m.CNS2014[m.CNS2014 != 0] <- 1
  class(m.CNS2014) <- "numeric"
  
  
  df.CNS_PCC2014 = read.table("A:/junkDNA.ai/MERIT/datasets/evaluation/CNS_PCC2014.txt", header = T, stringsAsFactors = F, sep = "\t")
  df.CNS_PCC2014 <- subset(df.CNS_PCC2014, df.CNS_PCC2014$Transcription.factor %in% v.tfs & df.CNS_PCC2014$Target.gene %in% v.genes)
  df.CNS_PCC2014.unique = unique(df.CNS_PCC2014[,1:2])
  df.CNS_PCC2014.unique["PCC"] = 0
  
  for(i in 1:nrow(df.CNS_PCC2014.unique)){
    
    TF.i = df.CNS_PCC2014.unique$Transcription.factor[i]
    TG.i = df.CNS_PCC2014.unique$Target.gene[i]
    
    df.CNS_PCC2014.i = subset(df.CNS_PCC2014, df.CNS_PCC2014$Transcription.factor %in% TF.i & df.CNS_PCC2014$Target.gene %in% TG.i)
    df.CNS_PCC2014.unique$PCC[i] = max(abs(df.CNS_PCC2014.i$PCC))
    
  }
  
  m.CNS_PCC2014 <- acast(df.CNS_PCC2014.unique, Transcription.factor~Target.gene, value.var = "PCC")
  m.CNS_PCC2014[is.na(m.CNS_PCC2014)] <- 0
  #m.CNS_PCC2014[m.CNS_PCC2014 != 0] <- 1
  class(m.CNS_PCC2014) <- "numeric"
  
  
  l.gn.comparison <- list(m.ATRM,
                          m.AGRIS,
                          m.dna,
                          m.Chip500,
                          m.DE)
  names(l.gn.comparison) <- c("ATRM", "AtRegNet", "DapSeq",  "TF2Network Benchmark 1", "TF2Network Benchmark 2")
  
  
  
  # subset to the genes in question 
  for(i in 1:length(l.gn.comparison)){
    tfs.i = rownames(l.gn.comparison[[i]])
    tgs.i = colnames(l.gn.comparison[[i]])
    l.gn.comparison[[i]] = l.gn.comparison[[i]][intersect(rownames(m.MERIT_grn), tfs.i), 
                                                intersect(colnames(m.MERIT_grn), tgs.i)]
  }
  
  
  v.methods <- c("CNS2014", "CNS_PCC2014", "Vermeissen2014",  "Barah 2015 PCC", "Barah 2015 NCA",
                 "Linear regression", "CLR", "Random Forest Regression", "MERIT")
  
  l.grn <- vector(mode = "list", length = length(v.methods))
  
  l.grn[[1]] <- m.CNS2014
  l.grn[[2]] <- m.CNS_PCC2014
  l.grn[[3]] <- m.Vermeissen2014
  l.grn[[4]] <- m.pcc_pre_NCA
  l.grn[[5]] <- m.NCA
  
  l.grn[[6]] = m.lr_grn
  l.grn[[7]] = m.clr_grn
  l.grn[[8]] = m.rf_grn
  l.grn[[9]] = m.MERIT_grn

  names(l.grn) = v.methods
  
  l.grn.selected = l.grn
  
  ####
  # 
  # l.grn <- vector(mode = "list", length = length(v.methods))
  # for(i in 1:length(l.grn.selected)){
  #   
  #   tgs = colnames(l.grn[[i]])
  #   tgs = tgs[!tgs %in% ""]
  #   
  #   l.grn[[i]] = l.grn[[i]]  # [, intersect(tg.selection, tgs)]
  #   l.grn.selected[[i]] <- l.grn[[i]] #[names(which(rowSums(l.grn[[i]]) > 0)), names(which(colSums(l.grn[[i]]) > 0))]#   * m.motifNet[names(which(rowSums(m.rf_w_treatments_w_motifs) > 0)), names(which(colSums(m.rf_w_treatments_w_motifs) > 0))]
  #   # l.grn.selected[[i]] <- l.grn[[i]][intersect(rownames(l.grn[[i]]), rownames(m.AGRIS_ATRM)),intersect(colnames(l.grn[[i]]),colnames(m.AGRIS_ATRM))]  #[names(which(rowSums(m.rf_w_treatments_w_motifs) > 0)), names(which(colSums(m.rf_w_treatments_w_motifs) > 0))]
  #   #l.grn.selected[[i]][l.grn.selected[[i]] > 0] <- 1
  # }
  # names(l.grn.selected) <- v.methods
  # 
  
  df.comparativeEvaluation.total <- c()
  for(d in 1:length(l.gn.comparison)){
    
    m.gn.comparison <- l.gn.comparison[[d]]
    m.gn.comparison[m.gn.comparison > 0] <- 1
    
    df.comparativeEvaluation <- data.frame(method = v.methods, 
                                           hits = rep(0,length(v.methods)),
                                           links = rep(0,length(v.methods)),
                                           foldchange = rep(0,length(v.methods)), 
                                           pval = rep(1,length(v.methods)),
                                           dataset = names(l.gn.comparison)[d]
    )
    
    
    for(i in 1:length(v.methods)){
      
      tfs <- intersect(rownames(l.grn.selected[[i]]), rownames(m.gn.comparison))
      tgs <- intersect(colnames(l.grn.selected[[i]]), colnames(m.gn.comparison))
      
      fc = -1
      p.val = -1
      hitInSample = -1
      sampleSize = -1
      
      if(length(tfs) >= 2 & length(tgs) >= 2){
        
        m.grn.select <- l.grn.selected[[i]][tfs, tgs]
        m.grn.select <- m.grn.select / max(m.grn.select)
        m.gn.comparison.select <- m.gn.comparison[tfs, tgs]
        
        m.grn.select[m.grn.select > 0] <- 1
        m.grn.select[m.grn.select <= 0] <- 0
        
        # optimize recovery or as in 
        n.links.total <- dim(m.grn.select)[1] * dim(m.grn.select)[2]
        fc <- (sum(m.grn.select * m.gn.comparison.select) / sum(m.grn.select)) / (sum(m.gn.comparison.select) / n.links.total)
        
        hitInSample <- sum(m.grn.select * m.gn.comparison.select) 
        sampleSize <- sum(m.grn.select)
        
        hitInPop <- sum(m.gn.comparison.select)
        failInPop <- (n.links.total - hitInPop)
        
        p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
        
        
      }
      ###
      
      df.comparativeEvaluation$foldchange[i] <- fc
      df.comparativeEvaluation$pval[i] <- p.val
      df.comparativeEvaluation$hits[i] <- hitInSample
      df.comparativeEvaluation$links[i] <- sampleSize
    }
    
    
    df.comparativeEvaluation.total <- rbind(df.comparativeEvaluation.total, df.comparativeEvaluation)
    
    
  }
  print(df.comparativeEvaluation.total)
  
  
  
  #n.min_hit_links = 5
  
  
  
  write.csv2(df.comparativeEvaluation.total, "Supp/df.GRN_Comparative_Evaluation_on_biological_datasets_subset_targets_and_tfs.csv",row.names = FALSE)
  
  
  
  
  
  
  
  ####
  
  # l.grn.selected <- vector(mode = "list", length = length(v.methods))
  # for(i in 1:length(l.grn.selected)){
  #   
  #   tgs = colnames(l.grn[[i]])
  #   tgs = tgs[!tgs %in% ""]
  #   
  #   l.grn[[i]] = l.grn[[i]]# [, intersect(tg.selection, tgs)]
  #   l.grn.selected[[i]] <- l.grn[[i]] #[names(which(rowSums(l.grn[[i]]) > 0)), names(which(colSums(l.grn[[i]]) > 0))]#   * m.motifNet[names(which(rowSums(m.rf_w_treatments_w_motifs) > 0)), names(which(colSums(m.rf_w_treatments_w_motifs) > 0))]
  #   # l.grn.selected[[i]] <- l.grn[[i]][intersect(rownames(l.grn[[i]]), rownames(m.AGRIS_ATRM)),intersect(colnames(l.grn[[i]]),colnames(m.AGRIS_ATRM))]  #[names(which(rowSums(m.rf_w_treatments_w_motifs) > 0)), names(which(colSums(m.rf_w_treatments_w_motifs) > 0))]
  #   #l.grn.selected[[i]][l.grn.selected[[i]] > 0] <- 1
  # }
  # names(l.grn.selected) <- v.methods
  # 
  # 
  
  library(plotROC)
  
  
  
  
  
  l.grn.selected <- l.grn.selected[c(2,3,4,6,7,8,9)]
  v.methods = names(l.grn.selected)
  
  dataset <- c()
  
  p.plots <- vector(mode = "list", length = length(l.gn.comparison))
  
  for(d in 1:length(l.gn.comparison)){
    
    m.gn.comparison <- l.gn.comparison[[d]]
    
    m.gn.comparison[m.gn.comparison > 0] <- 1
    df.idx.grn <- which(m.gn.comparison > 0, arr.ind = TRUE)
    df.gn.comparison <- data.frame(TF = rownames(m.gn.comparison)[df.idx.grn[,1]], TG = colnames(m.gn.comparison)[df.idx.grn[,2]], stringsAsFactors = FALSE)
    names(df.gn.comparison) <- c("TF","Target")
    
    n.selection = max(unlist(lapply(l.grn.selected, function(m) length(which(m > 0)))))
    
    l.predictions <- vector(mode = "list", length = length(l.grn.selected))
    l.performance <- vector(mode = "list", length = length(l.grn.selected))
    l.auc <- vector(mode = "list", length = length(l.grn.selected))
    
    l.D <- vector(mode = "list", length = length(l.grn.selected))
    l.M <- vector(mode = "list", length = length(l.grn.selected))
    
    for(i in 1:length(l.grn.selected)){      
      
      m.grn <- l.grn.selected[[i]]
      m.grn <- m.grn / max(m.grn) # map into range 0 - 1
      
      df.grn <- melt(m.grn, stringsAsFactors = FALSE)
      names(df.grn) <- c("TF", "Target","val")
      
      df.grn <- subset(df.grn, df.grn$val > 0)
      df.grn <- subset(df.grn, as.character(df.grn$TF) != as.character(df.grn$Target))
      df.grn <- df.grn[order(-df.grn$val),] # order by rank (highest first)
      
      n.links <- min(nrow(df.grn), n.selection) # select to top n hits (as predicted by MERIT)
      df.grn <- df.grn[1:n.links,]# nrow(df.regulatoryNetwork),]
      
      rownames(df.grn) <- seq(1:nrow(df.grn))
      
      df.grn.gs <- rbind(df.grn[,1:2], df.gn.comparison)
      i.set <- which(duplicated(df.grn.gs) == TRUE) # identify hits with GS
      
      n.add_links <- n.selection - n.links
      
      df.grn["class"] <- 0
      df.grn[df.grn$TF %in% df.grn.gs[i.set,1] & df.grn$Target %in% df.grn.gs[i.set,2], ]$class <- 1
      
      l.D[[i]] <- df.grn$class
      l.M[[i]] <- df.grn$val
      
    }
    
    rm(m.gn.comparison)
    
    
    
    D.rand <- rbinom(length(l.D[[length(v.methods)]]), 1, .5)
    M.rand <- rep(0, length(l.D[[length(v.methods)]]))
    Z.rand <- rep("random", length( l.D[[length(v.methods)]]))
    
    n_hits = numeric(length(v.methods))
    for(i in 1:length(v.methods)){
      n_hits[i] = sum(l.D[[i]])
    }
    idx = which(n_hits >= n.min_hit_links)
    
    n_hits.selection = n_hits[idx]
    v.methods.selection = v.methods[idx]
    v.colors <- c( "green", "blue", "pink", "yellow", "gray", "orange", "darkred")
    names(v.colors) <- v.methods
    
    v.colors.selection = c("black", v.colors[idx])
    names(v.colors.selection) = c("random", v.methods.selection)
    
    D <- D.rand
    M <- M.rand
    Z <- Z.rand
    
    for(i in 1:length(v.methods.selection)){
      D <- c(D, l.D[[i]])
      M <- c(M, l.M[[i]])
      Z <- c(Z, rep(v.methods.selection[i], length(l.D[[i]])))
      
      n_hits[i] = sum(l.D[[i]])
    }
    
    rocdata <- data.frame(D = D,
                          M = M,
                          method = Z)
    rocdata$method <- factor(rocdata$method, levels = c("random", v.methods.selection))
    rocdata["colors"] <- v.colors.selection[rocdata$method]
    
    ggroc5 <- ggplot(rocdata, aes(m = M, d = D, color = method)) + geom_roc(n.cuts = 0) + theme_bw() + scale_color_manual(values=v.colors.selection) + scale_size_manual(values = rep(0.5,6))  + scale_linetype_manual(values = c("dashed", rep("solid", 5))) +  ggtitle(names(l.gn.comparison)[d]) + labs(x="false positive rate",y="true positive rate") +
      theme(plot.margin=unit(c(0.5,0.5,1.5,1),"cm")) #+ geom_rocci()
    
    
    
    p.plots[[d]] <- ggroc5
    
    plot(ggroc5)
    print(calc_auc(ggroc5))    
    
    
    
    dataset = rbind(dataset , data.frame(dataset = rep(names(l.gn.comparison)[d], length(v.methods.selection) + 1 ), method = c("Random", v.methods.selection), hits = c(0, n_hits.selection), auroc  = calc_auc(ggroc5)[,3]))
    
    # pdf("figures_and_data_output/auroc_merit.pdf", 10, 4)
    # ggroc5
    # dev.off()
    
    # calc_auc(ggroc5)
    
  }
  
  p <- multiplot(p.plots[[1]], p.plots[[2]], p.plots[[3]], p.plots[[4]],p.plots[[5]], cols=3)
  filename <- paste("Supp//GRN_inference_performance_comparision_MIN_LINKS.pdf", sep = "")
  
  pdf(filename, width=20, height=10)
  p
  dev.off()  
  
  
  filename <- paste("Supp//df.GRN_inference_performance_comparision_AUROC_MIN_LINKS.csv", sep = "")
  write.csv2(dataset, filename, row.names = FALSE)
  
  
  
  # 
  
  # library
  
  
  df.auroc = dataset
  idx = which(df.auroc$hits <= n.min_hit_links & df.auroc$method != "Random")
  df.auroc$auroc[idx] = 0
  
  # create a dataset
  
  df.auroc$method <- factor(df.auroc$method, levels = c("Random", v.methods))
  
  # Grouped
  v.colors <- c("black", "cyan", "green", "blue", "pink", "yellow", "gray", "purple", "orange", "darkred")
  names(v.colors) <- c("random", v.methods)
  plot_auroc = ggplot(df.auroc, aes(fill=method, y=auroc, x=dataset)) + geom_bar(position="dodge", stat="identity") + theme_bw() + scale_fill_manual(values=v.colors)
  
  
  return(list(plot_auroc=plot_auroc, df.auroc=df.auroc, plot_auroc_curves = p, df.comparativeEvaluation.total=df.comparativeEvaluation.total ))
  
  # 
  # , fill = v.colors
  # 
  # dataset = read.csv2(filename, stringsAsFactors = F)
  # 
  # message("AUPR curve ")
  # 
  # 
  # 
  # df.comparativeEvaluation.total["AUROC"] = "-"
  # for(i in 1:nrow(df.comparativeEvaluation.total)){
  #   
  #   auroc = subset(dataset, dataset == as.character(df.comparativeEvaluation.total$dataset[i]) & 
  #                  method == as.character(df.comparativeEvaluation.total$method[i]))
  #   
  #   if(nrow(auroc) > 0)
  #     df.comparativeEvaluation.total$AUROC[i] = auroc$auroc[1]
  #   
  # }
  # 
  # df.comparativeEvaluation.total = df.comparativeEvaluation.total[,c(1,6,2,3,4,5,7)]
  # df.comparativeEvaluation.total$hits = as.integer(df.comparativeEvaluation.total$hits)
  # df.comparativeEvaluation.total$links = as.integer(df.comparativeEvaluation.total$links)
  # df.comparativeEvaluation.total$pval =  as.character(df.comparativeEvaluation.total$pval)
  # install.packages("xtable")
  #library(xtable)
  # return(list(auroc=calc_auc(ggroc5), p = ggroc5, df.gs_enrichment = df.gs_enrichment))
  # 
  
  #xtable(df.comparativeEvaluation.total)
  
  rm(l.grn)
  rm(l.grn.selected)
  rm(l.gn.comparison)
  
  
  
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

grn_dataframe_to_matrix <- function(df, val = "val"){
  m <- acast(df, TF~Target, value.var = val)
  m[is.na(m)] <- 0
  class(m) <- "numeric"
  return(m)
}

matrix_to_dataframe <- function(m){
  df <- melt(m, stringsAsFactors = FALSE)
  return(df)
}

prepare = function(){
  
  # transcription factor annotations
  
  tf.families = read.table("datasets/gene_annotation/Thale_cress-transcription factor.txt", header = FALSE, fill = TRUE, stringsAsFactors = FALSE)[,1:2] 
  tf.families$V1 <- gsub("\\..*", "",tf.families$V1)
  tf.families <- unique(tf.families)
  
  v.tf_families <- tf.families$V2
  names(v.tf_families) <- tf.families$V1
  
  v.tf_families[v.tf_families == "AP2/ERF-AP2"] <- "AP2_ERF-AP2"
  v.tf_families[v.tf_families == "AP2/ERF-ERF"] <- "AP2_ERF-ERF"
  v.tf_families[v.tf_families == "AP2/ERF-RAV"] <- "AP2_ERF-RAV"
  
  # add non present in thale cress - manual
  v.tffams.add <- c("C2H2", "zf-HD", "other", "other", "zf-HD", "other", "zf-HD")
  names(v.tffams.add) <- c("AT4G26030", "AT5G08750", "AT5G59430", "AT3G46590", "AT4G38170", "AT3G21890", "AT3G42860")
  
  v.tf_families <- (c(v.tf_families, v.tffams.add))
  # v.fams <- unique(tb.fams)
  v.tfs <- unique(names(v.tf_families))
  
  
  write.table(df.transcriptionFactorAnnotation, "data/df.transcriptionFactorAnnotation.txt", row.names = F, sep = "\t")
  
  
  
  
  # metabolic doains and enzymes
  df.domains <- read.table(filename.metabolicDomain_annotation, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  v.domains <- colnames(df.domains[,5:20])
  v.domains <- gsub(".Metabolism", "", v.domains)
  v.domains <- gsub("\\.", " ", v.domains)
  v.domains <- v.domains[!v.domains %in% c("Macromolecules", "Unclassified")]
  v.domains <- tolower(v.domains)
  
  df.domains <- df.domains[,c(1,5,6,7,8,9,10,11,12,13,15,16,17,18,19)]
  names(df.domains) <- c("Gene_ID",v.domains)
  tb.domains <- colSums(df.domains[,2:15])
  v.ath_coding_genes <- unique(df.domains$Gene_ID)
  v.enz <- df.domains$Gene_ID[which(rowSums(df.domains[,2:15]) > 0)]
  
  df.domains = subset(df.domains, df.domains$Gene_ID %in% v.enz)
  
  
  write.table(df.domains, "data/df.enzymes_w_metabolic_domains.txt", row.names = F,  sep = "\t")
  
  
}
  

compute_randomforest_based_GRN <- function(mat.expression, k="sqrt", nb.trees=10000, set.regulators = NULL, set.genes = NULL, seed=1234, importance.measure = "impurity", n.cpus = 5){
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  mat.expression.norm <- t(mat.expression)
  mat.expression.norm <- apply(mat.expression.norm, 2, function(x) { (x - mean(x)) / sd(x) } )
  n.genes <- dim(mat.expression.norm)[2]
  genes<- colnames(mat.expression.norm)
  n.samples <- dim(mat.expression.norm)[1]
  
  if(is.null(set.genes)){
    n.genes <- n.genes
    genes <- genes
    
  }else{
    n.genes <- length(set.genes)
    genes <- set.genes
    
  }
  
  #mat.expression.norm <- mat.expression.norm[,genes]
  
  if (is.null(set.regulators)) {
    n.regulators <- n.genes
    regulators <- genes
  } else {
    n.regulators <- length(set.regulators)
    # regulators provided by names
    if(is.character(set.regulators)){
      regulators <- set.regulators
      # genes.undefined <- setdiff(regulators, genes)
      # if (length(genes.undefined) != 0) {
      #   stop(paste("Error: genes: ", paste(genes.undefined, collapse = ",")," not represented in gene expression matrix \n", sep=""))
      # }
      # regulators provided by indices
    }else if (is.numeric(set.regulators)) {
      regulators <- genes[set.regulators]
    }else{
      stop("Error: invalid regulator format")
    }
  }
  # set mtry
  if (class(k) == "numeric") {
    mtry <- K
  } else if (k == "sqrt") {
    mtry <- round(sqrt(n.regulators))
  } else if (k == "all") {
    mtry <- n.regulators-1
  } else {
    stop("Error: invalid parameter k, options: \"sqrt\", \"all\", or integer")
  }
  
  print(paste("Performing random-forest regression based gene regulatory network inference (# of decision trees per gene is: ", nb.trees, ", # of regulators per decision tree node: ", mtry, sep=""))
  mat.rapid <- matrix(0.0, nrow=n.regulators, ncol=n.genes, dimnames = list(regulators, genes))
  
  strt<-Sys.time()
  for(i in 1:n.genes){
    cat("Processing... ", round(i/n.genes * 100, digits = 2) , "%", "\r"); flush.console() 
    gn.i <- genes[i]
    # remove target gene from set of regulators
    regulators.i <- setdiff(regulators, gn.i)
    x <- mat.expression.norm[,regulators.i, drop=FALSE] # avoid coercion to numeric
    y <- mat.expression.norm[,gn.i]
    
    
    
    df.xy <- cbind(as.data.frame(x),y)
    rf.model <- ranger(y ~ ., data = df.xy, mtry=mtry, num.trees=nb.trees, importance = "impurity", num.threads = n.cpus)
    imp_scores <- rf.model$variable.importance  
    imp_scores.names <- names(imp_scores)
    mat.rapid[imp_scores.names, gn.i] <- as.numeric(imp_scores)
  }
  print(Sys.time()-strt)
  print("..finished.")
  
  return(mat.rapid / n.samples)
}       


compute_linearRegressionWithStabilitySelection_based_GRN <- function(mat.expression, set.regulators = NULL, set.genes = NULL, nbootstrap = 100, nstepsLARS = 5, n.cpus = 5){
  
  mat.expression.norm <- t(mat.expression)
  mat.expression.norm <- apply(mat.expression.norm, 2, function(x) { (x - mean(x)) / sd(x) } )
  
  n.genes <- dim(mat.expression.norm)[2]
  genes<- colnames(mat.expression.norm)
  n.samples <- dim(mat.expression.norm)[1]
  
  if(is.null(set.genes)){
    n.genes <- n.genes
    genes <- genes
    
  }else{
    n.genes <- length(set.genes)
    genes <- set.genes
    
  }
  
  #mat.expression.norm <- mat.expression.norm[,genes]
  
  if (is.null(set.regulators)) {
    n.regulators <- n.genes
    regulators <- genes
  } else {
    n.regulators <- length(set.regulators)
    # regulators provided by names
    if(is.character(set.regulators)){
      regulators <- set.regulators
      # genes.undefined <- setdiff(regulators, genes)
      # if (length(genes.undefined) != 0) {
      #   stop(paste("Error: genes: ", paste(genes.undefined, collapse = ",")," not represented in gene expression matrix \n", sep=""))
      # }
      # regulators provided by indices
    }else if (is.numeric(set.regulators)) {
      regulators <- genes[set.regulators]
    }else{
      stop("Error: invalid regulator format")
    }
  }
  
  print(paste("Performing linear regression based gene regulatory network inference (# of bootstraps: ", nbootstrap, ", # of lars steps: ", nstepsLARS, sep=""))
  mat.rapid <- matrix(0.0, nrow=n.regulators, ncol=n.genes, dimnames = list(regulators, genes))
  
  strt<-Sys.time()
  
  cl<-makeCluster(n.cpus)
  registerDoParallel(cl)
  l.res <- foreach(i = 1:n.genes, .packages=c("lars", "MERIT")) %dopar% {
  
    #for(i in 1:n.genes){
    #  cat("Processing... ", round(i/n.genes * 100, digits = 2) , "%", "\r"); flush.console() 
    gn.i <- genes[i]
    # remove target gene from set of regulators
    regulators.i <- setdiff(regulators, gn.i)
    x <- mat.expression.norm[,regulators.i, drop=FALSE] # avoid coercion to numeric
    y <- mat.expression.norm[,gn.i]
    
    imp_scores <- stabilitySelection(x,y,nbootstrap=nbootstrap,nsteps=nstepsLARS,plotme=FALSE)
    
    imp_scores  
  }
  stopCluster(cl)
  
  print(Sys.time()-strt)
  print("..finished.")
  
  for(i in 1:n.genes){
    gn.i <- genes[i]
    imp_scores <- l.res[[i]]
    imp_scores.names <- names(imp_scores)
    mat.rapid[imp_scores.names, gn.i] <- as.numeric(imp_scores)
  }

  return(mat.rapid / n.samples)
}       




#' Step 1 - Gene regulatory network inference using ensemble regression with Monte Carlo based threshold selection
#'
#' 
#' @param mat.expression
#' @param k (default k="sqrt")
#' @param nb.trees
#' @param set.regulators
#' @param set.genes
#' @param seed
#' @param importance.measure
#' @param n.cpus
#' @param mat.expression
#' @param set.regulators
#' @param set.genes
#' @param nbootstrap
#' @param nstepsLARS
#' @param n.cpus
#' @keywords 
#' @export
#' @examples
compute_ensemble_regression_with_montecarlo_based_stability_selection <- function(m.foldChange_differentialExpression,
                                                                                  df.transcriptionFactorAnnotation,
                                                                                  df.geneGroups,
                                                                                  seed=1234,
                                                                                  importance.measure="impurity",
                                                                                  n.trees=1000,
                                                                                  n.lead_method_expression_shuffling = 3,
                                                                                  n.bootstrap=100,
                                                                                  n.stepsLARS=5,
                                                                                  n.cpus=5){
          
            df.transcriptionFactorAnnotation = subset(df.transcriptionFactorAnnotation, df.transcriptionFactorAnnotation$with_geneExpression == "yes")
            df.geneGroups = subset(df.geneGroups, df.geneGroups$with_geneExpression == "yes")
            v.tfs = unique(df.transcriptionFactorAnnotation$TF_ID)
            v.genes =  unique(c(v.tfs, rownames(df.geneGroups)))
            
            strt<-Sys.time()
            X <- m.foldChange_differentialExpression[v.genes,]
            colnames(X) = as.character(seq(1:dim(X)[2]))
            
            message("running lead method (random forest regression) with Monte Carlo based threshold selection")  
            m.rf_grn <- compute_randomforest_based_GRN(mat.expression=X, k="sqrt", nb.trees=n.trees, set.regulators = v.tfs, set.genes = v.genes, seed=seed, importance.measure = importance.measure, n.cpus = n.cpus)
            saveRDS(m.rf_grn, paste("tmp/m.grn.RF.rds", sep = ""))
            
            set.seed(seed)
            message("running linear regression") 
            m.lr_grn <- compute_linearRegressionWithStabilitySelection_based_GRN(mat.expression=X, set.regulators = v.tfs, set.genes = v.genes, nbootstrap = n.bootstrap, nstepsLARS = n.stepsLARS, n.cpus = n.cpus)
            saveRDS(m.lr_grn, paste("tmp/m.grn.LR.rds", sep = ""))
            
            message("running context likelihood of relatedness (CLR)")
            m.MI = knnmi.all(X)
            m.CLR = parmigene::clr(m.MI)
            m.CLR = m.CLR[v.tfs, v.genes]
            saveRDS(m.CLR, paste("tmp/m.grn.CLR.rds"))
            print(Sys.time()-strt)
            # message("running Pearson's correlation (PCC)")
            # m.PCC = cor(t(X), method = "pearson")
            # m.PCC = m.PCC[v.tfs, v.genes]
            # saveRDS(m.PCC, paste("tmp/m.grn.PCC.rds"))
            
            for(i in 1:n.lead_method_expression_shuffling){
              
              message("running background monte carlo run", i)
              strt<-Sys.time()
              set.seed(seed + 25 * i)
              X.shuffled <- t(apply (X, 1,  function(m) sample(m, length(m))))
              
              message("running random forest regression")
              m.rf.shuffled <- compute_randomforest_based_GRN(mat.expression=X.shuffled, k="sqrt", nb.trees=n.trees, set.regulators = v.tfs, set.genes = v.genes, seed= (seed + 25 * i), importance.measure = importance.measure, n.cpus = n.cpus)
              saveRDS(m.rf.shuffled, paste("tmp/m.grn.RF_bg_", i, ".rds", sep = ""))
        
              message("running linear regression") 
              set.seed(seed + 25 * i)
              m.lr_grn <- compute_linearRegressionWithStabilitySelection_based_GRN(mat.expression=X.shuffled, set.regulators = v.tfs, set.genes = v.genes, nbootstrap = n.bootstrap, nstepsLARS = n.stepsLARS, n.cpus = n.cpus)
              saveRDS(m.lr_grn, paste("tmp/m.grn.LR_bg_", i, ".rds", sep = ""))
              
              message("running context likelihood of relatedness (CLR)")
              m.MI = knnmi.all(X.shuffled)
              m.CLR = parmigene::clr(m.MI)
              m.CLR = m.CLR[v.tfs, v.genes]
              saveRDS(m.CLR, paste("tmp/m.grn.CLR_bg_", i, ".rds", sep = ""))
              print(Sys.time()-strt)
            }
            
           
            
            
            # for(i in 1:n.lead_method_expression_shuffling){
            #   
            #   message("running background monte carlo run", i)
            #   set.seed(seed + 25 * i)
            #   X.shuffled <- t(apply (X, 1,  function(m) sample(m, length(m))))
            #   
            #   # message("running random forest regression") 
            #   # strt<-Sys.time()
            #   # m.rf.shuffled <- compute_randomforest_based_GRN(mat.expression=X.shuffled, k="sqrt", nb.trees=n.trees, set.regulators = v.tfs, set.genes = v.genes, seed= (seed + 25 * i), importance.measure = importance.measure, n.cpus = n.cpus)
            #   # saveRDS(m.rf.shuffled, paste("tmp/m.grn.RF_bg_", i, ".rds", sep = ""))
            #   # print(Sys.time()-strt)
            #   
            #   # message("running linear regression") 
            #   # strt<-Sys.time()
            #   # set.seed(seed + 25 * i)
            #   # m.lr_grn <- compute_linearRegressionWithStabilitySelection_based_GRN(mat.expression=X.shuffled, set.regulators = v.tfs, set.genes = v.genes, nbootstrap = n.bootstrap, nstepsLARS = n.stepsLARS, n.cpus = n.cpus)
            #   # saveRDS(m.lr_grn, paste("tmp/m.grn.LR_bg_", i, ".rds", sep = ""))
            #   # print(Sys.time()-strt)
            #   # 
            #   message("running context likelihood of relatedness (CLR)")
            #   m.MI = knnmi.all(X.shuffled)
            #   m.CLR = parmigene::clr(m.MI)
            #   m.CLR = m.CLR[v.tfs, v.genes]
            #   saveRDS(m.CLR, paste("tmp/m.grn.CLR_bg_", i, ".rds", sep = ""))
            #   
            #   message("running Pearson's correlation (PCC)")
            #   m.PCC = cor(t(X.shuffled), method = "pearson")
            #   m.PCC = m.PCC[v.tfs, v.genes]
            #   saveRDS(m.PCC, paste("tmp/m.grn.PCC_bg_", i, ".rds", sep = ""))
            # 
            # }
            # 
            
            
}





#' load grns
#'
#' This function loads a computed grns 
#' @param df.transcriptionFactorAnnotation =l.data$df.transcriptionFactorAnnotation,
#' @param df.geneGroups = v.genes, 
#' @param th.lead_grn_method = 0.95,
#' @param th.support_grn_methods = 0.95,
#' @param n.lead_method_expression_shuffling = 3
#' @param ngrnSupport = 1
#' @keywords 
#' @export
load_lead_support_grn <- function(df.transcriptionFactorAnnotation,
                                   df.geneGroups, 
                                   th.lead_grn_method = 0.95,
                                   n.lead_method_expression_shuffling = 3,
                                   th.support_grn_methods = 0.95,
                                   n.grnSupport = 1){
  
  
  df.transcriptionFactorAnnotation = subset(df.transcriptionFactorAnnotation, df.transcriptionFactorAnnotation$with_geneExpression == "yes")
  df.geneGroups = subset(df.geneGroups, df.geneGroups$with_geneExpression == "yes")
  v.tfs = unique(df.transcriptionFactorAnnotation$TF_ID)
  v.genes =  unique(c(v.tfs, rownames(df.geneGroups)))

  

  # random forest regression - lead method
  m.rf_grn <- readRDS(paste("tmp/m.grn.RF.rds", sep = ""))[v.tfs,v.genes]
  v.th.grns <- numeric(n.lead_method_expression_shuffling)
  quantile(m.rf_grn, th.lead_grn_method)
  for(i in 1:n.lead_method_expression_shuffling){
    m.rf.shuffled = readRDS(paste("tmp/m.grn.RF_bg_", i, ".rds", sep = ""))[v.tfs,v.genes]
    v.th.grns[i] = quantile(m.rf.shuffled, th.lead_grn_method)
  }
  #print(v.th.grns)
  th.grns = sum(v.th.grns) / n.lead_method_expression_shuffling
  m.rf_grn[m.rf_grn < th.grns] <- 0

  # support method
  m.lr_grn = readRDS(paste("tmp/m.grn.LR.rds", sep = ""))[v.tfs,v.genes]
  quantile(m.lr_grn, th.lead_grn_method)
  v.th.grns <- numeric(n.lead_method_expression_shuffling)
  for(i in 1:n.lead_method_expression_shuffling){
    m.grn.shuffled = readRDS(paste("tmp/m.grn.LR_bg_", i, ".rds", sep = ""))[v.tfs,v.genes]
    v.th.grns[i] = quantile(m.grn.shuffled, th.support_grn_methods)
  }
  th.grns = sum(v.th.grns) / n.lead_method_expression_shuffling
  #print(v.th.grns)
  m.lr_grn[m.lr_grn < th.grns] <- 0
  m.lr_grn.classification = m.lr_grn
  m.lr_grn.classification[m.lr_grn.classification >= th.grns] = 1
  
  m.clr_grn = readRDS(paste("tmp/m.grn.CLR.rds"))[v.tfs,v.genes]
  quantile(m.clr_grn, th.lead_grn_method)
  v.th.grns <- numeric(n.lead_method_expression_shuffling)
  for(i in 1:n.lead_method_expression_shuffling){
    m.grn.shuffled = readRDS(paste("tmp/m.grn.CLR_bg_", i, ".rds", sep = ""))[v.tfs,v.genes]
    v.th.grns[i] = quantile(m.grn.shuffled, th.support_grn_methods)
  }
  th.grns = sum(v.th.grns) / n.lead_method_expression_shuffling
  #print(v.th.grns)
  m.clr_grn[m.clr_grn < th.grns] <- 0
  m.clr_grn.classification = m.clr_grn
  m.clr_grn.classification[m.clr_grn.classification >= th.grns] = 1
  
  # m.pcc_grn = abs(readRDS(paste("tmp/m.grn.PCC.rds", sep = "")))[v.tfs,v.genes]
  # quantile(m.pcc_grn, th.lead_grn_method)
  # v.th.grns <- numeric(n.lead_method_expression_shuffling)
  # for(i in 1:n.lead_method_expression_shuffling){
  #   m.grn.shuffled = abs(readRDS(paste("tmp/m.grn.PCC_bg_", i, ".rds", sep = ""))[v.tfs,v.genes])
  #   v.th.grns[i] = quantile(m.grn.shuffled, th.support_grn_methods)
  # }
  # th.grns = sum(v.th.grns) / n.lead_method_expression_shuffling
  # print(v.th.grns)
  # m.pcc_grn[m.pcc_grn < th.grns] <- 0
  # m.PCC.classification = m.pcc_grn
  # m.PCC.classification[m.PCC.classification >= th.grns] = 1
  # 
  sum(m.lr_grn.classification)
  sum(m.clr_grn.classification)
  #sum(m.PCC.classification)
  
  
  
  #####

  m.grn.support = m.lr_grn.classification + m.clr_grn.classification#  + m.PCC.classification
  m.grn.support[m.grn.support <  n.grnSupport] = 0
  m.grn.support[m.grn.support >= n.grnSupport] = 1
  
  m.lead_support.grn = m.rf_grn * m.grn.support 
  
  
  return(list(m.lead_support.grn = m.lead_support.grn, 
              m.rf_grn = m.rf_grn, 
              m.lr_grn = m.lr_grn, 
              m.clr_grn = m.clr_grn))
  
}




#' Load dataset function
#'
#' This function loads a datasets
#' @param 
#' @keywords 
#' @export
#' @examples
#' load_datasets()
load_datasets = function(filename.genes = "data/genes.txt",
                         filename.experiment_series_ids = "data/experiment_series_ids.txt",
                         filename.foldChange_differentialExpression = "data/m.foldChange_differentialExpression.txt",
                         filename.pvalue_differentialExpression =	"data/m.pvalue_differentialExpression.txt",
                         filename.experiment_condition_tissue_annotation =	"data/df.experiment_condition_annotation.txt",
                         filename.transcriptionfactor_annotation = "data/df.transcriptionFactorAnnotation.txt", 
                         filename.geneGroups = "data/df.enzymes_w_metabolic_domains.txt"){
  
            
              genes = read.table(filename.genes, header = F, sep = "\t", stringsAsFactors = F)[,1]
              experiment_series_ids = read.table(filename.experiment_series_ids, header = F, sep = "\t", stringsAsFactors = F)[,1]
              df.annotation = read.table(filename.experiment_condition_tissue_annotation, sep = "\t", stringsAsFactors = F, header = T)
              df.foldChange_differentialExpression = read.table(filename.foldChange_differentialExpression, header = F, sep = "\t", stringsAsFactors = F)
              df.pvalue_differentialExpression = read.table(filename.pvalue_differentialExpression, header = F, sep = "\t", stringsAsFactors = F)
     
              if(length(genes) == 0){
                stop("Error: no genes found")
              }
              if(length(experiment_series_ids) == 0){
                stop("Error: no experiments found")
              }
              if(nrow(df.annotation) == 0){
                stop("Error: no condition annotation found")
              }
              if(nrow(df.foldChange_differentialExpression) == 0){
                stop("Error: no differential expression foldchange found")
              }
              if(nrow(df.pvalue_differentialExpression) == 0){
                stop("Error: no differential expression pvalue found")
              }
          
              m.foldChange_differentialExpression = data.matrix(df.foldChange_differentialExpression, rownames.force = NA)
              m.pvalue_differentialExpression     = data.matrix(df.pvalue_differentialExpression, rownames.force = NA)
              rownames(m.foldChange_differentialExpression) = rownames(m.pvalue_differentialExpression) = genes
              colnames(m.foldChange_differentialExpression) = colnames(m.pvalue_differentialExpression) = experiment_series_ids
              
              v.tissues = unique(df.annotation$condition_tissue)
              v.treatments = unique(c(df.annotation$condition_treatment_1, df.annotation$condition_treatment_2))
              v.treatments = v.treatments[!v.treatments == ""]
              tb.experiment_series_ids = table(experiment_series_ids)
              
              df.annotation["number_series"] = 0
              for(i in 1:nrow(df.annotation)){
                df.annotation$number_series[i] = tb.experiment_series_ids[df.annotation$series_id[i]]
              }
              
              tb.tissues = numeric(length(v.tissues))
              names(tb.tissues) = v.tissues
              for(i in 1:length(v.tissues)){
                idx = which(df.annotation$condition_tissue == v.tissues[i])
                tissues.i = df.annotation$number_series[idx]
                tissues.i = tissues.i[!is.na(tissues.i)]
                tb.tissues[i] = sum(tissues.i)
              }
              tb.condition_tissues = tb.tissues
              
              
              tb.condition_treatments = numeric(length(v.treatments))
              names(tb.condition_treatments) = unique(v.treatments)
              for(i in 1:length(tb.condition_treatments)){
                idx_1 = which(df.annotation$condition_treatment_1 %in% names(tb.condition_treatments)[i])
                idx_2 = which(df.annotation$condition_treatment_2 %in% names(tb.condition_treatments)[i])
                tb.condition_treatments[i] = sum(df.annotation$number_series[idx_1]) + sum(df.annotation$number_series[idx_2])
              }
              
 
              df.transcriptionFactorAnnotation = read.table(filename.transcriptionfactor_annotation, header = T, sep = "\t", stringsAsFactors = F)
              if(nrow(df.transcriptionFactorAnnotation) == 0){
                stop("Error: no transcription factor annotation found")
              }
              
              df.transcriptionFactorAnnotation["with_geneExpression"] = "no"
              df.transcriptionFactorAnnotation$with_geneExpression[which(df.transcriptionFactorAnnotation$TF_ID %in% genes)] = "yes"
              
              
          
              df.geneGroups = read.table(filename.geneGroups, header = T,  sep = "\t", stringsAsFactors = F)
              if(nrow(df.geneGroups) == 0){
                stop("Error: no gene group annotation found")
              }
              
              rownames(df.geneGroups) = df.geneGroups$Gene_ID
              df.geneGroups <- df.geneGroups[,!names(df.geneGroups) %in% c("Gene_ID")]
              
              tb.geneGroups = colSums(df.geneGroups)
              v.geneGroups = colnames(df.geneGroups)
              
              l.geneGroups <- vector(mode = "list", length = length(v.geneGroups))
              names(l.geneGroups) <- v.geneGroups
              for(i in 1:length(v.geneGroups)){
                l.geneGroups[[i]] <- rownames(df.geneGroups)[which(df.geneGroups[,v.geneGroups[i]] == 1)]
                l.geneGroups[[i]] <- intersect(l.geneGroups[[i]], genes)
              }
              
              
              
              df.geneGroups["with_geneExpression"] = "no"
              df.geneGroups$with_geneExpression[which(rownames(df.geneGroups) %in% genes)] = "yes"
              
              
              return(list(m.foldChange_differentialExpression=m.foldChange_differentialExpression,
                          m.pvalue_differentialExpression=m.pvalue_differentialExpression,
                          df.experiment_condition_annotation=df.annotation,
                          tb.condition_treatments=tb.condition_treatments,
                          tb.condition_tissues=tb.condition_tissues,

                          df.transcriptionFactorAnnotation=df.transcriptionFactorAnnotation, 
                          df.geneGroups=df.geneGroups,
                          tb.geneGroups=tb.geneGroups,
                          v.geneGroups=v.geneGroups,
                          l.geneGroups=l.geneGroups
                          ))
                          
}
                        





#' Step 2 - Transcription factor direct target promoter binding based filtering of gene regulatory link predictions
#'
#' 
#' @param m.grn gene regulatory network 
#' @keywords 
#' @export
#' @examples
#' install_and_load_libraries()
transcriptionFactorBindingInference <- function(m.grn, 
                                                file.TF_to_Motif_IDs = "",
                                                file.TFBS_motifs = "",
                                                file.promoterSeq = "",
                                                file.geneSeq = "",
                                                th.pre_tss = 1000,
                                                th.post_tss = 200,
                                                th.min.score.motif = "80%",
                                                genome_nucleotide_distribution = c(0.3253439, 0.1746561, 0.1746561, 0.3253439 ),
                                                th.pval.known_motifs = 0.05,
                                                n.cpus = 2, 
                                                b.load = "no"){
  
  
  if(b.load == "no"){
  
    bg.genome = genome_nucleotide_distribution
    names(bg.genome) = c("A","C","G","T") 
  
    df.tfs_motifs <- read.table(file.TF_to_Motif_IDs, sep = "\t", header = T, stringsAsFactors = F)
    l.tfs_w_motifs = vector(mode = "list", length = nrow(df.tfs_motifs))
    names(l.tfs_w_motifs) = df.tfs_motifs$TF_ID
    for(i in 1:nrow(df.tfs_motifs)){
      l.tfs_w_motifs[[i]] = unlist(strsplit(df.tfs_motifs$Motif_ID[i], ","))
      if(length(l.tfs_w_motifs[[i]]) == 0){
        l.tfs_w_motifs[[i]] = NULL
      } 
    }
    
    v.tfs.w.open_motifs <- names(l.tfs_w_motifs)[(which(lapply(l.tfs_w_motifs, function(m) is.null(m)) == FALSE))]
    
    l.motifs_files <- read.fasta(file = file.TFBS_motifs, seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
    l.motifs <- vector(mode = "list", length = length(l.motifs_files))
    names(l.motifs) <- names(l.motifs_files)
    v.motif.ids <- names(l.motifs)
    for(m in 1:length(l.motifs_files)){
      
      motif <- l.motifs_files[[m]]
      motif <- unlist(strsplit(motif, "\t"))
      
      idx.c <- which(motif == "]c")
      idx.g <- which(motif == "]g")
      idx.t <- which(motif == "]t")
      
      v.a <- motif[2:(idx.c - 1)]
      v.c <- motif[(idx.c + 1):(idx.g - 1)]
      v.g <- motif[(idx.g + 1):(idx.t - 1)]
      v.t <- motif[(idx.t + 1):(length(motif) - 1)]
      
      v.a[1] <- gsub("\\[", "", v.a[1])
      v.c[1] <- gsub("\\[", "", v.c[1])
      v.g[1] <- gsub("\\[", "", v.g[1])
      v.t[1] <- gsub("\\[", "", v.t[1])
      
      m.motif <- matrix(0, nrow = 4, ncol = length(v.a), dimnames = list(c("A", "C", "G", "T"), seq(1:length(v.a))))
      m.motif[1,] <- as.numeric(v.a)
      m.motif[2,] <- as.numeric(v.c)
      m.motif[3,] <- as.numeric(v.g)
      m.motif[4,] <- as.numeric(v.t)
      
      l.motifs[[m]] <- m.motif
    }
    
    upstream_sequences <- read.fasta(file = file.promoterSeq, seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
    upstream_sequences <- lapply(upstream_sequences, function(m) {toupper(m)})
    df.promSequences <- DNAStringSet(unlist(upstream_sequences))
   
    gene_sequences <- read.fasta(file = file.geneSeq, seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
    gene_sequences <- lapply(gene_sequences, function(m) {toupper(m)})
    df.geneSequences <- DNAStringSet(unlist(gene_sequences))
    names(df.geneSequences) <- gsub("\\..*","", names(df.geneSequences)) # remove gene model 
    
    df.idx.grn <- which(m.grn > 0, arr.ind = TRUE)
    df.grn <- data.frame(TF = rownames(m.grn)[df.idx.grn[,1]], TG = colnames(m.grn)[df.idx.grn[,2]], stringsAsFactors = FALSE)
    tfs.grn <- unique(df.grn$TF)
    tgs.grn <- unique(df.grn$TG)
    tfs.grn <- intersect(tfs.grn, v.tfs.w.open_motifs)
    
    strt <- Sys.time()
    cl<-makeCluster(min(length(tfs.grn), n.cpus))
    registerDoParallel(cl)
    l.res <-  foreach(j = 1:length(tfs.grn), .packages = c("TFBSTools", "Biostrings")) %dopar% {
      
      v.pvals <- rep(1, length(tgs.grn))
      v.tgs <- rep(0,length(tgs.grn))
      names(v.pvals) <- names(v.tgs) <- tgs.grn
  
      tf.j <- tfs.grn[j]
      if(tf.j %in% names(l.tfs_w_motifs)){
        l.tfs_w_motifs.j <- l.tfs_w_motifs[tf.j]
        if(!is.null(unlist(l.tfs_w_motifs.j))){
          strt <- Sys.time()
          df.grn.r <- subset(df.grn, df.grn$TF == tf.j)
          if(nrow(df.grn.r) > 0){
            l.motifs.r <- l.motifs[unlist(l.tfs_w_motifs.j)]
            tgs.r <- df.grn.r$TG 
            tgs.r <- intersect(tgs.r, df.promSequences@ranges@NAMES)
            Sequences <- df.promSequences[tgs.r]
            
            l.pwms <- list()
            for(k in 1:length(l.motifs.r)){
              pwm.groundTruth <- l.motifs.r[[k]]
              pwm <- PWMatrix(ID = "", name = "", profileMatrix = pwm.groundTruth, bg = bg.genome)  
              l.pwms <- c(l.pwms, list(pwm))
            }
            pwmList <- do.call(PWMatrixList, l.pwms)
            
            #Sequences <- subseq(Sequences, 3001 - v.th.motif_binding_promoter[i], 3000) # pre startsite
            if(length(Sequences) > 0){
              sitesetList = searchSeq(pwmList, Sequences, min.score=th.min.score.motif, strand="*") #, mc.cores = n.cpus)
              l.pvalues_binding_per_sequence <- pvalues(sitesetList, type="sampling")
              pval.min.tgs <- unlist(lapply(l.pvalues_binding_per_sequence, function(m) min(unlist(m))))
              pval.min.tgs <- pval.min.tgs[pval.min.tgs != "Inf"]
              if(length(pval.min.tgs) > 0){
                tgs.pval_min <- unique(names(pval.min.tgs))
                pval.min <- numeric(length(tgs.pval_min))
                names(pval.min) <- tgs.pval_min
                for(k in 1:length(tgs.pval_min)){
                  idx.k <- which(names(pval.min.tgs) == tgs.pval_min[k])
                  pval.min[k] <- min(pval.min.tgs[idx.k])
                }
                pval.min.tgs <- pval.min
                tgs.bound <- names(pval.min.tgs) 
                v.pvals[tgs.bound] <- as.numeric(pval.min.tgs)
                v.tgs[tgs.bound] <- 1
              }
            }
            
            # post TSS 200 kb
            tgs.r <- intersect(tgs.r, df.geneSequences@ranges@NAMES)
            Sequences <- df.geneSequences[tgs.r]
            
            if(length(Sequences) > 0){
              
              for(k in 1:length(Sequences)){
                if(length(Sequences[[k]]) > th.post_tss){
                  Sequences[[k]] <-  subseq(Sequences[[k]], 1, th.post_tss)
                }
              }
              sitesetList = searchSeq(pwmList, Sequences, min.score=th.min.score.motif, strand="*") #, mc.cores = n.cpus)
              l.pvalues_binding_per_sequence <- pvalues(sitesetList, type="sampling")
              pval.min.tgs <- unlist(lapply(l.pvalues_binding_per_sequence, function(m) min(unlist(m))))
              pval.min.tgs <- pval.min.tgs[pval.min.tgs != "Inf"]
              
              if(length(pval.min.tgs) > 0){
                tgs.pval_min <- unique(names(pval.min.tgs))
                pval.min <- numeric(length(tgs.pval_min))
                names(pval.min) <- tgs.pval_min
                for(k in 1:length(tgs.pval_min)){
                  idx.k <- which(names(pval.min.tgs) == tgs.pval_min[k])
                  pval.min[k] <- min(pval.min.tgs[idx.k])
                }
                pval.min.tgs <- pval.min
                tgs.bound <- names(pval.min.tgs) #names(pval.min.tgs)[pval.min.tgs <= th.pval.known_motifs]
                
                v.pvals[tgs.bound] <- pmin(v.pvals[tgs.bound] , as.numeric(pval.min.tgs))
                v.tgs[tgs.bound] <- 1
                
              }
              
            }
          }
          print(Sys.time() - strt)
          
        }
      }
      
      l.tmp <- list(v.pvals=v.pvals, v.tgs = v.tgs)
      
      return(l.tmp)  
    }
    stopCluster(cl)
    print(Sys.time()-strt)
    
    #######
    
    m.motifNet.pval <- matrix(1, nrow = length(tfs.grn), ncol = length(tgs.grn), dimnames = list(tfs.grn, tgs.grn))
    m.motifNet.tgs <- matrix(0, nrow = length(tfs.grn), ncol = length(tgs.grn), dimnames = list(tfs.grn, tgs.grn))
    
    for(j in 1:length(tfs.grn)){
      m.motifNet.pval[tfs.grn[j],] <- l.res[[j]]$v.pvals
      m.motifNet.tgs[tfs.grn[j],] <- l.res[[j]]$v.tgs
    }
    
    saveRDS(m.motifNet.pval, "tmp/m.motifNet.pval.rds")
    saveRDS(m.motifNet.tgs, "tmp/m.motifNet.tgs.rds")
    
  }else{
    
    m.motifNet.pval = readRDS("tmp/m.motifNet.pval.rds")
    m.motifNet.tgs = readRDS("tmp/m.motifNet.tgs.rds")
    
    if(length(m.motifNet.pval) == 0){
      stop("Error: no m.motifNet.pval.rds")
    }
    
    if(length(m.motifNet.tgs) == 0){
      stop("Error: no m.motifNet.tgs.rds")
    }
    
  }
  
  m.motifNet = m.motifNet.pval
  m.motifNet[m.motifNet > th.pval.known_motifs] <- 10
  m.motifNet[m.motifNet <= th.pval.known_motifs] <- 1
  m.motifNet[m.motifNet == 10] <- 0
  
  tfs_w_motif_binding = intersect(rownames(m.motifNet), rownames(m.grn))
  tgs_w_motif_binding = intersect(colnames(m.motifNet), colnames(m.grn))
  m.lead_support_w_motif.grn <- m.motifNet[tfs_w_motif_binding, tgs_w_motif_binding] * m.grn[tfs_w_motif_binding, tgs_w_motif_binding]
  
  
  return(list(m.motifNet.pval=m.motifNet.pval, 
              m.motifNet.tgs=m.motifNet.tgs, 
              m.lead_support_w_motif.grn=m.lead_support_w_motif.grn))
}



do_treatment_filtering_all_links <- function(l.treatments_per_links){
  
  message("individual link based condition filtering based on hypergeometric test")
  
  pb <- txtProgressBar(min = 0, max = length(l.treatments_per_links), style = 3)
  for(l in 1:length(l.treatments_per_links)){
    setTxtProgressBar(pb, l)
    l.treatments.l <- l.treatments_per_links[[l]]
    
    for(j in 1:length(l.treatments.l)){ # root and shoot
      v.treatments.l <- l.treatments.l[[j]]
      p.prior <- l.treatments_per_links[[l]][[j]] / l.treatments.tissues[[j]][names(l.treatments_per_links[[l]][[j]])]
      if(length(v.treatments.l) > 0){
        p.treatment <- rep(1, length(v.treatments.l))
        names(p.treatment) <- names(v.treatments.l)
        for(i in 1:length(v.treatments.l)){
          hitInSample <- v.treatments.l[i] 
          sampleSize <- sum(v.treatments.l)
          hitInPop <- l.treatments.tissues[[j]][names(v.treatments.l)[i]]
          popSize <- sum(l.treatments.tissues[[j]])
          failInPop <- popSize - hitInPop
          fc <- ((hitInSample / sampleSize) / (hitInPop / popSize))
          if(fc > 1 & hitInSample >= th.min.samples)
            p.treatment[i] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
        }
        # p.treatment <- p.adjust(p.treatment,"bonferroni") 
        i.sets <- which(p.treatment <= th.prob)
        v.sets <- names(p.treatment)[i.sets]
        l.treatments.l[[j]] <- v.treatments.l[v.sets]
      }
    }
    l.treatments_per_links[[l]] <- l.treatments.l 
  }
  close(pb)
  
  return(l.treatments_per_links)
}


do_treatment_filtering_single_link <- function(tb.link_treatments, tb.condition_treatments, 
                                               th.pval.treatment = 0.05, th.min.samples = 1, 
                                               s.multipleTestCorrection = "none"){ 
  
  p.treatment <- rep(1, length(tb.link_treatments))
  names(p.treatment) <- names(tb.link_treatments)
  
  for(i in 1:length(tb.link_treatments)){
    
    hitInSample <- tb.link_treatments[i] 
    sampleSize <- sum(tb.link_treatments)
    hitInPop <- tb.condition_treatments[names(tb.link_treatments)[i]]
    popSize <- sum(tb.condition_treatments)
    failInPop <- popSize - hitInPop
    fc <- ((hitInSample / sampleSize) / (hitInPop / popSize))
    
    if(fc > 1 & hitInSample >= th.min.samples)
      p.treatment[i] <- phyper(hitInSample, hitInPop, failInPop, sampleSize, lower.tail = FALSE)
    
  }
  
  p.treatment <- p.adjust(p.treatment, s.multipleTestCorrection)
  i.sets <- which(p.treatment <= th.pval.treatment)
  tb.link_treatments <- tb.link_treatments[i.sets]
  
  return(tb.link_treatments)
}

evaluate_tissues_per_treatment <- function(tb.tissues_per_treatment_per_link,  tb.condition_tissues, th.pval.tissue = 0.05, th.min.samples = 1, s.multipleTestCorrection = "none"){ 
  
  p.treatment <- rep(1, length(tb.tissues_per_treatment_per_link))
  names(p.treatment) <- names(tb.tissues_per_treatment_per_link)
  
  for(i in 1:length(tb.tissues_per_treatment_per_link)){
    
    hitInSample <- tb.tissues_per_treatment_per_link[i] 
    sampleSize <- sum(tb.tissues_per_treatment_per_link)
    hitInPop <- tb.condition_tissues[names(tb.tissues_per_treatment_per_link)[i]]
    popSize <- sum(tb.condition_tissues)
    failInPop <- popSize - hitInPop
    
    fc <- ((hitInSample / sampleSize) / (hitInPop / popSize))
    
    if(fc > 1 & hitInSample >= th.min.samples)
      p.treatment[i] <- phyper(hitInSample, hitInPop, failInPop, sampleSize, lower.tail = FALSE)
    
  }
  
  p.treatment <- p.adjust(p.treatment, s.multipleTestCorrection)
  tb.tissues_per_treatment_per_link <- tb.tissues_per_treatment_per_link[which(p.treatment <= th.pval.tissue)]  
  return(tb.tissues_per_treatment_per_link)    
  
}




perform_treatment_and_tissue_filtering <- function(m.grn, 
                                                   m.de.bin, 
                                                   v.conditionGroups, 
                                                   v.tissueGroups,
                                                   df.annotation,
                                                   tb.condition_treatments,
                                                   tb.condition_tissues,
                                                   th.pval.treatment = 0.05, 
                                                   th.pval.tissue = 0.05,
                                                   th.min.samples = 1, 
                                                   s.multipleTestCorrection = "none"){  
  

  m.rf <- m.grn
  m.rf_w_treatments <- m.rf 
  v.tfs.rf <- rownames(m.rf)
  v.tgs.rf <- colnames(m.rf)
  
  idx.rf <- which(m.rf > 0, arr.ind = TRUE)
  
  # regulatory network 
  l.treatments_and_tissues <- vector(mode = "list", length = nrow(idx.rf))
  
  pb <- txtProgressBar(min = 0, max = nrow(idx.rf), style = 3)
  for(i in 1:nrow(idx.rf)){
    
    setTxtProgressBar(pb, i)
    
    tf <- v.tfs.rf[idx.rf[i, 1]]
    tg <- v.tgs.rf[idx.rf[i, 2]]
    
    sets.tf <- as.numeric(which(m.de.bin[tf,] == 1))
    sets.tg <- as.numeric(which(m.de.bin[tg,] == 1))
    
    sets.vals <- intersect(sets.tf, sets.tg)
    series.i <- colnames(m.de.bin)[sets.vals]
    
    tb.series.i <- table(series.i)
    series.i <- unique(series.i)
    
    df.annotation.i <- subset(df.annotation, df.annotation$series_id %in% series.i)
    
    ###
    
    tb.link_treatments <- numeric(length(v.conditionGroups))
    names(tb.link_treatments) = v.conditionGroups
    
    for(j in 1:length(series.i)){
      df.annotation.ij = subset(df.annotation.i, df.annotation.i$series_id == series.i[j])
      condition_treatment_1 = df.annotation.ij$condition_treatment_1
      condition_treatment_1 = condition_treatment_1[condition_treatment_1 != ""]
      if(length(condition_treatment_1) > 0){
        tb.link_treatments[v.conditionGroups[condition_treatment_1]] = tb.link_treatments[v.conditionGroups[condition_treatment_1]] + tb.series.i[series.i[j]]
      }
      condition_treatment_2 = df.annotation.ij$condition_treatment_2
      condition_treatment_2 = condition_treatment_2[condition_treatment_2 != ""]
      if(length(condition_treatment_2) > 0){
        tb.link_treatments[v.conditionGroups[condition_treatment_2]] = tb.link_treatments[v.conditionGroups[condition_treatment_2]] + tb.series.i[series.i[j]]
      }
    }
    tb.link_treatments = tb.link_treatments[tb.link_treatments > 0]
    
    if(length(tb.link_treatments) > 0){
      
      # identify significant treatments
      tb.link_treatments <- do_treatment_filtering_single_link(tb.link_treatments, 
                                                              tb.condition_treatments, 
                                                              th.pval.treatment, 
                                                              th.min.samples = th.min.samples, 
                                                              s.multipleTestCorrection = s.multipleTestCorrection)
          
      # remove all links without 
      if(length(tb.link_treatments) > 0){
        
        ## Step 2 - assign dominant tissues to these conditions 
        i.set <- which(v.conditionGroups[df.annotation.i$condition_treatment_1] %in% names(tb.link_treatments))
        i.set <- unique(c(i.set, which(v.conditionGroups[df.annotation.i$condition_treatment_2] %in% names(tb.link_treatments))))
        tb.tissues.i <- table(df.annotation.i$condition_tissue[i.set]) # tissue annotations all treatments 
        
        l.treatments <- vector(mode = "list", length = length(tb.link_treatments))
        names(l.treatments) <- names(tb.link_treatments)
        
        # dominant treatment filter per link
        for(t in 1:length(tb.link_treatments)){ 
          
          l.treatments[[t]] <- list()
          
          # which conditionset 
          i.set <- which(v.conditionGroups[df.annotation.i$condition_treatment_1] == names(tb.link_treatments)[t])
          i.set <- unique(c(i.set, which(v.conditionGroups[df.annotation.i$condition_treatment_2] == names(tb.link_treatments)[t])))
          df.annotation.it = df.annotation.i[i.set,]
          tissues.t = unique(df.annotation.it$condition_tissue)
          
          tb.tissues_link_treatments = numeric(length(tissues.t))
          names(tb.tissues_link_treatments) = tissues.t
          
          for(k in 1:length(tissues.t)){
            df.annotation.itk = subset(df.annotation.it, df.annotation.it$condition_tissue == tissues.t[k])
            tb.tissues_link_treatments[k] = sum(tb.series.i[df.annotation.itk$series_id])
          }
        
          tb.tissues_link_treatments <- evaluate_tissues_per_treatment(tb.tissues_link_treatments,  tb.condition_tissues, 
                                                           th.pval.tissue = th.pval.tissue, 
                                                           th.min.samples = th.min.samples, 
                                                           s.multipleTestCorrection = s.multipleTestCorrection)
          
          # if significant tissues for treatment are found
          if(length(tb.tissues_link_treatments) > 0){
            l.treatments[[t]] <- names(tb.tissues_link_treatments)
          }
          
        }
        
        idx.treatment <- which(sapply(l.treatments, function(m) length(m) > 0) == TRUE)
        # list with individual (multiple) tissue elements
        l.treatments <- l.treatments[idx.treatment]
        
        if(length(l.treatments) > 0){
          l.treatments_and_tissues[[i]] <- l.treatments
        }else{
          m.rf_w_treatments[tf,tg] <- 0
          l.treatments_and_tissues[[i]] <- list()
        }
      }else{
        m.rf_w_treatments[tf,tg] <- 0
        l.treatments_and_tissues[[i]] <- list()
      }
    }else{
      m.rf_w_treatments[tf,tg] <- 0
      l.treatments_and_tissues[[i]] <- list()
    }
  }
  close(pb)
  
  ###
  
  idx.treatment <- which(sapply(l.treatments_and_tissues, function(m) length(m) > 0) == TRUE)
  l.treatments_and_tissues <- l.treatments_and_tissues[idx.treatment]
  
  # list with individual (multiple) tissue elements
  # l.treatments <- l.treatments[idx.treatment]
  idx.rf <- idx.rf[idx.treatment,]

  # extract the link condition
  l.regulatoryNetwork_treatments_and_tissues <- vector(mode = "list", length= length(l.treatments_and_tissues))
  pb <- txtProgressBar(min = 0, max = length(l.treatments_and_tissues), style = 3)
  for(i in 1:length(l.treatments_and_tissues)){
    setTxtProgressBar(pb, i)
    idx_links.treatment_tissue_filter <- which(unlist(lapply(l.treatments_and_tissues[[i]], function(m) if(length(m) == 0){FALSE}else{TRUE})) == TRUE)
    l.regulatoryNetwork_treatments_and_tissues[[i]] <-  l.treatments_and_tissues[[i]][idx_links.treatment_tissue_filter]
  }
  close(pb)
  l.regulatoryNetwork_treatments_and_tissues <-  unlist(lapply(l.regulatoryNetwork_treatments_and_tissues, function(m) return(m)), recursive=FALSE)
  

  v.gns <- rownames(m.de.bin)
  v.condition_tissue_pairs <- sapply(unique(v.conditionGroups), function(m) paste(m, "-",  unique(v.tissueGroups)))
  m.gn_condition_tissue_differentialExpression <- matrix(0, nrow = length(v.gns), ncol = (length(v.condition_tissue_pairs)),dimnames = list(v.gns, v.condition_tissue_pairs))
  pb <- txtProgressBar(min = 0, max = length(v.gns), style = 3)
  for(i in 1:length(v.gns)){
    
    setTxtProgressBar(pb, i)
    
    tg <- v.gns[i]
    
    sets.tg <- as.numeric(which(m.de.bin[tg,] == 1))
    sets.vals <- unique(colnames(m.de.bin)[sets.tg])
    
    df.annotation.i <- subset(df.annotation, df.annotation$series_id %in% sets.vals)
    
    ## Step 1 - identify conditions (sublevel treatments) per link
    tb.treatments <- table(c(v.conditionGroups[df.annotation.i$condition_treatment_1], 
                             v.conditionGroups[df.annotation.i$condition_treatment_2]))
    tb.treatments <- tb.treatments[names(tb.treatments) != ""]
    
    if(length(tb.treatments) > 0){
      
      # identify significant treatments
      # tb.treatments <- do_treatment_filtering_single_link(tb.treatments, tb.treatments.total, th.pval.treatment, th.min.samples = th.min.samples)
      
      # # remove all links without 
      # if(length(tb.treatments) > 0){
      
      ## Step 2 - assign dominant tissues to these conditions 
      # tb.tissues.total <- table(df.annotation$tissue)
      
      i.set <- which(v.conditionGroups[df.annotation.i$condition_treatment_1] %in% names(tb.treatments))
      i.set <- unique(c(i.set, which(v.conditionGroups[df.annotation.i$condition_treatment_2] %in% names(tb.treatments))))
      
      # tb.tissues.j <- table(df.annotation.j$tissue)
      tb.tissues.i <- table(df.annotation.i$condition_tissue[i.set]) # tissue annotations all treatments 
      
      # dominant treatment filter per link
      for(t in 1:length(tb.treatments)){ 
        
        # which conditionset 
        i.set <- which(v.conditionGroups[df.annotation.i$condition_treatment_1] == names(tb.treatments)[t])
        i.set <- unique(c(i.set, which(v.conditionGroups[df.annotation.i$condition_treatment_2] == names(tb.treatments)[t])))
        
        tb.tissues.i.t <- table(df.annotation.i$condition_tissue[i.set])
        #tb.tissues.i.t <- evaluate_tissues_per_treatment(tb.tissues.i.t,  tb.tissues.total, th.pval.tissue = th.pval.tissue, th.min.samples = th.min.samples)
        
        # if significant tissues for treatment are found
        if(length(tb.tissues.i.t) > 0){
          v.cond_tis_pairs.i <- paste(names(tb.treatments)[t], "-",  names(tb.tissues.i.t))
          m.gn_condition_tissue_differentialExpression[tg , v.cond_tis_pairs.i] <- 1
        }
        
      }
    }
    # }
  }
  close(pb)
  
  
  
  # perform treatment filtering for the entire random forest set - integrate with motif late
  v.tfs.rf <- rownames(m.rf_w_treatments)
  v.tgs.rf <- colnames(m.rf_w_treatments)
  idx.rf <- which(m.rf_w_treatments > 0, arr.ind = TRUE)
  
  
  idx.treatment <- which(sapply(l.treatments_and_tissues, function(m) length(m) > 0) == TRUE)
  l.treatments_and_tissues <- l.treatments_and_tissues[idx.treatment]
  
  
  tb.condition_tissue_differentialExpression <- colSums(m.gn_condition_tissue_differentialExpression)
  #message(length(which(tb.condition_tissue_differentialExpression > 0)), " out of ", length(tb.condition_tissue_differentialExpression), " conditions and tissue pairs with expressed genes ")
  
  idx.pairs <- which(tb.condition_tissue_differentialExpression > 0)
  
  m.gn_condition_tissue_differentialExpression <- m.gn_condition_tissue_differentialExpression[ , idx.pairs]
  tb.condition_tissue_differentialExpression <- tb.condition_tissue_differentialExpression[idx.pairs]
  
  v.cond_tiss_pairs <- colnames(m.gn_condition_tissue_differentialExpression)
  
  l.grn_subnetworks <- vector(mode = "list", length = length(v.cond_tiss_pairs))
  names(l.grn_subnetworks) <- v.cond_tiss_pairs
  
  for(i in 1:length(l.grn_subnetworks)){
    l.grn_subnetworks[[i]] <- matrix(0, nrow = nrow(m.rf), ncol = ncol(m.rf), dimnames = list(rownames(m.rf), colnames(m.rf)))
  }
  
  pb <- txtProgressBar(min = 0, max = nrow(idx.rf), style = 3)
  for(i in 1:nrow(idx.rf)){
    
    setTxtProgressBar(pb, i)
    
    tf <- v.tfs.rf[idx.rf[i, 1]]
    tg <- v.tgs.rf[idx.rf[i, 2]]
    
    tb.treatments <- l.treatments_and_tissues[[i]]
    for(j in 1:length(tb.treatments)){
      v.ct_pair.ij <- paste(names(tb.treatments)[j], "-", tb.treatments[[j]]) 
      for(k in 1:length(v.ct_pair.ij)){
        idx <- which(names(l.grn_subnetworks) == v.ct_pair.ij[k])
        l.grn_subnetworks[[idx]][tf, tg] <- 1
      }
    }
  }
  close(pb)
  
  idx.subnetworks <- which( lapply(l.grn_subnetworks, sum) > 0)
  l.grn_subnetworks <- l.grn_subnetworks[idx.subnetworks]
  
  #message("number of inferred subnetworks ", length(idx.subnetworks))

  return(list(l.treatments_and_tissues=l.treatments_and_tissues, m.rf_w_treatments = m.rf_w_treatments, l.regulatoryNetwork_treatments_and_tissues = l.regulatoryNetwork_treatments_and_tissues, l.grn_subnetworks = l.grn_subnetworks,  m.gn_condition_tissue_differentialExpression = m.gn_condition_tissue_differentialExpression, tb.condition_tissue_differentialExpression = tb.condition_tissue_differentialExpression))
}








#' Step 3 - Context specific annotation and filtering of gene regulatory link predictions
#'
#' This function filter... 
#' @param m.grn = m.grn,
#' @param l.grn_subnetworks = l.grn_subnetworks, 
#' @param df.geneGroups,
#' @param tb.geneGroups,
#' @param v.geneGroups,
#' @param l.geneGroups,
#' @param th.min_number_targets = 2,
#' @param th.min_number_MR_targets = 2,
#' @param th.pval = 0.05
#' @keywords 
#' @export
#' @examples
annotate_links_with_treatments_and_tissues <- function(m.lead_support_w_motif.grn, 
                                                       m.pvalue_differentialExpression,
                                                       df.experiment_condition_annotation,
                                                       tb.condition_treatments,
                                                       tb.condition_tissues,
                                                       v.conditionGroups, 
                                                       v.tissueGroups,
                                                       th.diffexp = 0.05,
                                                       th.pval.treatment = 0.05, 
                                                       th.pval.tissue = 0.05,
                                                       th.min.samples = 1, 
                                                       s.multipleTestCorrection = "none",
                                                       b.load = "no"){
  
  
  v.conditionGroups=names(tb.condition_treatments)
  names(v.conditionGroups)=names(tb.condition_treatments)
  
  v.tissueGroups=names(tb.condition_tissues)
  names(v.tissueGroups)=names(tb.condition_tissues)
  
  subDir <- paste("tmp/", th.diffexp, sep ="")
  
  if(b.load == "no"){
    
    m.grn = m.lead_support_w_motif.grn
    m.grn[m.grn > 0] = 1
    m.de = m.pvalue_differentialExpression
    
    if (!file.exists(subDir)){
      dir.create(subDir)
    }
    
    m.de.bin <- m.de
    m.de.bin[m.de.bin <= th.diffexp] <- 10
    m.de.bin[m.de.bin <= 1] <- 0
    m.de.bin[m.de.bin > 0] <- 1
    
    strt<-Sys.time()
    l.res <- perform_treatment_and_tissue_filtering(m.grn = m.grn, 
                                                    m.de.bin = m.de.bin, 
                                                    v.conditionGroups = v.conditionGroups, 
                                                    v.tissueGroups = v.tissueGroups,
                                                    df.annotation = df.experiment_condition_annotation, 
                                                    tb.condition_treatments=tb.condition_treatments,
                                                    tb.condition_tissues=tb.condition_tissues,
                                                    th.pval.treatment = th.pval.treatment, 
                                                    th.pval.tissue = th.pval.tissue,
                                                    th.min.samples = th.min.samples, 
                                                    s.multipleTestCorrection = s.multipleTestCorrection)
    
    
    l.treatments_and_tissues <- l.res$l.treatments_and_tissues
    m.rf_w_treatments <- l.res$m.rf_w_treatments
    l.regulatoryNetwork_treatments_and_tissues <- l.res$l.regulatoryNetwork_treatments_and_tissues
    m.gn_condition_tissue_differentialExpression <- l.res$m.gn_condition_tissue_differentialExpression
    tb.condition_tissue_differentialExpression = l.res$tb.condition_tissue_differentialExpression
    l.grn_subnetworks <- l.res$l.grn_subnetworks
    
    filename <- paste(subDir,"/l.treatments_and_tissues.rds", sep ="")
    saveRDS(l.treatments_and_tissues, filename)
    
    filename <- paste(subDir,"/m.rf_w_treatments.rds", sep ="")
    saveRDS(m.rf_w_treatments, filename)
    
    filename <- paste(subDir,"/l.regulatoryNetwork_treatments_and_tissues.rds", sep ="")
    saveRDS(l.regulatoryNetwork_treatments_and_tissues, filename)
    
    filename <- paste(subDir,"/m.gn_condition_tissue_differentialExpression.rds", sep ="")
    saveRDS(m.gn_condition_tissue_differentialExpression, filename)
    
    filename <- paste(subDir,"/l.grn_subnetworks.rds", sep ="")
    saveRDS(l.grn_subnetworks, filename)
    
    filename <- paste(subDir,"/tb.condition_tissue_differentialExpression.rds", sep ="")
    saveRDS(tb.condition_tissue_differentialExpression, filename)
    
    print(Sys.time()-strt)
    
  }
  
  filename <- paste(subDir,"/m.rf_w_treatments.rds", sep ="")
  m.grn = readRDS(filename)
  filename <- paste(subDir,"/l.grn_subnetworks.rds", sep ="")
  l.grn_subnetworks = readRDS(filename)

  filename <- paste(subDir,"/tb.condition_tissue_differentialExpression.rds", sep ="")
  tb.condition_tissue_differentialExpression = readRDS(filename)
  
  df.grn_subnetworks_statistics =  data.frame(conditions = names(l.grn_subnetworks), n_links = unlist(lapply(l.grn_subnetworks, sum)))
  
  return(list(m.grn = m.grn, 
              l.grn_subnetworks=l.grn_subnetworks,
              df.grn_subnetworks_statistics = df.grn_subnetworks_statistics, 
              tb.condition_tissue_differentialExpression = tb.condition_tissue_differentialExpression))

  
}





# identify master regulators => save iniital for stability selection
identify_bottom_tier_masterRegulators <- function(m.grn,
                                                  l.grn_subnetworks,
                                                  v.tfs = rownames(m.grn),
                                                  v.conds = names(l.grn_subnetworks),
                                                  df.geneGroups,
                                                  tb.geneGroups,
                                                  v.geneGroups,
                                                  l.geneGroups,
                                                  th.min_number_targets = 2,
                                                  th.pval = 0.05,
                                                  b.include_under_represented = "yes"){
  
  
  
  
  
  v.genes = rownames(df.geneGroups)

  # metabolic enzymes
  m.tfs_vs_conditions <- matrix(NA, nrow = length(v.tfs), ncol = length(v.conds), 
                                dimnames = list(v.tfs, v.conds))
  
  m.tfs_vs_conditions.numbers <- matrix(NA, nrow = length(v.tfs), ncol = length(v.conds), 
                                        dimnames = list(v.tfs, v.conds))
  
  # metabolic domains 
  l.tfs_vs_domains_given_condition <- vector(mode = "list", length = length(v.conds))
  names(l.tfs_vs_domains_given_condition) <- v.conds
  
  l.tfs_vs_domains_given_condition.numbers <- vector(mode = "list", length = length(v.conds))
  names(l.tfs_vs_domains_given_condition.numbers) <- v.conds
  
  pb <- txtProgressBar(min = 0, max = length(v.conds), style = 3)
  for(i in 1:length(v.conds)){
    
    setTxtProgressBar(pb, i)
    
    ct.i <- v.conds[i]
    
    m.grn.i <- l.grn_subnetworks[[i]]
    df.idx.grn.i <- which(m.grn.i == 1, arr.ind = TRUE)
    df.grn.i <- data.frame(TF = rownames(m.grn.i)[df.idx.grn.i[,1]], TG = colnames(m.grn.i)[df.idx.grn.i[,2]], stringsAsFactors = FALSE)
    df.grn.i <- subset(df.grn.i, df.grn.i$TG %in% v.genes) # subset to metabolic enzumes
    
    
    l.tfs_vs_domains_given_condition[[i]] <- matrix(NA, nrow = length(v.tfs), ncol = length(v.geneGroups), 
                                                    dimnames = list(v.tfs, (v.geneGroups)))
    
    
    l.tfs_vs_domains_given_condition.numbers[[i]] <- matrix(NA, nrow = length(v.tfs), ncol = length(v.geneGroups), 
                                                            dimnames = list(v.tfs, (v.geneGroups)))
    
    if(nrow(df.grn.i) > 0){
      
      tfs.grn.i <- unique(df.grn.i$TF) # all transcritption factors in the network
      
      for(j in 1:length(tfs.grn.i)){
        
        tf.ij <- tfs.grn.i[j]
        
        # number of links per transcription factor family member in the subnetwork
        df.grn.ij <- subset(df.grn.i, df.grn.i$TF == tf.ij)  # number of targets in condition by TF
        
        # condition specific network 
        hitInSample = n_A_B = nrow(df.grn.ij) # links per transcriptoin factor (family, conditios) in the network
        sampleSize = n_A = nrow(df.grn.i) # n.links.i # number of links in the subnetwork
        
        # condition independent network
        hitInPop = n_B = sum(m.grn[tf.ij,]) #n.tgs.grn_w_treatments_w_motifs.r #sum(m.gn_condition_tissue_differentialExpression[tf.IDs.global,ct.i])
        popSize = n_C = sum(m.grn) #n.tgs.grn_w_treatments_w_motifs  # (n.grn_w_treatments_w_motifs) # total number of links
        failInPop = n_C - n_B
        
        if(hitInSample >= th.min_number_targets){
          
          pval <- 1
          fc <- (n_A_B / n_A) / (n_B / n_C)
          
          if(fc > 1){
            pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= FALSE)
          }else if(fc < 1 & b.include_under_represented == "yes"){
            pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= TRUE)
          }
          
          if(pval <= th.pval){
            m.tfs_vs_conditions[tf.ij,ct.i] <- fc # enriched or depleated over 
            m.tfs_vs_conditions.numbers[tf.ij,ct.i] = hitInSample
          }else{
            m.tfs_vs_conditions[tf.ij,ct.i] <- 1
            m.tfs_vs_conditions.numbers[tf.ij,ct.i] = 0
          }
          
          # B) per condition - TFs versus Domains (P(TF,D|C)) => also cumulative plot 
          for(d in 1:length(v.geneGroups)){ # add the domain level 
            
            ## domain genes in target genes - condition dependent 
            hitInSample = n_A_B = length(intersect((unique(df.grn.ij$TG)), l.geneGroups[[d]]))
            sampleSize = n_A = length(l.geneGroups[[d]])
            
            ## condition independent part 
            tgs.ij.condition_independent = names(which(m.grn[tf.ij,] > 0))
            hitInPop = n_B = length(intersect(tgs.ij.condition_independent, v.genes)) 
            #   hitInPop = n_B = sum(m.rf_w_treatments.stability_selection[tf.ij,]) #n.tgs.grn_w_treatments_w_motifs.r #sum(m.gn_condition_tissue_differentialExpression[tf.IDs.global,ct.i])
            # popSize = n_C = sum(m.rf_w_treatments.stability_selection) 
            #   
            # length(intersect((unique(df.grn.ij$TG)), v.enz)) 
            popSize = n_C = length(v.genes) 
            failInPop = n_C-n_B
            
            if(hitInSample >= th.min_number_targets){
              
              pval <- 1
              fc <- (n_A_B / n_A) / (n_B / n_C)
              
              if(fc > 1){
                pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= FALSE)
              }else if(fc < 1){
                pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= TRUE)
              }
              
              if(pval <= th.pval){
                l.tfs_vs_domains_given_condition[[i]][tf.ij,d] <- fc
                l.tfs_vs_domains_given_condition.numbers[[i]][tf.ij,d] <- fc
              }
              
            }
          }
          
        }
        
      }
    }
  }
  close(pb)
  
  
  return(list(m.MR_vs_conditions = m.tfs_vs_conditions, 
              l.MR_vs_geneGroups_given_condition = l.tfs_vs_domains_given_condition, 
              l.MR_vs_geneGroups_given_condition.numbers = l.tfs_vs_domains_given_condition.numbers))
  
}




# define the bottom layer
identify_regulatory_hierachy = function(m.MR_vs_conditions, 
                                        l.MR_vs_geneGroups_given_condition, 
                                        
                                        m.grn,
                                        l.grn_subnetworks,
                                        
                                        v.tfs = rownames(m.grn),
                                        v.conds = names(l.grn_subnetworks),
                                        
                                        df.geneGroups,
                                        tb.geneGroups,
                                        v.geneGroups,
                                        l.geneGroups,
                                        
                                        th.min_number_MR_targets = 2,
                                        mode = "geneGroups",
                                        th.pval = 0.05){ # or genes
  
  
 
  v.genes = rownames(df.geneGroups)
  conds = v.conds
  
  l.Hierarchy = vector(mode = "list", length = length(conds))
  names(l.Hierarchy) = conds
  
  v.number_tiers = numeric(length(conds))
  names(v.number_tiers) = conds
  
  l.Hierarchy_nb_tfs_per_tier = vector(mode = "list", length = length(conds))
  names(l.Hierarchy_nb_tfs_per_tier) = conds
  
  l.Hierarchy_tfs_per_tier = vector(mode = "list", length = length(conds))
  names(l.Hierarchy_tfs_per_tier) = conds
  
  for(i in 1:length(conds)){
    
    m.grn.i <- l.grn_subnetworks[[i]]
    
    df.idx.grn.i <- which(m.grn.i == 1, arr.ind = TRUE)
    df.grn.i <- data.frame(TF = rownames(m.grn.i)[df.idx.grn.i[,1]], TG = colnames(m.grn.i)[df.idx.grn.i[,2]], stringsAsFactors = FALSE)
    v.tfs.i = unique(df.grn.i$TF)
    v.tgs.i = unique(df.grn.i$TG)
    
    v.MRs = rownames(m.MR_vs_conditions)
    
    if(mode == "genes"){
      v.MR_level1 = rownames(m.MR_vs_conditions)[which(m.MR_vs_conditions[,i] > 0)]
      l.Hierarchy[[i]] = matrix(0, nrow = length(v.MRs), ncol = length(c(v.MRs, v.genes)), dimnames = list(v.MRs, c(v.MRs, v.genes)))
    }
    
    if(mode == "geneGroups"){
      m.tfs_vs_domains_given_condition = l.MR_vs_geneGroups_given_condition[[i]]
      m.tfs_vs_domains_given_condition[is.na(m.tfs_vs_domains_given_condition)] = 0
      m.tfs_vs_domains_given_condition[m.tfs_vs_domains_given_condition < 1] = 0 # only analyze enrichment
      v.MR_level1 = v.MR_level1_domains = names(which(rowSums(m.tfs_vs_domains_given_condition) > 0))
      l.Hierarchy[[i]] = matrix(0, nrow = length(v.MRs), ncol = length(c(v.MRs, v.geneGroups)), dimnames = list(v.MRs, c(v.MRs, v.geneGroups)))
    }
    
    
    v.tgs_level_prevs = v.MR_level1    
    v.tfs.putative = v.tfs.i
    
    tiers = 1
    v.tfs_known = c()
    
    l.Hierarchy_tfs_per_tier[[i]] = list(v.MR_level1)
    
    b.continue <- TRUE
    
    while(b.continue){
      
      v.tgs_level_next = c()
      
      for(j in 1:length(v.tfs.i)){
        
        tf.ij <- v.tfs.i[j]
        # number of links per transcription factor family member in the subnetwork
        df.grn.ij <- subset(df.grn.i, df.grn.i$TF == tf.ij)  # number of targets in condition by T
        tgs.ij <- unique(df.grn.ij$TG)
        tgs_MR.ij = intersect(tgs.ij, v.tgs_level_prevs) # 
        
        # number of level prev master regulator trgets over all level prev master regulator | condition
        hitInSample = n_A_B = length(tgs_MR.ij) # number of master regulators as targets of TF in condition
        sampleSize = n_A = length(v.tgs_level_prevs) # number of master regulators in condition (not just targets) - and level
        
        
        # number of condition and TF independent targets
        tgs.ij = names(which(m.grn[tf.ij,] > 0)) 
        
        # tgs_w_tfs.ij.condition_independent = intersect(v.tfs.w.open_motifs, tgs.ij) # number of 
        
        hitInPop = n_B = length(tgs.ij) # sum(m.rf_w_treatments.stability_selection[tf.ij,]) #n.tgs.grn_w_treatments_w_motifs.r #sum(m.gn_condition_tissue_differentialExpression[tf.IDs.global,ct.i])
        popSize = n_C = length(colnames(m.grn)) # sum(m.rf_w_treatments.stability_selection) 
        
        # length(tgs.ij) # all targets of TF in condition 
        #popSize = n_C = 
        
        # length(v.tgs.i) # all targets in condition (TFs and TGS)
        failInPop = n_C - n_B
        if(hitInSample >= th.min_number_MR_targets){
          pval <- 1
          fc <- (n_A_B / n_A) / (n_B / n_C)
          if(fc > 1){
            pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= FALSE)
          } #else if(fc < 1){
          #  pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= TRUE)
          #}
          if(pval <= th.pval){
            l.Hierarchy[[i]][tf.ij, tgs_MR.ij] = fc 
            v.tgs_level_next = c(v.tgs_level_next, tf.ij)
          }
        }
      }
      if(length(v.tgs_level_next) == 0 | all(v.tgs_level_next %in% v.tfs_known)){
        b.continue = FALSE
      }else{
        tiers = tiers + 1 
        v.tgs_level_prevs = v.tgs_level_next
        v.tfs_known = unique(c(v.tfs_known, v.tgs_level_next))
        l.Hierarchy_tfs_per_tier[[i]] = c(l.Hierarchy_tfs_per_tier[[i]], list(v.tgs_level_next))
      }
    } 
    
    nb.tfs_per_tier = numeric(tiers)
    if(tiers > 1){
      for(t in 1:(tiers - 1)){
        v.tfs.t = unlist(l.Hierarchy_tfs_per_tier[[i]][t])
        v.tfs.t_next = unlist(l.Hierarchy_tfs_per_tier[[i]][t + 1])
        idx = (which(v.tfs.t %in% v.tfs.t_next))
        v.tfs.t = v.tfs.t[!v.tfs.t %in% v.tfs.t[idx]]
        l.Hierarchy_tfs_per_tier[[i]][[t]] = v.tfs.t
        nb.tfs_per_tier[t] = length(v.tfs.t)
        nb.tfs_per_tier[t + 1] = length(v.tfs.t_next)
      }
    }else{
      l.Hierarchy_tfs_per_tier[[i]][[1]] = unlist(l.Hierarchy_tfs_per_tier[[i]][1])
      nb.tfs_per_tier = length(unlist(l.Hierarchy_tfs_per_tier[[i]][1]))
    }
    
    
    l.Hierarchy_nb_tfs_per_tier[[i]] = nb.tfs_per_tier
    v.number_tiers[i] = tiers
  }
  
  for(i in 1:length(conds)){
    idx = which(l.Hierarchy_nb_tfs_per_tier[[i]] != 0)
    l.Hierarchy_nb_tfs_per_tier[[i]] = l.Hierarchy_nb_tfs_per_tier[[i]][idx]
    v.number_tiers[i] = length(l.Hierarchy_nb_tfs_per_tier[[i]])
  }
  
  ### 
  
  layers <- paste("layer_", 1:max(v.number_tiers), sep="")
  df.masterRegulatorHierarchy <- as.data.frame(matrix(rep(NA, 2 + length(layers)), nrow=1))
  names(df.masterRegulatorHierarchy) <- c("condition", "number_of_tiers",  layers)
  
  for(i in 1:length(l.MR_vs_geneGroups_given_condition)){
    
    df.masterRegulatorHierarchy.i <- as.data.frame(matrix(rep(NA, 2 + length(layers)), nrow=1))
    names(df.masterRegulatorHierarchy.i) <- c("condition", "number_of_tiers",  layers)
    
    df.masterRegulatorHierarchy.i$condition[1] = names(v.number_tiers)[i]
    df.masterRegulatorHierarchy.i$number_of_tiers[1] = v.number_tiers[i]
    
    if( v.number_tiers[i] > 0){
      idx = 3
      for(j in 1:length(l.Hierarchy_nb_tfs_per_tier[[i]])){
        df.masterRegulatorHierarchy.i[1,idx] = l.Hierarchy_nb_tfs_per_tier[[i]][j]
        idx = idx + 1
      }
    }
    df.masterRegulatorHierarchy = rbind(df.masterRegulatorHierarchy, df.masterRegulatorHierarchy.i)
  }
  
  
  return(list(l.Hierarchy=l.Hierarchy, l.Hierarchy_tfs_per_tier=l.Hierarchy_tfs_per_tier, l.Hierarchy_nb_tfs_per_tier = l.Hierarchy_nb_tfs_per_tier, df.masterRegulatorHierarchy=df.masterRegulatorHierarchy, v.number_tiers=v.number_tiers))
  
}





#' Step 4 - Master regulator hierarchy inference
#'
#' This function identifies the master regulator hierarchy
#' @param m.grn = m.grn,
#' @param l.grn_subnetworks = l.grn_subnetworks, 
#' @param df.geneGroups,
#' @param tb.geneGroups,
#' @param v.geneGroups,
#' @param l.geneGroups,
#' @param th.min_number_targets = 2,
#' @param th.min_number_MR_targets = 2,
#' @param th.pval = 0.05
#' @keywords 
#' @export
#' @examples
do_master_regulator_hierarchy_inference = function(m.grn,
                                                   l.grn_subnetworks, 
                                                   df.transcriptionFactorAnnotation,
                                                   df.geneGroups,
                                                   tb.geneGroups,
                                                   v.geneGroups,
                                                   l.geneGroups,
                                                   th.min_number_targets = 2,
                                                   th.min_number_MR_targets = 2,
                                                   th.pval = 0.05, 
                                                   foldername.results = "results/"){
  
  
  df.transcriptionFactorAnnotation = subset(df.transcriptionFactorAnnotation, df.transcriptionFactorAnnotation$with_geneExpression == "yes")

  v.tfs = unique(df.transcriptionFactorAnnotation$TF_ID)
  v.tf_families = character(length(v.tfs))
  names(v.tf_families) = v.tfs
  
  for(i in 1:length(v.tfs)){
    df.transcriptionFactorAnnotation.i = subset(df.transcriptionFactorAnnotation, df.transcriptionFactorAnnotation$TF_ID == v.tfs[i])
    v.tf_families[i] = paste(df.transcriptionFactorAnnotation.i$TF_Fam, collapse = ", ")
  }
   
    
  message("----- Identify bottom tier master regulators -----")
  
  # identify level 1 master regulators (bottom level - tfs x domains)
  l.res.MR <- identify_bottom_tier_masterRegulators(m.grn,
                                                 l.grn_subnetworks,
                                                 v.tfs = rownames(m.grn),
                                                 v.conds = names(l.grn_subnetworks),
                                                 df.geneGroups,
                                                 tb.geneGroups,
                                                 v.geneGroups,
                                                 l.geneGroups,
                                                 th.min_number_targets,
                                                 th.pval,
                                                 b.include_under_represented = "yes")

  m.MR_vs_conditions <- l.res.MR$m.MR_vs_conditions  # A) TFs versus Conditions (Matrix plot) P(TF,C)
  l.MR_vs_geneGroups_given_condition <- l.res.MR$l.MR_vs_geneGroups_given_condition  # B) per condition - TFs versus Domains (P(TF,D|C)) => also cumulative plot 
  
  # check for 
  m.MR_vs_conditions[!is.na(m.MR_vs_conditions)] = 1
  m.MR_vs_conditions[is.na(m.MR_vs_conditions)] = 0
  number_of_conditions_per_master_regulator = rowSums(m.MR_vs_conditions)
  
  
  message("----- Infer master regulator regulatory hierarchy -----")
  
  l.Hierarchy = vector(mode = "list", length = 2)
  l.Hierarchy_tfs_per_tier = vector(mode = "list", length = 2)
  l.Hierarchy_nb_tfs_per_tier = vector(mode = "list", length = 2)
  l.df.masterRegulatorHierarchy = vector(mode = "list", length = 2)
  v.number_tiers = vector(mode = "list", length = 2)
  
  names(l.Hierarchy) = c("genes", "geneGroups")
  names(l.Hierarchy_tfs_per_tier) = c("genes", "geneGroups")
  names(l.Hierarchy_nb_tfs_per_tier) = c("genes", "geneGroups")
  names(l.df.masterRegulatorHierarchy) = c("genes", "geneGroups")
  names(v.number_tiers) = c("genes", "geneGroups")
  
  # per condition
  l.res = identify_regulatory_hierachy(m.MR_vs_conditions, 
                                       l.MR_vs_geneGroups_given_condition, 
                                       m.grn,
                                       l.grn_subnetworks,
                                       v.tfs = rownames(m.grn),
                                       v.conds = names(l.grn_subnetworks),
                                       df.geneGroups,
                                       tb.geneGroups,
                                       v.geneGroups,
                                       l.geneGroups,
                                       th.min_number_MR_targets = th.min_number_MR_targets,
                                       mode = "genes",
                                       th.pval = th.pval)
    
  l.Hierarchy[[1]]=l.res$l.Hierarchy
  l.Hierarchy_tfs_per_tier[[1]]=l.res$l.Hierarchy_tfs_per_tier
  l.Hierarchy_nb_tfs_per_tier[[1]] = l.res$l.Hierarchy_nb_tfs_per_tier
  l.df.masterRegulatorHierarchy[[1]]=l.res$df.masterRegulatorHierarchy
  v.number_tiers[[1]]=l.res$v.number_tiers
  
  
  l.res = identify_regulatory_hierachy(m.MR_vs_conditions, 
                                       l.MR_vs_geneGroups_given_condition, 
                                       m.grn,
                                       l.grn_subnetworks,
                                       v.tfs = rownames(m.grn),
                                       v.conds = names(l.grn_subnetworks),
                                       df.geneGroups,
                                       tb.geneGroups,
                                       v.geneGroups,
                                       l.geneGroups,
                                       th.min_number_MR_targets = th.min_number_MR_targets,
                                       mode = "geneGroups",
                                       th.pval = th.pval)
    

  l.Hierarchy[[2]]=l.res$l.Hierarchy
  l.Hierarchy_tfs_per_tier[[2]]=l.res$l.Hierarchy_tfs_per_tier
  l.Hierarchy_nb_tfs_per_tier[[2]] = l.res$l.Hierarchy_nb_tfs_per_tier
  l.df.masterRegulatorHierarchy[[2]]=l.res$df.masterRegulatorHierarchy
  v.number_tiers[[2]]=l.res$v.number_tiers
  
 
  files = character(2)
  files[1] = paste(foldername.results, "genes/", sep = "")
  files[2] = paste(foldername.results, "geneGroups_masterRegulatorHierarchies/", sep = "")
  # 
  # if(!file.exists(paste(foldername.results, "genes/", sep = ""))){
  #   dir.create(paste(foldername.results, "genes/", sep = ""))
  # }
  if(!file.exists(paste(foldername.results, "geneGroups_masterRegulatorHierarchies/", sep = ""))){
    dir.create(paste(foldername.results, "geneGroups_masterRegulatorHierarchies/", sep = ""))
  }

  s = 2
  for(i in 1:length(v.number_tiers[[s]])){
    
    if(v.number_tiers[[s]][i] > 0){

      ct.i = names(v.number_tiers[[s]])[i]
      tfs.i = rownames(l.MR_vs_geneGroups_given_condition[ct.i][[1]])
      tgs.i = colnames(l.MR_vs_geneGroups_given_condition[ct.i][[1]])
      l.Hierarchy[[s]][[i]][tfs.i, tgs.i] = l.MR_vs_geneGroups_given_condition[ct.i][[1]]
      
      tfs.i = paste(rownames(l.Hierarchy[[s]][[i]]), "(", v.tf_families[rownames(l.Hierarchy[[s]][[i]])], ")", sep ="")
      tgs.i = colnames(l.Hierarchy[[s]][[i]])
      tgs_tfs.i = tgs.i[!tgs.i %in% v.geneGroups]
      tgs_tfs_w_fams.i = paste(tgs_tfs.i, "(", v.tf_families[tgs_tfs.i], ")", sep ="")
      
      m.net = l.Hierarchy[[s]][[i]]
      rownames(m.net) = tfs.i
      
      idx = which(colnames(m.net) %in% v.geneGroups)
      v.geneGroups.i = colnames(m.net)[idx]
      
      colnames(m.net) = c(tgs_tfs_w_fams.i, v.geneGroups.i)
      
      m.net[is.na(m.net)] = 0
      m.net[m.net > 1] = 1
      
      #g <- graph_from_adjacency_matrix(m.net)
      
      df.idx.MR_hierarchy.i <- which(m.net > 0, arr.ind = TRUE)
      
      if(nrow(df.idx.MR_hierarchy.i) > 0){
        
        df.MR_hierarchy.i <- data.frame(TF = rownames(m.net)[df.idx.MR_hierarchy.i[,1]], TG = colnames(m.net)[df.idx.MR_hierarchy.i[,2]], stringsAsFactors = FALSE)
        
        g <- graph_from_data_frame(df.MR_hierarchy.i, directed = TRUE)
        wc <- cluster_walktrap(g)
        members <- membership(wc)
        
        #Convert to object suitable for networkD3
        net_d3  <- igraph_to_networkD3(g,group = members)
        
        # filename = paste(output_folder, "Condition_specific_analyses\\df.regulatory_hiearchy_MR_per_tier_bottom_genes_enzymes.csv", sep = "")
        # filename = paste(output_folder, "Condition_specific_analyses\\df.regulatory_hiearchy_MR_per_tier_bottom_genes_enzymes.csv", sep = "")
        # subDir <- files[s]#  paste(output_folder, "C:/Users/Michael/Documents/Computational_Biology/MERIT_V1.0/figures/Analysis_per_condition/masterRegulatorHierarchies/", sep ="")
        
        subDir = files[s]
        filename <- paste(getwd(), "/", subDir, ct.i,".html", sep ="")
        net_d3$links["value"] <- 1
        net_d3$nodes <- remove.factors(net_d3$nodes)
        net_d3$nodes$name[!net_d3$nodes$name %in% v.geneGroups] <- paste(net_d3$nodes$name[!net_d3$nodes$name %in% v.geneGroups], "(", v.tf_families[as.character(net_d3$nodes$name[!net_d3$nodes$name %in% v.geneGroups])], ")", sep ="")
        #
        sankeyNetwork(Links = net_d3$links, Nodes = net_d3$nodes,
                      Source = "source", Target = "target",
                      NodeID = "name", Value = "value", 
                      fontSize = 20, nodeWidth = 20, units = "Letter(s)", fontFamily = "sans-serif", iterations = 0,
                      nodePadding = 1, height = 5000, width = 10000, sinksRight =TRUE, margin = NULL)%>%
          htmlwidgets::prependContent(htmltools::tags$h1(paste("Master Regulator Hierarchy in condition: ",ct.i))) %>%
          saveNetwork(file = filename)
      }
    }
  }
  
  return(list(l.Hierarchy=l.Hierarchy, 
              l.Hierarchy_tfs_per_tier=l.Hierarchy_tfs_per_tier,
              l.Hierarchy_nb_tfs_per_tier=l.Hierarchy_nb_tfs_per_tier,
              l.df.masterRegulatorHierarchy=l.df.masterRegulatorHierarchy,
              v.number_tiers=v.number_tiers,
              m.MR_vs_conditions = l.res.MR$m.MR_vs_conditions,  # A) TFs versus Conditions (Matrix plot) P(TF,C)
              l.MR_vs_geneGroups_given_condition = l.res.MR$l.MR_vs_geneGroups_given_condition,  # B) per condition - TFs versus Domains (P(TF,D|C)) => also cumulative plot 
              number_of_conditions_per_master_regulator=number_of_conditions_per_master_regulator))
  
}



#' format_results
#'
#' This function formats all results
#' @keywords 
#' @export
format_results = function(l.grn_subnetworks,
                          tb.condition_tissue_differentialExpression,
                          l.Hierarchy, 
                          l.Hierarchy_tfs_per_tier,
                          l.Hierarchy_nb_tfs_per_tier,
                          l.df.masterRegulatorHierarchy,
                          v.number_tiers,
                          m.MR_vs_conditions,  # A) TFs versus Conditions (Matrix plot) P(TF,C)
                          l.MR_vs_geneGroups_given_condition,  # B) per condition - TFs versus Domains (P(TF,D|C)) => also cumulative plot 
                          number_of_conditions_per_master_regulator,
                          tb.condition_treatments,
                          tb.condition_tissues,
                          df.transcriptionFactorAnnotation, 
                          df.geneGroups,
                          tb.geneGroups,
                          v.geneGroups,
                          l.geneGroups,
                          th.pval = 0.05,
                          foldername.results = "results/"){
  
  df.transcriptionFactorAnnotation = subset(df.transcriptionFactorAnnotation, df.transcriptionFactorAnnotation$with_geneExpression == "yes")
  v.tfs = unique(df.transcriptionFactorAnnotation$TF_ID)
  v.tf_families = character(length(v.tfs))
  names(v.tf_families) = v.tfs
  for(i in 1:length(v.tfs)){
    df.transcriptionFactorAnnotation.i = subset(df.transcriptionFactorAnnotation, df.transcriptionFactorAnnotation$TF_ID == v.tfs[i])
    v.tf_families[i] = paste(df.transcriptionFactorAnnotation.i$TF_Fam, collapse = ", ")
  }
  
  tb.tf_families = table(v.tf_families)
  
  df.transcriptionFactorAnnotation = subset(df.transcriptionFactorAnnotation, df.transcriptionFactorAnnotation$with_geneExpression == "yes")
  df.geneGroups = subset(df.geneGroups, df.geneGroups$with_geneExpression == "yes")
  v.tfs = unique(df.transcriptionFactorAnnotation$TF_ID)
  v.genes =  unique(c(v.tfs, rownames(df.geneGroups)))
  v.gns_geneGroups = rownames(df.geneGroups)
  
  ####
  
  filename = paste(foldername.results, "cumulative_numbers_of_active_genes_in_condition_groups.pdf", sep = "")
  df.data <- data.frame(cond = names(tb.condition_tissue_differentialExpression), nr = as.numeric(tb.condition_tissue_differentialExpression))
  p <- ggplot(df.data, aes(cond, nr))
  p <- p + geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("cumulative number of active (expressed and repressed) genes at differential expression p-value < 0.05")
  pdf(filename)
  p
  dev.off()
  
  ### 
  
  v.fams <- unique(v.tf_families)
  df.data <- data.frame(TF = character(), Fam = character(), Cond = character())
  
  for(i in 1:length(l.grn_subnetworks)){
    m.grn <- l.grn_subnetworks[[i]]
    df.idx.grn <- which(m.grn == 1, arr.ind = TRUE)
    df.grn <- data.frame(TF = rownames(m.grn)[df.idx.grn[,1]], TG = colnames(m.grn)[df.idx.grn[,2]], stringsAsFactors = FALSE)
    tfs.grn <- unique(df.grn$TF)
    ct.i <- names(l.grn_subnetworks)[i]
    for(j in 1:length(v.fams)){
      # family tfs in the grn (diff exp)
      tf.IDs <- names(v.tf_families[tfs.grn])[which(v.tf_families[tfs.grn] == v.fams[j])]
      # tf.IDs.global <- intersect(rownames(m.rf), names(v.tf_families)[which(v.tf_families == v.fams[j])])
      if(length(tf.IDs) > 0){
        df.data <- rbind(df.data, 
                         data.frame(TF = tf.IDs, 
                                    Fam = rep(v.fams[j] , length(tf.IDs)), 
                                    Cond = rep(ct.i, length(tf.IDs))))
      }
    }
  }
  
  # how many conditions does a TF occur in 
  df.data <- data.frame(TF = as.numeric(table(df.data$TF)) / length(l.grn_subnetworks)  * 100, Fam = as.character(v.tf_families[names(table(df.data$TF))]))
  p1 <- ggplot(df.data, aes(x=Fam, y=TF, fill=Fam))
  p1 <- p1 + geom_boxplot() + scale_fill_grey() + theme_bw() + theme(axis.text.x = element_text(angle = 270, hjust = 0)) + guides(fill=FALSE) + xlab("Transcription Factor family") + ylab("Percentage of conditions of TF family activity (%)")
  
  
  filename <- paste(foldername.results, "number_of_conditions_per_member_of_TFfam_in_percentage.pdf", sep = "")
  pdf(filename, 10, 8)
  p1
  dev.off()
  
  ######

  
  # Condition specific networks
  subDir = paste(foldername.results, "condition_specific_networks/", sep = "")
  # filename <- paste(getwd(),"/",subDir, ct.i,".html", sep ="")
  
  if (!file.exists(subDir)){
    dir.create(subDir)
  }
  
  for(i in 1:length(l.grn_subnetworks)){
    m.grn <- l.grn_subnetworks[[i]]
    df.idx.grn <- which(m.grn> 0, arr.ind = TRUE)
    df.grn <- data.frame(TF = rownames(m.grn)[df.idx.grn[,1]], TG = colnames(m.grn)[df.idx.grn[,2]], stringsAsFactors = FALSE)
    filename <- paste(subDir, names(l.grn_subnetworks)[i], ".csv", sep ="")
    write.csv2(df.grn, filename, row.names = FALSE)
  }
  

  # Master Regulator
  if(TRUE){
    hist(number_of_conditions_per_master_regulator / dim(m.MR_vs_conditions)[2] * 100)
    
    m.heatmap <- t(m.MR_vs_conditions)
    df.heatmap <- melt(m.heatmap)
    df.heatmap$value <- log(df.heatmap$value)
    df.heatmap$value[df.heatmap$value == -Inf] <- 0
    
    df.MR_vs_conditions <- df.heatmap
    df.MR_vs_conditions <- subset(df.MR_vs_conditions, !is.na(df.MR_vs_conditions$value))
    
    df.heatmap$Var2 <- paste(df.heatmap$Var2, "(", v.tf_families[df.heatmap$Var2], ")", sep = "")
    
    df.heatmap$value[df.heatmap$value < log(0.1)] <- log(0.1)
    df.heatmap$value[df.heatmap$value > log(10)] <- log(10)
    
    v.breaks = c(log(0.1), log(0.25), log(0.5), log(1), log(2), log(4), log(10))
    v.lables <- c("< 0.1", "0.25", "0.5","1","2","4", "> 10")#as.character(exp(v.breaks))
    
    ggheatmap <- ggplot(df.heatmap, aes(Var1, Var2, fill = value)) + geom_tile(color = "black") + scale_fill_gradient2(low = "yellow", high = "blue", mid = "black", midpoint = 0, name = "fold change", breaks=v.breaks, labels=v.lables) + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1, color = "black")) + theme(axis.title.y=element_blank(), axis.text.y = element_text(angle = 0, vjust = , size = 10, hjust = 1, color = "black"))
    
    filename <- paste(foldername.results, "masterRegulators_vs_conditions.pdf", sep = "")
    pdf(filename, width=15, height=120)
    plot(ggheatmap)
    dev.off()
    
    
    ####
    
    v.MRs = rownames(m.MR_vs_conditions)
    v.Conds = colnames(m.MR_vs_conditions)

    m.MR_vs_geneGroups_across_conditions <- matrix(0, nrow = length(v.MRs), ncol = length(v.geneGroups), dimnames = list(v.MRs, v.geneGroups))
    for(i in 1:length(v.Conds)){
      m.heatmap <- l.MR_vs_geneGroups_given_condition[v.Conds[i]][[1]]
      m.heatmap[is.na(m.heatmap)] <- 0
      m.heatmap[m.heatmap < 0] <- 0
      m.heatmap[m.heatmap > 0] <- 1
      m.MR_vs_geneGroups_across_conditions[, colnames(m.heatmap)] <- m.MR_vs_geneGroups_across_conditions[, colnames(m.heatmap)] + m.heatmap[rownames(m.MR_vs_geneGroups_across_conditions),]
    }
    m.MR_vs_geneGroups_across_conditions <- round(m.MR_vs_geneGroups_across_conditions / length(v.Conds) * 100)
    df.heatmap <- melt(t(m.MR_vs_geneGroups_across_conditions))
    df.heatmap$Var2 <- paste(df.heatmap$Var2, "(", v.tf_families[df.heatmap$Var2], ")", sep = "")
    ggheatmap <- ggplot(df.heatmap, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red", name = "Percentage of total master regulators") + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1, color = "black")) + theme(axis.title.y=element_blank(), axis.text.y = element_text(angle = 0, vjust = , size = 8, hjust = 1, color = "black"))
  
    filename <- paste(foldername.results, "masterRegulator_vs_GeneGroups_across_Conditions.pdf", sep = "")
    pdf(filename, width=10, height=100)
    plot(ggheatmap)
    dev.off()
    
    
    filename <- paste(foldername.results, "masterRegulator_vs_GeneGroups_across_Conditions_clustered.pdf", sep = "")
    m.heatmap <- m.MR_vs_geneGroups_across_conditions
    
    rownames(m.heatmap) <- paste(rownames(m.heatmap), "(", v.tf_families[rownames(m.heatmap)], ")", sep = "")
    
    paletteLength <- 10
    myColor <- colorRampPalette(c("white", "red"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- seq(min(m.heatmap), max(m.heatmap), length.out=paletteLength)
    
    p <- pheatmap(m.heatmap,  color=myColor, breaks=myBreaks)
    save_pheatmap_pdf(p, filename, width=5.5, height=110)
  }
  ####
  
  if(TRUE){
    filename <- paste(foldername.results, "df.masterRegulatorHierarchy_genes.csv", sep = "")
    write.table(l.df.masterRegulatorHierarchy[[1]], filename, row.names = F, sep = ";" )
    
    filename <- paste(foldername.results, "df.masterRegulatorHierarchy_geneGroups.csv", sep = "")
    write.table(l.df.masterRegulatorHierarchy[[2]], filename, row.names = F, sep = ";" )
  }
  
  #####

  #####


  if(TRUE){
  
    v.tf.fams.selection = unique(v.tf_families[rownames(m.MR_vs_geneGroups_across_conditions)])
    l.fams_vs_geneGroups = vector(mode = "list", length = 3)
    for(i in 1:2){l.fams_vs_geneGroups[[i]] = matrix(0, nrow = length(v.tf.fams.selection), ncol = length(v.geneGroups), dimnames = list(v.tf.fams.selection, v.geneGroups))}
    l.fams_vs_geneGroups[[3]] = matrix(1, nrow = length(v.tf.fams.selection), ncol = length(v.geneGroups), dimnames = list(v.tf.fams.selection, v.geneGroups))
    names(l.fams_vs_geneGroups) =c ("numbers", "percentage", "enrichment")
    
    for(i in 1:length(v.tf.fams.selection)){
      
      fam.i = v.tf.fams.selection[i]
      idx = which(v.tf_families[rownames(m.MR_vs_geneGroups_across_conditions)] == fam.i)
      
      for(j in 1:length(v.geneGroups)){
        
        v.ij.all = m.MR_vs_geneGroups_across_conditions[, v.geneGroups[j]]
        v.ij = m.MR_vs_geneGroups_across_conditions[idx, v.geneGroups[j]]
        v.tfs.ij = names(which(v.ij > 0))
        
        l.fams_vs_geneGroups[[1]][i,j] = length(v.tfs.ij)
        l.fams_vs_geneGroups[[2]][i,j] = length(v.tfs.ij) / tb.tf_families[fam.i] * 100
        
        # enrichment of the particular activity
        
        ## domain genes in target genes - condition dependent 
        hitInSample = n_A_B = length(v.tfs.ij)
        sampleSize = n_A = length(which(v.ij.all > 0))
        
        ## condition independent part 
        hitInPop = n_B = tb.tf_families[fam.i]
        popSize = n_C = sum(tb.tf_families)
        failInPop = n_C-n_B
        
        if(hitInSample >= 1){
          pval <- 1
          fc <- (n_A_B / n_A) / (n_B / n_C)
          if(fc > 1){
            pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= FALSE)
          }else if(fc < 1){
            pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= TRUE)
          }
          if(pval <= th.pval){
            l.fams_vs_geneGroups[[3]][i,j] <- fc
          }
        }
      }
    }
    
  
    subDir = paste(foldername.results, "family_v_geneGroups/", sep = "")
  
    if (!file.exists(subDir)){
      dir.create(subDir)
    }
    
    m.heatmap <- l.fams_vs_geneGroups[[1]]
    df.heatmap <- melt(t(m.heatmap))
    ggheatmap <- ggplot(df.heatmap, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red", name = "Numbers") + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1, color = "black")) + theme(axis.title.y=element_blank(), axis.text.y = element_text(angle = 0, vjust = , size = 8, hjust = 1, color = "black"))
    filename <- paste(subDir, "numbers.pdf", sep ="")
    pdf(filename, width=10, height=10)
    plot(ggheatmap)
    dev.off()
    
    m.heatmap <- l.fams_vs_geneGroups[[2]]
    df.heatmap <- melt(t(m.heatmap))
    ggheatmap <- ggplot(df.heatmap, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red", name = "Percentage") + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1, color = "black")) + theme(axis.title.y=element_blank(), axis.text.y = element_text(angle = 0, vjust = , size = 8, hjust = 1, color = "black"))
    filename <- paste(subDir, "percentage.pdf", sep ="")
    pdf(filename, width=10, height=10)
    plot(ggheatmap)
    dev.off()
    
    m.heatmap <- l.fams_vs_geneGroups[[3]]
    m.heatmap = log(m.heatmap)
    m.heatmap[m.heatmap == -Inf] = 0
    m.heatmap[m.heatmap == 0] = NA
    df.heatmap <- melt(t(m.heatmap))
    ggheatmap <- ggplot(df.heatmap, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red", name = "log FC") + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1, color = "black")) + theme(axis.title.y=element_blank(), axis.text.y = element_text(angle = 0, vjust = , size = 8, hjust = 1, color = "black"))
    filename <- paste(subDir, "enrichment.pdf", sep ="")
    pdf(filename, width=10, height=10)
    plot(ggheatmap)
    dev.off()
    
    
   
    v.tf.fams.selection = unique(v.tf_families[rownames(m.MR_vs_conditions)])
    v.conds.selection = colnames(m.MR_vs_conditions)
  
    l.fams_vs_conds = vector(mode = "list", length = 3)
    for(i in 1:2)
      l.fams_vs_conds[[i]] = matrix(0, nrow = length(v.tf.fams.selection), ncol = length(v.conds.selection), dimnames = list(v.tf.fams.selection, v.conds.selection))
    l.fams_vs_conds[[3]] = matrix(1, nrow = length(v.tf.fams.selection), ncol = length(v.conds.selection), dimnames = list(v.tf.fams.selection, v.conds.selection))
    names(l.fams_vs_conds) = c("numbers", "percentage", "enrichment")
    
    for(i in 1:length(v.tf.fams.selection)){
      fam.i = v.tf.fams.selection[i]
      idx = which(v.tf_families[rownames(m.MR_vs_conditions)] == fam.i)
      for(j in 1:length(v.conds.selection)){
        v.ij.all = m.MR_vs_conditions[, v.conds.selection[j]]
        v.ij = m.MR_vs_conditions[idx, v.conds.selection[j]]
        v.tfs.ij = names(which(v.ij > 0))
        l.fams_vs_conds[[1]][i,j] = length(v.tfs.ij)
        l.fams_vs_conds[[2]][i,j] = length(v.tfs.ij) / tb.tf_families[fam.i] * 100
        
        # enrichment of the particular activity
        ## domain genes in target genes - condition dependent 
        hitInSample = n_A_B = length(v.tfs.ij)
        sampleSize = n_A = length(which(v.ij.all > 0))
        ## condition independent part 
        hitInPop = n_B = tb.tf_families[fam.i]
        popSize = n_C = sum(tb.tf_families)
        failInPop = n_C-n_B
        
        if(hitInSample >= 1){
          pval <- 1
          fc <- (n_A_B / n_A) / (n_B / n_C)
          if(fc > 1){
            pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= FALSE)
          }else if(fc < 1){
            pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= TRUE)
          }
          
          if(pval <= th.pval){
            l.fams_vs_conds[[3]][i,j] <- fc
            
          }
        }
        
      }
      
    }
    
    subDir = paste(foldername.results, "family_v_conditions/", sep = "")
    
    if (!file.exists(subDir)){
      dir.create(subDir)
    }
    
    m.heatmap <- l.fams_vs_conds[[1]]
    df.heatmap <- melt(t(m.heatmap))
    ggheatmap <- ggplot(df.heatmap, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red", name = "Numbers") + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1, color = "black")) + theme(axis.title.y=element_blank(), axis.text.y = element_text(angle = 0, vjust = , size = 8, hjust = 1, color = "black"))
    filename <- paste(subDir, "numbers.pdf", sep ="")
    pdf(filename, width=10, height=10)
    plot(ggheatmap)
    dev.off()
    
    m.heatmap <- l.fams_vs_conds[[2]]
    df.heatmap <- melt(t(m.heatmap))
    ggheatmap <- ggplot(df.heatmap, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red", name = "Percentage") + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1, color = "black")) + theme(axis.title.y=element_blank(), axis.text.y = element_text(angle = 0, vjust = , size = 8, hjust = 1, color = "black"))
    filename <- paste(subDir, "percentage.pdf", sep ="")
    pdf(filename, width=10, height=10)
    plot(ggheatmap)
    dev.off()
    
    m.heatmap <- l.fams_vs_conds[[3]]
    m.heatmap = log(m.heatmap)
    m.heatmap[m.heatmap == -Inf] = 0
    m.heatmap[m.heatmap == 0] = NA
    df.heatmap <- melt(t(m.heatmap))
    ggheatmap <- ggplot(df.heatmap, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red", name = "log FC") + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1, color = "black")) + theme(axis.title.y=element_blank(), axis.text.y = element_text(angle = 0, vjust = , size = 8, hjust = 1, color = "black"))
    filename <- paste(subDir, "enrichment.pdf", sep ="")
    pdf(filename, width=10, height=10)
    plot(ggheatmap)
    dev.off()
    
  
    tb.geneGroup_activity = rep(0, length(v.geneGroups))
    names(tb.geneGroup_activity) = v.geneGroups
    
    b.tfs.active = rep(0, length(rownames(m.MR_vs_conditions)))
    names(b.tfs.active) = rownames(m.MR_vs_conditions)
    
    for(j in 1:length(v.geneGroups)){
      for(i in 1:length(v.conds.selection)){
        m.MR_vs_geneGroups_given_condition = l.MR_vs_geneGroups_given_condition[[i]]
        idx = which(m.MR_vs_geneGroups_given_condition[,v.geneGroups[j]] > 0)
        b.tfs.active[idx] = 1
      }
      tb.geneGroup_activity[j] = sum(b.tfs.active)
    }
    
    v.tf.fams.selection = unique(v.tf_families[rownames(m.MR_vs_conditions)])
    v.conds.selection = names(l.MR_vs_geneGroups_given_condition)
    
    m.fams_vs_geneGroups = matrix(0, nrow = length(v.conds.selection), ncol = length(v.geneGroups), dimnames = list(v.conds.selection, v.geneGroups))
    
    l.geneGroups_vs_conds = vector(mode = "list", length = 3)
    for(i in 1:2) 
      l.geneGroups_vs_conds[[i]] = matrix(0, nrow = length(v.conds.selection), ncol = length(v.geneGroups), dimnames = list(v.conds.selection, v.geneGroups))
    l.geneGroups_vs_conds[[3]] = matrix(1, nrow = length(v.conds.selection), ncol = length(v.geneGroups), dimnames = list(v.conds.selection, v.geneGroups))
    names(l.geneGroups_vs_conds) =c ("numbers", "percentage", "enrichment")
    
    for(i in 1:length(v.conds.selection)){
      m.MR_vs_geneGroups_given_condition = l.MR_vs_geneGroups_given_condition[[i]]
  
      for(j in 1:length(v.geneGroups)){
        
        m.ij = m.MR_vs_geneGroups_given_condition
        m.ij[is.na(m.ij)] = 0
        
        v.tfs.ij = names(which(m.MR_vs_geneGroups_given_condition[,v.geneGroups[j]] > 0))
        
        # number and percentage 
        l.geneGroups_vs_conds[[1]][i,j] = length(v.tfs.ij)
        l.geneGroups_vs_conds[[2]][i,j] = length(v.tfs.ij) / length(v.tf_families)
        
        # enrichment 
        hitInSample = n_A_B = length(v.tfs.ij) # active in domain given condition
        sampleSize = n_A = length(which(rowSums(m.ij) > 0)) # active across domains given condition
        
        ## condition independent part 
        hitInPop = n_B =  tb.geneGroup_activity[v.geneGroups[j]] # active in domain across conditions
        popSize = n_C = sum(tb.geneGroup_activity) # active in all domains across conditions
        
        if(hitInSample >= 1){
          pval <- 1
          fc <- (n_A_B / n_A) / (n_B / n_C)
          if(fc > 1){
            pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= FALSE)
          }else if(fc < 1){
            pval <- phyper(n_A_B, n_B, n_C-n_B, n_A,lower.tail= TRUE)
          }
          
          if(pval <= th.pval){
            l.geneGroups_vs_conds[[3]][i,j] <- fc
            
          }
        }
        
      }
      
    }
  
    
    
    subDir = paste(foldername.results, "geneGroups_v_conditions/", sep = "")
    
    if (!file.exists(subDir)){
      dir.create(subDir)
    }
    
    m.heatmap <- l.geneGroups_vs_conds[[1]]
    df.heatmap <- melt(t(m.heatmap))
    ggheatmap <- ggplot(df.heatmap, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red", name = "Numbers") + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1, color = "black")) + theme(axis.title.y=element_blank(), axis.text.y = element_text(angle = 0, vjust = , size = 8, hjust = 1, color = "black"))
    filename <- paste(subDir, "numbers.pdf", sep ="")
    pdf(filename, width=10, height=10)
    plot(ggheatmap)
    dev.off()
    
    m.heatmap <- l.geneGroups_vs_conds[[2]]
    df.heatmap <- melt(t(m.heatmap))
    ggheatmap <- ggplot(df.heatmap, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red", name = "Percentage") + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1, color = "black")) + theme(axis.title.y=element_blank(), axis.text.y = element_text(angle = 0, vjust = , size = 8, hjust = 1, color = "black"))
    filename <- paste(subDir, "percentage.pdf", sep ="")
    pdf(filename, width=10, height=10)
    plot(ggheatmap)
    dev.off()
    
    m.heatmap <- l.geneGroups_vs_conds[[3]]
    m.heatmap = log(m.heatmap)
    m.heatmap[m.heatmap == -Inf] = 0
    m.heatmap[m.heatmap == 0] = NA
    df.heatmap <- melt(t(m.heatmap))
    ggheatmap <- ggplot(df.heatmap, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red", name = "log FC") + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1, color = "black")) + theme(axis.title.y=element_blank(), axis.text.y = element_text(angle = 0, vjust = , size = 8, hjust = 1, color = "black"))
    filename <- paste(subDir, "enrichment.pdf", sep ="")
    pdf(filename, width=10, height=10)
    plot(ggheatmap)
    dev.off()
    
    
  }
    
  if(TRUE){
    
    df.ratio = data.frame(condition= character(length(l.grn_subnetworks)),
                          grn_size = numeric(length(l.grn_subnetworks)), 
                          tmt_number = numeric(length(l.grn_subnetworks)),
                          stringsAsFactors = F)
    
    v.sub_net = unlist(lapply(l.grn_subnetworks, sum))
    
    for(i in 1:length(v.sub_net)){
      
      condition = unlist(strsplit(names(v.sub_net)[i],  " - "))
      df.ratio$condition[i] = condition[1]
      df.ratio$grn_size[i] = as.numeric(v.sub_net[i])
      df.ratio$tmt_number[i] = as.numeric(tb.condition_treatments[condition[1]]) 
      
    }
    
    
    
    df.ratio["ratio"] = round(df.ratio$grn_size / df.ratio$tmt_number,2)

    
    filename <- paste(foldername.results, "df.ratio.csv", sep ="")
    write.table(df.ratio, filename, sep = ";", row.names = F)
    
    
  }

  if(TRUE){
    
    #library(xtable)
    
    df.condition_network_statistics = data.frame(condition = names(l.grn_subnetworks))
    df.condition_network_statistics["TFs"] = unlist(lapply(l.grn_subnetworks, function(m) length(which(rowSums(m) > 0))))
    df.condition_network_statistics["TF families"] = unlist(lapply(l.grn_subnetworks, function(m) length(unique(v.tf_families[(which(rowSums(m) > 0))]))))
    
    df.condition_network_statistics["TGs"] = unlist(lapply(l.grn_subnetworks, function(m) length(which(colSums(m) > 0))))
    df.condition_network_statistics["Links"] = unlist(lapply(l.grn_subnetworks, function(m) length(which(m > 0))))
    
    df.condition_network_statistics["Ratio TFs vs Enzymes"] = 0
    df.condition_network_statistics["Enzymes"] = 0
    
    df.condition_network_enzymes_statistics = data.frame(condition = names(l.grn_subnetworks))
    df.condition_network_enzymes_statistics["Enzymes"] = 0
    
    # df.condition_network_enzymes_statistics["amines and polyamines"] = 0
    # df.condition_network_enzymes_statistics["amino acids"] = 0
    # df.condition_network_enzymes_statistics["carbohydrates"] = 0
    # df.condition_network_enzymes_statistics["cofactors"] = 0
    # df.condition_network_enzymes_statistics["detoxification"] = 0
    # df.condition_network_enzymes_statistics["energy"] = 0
    # df.condition_network_enzymes_statistics["fatty acids and lipids"] = 0
    # df.condition_network_enzymes_statistics["hormones"] = 0
    # df.condition_network_enzymes_statistics["inorganic nutrients"] = 0
    # df.condition_network_enzymes_statistics["nucleotides"] = 0
    # df.condition_network_enzymes_statistics["redox"] = 0
    # df.condition_network_enzymes_statistics["primary specialized interfaced"] = 0
    # df.condition_network_enzymes_statistics["specialized"] = 0
    # df.condition_network_enzymes_statistics["other"] = 0
    
    for(i in 1:length(l.grn_subnetworks)){
      
      m.grn = l.grn_subnetworks[[i]]
      
      df.idx.grn <- which(m.grn > 0, arr.ind = TRUE)
      df.grn <- data.frame(TF = rownames(m.grn)[df.idx.grn[,1]], TG = colnames(m.grn)[df.idx.grn[,2]], stringsAsFactors = FALSE)
      names(df.grn) <- c("TF","Target")
      
      n.gns_geneGroups.i = length(intersect(df.grn$Target, v.gns_geneGroups))
      
      df.condition_network_statistics[i,"Ratio TFs vs Enzymes"] = as.character(round((df.condition_network_statistics$TFs[i] / n.gns_geneGroups.i),2))
      
      df.condition_network_statistics$Enzymes[i] = as.character(n.gns_geneGroups.i)
      df.condition_network_enzymes_statistics$Enzymes[i] = as.character(n.gns_geneGroups.i)
      if(n.gns_geneGroups.i > 0){
        for(j in 1:length(v.geneGroups)){
          id = paste("", v.geneGroups[j], sep = "")
          n.enz_dom.j = length(intersect(df.grn$Target, unlist(l.geneGroups[id])))
          percentage = round(n.enz_dom.j / n.gns_geneGroups.i * 100,0)
          df.condition_network_enzymes_statistics[i, v.geneGroups[j]] =  paste(n.enz_dom.j, " (", percentage, "%)", sep = "")
        }
      }
    }
    
    # 
    # df.condition_network_pathways_statistics = data.frame(condition = names(l.grn_subnetworks))
    # df.condition_network_pathways_statistics["Pathways"] = 0
    # 
    # df.condition_network_pathways_statistics["amines and polyamines"] = 0
    # df.condition_network_pathways_statistics["amino acids"] = 0
    # df.condition_network_pathways_statistics["carbohydrates"] = 0
    # df.condition_network_pathways_statistics["cofactors"] = 0
    # df.condition_network_pathways_statistics["detoxification"] = 0
    # df.condition_network_pathways_statistics["energy"] = 0
    # df.condition_network_pathways_statistics["fatty acids and lipids"] = 0
    # df.condition_network_pathways_statistics["hormones"] = 0
    # df.condition_network_pathways_statistics["inorganic nutrients"] = 0
    # df.condition_network_pathways_statistics["nucleotides"] = 0
    # df.condition_network_pathways_statistics["redox"] = 0
    # df.condition_network_pathways_statistics["primary specialized interfaced"] = 0
    # df.condition_network_pathways_statistics["specialized"] = 0
    # df.condition_network_pathways_statistics["other"] = 0
    
    
    # for(i in 1:length(l.grn_subnetworks.stability_selection)){
    #   
    #   m.grn = l.grn_subnetworks.stability_selection[[i]]
    #   
    #   df.idx.grn <- which(m.grn > 0, arr.ind = TRUE)
    #   df.grn <- data.frame(TF = rownames(m.grn)[df.idx.grn[,1]], TG = colnames(m.grn)[df.idx.grn[,2]], stringsAsFactors = FALSE)
    #   names(df.grn) <- c("TF","Target")
    #   
    #   v.enz.i = intersect(df.grn$Target, v.enz)
    #   
    #   n.enz.i = length(intersect(df.grn$Target, v.enz))
    #   
    #   df.pwys.i = subset(df.enzymes.2, df.enzymes.2$Gene.id %in% v.enz.i)
    #   
    #   n.pwys.i = length(unique(df.pwys.i$Pathway.id))
    #   df.condition_network_statistics$Pathways[i] = as.character(n.pwys.i)
    #   df.condition_network_pathways_statistics$Pathways[i] = as.character(n.pwys.i)
    #   
    #   if(n.enz.i > 0){
    #     for(j in 1:length(v.domains)){
    #       id = paste("", v.domains[j], sep = "")
    #       v.enz_dom.j = intersect(df.grn$Target, unlist(l.domain_enz[id]))
    #       
    #       df.pwys.ij = subset(df.enzymes.2, df.enzymes.2$Gene.id %in% v.enz_dom.j)
    #       n.pwys.ij = length(unique(df.pwys.ij$Pathway.id))
    #       
    #       percentage = round(n.pwys.ij / n.pwys.i * 100,0)
    #       df.condition_network_pathways_statistics[i, v.domains[j]] =  paste(n.pwys.ij, " (", percentage, "%)", sep = "")
    #     }
    #   }
    # }
    
    
    #xtable(df.condition_network_statistics)
    #xtable(df.condition_network_enzymes_statistics)
    #xtable(df.condition_network_pathways_statistics)
    
    write.table(df.condition_network_statistics, "df.condition_network_statistics.csv", row.names = F, quote = F, sep = ";")
    write.table(df.condition_network_enzymes_statistics, "df.condition_network_enzymes_statistics.csv", row.names = F, quote = F, sep = ";")
    #write.table(df.condition_network_pathways_statistics, "df.condition_network_pathways_statistics.csv", row.names = F, quote = F, sep = ";")
    
    
    filename <- paste(foldername.results, "df.condition_network_statistics.csv", sep ="")
    write.table(df.condition_network_statistics, filename, sep = ";", row.names = F)
    
    filename <- paste(foldername.results, "df.condition_network_enzymes_statistics.csv", sep ="")
    write.table(df.condition_network_enzymes_statistics, filename, sep = ";", row.names = F)
    
    #filename <- paste(foldername.results, "df.condition_network_pathways_statistics.csv", sep ="")
    #write.table(df.condition_network_pathways_statistics.csv, filename, sep = ";", row.names = F)
    
  }
  
}


#' run algorithm 
#'
#' @param b.load_grn_inference ("yes","no")
#' @param b.load_TFBS_inference "yes","no")
#' @param b.load_treatment_tissue_inference ("yes","no")
#' @param m.foldChange_differentialExpression
#' @param m.pvalue_differentialExpression
#' @param df.experiment_condition_annotation
#' @param tb.condition_treatments
#' @param tb.condition_tissues
#' @param df.transcriptionFactorAnnotation 
#' @param df.geneGroups
#' @param tb.geneGroups
#' @param v.geneGroups
#' @param l.geneGroups
#' @param n.cpus
#' @param seed (default = 1234)
#' @param importance.measure  (default = "impurity")
#' @param n.trees (default  = 1000)
#' @param n.lead_method_expression_shuffling (default  = 3)
#' @param nbootstrap (default = =100)
#' @param nstepsLARS (default = 5)
#' @param th.lead_grn_method (default = 0.95)
#' @param th.support_grn_methods (default = 0.95)
#' @param n.grnSupport (default = 1)
#' @param file.TF_to_Motif_IDs (default = "data/TF_to_Motif_IDs.txt")
#' @param file.TFBS_motifs (default = "data/Transcription_factor_weight_matrix_Arabidopsis_thaliana.txt")
#' @param file.promoterSeq (default = "data/TAIR10_upstream_1000_20101104.txt")
#' @param file.geneSeq (default = "data/TAIR10_seq_20110103_representative_gene_model_updated.txt")
#' @param th.pre_tss (default= 1000)
#' @param th.post_tss (default = 200)
#' @param genome_nucleotide_distribution ACGT distribution (default = c(0.3253439, 0.1746561, 0.1746561, 0.3253439),
#' @param th.pval.known_motifs (default = 0.05)
#' @param th.diffexp (default = 0.05)
#' @param th.pval.treatment (default = 0.05) 
#' @param th.pval.tissue (default= 0.05)
#' @param th.min.samples (default = 1) 
#' @param s.multipleTestCorrection (default = "none")
#' @param th.min_number_targets (default = 2)
#' @param th.min_number_MR_targets (default = 2)
#' @param th.pval_masterRegulator (default = 0.05)
#' @param foldername.tmp (default = "tmp/") 
#' @param foldername.results = (default = "results/")

#' @keywords 
#' @export
#' @examples
#' 
#' 
#' setwd(...) # set to dataset
#' 
#' 
#' l.data  =  load_datasets(filename.genes = "data/genes.txt",
#'                          filename.experiment_series_ids = "data/experiment_series_ids.txt",
#'                          filename.foldChange_differentialExpression = "data/m.foldChange_differentialExpression.txt",
#'                          filename.pvalue_differentialExpression =	"data/m.pvalue_differentialExpression.txt",
#'                          filename.experiment_condition_tissue_annotation =	"data/df.experiment_condition_annotation.txt",
#'                          filename.transcriptionfactor_annotation = "data/df.transcriptionFactorAnnotation.txt", 
#'                          filename.geneGroups = "data/df.enzymes_w_metabolic_domains.txt")
#' 
#' 
#' 
#' l.results = run_MERIT(b.load_grn_inference = "yes",
#'                       b.load_TFBS_inference = "yes",
#'                       b.load_treatment_tissue_inference = "yes",
#'                       m.foldChange_differentialExpression=l.data$m.foldChange_differentialExpression,
#'                       m.pvalue_differentialExpression=l.data$m.pvalue_differentialExpression,
#'                       df.experiment_condition_annotation=l.data$df.experiment_condition_annotation,
#'                       tb.condition_treatments=l.data$tb.condition_treatments,
#'                       tb.condition_tissues=l.data$tb.condition_tissues,
#'                       df.transcriptionFactorAnnotation=l.data$df.transcriptionFactorAnnotation, 
#'                       df.geneGroups=l.data$df.geneGroups,
#'                       tb.geneGroups=l.data$tb.geneGroups,
#'                       v.geneGroups=l.data$v.geneGroups,
#'                       l.geneGroups=l.data$l.geneGroups, 
#'                       n.cpus = 3,
#'                       seed=1234,
#'                       importance.measure="impurity",
#'                       n.trees=1000,
#'                       n.lead_method_expression_shuffling = 3,
#'                       n.bootstrap=100,
#'                       n.stepsLARS=5,
#'                       th.lead_grn_method = 0.95,
#'                       th.support_grn_methods = 0.95,
#'                       n.grnSupport = 1,
#'                       file.TF_to_Motif_IDs = "data/TF_to_Motif_IDs.txt",
#'                       file.TFBS_motifs = "data/Transcription_factor_weight_matrix_Arabidopsis_thaliana.txt",
#'                       file.promoterSeq = "data/TAIR10_upstream_1000_20101104.txt",
#'                       file.geneSeq = "data/TAIR10_seq_20110103_representative_gene_model_updated.txt",
#'                       th.pre_tss = 1000,
#'                       th.post_tss = 200,
#'                       genome_nucleotide_distribution = c(0.3253439, 0.1746561, 0.1746561, 0.3253439),
#'                       th.pval.known_motifs = 0.05,
#'                       th.diffexp = 0.05,
#'                       th.pval.treatment = 0.05, 
#'                       th.pval.tissue = 0.05,
#'                       th.min.samples = 1, 
#'                       s.multipleTestCorrection = "none",
#'                       th.min_number_targets = 2,
#'                       th.min_number_MR_targets = 2,
#'                       th.pval_masterRegulator = 0.05, 
#'                       foldername.tmp = "tmp/", 
#'                       foldername.results = "results/")
#' 
#' 
#' l.res.grn = l.results$l.res.grn
#' l.res.grn_tfbs = l.results$l.res.grn_tfbs
#' l.res.link_annotation = l.results$l.res.link_annotation
#' l.res.MR_hierarchy = l.results$l.res.MR_hierarchy
#' 
#' 
#' format_results(l.grn_subnetworks = l.res.link_annotation$l.grn_subnetworks,
#'                tb.condition_tissue_differentialExpression = l.res.link_annotation$tb.condition_tissue_differentialExpression,
#'                l.Hierarchy=l.res.MR_hierarchy$l.Hierarchy, 
#'                l.Hierarchy_tfs_per_tier=l.res.MR_hierarchy$l.Hierarchy_tfs_per_tier,
#'                l.Hierarchy_nb_tfs_per_tier=l.res.MR_hierarchy$l.Hierarchy_nb_tfs_per_tier,
#'                l.df.masterRegulatorHierarchy=l.res.MR_hierarchy$l.df.masterRegulatorHierarchy,
#'                v.number_tiers=l.res.MR_hierarchy$v.number_tiers,
#'                m.MR_vs_conditions = l.res.MR_hierarchy$m.MR_vs_conditions,  # A) TFs versus Conditions (Matrix plot) P(TF,C)
#'                l.MR_vs_geneGroups_given_condition = l.res.MR_hierarchy$l.MR_vs_geneGroups_given_condition,  # B) per condition - TFs versus Domains (P(TF,D|C)) => also cumulative plot 
#'                number_of_conditions_per_master_regulator=l.res.MR_hierarchy$number_of_conditions_per_master_regulator,
#'                tb.condition_treatments=l.data$tb.condition_treatments,
#'                tb.condition_tissues=l.data$tb.condition_tissues,
#'                df.transcriptionFactorAnnotation=l.data$df.transcriptionFactorAnnotation, 
#'                df.geneGroups=l.data$df.geneGroups,
#'                tb.geneGroups=l.data$tb.geneGroups,
#'                v.geneGroups=l.data$v.geneGroups,
#'                l.geneGroups=l.data$l.geneGroups,
#'                th.pval = 0.05,
#'                foldername.results = "results/")
#'                
run_MERIT <- function(b.load_grn_inference = "yes",
                      b.load_TFBS_inference = "yes",
                      b.load_treatment_tissue_inference = "yes",
                      m.foldChange_differentialExpression=m.foldChange_differentialExpression,
                      m.pvalue_differentialExpression=m.pvalue_differentialExpression,
                      df.experiment_condition_annotation=df.experiment_condition_annotation,
                      tb.condition_treatments=tb.condition_treatments,
                      tb.condition_tissues=tb.condition_tissues,
                      df.transcriptionFactorAnnotation=df.transcriptionFactorAnnotation, 
                      df.geneGroups=df.geneGroups,
                      tb.geneGroups=tb.geneGroups,
                      v.geneGroups=v.geneGroups,
                      l.geneGroups=l.geneGroups, 
                      n.cpus = 3,
                      seed=1234,
                      importance.measure="impurity",
                      n.trees=1000,
                      n.lead_method_expression_shuffling = 3,
                      n.bootstrap=100,
                      n.stepsLARS=5,
                      th.lead_grn_method = 0.95,
                      th.support_grn_methods = 0.95,
                      n.grnSupport = 1,
                      file.TF_to_Motif_IDs = "data/TF_to_Motif_IDs.txt",
                      file.TFBS_motifs = "data/Transcription_factor_weight_matrix_Arabidopsis_thaliana.txt",
                      file.promoterSeq = "data/TAIR10_upstream_1000_20101104.txt",
                      file.geneSeq = "data/TAIR10_seq_20110103_representative_gene_model_updated.txt",
                      th.pre_tss = 1000,
                      th.post_tss = 200,
                      genome_nucleotide_distribution = c(0.3253439, 0.1746561, 0.1746561, 0.3253439),
                      th.pval.known_motifs = 0.05,
                      th.diffexp = 0.05,
                      th.pval.treatment = 0.05, 
                      th.pval.tissue = 0.05,
                      th.min.samples = 1, 
                      s.multipleTestCorrection = "none",
                      th.min_number_targets = 2,
                      th.min_number_MR_targets = 2,
                      th.pval_masterRegulator = 0.05, 
                      foldername.tmp = "tmp/", 
                      foldername.results = "results/"){
  
  
 
  
  list.of.packages <- c("Biostrings", "TFBSTools")#,"seqLogo", "PWMEnrich", "BSgenome.Athaliana.TAIR.TAIR9")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)){
    source("https://bioconductor.org/biocLite.R")
    biocLite(new.packages)
  } 
  
  if(!"Biostrings" %in% installed.packages()){
    stop("Error: Biostrings package could not be installed from Bioconductor")
  }
  require(Biostrings)
  
  if(!"TFBSTools" %in% installed.packages()){
    stop("Error: TFBSTools package could not be installed from Bioconductor")
  }
  require(TFBSTools)

  ###
  
  require(reshape2)
  require(lars)
  require(foreach)
  require(doParallel)
  require(pheatmap)
  require(ggplot2)
  require(igraph)
  require(seqinr)
  require(networkD3)
  require(taRifx)
  require(parmigene)
  
  
  if(!file.exists(foldername.tmp)){
    dir.create(foldername.tmp)
  }
  
  if(!file.exists(foldername.results)){
    dir.create(foldername.results)
  }
  
  message("--------------------------------------------------")
  message("Step 1 - Gene regulatory network inference using ensemble regression with Monte Carlo based threshold selection")
  
  if(b.load_grn_inference == "no"){
    
    l.res.grn = compute_ensemble_regression_with_montecarlo_based_stability_selection(m.foldChange_differentialExpression=m.foldChange_differentialExpression,
                                                                                      df.transcriptionFactorAnnotation=df.transcriptionFactorAnnotation,
                                                                                      df.geneGroups=df.geneGroups,
                                                                                      seed=seed,
                                                                                      importance.measure=importance.measure,
                                                                                      n.trees=n.trees,
                                                                                      n.lead_method_expression_shuffling = n.lead_method_expression_shuffling,
                                                                                      n.bootstrap=n.bootstrap,
                                                                                      n.stepsLARS=n.stepsLARS,
                                                                                      n.cpus=n.cpus)
    
  }
  
  l.res.grn = load_lead_support_grn(df.transcriptionFactorAnnotation=df.transcriptionFactorAnnotation,
                                    df.geneGroups = df.geneGroups, 
                                    th.lead_grn_method = th.lead_grn_method,
                                    n.lead_method_expression_shuffling = n.lead_method_expression_shuffling,
                                    th.support_grn_methods = th.support_grn_methods,
                                    n.grnSupport = n.grnSupport)            
  
  

  message("")
  message("--------------------------------------------------")
  message("Step 2 - Transcription factor direct target promoter binding based filtering of gene regulatory link predictions")
  
  l.res.grn_tfbs = transcriptionFactorBindingInference(m.grn = l.res.grn$m.lead_support.grn , 
                                                       file.TF_to_Motif_IDs = file.TF_to_Motif_IDs,
                                                       file.TFBS_motifs = file.TFBS_motifs,
                                                       file.promoterSeq = file.promoterSeq,
                                                       file.geneSeq = file.geneSeq,
                                                       th.pre_tss = th.pre_tss,
                                                       th.post_tss = th.post_tss,
                                                       genome_nucleotide_distribution = genome_nucleotide_distribution,
                                                       th.pval.known_motifs=th.pval.known_motifs,
                                                       n.cpus = n.cpus,
                                                       b.load = b.load_TFBS_inference)
  
  message("")
  message("--------------------------------------------------")
  message("Step 3 - Context specific annotation and filtering of gene regulatory link predictions")
  
  l.res.link_annotation = annotate_links_with_treatments_and_tissues(m.lead_support_w_motif.grn=l.res.grn_tfbs$m.lead_support_w_motif.grn, 
                                                                     m.pvalue_differentialExpression=m.pvalue_differentialExpression,
                                                                     df.experiment_condition_annotation=df.experiment_condition_annotation,
                                                                     tb.condition_treatments=tb.condition_treatments,
                                                                     tb.condition_tissues=tb.condition_tissues,
                                                                     th.diffexp = th.diffexp,
                                                                     th.pval.treatment = th.pval.treatment, 
                                                                     th.pval.tissue = th.pval.tissue,
                                                                     th.min.samples = th.min.samples, 
                                                                     s.multipleTestCorrection = "none",
                                                                     b.load = b.load_treatment_tissue_inference)
  
  message("")
  message("--------------------------------------------------")
  message("Step 4 - Master regulator hierarchy inference")
  message("")
  
  l.res.MR_hierarchy = do_master_regulator_hierarchy_inference(m.grn = l.res.link_annotation$m.grn,
                                                               l.grn_subnetworks = l.res.link_annotation$l.grn_subnetworks, 
                                                               df.transcriptionFactorAnnotation = df.transcriptionFactorAnnotation,
                                                               df.geneGroups,
                                                               tb.geneGroups,
                                                               v.geneGroups,
                                                               l.geneGroups,
                                                               th.min_number_targets = th.min_number_targets,
                                                               th.min_number_MR_targets = th.min_number_MR_targets,
                                                               th.pval = th.pval_masterRegulator, 
                                                               foldername.results = foldername.results)
  

  message("--------------------------------------------------")
  return(list(l.res.grn = l.res.grn, 
              l.res.grn_tfbs = l.res.grn_tfbs, 
              l.res.link_annotation = l.res.link_annotation,
              l.res.MR_hierarchy = l.res.MR_hierarchy))
  
}                       










#### OBSOLETE ####
run_grn_stability_selection <- function(m.lead_suppport_w_motif.grn=m.lead_suppport_w_motif.grn, 
                                        m.foldChange_differentialExpression=m.foldChange_differentialExpression,
                                        m.pvalue_differentialExpression=m.pvalue_differentialExpression,
                                        df.experiment_condition_annotation=df.annotation,
                                        tb.condition_treatments=tb.condition_treatments,
                                        tb.condition_tissues=tb.condition_tissues,
                                        v.conditionGroups=v.conditionGroups, 
                                        v.tissueGroups=v.tissueGroups,
                                        v.th.diffexp.stabilitySelection = c(0.01, 0.05, 0.065) ,
                                        th.pval.treatment = 0.05, 
                                        th.pval.tissue = 0.05,
                                        th.min.samples = 1, 
                                        s.multipleTestCorrection = "none"
                                        
                                        ){
  
  #### LINK CONDITION STABILITY ANALYSIS - (variation on the differential expression thresholds)
  m.grn = m.lead_suppport_w_motif.grn
  m.grn[m.grn > 0] = 1
  m.de = m.pvalue_differentialExpression

  if (!file.exists("tmp/grn_stability_selection/")){
    dir.create("tmp/grn_stability_selection/")
  }
  
  strt<-Sys.time() 
  cl<-makeCluster(min(length(v.th.diffexp.stabilitySelection), n.cpus))
  registerDoParallel(cl)
  l.res <-  foreach(p = 1:length(v.th.diffexp.stabilitySelection)) %dopar% { 
    
    th.diffexp.p <- v.th.diffexp.stabilitySelection[p]
   
    subDir <- paste("tmp/grn_stability_selection/", th.diffexp.p, sep ="")
    
    if (!file.exists(subDir)){
      dir.create(subDir)
    }
    
    m.de.bin <- m.de
    m.de.bin[m.de.bin <= th.diffexp.p] <- 10
    m.de.bin[m.de.bin <= 1] <- 0
    m.de.bin[m.de.bin > 0] <- 1
    
    strt<-Sys.time()
    l.res <- perform_treatment_and_tissue_filtering(m.grn = m.grn, 
                                                    m.de.bin = m.de.bin, 
                                                    v.conditionGroups = v.conditionGroups, 
                                                    v.tissueGroups = v.tissueGroups,
                                                    df.annotation = df.annotation, 
                                                    tb.condition_treatments=tb.condition_treatments,
                                                    tb.condition_tissues=tb.condition_tissues,
                                                    th.pval.treatment = th.pval.treatment, 
                                                    th.pval.tissue = th.pval.tissue,
                                                    th.min.samples = th.min.samples, 
                                                    s.multipleTestCorrection = s.multipleTestCorrection)
    
    
    l.treatments_and_tissues <- l.res$l.treatments_and_tissues
    m.rf_w_treatments <- l.res$m.rf_w_treatments
    l.regulatoryNetwork_treatments_and_tissues <- l.res$l.regulatoryNetwork_treatments_and_tissues
    m.gn_condition_tissue_differentialExpression <- l.res$m.gn_condition_tissue_differentialExpression
    tb.condition_tissue_differentialExpression = l.res$tb.condition_tissue_differentialExpression
    l.grn_subnetworks <- l.res$l.grn_subnetworks
    
    filename <- paste(subDir,"/l.treatments_and_tissues.rds", sep ="")
    saveRDS(l.treatments_and_tissues, filename)
    
    filename <- paste(subDir,"/m.rf_w_treatments.rds", sep ="")
    saveRDS(m.rf_w_treatments, filename)
    
    filename <- paste(subDir,"/l.regulatoryNetwork_treatments_and_tissues.rds", sep ="")
    saveRDS(l.regulatoryNetwork_treatments_and_tissues, filename)
    
    filename <- paste(subDir,"/m.gn_condition_tissue_differentialExpression.rds", sep ="")
    saveRDS(m.gn_condition_tissue_differentialExpression, filename)
    
    filename <- paste(subDir,"/l.grn_subnetworks.rds", sep ="")
    saveRDS(l.grn_subnetworks, filename)
    
    filename <- paste(subDir,"/tb.condition_tissue_differentialExpression.rds", sep ="")
    saveRDS(tb.condition_tissue_differentialExpression, filename)
    
    print(Sys.time()-strt)
    
    rm(m.de.bin)
    rm(m.rf_w_treatments)
    rm(l.res)
    rm(l.treatments_and_tissues)
    rm(l.regulatoryNetwork_treatments_and_tissues)
    rm(m.gn_condition_tissue_differentialExpression)
    rm(tb.condition_tissue_differentialExpression)
    rm(l.grn_subnetworks)
    
    return(0)
  }
  stopCluster(cl)
  print(Sys.time()-strt)
  
}


generate_ensemble_grn_based_on_stability_selection <- function(df.transcriptionFactorAnnotation=df.transcriptionFactorAnnotation,
                                                               df.geneGroups=df.geneGroups,
                                                               v.th.diffexp.stabilitySelection = c(0.01, 0.05, 0.1) , 
                                                               th.diffexp_lead = 0.05, 
                                                               th.prob = 0.55){
  
  
  df.transcriptionFactorAnnotation = subset(df.transcriptionFactorAnnotation, df.transcriptionFactorAnnotation$with_geneExpression == "yes")
  df.geneGroups = subset(df.geneGroups, df.geneGroups$with_geneExpression == "yes")
  v.tfs = unique(df.transcriptionFactorAnnotation$TF_ID)
  v.genes =  unique(c(v.tfs, rownames(df.geneGroups)))
  
  if(!th.diffexp_lead %in% v.th.diffexp.stabilitySelection){
    stop("lead diff exp.  not in selection")
  }

  # lead 
  filename.parameter <- th.diffexp_lead
  subDir <- paste("tmp/grn_stability_selection/", filename.parameter, sep ="")
  
  filename <- paste(subDir,"/l.grn_subnetworks.rds", sep ="")
  l.grn_subnetworks = readRDS(filename)
  v.conditions.selected = names(l.grn_subnetworks)
  
  l.grn_subnetworks.stability_selection = vector(mode = "list", length = length(v.conditions.selected))
  names(l.grn_subnetworks.stability_selection) = v.conditions.selected
  for(c in 1:length(v.conditions.selected)){
    l.grn_subnetworks.stability_selection[[c]] = matrix(0, nrow = length(v.tfs), ncol = length(v.genes), dimnames = list(v.tfs, v.genes))
  }
  
  n.max.grn_subnetworks.stability_selection = numeric(length(v.conditions.selected))
  names(n.max.grn_subnetworks.stability_selection) = v.conditions.selected
  
  # support 
  pb <- txtProgressBar(min = 0, max = length(v.th.diffexp.stabilitySelection), style = 3)
  for(p in 1:length(v.th.diffexp.stabilitySelection)){
    setTxtProgressBar(pb, p)
    filename.parameter <- v.th.diffexp.stabilitySelection[p]
    subDir <- paste("tmp/grn_stability_selection/", filename.parameter, sep ="")
    filename <- paste(subDir,"/m.rf_w_treatments.rds", sep ="")
    m.rf_w_treatments = readRDS(filename)
    
    filename <- paste(subDir,"/l.grn_subnetworks.rds", sep ="")
    l.grn_subnetworks = readRDS(filename)
    
    for(c in 1:length(v.conditions.selected)){
      if(v.conditions.selected[c] %in% names(l.grn_subnetworks)){
        m.grn_subnetwork = l.grn_subnetworks[v.conditions.selected[c]][[1]]
        m.grn_subnetwork[m.grn_subnetwork > 0] = 1
        l.grn_subnetworks.stability_selection[[c]][rownames(m.grn_subnetwork), colnames(m.grn_subnetwork)] = l.grn_subnetworks.stability_selection[[c]][rownames(m.grn_subnetwork), colnames(m.grn_subnetwork)] + m.grn_subnetwork
        n.max.grn_subnetworks.stability_selection[c] = n.max.grn_subnetworks.stability_selection[c] + 1
      }
    }
  }
  close(pb)
  
  filename <- "tmp/grn_stability_selection/l.grn_subnetworks.stability_selection.rds"
  saveRDS(l.grn_subnetworks.stability_selection, filename)

  filename <- "tmp/grn_stability_selection/n.max.grn_subnetworks.stability_selection.rds"
  saveRDS(n.max.grn_subnetworks.stability_selection, filename)

  
  #### 
  
  
  filename.parameter <- th.diffexp_lead
  subDir <- paste("tmp/grn_stability_selection/", filename.parameter, sep ="")
  filename <- paste(subDir,"/m.rf_w_treatments.rds", sep ="")
  m.rf_w_treatments = readRDS(filename)
  filename <- paste(subDir,"/l.grn_subnetworks.rds", sep ="")
  l.grn_subnetworks = readRDS(filename)
  
  th.stabilitySelection <- round(length(v.th.diffexp.stabilitySelection) * th.prob, 0)
  for(c in 1:length(v.conditions.selected)){
    if(sum(l.grn_subnetworks.stability_selection[[c]]) > 0){
      l.grn_subnetworks.stability_selection[[c]][l.grn_subnetworks.stability_selection[[c]] < th.stabilitySelection] = 0
      l.grn_subnetworks.stability_selection[[c]][l.grn_subnetworks.stability_selection[[c]] >= th.stabilitySelection] = 1
      l.grn_subnetworks.stability_selection[[c]] = l.grn_subnetworks.stability_selection[[c]][rownames(l.grn_subnetworks[[c]]),] * l.grn_subnetworks[[c]]
    }
  }
  
  
  #filename <- paste(subDir.selection,"/m.rf_w_treatments.rds", sep ="")
  m.rf_w_treatments.stability_selection = matrix(0, 
                                                 nrow = dim(l.grn_subnetworks.stability_selection[[1]])[1],
                                                 ncol = dim(l.grn_subnetworks.stability_selection[[1]])[2],
                                                 dimnames = list(rownames(l.grn_subnetworks.stability_selection[[1]]), 
                                                                 colnames(l.grn_subnetworks.stability_selection[[1]])))
  
  #m.rf_w_treatments.stability_selection[m.rf_w_treatments.stability_selection > 0] = 1
  for(c in 1:length(v.conditions.selected)){
    if(sum(l.grn_subnetworks.stability_selection[[c]]) > 0){
      m.rf_w_treatments.stability_selection = m.rf_w_treatments.stability_selection + l.grn_subnetworks.stability_selection[[c]]
    }
  }
  
  m.rf_w_treatments.stability_selection = m.rf_w_treatments.stability_selection * m.rf_w_treatments
  m.number_of_conditions_per_link.stability_selection = m.rf_w_treatments.stability_selection
  m.rf_w_treatments.stability_selection[m.rf_w_treatments.stability_selection > 0] = 1
  
  if (!file.exists("results/condition_specific_networks/")){
    dir.create("results/condition_specific_networks/")
  }
  
  for(i in 1:length(l.grn_subnetworks.stability_selection)){
    m.grn <- l.grn_subnetworks.stability_selection[[i]]
    df.idx.grn <- which(m.grn> 0, arr.ind = TRUE)
    df.grn <- data.frame(TF = rownames(m.grn)[df.idx.grn[,1]], TG = colnames(m.grn)[df.idx.grn[,2]], stringsAsFactors = FALSE)
    filename <- paste("results/condition_specific_networks/", names(l.grn_subnetworks.stability_selection)[i], ".csv", sep ="")
    write.csv2(df.grn, filename, row.names = FALSE)
  }
  
  df.grn_subnetworks_stability_selection_statistics =  data.frame(conditions = names(l.grn_subnetworks), n_before = unlist(lapply(l.grn_subnetworks, sum)), 
                                                                  n_stability_selection = unlist(lapply(l.grn_subnetworks.stability_selection, sum)))
  
  
  return(list(m.grn.stability_selection=m.rf_w_treatments.stability_selection, 
              l.grn_subnetworks.stability_selection=l.grn_subnetworks.stability_selection, 
              m.number_of_conditions_per_link.stability_selection=m.number_of_conditions_per_link.stability_selection,
              df.grn_subnetworks_stability_selection_statistics = df.grn_subnetworks_stability_selection_statistics))
  
}














