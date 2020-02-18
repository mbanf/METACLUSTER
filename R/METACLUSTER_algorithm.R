
do_treatment_filtering_single_gene_pair <- function(tb.treatments_gene_pair, tb.condition_treatments, th.pval.treatment = 0.05, th.min.samples = 1, s.multipleTestCorrection = "none"){ 
  
  p.treatments_gene_pair <- rep(1, length(tb.treatments_gene_pair))
  names(p.treatments_gene_pair) <- names(tb.treatments_gene_pair)
  
  for(i in 1:length(tb.treatments_gene_pair)){
    hitInSample <- tb.treatments_gene_pair[i] 
    sampleSize <- sum(tb.treatments_gene_pair)
    hitInPop <- tb.condition_treatments[names(tb.treatments_gene_pair)[i]]
    popSize <- sum(tb.condition_treatments)
    failInPop <- popSize - hitInPop
    fc <- ((hitInSample / sampleSize) / (hitInPop / popSize))
    
    if(fc > 1 & hitInSample >= th.min.samples)
      p.treatments_gene_pair[i] <- phyper(hitInSample, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
  }
  
  p.treatments_gene_pair <- p.adjust(p.treatments_gene_pair, s.multipleTestCorrection)
  i.sets <- which(p.treatments_gene_pair <= th.pval.treatment)
  
  if(length(i.sets) > 0){
    tb.treatments_gene_pair <- tb.treatments_gene_pair[i.sets]
    p.treatments_gene_pair <- p.treatments_gene_pair[i.sets]
  }else{ 
    tb.treatments_gene_pair <- c()
    p.treatments_gene_pair <- c()
  }
  
  return(list(tb.treatments_gene_pair = tb.treatments_gene_pair, p.treatments_gene_pair = p.treatments_gene_pair))
}


evaluate_tissues_per_treatment <- function(tb.tissues_per_treatment,  tb.tissues, th.pval.tissue = 0.05, th.min.samples = 1, s.multipleTestCorrection = "none"){ 
  
  p.treatment <- rep(1, length(tb.tissues_per_treatment))
  names(p.treatment) <- names(tb.tissues_per_treatment)
  
  for(i in 1:length(tb.tissues_per_treatment)){
    
    hitInSample <- tb.tissues_per_treatment[i] 
    sampleSize <- sum(tb.tissues_per_treatment)
    hitInPop <- tb.tissues[names(tb.tissues_per_treatment)[i]]
    popSize <- sum(tb.tissues)
    failInPop <- popSize - hitInPop
    
    fc <- ((hitInSample / sampleSize) / (hitInPop / popSize))
    
    if(fc > 1 & hitInSample >= th.min.samples)
      p.treatment[i] <- phyper(hitInSample, hitInPop, failInPop, sampleSize, lower.tail = FALSE)
    
  }
  
  p.treatment <- p.adjust(p.treatment, s.multipleTestCorrection)
  
  i.sets = which(p.treatment <= th.pval.tissue)
  if(length(i.sets) > 0){
    tb.tissues <- tb.tissues_per_treatment[i.sets]
    p.treatment <- p.treatment[i.sets]
  }else{
    tb.tissues <- c()
    p.treatment <- c()
  }
  
  return(list(tb.tissues=tb.tissues, p.treatment=p.treatment))    
  
}


estimate_cluster_coexpression <- function(m.fc, 
                                          m.de.ternary, 
                                          m.co_diff, 
                                          df.annotation,
                                          df.geneCluster,
                                          tb.condition_treatments,
                                          tb.condition_tissues,
                                          n.gene_pair_samples = 1,
                                          th.min_overlap = 1, 
                                          th.min.samples = 1,
                                          th.gns.min = 2, 
                                          
                                          th.p.cluster = 0.05,
                                          th.rho_prob = 0.05,
                                          th.pval.treatment=0.05,
                                          th.pval.tissue=0.05,
                                          
                                          s.multipleHypthesisCorrection = "none",
                                          seed = 1234,
                                          b.choseSimilarConditionInSimulation = TRUE,
                                          b.signatureEnzymeAnalysis = FALSE,
                                          b.prepare = TRUE, v.rho_random.gene_pairs = NULL,
                                          foldername.results=foldername.results,
                                          heatmap_width = 10, heatmap_height = 10){
  
  
  set.seed(seed)
  
  v.conditionGroups = names(tb.condition_treatments)
  names(v.conditionGroups) = v.conditionGroups
  v.tissueGroups = names(tb.condition_tissues)
  names(v.tissueGroups) = v.tissueGroups
  
  diag(m.co_diff) = 0
  hitInPop_codiff = hitInPop <- length(which(m.co_diff >= th.min_overlap))
  popSize_codiff = popSize <- dim(m.co_diff)[1] * dim(m.co_diff)[2]

  p.prior.codiff <- hitInPop / popSize

  if(b.signatureEnzymeAnalysis){
    p.bg <- length(which(df.geneCluster$Enzyme.classification == "")) / nrow(df.geneCluster)
    p.sig <- length(which(df.geneCluster$Enzyme.classification != "")) / nrow(df.geneCluster)
  }
  
  m.functionality <- matrix(1, nrow = length(v.conditionGroups), ncol = length(v.tissueGroups) + 1, dimnames = list(v.conditionGroups, c(v.tissueGroups, "nonspecific")))
  
  tb.series_ids <- colSums(abs(m.de.ternary))
  v.gns <- rownames(m.de.ternary)
  
  v.gcs <- unique(df.geneCluster$Cluster.ID)
  
  df.cluster_annotations <- c()
  
  s.series <- colSums(abs(m.de.ternary)) 
  n.gns.total <- dim(m.de.ternary)[1]
  
  v.rank.clusters <- numeric(length(v.gcs))
  names(v.rank.clusters) <- v.gcs
  
  v.n.genepairs.clusters <- numeric(length(v.gcs))
  names(v.n.genepairs.clusters) <- v.gcs
  
  v.n.bg.genes.clusters <- numeric(length(v.gcs))
  names(v.n.bg.genes.clusters) <- v.gcs
  
  v.n.sig.genes.clusters <- numeric(length(v.gcs))
  names(v.n.sig.genes.clusters) <- v.gcs
  
  p.sig.genes.clusters <- rep(1,length(v.gcs))
  names(p.sig.genes.clusters) <- v.gcs
  
  # compute gene expression significance cutoff
  if(!b.prepare){
    l.treatments <- list()
    # m.functionality <- matrix(0, nrow = length(tb.tissues), ncol = length(tb.treatments), dimnames = list(names(tb.tissues), names(tb.treatments)))
    p.rho <- ecdf(v.rho_random.gene_pairs)
    th.abs_rho <- quantile(v.rho_random.gene_pairs, 1 - th.rho_prob) 
    
    
    # th.abs_rho <- quantile(v.rho_random.gene_pairs, 1 - th.rho_prob) 
    # message("pcc co-expression threshold at ", th.rho_prob," : ", th.abs_rho)
  }
  
  pb <- txtProgressBar(min = 0, max = length(v.gcs), style = 3)
  
  if(is.null(v.rho_random.gene_pairs))
    v.rho_random.gene_pairs <- numeric()
  
  # estimate gene cluster rank 
  for(i in 1:length(v.gcs)){ # per gene cluster
    
    df.geneCluster.i <- subset(df.geneCluster, df.geneCluster$Cluster.ID == v.gcs[i])
    gns.i <- intersect(df.geneCluster.i$Gene.ID, v.gns) # genes on the microarray chip
    
    setTxtProgressBar(pb, i)
    
    # at least two genes in genomic cluster in microarray
    if(length(gns.i) >= th.gns.min){ # retain all above threshold
      
      if(b.signatureEnzymeAnalysis){
        
        v.n.sig.genes.clusters[i] <-  n.sig <- length(which(df.geneCluster.i$Gene.Name != ""))
        v.n.bg.genes.clusters[i] <- n.bg <- length(which(df.geneCluster.i$Gene.Name == ""))
        n.sig.total <- length(which(df.geneCluster$Gene.Name != ""))
        n.bg.total <- length(which(df.geneCluster$Gene.Name == ""))
        
        p.prior <- n.sig.total / (n.sig.total + n.bg.total)
        p.sig.genes.clusters[i] <- binom.test(n.sig, (n.sig + n.bg), p = p.prior)$p.value
        
      }
      
      # co-differential expression sets (across conditions)
      m.gc.i <- m.co_diff[gns.i, gns.i]
      diag(m.gc.i) <- 0
      m.gc.i[lower.tri(m.gc.i)] <- 0
      i.gc_pairs <- which(m.gc.i >= th.min_overlap, arr.ind = TRUE)
      
      if(nrow(i.gc_pairs) > 0){ 
        
        v.codiff_genes = unique(c(rownames(m.gc.i)[i.gc_pairs[,1]], colnames(m.gc.i)[i.gc_pairs[,2]]))
        
        n.codiff_pairs <-  nrow(i.gc_pairs)
        n.all_pairs <- length(gns.i) * (length(gns.i) - 1) / 2 # in cluster 
        p.cluster.codiff <- binom.test(n.codiff_pairs, n.all_pairs, p = p.prior.codiff, alternative ="greater")$p.value

        hitInSample <- n.codiff_pairs
        sampleSize <- n.all_pairs
        hitInPop = hitInPop_codiff
        popSize = popSize_codiff
        failInPop <- popSize - hitInPop
        fc <- ((hitInSample / sampleSize) / (hitInPop / popSize))
        
        p.cluster.codiff <- phyper(hitInSample, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
        
        if(p.cluster.codiff == 0)
          p.cluster.codiff = 1e-100
        
        # correlation analysis (across conditions)
        # correlation threshold on sub-condition groups
        v.rho.gene_pairs <- numeric(nrow(i.gc_pairs))
        for(j in 1:nrow(i.gc_pairs)){ # per gene pair
          
          g1 <- gns.i[i.gc_pairs[j,1]]
          g2 <- gns.i[i.gc_pairs[j,2]]
          # identify shared similarity => pairwise multiplication should be 1 (1 * 1 or -1 * -1)
          idx.vals <- which((m.de.ternary[g1,] * m.de.ternary[g1,]) == 1)
          #l.sets.gene_pairs[j] <- length(idx.vals)
          expr.g1 <- m.fc[g1, idx.vals]
          expr.g2 <- m.fc[g2, idx.vals]
          # rho <- cor.test(expr.g1,expr.g2) # 
          # if(rho$p.value <= th.rho_prob)
          rho <- (cor(expr.g1, expr.g2)) # pearson correlation 
          # if th.min_overlap smaller - use rho for additional statistical analysis ( distribution) 
          v.rho.gene_pairs[j] <- as.numeric(rho) # 
          
          if(b.prepare){
            ### generate pair specific background correlations
            # simulate same treatment sets different genes rather than 
            v.rho_random.gene_pairs.j <- numeric(n.gene_pair_samples)
            for(k in 1:n.gene_pair_samples){ 
              # sample gene pairs
              g.sim <- sample(dim(m.fc)[1], 2) 
              if(b.choseSimilarConditionInSimulation){
                s.sim <- idx.vals# use original conditions
              }else{ 
                s.sim <- sample(dim(m.fc)[2], length(idx.vals))   # sample conditions
              }
              expr.g1 <- m.fc[g.sim[1], s.sim]
              expr.g2 <- m.fc[g.sim[2], s.sim]
              rho <- cor(expr.g1, expr.g2)
              v.rho_random.gene_pairs.j[k] <- rho
            }
            v.rho_random.gene_pairs <- c(v.rho_random.gene_pairs, v.rho_random.gene_pairs.j)
          }
        }
        
        if(!b.prepare){
          
          # selection of gene coexpression support 
          df.codiff <- data.frame(gn.1 = gns.i[i.gc_pairs[,1]], gn.2 = gns.i[i.gc_pairs[,2]], stringsAsFactors = FALSE)
          df.rho <- data.frame(gn.1 = gns.i[i.gc_pairs[,1]], gn.2 = gns.i[i.gc_pairs[,2]], v.pcc = v.rho.gene_pairs, stringsAsFactors = FALSE)
          # add wilcox significance value - not used currently
          # df.rho["p.wilcox"] <- wilcox.test(v.rho.gene_pairs, v.rho_random.gene_pairs)$p.value
          # how many gene pairs make the threshold?
          i.rho <- which(v.rho.gene_pairs >= th.abs_rho)
          
          if(length(i.rho) > 0){ # initial pairs need to be above threshold 
            
            df.rho <- df.rho[i.rho,]
            v.gns.rho <- unique(c(df.rho$gn.1, df.rho$gn.2))
            
            # sort by the coexpression values
            v.gns.rho <- unique(c(df.rho$gn.1, df.rho$gn.2))
            df.geneCluster.rho <- subset(df.geneCluster.i, df.geneCluster.i$Gene.ID %in% v.gns.rho)
            
            # binomial likelihood of the gene cluster to occur by chance
            p.mu <- 1 -  mean(p.rho(df.rho$v.pcc))
            #p.mu <- th.rho_prob 
            
            n.coex_pairs <-  nrow(df.rho)
            n.all_pairs <- length(gns.i) * (length(gns.i) - 1) / 2
            # n.non_coex_pairs <- n.all_pairs - n.coex_pairs
            
            p.cluster.coex <- binom.test(n.coex_pairs, n.codiff_pairs, p = p.mu, alternative ="greater")$p.value
          
            if(p.cluster.coex == 0)
              p.cluster.coex = 1e-100
            
          }else{
            p.cluster.coex = 1 # no coexpression 
          }
            
        
          if(b.signatureEnzymeAnalysis){

            v.n.sig.genes.clusters[i] <-  n.sig <- length(which(df.geneCluster.rho$Gene.Name != ""))
            v.n.bg.genes.clusters[i] <- n.bg <- length(which(df.geneCluster.rho$Gene.Name == ""))
            n.sig.total <- length(which(df.geneCluster$Gene.Name != ""))
            n.bg.total <- length(which(df.geneCluster$Gene.Name == ""))
            p.prior <- n.sig.total / (n.sig.total + n.bg.total)
            p.sig.genes.clusters[i] <- binom.test(n.sig, (n.sig + n.bg), p = p.prior)$p.value

            # signature enzymes
            p.cluster.sig <- choose(n.bg, n.sig) * (p.sig^n.sig) * (p.bg^n.bg)
          }
          
          # identify most predominant treatments
          l.sets_per_gene <- vector(mode = "list", length = length(colnames(m.de.ternary)))
          names(l.sets_per_gene) <- colnames(m.de.ternary)
          
          for(j in 1:nrow(df.codiff)){ # per gene pair
            g1 <- df.codiff$gn.1[j]
            g2 <- df.codiff$gn.2[j]
            # identify shared similarity => pairwise multiplication should be 1 (1 * 1 or -1 * -1)
            idx.vals <- which((m.de.ternary[g1,] * m.de.ternary[g1,]) == 1)
            # annotate gene pairs to conditions 
            l.sets_per_gene[idx.vals] <- lapply(l.sets_per_gene[idx.vals], function(m) unique(c(m,c(g1,g2)))) 
          }
          
          i.set <- which(unlist(lapply(l.sets_per_gene, function(m) if(is.null(m)){ FALSE} else{ TRUE } ) ) == TRUE)
          tb.fc.sets <- unlist(lapply(l.sets_per_gene, length))
          
        
          # minimum of 2 genes to consider condition
          if(any(tb.fc.sets >= 2)){
          
            pvals.tmts <- rep(1, length(tb.fc.sets))
         
            # significance of an experiment based on the number of genes 
            for(j in 1:length(tb.fc.sets)){
              hitInSample <- tb.fc.sets[j]
              if(hitInSample >= 2){
                sampleSize <- length(v.codiff_genes) # codifferentiated genes
                
                hitInPop <- tb.series_ids[j]
                failInPop <- length(v.gns) - hitInPop # genome backgrounds
                
                # number of enzymes instead of genomes
                pvals.tmts[j] <- phyper(hitInSample, hitInPop, failInPop, sampleSize, lower.tail = FALSE)
              }
            }
            
            # identify significant treatment and tissue sets 
            pvals.tmts <- p.adjust(pvals.tmts, s.multipleHypthesisCorrection)
            i.set <- which(pvals.tmts <= th.p.cluster)
            
            if(length(i.set) > 0){
              
              m.functionality.i = matrix(1, nrow = length(v.conditionGroups), ncol = length(v.tissueGroups) + 1, dimnames = list(v.conditionGroups, c(v.tissueGroups, "non specific")) )
              df.cluster_annotations.i <- c()
              
              tb.fc.sets.significant <- tb.fc.sets[i.set]
              series.set = names(tb.fc.sets.significant)
              
              df.annotation.set = subset(df.annotation, df.annotation$unique_ID %in% series.set)
              
              v.gn.names <- df.geneCluster.i$Gene.Name
              names(v.gn.names) <- df.geneCluster.i$Gene.ID  
              
              
              
              number_of_codiff_expressed_genes = length(v.codiff_genes)
              percentage_of_codiff_expressed_genes =  round(length(v.codiff_genes) / length(gns.i) * 100,1)
              codiff_expressed_genes = paste(v.codiff_genes, collapse = " / ")
              names_of_codiff_expressed_genes = paste(v.gn.names[v.codiff_genes], collapse = " / ")
              
              # coexpressed genes
              if(length(i.rho) > 0){
                number_of_coexpressed_genes = length(unique(c(df.rho$gn.1, df.rho$gn.2)))
                percentage_of_coexpressed_genes =  round(length(unique(c(df.rho$gn.1, df.rho$gn.2))) / length(gns.i) * 100,1)
                coexpressed_genes =paste(unique(c(df.rho$gn.1, df.rho$gn.2)), collapse = " / ")
                names_of_coexpressed_genes = paste(v.gn.names[unique(c(df.rho$gn.1, df.rho$gn.2))], collapse = " / ")
              }else{
                number_of_coexpressed_genes = 0
                percentage_of_coexpressed_genes = 0
                coexpressed_genes = ""
                names_of_coexpressed_genes = ""
              }
              
              
              # TRANSFORM treatments, tissues to CONDITIONS (treatment super groups and tissue super groups)
              tb.treatments.j <- numeric(length(v.conditionGroups))
              names(tb.treatments.j) = v.conditionGroups
              for(j in 1:length(series.set)){
                df.annotation.j = subset(df.annotation, df.annotation$unique_ID == series.set[j])
                condition_treatment_1 = df.annotation.j$condition_treatment_1
                condition_treatment_1 = condition_treatment_1[condition_treatment_1 != ""]
                if(length(condition_treatment_1) > 0){
                  tb.treatments.j[v.conditionGroups[condition_treatment_1]] = tb.treatments.j[v.conditionGroups[condition_treatment_1]] + 1
                }
                condition_treatment_2 = df.annotation.j$condition_treatment_2
                condition_treatment_2 = condition_treatment_2[condition_treatment_2 != ""]
                if(length(condition_treatment_2) > 0){
                  tb.treatments.j[v.conditionGroups[condition_treatment_2]] = tb.treatments.j[v.conditionGroups[condition_treatment_2]] + 1
                }
              }
              tb.treatments.j = tb.treatments.j[tb.treatments.j > 0]
                
              if(length(tb.treatments.j) > 0){
          
                # identify significant treatments
                l.res <- do_treatment_filtering_single_gene_pair(tb.treatments.j, tb.condition_treatments, th.pval.treatment, th.min.samples = th.min.samples, s.multipleTestCorrection = "none")
                
                tb.conditions_treatments.significant = l.res$tb.treatments_gene_pair
                p.conditions_treatments.significant = l.res$p.treatments_gene_pair
                
                # remove all links without 
                if(length(tb.conditions_treatments.significant) > 0){
                  # dominant treatment filter per link
                  for(l in 1:length(tb.conditions_treatments.significant)){ 
                    
                    condition = names(tb.conditions_treatments.significant)[l]
                    
                    if(p.conditions_treatments.significant[l] == 0)
                      p.conditions_treatments.significant[l] = 1e-100
                    
                    # which conditionset 
                    i.set <- which(v.conditionGroups[df.annotation.set$condition_treatment_1] == condition)
                    i.set <- unique(c(i.set, which(v.conditionGroups[df.annotation.set$condition_treatment_2] == condition)))
                    
                    df.annotation.set.l = df.annotation.set[i.set,]
                    
                    tb.tissues_per_treatment = table(df.annotation.set.l$condition_tissue)
                   
                    l.res <- evaluate_tissues_per_treatment(tb.tissues_per_treatment,  tb.condition_tissues, 
                                                                               th.pval.tissue = th.pval.tissue, 
                                                                               th.min.samples = th.min.samples, 
                                                                               s.multipleTestCorrection = "none")
                    
                    
                    tb.tissues_per_treatment = l.res$tb.tissues
                    p.tissues_per_treatment = l.res$p.treatment
                    
                 
                    # if significant tissues for treatment are found
                    if(length(tb.tissues_per_treatment) > 0){
                  
                      for(s in 1:length(tb.tissues_per_treatment)){
                        
                        tissue = names(tb.tissues_per_treatment)[s]
                     
                        if(p.tissues_per_treatment[s] == 0)
                          p.tissues_per_treatment[s] = 1e-100
                        
                        p.cluster <- sumlog(c(p.cluster.codiff, p.cluster.coex, p.conditions_treatments.significant[l], p.tissues_per_treatment[s]))$p
                        
                        if(p.cluster <= th.p.cluster){
                        
                          m.functionality.i[condition,tissue] = p.cluster
                        
                          df.annotation.set.l.s = subset(df.annotation.set.l, df.annotation.set.l$condition_tissue == tissue)
                          
                          enzymes = unique(unlist(l.sets_per_gene[df.annotation.set.l.s$unique_ID]))
                
                       
                          # treatment and tissue annotations with at least th.min of genes 
                          df.cluster_annotations.j <- data.frame(cluster.ID = v.gcs[i],
                                                                 
                                                                 condition = condition,
                                                                 tissue = tissue,
                                                                 
                                                                 p.cluster = p.cluster,
                                                                 p.codiff = p.cluster.codiff, 
                                                                 p.pcc = p.cluster.coex, 
                                                                 p.condition = p.conditions_treatments.significant[l],
                                                                 p.tissue = p.tissues_per_treatment[s],
                                                                                                             
                                                                 condition_and_tissue_specific_genes = paste(enzymes, collapse = " / "),
                                                                 names_condition_and_tissue_specific_genes = paste(v.gn.names[enzymes], collapse = " / "),
                                                                 
                                                                 number_of_genes_in_cluster = length(gns.i),
                                                                 genes_in_cluster = paste(gns.i, collapse = " / "),
                                                                 names_genes_in_cluster = paste(v.gn.names[gns.i], collapse = " / "),
                                                                 
                                                                 number_of_codiff_expressed_genes = number_of_codiff_expressed_genes,
                                                                 percentage_of_codiff_expressed_genes = percentage_of_codiff_expressed_genes,
                                                                 codiff_expressed_genes = codiff_expressed_genes,
                                                                 names_of_codiff_expressed_genes = names_of_codiff_expressed_genes,
                                                                 
                                                                 number_of_coexpressed_genes = number_of_coexpressed_genes,
                                                                 percentage_of_coexpressed_genes = percentage_of_coexpressed_genes,
                                                                 coexpressed_genes = coexpressed_genes,
                                                                 names_of_coexpressed_genes = names_of_coexpressed_genes
                                                                 
                                                                 )
                          
                          df.cluster_annotations.i <- rbind(df.cluster_annotations.i, df.cluster_annotations.j)
                          
                        }
                        
                      }
                      
                        
                    }else{ # add non specific
                      
                      p.cluster <- sumlog(c(p.cluster.codiff, p.cluster.coex, p.conditions_treatments.significant[l]))$p
                      
                      
                      if(p.cluster <= th.p.cluster){
                        
                        m.functionality.i[condition,"non specific"] = p.cluster
             
                        enzymes = unique(unlist(l.sets_per_gene[df.annotation.set.l$unique_ID]))
                        
                        v.gn.names <- df.geneCluster.i$Gene.Name
                        names(v.gn.names) <- df.geneCluster.i$Gene.ID  
                        
                        
                        df.cluster_annotations.j <- data.frame(cluster.ID = v.gcs[i],
                                                               
                                                               condition = condition,
                                                               tissue = "non specific",
                                                               
                                                               p.cluster = p.cluster,
                                                               p.codiff = p.cluster.codiff, 
                                                               p.pcc = p.cluster.coex, 
                                                               p.condition = p.conditions_treatments.significant[l],
                                                               p.tissue = NA,
                                                               
                                                               condition_and_tissue_specific_genes = paste(enzymes, collapse = " / "),
                                                               names_condition_and_tissue_specific_genes = paste(v.gn.names[enzymes], collapse = " / "),
                                                               
                                                               number_of_genes_in_cluster = length(gns.i),
                                                               genes_in_cluster = paste(gns.i, collapse = " / "),
                                                               names_genes_in_cluster = paste(v.gn.names[gns.i], collapse = " / "),
                                                               
                                                               number_of_codiff_expressed_genes = number_of_codiff_expressed_genes,
                                                               percentage_of_codiff_expressed_genes = percentage_of_codiff_expressed_genes,
                                                               codiff_expressed_genes = codiff_expressed_genes,
                                                               names_of_codiff_expressed_genes = names_of_codiff_expressed_genes,
                                                               
                                                               number_of_coexpressed_genes = number_of_coexpressed_genes,
                                                               percentage_of_coexpressed_genes = percentage_of_coexpressed_genes,
                                                               coexpressed_genes = coexpressed_genes,
                                                               names_of_coexpressed_genes = names_of_coexpressed_genes
                                                               
                                                               )
                        
                        df.cluster_annotations.i <- rbind(df.cluster_annotations.i, df.cluster_annotations.j)
                        
                      } 
                      
                    }
                    
                  }
                }
                
              }
              
              if(length(df.cluster_annotations.i) > 0)
                df.cluster_annotations <- rbind(df.cluster_annotations, df.cluster_annotations.i)
              
                        
              if(FALSE){ # if(any(m.functionality.i > 0)){
              
                  df.cluster_annotations <- rbind(df.cluster_annotations, df.cluster_annotations.i)
                  
                  #m.functionality.i <- m.functionality.i[rowSums(m.functionality.i) > 0, colSums(m.functionality.i) > 0]
                  m.tmp = m.functionality.i
                  # m.tmp = round(m.tmp * 100,1)
         
                  m.availability = matrix(-1, nrow = nrow(m.functionality.i), ncol = ncol(m.functionality.i), dimnames = list(rownames(m.functionality.i), colnames(m.functionality.i)))
                  for(x in 1:nrow(df.annotation)){
                    if(df.annotation$condition_treatment_1[x] != ""){
                      m.availability[df.annotation$condition_treatment_1[x], df.annotation$condition_tissue[x]] = 1
                    }
                    if(df.annotation$condition_treatment_2[x] != ""){
                      m.availability[df.annotation$condition_treatment_2[x], df.annotation$condition_tissue[x]] = 1
                    }
                  }
                  
                  for(x in 1:nrow(m.availability)){
                    for(y in 1:ncol(m.availability)){
                      if(colnames(m.availability)[y] != "non specific"){
                        if(m.availability[x,y] == -1){
                          m.tmp[x,y] = -1
                        }  
                      }
                    }
                  }
                  m.tmp[m.tmp == -1] = NA
                
                  
                  m.heatmap <- m.tmp[sort(rownames(m.tmp)),]
                  m.heatmap <- t(m.heatmap) #  m.MR # * m.regulatorActivity.pvalue
                  
                  breaks = seq(0,floor(max(m.heatmap[!is.na(m.heatmap)]) + 1),1)
                  color <- c( "black",  "orange" ,colorRampPalette(c( "orange", "red"))(length(breaks) - 1))
                  
                  # plot pdf 3.5 - 9
                  p = pheatmap(m.heatmap, color = color, breaks = breaks, border_color = "orange",  show_rownames = T, show_colnames = T , cluster_rows = F, cluster_cols = F, treeheight_row = 0, treeheight_col = 0,
                               main = "Colors indicate the rank of the cluster to be active per condition and tissue. \n (Gray tiles indicate condition-tissue combinations not present in the expression dataset)",
                               fontsize = 7, fontsize_row = 9, fontsize_col = 9,
                               silent = T)
                  
                               # main = "Gene cluster condition functionality map \n (colors indicate the rank of the cluster active per condition and tissue)")#fontsize = 1)
                  
                  # dev.off()
              
                  save_pheatmap_pdf(p, paste(foldername.results, "/condition_activity_per_cluster_heatmaps/", v.gcs[i],".pdf", sep = ""),  width=heatmap_width, height=heatmap_height)
              
                 
                  #### enzyme coexpression map 
                  
                  # if(FALSE){
                  #   v.gns_plus_names <- paste(gns.i, " (" ,v.gn.names[gns.i], ")", sep = "")
                  #   
                  #   m.coexpression.i <- matrix(0, nrow = length(v.gns_plus_names), ncol = length(v.gns_plus_names), 
                  #                              dimnames = list(v.gns_plus_names, v.gns_plus_names))
                  #   
                  #   
                  #   for(k in 1:nrow(df.rho)){
                  #     k1 <- paste(df.rho$gn.1[k], " (" ,v.gn.names[df.rho$gn.1[k]], ")", sep = "")
                  #     k2 <- paste(df.rho$gn.2[k], " (" ,v.gn.names[df.rho$gn.2[k]], ")", sep = "")
                  #     m.coexpression.i[k1,k2] <- df.rho$v.pcc[k]
                  #   }
                  # }else{
                  # 
                  #   m.coexpression.i <- matrix(0, nrow = length(gns.i), ncol = length(gns.i), 
                  #                              dimnames = list(gns.i, gns.i))
                  #   
                  #   
                  #   for(k in 1:nrow(df.rho)){
                  #     k1 <- df.rho$gn.1[k]
                  #     k2 <- df.rho$gn.2[k] 
                  #     m.coexpression.i[k1,k2] <- df.rho$v.pcc[k]
                  #   }
                  # }
                  # 
                  # 
                  # m.coexpression.i <- m.coexpression.i + t(m.coexpression.i)
                  # 
                  # p <- pheatmap(m.coexpression.i, cluster_rows = FALSE,cluster_cols = FALSE, legend = TRUE, 
                  #               main = paste("Coexpression map of gene cluster ", v.gcs[i], "\n (colors indicate Pearson's correlation values)", sep = ""))
                  # save_pheatmap_pdf(p, paste(foldername.results, "/coexpression_heatmaps/", v.gcs[i],".pdf", sep = ""),  width=heatmap_width, height=heatmap_height)
                  # 
                  # 
                  # functionality map - pvalue 
                  m.functionality.i[m.functionality.i > 0] <- 1
                  m.functionality <- m.functionality + m.functionality.i # global cluster 
                  
              }
        
        
              
            }
          
          }

        }
      }
    }
  }
  
  close(pb)
  
  if(!b.prepare){
    df.cluster_annotations <- df.cluster_annotations[order(df.cluster_annotations$p.cluster),]
    return(df.cluster_annotations) #, m.functionality = m.functionality))
  }else{
    return(v.rho_random.gene_pairs)
  }
  
}



#' run Metacluster algorithm
#'
#' This function runs the condition specific gene cluster annotation
#' @param m.foldChange_differentialExpression differential expression foldchange matrix - rows are genes, cols are experiments
#' @param m.pvalue_differentialExpression differential expression pvalue matrix - rows are genes, cols are experiments
#' @param df.experiment_condition_annotation experiment condition annotation
#' @param df.geneCluster gene cluster dataset
#' @param tb.condition_treatments table of conditions
#' @param tb.condition_tissues table of tissues
#' @param v.conditionGroups treatment and condition map
#' @param v.tissueGroups tissue maps 
#' @param n.cpus number of cores used (default = 1)
#' @param b.load_codifferentialAnalysis_monteCarloSimulation load codifferential expression data ("yes", "no")
#' @param pvalue_DifferentialExpression pvalue treshold for differential expession (default = 0.05)
#' @param probability_codifferentialExpression_MonteCarloSimulation probability threshold codifferential expression (default = 0.05)
#' @param pvalue_coexpression_distribution pvalue treshold context specific coexpression (default = 0.05)
#' @param pvalue_geneClusterPrediction pvalue gene cluster inference enzyme presence (default = 0.05)
#' @param pvalue_geneClusterConsistency pvalue gene cluster enzyme condition consistency (default = 0.05)
#' @param pvalue_treatment_per_condition pvalue gene pair condition annotation (default = 0.05)
#' @param pvalue_tissue_per_condition pvalue gene pair tissue annotation (default = 0.05)
#' @param th.consistent_condition_presence_percentage percentage of gene cluster enyzmes that are expressed in each condition in order to annotate the condition to the cluster (default = 0.7
#' @param number_codifferentialExpression_MonteCarloSimulations number of codiffernetial expression background monte carlo simulations (default = 1)
#' @param number_conditionSpecificCoexpressionBackgroundGenePairs number of context specific coexpression simulation background gene pairs (default = 50)
#' @param min_number_condition_samples minimum number of condition samples for significance test (default 1)
#' @param foldername.tmp temp file folder name (default = /tmp)
#' @param foldername.results results file folder name (default = /results)
#' @param heatmap_width default = 10
#' @param heatmap_height default = 5
#' @return a list of results
#' @keywords 
#' @export
#' @examples
#' 
#' # install_and_load_libraries()
#' 
#' # set directory to dataset directory, e.g. /User/home/athaliana_schlapfer2017/
#' setwd(...)
#' 
#' 
#' message("load datasets")
#' l.data = load_datasets(input_format = "PCF2017_enzymes_only",
#'                        filename.genes = "data/genes.txt",
#'                        filename.experiment_ids = "data/experiment_ids.txt",
#'                        filename.geneCluster = "data/ath_geneInCluster_3_aracyc.txt-labeled_NoHypoGenes.txt",
#'                        filename.foldChange_differentialExpression = "data/m.foldChange_differentialExpression.txt",
#'                        filename.pvalue_differentialExpression =	"data/m.pvalue_differentialExpression.txt",
#'                        filename.experiment_condition_tissue_annotation ="data/experiment_annotation.txt")
#' 
#' message("run METACLUSTER")
#' df.cluster_annotations = run_METACLUSTER(m.foldChange_differentialExpression = l.data$m.foldChange_differentialExpression,
#'                                          m.pvalue_differentialExpression = l.data$m.pvalue_differentialExpression,
#'                                          df.experiment_condition_annotation = l.data$df.experiment_condition_annotation,
#'                                          df.geneCluster = l.data$df.geneCluster,
#'                                          tb.condition_treatments = l.data$tb.condition_treatments,
#'                                          tb.condition_tissues = l.data$tb.condition_tissues,
#'                                          n.cpus = 3,
#'                                          b.load_codifferentialAnalysis_monteCarloSimulation = "yes",
#'                                          pvalue_DifferentialExpression = 0.05,
#'                                          probability_codifferentialExpression_MonteCarloSimulation = 0.95,
#'                                          pvalue_coexpression_distribution = 0.05,
#'                                          pvalue_geneClusterPrediction = 0.05,
#'                                          pvalue_geneClusterConsistency = 0.05,
#'                                          pvalue_treatment_per_condition = 0.05,
#'                                          pvalue_tissue_per_condition = 0.05,
#'                                          number_codifferentialExpression_MonteCarloSimulations = 1,
#'                                          number_conditionSpecificCoexpressionBackgroundGenePairs = 100,
#'                                          min_number_condition_samples = 1,
#'                                          seed = 1234,
#'                                          heatmap_width = 10,
#'                                          heatmap_height = 5,
#'                                          foldername.results = "results/",
#'                                          foldername.tmp = "tmp/")
#' 
#' 
#' evaluate_and_store_results(df.cluster_annotations=df.cluster_annotations,
#'                            df.experiment_condition_annotation = l.data$df.experiment_condition_annotation,
#'                            tb.condition_treatments = l.data$tb.condition_treatments,
#'                            tb.condition_tissues = l.data$tb.condition_tissues,
#'                            min_number_of_genes = 3,
#'                            heatmap_width = 4, heatmap_height = 7, fontsize = 7, fontsize_row = 10, fontsize_col = 10,
#'                            foldername.results = "results/")

run_METACLUSTER = function(m.foldChange_differentialExpression,
                           m.pvalue_differentialExpression,
                           df.experiment_condition_annotation,
                           df.geneCluster,
                           tb.condition_treatments,
                           tb.condition_tissues,
                           n.cpus = 1,
                           b.load_codifferentialAnalysis_monteCarloSimulation = "yes",
                           pvalue_DifferentialExpression = 0.05,
                           probability_codifferentialExpression_MonteCarloSimulation = 0.95,
                           pvalue_coexpression_distribution = 0.05,
                           pvalue_geneClusterPrediction = 0.05,
                           pvalue_geneClusterConsistency = 0.05,
                           pvalue_treatment_per_condition = 0.05,
                           pvalue_tissue_per_condition = 0.05,
                          
                           number_codifferentialExpression_MonteCarloSimulations = 3,
                           number_conditionSpecificCoexpressionBackgroundGenePairs = 100,
                           min_number_condition_samples = 1,
                           seed = 1234,
                           heatmap_width = 10, heatmap_height = 6,
                           foldername.tmp = "tmp/",
                           foldername.results = "results/"){

  
  list.of.packages <- c("reshape2", "foreach", "doParallel", "pheatmap", "metap")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)){
    install.packages(new.packages)
  } 
  
  library(reshape2)
  library(foreach)
  library(doParallel)
  library(pheatmap)
  library(metap)
  
  if(!file.exists(foldername.tmp)){
    dir.create(foldername.tmp)
  }
  
  if(!file.exists(foldername.results)){
    dir.create(foldername.results)
  }
  
  if(!file.exists(paste(foldername.results, "/condition_activity_per_cluster_heatmaps/", sep = ""))){
    dir.create(paste(foldername.results, "/condition_activity_per_cluster_heatmaps/", sep = ""))
  }
  
  m.fc = m.foldChange_differentialExpression
  m.de = m.pvalue_differentialExpression
  df.annotation = df.experiment_condition_annotation
  
  th.diffexp = pvalue_DifferentialExpression
  th.prob.MC = probability_codifferentialExpression_MonteCarloSimulation
  th.rho_prob = pvalue_coexpression_distribution
  th.p.cluster = pvalue_geneClusterPrediction
  th.p.cluster_consistency = pvalue_geneClusterConsistency
  th.pval.treatment = pvalue_treatment_per_condition
  th.pval.tissue = pvalue_tissue_per_condition
  n.sim = number_codifferentialExpression_MonteCarloSimulations
  n.gene_pair_samples = number_conditionSpecificCoexpressionBackgroundGenePairs
  th.min.samples = min_number_condition_samples
  
  
  
  
  # make ternary - up and down regulation 
  m.fc.ternary <- m.fc
  m.fc.ternary[m.fc.ternary > 0] <- 1
  m.fc.ternary[m.fc.ternary < 0] <- -1
  
  m.de.bin <- m.de
  m.de.bin[m.de.bin <= th.diffexp] <- 10
  m.de.bin[m.de.bin <= 1] <- 0
  m.de.bin[m.de.bin > 0] <- 1
  
  m.de.ternary <- m.fc.ternary * m.de.bin
  
  #####
  
  # shuffle rowwise (per gene)
  if(b.load_codifferentialAnalysis_monteCarloSimulation == "no"){
    
    message("Computing co-differential expression matrix...")
    
    m.co_diff <- matrix(0, nrow = nrow(m.de.ternary), ncol = nrow(m.de.ternary), dimnames = list(rownames(m.de.ternary), rownames(m.de.ternary)))
    pb <- txtProgressBar(min = 0, max = nrow(m.de.ternary), style = 3)
    for(i in 1:nrow(m.de.ternary)){
      setTxtProgressBar(pb, i)
      for(j in i:nrow(m.de.ternary)){
        v.codif <- m.de.ternary[i,] * m.de.ternary[j,]
        m.co_diff[i,j] <- sum(v.codif[v.codif == 1])
      }
    }
    close(pb)
    saveRDS(m.co_diff, paste(foldername.tmp, "m.co_diff.", th.diffexp, ".rds", sep = ""))

    tb.series_ids <- colSums(m.de.bin)
    v.gns <- rownames(m.de.ternary)
    rm(m.co_diff)
    
    message("...finished")
    message("")
    
    
    message(paste("Computing significant minimum level of overlapping conditions based on", n.sim, "runs of Monte Carlo simulations for co-differential expression..."))
    # 
    strt<-Sys.time()
    cl<-makeCluster(min(n.cpus, n.sim))
    registerDoParallel(cl)
    l.res <- foreach(s = 1:n.sim) %dopar% {
      # for(s in 1:n.sim)
      set.seed((seed + 123 * s))
      m.shuffled <- t(apply(m.de.ternary, 1, sample))
      m.co_diff.shuffled <- matrix(0, nrow = nrow(m.de.ternary), ncol = nrow(m.de.ternary), dimnames = list(rownames(m.de.ternary), rownames(m.de.ternary)))

      pb <- txtProgressBar(min = 0, max = nrow(m.shuffled), style = 3)
      for(i in 1:nrow(m.shuffled)){
        setTxtProgressBar(pb, i)
        for(j in i:nrow(m.shuffled)){
          v.codif <- m.shuffled[i,] * m.shuffled[j,]
          m.co_diff.shuffled[i,j] <- sum(v.codif[v.codif == 1])
        }
      }
      close(pb)

      saveRDS(m.co_diff.shuffled, paste(foldername.tmp, "m.co_diff_shuffled.", s, ".", th.diffexp, ".rds", sep = ""))
      rm(m.co_diff.shuffled)
    }
    stopCluster(cl)
    print(Sys.time()-strt)
    
    v.th.likelihood <- numeric(n.sim)
    
    for(s in 1:n.sim){
      m.co_diff.shuffled <- readRDS(paste(foldername.tmp, "m.co_diff_shuffled.", s, ".", th.diffexp, ".rds", sep = ""))
      m.co_diff.shuffled <- m.co_diff.shuffled + t(m.co_diff.shuffled)
      diag(m.co_diff.shuffled) <- 0
      #th.likelihood <- quantile(m.co_diff.shuffled, 0.95)
      #print(th.likelihood)
      #th.likelihood <- quantile(m.co_diff.shuffled, 0.99)
      #print(th.likelihood)
      v.th.likelihood[s] <- quantile(m.co_diff.shuffled, th.prob.MC)
      rm(m.co_diff.shuffled)
    }
    
    print(Sys.time()-strt)
    
    saveRDS(v.th.likelihood, paste(foldername.tmp, "v.th.likelihood_", th.diffexp, "_" ,th.prob.MC, "_",n.sim, ".rds", sep = ""))
    
    message("...finished")
    message("")
  }
  
  if(!file.exists(paste(foldername.tmp, "v.th.likelihood_", th.diffexp, "_" ,th.prob.MC,"_",n.sim, ".rds", sep = ""))){
    stop("Error: Monte Carlo simulation and / or co-differential expression matrix does not exist!")
  }
  
  m.co_diff <- readRDS(paste(foldername.tmp, "m.co_diff.", th.diffexp, ".rds", sep = ""))
  v.th.likelihood <- readRDS(paste(foldername.tmp, "v.th.likelihood_", th.diffexp, "_" ,th.prob.MC, "_",n.sim, ".rds", sep = ""))
  
  th.min_overlap <- round(mean(v.th.likelihood),0)
  
  #####
  
  message("Establishing condition specific co-expression by chance distribution... ")
  
  v.rho_random.gene_pairs <- estimate_cluster_coexpression(m.fc = m.fc, 
                                                           m.de.ternary = m.de.ternary, 
                                                           m.co_diff = m.co_diff,
                                                           df.annotation = df.annotation,
                                                           df.geneCluster = df.geneCluster,
                                                           tb.condition_treatments = tb.condition_treatments,
                                                           tb.condition_tissues = tb.condition_tissues,
                                                           n.gene_pair_samples = n.gene_pair_samples,
                                                           th.min_overlap = th.min_overlap, 
                                                           th.min.samples = th.min.samples,
                                                           th.gns.min = 2,
                                                           th.p.cluster = th.p.cluster,
                                                           th.rho_prob = th.rho_prob,
                                                           th.pval.treatment=th.pval.treatment,
                                                           th.pval.tissue=th.pval.tissue,
                                                           s.multipleHypthesisCorrection = "none",
                                                           seed=seed,
                                                           b.choseSimilarConditionInSimulation = TRUE,
                                                           b.signatureEnzymeAnalysis = FALSE,
                                                           b.prepare = TRUE, v.rho_random.gene_pairs = NULL,
                                                           foldername.results=foldername.results,
                                                           heatmap_width = heatmap_width, heatmap_height = heatmap_height)
  
  message("...finished")
  message("")
  message("Starting gene cluster context specific transcriptional activity analysis...")
  
  
  # s.multipleHypthesisCorrection = "none"
  # b.choseSimilarConditionInSimulation = TRUE
  # b.signatureEnzymeAnalysis = FALSE
  # b.prepare = FALSE

  # estimate genomic cluster correlation - integrate automatic adaption 
  l.results <- estimate_cluster_coexpression(m.fc = m.fc, m.de.ternary = m.de.ternary, m.co_diff = m.co_diff, 
                                             df.annotation = df.annotation,
                                             df.geneCluster = df.geneCluster,
                                             tb.condition_treatments = tb.condition_treatments,
                                             tb.condition_tissues = tb.condition_tissues,
                                             n.gene_pair_samples = n.gene_pair_samples,
                                             th.min_overlap = th.min_overlap, 
                                             th.min.samples = th.min.samples,
                                             th.gns.min = 2,
                                             th.p.cluster = th.p.cluster,
                                             th.rho_prob = th.rho_prob,
                                             th.pval.treatment=th.pval.treatment,
                                             th.pval.tissue=th.pval.tissue,
                                             s.multipleHypthesisCorrection = "none",
                                             seed=seed,
                                             b.choseSimilarConditionInSimulation = TRUE,
                                             b.signatureEnzymeAnalysis = FALSE,
                                             b.prepare = FALSE, v.rho_random.gene_pairs = v.rho_random.gene_pairs,
                                             foldername.results=foldername.results,
                                             heatmap_width = heatmap_width, heatmap_height = heatmap_height)
  
  
  
  message("...finished")
  message("")
  
  return(l.results)
}
