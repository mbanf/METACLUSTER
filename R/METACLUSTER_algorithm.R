
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
  
  if(length(i.sets) > 0)
    tb.treatments_gene_pair <- tb.treatments_gene_pair[i.sets]
  else 
    tb.treatments_gene_pair <- c()
  
  return(tb.treatments_gene_pair)
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
  if(length(i.sets) > 0)
    tb.tissues_per_treatment <- tb.tissues_per_treatment[i.sets]  
  else 
    tb.tissues_per_treatment <- c()
  
  
  return(tb.tissues_per_treatment)    
  
}


estimate_cluster_coexpression <- function(m.fc = m.fc, 
                                          m.de.ternary = m.de.ternary, 
                                          m.co_diff = m.co_diff, 
                                          df.annotation = df.annotation,
                                          df.geneCluster = df.geneCluster,
                                          tb.condition_treatments = tb.condition_treatments,
                                          tb.condition_tissues = tb.condition_tissues,
                                          n.gene_pair_samples = n.gene_pair_samples,
                                          th.min_overlap = th.min_overlap, 
                                          th.min.samples = th.min.samples,
                                          th.gns.min = th.gns.min, 
                                          th.p.cluster = th.p.cluster,
                                          th.rho_prob = th.rho_prob,
                                          th.pval.treatment=th.pval.treatment,
                                          th.pval.tissue=th.pval.tissue,
                                          th.consistent_condition_presence_percentage = th.consistent_condition_presence_percentage,
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
  hitInPop <- length(which(m.co_diff >= th.min_overlap))
  popSize <- dim(m.co_diff)[1] * dim(m.co_diff)[2]
 
  p.prior.codiff <- hitInPop / popSize
  
  if(b.signatureEnzymeAnalysis){
    p.bg <- length(which(df.geneCluster$Enzyme.classification == "")) / nrow(df.geneCluster)
    p.sig <- length(which(df.geneCluster$Enzyme.classification != "")) / nrow(df.geneCluster)
  }
  
  m.functionality <- matrix(0, nrow = length(v.conditionGroups), ncol = length(v.tissueGroups), dimnames = list(v.conditionGroups, v.tissueGroups))
                          
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
    
    # at least three genes in genomic cluster in microarray
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
        
        n.codiff_pairs <-  nrow(i.gc_pairs)
        n.all_pairs <- length(gns.i) * (length(gns.i) - 1) / 2 # in cluster 
        p.cluster.codiff <- binom.test(n.codiff_pairs, n.all_pairs, p = p.prior.codiff, alternative ="greater")$p.value
          
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
            
            for(j in 1:nrow(df.rho)){ # per gene pair
              g1 <- df.rho$gn.1[j]
              g2 <- df.rho$gn.2[j]
              # identify shared similarity => pairwise multiplication should be 1 (1 * 1 or -1 * -1)
              idx.vals <- which((m.de.ternary[g1,] * m.de.ternary[g1,]) == 1)
              # annotate gene pairs to conditions 
              l.sets_per_gene[idx.vals] <- lapply(l.sets_per_gene[idx.vals], function(m) unique(c(m,c(g1,g2)))) 
            }
            
            i.set <- which(unlist(lapply(l.sets_per_gene, function(m) if(is.null(m)){ FALSE} else{ TRUE } ) ) == TRUE)
            # l.sets_per_gene <- lapply(l.sets_per_gene[sets.vals], function(m)  ) 
            #tb.fc.sets <- lapply(l.sets_per_gene[i.set], function(m) unique(c(m,c(g1,g2)))) 
            tb.fc.sets <- unlist(lapply(l.sets_per_gene[i.set], length))
            l.genes_per_sets <- l.sets_per_gene[i.set]
            
            # filter treatment and tissue annotations by min 50 % of coexpressed genes
            # tb.fc.sets <- tb.fc.sets[which(tb.fc.sets >= nrow(df.rho) * th.cond_min)]
            i.set <- which(tb.fc.sets >= th.gns.min)
            
            tb.fc.sets <- tb.fc.sets[i.set]
            l.genes_per_sets <- l.genes_per_sets[i.set]
            
            if(length(tb.fc.sets) > 0){
              
              pvals.tmts <- rep(1, length(tb.fc.sets))
              df.annotation.j <- subset(df.annotation, df.annotation$series_id %in% names(tb.fc.sets)) # requires two sets of curations (tissue and treatment)
              df.annotation.j <- unique(df.annotation.j[,c("series_id", "condition_treatment_1", "condition_treatment_2", "condition_tissue")])
              
              for(j in 1:length(tb.fc.sets)){
                
                hitInSample <- tb.fc.sets[j]
                sampleSize <- length(v.gns.rho)
                
                hitInPop <- tb.series_ids[names(tb.fc.sets)[j]]
                failInPop <- length(v.gns) - hitInPop # genome backgrounds
                
                # number of enzymes instead of genomes
                pvals.tmts[j] <- phyper(hitInSample, hitInPop, failInPop, sampleSize, lower.tail = FALSE)
              }
            
              # identify significant treatment and tissue sets 
              pvals.tmts <- p.adjust(pvals.tmts, s.multipleHypthesisCorrection)
              i.set <- which(pvals.tmts < th.p.cluster)
              
              if(length(i.set) > 0){
                
                tb.fc.sets.significant <- tb.fc.sets[i.set]
                l.genes_per_sets.significant <- l.genes_per_sets[i.set]
           
                p.prior <- s.series[names(tb.fc.sets.significant)] / n.gns.total
                p.conds <- sapply(1:length(tb.fc.sets.significant), 
                                  function(i) binom.test(as.integer(tb.fc.sets.significant[i]), 
                                                         length(v.gns.rho), 
                                                         p = p.prior[i])$p.value) 
                
                # Dimension 2 - all genes represented in each condition
                p.all_gns_in_each_condition <- p.adjust(p.conds, s.multipleHypthesisCorrection)
                
                # conditins to keep 
                i.set <- (which(p.all_gns_in_each_condition < th.p.cluster))
                perc.sets.significant <- tb.fc.sets.significant[i.set] / length(v.gns.rho)
                i.set <- which(perc.sets.significant >= th.consistent_condition_presence_percentage)
                
                if(length(i.set) > 0){
                  
                  tb.fc.sets.significant <- tb.fc.sets.significant[i.set]
                  p.all_gns_in_each_condition <- p.all_gns_in_each_condition[i.set]
                  
                  if(length(i.set) > 1){
                    p.all_gns_in_each_condition <- sumlog(p.all_gns_in_each_condition)$p
                  }
                  
             
                  p.cluster <- sumlog(c(p.cluster.codiff, p.cluster.coex, p.all_gns_in_each_condition))$p
                  
                    if(p.cluster <= th.p.cluster){
                  
                      tb.series.set = table(names(tb.fc.sets.significant))
                      series.set = names(tb.series.set)
                      
                      # TRANSFORM treatments, tissues to CONDITIONS (treatment super groups and tissue super groups)
                      df.annotation.j <- subset(df.annotation.j, df.annotation.j$series_id %in% series.set)
    
                      ####
                      
                      tb.treatments.j <- numeric(length(v.conditionGroups))
                      names(tb.treatments.j) = v.conditionGroups
                      for(j in 1:length(series.set)){
                        df.annotation.ij = subset(df.annotation.j, df.annotation.j$series_id == series.set[j])
                        condition_treatment_1 = df.annotation.ij$condition_treatment_1
                        condition_treatment_1 = condition_treatment_1[condition_treatment_1 != ""]
                        if(length(condition_treatment_1) > 0){
                          tb.treatments.j[v.conditionGroups[condition_treatment_1]] = tb.treatments.j[v.conditionGroups[condition_treatment_1]] + tb.series.set[series.set[j]]
                        }
                        condition_treatment_2 = df.annotation.ij$condition_treatment_2
                        condition_treatment_2 = condition_treatment_2[condition_treatment_2 != ""]
                        if(length(condition_treatment_2) > 0){
                          tb.treatments.j[v.conditionGroups[condition_treatment_2]] = tb.treatments.j[v.conditionGroups[condition_treatment_2]] + tb.series.set[series.set[j]]
                        }
                      }
                      tb.treatments.j = tb.treatments.j[tb.treatments.j > 0]
                     
                      ###
                      
                      if(length(tb.treatments.j) > 0){
                       
                        # identify significant treatments
                        tb.conditions_treatments.significant <- do_treatment_filtering_single_gene_pair(tb.treatments.j, tb.condition_treatments, th.pval.treatment, th.min.samples = th.min.samples, s.multipleTestCorrection = "none")
                        
                        # remove all links without 
                        if(length(tb.conditions_treatments.significant) > 0){
                          
                          ## Step 2 - assign dominant tissues to these conditions 
                          i.set <- which(v.conditionGroups[df.annotation.j$condition_treatment_1] %in% names(tb.conditions_treatments.significant))
                          i.set <- unique(c(i.set, which(v.conditionGroups[df.annotation.j$condition_treatment_2] %in% names(tb.conditions_treatments.significant))))
                          
                          # tb.tissues.j <- table(df.annotation.j$tissue)
                          tb.tissues.j <- table(df.annotation.j$condition_tissue[i.set]) # tissue annotations all treatments 
                          
                          l.conditions_treatment_and_tissue <- vector(mode = "list", length = length(tb.conditions_treatments.significant))
                          names(l.conditions_treatment_and_tissue) <- names(tb.conditions_treatments.significant)
                          
                          # dominant treatment filter per link
                          for(t in 1:length(tb.conditions_treatments.significant)){ 
                            
                            l.conditions_treatment_and_tissue[[t]] <- character()
                            
                            # which conditionset 
                            i.set <- which(v.conditionGroups[df.annotation.j$condition_treatment_1] == names(tb.conditions_treatments.significant)[t])
                            i.set <- unique(c(i.set, which(v.conditionGroups[df.annotation.j$condition_treatment_2] == names(tb.conditions_treatments.significant)[t])))
                      
                            df.annotation.jt = df.annotation.j[i.set,]
                            tissues.t = unique(df.annotation.jt$condition_tissue)
                            
                            tb.tissues_per_treatment = numeric(length(tissues.t))
                            names(tb.tissues_per_treatment) = tissues.t
                            
                            for(k in 1:length(tissues.t)){
                              df.annotation.jtk = subset(df.annotation.jt, df.annotation.jt$condition_tissue == tissues.t[k])
                              tb.tissues_per_treatment[k] = sum(tb.series.set[df.annotation.jtk$series_id])
                            }
                            
                            tb.tissues_per_treatment <- evaluate_tissues_per_treatment(tb.tissues_per_treatment,  tb.condition_tissues, 
                                                                                         th.pval.tissue = th.pval.tissue, 
                                                                                         th.min.samples = th.min.samples, 
                                                                                         s.multipleTestCorrection = "none")
                            

                            # if significant tissues for treatment are found
                            if(length(tb.tissues_per_treatment) > 0){
                              l.conditions_treatment_and_tissue[[t]] <- names(tb.tissues_per_treatment)
                            }
                          }
                          
                          idx.treatment <- which(sapply(l.conditions_treatment_and_tissue, function(m) length(m) > 0) == TRUE)
                          
                          if(length(idx.treatment) > 0){
                            
                          
                            rank.i <- - (log(p.cluster.codiff) + log(p.cluster.coex) + log(p.all_gns_in_each_condition)) # rank only if cluster is selected
                            
                            
                             
                            #v.rank.clusters[i] <- rank.i
                            # v.n.genepairs.clusters[i] <- nrow(df.rho) # totoal number of gene cluster enzymes
                            
                            
                            # list with individual (multiple) tissue elements
                            l.conditions_treatment_and_tissue <- l.conditions_treatment_and_tissue[idx.treatment]
                          
                            m.functionality.i <- matrix(0, nrow = length(v.conditionGroups), ncol = length(v.tissueGroups), 
                                                        dimnames = list(v.conditionGroups, v.tissueGroups))
                         
                            df.cluster_annotations.i <- c()
                            
                            # treatment and tissue annotations with at least th.min of genes 
                            for(t in 1:length(l.conditions_treatment_and_tissue)){
                              
                              t.set <- unique(c(which(v.conditionGroups[df.annotation.j$condition_treatment_1] == names(l.conditions_treatment_and_tissue)[t]),
                                                which(v.conditionGroups[df.annotation.j$condition_treatment_2] == names(l.conditions_treatment_and_tissue)[t])))
                              
                              for(k in 1:length(l.conditions_treatment_and_tissue[[t]])){
                                
                                k.set <- which(df.annotation.j$condition_tissue == l.conditions_treatment_and_tissue[[t]][k])
                                                
                                tk.set <- intersect(t.set, k.set) # treatment and tissue
                                series.id.tk <- df.annotation.j$series_id[tk.set]
                                
                                # identify all enzymes according 
                                v.gns.tk <- unique(unlist(l.genes_per_sets.significant[series.id.tk]))
                                m.functionality.i[names(l.conditions_treatment_and_tissue)[t], l.conditions_treatment_and_tissue[[t]][k]] <- 1 #length(v.gns.tk) # round(length(v.gns.tk) / length(gns.i) * 100,1)
                              
                                #df.cluster_annotations.j <- data.frame()
                                
                                v.gn.names <- df.geneCluster.i$Gene.Name
                                names(v.gn.names) <- df.geneCluster.i$Gene.ID  
                                
                                
                                # treatment and tissue annotations with at least th.min of genes 
                                df.cluster_annotations.j <- data.frame(cluster.ID = v.gcs[i],
                                                                       rank = rank.i,
                                                                       number_of_enzymes_total = length(gns.i),
                              
                                                                       number_of_coexpressed_enzymes = length(unique(c(df.rho$gn.1, df.rho$gn.2))),
                                                                       percentage_of_coexpressed_enzymes =  round(length(unique(c(df.rho$gn.1, df.rho$gn.2))) / length(gns.i) * 100,1),
                                
                                ###
                                                                        
                                                                        condition_treatment =  names(l.conditions_treatment_and_tissue)[t],
                                                                        condition_tissue =  l.conditions_treatment_and_tissue[[t]][k],
                                                                        
                                                                        number_of_coexpressed_enzymes_in_condition = length(v.gns.tk),
                                                                        percentage_of_coexpressed_enzymes_in_condition = round(length(v.gns.tk) / length(gns.i) * 100,1),
                                                                        ids_of_coexpressed_enzymes_in_condition =paste(v.gns.tk, collapse = " / "),
                                                                        
                                                                       
                                                                        names_of_enzymes_in_condition = paste(v.gn.names[v.gns.tk], collapse = " / "))
                                                                        
                                df.cluster_annotations.i <- rbind(df.cluster_annotations.i, df.cluster_annotations.j)
                                
                              }
                            }
                          
                            
                            df.cluster_annotations <- rbind(df.cluster_annotations, df.cluster_annotations.i)
                            
                            
                            #m.functionality.i <- m.functionality.i[rowSums(m.functionality.i) > 0, colSums(m.functionality.i) > 0]
                            m.heatmap <- t(m.functionality.i) #  m.MR # * m.regulatorActivity.pvalue
                            # colnames(m.heatmap)[21] <- "pathogen (pythopthera, botrytis)"
                            # colnames(m.heatmap)[7] <- "hormone (BRs, ethylene)"
                            # 
                            #m.heatmap <- m.functionality.i
                            #m.heatmap <- m.heatmap[which(rowSums(m.heatmap) > 0), which(colSums(m.heatmap) > 0)]
                            p <- pheatmap((m.heatmap), cluster_rows = FALSE,cluster_cols = FALSE, color = c("black", "white"), legend = FALSE,
                                          main = paste("Condition activity map of gene cluster ", v.gcs[i], "\n (white indicates activity in condition)", sep = ""))
                            save_pheatmap_pdf(p, paste(foldername.results, "/condition_activity_per_cluster_heatmaps/", v.gcs[i],".pdf", sep = ""),  width=heatmap_width, height=heatmap_height)
                      
                            #### enzyme coexpression map 
                            
                            v.gns_plus_names <- paste(gns.i, " (" ,v.gn.names[gns.i], ")", sep = "")
                            
                            m.coexpression.i <- matrix(0, nrow = length(v.gns_plus_names), ncol = length(v.gns_plus_names), 
                                                       dimnames = list(v.gns_plus_names, v.gns_plus_names))
                            
                            
                            for(k in 1:nrow(df.rho)){
                              k1 <- paste(df.rho$gn.1[k], " (" ,v.gn.names[df.rho$gn.1[k]], ")", sep = "")
                              k2 <- paste(df.rho$gn.2[k], " (" ,v.gn.names[df.rho$gn.2[k]], ")", sep = "")
                              m.coexpression.i[k1,k2] <- df.rho$v.pcc[k]
                            }
                            
                            m.coexpression.i <- m.coexpression.i + t(m.coexpression.i)
                            
                            p <- pheatmap(m.coexpression.i, cluster_rows = FALSE,cluster_cols = FALSE, legend = TRUE, 
                                          main = paste("Coexpression map of gene cluster ", v.gcs[i], "\n (colors indicate Pearson's correlation values)", sep = ""))
                            save_pheatmap_pdf(p, paste(foldername.results, "/coexpression_heatmaps/", v.gcs[i],".pdf", sep = ""),  width=heatmap_width, height=heatmap_height)
                                              

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
            }
          }
      }
    }
  }
  
  close(pb)

  if(!b.prepare){
    df.cluster_annotations <- df.cluster_annotations[order(-df.cluster_annotations$rank),]
    return(list(df.cluster_annotations=df.cluster_annotations, m.functionality = m.functionality))
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
#' @param n.cpus number of cores used
#' @param b.load_codifferentialAnalysis_monteCarloSimulation load codifferential expression data ("yes", "no")
#' @param pvalue_DifferentialExpression pvalue treshold for differential expession (default = 0.05)
#' @param probability_codifferentialExpression_MonteCarloSimulation probability threshold codifferential expression (default = 0.05)
#' @param pvalue_coexpression_distribution pvalue treshold context specific coexpression (default = 0.05)
#' @param pvalue_geneClusterPrediction pvalue gene cluster inference enzyme presence (default = 0.05)
#' @param pvalue_geneClusterConsistency pvalue gene cluster enzyme condition consistency (default = 0.05)
#' @param pvalue_treatment_per_condition pvalue gene pair condition annotation (default = 0.05)
#' @param pvalue_tissue_per_condition pvalue gene pair tissue annotation (default = 0.05)
#' @param th.consistent_condition_presence_percentage percentage of enyzmes with similar annotation per gene cluster (default = 0.8)
#' @param min_number_of_genes min number of enzymes per gene cluster (default = 3)
#' @param number_codifferentialExpression_MonteCarloSimulations number of codiffernetial expression background monte carlo simulations (default = 3)
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
#' l.data = load_datasets(input_format = "PCF2017",
#'                        filename.genes = "data/genes.txt",
#'                        filename.experiment_series_ids = "data/experiment_series_ids.txt",
#'                        filename.geneCluster = "data/ath_geneInCluster_3_aracyc.txt-labeled_NoHypoGenes.txt",
#'                        filename.foldChange_differentialExpression = "data/m.foldChange_differentialExpression.txt",
#'                        filename.pvalue_differentialExpression =	"data/m.pvalue_differentialExpression.txt",
#'                        filename.experiment_condition_tissue_annotation ="data/df.experiment_condition_annotation.txt")
#' 
#' message("run METACLUSTER")
#' l.results = run_METACLUSTER(m.foldChange_differentialExpression = l.data$m.foldChange_differentialExpression,
#'                             m.pvalue_differentialExpression = l.data$m.pvalue_differentialExpression,
#'                             df.experiment_condition_annotation = l.data$df.experiment_condition_annotation,
#'                             df.geneCluster = l.data$df.geneCluster,
#'                             tb.condition_treatments = l.data$tb.condition_treatments,
#'                             tb.condition_tissues = l.data$tb.condition_tissues,
#'                             n.cpus = 3,
#'                             b.load_codifferentialAnalysis_monteCarloSimulation = "yes",
#'                             pvalue_DifferentialExpression = 0.05,
#'                             probability_codifferentialExpression_MonteCarloSimulation = 0.95,
#'                             pvalue_coexpression_distribution = 0.05,
#'                             pvalue_geneClusterPrediction = 0.05,
#'                             pvalue_geneClusterConsistency = 0.05,
#'                             pvalue_treatment_per_condition = 0.05,
#'                             pvalue_tissue_per_condition = 0.05,
#'                             th.consistent_condition_presence_percentage = 0.95,
#'                             min_number_of_genes = 3,
#'                             number_codifferentialExpression_MonteCarloSimulations = 3,
#'                             number_conditionSpecificCoexpressionBackgroundGenePairs = 100,
#'                             min_number_condition_samples = 1,
#'                             seed = 1234,
#'                             heatmap_width = 10, 
#'                             heatmap_height = 5,
#'                             foldername.tmp = "tmp/", 
#'                             foldername.results = "results/")
#'                             
#' print(head(l.results$df.cluster_annotations))
#' 
#' evaluate_and_store_results(df.cluster_annotations=l.results$df.cluster_annotations,
#'                            df.experiment_condition_annotation = l.data$df.experiment_condition_annotation,
#'                            m.functionality=l.results$m.functionality, 
#'                            heatmap_width = 10, heatmap_height = 5,
#'                            foldername.results = "results/")
run_METACLUSTER = function(m.foldChange_differentialExpression,
                          m.pvalue_differentialExpression,
                          df.experiment_condition_annotation,
                          df.geneCluster,
                          tb.condition_treatments,
                          tb.condition_tissues,
                          n.cpus = 3,
                          b.load_codifferentialAnalysis_monteCarloSimulation = "yes",
                          pvalue_DifferentialExpression = 0.05,
                          probability_codifferentialExpression_MonteCarloSimulation = 0.95,
                          pvalue_coexpression_distribution = 0.05,
                          pvalue_geneClusterPrediction = 0.05,
                          pvalue_geneClusterConsistency = 0.05,
                          pvalue_treatment_per_condition = 0.05,
                          pvalue_tissue_per_condition = 0.05,
                          th.consistent_condition_presence_percentage = 0.95,
                          min_number_of_genes = 3,
                          number_codifferentialExpression_MonteCarloSimulations = 3,
                          number_conditionSpecificCoexpressionBackgroundGenePairs = 100,
                          min_number_condition_samples = 1,
                          seed = 1234,
                          heatmap_width = 10, heatmap_height = 6,
                          foldername.tmp = "tmp/",
                          foldername.results = "results/"){

  
  

  if(!file.exists(foldername.tmp)){
    dir.create(foldername.tmp)
  }
  
  if(!file.exists(foldername.results)){
    dir.create(foldername.results)
  }
  
  if(!file.exists(paste(foldername.results, "/condition_activity_per_cluster_heatmaps/", sep = ""))){
    dir.create(paste(foldername.results, "/condition_activity_per_cluster_heatmaps/", sep = ""))
  }
  if(!file.exists(paste(foldername.results, "/coexpression_heatmaps/", sep = ""))){
    dir.create(paste(foldername.results, "/coexpression_heatmaps/", sep = ""))
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
  th.gns.min = min_number_of_genes
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
    
    message("compute codifferential expression matrix")
    
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
    
    
    message(paste("Compute significant minimum level of overlapping conditions based on", n.sim, "runs of Monte Carlo simulations for weighted gene coordination"))
    
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
      rm(m.co_diff)
    }
    stopCluster(cl)
    print(Sys.time()-strt)
    
    v.th.likelihood <- numeric(n.sim)
    
    for(s in 1:n.sim){
      m.co_diff.shuffled <- readRDS(paste(foldername.tmp, "m.co_diff_shuffled.", s, ".", th.diffexp, ".rds", sep = ""))
      m.co_diff.shuffled <- m.co_diff.shuffled + t(m.co_diff.shuffled)
      diag(m.co_diff.shuffled) <- 0
      th.likelihood <- quantile(m.co_diff.shuffled, 0.95)
      print(th.likelihood)
      th.likelihood <- quantile(m.co_diff.shuffled, 0.99)
      print(th.likelihood)
      v.th.likelihood[s] <- quantile(m.co_diff.shuffled, th.prob.MC)
      rm(m.co_diff.shuffled)
    }
    
    print(Sys.time()-strt)

    saveRDS(v.th.likelihood, paste(foldername.tmp, "v.th.likelihood_", th.diffexp, "_" ,th.prob.MC, "_",n.sim, ".rds", sep = ""))
    
  }

  if(!file.exists(paste(foldername.tmp, "v.th.likelihood_", th.diffexp, "_" ,th.prob.MC,"_",n.sim, ".rds", sep = ""))){
    stop("Error: Monte Carlo simulation and / or co-differential expression matrix does not exist!")
  }
  
  m.co_diff <- readRDS(paste(foldername.tmp, "m.co_diff.", th.diffexp, ".rds", sep = ""))
  v.th.likelihood <- readRDS(paste(foldername.tmp, "v.th.likelihood_", th.diffexp, "_" ,th.prob.MC, "_",n.sim, ".rds", sep = ""))

  th.min_overlap <- round(mean(v.th.likelihood),0)

  #####
  
  message("compute sampled condition specific coexpression ... ")
  
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
                                                           th.consistent_condition_presence_percentage = th.consistent_condition_presence_percentage,
                                                           s.multipleHypthesisCorrection = "none",
                                                           seed=seed,
                                                           b.choseSimilarConditionInSimulation = TRUE,
                                                           b.signatureEnzymeAnalysis = FALSE,
                                                           b.prepare = TRUE, v.rho_random.gene_pairs = NULL,
                                                           foldername.results=foldername.results,
                                                           heatmap_width = heatmap_width, heatmap_height = heatmap_height)
  
  
  message("compute gene cluster condition specific coexpression ... ")
  
  
  # estimate genomic cluster correlation - integrate automatic adaption 
  l.results <- estimate_cluster_coexpression(m.fc = m.fc, m.de.ternary = m.de.ternary, m.co_diff = m.co_diff, 
                                         df.annotation = df.annotation,
                                         df.geneCluster = df.geneCluster,
                                         tb.condition_treatments = tb.condition_treatments,
                                         tb.condition_tissues = tb.condition_tissues,
                                         n.gene_pair_samples = n.gene_pair_samples,
                                         th.min_overlap = th.min_overlap, 
                                         th.min.samples = th.min.samples,
                                         th.gns.min = th.gns.min,
                                         th.p.cluster = th.p.cluster,
                                         th.rho_prob = th.rho_prob,
                                         th.pval.treatment=th.pval.treatment,
                                         th.pval.tissue=th.pval.tissue,
                                         th.consistent_condition_presence_percentage = th.consistent_condition_presence_percentage,
                                         s.multipleHypthesisCorrection = "none",
                                         seed=seed,
                                         b.choseSimilarConditionInSimulation = TRUE,
                                         b.signatureEnzymeAnalysis = FALSE,
                                         b.prepare = FALSE, v.rho_random.gene_pairs = v.rho_random.gene_pairs,
                                         foldername.results=foldername.results,
                                         heatmap_width = heatmap_width, heatmap_height = heatmap_height)

  
  return(l.results)
}
