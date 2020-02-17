
trim.leading <- function (x) sub("^\\s+", "", x)


trim.trailing <- function (x) sub("\\s+$", "", x)


trim <- function (x) gsub("^\\s+|\\s+$", "", x)



perform_significance_tests = function(df.cluster_annotations,
                                      df.geneCluster,
                                      min_number_of_genes = 2){
  
  df.cluster_annotations = subset(df.cluster_annotations, df.cluster_annotations$number_of_codiff_expressed_genes >= min_number_of_genes)
  v.enz.predicted = paste(as.character(df.cluster_annotations$condition_and_tissue_specific_genes), collapse = " / ")
  v.enz.predicted = unique(unlist(strsplit(v.enz.predicted, " / ")))
  
  message("signature enzyme enrichment - compare among enzymes")
  n.sigEnzymes <- length(which(df.geneCluster$Gene.Name != ""))
  n.enz <- nrow(df.geneCluster)
  
  df.geneCluster.predicted = subset(df.geneCluster, df.geneCluster$Gene.ID %in% v.enz.predicted)
  n.enz_predicted = length(v.enz.predicted)
  n.sig.enz.counter <- length(which(df.geneCluster.predicted$Gene.Name != ""))
  
  n.background <- n.enz
  hitInSample <- n.sig.enz.counter
  sampleSize <- n.enz_predicted
  hitInPop <-  n.sigEnzymes  
  failInPop <- n.background - hitInPop
  
  foldChange <- (hitInSample / sampleSize) / (hitInPop / n.background)
  p.val <- print(phyper(hitInSample, hitInPop, failInPop, sampleSize, lower.tail = FALSE))
  
  message(foldChange)
  message(p.val)
  
  message("Schlapfer et al. high confidence overlap")
  
  df.Schlapfer_hc = read.table("data/high_confidence_geneInCluster_3_aracyc.txt-labeled_NoHypoGenes.txt", header = T, sep ="\t", stringsAsFactors = F)
  n.sigEnzymes.schlapfer <- length(which(df.Schlapfer_hc$GeneName != ""))

  v.gcs_predicted = unique(df.cluster_annotations$cluster.ID)
  v.clusters.highConfidence = unique(df.Schlapfer_hc$ClusterID)
  
  n.background <- length(unique(df.geneCluster$Cluster.ID))
  hitInSample <- length(intersect(v.clusters.highConfidence, v.gcs_predicted))
  sampleSize <- length(v.gcs_predicted)
  hitInPop <-  length(v.clusters.highConfidence)  
  failInPop <- n.background - hitInPop
  
  foldChange <- (hitInSample / sampleSize) / (hitInPop / n.background)
  
  
  p.val <- print(phyper(hitInSample, hitInPop, failInPop, sampleSize, lower.tail= FALSE))
  
  message(hitInSample)
  message(p.val)
  message(foldChange)
  
  message("sig. enz")
  
  n.background <- nrow(df.Schlapfer_hc)
  hitInSample <- n.sig.enz.counter
  
  sampleSize <- n.enz_predicted
  hitInPop <-  n.sigEnzymes.schlapfer  
  
  
  failInPop <- n.background - hitInPop
  
  foldChange <- (hitInSample / sampleSize) / (hitInPop / n.background)
  # p.val <- print(phyper(hitInSample, hitInSample + hitInPop, failInPop + sampleSize - hitInSample, sampleSize, lower.tail = FALSE))
  
  p.val = fisher.test(matrix(c(hitInSample, hitInPop, sampleSize-hitInSample, failInPop), 2, 2), alternative='greater')$p.val;
  
  message(foldChange)
  message(p.val)
  
  

}


#' Load dependency function
#'
#' This function installs or loads dependencies needed
#' @param 
#' @keywords 
#' @export
#' @examples
#' install_and_load_libraries()
install_and_load_libraries <- function(){
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("multtest")
  
  list.of.packages <- c("metap", "reshape2","doParallel", "pheatmap", "foreach")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  # loading required libraries
  require(reshape2)
  require(foreach)
  require(doParallel)
  require(pheatmap)
  require(metap)
  
}



#' Load dataset function
#'
#' This function loads a datasets
#' @param input_format "custom" or "PCF2017" (default = "PCF2017")
#' @param filename.genes genes (rows of the expression datasets)
#' @param filename.sample_ids_differentialExpression experimental datasets (columns of the expression datasets)
#' @param filename.geneCluster filename gene clusters
#' @param filename.foldChange_differentialExpression differential expression data (fold changes)
#' @param filename.pvalue_differentialExpression  differential expression data (p-values)
#' @param filename.experiment_condition_tissue_annotation experiment to treatment and tissue annotation
#' @keywords 
#' @export
#' @examples
#' load_datasets()
load_datasets = function(input_format = "PCF2017",
                         filename.genes = "",
                         filename.sample_ids_differentialExpression  = "",
                         filename.geneCluster = "",
                         filename.foldChange_differentialExpression = "",
                         filename.pvalue_differentialExpression =	"",
                         filename.experiment_condition_tissue_annotation =	""){
  
  
  df.geneCluster = load_gene_cluster_data(filename.geneCluster=filename.geneCluster, input_format=input_format)
  
  genes = read.table(filename.genes, header = F, sep = "\t", stringsAsFactors = F)[,1]
  experiment_series_ids = read.table(filename.sample_ids_differentialExpression, header = F, sep = "\t", stringsAsFactors = F)[,1]
  experiment_series_ids = as.character(experiment_series_ids)
  
  df.annotation = read.table(filename.experiment_condition_tissue_annotation, sep = "\t", stringsAsFactors = F, header = T)
  df.annotation$unique_ID = as.character(df.annotation$unique_ID)
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
    df.annotation$number_series[i] = tb.experiment_series_ids[as.character(df.annotation$unique_ID[i])]
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
    tb.condition_treatments[i] = length(idx_1) + length(idx_2)# sum(df.annotation$number_series[idx_1]) + sum(df.annotation$number_series[idx_2])
  }
  
  
  return(list(m.foldChange_differentialExpression=m.foldChange_differentialExpression,
              m.pvalue_differentialExpression=m.pvalue_differentialExpression,
              df.experiment_condition_annotation=df.annotation,
              df.geneCluster=df.geneCluster,
              tb.condition_treatments=tb.condition_treatments,
              tb.condition_tissues=tb.condition_tissues))
  
}



#' Load gene cluster Function
#'
#' This function loads gene cluster information 
#' @param filename.geneCluster includes gene cluster filenames as well as the input  
#' @param input_format format of plant cluster finder (PCF2017_enzymes_only), PCF2017 or custom 
#' @keywords 
#' @export
#' @examples
#' load_gene_cluster_data()
load_gene_cluster_data = function(filename.geneCluster = "", 
                                  input_format = "PCF2017_enzymes_only"){
  
  format = input_format
  
  v.colnames_mandatory = c("Cluster.ID", "Gene.ID", "Gene.Name")
  if(format == "custom"){
    df.geneCluster <- read.table(filename.geneCluster, sep = "\t", header = T, fill = T, stringsAsFactors = FALSE, quote = "")
    if(nrow(df.geneCluster) == 0){
      stop("Error: no gene clusters found")
    }
    if(!all(v.colnames_mandatory %in% names(df.geneCluster))){
      stop(paste("error: could not find all mandatory columns in file:", paste(v.colnames_mandatory, collapse = ", ")))
    }
  }else if(format == "PCF2017" | format == "PCF2017_enzymes_only"){
    # FORMAT the gene cluster file
    df.geneCluster <- read.table(filename.geneCluster, sep = "\t", header = FALSE, skip = 1, fill = T, stringsAsFactors = FALSE, quote = "")
    if(nrow(df.geneCluster) == 0){
      stop("Error: no gene clusters found")
    }
    names(df.geneCluster) <- c("id", "Gene.Name", "Gene.ID", "Cluster.ID", "rxn.id", "ec", "pwy.id","pwy_of_gene", "Main_functional_domains", "SM_functional_domains")
    # df.geneCuster$Cluster.ID <- gsub(" ", "", df.geneCuster$Cluster.ID)
    for(j in 1:ncol(df.geneCluster)){
      df.geneCluster[,j]<- sapply(df.geneCluster[,j], trim)
    }
    
    if(format == "PCF2017_enzymes_only"){
      df.geneCluster = subset(df.geneCluster, !df.geneCluster$rxn.id %in% c("[]", "[]/[]", "[]/[]/[]", "[]/[]/[]/[]", "[]/[]/[]/[]/[]", "[]/[]/[]/[]/[]/[]"))
    }
    
  }else{
    stop(paste("error: gene cluster format needs to be \"PCF2017\" or \"custom\" or \"PCF2017_enzymes_only\" ", sep = ""))
  }
  
  return(df.geneCluster)
}


#' save heatmapt
#'
#' This function saves a heatmap plot
#' @param x matrix 
#' @keywords 
#' @examples
#' save_pheatmap_pdf()
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


#' A result statistics
#'
#' This function analysis the results
#' @param df.cluster_annotations inferred cluster condition annotations 
#' @param df.experiment_condition_annotation gene cluster acitivity heatmap
#' @param min_number_of_genes
#' @param v.gc_validated
#' @param heatmap_height
#' @param heatmap_height
#' @keywords 
#' @export
#' @examples
#' cat_function()
evaluate_and_store_results = function(df.cluster_annotations,
                                      df.experiment_condition_annotation,
                                      tb.condition_treatments = l.data$tb.condition_treatments,
                                      tb.condition_tissues = l.data$tb.condition_tissues,
                                      
                          
                                      min_number_of_genes = 2,
                                      v.gc_validated = c("C628_3", "C463_3", "C615_3", "C641_3"),
                                      heatmap_width = 10, heatmap_height = 6, fontsize = 7, fontsize_row = 9, fontsize_col = 9,
                                      foldername.results = "results/"){
  


  if(!file.exists(foldername.results)){
    dir.create(foldername.results)
  }
  
  if(!file.exists(paste(foldername.results, "/condition_activity_per_cluster_heatmaps/", sep = ""))){
    dir.create(paste(foldername.results, "/condition_activity_per_cluster_heatmaps/", sep = ""))
  }
  
  df.cluster_annotations = subset(df.cluster_annotations, df.cluster_annotations$number_of_codiff_expressed_genes >= min_number_of_genes)
  
  
  write.table(df.cluster_annotations, paste(foldername.results, "/df.cluster_annotations.txt", sep = ""), row.names = F, sep = "\t")
  
  v.enz.predicted = paste(as.character(df.cluster_annotations$condition_and_tissue_specific_genes), collapse = " / ")
  v.enz.predicted = unique(unlist(strsplit(v.enz.predicted, " / ")))
  
  message("Statistics")
  v.gcs_predicted <- as.character(unique(df.cluster_annotations$cluster.ID))
  # v.enz_predicted <- as.character()
  # sig.enz.counter <- 0
  # for(i in 1:length(v.gcs_predicted)){
  #   df.cluster_annotations.i <- subset(df.cluster_annotations, df.cluster_annotations$cluster.ID == as.character(v.gcs_predicted[i]))
  #   v.gn.names.i <- unlist(strsplit(as.character(df.cluster_annotations.i$number_of_codiff_expressed_genes)[1], " / "))
  #   sig.enz.counter <- sig.enz.counter + length(which(v.gn.names.i != ""))
  #   v.enz_predicted <- c(v.enz_predicted, unlist(strsplit(as.character(df.cluster_annotations.i$codiff_expressed_genes)[1], " / ")))
  #   # number of coexpressed signature enzymes
  # }
  message("# enzymes predicted: ", length(unique(v.enz.predicted)))
  message("# gene clusters with context predicted: ", length(v.gcs_predicted))
  #message("# signature enzymes predicted: ", sig.enz.counter)
  
  message("known clusters found: ", paste(v.gc_validated[which(v.gc_validated %in% v.gcs_predicted)],collapse = ", "))
  
  
  message("plot functionality maps")
  
  v.conditionGroups = names(tb.condition_treatments)
  v.tissueGroups = names(tb.condition_tissues)
  
  m.functionality = matrix(0, nrow = length(v.conditionGroups), ncol = length(v.tissueGroups) + 1, dimnames = list(v.conditionGroups, c(v.tissueGroups, "non specific")) )
  m.functionality <- m.functionality[sort(rownames(m.functionality)),]
  
  m.availability = matrix(-1, nrow = nrow(m.functionality), ncol = ncol(m.functionality), dimnames = list(rownames(m.functionality), colnames(m.functionality)))
  for(x in 1:nrow(df.experiment_condition_annotation)){
    if(df.experiment_condition_annotation$condition_treatment_1[x] != ""){
      m.availability[df.experiment_condition_annotation$condition_treatment_1[x], df.experiment_condition_annotation$condition_tissue[x]] = 1
    }
    if(df.experiment_condition_annotation$condition_treatment_2[x] != ""){
      m.availability[df.experiment_condition_annotation$condition_treatment_2[x], df.experiment_condition_annotation$condition_tissue[x]] = 1
    }
  }
  
  
  if(length(v.gcs_predicted) > 0){
    
    for(i in 1:length(v.gcs_predicted)){
      
      df.cluster_annotations.i = subset(df.cluster_annotations, df.cluster_annotations$cluster.ID == v.gcs_predicted[i])
 
      m.functionality.i = matrix(1, nrow = length(v.conditionGroups), ncol = length(v.tissueGroups) + 1, dimnames = list(v.conditionGroups, c(v.tissueGroups, "non specific")) )
      m.functionality.i <- m.functionality.i[sort(rownames(m.functionality.i)),]
      
      for(j in 1:nrow(df.cluster_annotations.i)){
        
        condition = as.character(df.cluster_annotations.i$condition[j])
        tissue = as.character(df.cluster_annotations.i$tissue[j])
        pval = as.numeric(df.cluster_annotations.i$p.cluster[j])
        
        m.functionality.i[condition,tissue] = pval
      }  
      
      
      m.tmp = m.functionality.i
      
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
      
      m.rank <- - log(m.heatmap)
      
      # breaks = c( 1, seq(0.05, min(m.heatmap[!is.na(m.heatmap)]),-0.0001))
      
 
      
      # breaks = c( 1, seq(0.05, min(m.heatmap[!is.na(m.heatmap)]),-0.0001))
      # 
      # break_map = unique(as.numeric(m.heatmap))
      # names(break_map) = unique(as.numeric(m.rank))
      # break_map = break_map[!is.na(break_map)]
      # 
      # break_map = break_map[order(-break_map)]
      # 
      # breaks = as.numeric(names(break_map))
      # legend_labels = as.character(break_map)
      # 
  
      legend_breaks = c(0, 2, -log(0.001), -log(1e-10))
      legend_labels = c(">0.05","<0.05", "<0.001", "<1e-10")
      
      
      inc = 1

      breaks = c(seq(0, floor(max(m.rank[!is.na(m.rank)]) + 1), inc))
      color <- c( "black",  "yellow" ,colorRampPalette(c( "yellow", "yellow", "blue",  "blue"))(length(breaks) - 1))
      
      # plot pdf 3.5 - 9
      p = pheatmap(m.rank, color = color, breaks = breaks, legend_breaks = legend_breaks, legend_labels = legend_labels, border_color = "black",  show_rownames = T, show_colnames = T , cluster_rows = F, cluster_cols = F, treeheight_row = 0, treeheight_col = 0,
                   main = "Colors indicate the significance of the cluster to be active per condition and tissue. \n (Gray tiles indicate condition-tissue combinations not present in the expression dataset)",
                   fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
                   silent = T)
      
      save_pheatmap_pdf(p, paste(foldername.results, "/condition_activity_per_cluster_heatmaps/", v.gcs_predicted[i],".pdf", sep = ""),  width=heatmap_height, height= heatmap_width )
      
    }
    
    
    # global
    for(j in 1:nrow(df.cluster_annotations)){
      
      condition = as.character(df.cluster_annotations$condition[j])
      tissue = as.character(df.cluster_annotations$tissue[j])
      
      m.functionality[condition,tissue] = m.functionality[condition,tissue] + 1
    }  
    
    
    for(i in 1:nrow(m.availability)){
      for(j in 1:ncol(m.availability)){
        if(colnames(m.functionality)[j] != "non specific"){
          if(m.availability[i,j] == -1){
            m.functionality[i,j] = -1
          }  
        }
      }
    }
    
    m.functionality[m.functionality == -1] = NA
    
    m.heatmap <- m.functionality # t(m.functionality) #  m.MR # * m.regulatorActivity.pvalue
    
    breaks = c( 0, 1, seq(2,max(m.heatmap[!is.na(m.heatmap)]),1))
    color <- c( "black",  "orange" ,colorRampPalette(c( "orange", "red"))(length(breaks) - 1))
    
    # plot pdf 3.5 - 9
    p = pheatmap(m.heatmap, color = color, breaks = breaks, border_color = "orange",  show_rownames = T, show_colnames = T , cluster_rows = F, cluster_cols = F, treeheight_row = 0, treeheight_col = 0,
                 main = "Gene cluster condition functionality map \n (colors indicate numbers of clusters active per condition)",
                 fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
                 silent = T)
    
    
    save_pheatmap_pdf(p, paste(foldername.results, "/condition_activity_number_of_clusters.pdf", sep = ""), width=heatmap_width, height=heatmap_height)
    
    
  }
  

                       
}



