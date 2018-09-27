
trim.leading <- function (x) sub("^\\s+", "", x)


trim.trailing <- function (x) sub("\\s+$", "", x)


trim <- function (x) gsub("^\\s+|\\s+$", "", x)


#' Load dependency function
#'
#' This function installs or loads dependencies needed
#' @param 
#' @keywords 
#' @export
#' @examples
#' install_and_load_libraries()
install_and_load_libraries <- function(){
  
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
#' @param filename.experiment_series_ids experimental datasets (columns of the expression datasets)
#' @param filename.condition_groups treatment and tissue to condition maps
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
                         filename.experiment_series_ids = "",
                         filename.condition_groups = "",
                         filename.geneCluster = "",
                         filename.foldChange_differentialExpression = "",
                         filename.pvalue_differentialExpression =	"",
                         filename.experiment_condition_tissue_annotation =	""){
  

  df.geneCluster = load_gene_cluster_data(filename.geneCluster=filename.geneCluster, input_format=input_format)

  genes = read.table(filename.genes, header = F, sep = "\t", stringsAsFactors = F)[,1]
  experiment_series_ids = read.table(filename.experiment_series_ids, header = F, sep = "\t", stringsAsFactors = F)[,1]
  df.annotation = read.table(filename.experiment_condition_tissue_annotation, sep = "\t", stringsAsFactors = F, header = T)
  df.foldChange_differentialExpression = read.table(filename.foldChange_differentialExpression, header = F, sep = "\t", stringsAsFactors = F)
  df.pvalue_differentialExpression = read.table(filename.pvalue_differentialExpression, header = F, sep = "\t", stringsAsFactors = F)
  d.conditionGroups = read.table(filename.condition_groups, header = T, sep = "\t", stringsAsFactors = F)

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
  if(nrow(d.conditionGroups) == 0){
    stop("Error: no condition group maps found")
  }
  
  
  v.conditionGroups = d.conditionGroups$condition
  names(v.conditionGroups) = d.conditionGroups$treatment
  
  
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
  
  
  tb.condition_treatments = numeric(length(unique(v.conditionGroups)))
  names(tb.condition_treatments) = unique(v.conditionGroups)
  for(i in 1:length(tb.condition_treatments)){
    idx.i = which(names(tb.condition_treatments)[i] == v.conditionGroups)
    v.treatments.i = names(v.conditionGroups)[idx.i]
    idx_1 = which(df.annotation$condition_treatment_1 %in% v.treatments.i)
    idx_2 = which(df.annotation$condition_treatment_2 %in% v.treatments.i)
    tb.condition_treatments[i] = sum(df.annotation$number_series[idx_1]) + sum(df.annotation$number_series[idx_2])
  }
  
  
  return(list(m.foldChange_differentialExpression=m.foldChange_differentialExpression,
              m.pvalue_differentialExpression=m.pvalue_differentialExpression,
              df.experiment_condition_annotation=df.annotation,
              df.geneCluster=df.geneCluster,
              tb.condition_treatments=tb.condition_treatments,
              tb.condition_tissues=tb.condition_tissues,
              v.conditionGroups=v.conditionGroups, 
              v.tissueGroups=v.tissues))
  
}
  


#' Load gene cluster Function
#'
#' This function loads gene cluster information 
#' @param filename.geneCluster includes gene cluster filenames as well as the input  
#' @param input_format format of plant cluster finder (PCF2017) or custom 
#' @keywords 
#' @export
#' @examples
#' load_gene_cluster_data()
load_gene_cluster_data = function(filename.geneCluster = "", 
                                  input_format = "PCF2017"){
  
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
  }
  if(format == "PCF2017"){
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
  }else{
    stop(paste("error: gene cluster format needs to be \"PCF2017\" or \"custom\" ", sep = ""))
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
#' @param m.functionality gene cluster acitivity heatmap
#' @param width.heatmap
#' @param height.heatmap
#' @keywords 
#' @export
#' @examples
#' cat_function()
evaluate_and_store_results = function(df.cluster_annotations=df.cluster_annotations,
                                      m.functionality=m.functionality,
                                      foldername.results = "results/",
                                      width.heatmap = 7,
                                      height.heatmap = 7){
              
  
  write.table(df.cluster_annotations, paste(foldername.results, "/df.cluster_annotations.txt", sep = ""), row.names = F, sep = "\t")
  
  
  message("Statistics")
  v.gcs_predicted <- as.character(unique(df.cluster_annotations$cluster.ID))
  v.enz_predicted <- as.character()
  sig.enz.counter <- 0
  for(i in 1:length(v.gcs_predicted)){
    df.cluster_annotations.i <- subset(df.cluster_annotations, df.cluster_annotations$cluster.ID == as.character(v.gcs_predicted[i]))
    v.gn.names.i <- unlist(strsplit(as.character(df.cluster_annotations.i$names_of_enzymes_in_condition)[1], " / "))
    sig.enz.counter <- sig.enz.counter + length(which(v.gn.names.i != ""))
    v.enz_predicted <- c(v.enz_predicted, unlist(strsplit(as.character(df.cluster_annotations.i$ids_of_coexpressed_enzymes_in_condition)[1], " / ")))
    # number of coexpressed signature enzymes
  }
  message("# enzymes predicted: ", length(unique(v.enz_predicted)))
  message("# gene clusters with context predicted: ", length(v.gcs_predicted))
  message("# signature enzymes predicted: ", sig.enz.counter)
  
  
  message("plot functionality map")
  
  m.functionality <- m.functionality[rowSums(m.functionality) > 0, colSums(m.functionality) > 0]
  m.heatmap <- t(m.functionality) #  m.MR # * m.regulatorActivity.pvalue
  
  breaks = c(0, 1, seq(2,max(m.heatmap),1))
  color <- c("black",  "orange" ,colorRampPalette(c( "orange", "red"))(length(breaks) - 1))
  
  # plot pdf 3.5 - 9
  p = pheatmap(m.heatmap, color = color, breaks = breaks, border_color = "orange",  show_rownames = T, show_colnames = T , cluster_rows = F, cluster_cols = F, treeheight_row = 0, treeheight_col = 0)#fontsize = 1)
  save_pheatmap_pdf(p, paste(foldername.results, "/condition_activity_number_of_clusters.pdf", sep = ""), width=width.heatmap, height=height.heatmap)

}