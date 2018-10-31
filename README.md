# MERIT

## About
MERIT

[Contact the author](mailto:mbanf.research@gmail.com) in case you've found a bug. 

## Installation
The easiest way to install `MERIT` is through `devtools` (see OS specific notes on installing devtools at the end). 

```
# install.packages("devtools")
library(devtools)
install_github("MERIT","mbanf")

```

## Usage

To run the MERIT with the A. thaliana data you can download all neccessary datasets from onedrive: [datasets_athaliana_w_PMN_2017]


```
library(MERIT) # load package

setwd("/User/home/athaliana_PMN_2017") # set working directory to the dataset files


```

Load datasets parameters:

* `filename.genes` genes (rows of the expression datasets)
* `filename.experiment_series_ids` experimental datasets (columns of the expression datasets)
* `filename.foldChange_differentialExpression` differential expression data (fold changes)
* `filename.pvalue_differentialExpression`  differential expression data (p-values)
* `filename.experiment_condition_tissue_annotation` experiment to treatment and tissue annotation
* `filename.transcriptionfactor_annotation` transcription factor with family annotation 
* `filename.geneGroups` gene group dataset (here metabolic domains)

```
l.data  =  load_datasets(filename.genes = "data/genes.txt",
                         filename.experiment_series_ids = "data/experiment_series_ids.txt",
                         filename.foldChange_differentialExpression = "data/m.foldChange_differentialExpression.txt",
                         filename.pvalue_differentialExpression =	"data/m.pvalue_differentialExpression.txt",
                         filename.experiment_condition_tissue_annotation =	"data/df.experiment_condition_annotation.txt",
                         filename.transcriptionfactor_annotation = "data/df.transcriptionFactorAnnotation.txt", 
                         filename.geneGroups = "data/df.enzymes_w_metabolic_domains.txt")
```

MERIT Parameter sets:

!We set b.load_grn_inference = "yes", b.load_TFBS_inference = "yes, and b.load_treatment_tissue_inference = "yes for the PMN 2017 A.thaliana gene cluster predictions data, as we have pre-computed and provided all co-differential expression datasets  - for other datasets, set to "no"! - modular structrucar


* `b.load_grn_inference` load precomputed grns ("yes", "no")
* `b.load_TFBS_inference` load precomputed TFBS with GRN ("yes","no")
* `b.load_treatment_tissue_inference` load precomputed annotation filgerung ("yes","no")

* `m.foldChange_differentialExpression` differential expression foldchange matrix - rows are genes, cols are experiments
* `m.pvalue_differentialExpression` differential expression pvalue matrix - rows are genes, cols are experiments
* `df.experiment_condition_annotation` experiment condition annotation
* `df.geneGroups` 
* `tb.geneGroups`
* `v.geneGroups` 
* `l.geneGroups`
* `tb.condition_treatments` table of conditions
* `tb.condition_tissues` table of tissues
* `n.cpus` number of cores used

* `th.lead_grn_method` (default = 0.95)
* `th.support_grn_methods` (default = 0.95)
* `n.grnSupport` (default = 1)

* `file.TF_to_Motif_IDs` = "data/TF_to_Motif_IDs.txt",
* `file.TFBS_motifs` = "data/Transcription_factor_weight_matrix_Arabidopsis_thaliana.txt")
* `file.promoterSeq` = "data/TAIR10_upstream_1000_20101104.txt")
* `file.geneSeq` = "data/TAIR10_seq_20110103_representative_gene_model_updated.txt")
* `th.pre_tss` (default = 1000)
* `th.post_tss` (default = 200)
* `genome_nucleotide_distribution` A,C,G,T nucleotide distritbution (default = c(0.3253439, 0.1746561, 0.1746561, 0.3253439))
* `th.pval.known_motifs` = 0.05)

* `th.diffexp` = 0.05)
* `th.pval.treatment` = 0.05) 
* `th.pval.tissue` = 0.05)
* `th.min.samples` = 1)
* `s.multipleTestCorrection` = "none")
* `seed` (default = 1234)

* `th.min_number_targets` (default = 2),
* `th.min_number_MR_targets` (default = 2),
* `th.pval_masterRegulator` (default = 0.05)

* `importance.measure` (default = "impurity")
* `n.trees` (default = 1000)
* `n.lead_method_expression_shuffling` (default  = 3)
* `nbootstrap` (default = 100)
* `nstepsLARS` (default = 5)

* `foldername.tmp` temp file folder name (default = tmp/)
* `foldername.results` results file folder name (default = results/)

```
l.results = run_MERIT(b.load_grn_inference = "yes",
                      b.load_TFBS_inference = "yes",
                      b.load_treatment_tissue_inference = "yes",
                      m.foldChange_differentialExpression=l.data$m.foldChange_differentialExpression,
                      m.pvalue_differentialExpression=l.data$m.pvalue_differentialExpression,
                      df.experiment_condition_annotation=l.data$df.experiment_condition_annotation,
                      tb.condition_treatments=l.data$tb.condition_treatments,
                      tb.condition_tissues=l.data$tb.condition_tissues,
                      df.transcriptionFactorAnnotation=l.data$df.transcriptionFactorAnnotation, 
                      df.geneGroups=l.data$df.geneGroups,
                      tb.geneGroups=l.data$tb.geneGroups,
                      v.geneGroups=l.data$v.geneGroups,
                      l.geneGroups=l.data$l.geneGroups, 
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
                      foldername.results = "results/")
```

Next evaluate and store the results
```

# Results
# Step 1 - Gene regulatory network inference using ensemble regression with Monte Carlo based threshold selection 
l.res.grn = l.results$l.res.grn

# Step 2 - Transcription factor direct target promoter binding based filtering of gene regulatory link predictions
l.res.grn_tfbs = l.results$l.res.grn_tfbs

# Step3 - Context specific annotation and filtering of gene regulatory link predictions
l.res.link_annotation = l.results$l.res.link_annotation

# Step 4 - Master regulator hierarchy inference
l.res.MR_hierarchy = l.results$l.res.MR_hierarchy


format_results(l.grn_subnetworks = l.res.link_annotation$l.grn_subnetworks,
               tb.condition_tissue_differentialExpression = l.res.link_annotation$tb.condition_tissue_differentialExpression,
               l.Hierarchy=l.res.MR_hierarchy$l.Hierarchy, 
               l.Hierarchy_tfs_per_tier=l.res.MR_hierarchy$l.Hierarchy_tfs_per_tier,
               l.Hierarchy_nb_tfs_per_tier=l.res.MR_hierarchy$l.Hierarchy_nb_tfs_per_tier,
               l.df.masterRegulatorHierarchy=l.res.MR_hierarchy$l.df.masterRegulatorHierarchy,
               v.number_tiers=l.res.MR_hierarchy$v.number_tiers,
               m.MR_vs_conditions = l.res.MR_hierarchy$m.MR_vs_conditions,  
               l.MR_vs_geneGroups_given_condition = l.res.MR_hierarchy$l.MR_vs_geneGroups_given_condition, 
               number_of_conditions_per_master_regulator=l.res.MR_hierarchy$number_of_conditions_per_master_regulator,
               tb.condition_treatments=l.data$tb.condition_treatments,
               tb.condition_tissues=l.data$tb.condition_tissues,
               df.transcriptionFactorAnnotation=l.data$df.transcriptionFactorAnnotation, 
               df.geneGroups=l.data$df.geneGroups,
               tb.geneGroups=l.data$tb.geneGroups,
               v.geneGroups=l.data$v.geneGroups,
               l.geneGroups=l.data$l.geneGroups,
               th.pval = 0.05,
               foldername.results = "results/")

```


A) Overview of the MERIT framework.
![Alt text](/figure1.jpg?raw=true "functionality map")



## Notes

Installation of devtools dependencies under Ubuntu (prior to installing devtools):
sudo apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

Subsequently, install devtools in R:
install.packages("devtools")

