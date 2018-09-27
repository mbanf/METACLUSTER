# METACLUSTER - an R package for context-specific functionality analysis of metabolic gene clusters

## About
METACLUSTER facilitates comprehensive condition and tissue-specific functionality analysis of metabolic gene clusters based on a probabilistic framework for characterizing metabolic gene clusters using context-specific gene expression information

[Contact the author](mailto:mbanf.research@gmail.com) in case you've found a bug. 

## Installation
The easiest way to install `METACLUSTER` is through `devtools` (see OS specific notes on installing devtools at the end)

```
library(devtools)
install_github("METACLUSTER","mbanf")
```

## Usage

To run the METACLUSTER with the Schlapfer et al. 2017 A.thaliana gene cluster predictions data you can download all neccessary datasets from onedrive: [datasets_athaliana_w_gencluster_schlapfer2017](https://1drv.ms/u/s!Avm82Xhe9EZj1hRSMsADlU4T6EB3)


```
library(METACLUSTER) # load package

setwd("/User/home/athaliana_schlapfer2017") # set working directory to the dataset files

# install_and_load_libraries()


```

Load datasets parameters:

* `input_format` "custom" or "PCF2017" (default = "PCF2017")
* `filename.genes` genes (rows of the expression datasets)
* `filename.experiment_series_ids` experimental datasets (columns of the expression datasets)
* `filename.condition_groups` treatment and tissue to condition maps
* `filename.geneCluster` filename gene clusters
* `filename.foldChange_differentialExpression` differential expression data (fold changes)
* `filename.pvalue_differentialExpression`  differential expression data (p-values)
* `filename.experiment_condition_tissue_annotation` experiment to treatment and tissue annotation

```
l.data = load_datasets(input_format = "PCF2017",
                       filename.genes = "/data/genes.txt",
                       filename.experiment_series_ids = "/data/experiment_series_ids.txt",
                       filename.condition_groups = "/data/conditionGroups.txt",
                       filename.geneCluster = "/data/ath_geneInCluster_3_aracyc.txt-labeled_NoHypoGenes.txt",
                       filename.foldChange_differentialExpression =      "/data/m.foldChange_differentialExpression.txt",
                       filename.pvalue_differentialExpression =	"/data/m.pvalue_differentialExpression.txt",
                       filename.experiment_condition_tissue_annotation =	"/data/experiment_series_annotation_He_et_al_2015.txt")
```


METACLUSTER Parameter sets:

* `m.foldChange_differentialExpression` differential expression foldchange matrix - rows are genes, cols are experiments
* `m.pvalue_differentialExpression` differential expression pvalue matrix - rows are genes, cols are experiments
* `df.experiment_condition_annotation` experiment condition annotation
* `df.geneCluster` gene cluster dataset
* `tb.condition_treatments` table of conditions
* `tb.condition_tissues` table of tissues
* `v.conditionGroups` treatment and condition map
* `v.tissueGroups` tissue maps 
* `n.cpus` number of cores used
* `b.load_codifferentialAnalysis_monteCarloSimulation` load codifferential expression data ("yes", "no")
* `pvalue_DifferentialExpression` pvalue treshold for differential expession (default = 0.05)
* `probability_codifferentialExpression_MonteCarloSimulation` probability threshold codifferential expression (default = 0.05)
* `pvalue_coexpression_distribution` pvalue treshold context specific coexpression (default = 0.05)
* `pvalue_geneClusterPrediction` pvalue gene cluster inference enzyme presence (default = 0.05)
* `pvalue_geneClusterConsistency` pvalue gene cluster enzyme condition consistency (default = 0.05)
* `pvalue_treatment_per_condition` pvalue gene pair condition annotation (default = 0.05)
* `pvalue_tissue_per_condition` pvalue gene pair tissue annotation (default = 0.05)
* `th.consistent_condition_presence_percentage` percentage of enyzmes with similar annotation per gene cluster (default = 0.95)
* `min_number_of_genes` min number of enzymes per gene cluster (default = 3)
* `number_codifferentialExpression_MonteCarloSimulations` number of codiffernetial expression background monte carlo simulations (default = 3)
* `number_conditionSpecificCoexpressionBackgroundGenePairs` number of context specific coexpression simulation background gene pairs (default = 50)
* `min_number_condition_samples` minimum number of condition samples for significance test (default 1)
* `foldername.tmp` temp file folder name (default = /tmp)
* `foldername.results` results file folder name (default = /results)

```
l.results = run_METACLUSTER(m.foldChange_differentialExpression = l.data$m.foldChange_differentialExpression,
                            m.pvalue_differentialExpression = l.data$m.pvalue_differentialExpression,
                            df.experiment_condition_annotation = l.data$df.experiment_condition_annotation,
                            df.geneCluster = l.data$df.geneCluster,
                            tb.condition_treatments = l.data$tb.condition_treatments,
                            tb.condition_tissues = l.data$tb.condition_tissues,
                            v.conditionGroups = l.data$v.conditionGroups,
                            v.tissueGroups = l.data$v.tissueGroups,
                            n.cpus = 3,
                            b.load_codifferentialAnalysis_monteCarloSimulation = "no",
                            pvalue_DifferentialExpression = 0.05,
                            probability_codifferentialExpression_MonteCarloSimulation = 0.05,
                            pvalue_coexpression_distribution = 0.05,
                            pvalue_geneClusterPrediction = 0.05,
                            pvalue_geneClusterConsistency = 0.05,
                            pvalue_treatment_per_condition = 0.05,
                            pvalue_tissue_per_condition = 0.05,
                            th.consistent_condition_presence_percentage = 0.95,
                            min_number_of_genes = 3,
                            number_codifferentialExpression_MonteCarloSimulations = 3,
                            number_conditionSpecificCoexpressionBackgroundGenePairs = 50,
                            min_number_condition_samples = 1,
                            foldername.results = "results/",
                            foldername.tmp = "tmp/")
```

Next evaluate and store the results
```
print(head(l.results$df.cluster_annotations))
evaluate_and_store_results(df.cluster_annotations=l.results$df.cluster_annotations,
                           m.functionality=l.results$m.functionality, 
                           foldername.results = "results/")
```


Gene cluster functionality annotation inferred by METACLUSTER for the Schlapfer et al. 2017 A.thaliana gene cluster predictions data

![Alt text](/C666_3.jpg?raw=true "coexpression map")

Gene cluster context specific co-expression heatmap inferred by METACLUSTER of the C666 from Schlapfer et al. 2017
![Alt text](/functionality.jpg?raw=true "functionality map")

## Notes

We are currently compiling a set of helper functions to assist with compiling differential expression matrices (fold change and p-value)


Installation of devtools dependencies under Ubuntu (prior to installing devtools):
sudo apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

