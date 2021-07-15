library(GEOquery)
library(genefilter)
library(limma)
library(dplyr)
library(caTools)
library(beadarray)
library(illuminaHumanv3.db)
library(CEMiTool)
library(ggplot2)
library(randomForest)
library(e1071)
library(neuralnet)
source("C:/Users/rajra/OneDrive/Desktop/dge_full/alzheimers/functions.R")
gmt_in = read_gmt("C:/Users/rajra/OneDrive/Desktop/dge_full/alzheimers/gsea_reactome_pathway.xls")
int_df = read.table(file = "C:/Users/rajra/OneDrive/Desktop/dge_full/alzheimers/string_interactions.tsv", sep = '\t', header = TRUE)
int_df = int_df[,1:2]


#For temporal samples 
#gse = getGEO(filename = system.file("GSE118553_series_matrix.txt.gz", package = "GEOquery"))
gse <- getGEO(filename= "C:/Users/rajra/OneDrive/Desktop/dge_full/microarray_geo_dge_alzheimers_full/GSE118553_series_matrix.txt.gz")
gse_filt = pre_process_dataset(gse)
dataframes = get_dataframes(gse_filt, "Temporal_Cortex")
expression_df = dataframes[[1]]
pheno_df = dataframes[[2]]
feature_df = dataframes[[3]]

limma_results = run_limma(expression_df, pheno_df)
annotated_limma_results = lapply(limma_results, annotate_results, feature_df)
analyze_limma = lapply(annotated_limma_results, analyze_limma_results)

all_exprs_dfs = lapply(annotated_limma_results, get_all_genes, expression_df, feature_df)
new_exprs_dfs = lapply(annotated_limma_results, get_DE_genes, expression_df, feature_df)

cemi_results_all = run_CEMitool(all_exprs_dfs, pheno_df, c(16,16), c(TRUE,TRUE,FALSE))
cemi_results_dge = run_CEMitool(new_exprs_dfs, pheno_df, c(20,16, 8), c(TRUE,TRUE,TRUE))

cemi_results_all = lapply(cemi_results_all, run_enrichment_and_ora ,gmt_in)
cemi_results_dge = lapply(cemi_results_dge, run_enrichment_and_ora, gmt_in)

#Show protein-protein interactions
add_interactions(cemi_results_all[[1]],int_df)

#Get significant pathways 
pathways_all = lapply(cemi_results_all, get_pathways)
pathways_dge = lapply(cemi_results_dge[1:2], get_pathways)

#Get genes in the selected pathway
genes_list = pathways_dge[[1]][pathways_dge[[1]]$ID == "REACTOME_NEURONAL_SYSTEM",]$geneID
genes = strsplit(genes_list,"/")

#Run predictions on full dataset using selected genes
predictions = run_predictions(all_exprs_dfs[[1]],pheno_df, genes, "disease state: control")
