pre_process_dataset = function(geo_dataset){
  geo_norm <- normaliseIllumina(geo_dataset)
  ids <- as.character(featureNames(geo_norm))
  qual <- unlist(mget(ids, illuminaHumanv3PROBEQUALITY, ifnotfound=NA))
  rem <- qual == "No match" | qual == "Bad" | is.na(qual)
  geo_filt_1 <- geo_norm[!rem,]
  geo_filt <- varFilter(geo_filt_1)
  return(geo_filt)
}

get_dataframes = function(f_gse, source_name){
  expression_df = exprs(f_gse)
  pheno_df = pData(f_gse)
  feature_df = fData(f_gse)
  indices_source = which(pheno_df$source_name_ch1 == source_name)
  pheno_df = pheno_df[indices_source,]
  expression_df = data.frame(exprs(gse_filt) [,indices_source])
  return(list(expression_df, pheno_df, feature_df))
}

run_limma = function(l_expression_df, l_pheno_df){
  design_df = data.frame(row.names = colnames(l_expression_df))
  design_df$sample = l_pheno_df[,"characteristics_ch1.4"]
  levels(design_df$sample) = c('AD','AsymAD','control')
  design = model.matrix(~ 0 + design_df$sample)
  colnames(design) <- c("AD","AsymAD","control")
  rownames(design) = rownames(l_pheno_df)
  
  indices_AD_vs_AsymAD = which(pheno_df$characteristics_ch1.4 != "disease state: control")
  indices_AD_vs_control = which(pheno_df$characteristics_ch1.4 != "disease state: AsymAD")
  indices_AsymAD_vs_control = which(pheno_df$characteristics_ch1.4 != "disease state: AD")
  
  sample_AD_vs_AsymAD = lmFit(l_expression_df[,indices_AD_vs_AsymAD], design[(which(design[,3] == 0)),-3])
  sample_AD_vs_control = lmFit(l_expression_df[,indices_AD_vs_control], design[(which(design[,2] == 0)),-2])
  sample_AsymAD_vs_control = lmFit(l_expression_df[,indices_AsymAD_vs_control], design[(which(design[,1] == 0)),-1])
  
  contrasts_AD_vs_AsymAD <- makeContrasts( AD - AsymAD , levels=design[,-3])
  contrasts_AD_vs_control <- makeContrasts( AD - control , levels=design[,-2])
  contrasts_AsymAD_vs_control <- makeContrasts( AsymAD - control , levels=design[,-1])
  
  sample_AD_vs_AsymAD_fit <- contrasts.fit(sample_AD_vs_AsymAD, contrasts_AD_vs_AsymAD )
  sample_AD_vs_control_fit <- contrasts.fit(sample_AD_vs_control, contrasts_AD_vs_control )
  sample_AsymAD_vs_control_fit <- contrasts.fit(sample_AsymAD_vs_control, contrasts_AsymAD_vs_control )
  
  sample_AD_vs_AsymAD_final <- eBayes(sample_AD_vs_AsymAD_fit)
  sample_AD_vs_control_final <- eBayes(sample_AD_vs_control_fit)
  sample_AsymAD_vs_control_final <- eBayes(sample_AsymAD_vs_control_fit)
  
  return(list(sample_AD_vs_AsymAD_final, sample_AD_vs_control_final, sample_AsymAD_vs_control_final))
}

annotate_results = function(limma_result, feature_df){
  anno = feature_df
  anno <- anno[,c("Symbol","Entrez_Gene_ID","Chromosome")]
  limma_result$genes = anno
  return(limma_result)
}

analyze_limma_results = function(a_limma_results){
  total_genes = length(subset(topTable(a_limma_results, number = Inf), P.Value<0.001)$Symbol)
  
  #volcano_names = ifelse(abs(a_limma_results$coefficients) >= 1, a_limma_results$genes$Symbol, NA)
  
  volcanoplot(a_limma_results, names = volcano_names, xlab = "Log2 fold change", ylab = NULL, pch = 16, cex = 0.35, main = colnames(a_limma_results$contrasts)[1])
  
  return(total_genes)
}

get_all_genes = function(limma_results,l_expression_df, l_feature_df){
  top_probes = topTable(limma_results, number = Inf)
  top_probes = top_probes[match(unique(top_probes$Symbol),top_probes$Symbol),]
  probe_list = rownames(top_probes)
  l_expression_df = l_expression_df[probe_list,]
  l_feature_df = l_feature_df[probe_list,]
  
  ids = l_feature_df$ID
  l_expression_df = l_expression_df[ids,]
  rownames(l_expression_df) = l_feature_df$Symbol
  return(l_expression_df)
}

get_DE_genes = function(limma_results,l_expression_df, l_feature_df){
  top_probes = topTable(limma_results, number = Inf)
  top_probes = top_probes[top_probes$adj.P.Val < 0.05,]
  probe_list = rownames(top_probes)
  l_expression_df = l_expression_df[probe_list,]
  l_feature_df = l_feature_df[probe_list,]
  
  l_feature_df = l_feature_df[!duplicated(l_feature_df$Symbol),]
  ids = l_feature_df$ID
  l_expression_df = l_expression_df[ids,]
  rownames(l_expression_df) = l_feature_df$Symbol
  return(l_expression_df)
}



run_CEMitool = function(c_exprs_dfs,c_pheno_df, beta, set_beta_){

  exprs_AD_vs_AsymAD = c_exprs_dfs[[1]]
  exprs_AD_vs_control = c_exprs_dfs[[2]]
  exprs_AsymAD_vs_control = c_exprs_dfs[[3]]
  
  indices_AD_vs_AsymAD = which(c_pheno_df$characteristics_ch1.4 != "disease state: control")
  indices_AD_vs_control = which(c_pheno_df$characteristics_ch1.4 != "disease state: AsymAD")
  indices_AsymAD_vs_control = which(c_pheno_df$characteristics_ch1.4 != "disease state: AD")
  
  exprs_AD_vs_AsymAD = exprs_AD_vs_AsymAD[,indices_AD_vs_AsymAD]
  exprs_AD_vs_control = exprs_AD_vs_control[,indices_AD_vs_control]
  exprs_AsymAD_vs_control = exprs_AsymAD_vs_control[,indices_AsymAD_vs_control]
  
  #AD vs AsymAD
  
  samples_anno_1 = data.frame(row.names = c(1:ncol(exprs_AD_vs_AsymAD)))
  samples_anno_1$SampleName = colnames(exprs_AD_vs_AsymAD)
  samples_anno_1$Class = c_pheno_df[colnames(exprs_AD_vs_AsymAD),"characteristics_ch1.4"]
  if(set_beta_[1] == TRUE){
    cem_AD_vs_AsymAD = cemitool(exprs_AD_vs_AsymAD, samples_anno_1, set_beta = beta[1])
  }
  else{
    cem_AD_vs_AsymAD = cemitool(exprs_AD_vs_AsymAD, samples_anno_1)
  }
  
  #AD vs control
  
  samples_anno_2 = data.frame(row.names = c(1:ncol(exprs_AD_vs_control)))
  samples_anno_2$SampleName = colnames(exprs_AD_vs_control)
  samples_anno_2$Class = c_pheno_df[colnames(exprs_AD_vs_control),"characteristics_ch1.4"]
  if(set_beta_[2] == TRUE){
    cem_AD_vs_control = cemitool(exprs_AD_vs_control, samples_anno_2, set_beta = beta[2])
  }
  else{
    cem_AD_vs_control = cemitool(exprs_AD_vs_control, samples_anno_2)
  }
  
  #AsymAD vs control
  
  samples_anno_3 = data.frame(row.names = c(1:ncol(exprs_AsymAD_vs_control )))
  samples_anno_3$SampleName = colnames(exprs_AsymAD_vs_control )
  samples_anno_3$Class = c_pheno_df[colnames(exprs_AsymAD_vs_control ),"characteristics_ch1.4"]
  if(set_beta_[3] == TRUE){
    cem_AsymAD_vs_control = cemitool(exprs_AsymAD_vs_control, samples_anno_3, set_beta = beta[3])
  }
  else{
    cem_AsymAD_vs_control = cemitool(exprs_AsymAD_vs_control, samples_anno_3)
  }
  return(list(cem_AD_vs_AsymAD, cem_AD_vs_control, cem_AsymAD_vs_control))
}



run_enrichment_and_ora = function(cemi_results, c_gmt){
  
  if (nmodules(cemi_results) == 0){
    return_statement = "No modules found"
  }
  else{
  #Module enrichment 
  cemi_results <- mod_gsea(cemi_results)
  cemi_results <- plot_gsea(cemi_results)
  show_plot(cemi_results, "gsea")
  
  #Over representation analysis 
  cemi_results = mod_ora(cemi_results, c_gmt)
  cemi_results <- plot_ora(cemi_results)
  plots = show_plot(cemi_results, "ora")
  return_statement = cemi_results
  }
  return(return_statement)
}

add_interactions = function(cemi_results, int_df){
  interactions_data(cemi_results) <- int_df 
  cemi_results <- plot_interactions(cemi_results) 
  show_plot(cemi_results, "interaction") 
}

get_pathways = function(cemi_results){
  modules_indices = which(gsea_data(cemi_results)$padj[2] < 0.05 & gsea_data(cemi_results)$padj[3] < 0.05)
  if(length(modules_indices)==0){
    final_return = "No significant modules found"
  }
  else{
    modules_indices = which(gsea_data(cemi_results)$padj[2] < 0.05 & gsea_data(cemi_results)$padj[3] < 0.05)
    modules = gsea_data(cemi_results)$padj[modules_indices, "pathway"]
    modules = modules[modules != "Not.Correlated"]
    if(length(modules) == 0){
      final_return = "Only non correlated modules significant"
    }
    else{
    p_vals = which(ora_data(cemi_results)$p.adjust < 0.05)
    if(length(p_vals) == 0){
      final_return = "No significant genes found"
    }
    else{
    new = ora_data(cemi_results)[p_vals,]
    new = new[new$Module != "Not.Correlated",]
    indices = vector('list', length = length(unique(new$Module)))
    j = 1
    for (i in 1:length(modules)){
      find_match = match(modules[i],new$Module)
      if (is.na(find_match)== FALSE){
        indices[[j]] = which(new$Module == modules[i])
        j = j+1
      }
    }
  pathways = new[unlist(indices),]
  pathways = pathways[order(pathways$p.adjust),]
  final_return = pathways[,c("Module","ID","geneID","p.adjust")]
  }
  }
  }
  return(final_return)
}

run_predictions = function(gene_exprs_df, pheno_df, genes, not_comp_cond){

  gene_exprs_df_t = as.data.frame(t(as.matrix(gene_exprs_df)))
  gene_exprs_df_t = gene_exprs_df_t[,genes]
  gene_exprs_df_t$target = pheno_df$characteristics_ch1.4
  gene_exprs_df_t = gene_exprs_df_t[gene_exprs_df_t$target != not_comp_cond,]
  gene_exprs_df_t$target = factor(gene_exprs_df_t$target)
  
  if(not_comp_cond == "disease state: control"){
    levels(gene_exprs_df_t$target) = c("AD", "AsymAD")
  }
  
  if(not_comp_cond == "disease state: AsymAD"){
    levels(gene_exprs_df_t$target) = c("AD", "control")
  }
  
  if(not_comp_cond == "disease state: AD"){
    levels(gene_exprs_df_t$target) = c("AsymAD", "control")
  }
  
  split = sample.split(gene_exprs_df_t$target, SplitRatio = 0.60)
  training_set = subset(gene_exprs_df_t, split == TRUE)
  test_set = subset(gene_exprs_df_t, split == FALSE )
  training_set = training_set %>% mutate_if(is.numeric, scale)
  test_set = test_set %>% mutate_if(is.numeric, scale)

  #SVM
  svm_obj = svm(formula = target ~., data = training_set, type = 'C-classification', kernel = 'linear')
  svm_pred =  predict(svm_obj, newdata = test_set[,genes])
  svm_cm = table(test_set$target,svm_pred)
  svm_accuracy = (sum((svm_pred == test_set$target)== TRUE)*100)/nrow(test_set)
  
  #Random Forests
  rf <- randomForest(target ~ .,data=training_set)
  rf_pred = predict(rf, newdata = test_set[,genes])
  rf_cm = table(test_set$target, rf_pred)
  rf_accuracy = (sum((rf_pred == test_set$target)== TRUE)*100)/nrow(test_set)
  
  #Neural Network
  nn=neuralnet(target ~ .,data= training_set, hidden=4, act.fct = "logistic",linear.output = FALSE)
  nn_pred = predict(nn, newdata = test_set[,genes], all.units = FALSE)
  results = vector("list", length = nrow(nn_pred))
  results[which(nn_pred[,1] > nn_pred[,2])] = levels(test_set$target)[1]
  results[which(nn_pred[,1] < nn_pred[,2])] = levels(test_set$target)[2]
  nn_cm = table(test_set$target, unlist(results))
  nn_accuracy = (sum((unlist(results) == test_set$target)== TRUE)*100)/nrow(test_set)
  
  return(list(svm_accuracy, rf_accuracy, nn_accuracy))
}
