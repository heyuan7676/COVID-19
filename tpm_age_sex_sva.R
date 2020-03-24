library(sva)
source("utils.R")

test_association <- function(tissue){
  # sample covariates
  sample_in_the_tissue = samples %>% filter(SMTSD == tissue)
  sample_in_the_tissue = as.data.frame(sample_in_the_tissue)
  rownames(sample_in_the_tissue) = sample_in_the_tissue$SAMPID

  tis = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
  
  # gene TPM
  gene_tpm_in_the_tissue = read.table(paste0(datadir, 'tissue_tpm/', tis, '_gene_TPM.txt'),
                                      sep = '\t', header = T, stringsAsFactors = F, row.names = 1)
  gene_tpm_in_the_tissue = gene_tpm_in_the_tissue[apply(gene_tpm_in_the_tissue, 1, function(x) sum(x>0.1)) > 0.2 * ncol(gene_tpm_in_the_tissue),]
  gene_tpm_in_the_tissue = log10(gene_tpm_in_the_tissue + 1)

  sample_in_the_tissue = sample_in_the_tissue[gsub("\\.", "-", colnames(gene_tpm_in_the_tissue)), ]

  ## in generate_tissue_wise_TPM.R: make sure that sample orders in sample_in_the_tissue is the same as in gene_tpm_in_the_tissue 
  save_df = NULL
  save_df = cbind(save_df, sample_in_the_tissue$SAMPID)
  save_df = cbind(save_df, sample_in_the_tissue$AGE)
  save_df = cbind(save_df, sample_in_the_tissue$Gender)
  save_df_cols = c('SAMPID', 'AGE', 'SEX')
 
  ### 1. Test association with AGE_GROUP
  
  ## 1). keep AGE_GROUP when estimating SVs
  mod = model.matrix(~AGE_GROUP,data=sample_in_the_tissue)
  mod0 = model.matrix(~1, data=sample_in_the_tissue)
  n.sv = num.sv(as.matrix(gene_tpm_in_the_tissue), mod, method = 'be')
  sva1 = sva(as.matrix(gene_tpm_in_the_tissue),mod,mod0,n.sv=n.sv)
  print(paste0('Keep AGE, there are ', n.sv, ' SVs removed'))
  
  ## 2). combine AGE_GROUP and SVs 
  cov = data.frame("AGE_GROUP" = sample_in_the_tissue$AGE_GROUP)
  cov = cbind(cov, sva1$sv)
 
  ## 3). fit linear model

  # ACE2
  if('ENSG00000130234.10' %in% rownames(gene_tpm_in_the_tissue)){
    ace2 = gene_tpm_in_the_tissue['ENSG00000130234.10',]
    fitsv = lm(t(ace2)~as.matrix(cov))
    coef = summary(fitsv)$coefficients[2, 3]
    pvalue = summary(fitsv)$coefficients[2, 4]
    ace2_result = c(coef, pvalue, median(as.numeric(ace2)))
    # save
    control_model = lm(t(ace2)~as.matrix(cov)[,2:ncol(cov)])
    save_df = cbind(save_df, control_model$residuals)
    save_df_cols = c(save_df_cols, 'ACE2_AGE')
  }else{
    ace2_result = c(0, -1, 0)
  }
  
  # TMPRSS2
  if('ENSG00000184012.11' %in% rownames(gene_tpm_in_the_tissue)){
    TMPRSS2 = gene_tpm_in_the_tissue['ENSG00000184012.11', ]
    fitsv = lm(t(TMPRSS2)~as.matrix(cov))
    coef = summary(fitsv)$coefficients[2, 3]
    pvalue = summary(fitsv)$coefficients[2, 4]
    TMPRSS2_result = c(coef, pvalue, median(as.numeric(TMPRSS2)))
    # save
    control_model = lm(t(TMPRSS2)~as.matrix(cov)[,2:ncol(cov)])
    save_df = cbind(save_df, control_model$residuals)
        save_df_cols = c(save_df_cols, 'TMPRSS2_AGE')
  }else{
    TMPRSS2_result = c(0, -1, 0)
  }
  

  ### Test association with SEX
  # if only one sex, skip this step
  # if the tissue has less than 10 samples in either gender group, skip this step
  if(sum(table(exp_for_tiss$SEX) < 10) > 0){
    ace2_result = rbind(ace2_result, c(0, -1, 0))
    TMPRSS2_result = rbind(TMPRSS2_result, c(0, -1, 0))
  }else if(length(table(sample_in_the_tissue$SEX)) == 1){
    ace2_result = rbind(ace2_result, c(0, -1, 0))
    TMPRSS2_result = rbind(TMPRSS2_result, c(0, -1, 0))
  }else{     

  ## keep SEX when estimating SVs
  mod = model.matrix(~SEX,data=sample_in_the_tissue)
  mod0 = model.matrix(~1, data=sample_in_the_tissue)
  n.sv = num.sv(as.matrix(gene_tpm_in_the_tissue), mod, method = 'be')
  sva1 = sva(as.matrix(gene_tpm_in_the_tissue),mod,mod0,n.sv=n.sv)
  print(paste0('Keep SEX, there are ', n.sv, ' SVs removed'))
  
  ## combine SEX and SVs 
  cov = data.frame("SEX" = sample_in_the_tissue$SEX)
  cov = cbind(cov, sva1$sv)
  
  ## fit linear model
  # ACE2
  if('ENSG00000130234.10' %in% rownames(gene_tpm_in_the_tissue)){
    ace2 = gene_tpm_in_the_tissue['ENSG00000130234.10',]
    fitsv = lm(t(ace2)~as.matrix(cov))
    coef = summary(fitsv)$coefficients[2, 3]
    pvalue = summary(fitsv)$coefficients[2, 4]
    ace2_result = rbind(ace2_result, c(coef, pvalue, median(as.numeric(ace2))))
    # save
    control_model = lm(t(ace2)~as.matrix(cov)[,2:ncol(cov)])
    save_df = cbind(save_df, control_model$residuals)
        save_df_cols = c(save_df_cols, 'ACE2_SEX')
  }else{
    ace2_result = rbind(ace2_result, c(0, -1, 0))
  }

  
  # TMPRSS2
  if('ENSG00000184012.11' %in% rownames(gene_tpm_in_the_tissue)){
    TMPRSS2 = gene_tpm_in_the_tissue['ENSG00000184012.11',]
    fitsv = lm(t(TMPRSS2)~as.matrix(cov))
    coef = summary(fitsv)$coefficients[2, 3]
    pvalue = summary(fitsv)$coefficients[2, 4]
    TMPRSS2_result = rbind(TMPRSS2_result, c(coef, pvalue, median(as.numeric(TMPRSS2))))
    # save
    control_model = lm(t(TMPRSS2)~as.matrix(cov)[,2:ncol(cov)])
    save_df = cbind(save_df, control_model$residuals)
        save_df_cols = c(save_df_cols, 'TMPRSS2_SEX')
  }else{
    TMPRSS2_result = rbind(TMPRSS2_result, c(0, -1, 0))
  }
  }

  save_df = as.data.frame(save_df)
  colnames(save_df) = save_df_cols 
  write.table(save_df, paste0(datadir, 'SVA_corrected/', tis, '_SV_removed.txt'), sep = '\t', quote = FALSE)

 
  result = as.data.frame(rbind(ace2_result, TMPRSS2_result))
  colnames(result) = c("coefficient", "p-value", "median_TPM")
  result$Variable = c("AGE", "SEX", "AGE", "SEX")
  result$Gene = c("ACE2", "ACE2", "TMPRSS2", "TMPRSS2")
  result$Tissue = tissue

 dir.create(file.path(outdir, "Assoc_results_SVs"), showWarnings = FALSE)
 
    write.table(result, paste0(outdir, 'Assoc_results_SVs/Association_test_',tis,'.txt'), sep='\t', row.names=F, quote = FALSE)

  return(result)
}



args <- commandArgs(TRUE)


## Among these, samples were selected based on donor genotype
## availability and a threshold of at least 70 samples per tissue
tissue_idx = as.numeric(args[1])
tissue = sort(unique(samples$SMTSD))[tissue_idx]
print(tissue)
tis_df = test_association(tissue)
