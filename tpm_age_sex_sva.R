library(sva)
source("utils.R")

test_association <- function(tissue){
  # sample covariates
  sample_in_the_tissue = samples %>% filter(SMTSD == tissue)
  
  # gene TPM
  tis = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
  gene_tpm_in_the_tissue = read.table(paste0(datadir, 'tissue_tpm/', tis, '_gene_TPM.txt'),
                                      sep = '\t', header = T, stringsAsFactors = F, row.names = 1)
  gene_tpm_in_the_tissue = gene_tpm_in_the_tissue[apply(gene_tpm_in_the_tissue, 1, function(x) sum(x>0.1)) > 0.2 * ncol(gene_tpm_in_the_tissue),]
  gene_tpm_in_the_tissue = log10(gene_tpm_in_the_tissue + 1)
  
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
  	ace2_result = c(coef, pvalue)
  	# save
  	control_model = lm(t(ace2)~as.matrix(cov)[,2:ncol(cov)])
  	save_df = NULL
  	save_df = cbind(save_df, control_model$residuals)
  }else{
	ace2_result = c(0, -1)
  }
  
  # TMPRSS2
  if('ENSG00000184012.11' %in% rownames(gene_tpm_in_the_tissue)){
  	TMPRSS2 = gene_tpm_in_the_tissue['ENSG00000184012.11', ]
  	fitsv = lm(t(TMPRSS2)~as.matrix(cov))
  	coef = summary(fitsv)$coefficients[2, 3]
  	pvalue = summary(fitsv)$coefficients[2, 4]
  	TMPRSS2_result = c(coef, pvalue)
  	# save
  	control_model = lm(t(TMPRSS2)~as.matrix(cov)[,2:ncol(cov)])
  	save_df = cbind(save_df, control_model$residuals)
  }else{
	TMPRSS2_result = c(0, -1)
  }
  

  ### Test association with SEX
  
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
  	ace2_result = rbind(ace2_result, c(coef, pvalue))
  	# save
  	control_model = lm(t(ace2)~as.matrix(cov)[,2:ncol(cov)])
  	save_df = cbind(save_df, control_model$residuals)
  }else{
	ace2_result = rbind(ace2_result, c(0, -1))
  }

  
  # TMPRSS2
  if('ENSG00000184012.11' %in% rownames(gene_tpm_in_the_tissue)){
  	TMPRSS2 = gene_tpm_in_the_tissue['ENSG00000184012.11',]
  	fitsv = lm(t(TMPRSS2)~as.matrix(cov))
  	coef = summary(fitsv)$coefficients[2, 3]
  	pvalue = summary(fitsv)$coefficients[2, 4]
  	TMPRSS2_result = rbind(TMPRSS2_result, c(coef, pvalue))
  	# save
  	control_model = lm(t(TMPRSS2)~as.matrix(cov)[,2:ncol(cov)])
  	save_df = cbind(save_df, control_model$residuals)
  }else{
	TMPRSS2_result = rbind(TMPRSS2_result, c(0, -1))
  }
  
  save_df = as.data.frame(save_df)
  colnames(save_df) = c('ACE2_AGE', 'TMPRSS2_AGE', 'ACE2_SEX', 'TMPRSS2_SEX') 
  write.table(save_df, paste0(datadir, 'SVA_corrected/', tissue, '_SV_removed.txt'))

 
  result = as.data.frame(rbind(ace2_result, TMPRSS2_result))
  colnames(result) = c("coefficient", "p-value")
  result$Variable = c("AGE", "SEX", "AGE", "SEX")
  result$Gene = c("ACE2", "ACE2", "TMPRSS2", "TMPRSS2")
  result$Tissue = tissue

  return(result)
}


result_all = NULL
for (tissue in sort(unique(samples$SMTSD))){
  print(tissue)
  tis_df = test_association(tissue)
  result_all = rbind(result_all, tis_df)
}

write.table(result_all, paste0(outdir, 'Association_test.txt'), sep='\t', row.names = F)
