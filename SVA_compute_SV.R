library(sva)
source("utils.R")

dir.create(file.path(datadir, "SVs/"), showWarnings = FALSE)

estimate_SVs <- function(tissue){
  # sample covariates
  sample_in_the_tissue = samples %>% filter(SMTSD == tissue)
  sample_in_the_tissue = as.data.frame(sample_in_the_tissue)
  rownames(sample_in_the_tissue) = sample_in_the_tissue$SAMPID

  tis = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
  print(paste0('Tissue: ', tis))  
  # gene TPM
  gene_tpm_in_the_tissue = read.table(paste0(datadir, 'tissue_tpm/', tis, '_gene_TPM.txt'),
                                      sep = '\t', header = T, stringsAsFactors = F, row.names = 1)
  gene_tpm_in_the_tissue = gene_tpm_in_the_tissue[apply(gene_tpm_in_the_tissue, 1, function(x) sum(x>0.1)) > 0.2 * ncol(gene_tpm_in_the_tissue),]
  gene_tpm_in_the_tissue = log10(gene_tpm_in_the_tissue + 1)

  sample_in_the_tissue = sample_in_the_tissue[gsub("\\.", "-", colnames(gene_tpm_in_the_tissue)), ]

  ## 1). keep AGE_GROUP when estimating SVs
  mod = model.matrix(~AGE_GROUP,data=sample_in_the_tissue)
  mod0 = model.matrix(~1, data=sample_in_the_tissue)
  n.sv = num.sv(as.matrix(gene_tpm_in_the_tissue), mod, method = 'be')
  sva1 = sva(as.matrix(gene_tpm_in_the_tissue),mod,mod0,n.sv=n.sv)
  print(paste0('Keep AGE, there are ', n.sv, ' SVs removed'))
  
  ## 2). combine AGE_GROUP and SVs 
  cov = data.frame(sva1$sv)
  colnames(cov) = paste0('SV', seq(1, ncol(cov)))
  cov$AGE_GROUP = sample_in_the_tissue$AGE_GROUP
  write.table(cov, paste0(datadir, 'SVs/', tis, '_SVs_AGE.txt'), sep = '\t', row.names = F) 


  ### Estimate SVs preserving SEX
  # if only one sex, skip this step
  # if the tissue has less than 10 samples in either gender group, skip this step
  if(sum(table(sample_in_the_tissue$SEX) < 10) > 0){
	return
  }else if(length(table(sample_in_the_tissue$SEX)) == 1){
	return
  }else{
  	## keep SEX when estimating SVs
  	mod = model.matrix(~SEX,data=sample_in_the_tissue)
  	mod0 = model.matrix(~1, data=sample_in_the_tissue)
  	n.sv = num.sv(as.matrix(gene_tpm_in_the_tissue), mod, method = 'be')
  	sva1 = sva(as.matrix(gene_tpm_in_the_tissue),mod,mod0,n.sv=n.sv)
  	print(paste0('Keep SEX, there are ', n.sv, ' SVs removed'))

  	## combine SEX and SVs 
	cov = data.frame(sva1$sv)
  	colnames(cov) = paste0('SV', seq(1, ncol(cov)))
  	cov$SEX = sample_in_the_tissue$SEX
  	write.table(cov, paste0(datadir, 'SVs/', tis, '_SVs_SEX.txt'), sep = '\t', row.names = F) 
  }	
}



## Among these, samples were selected based on donor genotype
## availability and a threshold of at least 70 samples per tissue
for (tissue in sort(unique(samples$SMTSD))){
	estimate_SVs(tissue)
}
