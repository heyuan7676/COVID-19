source("utils.R")

tpm.whole = fread(paste0(datadir, 'gene_pc_lc_tpm.txt'))

for (tissue in sort(unique(samples$SMTSD))){
  print(tissue)
  # sample covariates
  sample_in_the_tissue = samples %>% filter(SMTSD == tissue)
  # restrict to the samples used in GTEx
  tis = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
   genotype_PCs  = tryCatch(read.table(paste0(datadir, 'GTEx_Analysis_v8_eQTL_covariates/',tis,'.v8.covariates.txt'),
                                 sep='\t', header = T, stringsAsFactors = F, row.names = 1), warning = function (w) {print(paste("No data available for tissue type", tis))}, error = function(f) {return("failed")})
  if(inherits(genotype_PCs, "character")){
    print(paste(" ", "Skipping tissue", tis))
    next
  }

  samples_used = gsub("\\.", "-", colnames(genotype_PCs))
  sample_in_the_tissue = sample_in_the_tissue %>% filter(SUBJID %in% samples_used) 

 
  # gene TPM
  gene_tpm_in_the_tissue = tpm.whole[,.SD, .SDcols=sample_in_the_tissue$SAMPID] 
  gene_tpm_in_the_tissue = as.data.frame(gene_tpm_in_the_tissue)
  rownames(gene_tpm_in_the_tissue) = tpm.whole$Name
  write.table(gene_tpm_in_the_tissue, 
              paste0(datadir, 'tissue_tpm/', tis, '_gene_TPM.txt'), 
              sep='\t')
  rm(gene_tpm_in_the_tissue)
}