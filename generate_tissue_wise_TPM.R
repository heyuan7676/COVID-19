source("utils.R")

tpm.whole = fread(paste0(datadir, 'gene_pc_lc_tpm.txt'))

for (tissue in sort(unique(samples$SMTSD))[40:55]){
  print(tissue)
  # sample covariates
  sample_in_the_tissue = samples %>% filter(SMTSD == tissue)
  
  # gene TPM
  gene_tpm_in_the_tissue = tpm.whole[,.SD, .SDcols=sample_in_the_tissue$SAMPID] 
  
  # save
  tis = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
  rownames(gene_tpm_in_the_tissue) = tpm.whole$Name
  write.table(gene_tpm_in_the_tissue, 
              paste0(datadir, 'tissue_tpm/', tis, '_gene_TPM.txt'), 
              sep='\t')
  rm(gene_tpm_in_the_tissue)
}




