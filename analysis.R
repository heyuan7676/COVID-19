source("utils.R")

tissue = 'Lung'
tissue_name = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
genotype_PCs  = tryCatch(read.table(paste0(datadir, 'GTEx_Analysis_v8_eQTL_covariates/',tissue_name,'.v8.covariates.txt'), 
                                    sep='\t', header = T, stringsAsFactors = F, row.names = 1), 
                         warning = function (w) {print(paste("No data available for tissue type", tis))}, 
                         error = function(f) {return("failed")}
                         )
sample_in_the_tissue = samples %>% filter(SMTSD == tissue)
sample_in_the_tissue$SUBJID = sapply(sample_in_the_tissue$SUBJID, function(x) gsub("-","\\.", x))

sample_in_the_tissue = sample_in_the_tissue%>% filter(SUBJID %in% colnames(genotype_PCs))
