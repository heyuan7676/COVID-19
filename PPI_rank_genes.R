source("utils.R")


### Didn't work out!

dir.create(file.path(datadir, "SVs/"), showWarnings = FALSE)

samples_used_in_GTEx <- function(){
  samples_used = c()
  for (tissue in sort(unique(samples$SMTSD))){
  	sample_in_the_tissue = samples %>% filter(SMTSD == tissue)
  	sample_in_the_tissue = as.data.frame(sample_in_the_tissue)
  	# restrict to the samples used in GTEx
  	tissue_name = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
   	genotype_PCs  = tryCatch(read.table(paste0(datadir, 'GTEx_Analysis_v8_eQTL_covariates/',tissue_name,'.v8.covariates.txt'),
                                 sep='\t', header = T, stringsAsFactors = F, row.names = 1), 
				 warning = function (w) {print(paste("No data available for tissue type", tissue_name))}, error = function(f) {return("failed")})
  	if(inherits(genotype_PCs, "character")){
    		print(paste(" ", "Skipping tissue", tissue_name))
    	next
  	}

  	sample_in_the_tissue = sample_in_the_tissue %>% filter(SUBJID %in% gsub("\\.", "-", colnames(genotype_PCs)))
	samples_used = c(samples_used, sample_in_the_tissue$SAMPID)
   }

   return(samples_used)
}


readin_tpm <- function(){
	tpm.whole = fread(paste0(datadir, 'gene_pc_lc_tpm.txt'))
	all.genes = as.character(tpm.whole %>% pull(Name))
	all.genes = sapply(all.genes, function(x) strsplit(x, '\\.')[[1]][1])

	## subset columns
	samples.used = samples_used_in_GTEx()
	samples = samples %>% filter(SAMPID %in% samples.used)
	tpm.exp = t(tpm.whole[, ..samples.used])

	## subset rows
	colnames(tpm.exp) = all.genes
	idxxx = apply(tpm.exp, 2, function(x) sum(x>0.1) > 0.2 * nrow(tpm.exp))
	print(paste0('Remove ', ncol(tpm.exp) - sum(idxxx), ' genes with low expression'))
	tpm.exp = tpm.exp[, idxxx]

	## log transformation
	tpm.exp = log10(tpm.exp + 1)
	return(tpm.exp)
}


#This outputs a table listing each "Tissue", "Gene", "Variable", "Median_TPM","coefficient", "pvalue", FDR
rank_TestGene <- function(Test_gene, Test_gene_name){ 
    Test_gene_tpm = data.frame("SAMPID" = rownames(tpm.exp.matrix), "geneExp" = tpm.exp.matrix[, Test_gene])
    Test_gene_tpm_cov = merge(Test_gene_tpm, samples, by = 'SAMPID')
    ### remove missing data
    cols = c("SAMPID", "AGE", "SEX", "DTHHRDY", "SMRIN", "SMTSISCH", "SMEXNCRT", "COHORT", "RACE", "SMTSD", "geneExp")
    Test_gene_tpm_cov = Test_gene_tpm_cov[, cols]
    Test_gene_tpm_cov.complete = Test_gene_tpm_cov[complete.cases(Test_gene_tpm_cov), ]
    #print(paste0("Removed ", nrow(Test_gene_tpm_cov)-nrow(Test_gene_tpm_cov.complete), " data points with missing data"))
    #print(paste0("Test with ", nrow(Test_gene_tpm_cov.complete), " data points"))

    ### fit on geneEXP
    #model   = lm(geneExp~AGE_GROUP+SEX+factor(DTHHRDY)+SMRIN+SMTSISCH+SMEXNCRT,
    #             data = Test_gene_tpm_cov.complete)
    model   = lm(geneExp~AGE+SEX+factor(DTHHRDY)+SMRIN+SMTSISCH+SMEXNCRT+factor(COHORT)+factor(RACE),
                 data = Test_gene_tpm_cov.complete)
    Test_gene_tpm_cov.complete$corrected_geneEXP = residuals(model)
    write.table(Test_gene_tpm_cov.complete, paste0(datadir, 'corrected_gene_exp/', Test_gene_name, '.csv'), sep=',', row.names = F, quote = F)


    ### compare the median of corrected expression
    corrected_exp = Test_gene_tpm_cov.complete %>% group_by(SMTSD) %>% summarize(Median = median(corrected_geneEXP)) 
    ranks = corrected_exp$SMTSD[rev(order(corrected_exp$Median))]

    return(ranks)
}


datadir_PPI = './PPI_data/'
dir.create(file.path(outdir, 'PPI_results'), showWarnings = FALSE)
dir.create(file.path(datadir, 'corrected_gene_exp'), showWarnings = FALSE)


## read in gene TPM
tpm.exp.matrix = readin_tpm()

PPI_gene_lists = read.table(paste0(datadir_PPI, 'PPI_interactions.csv'), sep = ',',
                            header = T, stringsAsFactors = F)

tissue_ranks = list()
for(tissue in sort(unique(samples$SMTSD))){
  tissue_ranks[[tissue]] = c()
}

columns_names = c()
for(i in seq(1, nrow(PPI_gene_lists))){
  rowi = PPI_gene_lists[i, ]
  print(c(as.character(rowi[1]), as.character(rowi[2])))
  rank_tis = rank_TestGene(as.character(rowi[2]), as.character(rowi[1]))
  if(is.null(rank_tis)){
    next
  }
  for(k in seq(1, length(rank_tis))){
    tis = rank_tis[k]
    tissue_ranks[[tis]] = c(tissue_ranks[[tis]], k)
  }
  columns_names = c(columns_names, as.character(rowi[1]))
}

tissue_ranks_df = t(as.data.frame(tissue_ranks))
colnames(tissue_ranks_df) = columns_names
write.table(tissue_ranks_df, paste0(file.path(outdir, 'PPI_results/'), 'Rank_median_in_tissues_corrected_for_COVs.csv'),
            sep=',', quote = F, row.names = T)








 
