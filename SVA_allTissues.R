library(sva)
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


estimate_SVs_all_tissues <- function(){

  tpm.exp.tissues = data.frame("Tissue" = rep("0", nrow(tpm.exp)), "samples" = rownames(tpm.exp))
  tpm.exp.tissues$Tissue = as.character(tpm.exp.tissues$Tissue)
  rownames(tpm.exp.tissues) = rownames(tpm.exp)
  for(tissue in sort(unique(samples$SMTSD))){
	samples_tis = as.character(samples %>% filter(SMTSD == tissue) %>% pull(SAMPID))
	tissue_name = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
	tpm.exp.tissues[samples_tis, "Tissue"] = tissue_name
  }
	
  ## 1). keep TISSUE_GROUP when estimating SVs
  mod = model.matrix(~factor(Tissue),data=tpm.exp.tissues)
  mod0 = model.matrix(~1, data=tpm.exp.tissues)
  n.sv = num.sv(as.matrix(tpm.exp), mod)
  sva1 = sva(as.matrix(tpm.exp),mod,mod0,n.sv=n.sv)
  print(paste0('Keep AGE, there are ', n.sv, ' SVs removed'))
  
  ## 2). combine samples, tissues and SVs 
  cov = data.frame(sva1$sv)
  colnames(cov) = paste0('SV', seq(1, ncol(cov)))
  cov$Tissue = tpm.exp.tissues$Tissue
  cov$SAMPID = tpm.exp.tissues$samples
  write.table(cov, paste0(datadir, 'SVs/All_samples_control_for_Tissue.txt'), sep = '\t', row.names = F, quote = F) 

}



tpm.whole = fread(paste0(datadir, 'gene_pc_lc_tpm.txt'))
all.genes = as.character(tpm.whole %>% pull(Name))

samples.used = samples_used_in_GTEx()
samples = samples %>% filter(SAMPID %in% samples.used)
tpm.exp = t(tpm.whole[, ..samples.used])

idxxx = apply(tpm.exp, 2, function(x) sum(x>0.1) > 0.2 * nrow(tpm.exp))
print(paste0('Remove ', ncol(tpm.exp) - sum(idxxx), ' genes with low expression'))

tpm.exp = tpm.exp[, idxxx]
tpm.exp = log10(tpm.exp + 1)


for (tissue in sort(unique(samples$SMTSD))){
	estimate_SVs(tissue)
}
