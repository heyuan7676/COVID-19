library(sva)
source("utils.R")

dir.create(file.path(datadir, "SVA_corrected/"), showWarnings = FALSE)

test_association <- function(tissue){
  # sample covariates
  sample_in_the_tissue = samples %>% filter(SMTSD == tissue)

  ## remove tissues with less than 70 samples
  if(nrow(sample_in_the_tissue) < 70){ return ()}

  sample_in_the_tissue = as.data.frame(sample_in_the_tissue)
  rownames(sample_in_the_tissue) = sample_in_the_tissue$SAMPID

  tissue_name = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
  
  # gene TPM
  gene_tpm_in_the_tissue = read.table(paste0(datadir, 'tissue_tpm/', tissue_name, '_gene_TPM.txt'),
                                      sep = '\t', header = T, stringsAsFactors = F, row.names = 1)
  gene_tpm_in_the_tissue = log10(gene_tpm_in_the_tissue + 1)

  sample_in_the_tissue = sample_in_the_tissue[gsub("\\.", "-", colnames(gene_tpm_in_the_tissue)), ]

  ## in generate_tissue_wise_TPM.R: make sure that sample orders in sample_in_the_tissue is the same as in gene_tpm_in_the_tissue 
  save_df = NULL
  save_df = cbind(save_df, sample_in_the_tissue$SAMPID)
  save_df = cbind(save_df, sample_in_the_tissue$AGE)
  save_df = cbind(save_df, sample_in_the_tissue$Gender)
  save_df_cols = c('SAMPID', 'AGE', 'SEX')
 
  ### 1. Test association with AGE_GROUP
  ## fit linear model
  if(Test_gene %in% rownames(gene_tpm_in_the_tissue)){
    print(paste0('Perform LR using SVs for ', tissue))
    gene_y = gene_tpm_in_the_tissue[Test_gene,]
    # AGE
    cov = read.table(paste0(datadir, 'SVs/', tissue_name, '_SVs_AGE.txt'), sep = '\t', header = T)
	  cov = as.matrix(cov)
    fitsv = lm(t(gene_y)~cov)
    coef = summary(fitsv)$coefficients['covAGE_GROUP', 3]
    pvalue = summary(fitsv)$coefficients['covAGE_GROUP', 4]
    result_tested_gene = c(coef, pvalue, median(as.numeric(gene_y)))
    # save
    sv_cols = paste0("SV", seq(1, ncol(cov)-1))
    control_model = lm(t(gene_y)~as.matrix(cov)[,sv_cols])
    save_df = cbind(save_df, control_model$residuals)
    save_df_cols = c(save_df_cols, 'SVA_AGE')
 
    # SEX
    cov = tryCatch(read.table(paste0(datadir, 'SVs/',tissue_name,'_SVs_SEX.txt'),
                                 sep='\t', header = T, stringsAsFactors = F), 
		   warning = function (w) {print(paste("No data available for tissue type", tis))}, 
		   error = function(f) {return("failed")})
    if(inherits(cov, "character")){
	print(paste0(" ", "Skipping ", tissue_name, " for sex"))
	result_tested_gene = rbind(result_tested_gene, c(0, -1, 0))
    }else{
		cov = as.matrix(cov)
    	fitsv = lm(t(gene_y)~cov)
    	coef = summary(fitsv)$coefficients['covSEX', 3]
    	pvalue = summary(fitsv)$coefficients['covSEX', 4]
    	result_tested_gene = rbind(result_tested_gene, c(coef, pvalue, median(as.numeric(gene_y))))
    	# save
    	control_model = lm(t(gene_y)~as.matrix(cov)[,2:ncol(cov)])
    	save_df = cbind(save_df, control_model$residuals)
    	save_df_cols = c(save_df_cols, 'SVA_SEX')
    }
  }else{
    print(paste(Test_gene, " has no expression in ", tissue_name))
    result_tested_gene = rbind(c(0,-1,0), c(0, -1, 0))
  }

  save_df = as.data.frame(save_df)
  colnames(save_df) = save_df_cols 
  write.table(save_df, paste0(datadir, 'SVA_corrected/', Test_gene_name, '/', tissue_name, '_SV_removed.txt'), 
              sep = '\t', quote = FALSE, row.names = F)

  result_tested_gene = as.data.frame(result_tested_gene)
  colnames(result_tested_gene) = c("coefficient", "p-value", "median_TPM")
  result_tested_gene$Variable  = c("AGE", "SEX")
  result_tested_gene$Gene      = Test_gene
  result_tested_gene$Tissue    = tissue

  return(result_tested_gene)
}


check_Test_gene_SVA <- function(Test_gene){
  dir.create(file.path(datadir, "SVA_corrected/", Test_gene_name, '/'), showWarnings = FALSE)
  result = data.frame()
  for(tissue in sort(unique(samples$SMTSD))){
    result_tested_gene = tryCatch(test_association(tissue), warning = function (w) {print(paste("No data available for tissue type", tis))},
                                  error = function(f) {return("failed")})
    if(!inherits(result_tested_gene, "character")){
      result = rbind(result, result_tested_gene)
    }else{
      print(paste0("Skipping ", tissue))
    }
  }
  
  ## remove null test
  result = result[result$"p-value" > 0, ]
  
  ## remove tissues with low gene expression (median TPM < 1)
  result$median_TPM = 10^(result$median_TPM) - 1
  result = result[result$median_TPM > 1, ]
  
  ## compute FDR for each gene seperately
  result$FDR = p.adjust(result$"p-value", method = 'BH')
  result = result[,c("Tissue", "Gene", "Variable", "median_TPM", "coefficient", "p-value", "FDR")]
  result = result[order(result$"p-value"), ]
  write.table(result, paste0(outdir, 'Association_tests_', Test_gene_name,'_SVA.csv'), sep=',', row.names=F, quote = FALSE)
  
  return(result)
}


## plot
plot_one_row <- function(rowi){
  tissue = rowi['Tissue']
  tissue_name = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
  
  # read in 
  corrected_df = try(read.table(paste0(datadir, 'SVA_corrected/', Test_gene_name, '/', tissue_name, '_SV_removed.txt'),
                                sep = '\t', header = T, stringsAsFactors = F))
  if(inherits(corrected_df, "try-error")){
    return()
  }
  y = corrected_df[, paste('SVA', rowi['Variable'], sep = '_')]
  x = corrected_df[, as.character(rowi['Variable'])]
  variable = as.character(rowi['Variable'])
  if(variable == 'AGE'){
    x = factor(x, levels = c("20-29", "30-39", "40-49", "50-59",
                             "60-69", "70-79"))
    color_p = 'Greens'
  }else{
    x = factor(x , levels = c(1,2), labels = c("Male", "Female"))
    color_p = 'Set1'
  }
  ggtitle_text = paste0(tissue,
                        ':\n coef = ', round(rowi$coefficient,3),
                        ':\n median TPM = ', round(rowi$median_TPM, 3))
  df_for_plot = data.frame("Covariate"=x, "y"=y)
  
  if(variable == 'AGE'){
    xlabs = paste(levels(df_for_plot$Covariate),"yr\n(N=",table(df_for_plot$Covariate),")",sep="")
  }else{
    xlabs = paste(levels(df_for_plot$Covariate),"\n(N=",table(df_for_plot$Covariate),")",sep="")
  }
  gg = ggplot(aes(x = Covariate, y = y), data = df_for_plot) +
    geom_boxplot(aes(fill = Covariate)) +
    theme_bw() +
    scale_x_discrete(labels=xlabs) + 
    xlab("") +
    ylab(paste0("Corrected expression of ", Test_gene_name)) +
    ggtitle(ggtitle_text) +
    scale_fill_brewer(palette = color_p) + 
    theme(legend.position = 'none')
  
  png(paste0(outdir, 'plots/', Test_gene_name, '_',tissue_name,'_',as.character(rowi['Variable']),'_SVA.png'),
      res = 130, height = 500, width = 600)
  print(gg)
  dev.off()
}




args <- commandArgs(TRUE)
Test_gene = args[1]
Test_gene_name = args[2]

reg_result = check_Test_gene_SVA(Test_gene)

# plot
#reg_result = read.table(paste0(outdir, 'Association_tests_', Test_gene_name,'_SVA.csv'), 
#                        sep=',', header = T, stringsAsFactors = F)
sig = reg_result[reg_result$FDR < 0.05, ]
for(i in seq(1, nrow(sig))){
	rowI = sig[i, ]
	plot_one_row(rowI)
}

