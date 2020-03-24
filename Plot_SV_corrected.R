source('utils.R')


collect_result = NULL
for(tissue in sort(unique(samples$SMTSD))){
  tis =  gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
  

  df = tryCatch(read.table(paste0(outdir, 'Assoc_results_SVs/Association_test_',tis,'.txt'), 
                  sep='\t', header = T, stringsAsFactors = F),warning = function (w) {print(paste("No data available for tissue type", tis))}, error = function(f) {return("failed")}
 )
  if(inherits(df, "character")){
    print(paste(" ", "No SVA results for ", tis))
    next
  }

## remove tissues with less than 70 samples

  sample_in_the_tissue = samples %>% filter(SMTSD == tissue)
  if(nrow(sample_in_the_tissue) < 70){ next }
  collect_result = rbind(collect_result, df)
}

## remove tests with < 10 samples in either gender group
collect_result = collect_result[collect_result$p.value > -1, ]

## remove tissues with low gene expression (median TPM < 1)
collect_result$median_TPM = 10^(collect_result$median_TPM) - 1
collect_result = collect_result[collect_result$median_TPM > 1, ]

## compute FDR for each gene seperately
ace2 = collect_result[collect_result$Gene == 'ACE2', ]
TMPRSS2 = collect_result[collect_result$Gene == 'TMPRSS2', ]
ace2$FDR = p.adjust(ace2$p.value, method = 'BH')
TMPRSS2$FDR = p.adjust(TMPRSS2$p.value, method = 'BH')
collect_result = rbind(ace2, TMPRSS2)
collect_result = collect_result[,c("Tissue", "Gene", "Variable",
                                   "median_TPM", "coefficient", "p.value", "FDR")]
collect_result = collect_result[order(collect_result$p.value), ]
write.table(collect_result, paste0(outdir, 'gene_cov_correlations_SVA.csv'),sep = ',',row.names = F)



collect_result = collect_result[collect_result$FDR < 0.05, ]

# plot
for(i in seq(1, nrow(collect_result))){
  rowi = collect_result[i, ]
  tissue = rowi['Tissue']
  tis =  gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
  
  # read in 
  corrected_df = try(read.table(paste0(datadir, 'SVA_corrected/', tis, '_SV_removed.txt'),
                            sep = '\t', header = T, stringsAsFactors = F))
  if(inherits(corrected_df, "try-error")){
    next
  }
  
  y = corrected_df[, paste(rowi['Gene'], rowi['Variable'], sep = '_')]
  x = corrected_df[, as.character(rowi['Variable'])]
  if(as.character(rowi['Variable']) == 'AGE'){
    x = factor(x, levels = c("20-29", "30-39", "40-49", "50-59",
                             "60-69", "70-79"))
    color_p = 'Greens'
  }else{
    x = factor(x , levels = c(1,2), labels = c("Female", "Male"))
    color_p = 'Set1'
  }
  ggtitle_text = paste0(tissue, 
                   ':\n coef = ', round(rowi$coefficient,3),
                   ':\n median TPM = ', round(rowi$median_TPM, 3))
  df_for_plot = data.frame("Covariate"=x, "y"=y)
  gg = ggplot(aes(x = Covariate, y = y), data = df_for_plot) + 
    geom_boxplot(aes(fill = Covariate)) + 
    theme_bw() + 
    theme(axis.text.x = element_blank()) + 
    xlab("") + 
    ylab(paste0("Corrected expression of ", rowi['Gene'])) + 
    ggtitle(ggtitle_text) + 
    scale_fill_brewer(palette = color_p)
  
  tis =  gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', as.character(rowi['Tissue'])))))
  png(paste0(outdir, rowi['Gene'], '_',tis,'_',as.character(rowi['Variable']),'_SVA.png'), 
      res = 130, height = 500, width = 600)
  print(gg)
  dev.off()
  
}



