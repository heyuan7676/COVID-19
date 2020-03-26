source("utils.R")

###
##Input: gene x TPM matrix, Donor attributes, Sample attributes, covariates (including the genotype PC2)
##All are available from GTEx portal

##Output: "gene_cov_correlations.txt"
##Each row is the test for one gene in one tissue, for either AGE or SEX
###


#Returns table of TPM and covariates for patients with given tissue type
readin_data_in_tissue <- function(tissue, Test_gene){
  # sample covariates
  sample_in_the_tissue = samples %>% filter(SMTSD == tissue)
  
  # read in genotype PCs
  tissue_name = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))
  #print(paste0("Tissue: ", tissue, "; Tis: ",tis))
  genotype_PCs  = tryCatch(read.table(paste0(datadir, 'GTEx_Analysis_v8_eQTL_covariates/',tissue_name,'.v8.covariates.txt'), 
                                 sep='\t', header = T, stringsAsFactors = F, row.names = 1), warning = function (w) {print(paste("No data available for tissue type", tis))}, error = function(f) {return("failed")}
                                  )
  if(inherits(genotype_PCs, "character")){
    print(paste(" ", "Skipping tissue", tissue_name))
    return()
  }
  genotype_PCs = genotype_PCs[1:5, ]
  genotype_PCs = as.data.frame(t(genotype_PCs))
  genotype_PCs$SUBJID = rownames(genotype_PCs)
  samples_used = rownames(genotype_PCs)
  
  # gene TPM
  gene_tpm_in_the_tissue = read.table(paste0(datadir, 'tissue_tpm/', tissue_name, '_gene_TPM.txt'),
                                      sep = '\t', header = T, stringsAsFactors = F, row.names = 1)
  gene_tpm_in_the_tissue = log10(gene_tpm_in_the_tissue + 1)
  colnames(gene_tpm_in_the_tissue) = sapply(colnames(gene_tpm_in_the_tissue), 
                                            function(x) paste(strsplit(x, '\\.')[[1]][1][1], strsplit(x, '\\.')[[1]][2], sep = '.'))
  
  if(Test_gene %in% rownames(gene_tpm_in_the_tissue)){
   Test_gene_tpm = data.frame("SUBJID" = samples_used)
   Test_gene_tpm$geneEXP = as.numeric(gene_tpm_in_the_tissue[Test_gene, samples_used])
  }
  
  # merge
  sample_in_the_tissue$SUBJID = sapply(sample_in_the_tissue$SUBJID, function(x) gsub("-","\\.", x))
  df_test = merge(sample_in_the_tissue, Test_gene_tpm, by = 'SUBJID')
  df_test = merge(df_test, genotype_PCs, by = 'SUBJID')
  
  return(df_test)
}


#This outputs a table listing each "Tissue", "Gene", "Variable", "Median_TPM","coefficient", "pvalue", FDR
check_Test_gene_LR <- function(Test_gene){
  collect_result = NULL
  for(tissue in sort(unique(samples$SMTSD))){
    ### read in
    exp_for_tiss = readin_data_in_tissue(tissue, Test_gene)
    if(is.null(exp_for_tiss)){
      next
    }
    
    print(paste0('Perform LR using confounders for ', tissue))
    ### remove missing data
    exp_for_tiss.complete = exp_for_tiss[complete.cases(exp_for_tiss), ]
    #print(paste0("Removed ", dim(exp_for_tiss)[1]-dim(exp_for_tiss.complete)[1], " data points with missing data"))
    #print(paste0("Test with ", dim(exp_for_tiss)[1], " data points"))

    ### fit on geneEXP
    model   = lm(geneEXP~PC1+PC2+PC3+PC4+PC5+AGE_GROUP+SEX+factor(DTHHRDY)+SMRIN+SMTSISCH+SMEXNCRT, 
                 data = exp_for_tiss.complete)
    AGE_GROUP_test = c(tissue, Test_gene, "AGE", median(as.numeric(exp_for_tiss$geneEXP)), 
                       summary(model)$coefficients[,3]["AGE_GROUP"], 
                       summary(model)$coefficients[,4]["AGE_GROUP"])
    collect_result = rbind(collect_result, AGE_GROUP_test)
    
    
    # if only one sex, skip this step
    # if the tissue has less than 10 samples in either gender group, skip this step
    if(sum(table(exp_for_tiss$SEX) < 10) > 0){
      SEX_test = c(tissue, Test_gene,"SEX",median(as.numeric(exp_for_tiss$geneEXP)), 0, -1)
    }else if(length(unique(exp_for_tiss$SEX)) == 2){
      SEX_test = c(tissue, Test_gene, "SEX", median(as.numeric(exp_for_tiss$geneEXP)),
                   summary(model)$coefficients[,3]["SEX"], 
                   summary(model)$coefficients[,4]["SEX"])
    }else{
      SEX_test = c(tissue, Test_gene,"SEX",median(as.numeric(exp_for_tiss$geneEXP)), 0, -1)
    }
    collect_result = rbind(collect_result, SEX_test)
  }
  
  
  collect_result = as.data.frame(collect_result)
  colnames(collect_result) = c("Tissue", "Gene", "Variable", "Median_TPM","coefficient", "pvalue")
  collect_result$coefficient = as.numeric(as.character(collect_result$coefficient))
  collect_result$pvalue = as.numeric(as.character(collect_result$pvalue))
  collect_result = collect_result[collect_result$pvalue > -1, ]
  collect_result = collect_result[order(collect_result$pvalue), ]
  
  collect_result$Median_TPM = as.numeric(as.character(collect_result$Median_TPM))
  collect_result$Median_TPM = 10^(collect_result$Median_TPM) - 1
  collect_result = collect_result[collect_result$Median_TPM > 1, ]
  collect_result$FDR = p.adjust(collect_result$pvalue, method = 'BH')
  collect_result = collect_result[order(collect_result$pvalue), ]
  
  collect_result$Tissue = as.character(collect_result$Tissue)
  collect_result$Variable = as.character(collect_result$Variable)

  write.table(collect_result, paste0(outdir, 'Association_tests_',Test_gene_name, '_LR.csv'), sep=',', row.names = F)
  
  return(collect_result)
}



### Plot: Gene - SEX
plot_gene_sex <- function(Test_gene, df){
  Gene_SEX = df[df$Variable == 'SEX', ]
  if(dim(Gene_SEX)[1] == 0){
    return ()
  }
  df_for_plot = NULL
  
  for(i in seq(1, dim(Gene_SEX)[1])){
    rowi = Gene_SEX[i, ] 
    tissue = as.character(rowi['Tissue'])
    
    ### read in
    exp_for_tiss = readin_data_in_tissue(tissue,Test_gene)
    exp_for_tiss.complete = exp_for_tiss[complete.cases(exp_for_tiss), ]
    
    ### fit the model
    model   = lm(geneEXP~PC1+PC2+PC3+PC4+PC5+AGE_GROUP+factor(DTHHRDY)+SMRIN+SMTSISCH+SMEXNCRT, 
                 data = exp_for_tiss.complete)
    exp_for_tiss.complete$corrected_expression = model$residuals
    
    df_for_plot = exp_for_tiss.complete[,c("SMTSD", "corrected_expression", "Gender")]
    df_for_plot$coefficient = rowi$coefficient
    df_for_plot$Median_TPM = rowi$Median_TPM
    #df_for_plot = rbind(df_for_plot, df_for_plot_i)
    
    ggtitle_text = paste0(df_for_plot$SMTSD, 
                                      ":\n coef = ", round(df_for_plot$coefficient, 3),
                                      ":\n median TPM = ", round(df_for_plot$Median_TPM,3))
    xlabs = paste(levels(df_for_plot$Gender),"\n(N=",table(df_for_plot$Gender),")",sep="")
    g_sex = ggplot(aes(x = Gender, y = corrected_expression), data = df_for_plot) + 
      geom_boxplot(aes(fill = Gender)) + 
      ggtitle(ggtitle_text) + 
      theme_bw() + 
      scale_x_discrete(labels=xlabs) + 
      xlab("") + 
      theme(legend.position = 'none') + 
      ylab(paste0("Corrected expression of ", Test_gene_name)) + 
      scale_fill_brewer(palette = 'Set1')
    
    tis_name = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue)))) 
    png(paste0(outdir,'plots/', Test_gene_name, '_',tis_name,'_SEX_LR.png'), 
        res = 130, height = 500, width = 600)
    print(g_sex)
    dev.off()
  }
  
}

### Plot: Gene - AGE
plot_gene_age <- function(Test_gene, df){
  df = df[df$Gene == Test_gene, ]
  Gene_AGE = df[df$Variable == 'AGE', ]
  if(dim(Gene_AGE)[1] == 0){
    return ()
  }
  df_for_plot = NULL
  
  for(i in seq(1, dim(Gene_AGE)[1])){
    rowi = Gene_AGE[i, ]
    tissue = as.character(rowi['Tissue'])
    
    ### read in
    exp_for_tiss = readin_data_in_tissue(tissue, Test_gene)
    exp_for_tiss.complete = exp_for_tiss[complete.cases(exp_for_tiss), ]
    
    ### fit the model
    model   = lm(geneEXP~PC1+PC2+PC3+PC4+PC5+SEX+factor(DTHHRDY)+SMRIN+SMTSISCH+SMEXNCRT, 
                 data = exp_for_tiss.complete)
    exp_for_tiss.complete$corrected_expression = model$residuals
    
    df_for_plot = exp_for_tiss.complete[,c("SMTSD","corrected_expression", "AGE")]
    df_for_plot$coefficient = rowi$coefficient
    df_for_plot$Median_TPM = rowi$Median_TPM
    df_for_plot = rbind(df_for_plot, df_for_plot)
    
    ggtitle_text = paste0(df_for_plot$SMTSD, 
                                      ":\n coef = ", round(df_for_plot$coefficient, 3),
                                      ":\n median TPM = ", round(df_for_plot$Median_TPM,3))
    xlabs = paste(names(table(df_for_plot$AGE)),"yr\n(N=",table(df_for_plot$AGE),")",sep="")
    g_AGE = ggplot(aes(x = AGE, y = corrected_expression), data = df_for_plot) + 
      geom_boxplot(aes(fill = AGE)) + 
      ggtitle(ggtitle_text) + 
      theme_bw() + 
      theme(legend.position = 'none') + 
      scale_x_discrete(labels=xlabs) + 
      xlab("") + 
      ylab(paste0("Corrected expression of ", Test_gene_name)) + 
      scale_fill_brewer(palette = 'Greens')
    
    tis_name = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue)))) 
    png(paste0(outdir, 'plots/', Test_gene_name, '_',tis_name,'_AGE_LR.png'), 
        res = 130, height = 500, width = 600)
    print(g_AGE)
    dev.off()
  }
  

}



args <- commandArgs(TRUE)
Test_gene = args[1]
Test_gene_name = args[2]

reg_result = check_Test_gene_LR(Test_gene)

#### Plot
reg_result = read.table(paste0(outdir, 'Association_tests_',Test_gene_name, '_LR.csv'),
                        sep= ',', header = T, stringsAsFactors = F)
sig = reg_result[reg_result$FDR < 0.1, ]
plot_gene_sex(Test_gene, sig)
plot_gene_age(Test_gene, sig)

