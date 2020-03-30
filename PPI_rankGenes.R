setwd('Desktop/COVID_19/scripts/Github/')
library(readr)
library(dplyr)

source('utils.R')



capture_tissue_ranks <- function(Test_gene, Test_gene_name){
  geneExp_across_tissue = data.frame()
  for(tissue in sort(unique(samples$SMTSD))){
    tissue_name = gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', tissue))))   
    gene_tpm_in_the_tissue = tryCatch(read.table(paste0(datadir, 'tissue_tpm/', tissue_name, '_gene_TPM.txt'),
                                                 sep = '\t', header = T, stringsAsFactors = F, row.names = 1),
                                      warning = function (w) {print(paste("No data available for tissue type", tissue_name))}, 
                                      error = function(f) {return ("failed")}
    )
    if(inherits(gene_tpm_in_the_tissue, "character")){
      print(paste(" ", "Skipping tissue", tissue_name))
      next
    }
    
    gene_tpm_in_the_tissue = log10(gene_tpm_in_the_tissue + 1)
    samples_used = colnames(gene_tpm_in_the_tissue)
    
    if(Test_gene %in% rownames(gene_tpm_in_the_tissue)){
      geneEXP = as.numeric(gene_tpm_in_the_tissue[Test_gene, samples_used])
    }else if(Test_gene %in% sapply(rownames(gene_tpm_in_the_tissue), function(x) strsplit(x, '\\.')[[1]][1])){
      tpm_matrix_gene_names = sapply(rownames(gene_tpm_in_the_tissue), function(x) strsplit(x, '\\.')[[1]][1])
      geneEXP = as.numeric(gene_tpm_in_the_tissue[which(tpm_matrix_gene_names == Test_gene), samples_used])
    }else{
      geneEXP = rep(0, length(samples_used))
    }
    geneEXP_in_tissue = data.frame("geneExp" = geneEXP, "tissue" = rep(tissue_name, length(geneEXP)))
    geneExp_across_tissue = rbind(geneExp_across_tissue, geneEXP_in_tissue)
  }
  
  summary_across_tissues = geneExp_across_tissue %>% 
    group_by(tissue) %>%
    summarize(Median=median(geneExp), Mean = mean(geneExp), Std = sd(geneExp),
              Min = min(geneExp), Max = max(geneExp))
  
  if(sum(geneExp_across_tissue$geneEXP != 0) == 0){
    return ()
  }
  
  ### rank the tissues
  ranks = summary_across_tissues[rev(order(summary_across_tissues$Median)), ] %>%
    pull(tissue)
  ranked_tissues = as.character(ranks)
  
  return (ranked_tissues)
}


datadir_PPI = './PPI_data/'
dir.create(file.path(outdir, 'PPI_results'), showWarnings = FALSE)

# gene TPM
PPI_gene_lists = read.table(paste0(datadir_PPI, 'PPI_interactions.csv'), sep = ',', 
                            header = T, stringsAsFactors = F)

tissue_ranks = list()
for(tissue in sort(unique(samples$SMTSD))){
  tissue_ranks[[tissue]] = c()
}

columns_names = c()
for(i in seq(1, nrow(PPI_gene_lists))[1:2]){
  rowi = PPI_gene_lists[i, ]
  rank_tis = capture_tissue_ranks(as.character(rowi[2]), as.character(rowi[1]))
  if(is.null(rank_tis)){
    next
  }
  for(k in seq(1, length(rank_tis))){
    tis = rank_tis[k]
    tissue_ranks[[tis]] = c(tissue_ranks[[tis]], k)
    columns_names = c(columns_names, as.character(rowi[1]))
  }
}

tissue_ranks_df = t(as.data.frame(tissue_ranks))
colnames(tissue_ranks_df) = columns_names
write.table(tissue_ranks_df, paste0(file.path(outdir, 'PPI_results'), 'Rank_median_in_tissues.csv'), 
            sep=',', quote = F, row.names = T)




