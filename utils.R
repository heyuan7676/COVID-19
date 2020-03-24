library(dplyr)
library(data.table)
library(ggplot2)

readin_samples <- function(){
  samples = fread(paste0(datadir, 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'))
  # keep the covaraites for RNA-seq
  samples = samples %>% filter(SMAFRZE == 'RNASEQ')
  
  samples = samples[, c("SAMPID", "SMTSD", "SMRIN", "SMTSISCH", "SMEXNCRT", "SMUNMPRT")]
  samples$SUBJID = sapply(samples$SAMPID, function(x) paste(strsplit(x, '-')[[1]][1], strsplit(x, '-')[[1]][2], sep='-'))
  samples = merge(samples, donors, by = 'SUBJID')
  
  ## make age an ordinal variable
  samples$AGE_GROUP = 0
  samples[samples$AGE == '20-29', "AGE_GROUP"] = 1
  samples[samples$AGE == '30-39', "AGE_GROUP"] = 2
  samples[samples$AGE == '40-49', "AGE_GROUP"] = 3
  samples[samples$AGE == '50-59', "AGE_GROUP"] = 4
  samples[samples$AGE == '60-69', "AGE_GROUP"] = 5
  samples[samples$AGE == '70-79', "AGE_GROUP"] = 6
  
  samples$Gender = '0'
  samples[samples$SEX == 1, "Gender"] = 'Female'
  samples[samples$SEX == 2, "Gender"] = 'Male'
  samples$Gender = factor(samples$Gender, levels = c("Female", "Male"))
  

  return(samples)
}



datadir = './GTEx_data/'
outdir = './results/'
donors = fread(paste0(datadir, 'GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt'))
samples = readin_samples()

