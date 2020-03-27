setwd('Desktop/COVID_19/scripts/Github/')
library(readr)

datadir = './PPI_data/'

PPI_dat = read_csv(paste0(datadir, 'PPI_TableS2.csv'))

## map gene names to ensemble IDs
library(biomaRt)
geneIDs <- function(geneName){
  ensembl=useMart("ensembl")
  ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
  geneNames=getBM(attributes=c("ensembl_gene_id", "external_gene_name", 
                               "gene_biotype", "entrezgene_id", "chromosome_name"),filters= 'external_gene_name', 
                  values=geneName, mart = ensembl, verbose=F)
  
  ## restrict to chr1-22, chrX, chrY
  geneNames = geneNames[!grepl("CHR", geneNames$chromosome_name), ]
  return(geneNames)
}


gene_Info = geneIDs(PPI_dat$PreyGene)

## attach virus protein information
dat_save = merge(gene_Info, PPI_dat[,c("Bait", "PreyGene")], 
                 by.x = 'external_gene_name', by.y = 'PreyGene')

write.table(dat_save, paste0(datadir, 'PPI_interactions.csv'), sep = ',', row.names = F)



