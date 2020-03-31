source('utils.R')

gtex_col = read.table('../../data/gtex_colors.txt', sep='\t', header = T, stringsAsFactors = F)
gtex_col = gtex_col[order(gtex_col$tissue_id), ]
gtex_col$tissue_color_hex = paste0("#", gtex_col$tissue_color_hex)


### read in ACE2 and TMPRSS2 gene
plot_corrected_TPM <- function(Test_gene){
  df = read.table(paste0(datadir, 'corrected_gene_exp/',Test_gene,'.csv'), sep=',', header = T, stringsAsFactors = F)
  medians = df %>% group_by(SMTSD) %>% summarize(Median = median(corrected_geneEXP))
  
  colors = gtex_col$tissue_color_hex[rev(order(medians$Median))]
  df$SMTSD = factor(df$SMTSD, levels = medians$SMTSD[rev(order(medians$Median))])
  
  g = ggplot(aes(x = SMTSD, y = corrected_geneEXP), data = df) + 
    geom_boxplot(aes(fill  = SMTSD)) + 
    scale_fill_manual(values = colors) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.2),
          legend.position =  'none') + 
    ylab("Corrected expression values") + 
    xlab("")
  
  png(paste0('results/plots_corrected_exp/',Test_gene,'_corrected_for_cov.png'), res = 120, width = 850, height = 800)
  print(g)
  dev.off()
}

for(gene in c("AGPS", "AP3B1", "ERC1", "FBN1","LOX",
              "MOV10", "PABPC4", "QSOX2", "SLC25A21")){
  plot_corrected_TPM(gene)
}
plot_corrected_TPM('ACE2')
plot_corrected_TPM("TMPRSS2")
