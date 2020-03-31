source('utils.R')

gtex_col = read.table('../../data/gtex_colors.txt', sep='\t', header = T, stringsAsFactors = F)
gtex_col = gtex_col[order(gtex_col$tissue_id), ]
gtex_col$tissue_color_hex = paste0("#", gtex_col$tissue_color_hex)

get_tissue_above_lung_ace2 <- function(){
  ace2 = read.table(paste0(datadir, 'corrected_gene_exp/ACE2.csv'), sep=',', header = T, stringsAsFactors = F)
  TMPRSS2 = read.table(paste0(datadir, 'corrected_gene_exp/TMPRSS2.csv'), sep=',', header = T, stringsAsFactors = F)
  
  ace2_med = ace2 %>% group_by(SMTSD) %>% summarize(Median = 10^median(corrected_geneEXP) - 1)
  lung_ace2 = ace2_med %>% filter(SMTSD == 'Lung') %>% pull(Median)
  qualify_tissues = ace2_med[ace2_med$Median >= lung_ace2, ] %>% pull(SMTSD)
  qualify_tissues = sapply(qualify_tissues, 
                           function(x)gsub(" ", "_", gsub('\\)', '', gsub(' \\(', '_', gsub(' - ', '_', x)))))
  return(qualify_tissues) 
}

PPI_gene_ranks = read.table('results/PPI_results/Rank_median_in_tissues_corrected_for_COVs.csv', sep = ',',
                            header = T, stringsAsFactors = F)
PPI_gene_ranks = PPI_gene_ranks[order(rownames(PPI_gene_ranks)), ]
PPI_gene_ranks$tissue = gtex_col$tissue_id
  
tissues_pass_lung = get_tissue_above_lung_ace2()
df_for_plot = reshape2::melt(PPI_gene_ranks)
#df_for_plot = df_for_plot %>% filter(tissue %in% tissues_pass_lung)
medians = df_for_plot %>% group_by(tissue) %>% summarize(Median = median(value))
df_for_plot$tissue = factor(df_for_plot$tissue,
                            levels = medians$tissue[order(medians$Median)])

rownames(gtex_col) = gtex_col$tissue_id
colors = gtex_col[medians$tissue[order(medians$Median)], "tissue_color_hex"]


g = ggplot(aes(x = tissue, y = value), data = df_for_plot) + 
  geom_boxplot(aes(fill  = tissue)) + 
  scale_fill_manual(values = colors) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.2),
        legend.position =  'none') + 
  scale_y_reverse() + 
  ylab("Rank of expression across tissues") + 
  xlab("")

png('results/PPI_results/PPI_ranks_corrected_for_cov.png', res = 120, width = 850, height = 800)
print(g)
dev.off()





############### Heatmap

library(ComplexHeatmap)
library(circlize)
library(readr)

assign_colors_to_baits <- function(){
  library(colortools)
  unique_baits = c(paste0("nsp", seq(1,5)), "nsp5_C145A",paste0("nsp", seq(6,15)),
                   "E", "M", "N","Spike",
                   "orf3a", "orf3b", "orf6", "orf7a", "orf7b", "orf8",
                   "orf9b", "orf9c", "orf10")
  unique_baits = paste('SARS-CoV2', unique_baits)
  unique_baits = unique_baits[unique_baits%in%PPI_dat$Bait ]
  unique_baits_cols = c()
  i = 1
  for(b in unique_baits[1:15]){
    unique_baits_cols[b] = sequential("blue", percentage = 6, plot=F)[i+2]
    i = i + 1
  }
  i = 1
  for(b in unique_baits[16:19]){
    unique_baits_cols[b] = sequential("red", percentage = 10, plot=F)[5+i]
    i = i + 1
  }
  
  i = 1
  for(b in unique_baits[20:27]){
    unique_baits_cols[b] = sequential("green", percentage=10, plot=F)[2+i]
    i = i+1
  }
  return(unique_baits_cols)
}

## read in PPI interactions
datadir_PPI = './PPI_data/'
PPI_dat = read_csv(paste0(datadir_PPI, 'PPI_TableS2.csv'))
unique_baits_cols = assign_colors_to_baits()

## read in gene expression ranks
PPI_gene_ranks = read.table('results/PPI_results/Rank_median_in_tissues_corrected_for_COVs.csv', sep = ',',
                            header = T, stringsAsFactors = F)
tissues_pass_lung = get_tissue_above_lung_ace2()
PPI_gene_ranks = PPI_gene_ranks[order(rownames(PPI_gene_ranks)), ]
rownames(PPI_gene_ranks) = gtex_col$tissue_id
PPI_gene_ranks = PPI_gene_ranks[tissues_pass_lung,]


## merge the two dataframes
baits_withGene_ranks = intersect(as.character(PPI_dat %>% pull(PreyGene)), colnames(PPI_gene_ranks))

PPI_dat = PPI_dat %>% filter(PreyGene %in% baits_withGene_ranks)
PPI_gene_ranks = PPI_gene_ranks[,as.character(PPI_dat %>% pull(PreyGene))]
PPI_gene_ranks_baits = as.character(PPI_dat %>% pull(Bait))


ht_opt(
  legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
  legend_labels_gp = gpar(fontsize = 8), 
  heatmap_column_names_gp = gpar(fontsize = 8),
  heatmap_column_title_gp = gpar(fontsize = 10),
  heatmap_row_title_gp = gpar(fontsize = 8)
)


bait_annotation = HeatmapAnnotation(bait = PPI_gene_ranks_baits, 
                       col = list(bait = unique_baits_cols),
                       annotation_name_side = "left",
                       annotation_legend_param = list(
                         bait = list(
                           direction = 'horizontal',
                           title = 'Bait', 
                           at = names(unique_baits_cols),
                           labels = names(unique_baits_cols)
                         )
                       ))


col_fun = colorRamp2(c(1, 49), c('#53bdf8','gray10'))
ht_cluster_genes = Heatmap(as.matrix(PPI_gene_ranks), 
             show_row_names =T, 
             show_row_dend = T, 
             show_column_dend = F,
             show_column_names = F,
             cluster_rows = T,
             cluster_columns = T,
             col = col_fun, 
             row_dend_side = "left",
             row_km = 6,   ## split by kmeans cluster 
             row_km_repeats = 100, 
             cluster_row_slices = F,
             bottom_annotation  = bait_annotation, 
             column_title = "SARS-CoV-2 genome components", 
             column_title_side = "bottom",
             heatmap_legend_param = list(
               title = "Gene expression Rank",
               at = c(1, 49),
               labels = c("1", "49"),
               direction = "horizontal",
               legend_width = unit(10, "cm"))
)



setEPS()
postscript("results/PPI_results/PPI_ranks_corrected_for_cov_heatmap_clusterGenes.eps", height = 5, width = 16)
draw(ht_cluster_genes, heatmap_legend_side = 'bottom') 
dev.off()



### annotate drugable genes
druggable_genes = read.table('PPI_data/druggable_genes.txt', 
                             sep='\t', header = T,stringsAsFactors = F)
druggable_genes = druggable_genes$Gene_Name

### interesting gene
druggable_genes = c("AGPS", "AP3B1", "ERC1", "FBN1", 
                    "LOX", "MOV10", "PABPC4", "QSOX2", "SLC25A21")
show_idx = c()
for(g in druggable_genes){
  show_idx = c(show_idx, which(colnames(PPI_gene_ranks) == g))
}


anno = columnAnnotation(genes = anno_mark(at = show_idx, 
                 labels = druggable_genes, 
                 labels_gp = gpar(fontface = 'italic')))

ht_cluster_genes = Heatmap(as.matrix(PPI_gene_ranks), 
                           show_row_names =T, 
                           show_row_dend = T, 
                           show_column_dend = F,
                           show_column_names = F,
                           cluster_rows = T,
                           cluster_columns = T,
                           col = col_fun, 
                           row_dend_side = "left",
                           row_km = 6,   ## split by kmeans cluster 
                           row_km_repeats = 100, 
                           cluster_row_slices = F,
                           bottom_annotation = bait_annotation, 
                           top_annotation  = anno,
                           column_title = "SARS-CoV-2 genome components", 
                           column_title_side = "bottom",
                           heatmap_legend_param = list(
                             title = "Gene expression Rank",
                             at = c(1, 49),
                             labels = c("1", "49"),
                             direction = "horizontal",
                             legend_width = unit(10, "cm"))
)




setEPS()
postscript("results/PPI_results/PPI_ranks_corrected_for_cov_heatmap_clusterGenes_v3.eps", height = 5, width = 16)
draw(ht_cluster_genes, heatmap_legend_side = 'bottom') 
dev.off()




### only the 49 genes
baits_withGene_ranks = intersect(as.character(PPI_dat %>% pull(PreyGene)), druggable_genes)

PPI_dat = PPI_dat %>% filter(PreyGene %in% baits_withGene_ranks)
PPI_gene_ranks = PPI_gene_ranks[,as.character(PPI_dat %>% pull(PreyGene))]
PPI_gene_ranks_baits = as.character(PPI_dat %>% pull(Bait))


ht_opt(
  legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
  legend_labels_gp = gpar(fontsize = 8), 
  heatmap_column_names_gp = gpar(fontsize = 8),
  heatmap_column_title_gp = gpar(fontsize = 10),
  heatmap_row_title_gp = gpar(fontsize = 8)
)


bait_annotation = HeatmapAnnotation(bait = PPI_gene_ranks_baits, 
                                    col = list(bait = unique_baits_cols),
                                    annotation_name_side = "left",
                                    annotation_legend_param = list(
                                      bait = list(
                                        direction = 'horizontal',
                                        title = 'Bait', 
                                        at = names(unique_baits_cols),
                                        labels = names(unique_baits_cols)
                                      )
                                    ))


col_fun = colorRamp2(c(1, 49), c('#53bdf8','gray10'))
show_idx = c()
for(g in druggable_genes){
  show_idx = c(show_idx, which(colnames(PPI_gene_ranks) == g))
}


anno = columnAnnotation(genes = anno_mark(at = show_idx, 
                                          labels = druggable_genes, 
                                          labels_gp = gpar(fontface = 'italic')))

ht_cluster_genes = Heatmap(as.matrix(PPI_gene_ranks), 
                           show_row_names =T, 
                           show_row_dend = T, 
                           show_column_dend = F,
                           show_column_names = F,
                           cluster_rows = T,
                           cluster_columns = T,
                           col = col_fun, 
                           row_dend_side = "left",
                           row_km = 6,   ## split by kmeans cluster 
                           row_km_repeats = 100, 
                           cluster_row_slices = F,
                           bottom_annotation = bait_annotation, 
                           top_annotation  = anno,
                           column_title = "SARS-CoV-2 genome components", 
                           column_title_side = "bottom",
                           heatmap_legend_param = list(
                             title = "Gene expression Rank",
                             at = c(1, 49),
                             labels = c("1", "49"),
                             direction = "horizontal",
                             legend_width = unit(10, "cm"))
)



setEPS()
postscript("results/PPI_results/PPI_ranks_corrected_for_cov_heatmap_clusterGenes_v3.eps", height = 5, width = 16)
draw(ht_cluster_genes, heatmap_legend_side = 'bottom') 
dev.off()


