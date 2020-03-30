library(Seurat)
library(dplyr)

### process the barcodes - cell types 
cell_identity = read.table('sc_data/GSE92332_atlas_UMIcounts.txt', nrows = 1, stringsAsFactors = F)
cell_identity = as.character(cell_identity[1,])

batch_id = sapply(cell_identity, function(x) strsplit(x, '_')[[1]][1])
barcode_id = sapply(cell_identity, function(x) strsplit(x, '_')[[1]][2])
cell_type = sapply(cell_identity, function(x) strsplit(x, '_')[[1]][3])
cell_identity = data.frame("batch_id" = batch_id, "barcode_id" = barcode_id, "cell_type" = cell_type)
print(nrow(cell_identity))

cell_identity = cell_identity[!duplicated(cell_identity$barcode_id), ]
rownames(cell_identity) = cell_identity$barcode_id
print(nrow(cell_identity))


### signature genes
signature_genes = c()
celltypes = c("Enterocyte_Mature_Distal", "Enterocyte_Mature_Proximal",
              "Enteroendocrine", "goblet", "paneth", "tuft",
              "Enterocyte_Immature_Distal", "Enterocyte_Immature_Proximal",
              "stem", "TA.G2")
for(celltype in celltypes){
  sg = read.table(paste0('sc_data/', celltype, '_signature_genes.txt'), sep='\t', stringsAsFactors = F)
  signature_genes = c(signature_genes, sg$V1)
}


correct_count <- function(samples_name){
  sc.data = Read10X(data.dir = paste0('sc_data/GSE92332_RAW/', samples_name))
  sc = CreateSeuratObject(counts = sc.data, min.cells = 3, min.features = 800)
  sc = NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
  exp = sc[["RNA"]]@data
  
  cells_annotated = intersect(colnames(exp), barcode_id)
  exp = exp[,cells_annotated]
  colnames(exp) = as.character(cell_identity[cells_annotated, "cell_type"])
  exp = as.data.frame(exp)
  
  exp = exp[signature_genes, ]
  
  return(exp)
}


df_batchs = correct_count('Atlas1')
for(batch in c(2,5,6,9,10)){
  df_count = correct_count(paste0('Atlas', batch))
  df_batchs = cbind(df_batchs, df_count)
}
df_batchs[is.na(df_batchs)] = 0

df_batchs_save = as.data.frame(t(df_batchs))
df_batchs_save$celltype = colnames(df_batchs)
df_batchs_save_mean = aggregate(.~celltype, df_batchs_save, mean)

rownames(df_batchs_save_mean) = df_batchs_save_mean$celltype
df_batchs_save_mean = df_batchs_save_mean[,seq(2, ncol(df_batchs_save_mean))]
df_batchs_save_mean = as.data.frame(t(df_batchs_save_mean))
write.table(df_batchs_save_mean, 'sc_data/process_signature_genes.csv', sep=',', quote = F)



### check expression for signature genes
signature_genes_exp = list()
for(celltype in celltypes){
  sg = read.table(paste0('sc_data/', celltype, '_signature_genes.txt'), sep='\t', stringsAsFactors = F)
  sg = sg$V1
  sg_exp = df_batchs[sg, ]
  
  sg_exp = reshape2::melt(sg_exp)
  sg_exp$variable = as.character(sg_exp$variable)
  mean_value = sg_exp %>% group_by(variable) %>% summarise(mean = mean(value))
  signature_genes_exp[[celltype]] = mean_value[order(mean_value$variable), ]
}

check_celltype_signature = as.data.frame(signature_genes_exp)
mean_col = colnames(check_celltype_signature)[grep('mean', colnames(check_celltype_signature))]

df_plot = check_celltype_signature[,mean_col]
df_plot$cell_type = check_celltype_signature$tuft.variable
df_plot = reshape2::melt(df_plot, id_var = 'cell_type')

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
g = ggplot(aes(x = variable, y = value, fill = cell_type), data = df_plot) + 
  geom_bar(stat= 'identity', position='dodge') + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_fill_manual(values = getPalette(15)) + 
  xlab("Celltype of the signature genes") + 
  ylab("Mean value of genes in the signature geneset")



