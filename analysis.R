### Compare results between the two methods

args <- commandArgs(TRUE)
Test_gene = args[1]
Test_gene_name = args[2]

#Test_gene = 'ENSG00000130234.10'
#Test_gene_name = 'ACE2'
Test_gene = 'ENSG00000184012.11'
Test_gene_name = 'TMPRSS2'

sva_result = read.table(paste0(outdir, 'Association_tests_',Test_gene, '_SVA.csv'),
                       sep= ',', header = T, stringsAsFactors = F)

lr_result = read.table(paste0(outdir, 'Association_tests_',Test_gene, '_LR.csv'),
                        sep= ',', header = T, stringsAsFactors = F)

compare_result = merge(sva_result, lr_result,
                      by = c("Tissue", "Gene", "Variable"))

compare_result$sva_log10P = -log10(compare_result$p.value)
compare_result$lr_log10P  = -log10(compare_result$pvalue)

# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


df = data.frame(x = compare_result$sva_log10P, y = compare_result$lr_log10P, z = compare_result$FDR.x < 0.05)
g = ggplot(aes(x=x, y = y), data = df) + 
    geom_point(aes(colour=factor(z))) + 
    geom_smooth(method = "lm", se = T, size = 0.5) + 
    xlab("Control for SVs") + 
    ylab("Control for known confounders") + 
    ggtitle(paste0('Compare -log10(p-value): ', Test_gene_name)) + 
    geom_text(x = 25, y = 23, label = lm_eqn(df), parse = TRUE, size = 3) + 
    labs(colour = "Significant hits") +
    scale_color_manual(values=c("#999999", "red")) +
    theme_bw()

png(paste0("results/Compare_LR_SVA_",Test_gene_name,".png"),res = 130, height = 400, width = 600)
print(g)
dev.off()




