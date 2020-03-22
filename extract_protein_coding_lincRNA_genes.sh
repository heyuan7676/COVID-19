#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH -p shared



### Extract protein-coding genes and lincRNA from the TPM gene matrix

datadir=/Users/Yuan/Desktop/GTEx_data/blog_post_data
cd  ${datadir}
cat gencode.v26.annotation.gene.txt  | grep protein_coding  | awk '{print $1}' >> pc.lc.genes
cat gencode.v26.annotation.gene.txt  | grep lincRNA | awk '{print $1}' >> pc.lc.genes

sort pc.lc.genes > temp
mv temp pc.lc.genes


split -l 5000 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct  
for fn in x*
do
	grep -w -F -f pc.lc.genes ${fn} > ${fn}.pc.lc
done

head -n3 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | sed '1d' > gene_pc_lc_tpm.txt
cat x*pc.lc >> gene_pc_lc_tpm.txt


rm x*
