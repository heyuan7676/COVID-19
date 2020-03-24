#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p shared

ml R/3.5.1
ml parallel

## download data
bash download_data.sh

## data formatting
bash extract_protein_coding_lincRNA_genes.sh
Rscript generate_tissue_wise_TPM.R

## Learn SVs, test associations controlling for SVs

# Learn SVs - can be run in parallel
command_fn=run_SVA_compute_SV_parallel.txt
rm -f ${command_fn}
for i in {1..54}
do
	echo "Rscript SVA_compute_SV_parallel.R ${i}" >> ${command_fn}
done

parallel -j 10 :::: ${command_fn}

#Rscript SVA_compute_SV.R

# Linear regression controlling for SVs
Rscript SVA_followedby_LR.R ENSG00000130234.10
Rscript SVA_followedby_LR.R ENSG00000184012.11


## Linear regerssion controlling for known confounders
Rscript LR_confounders.R ENSG00000130234.10
Rscript LR_confounders.R ENSG00000184012.11
