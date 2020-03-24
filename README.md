# COVID-19: association between expression and AGE/SEX 
Script dependencies include the following R packages:
* dplyr
* data.table
* ggplot2
* sva

`download_data.sh`: download data from GTEX and Gencode for downstream analysis <br/>

Approach I: Perform SVA to learn about the SVs, and test asscoations with AGE and SEX controlling for the SVs <br/>

`extract_protein_coding_lincRNA_genes.sh`: Extract genes that are protein coding or lincRNA <br/>
`generate_tissue_wise_TPM.R`: Generate a matrix of genes by samples for each tissue <br/>
`SVA_compute_SV.R`: Infer about surrogate variables (SVs) using SVA <br/>
`SVA_followedby_LR.R`: Perform association tests controlling for SVs <br/>


Pipeline:<br/>
`bash download_data.sh` <br/>
`bash extract_protein_coding_lincRNA_genes.sh`<br/>
`Rscript generate_tissue_wise_TPM.R`<br/>
`Rscript SVA_compute_SV.R`<br/>
`Rscript SVA_followedby_LR.R ENSG00000130234.10`<br/>
`Rscript SVA_followedby_LR.R ENSG00000184012.11`<br/>


Approach II: Test asscoations with AGE and SEX controlling for known confounders <br/>
`Rscript LR_confounders.R`

* DTHHRDY  - Death Circumstances<br/>
* SMRIN  - RIN number <br/>
* SMTSISCH - Total Ischemic time <br/>
* SMEXNCRT - Exonic Rate <br/>
