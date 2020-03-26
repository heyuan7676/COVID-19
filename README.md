# COVID-19: association between expression and AGE/SEX 

## Script structure
###### Pipeline: `bash run.sh` <br/>
	`download_data.sh`: download data from GTEX and Gencode for downstream analysis


###### Script dependencies include the following R packages:<br/>
* dplyr
* data.table
* ggplot2
* sva
	
###### Approach I: Perform SVA to learn about the SVs, and test asscoations with AGE and SEX controlling for the SVs <br/>

	`extract_protein_coding_lincRNA_genes.sh`: Extract genes that are protein coding or lincRNA
	`generate_tissue_wise_TPM.R`: Generate a matrix of genes by samples for each tissue
	`SVA_compute_SV.R`: Infer about surrogate variables (SVs) using SVA
	`SVA_followedby_LR.R`: Perform association tests controlling for SVs


###### Approach II: Test asscoations with AGE and SEX controlling for known confounders <br/>
	`Rscript LR_confounders.R ENSG00000130234.10`
	`Rscript LR_confounders.R ENSG00000184012.11`
		* DTHHRDY  - Death Circumstances
		* SMRIN  - RIN number
		* SMTSISCH - Total Ischemic time
		* SMEXNCRT - Exonic Rate
