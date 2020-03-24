# COVID-19: association between expression and AGE/SEX 
Script dependencies include the following R packages:
* dplyr
* data.table
* ggplot2
* sva

`download_data.sh`: download data from GTEX and Gencode for downstream analysis

Approach I: <br/>
Scripts to be run in the following order:<br/>
`tpm_age_sex.R`: controlling for covariates below by regressing them out:<br/>
* DTHHRDY  - Death Circumstances<br/>
* SMRIN  - RIN number <br/>
* SMTSISCH - Total Ischemic time <br/>
* SMEXNCRT - Exonic Rate <br/>


Approach II:<br/>
`extract_protein_coding_lincRNA_genes.sh`: Extract genes that are protein coding or lincRNA<br/>

`generate_tissue_wise_TPM.R`: Generate a matrix of genes by samples for each tissue<br/>

`tpm_age_sex_sva.R`: infer about surrogate variables (SVs) using SVA and regress the SVs out<br/>
                    When testing for AGE, keep AGE in SVA; when testing for SEX, keep SEX in SVA.<br/>
                    Pass in the tissue of interest's number (number can be found by:
                    `cut -f7 GTEx_data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | sort -u | grep -n "tissue_of_interest"` - 1 

`Plot_SV_corrected.R`: plot out the significant associations
