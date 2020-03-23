# COVID-19: association between expression and AGE/SEX 
Script dependencies include the following R packages:
* dplyr
* data.table
* ggplot2
* sva

<<<<<<< HEAD
`download_data.sh`: download data from GTEX and Gencode for downstream analysis
=======
Approach I: <br/>
>>>>>>> 596cf2f6b407eab661b3c79982429367d860b24f
`tpm_age_sex.R`: controlling for covariates below by regressing them out:<br/>
                    # DTHHRDY  - Death Circumstances<br/>
                    # SMRIN  - RIN number <br/>
                    # SMTSISCH - Total Ischemic time <br/>
                    # SMEXNCRT - Exonic Rate <br/>
                 
<<<<<<< HEAD
                 
`tpm_age_sex_sva.R`: infer about surrogate variables (SVs) using SVA and regress the SVs out<br/>
                    When testing for AGE, keep AGE in SVA; when testing for SEX, keep SEX in SVA.<br/>
                    As an argument, takes in a tissue number (tissue numbers are the order listed by `cut -f7 GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | uniq`, 1-indexed).
                    
                    
                    
`extract_protein_coding_lincRNA_genes.sh`: Extract genes that are protein coding or lincRNA
`generate_tissue_wise_TPM.R`: Generate a matrix of genes by samples for each tissue. 
=======
Approach II:<br/>
`extract_protein_coding_lincRNA_genes.sh`: Extract genes that are protein coding or lincRNA

`generate_tissue_wise_TPM.R`: Generate a matrix of genes by samples for each tissue

`tpm_age_sex_sva.R`: infer about surrogate variables (SVs) using SVA and regress the SVs out<br/>
                    When testing for AGE, keep AGE in SVA; when testing for SEX, keep SEX in SVA.


`Plot_SV_corrected.R`: plot out the significant associations
>>>>>>> 596cf2f6b407eab661b3c79982429367d860b24f
