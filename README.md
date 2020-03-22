# COVID-19: association between expression and AGE/SEX 

Approach I:
`tpm_age_sex.R`: controlling for covariates below by regressing them out:<br/>
                    # DTHHRDY  - Death Circumstances<br/>
                    # SMRIN  - RIN number <br/>
                    # SMTSISCH - Total Ischemic time <br/>
                    # SMEXNCRT - Exonic Rate <br/>
                 
Approach II:
`extract_protein_coding_lincRNA_genes.sh`: Extract genes that are protein coding or lincRNA

`generate_tissue_wise_TPM.R`: Generate a matrix of genes by samples for each tissue

`tpm_age_sex_sva.R`: infer about surrogate variables (SVs) using SVA and regress the SVs out<br/>
                    When testing for AGE, keep AGE in SVA; when testing for SEX, keep SEX in SVA.

