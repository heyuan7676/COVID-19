#Download the data from GTEX:
#TPM data
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
gunzip  GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
TPM="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"

#Extract just the genes we are interested in.
grep -m1 "Description" $TPM > gene_tpm.gct ; sed -n 53469p $TPM >> gene_tpm.gct ; sed -n 51843p $TPM >> gene_tpm.gct

#Subject phenotypes:
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
#Sample annotations:
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
#Covariates data
wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_covariates.tar.gz
tar -xvf GTEx_Analysis_v8_eQTL_covariates.tar.gz
rm GTEx_Analysis_v8_eQTL_covariates.tar.gz



#Set up directories
mkdir -p GTEx_data
mv *.gct GTEx_data
mv *.txt GTEx_data
mv GTEx_Analysis_v8_eQTL_covariates ./GTEx_data/
mv GTEx_Analysis_v8_eQTL_expression_matrices GTEx_data/

mkdir -p results
