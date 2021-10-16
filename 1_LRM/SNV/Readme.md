# Create Logistic Regression model (LRM) for SNV calling

The SNV model has not been used in the creation of GCAT|panel, due to bias decisions. The model consider all Deepvaraint findings as True positives, and Haplotype Caller and Strelka2 findings as false positives. Nevertheless, here we upload the code to execute the script of SNV model creation, using for training the GIAB sample (NA12878) and test the in-silico sample.  

To execute the script "SNVs_GIAB_model.R" the following R libraries are requiered:  
library(data.table)  
library(dplyr)  
library(xtable)  
library(stargazer)  
All files to train the model are included in GIAB folder, and files to test the model in insilico3 folder.  

The five columns of input files from GIAB and insilico3 folders corresponds with: Chromosome, Position, Reference allele, Alternate allele, and genotype. These files were generated with outputs from Haplotype Caller, Deepvariant and Strelka2 (Commands in GIAB folder).  
In order to evaluate correctly the precision and recall for each variant caller in real sample (NA12878) we just used the variants in conservative genome regions provided by Genome In A Bottle Consortium (GIAB)(https://www.nist.gov/programs-projects/genome-bottle).

After the execution of the SNVs_GIAB_model.R script different outputs were generated:  
1- glm_model.rds. The model to use in real samples (folder 2_merge_callers)  
2- Tables with precision recall and genotype error of three variant callers and model for NA12878 (train) and in-silico (test): table_snvs_model_giab_Lrm.csv and table_snvs_model_insilico3.csv  
3- The outputs from GIAB and insilico3 folders save intermidiate files.  
