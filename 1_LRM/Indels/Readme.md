# Create Logistic Regression Model (LRM) for Indel calling  
Here we upload the script (indels_GIAB_model.R) to create the LRM for Indel calling. The model is created using has input the Deepvariant, Haplotype Caller and Strelka2 algorithms.  
To execute the script "indels_GIAB_model.R" the following R libraries are requiered:  
library(data.table)  
library(dplyr)  
library(xtable)  
library(stargazer)  

All files to train the model are included in GIAB folder, and files to test the model in insilico3 folder.  
The five columns of input files from GIAB and insilico3 folders corresponds with: Chromosome, Position, Reference allele, Alternate allele, and genotype. These files were generated with outputs from Haplotype Caller, Deepvariant and Strelka2 (Commands in GIAB folder).  
In order to evaluate correctly the precision and recall for each variant caller in real sample (NA12878) we just used the variants in conservative genome regions provided by Genome In A Bottle Consortium (GIAB)(https://www.nist.gov/programs-projects/genome-bottle). Further description of filters applied and model creation can be found in supplementary materials of [original paper](https://www.biorxiv.org/content/10.1101/2021.07.20.453041v1).  

After the execution of the SNVs_GIAB_model.R script different outputs were generated:  
1- glm_model.rds. The model to use in real samples (folder 2_merge_callers)  
2- Tables with precision recall and genotype error of three variant callers and model for NA12878 (train) and in-silico (test): table_indel_model_giab_Lrm.csv and table_indel_model_insilico3.csv  
3- The outputs from GIAB and insilico3 folders save intermidiate files.  

To execute the script, with other input data, you have to change manually the input files. Besides, the script must be executed in the cloned folder.
