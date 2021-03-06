# Merge indels detected in the GCAT cohort  

After combine indels detected by Deepvariant, Haplotype Caller and Strelka2 per sample and predict the accuracy of their detections ([section 2_merge_callers](https://github.com/gcatbiobank/GCAT_panel/tree/main/2_merge_callers/indels)), the following R script combine all indels detected in the cohort into a single VCF file, providing the final quality of each variant detected in our cohort. An indel is considered PASS if >= 0.5 of the indel evaluated into the GCAT cohort was predicted as true positive, otherwise was NO_PASS. Further indel metrics were calculated, such as  Minor Allele Frequency (MAF), Population variantion (POPVAR), Allele Frequency (AF), among others. Detailed description can be found in [Supplementary materials](https://www.biorxiv.org/content/10.1101/2021.07.20.453041v1).  

Before to execute the samples_merge_indels.R script you have to full fill the following requirements:  
1- Install the following libraries:  
library(data.table)  
library(dplyr)  
library(R.utils)  

2- To create subfolders for each chromosome in the outputs folder ( i.e: chr_10) and empty_batches folder too.  

3- Before to execute samples_merge_indels.R script, you must be sure that merge_callers_indels.R is executed for each chromosome and sample. samples_merge_indels.R script cut by batches of 50000 the chromosome 1 to X, for this reason is recommended obtain all chromosomes merged by caller. Otherwise, the script will breaks, because will miss some chromosome.  

4- To execute the R script in cloned folder like this:

Rscript samples_merge_indels.R 33608 2 path_input files

First argument is a batch to merge, in this case, the 33608 batch corresponds to chromosome 10 at position ranges between 150001-200000.  
Second argument is the sample number to merge (this is a cohort of two individuals).  
Third is the path of outputs obtained from merge_callers step.

The file length_chromosomes_for_merge.txt contains the length of chromosomes from reference genome hs37d5.
The file samplesok is a list of all samples to merge in the VCF file.
The header includes the description of INFO column from VCF file.

The result is a VCF file with all samples merged, including their genotypes and relevant population information.
