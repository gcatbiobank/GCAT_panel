# Merge indel outputs by sample  

In this section we combine by sample the outputs of Deepvariant, Haplotype Caller and Strelka2, filtering the potential false positives with a Logistic Regression Model (created in [1_LRM/indels](https://github.com/gcatbiobank/GCAT_panel/tree/main/1_LRM/Indels) section). The model has been generated using the indels validated in NA12878 sample from Genome In A Bottle Consortium (GIAB) and tested with indels from in-silico sample. The strategy to merge callers by sample and filtering applied can be found in Supplementary material [sections 8.2 and 8.4](https://www.biorxiv.org/content/10.1101/2021.07.20.453041v1) respectively.  

The script "merge_callers_indels.R" is used to merge Deepvariant, Haplotype Caller and Strelka2 VCF outputs. This script requires the following libraries:  
library(data.table)  
library(dplyr)  
library(rlang)  
library(caret)  
library(e1071)  
library(R.utils)  

The script requires different inputs:  
#### 1- The strategies used by variant callers to detect the genetic variability (example found in input/strategies.csv)  
#### 2- A list of samples which we performed the variant calling (example found in input/samples)  
#### 3- The variant caller outputs for each sample (example of each variant caller output AFTER preprocessing can be found in inputs/Deepvaraint, GATK, Strelka2 respectively)  

The script is executed as follows:  
##### Rscript merge_callers_indels.R 10 1  

Where the first argument is the chromosome (10), and second argument (1) is the first sample of ./inputs/samples list file. It is important create the same name file as example presented in this demo, otherwise the merge_callers_indels.R code have to be modified. These strategies allow to parallelise the merge step by sample and chromosome. We recommended follow the same folder structure presented here to execute the R script, otherwise, the R code should be modified.  
