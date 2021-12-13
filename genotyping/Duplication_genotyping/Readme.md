# Duplications re-genotyping  

To re-genotype all duplications obtained after [2_merge_callers](https://github.com/gcatbiobank/GCAT_panel/tree/main/2_merge_callers) step, we need to perform three steps.

1- Create the files to count reads in the BAM file  
2- Count the reads in the BAM file  
3- Re-genotype duplications for each sample  

The pipline is a combination of python and bash scripts.

## 1- Create the files to count reads in the BAM file

To execute "1_generate_files_to_count_reads.py" we need to import the following python2 modules:

import re  
import sys  
import decimal  

Then, as an input, we require the output obtained after the 2_merge_callers step (here provided in inputs folder under the name "insilico_merge_callers_input"). To execute "1_generate_files_to_count_reads.py" we need three arguments":

python generate_files_to_count_reads.py inputs/insilico_merge_callers_input insilico 1 outputs/

Were the first argument is the 2_merge_callers output for duplications, second sample name, third the chromosome and fourth the folder where save the output. The result is the file to obtain read information from BAM file  

## 2- Count the reads in the BAM file
