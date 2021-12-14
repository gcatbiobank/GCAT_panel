# Duplications re-genotyping  

To re-genotype all duplications obtained after [2_merge_callers](https://github.com/gcatbiobank/GCAT_panel/tree/main/2_merge_callers), we need to perform three steps.

1- Selecting the coordinates to count reads in the BAM file  
2- Determining the median coverage of duplications using the BAM file  
3- Re-genotype duplications for each sample  

The pipline is a combination of python and bash scripts.

## 1- Selecting the coordinates to count reads in the BAM file

The first python script consist to select the coordinates to find the median coverage around duplications.

To execute "1_generate_files_to_count_reads.py" we need to import the following python2 modules:

import re  
import sys  
import decimal  

As input, we require the output obtained after the 2_merge_callers step (here provided in the inputs folder under the name "insilico_merge_callers_input"). To execute "1_generate_files_to_count_reads.py" we need three arguments":

**python 1_generate_files_to_count_reads.py inputs/insilico_merge_callers_input insilico 10 outputs/**

Where the first argument is the 2_merge_callers output for duplications, second is the sample name, third is the chromosome and fourth is the folder where save the output. The result is a file with following information separated by 7 columns:  

Chromosome, Start_dup, End_dup, Size_dup, distance between current and previous duplication, Breakposition error, and genotype obtanied after merge_callers step.

## 2- Determining the median coverage of duplications using the BAM file

To determine the duplications genotype state (0/1 or 1/1), we need the median coverage of that region. Determining the median coverage of a duplication is a challenge due to two main reasons: 1) the sequencing coverage of a duplication is at least the double of sample coverage and 2) the sequencing coverage fluctuates around the genome. For this reason, the objective of the second bash script consist to find the median coverage of a region where the duplication is detected.   

To execute 2_check_cx.sh script, we need the get_cx.sh script present in this folder. Besides, we need to install [Samtools](https://github.com/samtools/samtools) in this folder, or change the path where the samtools is installed in 2_check_cx.sh and get_cx.sh scripts. Finally, we need three arguments:

./2_check_cx.sh outputs/insilico_chr10.txt inputs/insilico3_test.bam insilico

First argument is the output obtained in the previous step, second argument is the BAM file, and finally, the sample name.

The result, is the same file as previous step, with the median coverage for each duplication.  

## 3- Re-genotype duplications for each sample

Finally, we have to execute 3_regenotyping_dups.py script, to re-genotype all duplications for each sample. This script requires the following python modules:  

from decimal import *  
import re  
import sys  
import pysam  
import gzip  

To execute 3_regenotyping_dups python script, we require five arguments:

python 3_regenotyping_dups.py inputs/insilico3_test.bam inputs/insilico_merge_callers_input outputs/insilico_10_readinfo.txt outputs/insilico_merge_callers_input_10_regentoype 100M

First is the BAM file, second is the 2_merge_callers output, third is the previous step output, fourth is the output and fifith the read lenght all matched (CIGAR).

The resulting output is the 2_merge_callers adding three columns. Our genotype evaluated, the median coverage and split reads from duplication breakpoint.

Further filter descriptions of the duplication genotyper can be found in [supplementary materials](https://www.biorxiv.org/content/10.1101/2021.07.20.453041v1).
