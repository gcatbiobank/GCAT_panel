# Translocation re-genotyping

Here we provided the code to re-genotype translocations after execute 2_merge_callers step.

To execute "regenotyping_translocation_bamok.py" script, we have to import the following libraries:  

from decimal import *  
import re  
import sys  
import pysam  
import gzip  
from bam_count_readsok import *  
import subprocess  

Besides we have to load the function bam_count_readsok.py script located in this folder. Is recommended install [Samtools](https://github.com/samtools/samtools) in this folder, otherwise you will have to modify the code of bam_count_readsok.py.

To execute the regenotyping_translocation_bamok.py script you have to include 7 arguments:

python regenotyping_translocation_bamok.py inputs/input_translocations_merge_callers.txt outputs/merge_callers_with_read_info.txt outputs/merge_callers_regenotype.txt inputs/insilico_small_set.bam insilico 100M outputs/

Where the first argument is the 2_merge_callers output, second argument is the 2_merge_callers output with additional column called "our_genotyper" with the new genotype:split read info:total coverage. Third argument is the 2_merge_callers output with the new genotype calculated. This file will be used for 3_merge_samples step. Fourth argument is the BAM file; Fifth argument is the sample name; Sixth argument is the CIGAR read length (in this case is 100M); and Seventh the output path.

With this script we can detect all reads that point to the translocation (just split reads). Additional information about filters applied in regenotyping_translocation_bamok.py, can be found in [supplementary materials](https://www.biorxiv.org/content/10.1101/2021.07.20.453041v1).
