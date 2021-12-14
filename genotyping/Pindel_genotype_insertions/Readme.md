# Genotype Long insertions detected by Pindel  

The genotyper for the de novo insertions detected by Pindel has been done in Python2. To execute "regenotyping_pindel_bam.py" script you need to import the following modules:

from decimal import *  
import re  
import sys  
import pysam

Then, the script requires 4 arguments:

**python regenotyping_pindel_bam.py input/insilico_pindel_LI.vcf output/insilico_pindel_LI_genotype.vcf input/insilico_small.bam 100M**

First argument is the Pindel VCF, second argument is the script output, third argument is the BAM file and fourth the read length (CIGAR nomenclature).
