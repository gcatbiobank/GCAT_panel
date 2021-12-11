# Genotyping Duplications, Translocations and de novo insertions detected by Pindel

Due to their high genotype error produced by variant callers, here we deposited the code to re-genotype Duplications and Tranlocations detected after [section 2_merge_callers](https://github.com/gcatbiobank/GCAT_panel/tree/main/2_merge_callers/). However, Pindel does not provide the genotypes for Long insertion outputs. For this reason, after execute Pindel for all 785 samples, we appiled a custom genotyper for this particular SV type.  

## BAM file creation  

All  scripts provided requires a BAM file for each sample. Here we show an scheme followed to construct the BAM files for this project. Further documentation to create BAM files can be found in [Supplementary materials](https://www.biorxiv.org/content/10.1101/2021.07.20.453041v1)

![supplementary_fig5](https://user-images.githubusercontent.com/28949802/145678443-242c71cb-75b2-4acb-9ce5-6fd37a289560.jpg)

All information to re-genotype Duplications, Translocations and Pindel de novo insertions are described in each folder.
