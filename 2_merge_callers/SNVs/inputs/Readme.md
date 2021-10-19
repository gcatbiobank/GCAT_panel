# Structure of input files  
1- The samples file is a list of samples with variant calling applied.  
2- The strategies.csv file consist in the name of callers used to perform the variant caller, and an integer, referring the strategy used to detect the SNVs. If any variant caller follow the same strategy to detect SNVs (Split-reads, Paired-end, Read-depth, combination) the number have to be the same.  
3- The inputs of each variant caller after VCF pre-processing.
  - Each file contains 9 columns, which are:
    chromosome, position, Reference allele, Alternative allel, Genotype, AD (Allelic Depth), DP (Read depth),GQ (Genome Quality),Highest PL (Phred-scaled genotype likelihoods) or GL (Genotype likelihood)

