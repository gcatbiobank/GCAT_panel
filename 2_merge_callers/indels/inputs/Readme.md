# Structure of input files  

1- The samples file is a list of samples with variant calling applied.  
2- The strategies.csv file is the name of callers used to perform the variant caller, and an integer, referring the strategy used to detect the indels. If any variant caller follow the same strategy to detect indels (Split-reads, Paired-end, Read-depth, combination) the number have to be the same.  
3- The inputs of each variant caller (Deepvariant, Haplotype Caller and Strelka2) after VCF pre-processing.  

Each file contains 9 columns:  
Chromosome, Position, Reference allele, Alternative allel, Genotype, AD (Allelic Depth), DP (Read depth), GQ (Genome Quality),Highest PL (Phred-scaled genotype likelihoods) or GL (Genotype likelihood)  
