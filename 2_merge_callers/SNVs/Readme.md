# Merge Variant caller outputs of SNVs by sample  

In this section we combine by sample the outputs of Deepvariant, Haplotype Caller and Strelka2 without applying the Logistic Regression Model. Here we applied co-occurrence method, considering an SNV as True Positive (PASS), if it is detected at least by two variant callers. The strategy to merge callers by sample and filtering applied can be found in Supplementary material [sections 8.2 and 8.4](https://www.biorxiv.org/content/10.1101/2021.07.20.453041v1) respectively.
