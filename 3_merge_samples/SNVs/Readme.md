# Merge SNVs detected by sample in a single file  

After combine SNVs detected by Deepvariant, Haplotype Caller and Strelka2 per sample and predict the accuracy of their detections ([section 2_merge_callers](https://github.com/gcatbiobank/GCAT_panel/tree/main/2_merge_callers/SNVs)), the following R script combine the SNVs detected by sample in a single VCF file, providing the final quality of each variant detected in our cohort. A SNV is considered PASS if >= 0.5 of the SNV evaluated into the GCAT cohort was predicted as true positive, otherwise was NO PASS. Further SNV metrics were calculated, such as  Minor Allele Frequency (MAF), Population variantion (POPVAR), Allele Frequency (AF), among others. Detailed description can be found in [Supplementary materials](https://www.biorxiv.org/content/10.1101/2021.07.20.453041v1).  


