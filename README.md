---
title: "Transmission analysis"
author: "Karim"
  
---
# **Studying malaria transmission dynamic in the Senegambia area**
Given the whole genome variant call file (VCF) that contains SNPs loci genotyped across several samples which belong to different population, we proposed an approach for the study of malaria transmission dynamic in the geographical space covered by the different sampling locations.    
In this study, isolates are from 9 Gambian and 2 Senegalese locations, we aimed at determining the transmission dynamic between Senegal and Gambia and within Gambia. Our method covers the following points:   
&nbsp;&nbsp;&nbsp;&nbsp;**1.** Data QC. This involves the filtration of SNPs and isolates based on their percentage of missing genotypes, and the filtration of SNPs based on the minor allele frequency (MAF) provided by the user.   
&nbsp;&nbsp;&nbsp;&nbsp;**2.** Estimation of the isolates's within-host diversity    
&nbsp;&nbsp;&nbsp;&nbsp;**3.** Construction of the input files for relatedness inference. Here the user can choose to recode the mixed genotypes using 3 methods: `raw` for considering the mixed allele as a third allele and recode it as 2, `minor_David` for recoding based on the minor allele on the given SNP, `minor_Karim` for recoding based on the minor allele obtained from the number of reads supporting each allele on the given site for the given sample.    
&nbsp;&nbsp;&nbsp;&nbsp;**4.** Inferrence of relatedness bet isolates from every pairwise population present in the data. This is based on the model developed by Aimee R. Taylor and co-authors. https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009101. We recommend that you modify the `runSim.sh` script befroe renning this step so that the instructions can suit your HPC.     
&nbsp;&nbsp;&nbsp;&nbsp;**5.** Concatenate the relatedness files created from `4`      
&nbsp;&nbsp;&nbsp;&nbsp;**6.** Summarize the relatedness data for each pair of population     
&nbsp;&nbsp;&nbsp;&nbsp;**7.** Combine the different pairwise-population relatedness files into one single file    
&nbsp;&nbsp;&nbsp;&nbsp;**8.** We propose an analysis approach


***

###  Depencies   
    1. softwares: _tabix, R_  
    2. R packages: _data.table, tictoc, doParallel, ggplot2, tidyr, dplyr, moimix, SeqArray, ggplot2, statip, Rcpp, Rfast, yaml, fst_     
    
***  


### ** cloning into the repository and preparing for the analysis**     
```{bash eval = FALSE}
git clone https://gitlab.internal.sanger.ac.uk/km28/transmission_analysis.git          
cd transmission_analysis          
mkdir Results/ Data/       
cp your/input/vcf/file.vcf.gz Data/    
```


***  

See the README.Rmd file for mor details.
