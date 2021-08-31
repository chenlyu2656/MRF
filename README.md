# Methylation Random Field Method
The Methylation Random Field (MRF) method is developed to detect methylation quantitative trait loci (mQTLs) by testing the association between the methylation level of a CpG site and a set of genetic variants within a genomic region.

The proposed MRF has two major advantages:
1) It uses a beta distribution to characterize the bimodal and interval properties of the methylation trait at a CpG site;
2) It considers multiple common and rare genetic variants within a genomic region to identify mQTLs. 

# MRF Function
Available [here](./R/MRF_fun.R)

# Tutorial
Sample data availability
- [Methylation trait](./Example/Trait.txt)
- [Genotypes](./Example/Genotype.txt)
- [Covariates](./Example/Covariate.txt)

Sample code availability
- [Example](./Example/Example.R)

Reference:
Lyu, Chen, Manyan Huang, Nianjun Liu, Zhongxue Chen, Philip J. Lupo, Benjamin Tycko, John S. Witte, Charlotte A. Hobbs, and Ming Li. "Detecting methylation quantitative trait loci using a methylation random field method." Briefings in Bioinformatics (2021)
