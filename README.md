**Version 1.0:** 12.22.17

This suite of modules is meant to ease the babysitting load for routinely processing genotype files to ADMIXTURE or imputation.  
* Requires Python 3.6.3, but will identify Python version and update if needed.  
* Starts from plink bed/bim/fam files. I'll eventually add other file formats.  

# Modules
## genoprocess
**Sees what version of Python you're running.**  
  Downloads Python 3.6.3 if you're on linux and exits telling you to download it if on Mac or Windows.  

**Asks the user if they would like to do the following:**  

genoprocess task | Short description | System requirements | Notes
----- | ----- | ----- | -----
1 | Download reference files and programs | Most any, some Mac or Linux | Some programs will not download if you aren't on Linux because they are Linux based  
2 | Update sex | Any |   
3 | Produce new dataset of people and SNPs with missince all rate < 10% | Any |   
4 | Run IBD matrix through Plink | Any |   
5 | Update family (FID) and individual IDs (IID) | Any |   
6 | Update maternal and paternal IDs | Any |   
7 | Harmonize with 1000G Phase 3 | Any |   
8 | Filter for extreme (+-3SD) heterozygosity values | Any |   
9 | Merge with 1000G Phase 3 | Any |   
10 | Prepare for ADMIXTURE with k = 3...9 | Mac or Linux | I will submit the job for you if you're on the Penn State ACI-B cluster  
11 | Run a phasing check and prepare files for phasing using SHAPEIT 2 | Mac or Linux | I will submit the job for you if you're on the Penn State ACI-B cluster  
12 | Prepares files for imputation using the Sanger Imputation Server | Mac or Linux |   
13 | Nothing | |   

## genodownload  
**Downloads the following:**  
  1. Python 3.6.3
  2. Plink 1.9
  3. 1000G Phase 3 VCF
  4. 1000G Phase 3 Hap/Legend/Sample
  5. GRCh37/hg19 1000G FASTA file
  6. Genotype Harmonizer
  7. pip
  8. snpflip
  9. shapeit
  10. vcftools
  11. bcftools
  12. htslib
  13. samtools

## genoqc  
Update sex  
Missing call rate  
Heterozygosity check  

## genorelatives  
Run IBD to identify relatives  
Update FID & IID information  
Update parental IDs  

## genoharmonize  
Harmonize with 1000G

## genomerge  
Merge with 1000G

## genoadmixture  
Prepare for admixture, submit job if on Penn Sate ACI-B cluster  

Should be done before admixture | genoprocess number | Module used  
------------------------------- | --------------------- | ------------  
Keep only SNPs and people with missing call rate < 10% | #3 | `genoqc.missing_call_rate()`  
Run IBD to identify relatives | #4 | `genorelatives.ibd()`
Update FID and/or IID information | #5 | `genorelatives.update_id()`  
Update parental IDs | #6 | `genorelatives.update_parental()`  
Harmonize with 1000G Phase3 | #7 | `genoharmonize.harmonize_with_1oooG()`  
Filter out individuals with extreme heterozygosity values | #8 | `genoqc.het()`  
Merge with 1000G | #9 | `genomerge.merge1000g()`  

## genophaseimpute
Checks data and prepares for phasing using shapeit.  
Checks phased data for imputation  

Should be done before phasing | genoprocess number | Module used  
----------------------------- | ------------------ | -----------
Genotypes must be on GRCh37/hg19 | User should check this | None
Keep only SNPs and people with missing call rate < 10% | #3 | `genoqc.missing_call_rate()`  
Filter out individuals with extreme heterozygosity values | #8 | `genoqc.het()`  
Filter out SNPs with MAF < 5% | #7 | `genoharmonize.harmonize_with_1000G()`
Filter out SNPs with HWE p-value < 0.05 | #7 | `genoharmonize.harmonize_with_1000G()`
Set haploid genotypes (male ChrX) as missing | User should check this if any warnings of het haploid genotypes appear | None
Check for gender mismatches | #2 | `genoqc.update_sex()`
Harmonize with 1000G Phase3 | #7 | `genoharmonize.harmonize_with_1000G()`
