# Julie's Processing script for Shriver Lab genotype data
# Date 9.5.17

# WHAT THIS DOES #
# -Run ibs matrix through plink and use R to spit out list of relatives as well as relatedness histogram
# -Update family and individual IDs
# -Update maternal and paternal IDs
# -Update sex
# -Using 1000 Genomes as a reference (this part based off Perl script by W. Rayner, 2015, wrayner@well.ox.ac.uk
#   -Removes SNPs not in 1000 Genomes
#   -Removes all A/T G/C SNPs with MAF > 40% in the reference data set
#   -Removes all SNPs with an AF difference >0.2, between reference and dataset frequency file, frequency file is
#    expected to be a plink frequency file with the same number of SNPs as the bim file
#   -Removes duplicates that may be introduced with the position update
#   -Removes indels
# -Merges your data with 1000 Genomes
# -Runs ADMIXTURE with k = 3...9
# -Phases using SHAPEIT2

# REQUIREMENTS #
# You must have R installed on your machine to run this script.
# You must have plink 1.9 in the same directory as your genotype files and this script, and it should be called 'plink'.
# Using this version of the script requires that you have rpy2 installed on your machine. If you have Anaconda, then you
#   can install rpy2 using: conda install -c r rpy2
#   Otherwise, install using directions at https://rpy2.readthedocs.io/en/version_2.8.x/index.html

# Getting the needed modules.
import os

# Identity-by-descent in Plink
geno_name = input('Enter the name of the genotype files: ')
os.system('plink --bfile ' + geno_name + ' --indep 50 5 2 --out ' + geno_name)
os.system('plink --bfile ' + geno_name + ' --exclude ' + geno_name + '.prune.out --genome --min 0.2 --out ' + geno_name)
os.system('sed -r "s/\s+/\t/g" ' + geno_name + '.genome > ' + geno_name + '.tab.genome')

# Important values of Pi-hat
#   -First-degree relative = 0.5
#   -Second-degree relative = 0.25
#   -Third-degree relative = 0.125
#   -Fourth-degree relative = 0.0625
#   -Fifth-degree relative = 0.03125

# Now your job is to use the .tab.genome file to investigate relatives and create files with the FIDs and IIDs