# Julie's Processing script for Shriver Lab genotype data
# Date 9.5.17

# WHAT THIS DOES #
# -Run IBD matrix through plink
# -Update family and individual IDs
# -Update maternal and paternal IDs
# -Update sex
# -Using 1000 Genomes as a reference (this part based off Perl script by W. Rayner, 2015, wrayner@well.ox.ac.uk)
#   -Removes SNPs with MAF < 5% in study dataset
#   -Removes SNPs not in 1000 Genomes
#   -Removes all A/T G/C SNPs with MAF > 40% in the reference data set
#   -Removes all SNPs with an AF difference >0.2, between reference and dataset frequency file, frequency file is
#       expected to be a plink frequency file with the same number of SNPs as the bim file
#   -Removes duplicates that may be introduced with the position update
#   -Removes indels
# -Merges your data with 1000 Genomes
# -Runs ADMIXTURE with k = 3...9
# -Phases using SHAPEIT2

# REQUIREMENTS #
# You must have plink 1.9 in the same directory as your genotype files and this script, and it should be called 'plink'.

# Getting the needed modules.
import os

to_do = input('What would you like to do?\n'
              '1) Run an IBD analysis to identify relatives. All you need are plink bed/bim/fam files.\n'
              '2) Update FID or IID information. You need a file with the following information Old FID, Old IID, '
              'New FID, New IID.\n'
              '3) Update parental IDs. You need a file with FID, IID, Paternal IID, and Maternal IID.\n'
              '4) Update sex. You need a file with FID, IID, Sex (M = 1, F = 2, Unknown = 0)\n'
              '5) Merge with 1000 Genomes Phase 3\n'
              '6) Nothing\n'
              'Please enter a number (i.e. 2): ')

if to_do == '1':

    # Identity-by-descent in Plink
    # This part of the script will prune for LD, calculate IBD, and exclude individuals who have IBD < 0.2
    # The IBD results will have .genome appended to your file name. I have also included a line to convert the IBD results
    #   from whitespace to tab delimited. This will have .tab.genome appended to your filename.

    # Important values of Pi-hat
    #   -First-degree relative = 0.5
    #   -Second-degree relative = 0.25
    #   -Third-degree relative = 0.125
    #   -Fourth-degree relative = 0.0625
    #   -Fifth-degree relative = 0.03125

    geno_name = input('Please enter the name of the genotype files (without bed/bim/fam extension: ')
    print("Your IBD results in a tab delimited file will have the name " + geno_name + ".tab.genome. You should use "
                                                                                       "this file to investigate your "
                                                                                       "relatives and possibly update "
                                                                                       "the FID and IIDs in your file.")
    os.system('plink --bfile ' + geno_name + ' --indep 50 5 2 --out ' + geno_name)
    os.system('plink --bfile ' + geno_name + ' --exclude ' + geno_name + '.prune.out --genome --min 0.2 --out ' + geno_name)
    os.system('sed -r "s/\s+/\t/g" ' + geno_name + '.genome > ' + geno_name + '.tab.genome')
        # Comment out this line if you prefer whitespace delimited files
    print("Analysis finished")
    # Now your job is to use the .tab.genome file to investigate relatives and possibly update FID/IID and parents.

elif to_do == '2':
    # File for updating FID should have four fields
    #  1) Old FID
    #  2) Old IID
    #  3) New FID
    #  4) New IID

    geno_name = input('Please enter the name of your genotype files (without bed/bim/fam extension): ')
    update_fid_filename = input('Please enter the name of your text file for updating FID or IID (with file extension): ')
    print("Your genotype files with the FID updated will have the name " + geno_name + "_FIDUpdated")
    os.system('plink --bfile ' + geno_name + ' --update-ids ' + update_fid_filename + ' --make-bed --out ' + geno_name +
            '_FIDUpdated')

elif to_do == '3':
    # File for updating parents should have four fields:
    #   1) FID
    #   2) IID
    #   3) New paternal IID
    #   4) New maternal IID

    geno_name = input('Please enter the name of your genotype files (without bed/bim/fam extension). Remember, if you '
                      'just updated FIDs, then your genotype name should have FIDUpdated at the end of it.: ')
    update_parents_filename = input('Please enter the name of your text file for updating parents (with file extension): ')
    print("Your genotype files with parents updated will have the name " + geno_name + "_ParentsUpdated")
    os.system('plink --bfile ' + geno_name + ' --update-parents ' + update_parents_filename + ' --make-bed --out ' +
              geno_name + '_ParentsUpdated')

elif to_do == '4':
    # File for updating sex should have:
    #   1) FID
    #   2) IID
    #   3) Sex (1 = M, 2 = F, 0 = missing)

    geno_name = input('Please enter the name of your genotype files (without bed/bim/fam extension). Remember, if you '
                      'just updated FIDs, then your genotype name should have _FIDUpdated at the end of it. If you just '
                      'updated parents, then your genotype name should have _ParentsUpdated at the end of it.: ')
    update_sex_filename = input('Please enter the name of your text file for updating sex (with file extension): ')
    print("Your genotype files with sex updated will have the name " + geno_name + "_SexUpdated")
    os.system('plink --bfile ' + geno_name + ' --update-sex ' + update_sex_filename + ' --make-bed --out ' + geno_name
              + '_SexUpdated')

elif to_do == '5':
    # Remove SNPs with MAF < 0.05, because they are often not imputed well
    geno_name = input('Please enter the name of the genotype files (without bed/bim/fam extension: ')
    os.system('plink --bfile ' + geno_name + ' -maf 0.05 --make-bed --out ' + geno_name + '_MAF0.05')

    # Remove SNPs that are not in 1000 G Phase 3
    # This is where I got the reference files: https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html
    import pandas as pd
    import numpy as np

    bim_file = pd.read_csv(geno_name + '_MAF0.05.bim', sep="\t", header=None, usecols = [0,1,3,4,5],
                           names=['chr', 'dataset_id', 'position', 'dataset_A1', 'dataset_A2'])

    snps_by_chr = ['chr%d_snps' % x for x in range(1, 23)]
    legend_file_names = ['1000GP_Phase3_chr%d.legend' % x for x in range(1, 23)]
    maf_AT_filter = ['chr%d_maf_AT_filter' % x for x in range(1, 23)]
    maf_TA_filter = ['chr%d_maf_AT_filter' % x for x in range(1, 23)]
    maf_GC_filter = ['chr%d_maf_AT_filter' % x for x in range(1, 23)]
    maf_CG_filter = ['chr%d_maf_AT_filter' % x for x in range(1, 23)]
    maf_per_chr_filter = ['chr%d_maf_filter' % x for x in range(1,23)]

    # Match the position in 1000 genomes with the position in our genotype file, on a chromosome by chromosome basis.
    # Removes all A/T G/C SNPs with MAF > 40% in the reference data set
    for i in range(0, len(snps_by_chr)):
        current_legend_file = pd.read_csv(legend_file_names[i], sep=" ", header = 0,
                                          dtype={'id': str, 'position': int, 'a0': str,'a1': str, 'TYPE': str,
                                                 'AFR': float, 'AMR': float, 'EAS': float, 'EUR': float, 'SAS': float,
                                                 'ALL': float})
        current_legend_file.rename(columns={'id': 'reference_id', 'a0': 'REF', 'a1': 'ALT', 'TYPE': 'type'}, inplace=True)
        print('Successfully read in chr' + str(i + 1) + ' legend file')
        current_legend_file['reference_MAF'] = np.where(current_legend_file['ALL'] <= 0.5, current_legend_file['ALL'],
                                                        1 - current_legend_file['ALL'])
        current_legend_file = current_legend_file[current_legend_file.type == 'Biallelic_SNP']

        snps_by_chr[i] = pd.merge(left=bim_file.loc[bim_file['chr'] == i + 1], right=current_legend_file, how='inner', on='position')
        print('chr' + str(i + 1) + ' overlap with 1000G complete')
        maf_AT_filter[i] = snps_by_chr[i].loc[((snps_by_chr[i]['REF'] == 'A') & (snps_by_chr[i]['ALT'] == 'T')) &
                                              (snps_by_chr[i]['reference_MAF'] >= 0.4)]
        maf_TA_filter[i] = snps_by_chr[i].loc[((snps_by_chr[i]['REF'] == 'T') & (snps_by_chr[i]['ALT'] == 'A')) &
                                              (snps_by_chr[i]['reference_MAF'] >= 0.4)]
        maf_GC_filter[i] = snps_by_chr[i].loc[((snps_by_chr[i]['REF'] == 'G') & (snps_by_chr[i]['ALT'] == 'C')) &
                                              (snps_by_chr[i]['reference_MAF'] >= 0.4)]
        maf_CG_filter[i] = snps_by_chr[i].loc[((snps_by_chr[i]['REF'] == 'C') & (snps_by_chr[i]['ALT'] == 'G')) &
                                              (snps_by_chr[i]['reference_MAF'] >= 0.4)]
        maf_per_chr_filter[i] = pd.concat([maf_AT_filter[i], maf_TA_filter[i], maf_GC_filter[i], maf_CG_filter[i]])
        print('Filtered out A/T G/C SNPs with MAF > 40% in chr' + str(i + 1))

    # Chromosome X
    # Plink codes the pseudoautosomal region as 23, for now I won't include the non-pseudoautosomal regions since they
    #   need to be treated differently.
    current_legend_file = pd.read_csv('1000GP_Phase3_chrX_PAR1.legend', sep = " ",
                                      dtype={'id': str, 'position': int, 'a0': str, 'a1': str, 'TYPE': str, 'AFR': float,
                                             'AMR': float, 'EAS': float, 'EUR': float, 'SAS': float, 'ALL': float})
    print('Successfully read in chrX pseudoautosomal region 1 legend file')
    current_legend_file.rename(
        columns={'id': 'reference_id', 'a0': 'REF', 'a1': 'ALT', 'TYPE': 'type'}, inplace=True)
    current_legend_file['reference_MAF'] = np.where(current_legend_file['ALL'] <= 0.5, current_legend_file['ALL'],
                                                    1 - current_legend_file['ALL'])
    current_legend_file = current_legend_file[current_legend_file.type == 'Biallelic_SNP']
    chrX_PAR1 = pd.merge(left=bim_file.loc[bim_file['chr'] == 23], right = current_legend_file, how='inner', on='position')
    print('chrX pseudoautosomal region 1 overlap with 1000G done')

    chrX_PAR1_AT_filter = chrX_PAR1.loc[
        ((chrX_PAR1['REF'] == 'A') & (chrX_PAR1['ALT'] == 'T')) & (chrX_PAR1['reference_MAF'] >= 0.4)]
    chrX_PAR1_TA_filter = chrX_PAR1.loc[
        ((chrX_PAR1['REF'] == 'T') & (chrX_PAR1['ALT'] == 'A')) & (chrX_PAR1['reference_MAF'] >= 0.4)]
    chrX_PAR1_GC_filter = chrX_PAR1.loc[
        ((chrX_PAR1['REF'] == 'G') & (chrX_PAR1['ALT'] == 'C')) & (chrX_PAR1['reference_MAF'] >= 0.4)]
    chrX_PAR1_CG_filter = chrX_PAR1.loc[
        ((chrX_PAR1['REF'] == 'C') & (chrX_PAR1['ALT'] == 'G')) & (chrX_PAR1['reference_MAF'] >= 0.4)]
    chrX_PAR1_maf_filter = pd.concat([chrX_PAR1_AT_filter, chrX_PAR1_TA_filter, chrX_PAR1_GC_filter, chrX_PAR1_CG_filter])
    print('Filtered out A/T G/C SNPs by MAF > 40% for chrX pseudoautosomal region 1')

    current_legend_file = pd.read_csv('1000GP_Phase3_chrX_PAR2.legend', sep=" ", header=0,
                                      dtype={'id': str, 'position': int, 'a0': str, 'a1': str, 'TYPE': str, 'AFR': float,
                                             'AMR': float, 'EAS': float, 'EUR': float, 'SAS': float, 'ALL': float})
    current_legend_file.rename(columns={'id': 'reference_id', 'a0': 'REF', 'a1': 'ALT', 'TYPE': 'type'}, inplace=True)
    print('Successfully read in chrX pseudoautosomal region 2 legend file')
    current_legend_file['reference_MAF'] = np.where(current_legend_file['ALL'] <= 0.5, current_legend_file['ALL'],
                                                    1 - current_legend_file['ALL'])
    current_legend_file = current_legend_file[current_legend_file.type == 'Biallelic_SNP']
    chrX_PAR2 = pd.merge(left=bim_file.loc[bim_file['chr'] == 23], right=current_legend_file, how='inner', on='position')
    print('chrX pseudoautosomal region 2 overlap with 1000G done')
    chrX_PAR2_AT_filter = chrX_PAR2.loc[
        ((chrX_PAR2['REF'] == 'A') & (chrX_PAR2['ALT'] == 'T')) & (chrX_PAR2['reference_MAF'] >= 0.4)]
    chrX_PAR2_TA_filter = chrX_PAR2.loc[
        ((chrX_PAR2['REF'] == 'T') & (chrX_PAR2['ALT'] == 'A')) & (chrX_PAR2['reference_MAF'] >= 0.4)]
    chrX_PAR2_GC_filter = chrX_PAR2.loc[
        ((chrX_PAR2['REF'] == 'G') & (chrX_PAR2['ALT'] == 'C')) & (chrX_PAR2['reference_MAF'] >= 0.4)]
    chrX_PAR2_CG_filter = chrX_PAR2.loc[
        ((chrX_PAR2['REF'] == 'C') & (chrX_PAR2['ALT'] == 'G')) & (chrX_PAR2['reference_MAF'] >= 0.4)]
    chrX_PAR2_maf_filter = pd.concat([chrX_PAR2_AT_filter, chrX_PAR2_TA_filter, chrX_PAR2_GC_filter, chrX_PAR2_CG_filter])
    print('Filtered out A/T G/C SNPs by MAF > 40% for chrX pseudoautosomal region 2')

    # Remove A/T G/C SNPs with MAF > 40% from list of snps in our dataset with matches in 1000 Genomes. 
    overlap_with_1000G = pd.concat([snps_by_chr[0], snps_by_chr[1], snps_by_chr[2], snps_by_chr[3], snps_by_chr[4],
                 snps_by_chr[5], snps_by_chr[6], snps_by_chr[7], snps_by_chr[8], snps_by_chr[9],
                 snps_by_chr[10], snps_by_chr[11], snps_by_chr[12], snps_by_chr[13], snps_by_chr[14],
                 snps_by_chr[15], snps_by_chr[16], snps_by_chr[17], snps_by_chr[18], snps_by_chr[19],
                 snps_by_chr[20], snps_by_chr[21], chrX_PAR1, chrX_PAR2])
    palindromic_MAF_filter = pd.concat([maf_per_chr_filter[0], maf_per_chr_filter[1], maf_per_chr_filter[2], maf_per_chr_filter[3],
                                  maf_per_chr_filter[4],maf_per_chr_filter[5], maf_per_chr_filter[6], maf_per_chr_filter[7],
                                  maf_per_chr_filter[8], maf_per_chr_filter[9], maf_per_chr_filter[10], maf_per_chr_filter[11],
                                  maf_per_chr_filter[12], maf_per_chr_filter[13], maf_per_chr_filter[14], maf_per_chr_filter[15],
                                  maf_per_chr_filter[16], maf_per_chr_filter[17], maf_per_chr_filter[18], maf_per_chr_filter[19],
                                  maf_per_chr_filter[20], maf_per_chr_filter[21], chrX_PAR1_maf_filter, chrX_PAR2_maf_filter])
    common_snps = pd.merge(left = overlap_with_1000G, right = palindromic_MAF_filter, how = 'inner')
    snps_to_keep = overlap_with_1000G[~overlap_with_1000G.reference_id.isin(common_snps.reference_id)]
    snps_to_keep = snps_to_keep.drop_duplicates(subset = 'dataset_id', keep = 'first', inplace = False)
    snps_to_keep.to_csv('SNPs_to_keep_info.txt', sep='\t', header=True, index=False)
    snps_to_keep['dataset_id'].to_csv('SNPs_to_keep.txt', sep = '\t', header = False, index = False)
    snps_to_keep['ChangeAlleleOrder'] = np.where((snps_to_keep['dataset_A1'] == snps_to_keep['REF']) &
                                                 (snps_to_keep['dataset_A2'] == snps_to_keep['ALT']),
                                                 'Change', 'KeepForNow')
    SNPsToChange = snps_to_keep[snps_to_keep['ChangeAlleleOrder'] == 'Change']
    SNPsToChange[['dataset_id', 'ALT']].to_csv('ForceA1Alleles.txt', sep='\t', header=False, index=False)
    os.system(
        'plink --bfile ' + geno_name + '_MAF0.05 --extract SNPs_to_keep.txt --reference-allele ForceA1Alleles.txt --freq --make-bed --out '
        + geno_name + '_MAF0.05_FilteredSNPs_ChangeA1Test')
    snps_to_keep.to_csv('Info_on_SNPs_to_keep.txt', sep='\t', header=True, index=False)

    #Need to check strand as well.

elif to_do == '6':
    print("You go, couch potato")

else:
    print("Please enter a number 1-6.")



