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
                      'just updated FIDs, then your genotype name should have "_FIDUpdated" at the end of it.: ')
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
                      'just updated FIDs, then your genotype name should have "_FIDUpdated" at the end of it. If you just '
                      'updated parents, then your genotype name should have "_ParentsUpdated" at the end of it.: ')
    update_sex_filename = input('Please enter the name of your text file for updating sex (with file extension): ')
    print("Your genotype files with sex updated will have the name " + geno_name + "_SexUpdated")
    os.system('plink --bfile ' + geno_name + ' --update-sex ' + update_sex_filename + ' --make-bed --out ' + geno_name
              + '_SexUpdated')

elif to_do == '5':
    # Removes SNPs not in the reference 1000G Phase 3
    # Removes SNPs with HWE p-value < 0.01
    # Removes SNPs with MAF < 5%
    # Update IDs to match the reference
    # Updates the reference allele to match 1000G
    # Outputs new files per chromosome.

    import pandas as pd
    import numpy as np

    geno_name = input('Please enter the name of the genotype files (without bed/bim/fam extension: ')

    vcf_file_names = ['ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes' % x for x in range(1,23)]
    freq_file_names = [geno_name + '_chr%d_harmonized.frq' % x for x in range(1,23)]
    bim_file_names = [geno_name + '_chr%d_harmonized.bim' % x for x in range(1,23)]
    legend_file_names = ['1000GP_Phase3_chr%d.legend.gz' % x for x in range(1,23)]
    harmonized_geno_names = [geno_name + '_chr%d_harmonized' % x for x in range(1,23)]
    AF_diff_removed_by_chr = ['chr%d_SNPsRemoved_AFDiff' % x for x in range(1,23)]
    final_snp_lists_by_chr = ['chr%d_SNPsKept' % x for x in range(1,23)]
    final_snp_list_names = ['chr%d_SNPsKept.txt' % x for x in range(1,23)]
    AF_checked_names = [geno_name + '_chr%dharmonized_AFChecked' % x for x in range(1,23)]

    for i in range(0, len(vcf_file_names)):
        os.system('java -Xmx1g -jar GenotypeHarmonizer.jar $* '
                  '--input ' + geno_name +
                  ' --ref ' + vcf_file_names[i] + '.vcf.gz'
                  ' --refType VCF '
                  '--chrFilter ' + str(i + 1) +
                  ' --hweFilter 0.01 '
                  '--mafFilter 0.05 '
                  '--update-id '
                  '--debug '
                  '--mafAlign 0 '
                  '--update-reference-allele '
                  '--outputType PLINK_BED '
                  '--output ' + harmonized_geno_names[i])

    # Since we are going to use a global reference population, removes all SNPs with an AF difference > 0.2
    #  between study dataset and all superpopulation allele frequencies. IF within 0.2 of any superpopulation frequency,
    #  keep variant. Frequency file is expected to be a Plink frequency file with the same number of variants as the bim file.

        os.system('plink --bfile ' + harmonized_geno_names[i] + ' --freq --out ' + harmonized_geno_names[i])

        freq_file = pd.read_csv(freq_file_names[i], sep='\s+', header = 0, usecols = [0,1,2,3,4],
                                dtype = {'CHR': int, 'SNP': str, 'A1': str, 'A2': str, 'MAF':float})
        freq_file.rename(columns={'A1': 'dataset_a1', 'A2': 'dataset_a2', 'MAF': 'dataset_a1_frq'}, inplace=True)
        freq_file['dataset_a2_frq'] = 1-freq_file['dataset_a1_frq']

        bim_file = pd.read_csv(bim_file_names[i], sep='\s+', header=None, usecols=[0, 1, 3],
                               names=['CHR', 'SNP', 'position'])

        freq_file_with_position = pd.merge(left = freq_file, right = bim_file, how = 'inner', on=['CHR','SNP'])

        legend_file = pd.read_csv(legend_file_names[i], compression="gzip", sep=" ", header = 0,
                                  dtype={'id': str, 'position': int, 'a0': str,'a1': str, 'TYPE': str, 'AFR': float,
                                         'AMR': float, 'EAS': float, 'EUR': float, 'SAS': float, 'ALL': float})

        legend_file.rename(columns={'id': 'reference_id', 'a0': 'reference_a0', 'a1': 'reference_a1'}, inplace=True)

        merged_file = pd.merge(left=freq_file_with_position, right=legend_file, how='inner', on='position')

        #The MAF column in the frq file is the allele frequency of the A1 allele (which is usually minor, but
        # possibly not in my case because I just updated the reference to match 1000G.
        #The AF columns in the legend file are the allele frequencies of the a1 allele in that file.
        #Need to match based on A1 alleles (hopefully this has already been done in the harmonization step, then
        # calculate the allele frequency difference.

        merged_file['AFR_Match_Diff'] = np.where((merged_file['dataset_a1']==merged_file['reference_a1']) &
                                                 (merged_file['dataset_a2'] == merged_file['reference_a0']),
                                                 abs(merged_file['dataset_a1_frq'] - merged_file['AFR']),'')
        merged_file['AMR_Match_Diff'] = np.where((merged_file['dataset_a1']==merged_file['reference_a1']) &
                                                 (merged_file['dataset_a2'] == merged_file['reference_a0']),
                                                 abs(merged_file['dataset_a1_frq'] - merged_file['AMR']), '')
        merged_file['EAS_Match_Diff'] = np.where((merged_file['dataset_a1']==merged_file['reference_a1']) &
                                                  (merged_file['dataset_a2'] == merged_file['reference_a0']),
                                                 abs(merged_file['dataset_a1_frq'] - merged_file['EAS']), '')
        merged_file['EUR_Match_Diff'] = np.where((merged_file['dataset_a1']==merged_file['reference_a1']) &
                                                  (merged_file['dataset_a2'] == merged_file['reference_a0']),
                                                 abs(merged_file['dataset_a1_frq'] - merged_file['EUR']), '')
        merged_file['SAS_Match_Diff'] = np.where((merged_file['dataset_a1']==merged_file['reference_a1']) &
                                                  (merged_file['dataset_a2'] == merged_file['reference_a0']),
                                                 abs(merged_file['dataset_a1_frq'] - merged_file['SAS']), '')
        merged_file['AFR_FlipMatch_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a0']) &
                                                      (merged_file['dataset_a2'] == merged_file['reference_a1']),
                                                     abs(merged_file['dataset_a2_frq'] - merged_file['AFR']), '')
        merged_file['AMR_FlipMatch_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a0']) &
                                                      (merged_file['dataset_a2'] == merged_file['reference_a1']),
                                                     abs(merged_file['dataset_a2_frq'] - merged_file['AMR']), '')
        merged_file['EAS_FlipMatch_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a0']) &
                                                      (merged_file['dataset_a2'] == merged_file['reference_a1']),
                                                     abs(merged_file['dataset_a2_frq'] - merged_file['EAS']), '')
        merged_file['EUR_FlipMatch_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a0']) &
                                                      (merged_file['dataset_a2'] == merged_file['reference_a1']),
                                                     abs(merged_file['dataset_a2_frq'] - merged_file['EUR']), '')
        merged_file['SAS_FlipMatch_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a0']) &
                                                      (merged_file['dataset_a2'] == merged_file['reference_a1']),
                                                     abs(merged_file['dataset_a2_frq'] - merged_file['SAS']), '')
        merged_file['AFR_Diff'] = merged_file['AFR_Match_Diff']+merged_file['AFR_FlipMatch_Diff']
        merged_file['AMR_Diff'] = merged_file['AMR_Match_Diff']+merged_file['AMR_FlipMatch_Diff']
        merged_file['EAS_Diff'] = merged_file['EAS_Match_Diff']+merged_file['EAS_FlipMatch_Diff']
        merged_file['EUR_Diff'] = merged_file['EUR_Match_Diff']+merged_file['EUR_FlipMatch_Diff']
        merged_file['SAS_Diff'] = merged_file['SAS_Match_Diff']+merged_file['SAS_FlipMatch_Diff']

        merged_file.drop(labels = ['AFR_Match_Diff', 'AFR_FlipMatch_Diff', 'AMR_Match_Diff', 'AMR_FlipMatch_Diff',
                                   'EAS_Match_Diff', 'EAS_FlipMatch_Diff', 'EUR_Match_Diff', 'EUR_FlipMatch_Diff',
                                   'SAS_Match_Diff', 'SAS_FlipMatch_Diff'],axis=1, inplace = True)

        merged_file[['AFR_Diff', 'AMR_Diff', 'EAS_Diff', 'EUR_Diff', 'SAS_Diff']] = \
            merged_file[['AFR_Diff', 'AMR_Diff', 'EAS_Diff', 'EUR_Diff', 'SAS_Diff']].apply(pd.to_numeric, errors = 'coerce')

        merged_file['AF_Decision'] = np.where((merged_file['AFR_Diff'] > 0.2) & (merged_file['AMR_Diff'] > 0.2) &
                                              (merged_file['EAS_Diff'] > 0.2) & (merged_file['EUR_Diff'] > 0.2) &
                                              (merged_file['SAS_Diff'] > 0.2), 'Remove', 'Keep')

        merged_file.drop_duplicates(subset=['SNP'], keep = False, inplace=True)

        AF_diff_removed_by_chr[i] = merged_file[merged_file['AF_Decision'] == 'Remove']

        final_snp_lists_by_chr[i] = merged_file[merged_file['AF_Decision'] == 'Keep']
        final_snp_lists_by_chr[i]['SNP'].to_csv(final_snp_list_names[i], sep='\t', header=False, index=False)

        os.system('plink --bfile ' + harmonized_geno_names[i] + ' --extract ' + final_snp_list_names[i] +
                  ' --recode vcf-iid bgz --make-bed --out ' + harmonized_geno_names[i] + '_AFChecked')

        print('Finished with chr' + str(i + 1))

# Now for chromosome X
    os.system('java -Xmx1g -jar GenotypeHarmonizer.jar $* '
              '--input ' + geno_name +
              '--ref ALL.chrX.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
              '--refType VCF '
              '--chrFilter 23 '
              '--hweFilter 0.01 '
              '--mafFilter 0.05 '
              '--update-id '
              '--debug '
              '--mafAlign 0 '
              '--update-reference-allele '
              '--outputType PLINK_BED '
              '--output ' + geno_name + '_chrX_harmonized')

    # Since we are going to use a global reference population, removes all SNPs with an AF difference > 0.2
    #  between study dataset and all superpopulation allele frequencies. IF within 0.2 of any superpopulation frequency,
    #  keep variant. Frequency file is expected to be a Plink frequency file with the same number of variants as the bim file.

    os.system('plink --bfile ' + geno_name + '_chrX_harmonized --freq --out ' + geno_name + '_chrX_harmonized')

    freq_file = pd.read_csv(geno_name + '_chrX_harmonized.frq', sep='\s+', header=0, usecols=[0, 1, 2, 3, 4],
                            dtype={'CHR': int, 'SNP': str, 'A1': str, 'A2': str, 'MAF': float})
    freq_file.rename(columns={'A1': 'dataset_a1', 'A2': 'dataset_a2', 'MAF': 'dataset_a1_frq'}, inplace=True)
    freq_file['dataset_a2_frq'] = 1 - freq_file['dataset_a1_frq']

    bim_file = pd.read_csv(geno_name + '_chrX_harmonized.bim', sep='\s+', header=None, usecols=[0, 1, 3],
                           names=['CHR', 'SNP', 'position'])

    freq_file_with_position = pd.merge(left=freq_file, right=bim_file, how='inner', on=['CHR', 'SNP'])

    PAR1_legend_file = pd.read_csv('1000GP_Phase3_chrX_PAR1.legend.gz', compression="gzip", sep=" ", header=0,
                              dtype={'id': str, 'position': int, 'a0': str, 'a1': str, 'TYPE': str, 'AFR': float,
                                     'AMR': float, 'EAS': float, 'EUR': float, 'SAS': float, 'ALL': float})
    PAR1_legend_file.rename(columns={'id': 'reference_id', 'a0': 'reference_a0', 'a1': 'reference_a1'}, inplace=True)

    PAR2_legend_file = pd.read_csv('1000GP_Phase3_chrX_PAR2.legend.gz', compression="gzip", sep=" ", header=0,
                              dtype={'id': str, 'position': int, 'a0': str, 'a1': str, 'TYPE': str, 'AFR': float,
                                     'AMR': float, 'EAS': float, 'EUR': float, 'SAS': float, 'ALL': float})
    PAR2_legend_file.rename(columns={'id': 'reference_id', 'a0': 'reference_a0', 'a1': 'reference_a1'}, inplace=True)

    legend_file = pd.concat(['PAR1_legend_file', 'PAR2_legend_file'])

    merged_file = pd.merge(left=freq_file_with_position, right=legend_file, how='inner', on='position')

    # The MAF column in the frq file is the allele frequency of the A1 allele (which is usually minor, but
    # possibly not in my case because I just updated the reference to match 1000G.
    # The AF columns in the legend file are the allele frequencies of the a1 allele in that file.
    # Need to match based on A1 alleles (hopefully this has already been done in the harmonization step, then
    # calculate the allele frequency difference.

    merged_file['AFR_Match_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a1']) &
                                             (merged_file['dataset_a2'] == merged_file['reference_a0']),
                                             abs(merged_file['dataset_a1_frq'] - merged_file['AFR']), '')
    merged_file['AMR_Match_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a1']) &
                                             (merged_file['dataset_a2'] == merged_file['reference_a0']),
                                             abs(merged_file['dataset_a1_frq'] - merged_file['AMR']), '')
    merged_file['EAS_Match_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a1']) &
                                             (merged_file['dataset_a2'] == merged_file['reference_a0']),
                                             abs(merged_file['dataset_a1_frq'] - merged_file['EAS']), '')
    merged_file['EUR_Match_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a1']) &
                                             (merged_file['dataset_a2'] == merged_file['reference_a0']),
                                             abs(merged_file['dataset_a1_frq'] - merged_file['EUR']), '')
    merged_file['SAS_Match_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a1']) &
                                             (merged_file['dataset_a2'] == merged_file['reference_a0']),
                                             abs(merged_file['dataset_a1_frq'] - merged_file['SAS']), '')
    merged_file['AFR_FlipMatch_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a0']) &
                                                 (merged_file['dataset_a2'] == merged_file['reference_a1']),
                                                 abs(merged_file['dataset_a2_frq'] - merged_file['AFR']), '')
    merged_file['AMR_FlipMatch_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a0']) &
                                                 (merged_file['dataset_a2'] == merged_file['reference_a1']),
                                                 abs(merged_file['dataset_a2_frq'] - merged_file['AMR']), '')
    merged_file['EAS_FlipMatch_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a0']) &
                                                 (merged_file['dataset_a2'] == merged_file['reference_a1']),
                                                 abs(merged_file['dataset_a2_frq'] - merged_file['EAS']), '')
    merged_file['EUR_FlipMatch_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a0']) &
                                                 (merged_file['dataset_a2'] == merged_file['reference_a1']),
                                                 abs(merged_file['dataset_a2_frq'] - merged_file['EUR']), '')
    merged_file['SAS_FlipMatch_Diff'] = np.where((merged_file['dataset_a1'] == merged_file['reference_a0']) &
                                                 (merged_file['dataset_a2'] == merged_file['reference_a1']),
                                                 abs(merged_file['dataset_a2_frq'] - merged_file['SAS']), '')
    merged_file['AFR_Diff'] = merged_file['AFR_Match_Diff'] + merged_file['AFR_FlipMatch_Diff']
    merged_file['AMR_Diff'] = merged_file['AMR_Match_Diff'] + merged_file['AMR_FlipMatch_Diff']
    merged_file['EAS_Diff'] = merged_file['EAS_Match_Diff'] + merged_file['EAS_FlipMatch_Diff']
    merged_file['EUR_Diff'] = merged_file['EUR_Match_Diff'] + merged_file['EUR_FlipMatch_Diff']
    merged_file['SAS_Diff'] = merged_file['SAS_Match_Diff'] + merged_file['SAS_FlipMatch_Diff']

    merged_file.drop(labels=['AFR_Match_Diff', 'AFR_FlipMatch_Diff', 'AMR_Match_Diff', 'AMR_FlipMatch_Diff',
                             'EAS_Match_Diff', 'EAS_FlipMatch_Diff', 'EUR_Match_Diff', 'EUR_FlipMatch_Diff',
                             'SAS_Match_Diff', 'SAS_FlipMatch_Diff'], axis=1, inplace=True)

    merged_file = merged_file[~np.isnan(merged_file['AFR_Diff', 'AMR_Diff', 'EAS_Diff', 'EUR_Diff', 'SAS_Diff'])]

    merged_file[['AFR_Diff', 'AMR_Diff', 'EAS_Diff', 'EUR_Diff', 'SAS_Diff']] = \
        merged_file[['AFR_Diff', 'AMR_Diff', 'EAS_Diff', 'EUR_Diff', 'SAS_Diff']].apply(pd.to_numeric, errors = 'coerce')

    merged_file['AF_Decision'] = np.where((merged_file['AFR_Diff'] > 0.2) & (merged_file['AMR_Diff'] > 0.2) &
                                          (merged_file['EAS_Diff'] > 0.2) & (merged_file['EUR_Diff'] > 0.2) &
                                          (merged_file['SAS_Diff'] > 0.2), 'Remove', 'Keep')

    merged_file.drop_duplicates(subset=['SNP'], keep=False, inplace=True)

    ChrX_SNPs_Removed = merged_file[merged_file['AF_Decision'] == 'Remove']

    ChrX_SNPs_Kept = merged_file[merged_file['AF_Decision'] == 'Keep']
    ChrX_SNPs_Kept['SNP'].to_csv('chrX_SNPsKept_List.txt', sep='\t', header=False, index=False)

    os.system('plink --bfile ' + geno_name + '_chrX_harmonized --extract chrX_SNPsKept_List.txt --recode vcf-iid bgz '
                                             '--make-bed --out ' + geno_name + '_chrX_harmonized_AFChecked')

    All_SNPs_Removed = pd.concat(AF_diff_removed_by_chr[0], AF_diff_removed_by_chr[1], AF_diff_removed_by_chr[2],
                                 AF_diff_removed_by_chr[3], AF_diff_removed_by_chr[4], AF_diff_removed_by_chr[5],
                                 AF_diff_removed_by_chr[6], AF_diff_removed_by_chr[7], AF_diff_removed_by_chr[8],
                                 AF_diff_removed_by_chr[9], AF_diff_removed_by_chr[10], AF_diff_removed_by_chr[11],
                                 AF_diff_removed_by_chr[12], AF_diff_removed_by_chr[13], AF_diff_removed_by_chr[14],
                                 AF_diff_removed_by_chr[15], AF_diff_removed_by_chr[16], AF_diff_removed_by_chr[17],
                                 AF_diff_removed_by_chrp[18], AF_diff_removed_by_chr[19], AF_diff_removed_by_chr[20],
                                 AF_diff_removed_by_chr[21], ChrX_SNPs_Removed)
    All_SNPs_Removed.to_csv('SNPs_Removed_AFDiff.txt', sep='\t', header = True, index = False)

    All_SNPs_Kept = pd.concat(final_snp_lists_by_chr[0], final_snp_lists_by_chr[1], final_snp_lists_by_chr[2],
                              final_snp_lists_by_chr[3], final_snp_lists_by_chr[4], final_snp_lists_by_chr[5],
                              final_snp_lists_by_chr[6], final_snp_lists_by_chr[7], final_snp_lists_by_chr[8],
                              final_snp_lists_by_chr[9], final_snp_lists_by_chr[10], final_snp_lists_by_chr[11],
                              final_snp_lists_by_chr[12], final_snp_lists_by_chr[13], final_snp_lists_by_chr[14],
                              final_snp_lists_by_chr[15], final_snp_lists_by_chr[16], final_snp_lists_by_chr[17],
                              final_snp_lists_by_chr[18], final_snp_lists_by_chr[19], final_snp_lists_by_chr[20],
                              final_snp_lists_by_chr[21], ChrX_SNPs_Kept)

    All_SNPs_Kept.to_csv('SNPs_Kept.txt', sep='\t', header = True, index = False)

    os.system('plink --vcf ' + vcf_file_names[i] + '.vcf.gz --double-id --biallelic-only --vcf-require-gt --make-bed --out '
              + vcf_file_names[i])
    os.system('plink --vcf All.chrX.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --double-id '
              '--biallelic-only --vcf-require-gt --make-bed --out All.chrX.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes')

    merge_list = AF_checked_names + vcf_file_names
    merge_list.extend([geno_name + '__chrX_harmonized_AFChecked', 'All.chrX.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes'])

    write_merge_list = open("MergeList.txt", "w")
    print >> write_merge_list, "\n".join(str(i) for i in merge_list)
    write_merge_list.close()

    os.system('plink --bmerge ' + harmonized_geno_names[0] + '_AFChecked --merge-list MergeList.txt --make-bed --out ' + geno_name + '_1000G' )

elif to_do == '6':
    print("You go, couch potato")

else:
    print("Please enter a number 1-6.")



