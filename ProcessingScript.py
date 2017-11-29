# Julie's Processing script for Shriver Lab genotype data
# Date 9.5.17

# WHAT THIS DOES #
# 1) Run IBD matrix through plink
# 2) Update family and individual IDs
# 3) Update maternal and paternal IDs
# 4) Update sex
# 5) Using 1000 Genomes as a reference (this part based off Perl script by W. Rayner, 2015, wrayner@well.ox.ac.uk)
    #   -Removes SNPs with MAF < 5% in study dataset
    #   -Removes SNPs not in 1000 Genomes
    #   -Removes all A/T G/C SNPs with MAF > 40% in the reference data set
    #   -Removes all SNPs with an AF difference >0.2, between reference and dataset frequency file, frequency file is
    #       expected to be a plink frequency file with the same number of SNPs as the bim file
    #   -Removes duplicates that may be introduced with the position update
    #   -Removes indels
# 6) Merges your data with 1000 Genomes
# 7) Prepares files for ADMIXTURE with k = 3...9
# 8) Prepares files for phasing using SHAPEIT2
# 9) Prepares files for imputation using the Sanger Imputation Server

# REQUIREMENTS #
# You must have plink 1.9 (https://www.cog-genomics.org/plink2) in the same directory as your genotype files and this script,
#   and it should be called 'plink' Alternatively, you could put plink in your PATH.
# You must have downloaded GenotypeHarmonizer (https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer-Download)
#   and either have it in the same directory as your genotype files and this script. Or, put it in your PATH
# You must have downloaded the 1000G Phase 3 legend files from http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/.
#   These can be anywhere you want, you'll tell me the path later
# You must have downloaded the 1000G Phase 3 VCF files from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/.
#   These can be anywhere you want, you'll tell me the path later

# Getting the needed modules.
import os

to_do = input('What would you like to do?\n'
              '1) Run an IBD analysis to identify relatives. All you need are plink bed/bim/fam files.\n'
              '2) Update FID or IID information. You need a file with the following information Old FID, Old IID, '
              'New FID, New IID.\n'
              '3) Update parental IDs. You need a file with FID, IID, Paternal IID, and Maternal IID.\n'
              '4) Update sex. You need a file with FID, IID, Sex (M = 1, F = 2, Unknown = 0)\n'
              '5) Harmonize with 1000 Genomes Phase 3 (need to do this before merging or phasing)\n'
              '6) Merge your data with 1000G\n'
              '7) Prepares files for ADMIXTURE with k = 3...9\n'
              '8) Nothing\n'
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
    # Outputs new files per chromosome, in vcf and plink bed/bim/fam format.

    import pandas as pd
    import numpy as np

    geno_name = input('Please enter the name of the genotype files (without bed/bim/fam extension: ')
    vcf_path = input('Please enter the pathname of where your 1000G vcf files are (i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\ etc.):')
    legend_path = input('Please enter the pathname of where your 1000G legend files are (i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\ etc.):')

    ref_file_names = ['ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' % x for x in range(1,23)]
    ref_file_names.extend(['ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'])
    harmonized_geno_names = [geno_name + '_chr%d_harmonized' % x for x in range(1, 24)]
    freq_file_names = [geno_name + '_chr%d_harmonized.frq' % x for x in range(1,24)]
    legend_file_names = ['1000GP_Phase3_chr%d.legend.gz' % x for x in range(1,23)]
    legend_file_names.extend(['1000GP_Phase3_chrX_NONPAR.legend'])
    AF_diff_removed_by_chr = ['chr%d_SNPsRemoved_AFDiff' % x for x in range(1,24)]
    final_snp_lists_by_chr = ['chr%d_SNPsKept' % x for x in range(1,24)]
    AF_checked_names = [geno_name + '_chr%d_HarmonizedTo1000G' % x for x in range(1,24)]

    # Harmonize for each chromosome
    for i in range(0, len(ref_file_names)):
        if i < 22:
            os.system('java -Xmx1g -jar GenotypeHarmonizer.jar $* '
                      '--input ' + geno_name +
                      ' --ref "' + os.path.join(vcf_path,ref_file_names[i]) +
                      '" --refType VCF '
                      '--chrFilter ' + str(i + 1) +
                      ' --hweFilter 0.01 '
                      '--mafFilter 0.05 '
                      '--update-id '
                      '--debug '
                      '--mafAlign 0 '
                      '--update-reference-allele '
                      '--outputType PLINK_BED '
                      '--output ' + harmonized_geno_names[i])

        elif i == 22:
            # Since chr X is labeled as 23 in the plink files and X in the vcf files, we need to separate it out and convert the 23 to X before harmonizing
            os.system('plink --bfile ' + geno_name + ' --chr X --make-bed --out ' + geno_name + '_chr23')
            bim_file = pd.read_csv(geno_name + '_chr23.bim', sep='\t', header=None)
            bim_file.iloc[:, 0].replace(23, 'X', inplace=True)
            bim_file.to_csv(geno_name + '_chr23.bim', sep='\t', header=False, index=False, na_rep='NA')

            os.system('java -Xmx1g -jar GenotypeHarmonizer.jar $* '
                      '--input ' + geno_name + '_chr23'
                      ' --ref "' + os.path.join(vcf_path,ref_file_names[i]) +
                      '" --refType VCF '
                      '--hweFilter 0.01 '
                      '--mafFilter 0.05 '
                      '--update-id '
                      '--debug '
                      '--mafAlign 0 '
                      '--update-reference-allele '
                      '--outputType PLINK_BED '
                      '--output ' + harmonize_geno_names[i])
        else:
            print("Something is wrong with the number/name of reference files")

    # Now for each that was just harmonized, remove all SNPs with an allele (AF) difference > 0.2 since we are going to use a global reference population
    # between study dataset and all superpopulation allele frequencies. IF within 0.2 of any superpopulation frequency,
    # keep variant. Frequency file is expected to be a Plink frequency file with the same number of variants as the bim file.
    for i in range(0, len(harmonized_geno_names)):
        os.system('plink --bfile ' + harmonized_geno_names[i] + ' --freq --out ' + harmonized_geno_names[i])
        freq_file = pd.read_csv(freq_file_names[i], sep='\s+', header = 0, usecols = [0,1,2,3,4],
                                dtype = {'CHR': int, 'SNP': str, 'A1': str, 'A2': str, 'MAF':float})
        freq_file.rename(columns={'A1': 'dataset_a1', 'A2': 'dataset_a2', 'MAF': 'dataset_a1_frq'}, inplace=True)
        freq_file['dataset_a2_frq'] = 1-freq_file['dataset_a1_frq']

        bim_file = pd.read_csv(harmonized_geno_names[i] + '.bim', sep='\s+', header=None, usecols=[0, 1, 3],
                               names=['CHR', 'SNP', 'position'])

        freq_file_with_position = pd.merge(left = freq_file, right = bim_file, how = 'inner', on=['CHR','SNP'])


        legend_file = pd.read_csv(os.path.join(legend_path,legend_file_names[i]), compression="gzip", sep=" ", header = 0,
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
        final_snp_lists_by_chr[i]['SNP'].to_csv(final_snp_lists_by_chr[i] + '.txt', sep='\t', header=False, index=False)

        # Make both bed and vcf files for each chromosomes. Need bed file for merging, and need vcf file for phasing
        os.system('plink --bfile ' + harmonized_geno_names[i] + ' --extract ' + final_snp_lists_by_chr[i] +
                  '.txt --recode vcf-iid bgz --make-bed --out ' + AF_checked_names[i])

        if os.path.getsize(AF_checked_names[i] + '.bim') > 0:
            os.system('rm ' + final_snp_lists_by_chr[i] + '.txt')
            os.system('rm ' + harmonized_geno_names[i] + '.bed')
            os.system('rm ' + harmonized_geno_names[i] + '.bim')
            os.system('rm ' + harmonized_geno_names[i] + '.fam')
            os.system('rm ' + harmonized_geno_names[i] + '.log')
            os.system('rm ' + harmonized_geno_names[i] + '.frq')
            os.system('rm ' + harmonized_geno_names[i] + '.nosex')
            os.system('rm ' + harmonized_geno_names[i] + '.hh')

        print('Finished with chr' + str(i + 1))

    # Write a list of all SNPs removed and all SNPs kept just for reference purposes.
    All_SNPs_Removed = pd.concat([AF_diff_removed_by_chr[0], AF_diff_removed_by_chr[1], AF_diff_removed_by_chr[2],
                                  AF_diff_removed_by_chr[3], AF_diff_removed_by_chr[4], AF_diff_removed_by_chr[5],
                                  AF_diff_removed_by_chr[6], AF_diff_removed_by_chr[7], AF_diff_removed_by_chr[8],
                                  AF_diff_removed_by_chr[9], AF_diff_removed_by_chr[10], AF_diff_removed_by_chr[11],
                                  AF_diff_removed_by_chr[12], AF_diff_removed_by_chr[13], AF_diff_removed_by_chr[14],
                                  AF_diff_removed_by_chr[15], AF_diff_removed_by_chr[16], AF_diff_removed_by_chr[17],
                                  AF_diff_removed_by_chr[18], AF_diff_removed_by_chr[19], AF_diff_removed_by_chr[20],
                                  AF_diff_removed_by_chr[21], AF_diff_removed_by_chr[22]])
    All_SNPs_Removed.to_csv('SNPs_Removed_AFDiff.txt', sep='\t', header = True, index = False)

    All_SNPs_Kept = pd.concat([final_snp_lists_by_chr[0], final_snp_lists_by_chr[1], final_snp_lists_by_chr[2],
                               final_snp_lists_by_chr[3], final_snp_lists_by_chr[4], final_snp_lists_by_chr[5],
                               final_snp_lists_by_chr[6], final_snp_lists_by_chr[7], final_snp_lists_by_chr[8],
                               final_snp_lists_by_chr[9], final_snp_lists_by_chr[10], final_snp_lists_by_chr[11],
                               final_snp_lists_by_chr[12], final_snp_lists_by_chr[13], final_snp_lists_by_chr[14],
                               final_snp_lists_by_chr[15], final_snp_lists_by_chr[16], final_snp_lists_by_chr[17],
                               final_snp_lists_by_chr[18], final_snp_lists_by_chr[19], final_snp_lists_by_chr[20],
                               final_snp_lists_by_chr[21], final_snp_lists_by_chr[22]])
    All_SNPs_Kept.to_csv('SNPs_Kept.txt', sep='\t', header = True, index = False)

elif to_do == '6':
    #Merges the 1000G with house dataset

    import csv
    import pandas as pd
    import numpy as np

    geno_name = input('Please enter the name of the genotype files (without bed/bim/fam extension: ')
    vcf_path = input('Please enter the pathname of where your 1000G vcf files are (i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\ etc.):')

    ref_file_names = ['ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' % x for x in
                      range(1, 23)]
    ref_file_names.extend(['ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'])
    AF_checked_names = [geno_name + '_chr%d_HarmonizedTo1000G' % x for x in range(1, 24)]
    chr_1000G_Phase3_names = ['chr%d_1000G_Phase3' % x for x in range(1, 24)]
    unique_rsid_per_chr = ['chr%d_unique_rsid' % x for x in range(1, 24)]

    with open("HouseMergeList.txt", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(AF_checked_names)
    os.system('plink --merge-list HouseMergeList.txt --geno 0.01 --make-bed --out ' + geno_name + '_HarmonizedTo1000G')

    if os.path.getsize(geno_name + '_HarmonizedTo1000G.bim') > 0:
        for i in range(0,len(AF_checked_names)):
            os.system('rm ' + AF_checked_names[i] + '.bed')
            os.system('rm ' + AF_checked_names[i] + '.bim')
            os.system('rm ' + AF_checked_names[i] + '.fam')
            os.system('rm ' + AF_checked_names[i] + '.log')
            os.system('rm ' + AF_checked_names[i] + '.nosex')
            os.system('rm ' + AF_checked_names[i] + '.hh')
    
    #Merge 1000G chr data into one plink formatted file, need to convert from vcf files - but only taking the snps that are in the house dataset
    house_snps_kept = pd.read_csv('SNPs_Kept.txt', header=0, sep='\t')
    house_snps_kept = house_snps_kept.loc[:, ['SNP']]
    house_snps_kept.to_csv('SNPs_Kept_List.txt', sep='\t', header=False, index=False)

    for i in range(0, len(ref_file_names)):
        os.system('plink --vcf "' + os.path.join(vcf_path,ref_file_names[i]) + '" --double-id --biallelic-only strict '
                                                                              '--vcf-require-gt --extract SNPs_Kept_List.txt '
                                                                              '--make-bed '
                                                                              '--out ' + chr_1000G_Phase3_names[i])
    os.system('rm *~')

    with open("1000GMergeList.txt", "w") as f:
        wr = csv.writer(f, delimiter = "\n")
        wr.writerow(chr_1000G_Phase3_names)
    os.system('plink --merge-list 1000GMergeList.txt --geno 0.01 --make-bed --out 1000G_Phase3')

    logfile=pd.DataFrame()
    with open('1000G_Phase3.log', 'r') as f:
        for line in f:
            logfile = pd.concat([logfile, pd.DataFrame([tuple(line.strip().split(" "))])], ignore_index=True)

    if logfile[0].str.contains('Warning:').any():
        rsid_warnings = logfile.loc[logfile[0] == 'Warning:', 6].str.split("'", expand=True)
        rsid_warnings[1].dropna(how='any').to_csv('1000G_MergeWarnings.txt', sep='\t', header=False, index=False)

    if os.path.exists('1000G_MergeWarnings.txt'):
        if os.path.exists('1000G_Phase3-merge.missnp'):
            # If merge warnings and missnps exist, exclude both from 1000G completely (there are plenty of other snps)
            missnp = pd.read_csv('1000G_Phase3-merge.missnp', sep='\t', header=None)
            warnings_missnp = pd.concat([rsid_warnings, missnp], axis=0)
            warnings_missnp[1].dropna(how='any').to_csv('1000G_warnings_missnp.txt', sep='\t',
                                                        header=False, index=False)
            for i in range(0,len(chr_1000G_Phase3_names)):
                os.system('plink --bfile ' + chr_1000G_Phase3_names[i] + ' --exclude 1000G_warnings_missnp.txt --geno 0.01 --make-bed --out '+ chr_1000G_Phase3_names[i])
            os.system('rm *~')
            os.system('plink --merge-list 1000GMergeList.txt --geno 0.01 --make-bed --out 1000G_Phase3')
        else:  # If only merge warnings exist, exclude from 1000G completely
            for i in range(0,len(chr_1000G_Phase3_names)):
                os.system('plink --bfile ' + chr_1000G_Phase3_names[i] + ' --exclude 1000G_MergeWarnings.txt --geno 0.01 --make-bed --out '+ chr_1000G_Phase3_names[i])
            os.system('rm *~')
            os.system('plink --merge-list 1000GMergeList.txt --geno 0.01 --make-bed --out 1000G_Phase3')
    elif os.path.exists('1000G_Phase3-merge.missnp') and os.path.getsize('1000G_MergeWarnings.txt') == 0:
        # If only the missnps still exist, remove them in 1000G.
        for i in range(0, len(chr_1000G_Phase3_names)):
            os.system('plink --bfile ' + chr_1000G_Phase3_names[i] + ' --exclude 1000G_Phase3-merge.missnp --geno 0.01 --make-bed --out ' + chr_1000G_Phase3_names[i])
        os.system('rm *~')
        os.system('plink --merge-list 1000GMergeList.txt --geno 0.01 --make-bed --out 1000G_Phase3')
    elif os.path.exists('1000G_Phase3.bim'):
        print("Successfully merged 1000G")
        for i in range(0,len(chr_1000G_Phase3_names)):
            os.system('rm ' + chr_1000G_Phase3_names[i] + '.*')
    else:
        print("Unable to merge 1000G chromosome files.")

    #Perform initial merge of house data and 1000G
    os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --bmerge 1000G_Phase3 --geno 0.01 --make-bed --out ' + geno_name + '_1000G')

    #Read in log file to see if anything went wrong
    logfile = pd.DataFrame()
    with open(geno_name + '_1000G.log', 'r') as f:
        for line in f:
            logfile = pd.concat([logfile, pd.DataFrame([tuple(line.strip().split(" "))])], ignore_index=True)

    #Check to see if the logfile has warnings or if there are triallelic snps that need to be flipped (missnp)
    if logfile[0].str.contains('Warning:').any():
        rsid_warnings = logfile.loc[logfile[0] == 'Warning:', 6].str.split("'", expand=True)
        rsid_warnings[1].dropna(how='any').to_csv(geno_name + '_1000G_MergeWarnings.txt', sep='\t', header=False,
                                                  index=False)

    if os.path.exists(geno_name + '_1000G_MergeWarnings.txt'):
        if os.path.exists(geno_name + '_1000G-merge.missnp'):
            #If merge warnings and missnps exist, exclude warning snps from both 1000G and house dataset, flip snps in house dataset
            os.system('plink --bfile 1000G_Phase3' + ' --exclude ' + geno_name + '_1000G_MergeWarnings.txt --geno 0.01 '
                                                                                 '--make-bed --out 1000G_Phase3')
            os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --exclude ' + geno_name +
                      '_1000G_MergeWarnings.txt --flip ' + geno_name + '_1000G-merge.missnp --geno 0.01 --make-bed --out '
                      + geno_name + '_HarmonizedTo1000G')
            os.system('rm *~')
            os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --bmerge 1000G_Phase3 --geno 0.01 --make-bed --out ' + geno_name + '_1000G_merge2')
        else: #If only mergewarnings exists, exclude warning snps from both 100)G and house dataset.
            os.system('plink --bfile 1000G_Phase3' + ' --exclude ' + geno_name + '_1000G_MergeWarnings.txt --geno 0.01 '
                                                                                 '--make-bed --out 1000G_Phase3')
            os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --exclude ' + geno_name +
                      '_1000G_MergeWarnings.txt --geno 0.01 --make-bed --out ' + geno_name + '_HarmonizedTo1000G')
            os.system('rm *~')
            os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --bmerge 1000G_Phase3 --geno 0.01 --make-bed --out ' + geno_name + '_1000G_merge2')
    elif os.path.exists(geno_name + '_1000G-merge.missnp') and not os.path.exists(geno_name + '_1000G_MergeWarnings.txt'):
        #If only the missnps exist, flip them in the house dataset.
        os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --flip ' + geno_name + '_1000G-merge.missnp --geno 0.01 --make-bed --out '
                  + geno_name + '_HarmonizedTo1000G')
        os.system('rm *~')
        os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --bmerge 1000G_Phase3 --geno 0.01 --make-bed --out ' + geno_name + '_1000G_merge2')
    elif os.path.exists(geno_name + '_1000G.bim'):
        print("Successfully merged house dataset with 1000G, though you should double check. I can't predict every error.")
    else:
        print("The house dataset did not merge properly with 1000G, but I'm not sure why. I'm sorry.")

    if os.path.exists(geno_name + '_1000G_merge2.log'):
        logfile = pd.DataFrame()
        with open(geno_name + '_1000G_merge2.log', 'r') as f:
            for line in f:
                logfile = pd.concat([logfile, pd.DataFrame([tuple(line.strip().split(" "))])], ignore_index=True)

        # Check to see if the logfile has warnings or if there are triallelic snps that need to be flipped (missnp)
        if logfile[0].str.contains('Warning:').any():
            rsid_warnings = logfile.loc[logfile[0] == 'Warning:', 6].str.split("'", expand=True)
            rsid_warnings[1].dropna(how='any').to_csv(geno_name + '_1000G_merge2_warnings.txt', sep='\t', header=False,
                                                      index=False)

        if os.path.exists(geno_name + '_1000G_merge2_warnings.txt'):
            if os.path.exists(geno_name + '_1000G_merge2-merge.missnp'):
                #If merge warnings and missnps still exist, exclude from both 1000G and house dataset.
                missnp = pd.read_csv(geno_name + '_1000G_merge2-merge.missnp', sep = '\t', header=None)
                warnings_missnp = pd.concat([rsid_warnings, missnp], axis=0)
                warnings_missnp[1].dropna(how='any').to_csv(geno_name + '_1000G_warnings_missnp.txt', sep='\t', header=False, index=False)
                os.system('plink --bfile 1000G_Phase3 --exclude ' + geno_name + '_1000G_merge2_warnings_missnp.txt --geno 0.01 --make-bed --out 1000G_Phase3')
                os.system('plink --bfile ' +geno_name + '_HarmonizedTo1000G --exclude ' + geno_name + '_1000G_merge2_warnings_missnp.txt --geno 0.01 --make-bed --out 1000G_Phase3')
                os.system('rm *~')
                os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --bmerge 1000G_Phase3 --geno 0.01 '
                                                         '--make-bed --out ' + geno_name + '_1000G_merge3')
            else:  # If only merge warnings still exist, exclude from both 1000G and house dataset.
                os.system('plink --bfile 1000G_Phase3' + ' --exclude ' + geno_name + '_1000G_merge2_warnings.txt '
                                                                                     '--geno 0.01 --make-bed '
                                                                                     '--out 1000G_Phase3')
                os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --exclude ' + geno_name +
                          '_1000G_merge2_warnings.txt --geno 0.01 --make-bed --out ' + geno_name + '_HarmonizedTo1000G')
                os.system('rm *~')
                os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --bmerge 1000G_Phase3 --geno 0.01 '
                                                         '--make-bed --out ' + geno_name + '_1000G_merge3')
        elif os.path.exists(geno_name + '_1000G_merge2-merge.missnp') and not os.path.exists(geno_name + '_1000G_merge2_warnings.txt'):
            # If only the missnps still exist, remove them in the both datasets.
            os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --exclude ' + geno_name +
                      '_1000G_merge2-merge.missnp --geno 0.01 --make-bed --out ' + geno_name + '_HarmonizedTo1000G')
            os.system('plink --bfile 1000G_Phase3 --exclude ' + geno_name + '_1000G_merge2-merge.missnp --geno 0.01 --make-bed --out 1000G_Phase3')
            os.system('rm *~')
            os.system('plink --bfile ' + geno_name + '_HarmonizedTo1000G --bmerge 1000G_Phase3 --geno 0.01 --make-bed '
                                                     '--out ' + geno_name + '_1000G_merge3')
        elif os.path.exists(geno_name + '_1000G_merge2.bim'):
            print("Successfully merged house dataset with 1000G on 2nd try, though you should double check. I can't predict every error.")
        else:
            print("House dataset and 1000G did not successfully merge 2nd time, but I'm not sure why. I'm sorry!")

    if os.path.exists(geno_name + '_1000G_merge3.log'):
        logfile = pd.DataFrame()
        with open(geno_name + '_1000G_merge3.log', 'r') as f:
            for line in f:
                logfile = pd.concat([logfile, pd.DataFrame([tuple(line.strip().split(" "))])], ignore_index=True)

        # Check to see if the logfile still has warnings, if so, the user will need to identify them and take care of them manually
        if not os.path.exists(geno_name + '_1000G_merge3.bim'):
            if logfile[0].str.contains('Warning:').any():
                print("I'm sorry, the logfile still has warnings, even after removing snps that threw errors in the "
                      "first two tries, you'll have to manually deal with these using the merge3 log file.")
            if os.path.exists(geno_name + '_1000G_merge3-merge.missnp'):
                print("I'm sorry, there are still snps with 3+ variants present, even after flipping some and removing "
                      "the ones that the flip didn't solve you'll have to deal with these manually using the merge3 log file")
        if os.path.exists(geno_name + '_1000G_merge3.bim'):
            print("House dataset and 1000G merged correctly on 3rd try. Though you should double-check, I can't predict every error.")

elif to_do == '7':
    #Prepares files for an admixture run k = 3...9
    admixture_proceed_check = input("Some cautions/notes before you perform this step:\n"
                                    "1) If you'd like to compare your population to 1000G, then you should perform step "
                                    "6 before this step. If not, then go right ahead\n"
                                    "2) This will perform admixture runs from k = 3 - 9. If you'd like other admixture "
                                    "runs performed, then you should change this code to reflect that. Or, submit a "
                                    "request to change the code and I'll get around to it.\n"
                                    "3) You should have an ACI-B cluster allocation at Penn State to perform this step.\n"
                                    "4) This will write the files that you need, but you are responsible for the memory, node, and "
                                    "time usage (walltime = 150 hrs, nodes 1, ppn = 8, pmem = 8gb) and for putting them "
                                    "on the cluster and submitting them to admixture \n"
                                    "5) You will need the admixture program either on your path or in the same folder where "
                                    "you will submit this job.\n"
                                    "6) You will need to transfer the pbs files and genotype bed/bim/fam files to your cluster before running.\n"
                                    "7) Are you okay with all of this? (y/n): ")
    if admixture_proceed_check in ('y', 'Y', 'Yes', 'yes'):
        geno_name = input('Please enter the name of the genotype files that you would like to perform admixture on (without bed/bim/fam extension: ')
        allocation_name = input('Please enter the name of your cluster allocation: ')
        os.system('plink --bfile ' + geno_name + ' --indep-pairwise 50 10 0.1 --out ' + geno_name)
        os.system('plink --bfile ' + geno_name + ' --extract ' + geno_name + '.prune.in --make-bed --out ' + geno_name + '_LDpruned')
        with open (geno_name + '_Admixture_k3to6.pbs', 'a') as file:
            file.write('#PBS -l walltime=150:00:00\n'
                       '#PBS -l nodes=1:ppn=8\n'
                       '#PBS -l pmem=8gb\n'
                       '#PBS -A ' + allocation_name + '\n'
                       '#PBS -j oe\n'
                       'cd $PBS_O_WORKDIR\n'
                       'for K in {3..6}; do ./admixture --cv ' + geno_name + '.bed $K | tee ' + geno_name + '.log${K}.out; done')

        with open(geno_name + '_Admixture_k7to9.pbs', 'a') as file:
            file.write('#PBS -l walltime=150:00:00\n'
                       '#PBS -l nodes=1:ppn=8\n'
                       '#PBS -l pmem=8gb\n'
                       '#PBS -A ' + allocation_name + '\n'
                       '#PBS -j oe\n'
                       'cd $PBS_O_WORKDIR\n'
                       'for K in {7..9}; do ./admixture --cv ' + geno_name + '.bed $K | tee ' + geno_name + '.log${K}.out; done')
        print("Transfer " + geno_name + "_LDpruned bed/bim/fam files, " + geno_name + "_Admixture_k3to6.pbs, and "
              + geno_name + "_Admixture_7to9.pbs files to the cluster.\n"
              "Submit them using qsub " + geno_name + "Admixture_k3to6.pbs and qsub " + geno_name + "Admixture_k7to9.pbs")
        print("When you get your results, you should evaluate them to see which makes sense given your study population and which has the lowest CV value.")

    elif admixture_proceed_check in ('n', 'N', 'No', 'no'):
        print("Okay, we will not perform admixture at this time.")

    else:
        print("Please answer 'yes' or 'no.'")

elif to_do == '8':
    print("You go, couch potato")

else:
    print("Please enter a number 1-8.")



