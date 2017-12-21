def harmonize_with_1000g(geno_name):
    # Using 1000 Genomes as a reference(based off Perl script by W.Rayner, 2015, wrayner @ well.ox.ac.uk)
    #   -Removes SNPs with MAF < 5% in study dataset
    #   -Removes SNPs not in 1000 Genomes Phase 3
    #   -Removes all A/T G/C SNPs with MAF > 40% in the reference data set
    #   -Removes all SNPs with an AF difference >0.2, between reference and dataset frequency file, frequency file is
    #         expected to be a plink frequency file with the same number of SNPs as the bim file
    #   -Removes duplicates that may be introduced with the position update
    #   -Removes indels #Need to figure out how to do this. Or even if it is necessary with plink files.
    #   -Removes SNPs with HWE p-value < 0.01
    #   -Updates the reference allele to match 1000G
    #   -Outputs new files per chromosome, in plink bed/bim/fam format.

    # Needed modules
    import os
    import subprocess
    import shutil
    import glob
    import sys
    import pandas as pd
    import numpy as np
    import csv
    import gzip

    # Get current working directory.
    orig_wd = os.getcwd()
    
    # Getting required reference files and genotype harmonizer program
    # Ask the user if they already have the 1000G Phase 3 vcf files.
    vcf_exists = input('\u001b[34;1m Have you already downloaded the 1000G Phase3 VCF files? (y/n): \u001b[0m').lower()

    if vcf_exists in ('y', 'yes'):
        # Getting user's path to VCF files
        vcf_path = input('\u001b[35;1m Please enter the pathname of where your 1000G VCF files are '
                         '(i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\VCF\\ etc.): \u001b[0m')

    elif vcf_exists in ('n', 'no'):
        # Get module where downloading instructions are.
        import genodownload
        # From that module, call download 1000G Phase 3 VCF
        genodownload.vcf_1000g_phase3()
        # Saving VCF path
        vcf_path = os.path.join(os.getcwd(), '1000G_Phase3_VCF')
    else:
        sys.exit("Please answer yes or no. Quitting now because no VCF files.")

    # Ask the user if they already have the 1000G Phase 3 Hap/Legend/Sample files.
    legend_exists = input('\u001b[36;1m Have you already downloaded the 1000G Phase 3 Hap/Legend/Sample files? '
                          '(y/n): \u001b[0m').lower()

    if legend_exists in ('y', 'yes'):
        # Ask the user where the Legend files are.
        legend_path = input('\u001b[31;1m Please enter the pathname of where your 1000G legend files are '
                            '(i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\ etc.): \u001b[0m')

    elif legend_exists in ('n', 'no'):
        # Get genodownload module
        import genodownload
        # Call download HapLegendSample command
        genodownload.hls_1000g_phase3()
        # Saving legend path
        legend_path = os.path.join(os.getcwd(), '1000G_Phase3_HapLegendSample')
    else:
        sys.exit('Please answer yes or no. Quitting now because no legend files.')
    
    # Ask if they have the hg19 fasta files.
    fasta_exists = input('\u001b[32;1m Have you already downloaded the 1000G hg19 fasta file? (y/n): \u001b[0m').lower()

    if fasta_exists in ('y', 'yes'):
        # Ask the user where the fasta file is.
        fasta_path = input('\u001b[33;1m Please enter the pathname of where the your 1000G hg19 fasta file is '
                           '(i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\Fasta\\ etc.): '
                           '\u001b[0m')
    elif fasta_exists in ('n', 'no'):
        # Get geno download module
        import genodownload
        # Call download fasta command
        genodownload.fasta_1000G_hg19()
        # Saving fasta path
        fasta_path = os.path.join(os.getcwd(), '1000G_hg19_fasta')
    else:
        sys.exit('Please answer yes or no. Quitting now because no fasta file.')
    
    # Ask if they have genotype harmonizer.
    harmonizer_exists = input('\u001b[34;1m Have you already downloaded Genotype Harmonizer? (y/n): \u001b[0m').lower()

    if harmonizer_exists in ('y', 'yes'):
        # Ask the user where genotype harmonizer is.
        harmonizer_path = input('\u001b[35;1m Please enter the pathname of where the Genotype Harmonizer folder is '
                                '(i.e. C:\\Users\\Julie White\\Box Sync\\Software\\GenotypeHarmonizer-1.4.20\\): '
                                '\u001b[0m')

    elif harmonizer_exists in ('n', 'no'):
        import genodownload
        genodownload.genotype_harmonizer()
        # Harmonize path now that we've downloaded it.
        harmonizer_path = os.path.join(os.getcwd(), 'GenotypeHarmonizer-1.4.20/GenotypeHarmonizer-1.4.20-SNAPSHOT/')

    else:
        sys.exit('\u001b[35;1m Please write yes or no. Quitting now because no Genotype Harmonizer. \u001b[0m')

    # Start harmonization
    # Make new folder where the harmonized files will be located.
    if not os.path.exists('Harmonized_To_1000G'):
        os.makedirs('Harmonized_To_1000G')

    # Copy genotype files to new folder.
    shutil.copy2(geno_name + '.bed', 'Harmonized_To_1000G')
    shutil.copy2(geno_name + '.bim', 'Harmonized_To_1000G')
    shutil.copy2(geno_name + '.fam', 'Harmonized_To_1000G')

    # Copy plink to new folder.
    plink_files = glob.glob(r'plink')
    plink_files.extend(glob.glob(r'plink.exe'))

    if not plink_files:
        sys.exit("We need plink to run this part of the script. Please run step 1 if you do not have plink, or make "
                 "sure that it is in this directory.")
    else:
        for file in plink_files:
            print(file)
            shutil.copy(file, 'Harmonized_To_1000G')
    
    # Switch to this directory.
    os.chdir('Harmonized_To_1000G')
    
    # Make the lists that we're going to need, since this is on a per chromosome basis.
    vcf_file_names = ['ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' % x for x in
                      range(1, 23)]
    vcf_file_names.extend(['ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'])
    legend_file_names = ['1000GP_Phase3_chr%d.legend.gz' % x for x in range(1, 23)]
    legend_file_names.extend(['1000GP_Phase3_chrX_NONPAR.legend.gz'])
    harmonized_geno_names = [geno_name + '_chr%d_Harmonized' % x for x in range(1, 24)]
    id_updates = [s + '_idUpdates.txt' for s in harmonized_geno_names]
    id_update_names = [s + '_idUpdates.txt' for s in harmonized_geno_names]
    snp_logs = [s + '_snpLog.log' for s in harmonized_geno_names]
    snp_log_names = [s + '_snpLog.log' for s in harmonized_geno_names]
    freq_file_names = [geno_name + '_chr%d_Harmonized.frq' % x for x in range(1, 24)]
    af_diff_removed_by_chr = ['chr%d_SNPsRemoved_AFDiff' % x for x in range(1, 24)]
    final_snps_by_chr = ['chr%d_SNPsKept' % x for x in range(1, 24)]
    final_snp_lists = ['chr%d_SNPsKept.txt' % x for x in range(1, 24)]
    af_checked_names = [geno_name + '_chr%d_HarmonizedTo1000G' % x for x in range(1, 24)]

    
    # Harmonize for each chromosome
    for i in range(0, len(vcf_file_names)):
        # Call genotype harmonizer for autosomes
        if i < 22:
            subprocess.check_output('java -Xmx1g -jar "' + harmonizer_path + '/GenotypeHarmonizer.jar" $* --input '
                                    + geno_name + ' --ref "' + os.path.join(vcf_path, vcf_file_names[i])
                                    + '" --refType VCF --chrFilter ' + str(i+1)
                                    + ' --hweFilter 0.01 --mafFilter 0.05 --update-id --debug --mafAlign 0 '
                                      '--update-reference-allele --outputType PLINK_BED --output '
                                    + harmonized_geno_names[i], shell=True)

            id_updates[i] = pd.read_csv(id_update_names[i], sep='\t', header=0,
                                        dtype={'chr':str, 'pos':int, 'originalId':str, 'newId':str})
            snp_logs[i] = pd.read_table(snp_log_names[i], sep='\t', header=0,
                                      dtype={'chr': str, 'pos': int, 'id': str, 'alleles': str, 'action':str,
                                             'message':str})

        elif i == 22:
            # Since chr X is labeled as 23 in the plink files and X in the vcf files, we need to separate it out and
            # convert the 23 to X before harmonizing
            subprocess.check_output(['plink', '--bfile', geno_name, '--chr', 'X', '--make-bed', '--out',
                                     geno_name + '_chr23'])
            # Read chrX file into pandas
            bim_file = pd.read_csv(geno_name + '_chr23.bim', sep='\t', header=None)
            # Replace '23' with 'X', which is how genotype harmonizer calls X
            bim_file.iloc[:, 0].replace(23, 'X', inplace=True)
            # Write new genotype
            bim_file.to_csv(geno_name + '_chr23.bim', sep='\t', header=False, index=False, na_rep='NA')
            # Call genotype harmonizer for X chromosome.
            subprocess.check_output('java -Xmx1g -jar "' + harmonizer_path + '/GenotypeHarmonizer.jar" $* --input '
                                    + geno_name + '_chr23 --ref "'
                                    + os.path.join(vcf_path, vcf_file_names[i])
                                    + '" --refType VCF --hweFilter 0.01 --mafFilter 0.05 --update-id --debug '
                                      '--mafAlign 0 --update-reference-allele --outputType PLINK_BED --output '
                                    + harmonized_geno_names[i], shell=True)
            subprocess.call('rm ' + geno_name + '_chr23.*', shell=True)
            
            id_updates[i] = pd.read_csv(id_update_names[i], sep='\t', header=0,
                                        dtype={'chr': str, 'pos': int, 'originalId': str, 'newId': str})
            snp_logs[i] = pd.read_csv(snp_log_names[i], sep='\t', header=0,
                                   dtype={'chr': str, 'pos': int, 'id': str, 'alleles': str, 'action': str,
                                          'message': str})
        else:
            sys.exit("\u001b[36;1m Something is wrong with the number/name of reference files \u001b[0m")

    # Concatenate all of the id updates into one file.
    all_id_updates = pd.concat([id_updates[0], id_updates[1], id_updates[2], id_updates[3], id_updates[4],
                                id_updates[5], id_updates[6], id_updates[7], id_updates[8], id_updates[9],
                                id_updates[10], id_updates[11], id_updates[12], id_updates[13], id_updates[14],
                                id_updates[15], id_updates[16], id_updates[17], id_updates[18], id_updates[19],
                                id_updates[20], id_updates[21], id_updates[22]])
    # Write list to text file.
    all_id_updates.to_csv('Harmonization_ID_Updates.txt', sep='\t', header=True, index=False)

    # Remove the clutter
    if os.path.getsize('Harmonization_ID_Updates.txt') > 0:
        for f in id_update_names:
            subprocess.call(['rm', f])

    all_snp_logs = pd.concat([snp_logs[0], snp_logs[1], snp_logs[2], snp_logs[3], snp_logs[4], snp_logs[5], snp_logs[6],
                              snp_logs[7], snp_logs[8], snp_logs[9], snp_logs[10], snp_logs[11], snp_logs[12],
                              snp_logs[13], snp_logs[14], snp_logs[15], snp_logs[16], snp_logs[17], snp_logs[18],
                              snp_logs[19], snp_logs[20], snp_logs[21], snp_logs[22]])
    # Write list to text file.
    all_snp_logs.to_csv('Harmonization_SNP_Logs.txt', sep='\t', header=True, index=False)

    # Remove the clutter
    if os.path.getsize('Harmonization_SNP_Logs.txt') > 0:
        for f in snp_log_names:
            subprocess.call(['rm', f])

    # Now for each that was just harmonized, remove all SNPs with an allele (AF) difference > 0.2 since we are going to
    # use a global reference population between study dataset and all superpopulation allele frequencies. IF within 0.2
    # of any superpopulation frequency, keep variant. Frequency file is expected to be a Plink frequency file with the
    # same number of variants as the bim file.
    for i in range(0, len(harmonized_geno_names)):
        # Create plink file
        subprocess.check_output(['plink', '--bfile', harmonized_geno_names[i], '--freq', '--out',
                                 harmonized_geno_names[i]])
        # Read freq file into python
        freq_file = pd.read_csv(freq_file_names[i], sep='\s+', header=0, usecols=[0, 1, 2, 3, 4],
                                dtype={'CHR': int, 'SNP': str, 'A1': str, 'A2': str, 'MAF': float})
        # Rename columns of freq file.
        freq_file.rename(columns={'A1': 'dataset_a1', 'A2': 'dataset_a2', 'MAF': 'dataset_a1_frq'}, inplace=True)
        # Calculate the frequency of the second allele, since it's not given in the freq file.
        freq_file['dataset_a2_frq'] = 1 - freq_file['dataset_a1_frq']
        # Read in bim file.
        bim_file = pd.read_csv(harmonized_geno_names[i] + '.bim', sep='\s+', header=None, usecols=[0, 1, 3],
                               names=['CHR', 'SNP', 'position'])
        # Merge frequency file with bim file to get position for each SNP
        freq_file_with_position = pd.merge(left=freq_file, right=bim_file, how='inner', on=['CHR', 'SNP'])
        # Read in legend file.
        legend_file = pd.read_csv(os.path.join(legend_path, legend_file_names[i]), compression="gzip", sep=" ",
                                  header=0, dtype={'id': str, 'position': int, 'a0': str, 'a1': str, 'TYPE': str,
                                                   'AFR': float, 'AMR': float, 'EAS': float, 'EUR': float, 'SAS': float,
                                                   'ALL': float})
        # Rename columns of legend file
        legend_file.rename(columns={'id': 'reference_id', 'a0': 'reference_a0', 'a1': 'reference_a1'}, inplace=True)

        # To Remove A/T or G/C SNPs in reference file that have an MAF > 40%, first identify which are AT/GC SNPs.
        legend_file['ATGC_SNP'] = np.where(
            ((legend_file['reference_a0'] == 'A') & (legend_file['reference_a1'] == 'T')) |
            ((legend_file['reference_a0'] == 'T') & (legend_file['reference_a1'] == 'A')) |
            ((legend_file['reference_a0'] == 'C') & (legend_file['reference_a1'] == 'G')) |
            ((legend_file['reference_a0'] == 'G') & (legend_file['reference_a1'] == 'C')), 'ATGC', 'Fine')
        # For each superpopulation, create a new column with the MAF. If the AF in that column is less than 0.5, then
        # that is the MAF, if not then 1-AF is MAF
        legend_file['AFR_MAF'] = np.where(legend_file['AFR'] < 0.5, legend_file['AFR'], 1 - legend_file['AFR'])
        legend_file['AMR_MAF'] = np.where(legend_file['AMR'] < 0.5, legend_file['AMR'], 1 - legend_file['AMR'])
        legend_file['EAS_MAF'] = np.where(legend_file['EAS'] < 0.5, legend_file['EAS'], 1 - legend_file['EAS'])
        legend_file['EUR_MAF'] = np.where(legend_file['EUR'] < 0.5, legend_file['EUR'], 1 - legend_file['EUR'])
        legend_file['SAS_MAF'] = np.where(legend_file['SAS'] < 0.5, legend_file['SAS'], 1 - legend_file['SAS'])

        # Create column in legend file with decision about whether to keep or remove SNP, if it is an ATGC SNP where
        # the MAF in all superpopulations is greater than 40%, then remove it.
        legend_file['MAF_Decision'] = np.where((legend_file['ATGC_SNP'] == 'ATGC') & (legend_file['AFR_MAF'] > 0.4) &
                                               (legend_file['AMR_MAF'] > 0.4) & (legend_file['EAS_MAF'] > 0.4) &
                                               (legend_file['EUR_MAF'] > 0.4) & (legend_file['SAS_MAF'] > 0.4),
                                               'Remove', 'Keep')

        # Make new legend file with just the SNPs that pass this threshold.
        legend_file = legend_file[legend_file['MAF_Decision'] == 'Keep']

        legend_file.drop(labels=['ATGC_SNP', 'AFR_MAF', 'AMR_MAF', 'EAS_MAF', 'EUR_MAF', 'SAS_MAF', 'MAF_Decision'],
                         axis=1, inplace=True)

        # Merge freq file with positions with legend file to get overlap. This file contains only SNPs with matches in
        # 1000G
        merged_file = pd.merge(left=freq_file_with_position, right=legend_file, how='inner', on='position')
        merged_file = merged_file[merged_file['TYPE'] == 'Biallelic_SNP']

        # The MAF column in the frq file is the allele frequency of the A1 allele (which is usually minor, but
        # possibly not in my case because I just updated the reference to match 1000G.
        # The AF columns in the legend file are the allele frequencies of the a1 allele in that file.
        # Need to match based on A1 alleles (hopefully this has already been done in the harmonization step, then
        # calculate the allele frequency difference.
        # For each population group, match the alleles where they aren't flipped (i.e. in the same order in house and
        # ref datasets..
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
        # For each population group, match the flipped alleles.
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
        # Add the nonflipped and flipped columns together to find out which files have an allele frequency
        # difference > 0.2
        merged_file['AFR_Diff'] = merged_file['AFR_Match_Diff'] + merged_file['AFR_FlipMatch_Diff']
        merged_file['AMR_Diff'] = merged_file['AMR_Match_Diff'] + merged_file['AMR_FlipMatch_Diff']
        merged_file['EAS_Diff'] = merged_file['EAS_Match_Diff'] + merged_file['EAS_FlipMatch_Diff']
        merged_file['EUR_Diff'] = merged_file['EUR_Match_Diff'] + merged_file['EUR_FlipMatch_Diff']
        merged_file['SAS_Diff'] = merged_file['SAS_Match_Diff'] + merged_file['SAS_FlipMatch_Diff']

        # Delete these columns because we don't need them anymore
        merged_file.drop(labels=['AFR_Match_Diff', 'AFR_FlipMatch_Diff', 'AMR_Match_Diff', 'AMR_FlipMatch_Diff',
                                 'EAS_Match_Diff', 'EAS_FlipMatch_Diff', 'EUR_Match_Diff', 'EUR_FlipMatch_Diff',
                                 'SAS_Match_Diff', 'SAS_FlipMatch_Diff'], axis=1, inplace=True)

        # Make allele freuqency columns numeric.
        merged_file[['AFR_Diff', 'AMR_Diff', 'EAS_Diff', 'EUR_Diff', 'SAS_Diff']] = \
            merged_file[['AFR_Diff', 'AMR_Diff', 'EAS_Diff', 'EUR_Diff', 'SAS_Diff']].apply(pd.to_numeric,
                                                                                            errors='coerce')

        # Make new column 'AF_Decision' where you remove alleles that have allele frequency differences > 0.2 from all
        # population groups.
        merged_file['AF_Decision'] = np.where((merged_file['AFR_Diff'] > 0.2) & (merged_file['AMR_Diff'] > 0.2) &
                                              (merged_file['EAS_Diff'] > 0.2) & (merged_file['EUR_Diff'] > 0.2) &
                                              (merged_file['SAS_Diff'] > 0.2), 'Remove', 'Keep')

        # Drop duplicate SNPs
        merged_file.drop_duplicates(subset=['SNP'], keep=False, inplace=True)

        # Write file for each chromosome of the SNPs that we've removed in this step.
        af_diff_removed_by_chr[i] = merged_file[merged_file['AF_Decision'] == 'Remove']

        # Write list for each chromosome of final SNPs that we are keeping.
        final_snps_by_chr[i] = merged_file[merged_file['AF_Decision'] == 'Keep']
        # Write list for each chromosome, because we're going to use it to filter the chromosomes to create new files.
        final_snps_by_chr[i]['SNP'].to_csv(final_snp_lists[i], sep='\t', header=False, index=False)

        # Make plink files for each chromosomes. Need bed file for merging
        subprocess.check_output(['plink', '--bfile', harmonized_geno_names[i], '--extract', final_snp_lists[i],
                                 '--make-bed', '--out', af_checked_names[i]])

        # Remove extra files that we don't need anymore. These were files that were harmonized, but not checked for
        # allele frequency differences.
        if os.path.getsize(af_checked_names[i] + '.bim') > 0:
            subprocess.call('rm ' + final_snp_lists[i], shell=True)
            subprocess.call('rm ' + harmonized_geno_names[i] + '.*', shell=True)
            
        # Done with one chromosome.

        print('\u001b[36;1m Finished with chr' + str(i + 1) + '\u001b[0m')
    
    # Make a big list of all SNPs removed and all SNPs kept just for reference purposes.
    all_snps_removed = pd.concat([af_diff_removed_by_chr[0], af_diff_removed_by_chr[1], af_diff_removed_by_chr[2],
                                  af_diff_removed_by_chr[3], af_diff_removed_by_chr[4], af_diff_removed_by_chr[5],
                                  af_diff_removed_by_chr[6], af_diff_removed_by_chr[7], af_diff_removed_by_chr[8],
                                  af_diff_removed_by_chr[9], af_diff_removed_by_chr[10], af_diff_removed_by_chr[11],
                                  af_diff_removed_by_chr[12], af_diff_removed_by_chr[13], af_diff_removed_by_chr[14],
                                  af_diff_removed_by_chr[15], af_diff_removed_by_chr[16], af_diff_removed_by_chr[17],
                                  af_diff_removed_by_chr[18], af_diff_removed_by_chr[19], af_diff_removed_by_chr[20],
                                  af_diff_removed_by_chr[21], af_diff_removed_by_chr[22]])
    # Write list to text file.
    all_snps_removed.to_csv('SNPs_Removed_AFCheck.txt', sep='\t', header=True, index=False)

    # Make one big list of all SNPs kept
    all_snps_kept = pd.concat([final_snps_by_chr[0], final_snps_by_chr[1], final_snps_by_chr[2],
                               final_snps_by_chr[3], final_snps_by_chr[4], final_snps_by_chr[5],
                               final_snps_by_chr[6], final_snps_by_chr[7], final_snps_by_chr[8],
                               final_snps_by_chr[9], final_snps_by_chr[10], final_snps_by_chr[11],
                               final_snps_by_chr[12], final_snps_by_chr[13], final_snps_by_chr[14],
                               final_snps_by_chr[15], final_snps_by_chr[16], final_snps_by_chr[17],
                               final_snps_by_chr[18], final_snps_by_chr[19], final_snps_by_chr[20],
                               final_snps_by_chr[21], final_snps_by_chr[22]])
    # Write this list to a text file.
    all_snps_kept.to_csv('SNPs_Kept_AFCheck.txt', sep='\t', header=True, index=False)
    
    # Merge harmonized dataset genotypes
    with open("HouseMergeList.txt", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(af_checked_names)
    subprocess.check_output(['plink', '--merge-list', 'HouseMergeList.txt', '--geno', '0.01', '--make-bed', '--out',
                             geno_name + '_HarmonizedTo1000G'])

    if os.path.getsize(geno_name + '_HarmonizedTo1000G.bim') > 0:
        for i in range(0, len(af_checked_names)):
            subprocess.call('rm ' + str(af_checked_names[i]) + '.*', shell=True)

    else:
        sys.exit("\u001b[36;1m For some reason the house gentoypes did not merge. You should try it manually \u001b[0m")
    
    #Check to make sure the snps are on the same strand as the reference
    # First need to change the chromosome names to match the fasta file so they can match.
    # Read chrX file into pandas
    bim_file = pd.read_csv(geno_name + '_HarmonizedTo1000G.bim', sep='\t', header=None,
                           dtype={0:str, 1:str, 2:int, 3:int, 4: str, 5:str})
    # Replace '23' with 'X', which is how the fasta file calls X
    bim_file.iloc[:, 0].replace('23', 'X', inplace=True)
    # Replace '24' with 'Y', which is how the fasta file calls Y
    bim_file.iloc[:, 0].replace('24', 'Y', inplace=True)
    # Replace '26' with 'MT' which is how the fasta file calls the mitochondrial DNA
    bim_file.iloc[:, 0].replace('26', 'MT', inplace=True)
    # Write new genotype
    bim_file.to_csv(geno_name + '_HarmonizedTo1000G.bim', sep='\t', header=False, index=False, na_rep='NA')

    # Fasta file needs to be unzipped for snpflip to work.
    if os.path.exists(os.path.join(fasta_path, 'human_g1k_v37.fasta')):
        pass
    elif os.path.exists(os.path.join(fasta_path, 'human_g1k_v37.fasta.gz')):
        with gzip.open(os.path.join(fasta_path, 'human_g1k_v37.fasta.gz'), 'rb') as f_in, \
                open(os.path.join(fasta_path, 'human_g1k_v37.fasta'), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        sys.exit("Quitting because I cannot find the fasta file.")

    try:
        # Find where snpflip is.
        for r, d, f in os.walk("C:\\"):
            for files in f:
                if files == "snpflip":
                    snpflip_path = (os.path.join(r, files))
        # Perform flip check.
        subprocess.check_output('python ' + snpflip_path + ' --fasta-genome "'
                                + os.path.join(fasta_path, 'human_g1k_v37.fasta')
                                + '" --bim-file ' + geno_name + '_HarmonizedTo1000G.bim --output-prefix ' + geno_name
                                + '_HarmonizedTo1000G', shell=True)

    except:
        # Import module where I have snpflip
        import genodownload
        # Download snpflip
        genodownload.snpflip()
        # Find where snpflip is:
        for r, d, f in os.walk("C:\\"):
            for files in f:
                if files == "snpflip":
                    snpflip_path = (os.path.join(r, files))
        # Re do
        subprocess.check_output('python ' + snpflip_path + ' --fasta-genome "'
                                + os.path.join(fasta_path, 'human_g1k_v37.fasta') + '" --bim-file ' + geno_name
                                + '_HarmonizedTo1000G.bim --output-prefix ' + geno_name + '_HarmonizedTo1000G',
                                shell=True)

    # If SNPs exist that are on the reverse strand, then flip them.
    # Currently ignores snps that are ambiguous, since I already removed those that would be hard to phase. Could change
    # this later.
    if os.path.getsize(geno_name + '_HarmonizedTo1000G.reverse') > 0:
        subprocess.check_output(['plink', '--bfile', geno_name + '_HarmonizedTo1000G', '--flip',
                                 geno_name + '_HarmonizedTo1000G.reverse', '--make-bed', '--out',
                                 geno_name + '_HarmonizedTo1000G_StrandChecked'])
    # If .reverve doesn't exist, still make a new file to signify that you've done the strand check.
    else:
        subprocess.check_output(['plink', '--bfile', geno_name + '_HarmonizedTo1000G', '--make-bed', '--out',
                                 geno_name + '_HarmonizedTo1000G_StrandChecked'])

    # Finished
    shutil.copy2(geno_name + '_HarmonizedTo1000G_StrandChecked.bed', orig_wd)
    shutil.copy2(geno_name + '_HarmonizedTo1000G_StrandChecked.bim', orig_wd)
    shutil.copy2(geno_name + '_HarmonizedTo1000G_StrandChecked.fam', orig_wd)

    # Change back to original working directory.
    os.chdir(orig_wd)


