import platform

try:
    import colorama
except ImportError:
    import genodownload
    genodownload.getcolorama()

from colorama import init, Fore, Style
init()

# Since we use plink a lot, I'm going to go ahead and set a plink variable with the system-specific plink name.
system_check = platform.system()
if system_check in ("Linux", "Darwin"):
    plink = "./plink"
    rm = "rm "
elif system_check == "Windows":
    plink = 'plink.exe'
    rm = "del "


def merge1000g(geno_name, harmonized_path):
    import os
    import sys
    import subprocess

    # Get original working directory
    orig_wd = os.getcwd()

    # Make sure user has harmonized first.
    print(Fore.BLUE + Style.BRIGHT)
    merge_proceed = input("You must harmonize your data with 1000G before this step. Have you already done this? "
                          "(y/n): ").lower()
    print(Style.RESET_ALL)

    if merge_proceed in ("y", "yes"):
        import csv
        import shutil
        import glob

        try:
            import pandas as pd
        except (ImportError, ModuleNotFoundError):
            import genodownload
            genodownload.getpandas()
            import pandas as pd

        # Create new directory for storing these files.
        if not os.path.exists('Merged_With_1000G'):
            os.makedirs('Merged_With_1000G')

        # Copy plink to new folder.
        if not glob.glob(r'plink*'):
            sys.exit("We need plink to run this part of the script. Please run step 1 if you do not have plink, or "
                     "make sure that it is in this directory.")
        else:
            for file in glob.glob(r'plink*'):
                print(file)
                shutil.copy(file, 'Merged_With_1000G')

        # Making sure they still have the 1000G vcf files. They should if they have just harmonized, but they might
        # have deleted them or something.
        print(Fore.MAGENTA + Style.BRIGHT)
        vcf_exists = input('Have you already downloaded the 1000G Phase3 VCF files? (y/n): ').lower()
        print(Style.RESET_ALL)

        if vcf_exists in ('y', 'yes'):
            # Getting user's path to VCF files
            print(Fore.GREEN)
            vcf_path = input('Please enter the pathname of where your 1000G vcf files are '
                             '(i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\VCF\\ etc.): ')
            print(Style.RESET_ALL)

        elif vcf_exists in ('n', 'no'):
            # Get module where downloading instructions are.
            import genodownload
            # From that module, call download 1000G Phase 3 VCF
            genodownload.vcf_1000g_phase3()
            # Saving VCF path
            vcf_path = os.path.join(os.getcwd(), '1000G_Phase3_VCF')
        else:
            sys.exit("Please answer yes or no. Quitting now because no VCF files.")

        # Names of per chromosome files that we're going to merge into one big file.
        ref_file_names = ['ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' % x for x in
                          range(1, 23)]
        ref_file_names.extend(['ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'])
        chr_1000G_Phase3_names = ['chr%d_1000G_Phase3' % x for x in range(1, 24)]

        # Merge 1000G chr data into one plink formatted file, need to convert from vcf files - but only taking the snps
        # that are in the house dataset
        # Read in SNPs_Kept file from harmonization process
        if os.path.exists(os.path.join(harmonized_path,'SNPs_Kept.txt')):
            house_snps_kept = pd.read_csv(os.path.join(harmonized_path,'SNPs_Kept.txt'), header=0, sep='\t')
            # Keep only the 'SNP' column
            house_snps_kept = house_snps_kept.loc[:, ['SNP']]
            # Write that column to a file to be used by plink
            house_snps_kept.to_csv('Merged_With_1000G/SNPs_Kept_List.txt', sep='\t', header=False, index=False)
        else:
            sys.exit("Quitting because I cannot find a file called 'SNPs_Kept.txt' at "
                     + harmonized_path + ". This is a product of the harmonization process and is necessary for "
                                         "merging with 1000G." )

        # Change to directory where we're going to merge the files.
        os.chdir('Merged_With_1000G')

        # Convert vcf files to plink format.
        for i in range(0, len(ref_file_names)):
            subprocess.check_output([plink,'--vcf', os.path.join(vcf_path, ref_file_names[i]), '--double-id',
                                     '--biallelic-only', 'strict', '--vcf-require-gt', '--extract',
                                     'SNPs_Kept_List.txt', '--make-bed', '--out', chr_1000G_Phase3_names[i]])
        subprocess.call(rm + '*~', shell=True)

        # Create list of files to be merged into one large file.
        with open("1000GMergeList.txt", "w") as f:
            wr = csv.writer(f, delimiter="\n")
            wr.writerow(chr_1000G_Phase3_names)

        # Use plink to merge those files into one large file.
        subprocess.check_output([plink, '--merge-list', '1000GMergeList.txt', '--geno', '0.01', '--make-bed', '--out',
                                 '1000G_Phase3'])

        # Read in log file from merge.
        logfile = pd.DataFrame()
        # Format the log file into a pandas dataframe.
        with open('1000G_Phase3.log', 'r') as f:
            for line in f:
                logfile = pd.concat([logfile, pd.DataFrame([tuple(line.strip().split(" "))])], ignore_index=True)

        # If the logfile contains warnings, write a text file '1000G_MergeWarnings.txt' with the SNPs that threw
        # warnings.
        if logfile[0].str.contains('Warning:').any():
            rsid_warnings = logfile.loc[logfile[0] == 'Warning:', 6].str.split("'", expand=True)
            rsid_warnings[1].dropna(how='any').to_csv('1000G_MergeWarnings.txt', sep='\t', header=False, index=False)

        # If the text file 1000G_MergeWarnings exists...
        if os.path.exists('1000G_MergeWarnings.txt'):
            # If merge warnings and missnps exist, exclude both from 1000G completely (there are plenty of other snps)
            if os.path.exists('1000G_Phase3-merge.missnp'):
               # Read in missnp file
                missnp = pd.read_csv('1000G_Phase3-merge.missnp', sep='\t', header=None)
                # Merge the warning snps with the missnps
                warnings_missnp = pd.concat([rsid_warnings, missnp], axis=0)
                # Drop duplicates and write to a file to be used in plink.
                warnings_missnp[1].dropna(how='any').to_csv('1000G_warnings_missnp.txt', sep='\t',
                                                            header=False, index=False)
                # Remove these snps from plink files and create new plink files.
                for i in range(0, len(chr_1000G_Phase3_names)):
                    subprocess.check_output([plink, '--bfile', chr_1000G_Phase3_names[i], '--exclude',
                                             '1000G_warnings_missnp.txt', '--geno', '0.01', '--make-bed', '--out',
                                             chr_1000G_Phase3_names[i]])
                # Remove old plink files.
                subprocess.call(rm + '*~', shell=True)
                # Retry the merge
                subprocess.check_output([plink, '--merge-list', '1000GMergeList.txt', '--geno', '0.01', '--make-bed',
                                         '--out', '1000G_Phase3'])
                # The merge should be successful this time, but the user should double check.
                print("Successfully merged 1000G, though you should double-check the log file to be sure.")

            else:  # If only merge warnings exist, exclude from 1000G completely
                # Use plink to exclude the merge warning snps.
                for i in range(0, len(chr_1000G_Phase3_names)):
                    subprocess.check_output([plink, '--bfile', chr_1000G_Phase3_names[i], '--exclude',
                                             '1000G_MergeWarnings.txt', '--geno', '0.01', '--make-bed', '--out',
                                             chr_1000G_Phase3_names[i]])
                # Remove old plink files.
                subprocess.call(rm + '*~', shell=True)
                # Try merge again.
                subprocess.check_output([plink, '--merge-list', '1000GMergeList.txt', '--geno', '0.01', '--make-bed',
                                         '--out', '1000G_Phase3'])
                # Merge should be successful this time, but the user should double check.
                print("Successfully merged 1000G, though you should double-check the log file to be sure.")

        # If only the missnps exist, remove them in 1000G.
        elif os.path.exists('1000G_Phase3-merge.missnp') and os.path.getsize('1000G_MergeWarnings.txt') == 0:
            # Use plink to remove the missnps
            for i in range(0, len(chr_1000G_Phase3_names)):
                subprocess.check_output([plink, '--bfile', chr_1000G_Phase3_names[i], '--exclude',
                                         '1000G_Phase3-merge.missnp', '--geno', '0.01', '--make-bed', '--out',
                                         chr_1000G_Phase3_names[i]])
            # Remove old plink files
            subprocess.call(rm + '*~', shell=True)
            # Retry the merge
            subprocess.check_output([plink, '--merge-list', '1000GMergeList.txt', '--geno', '0.01', '--make-bed',
                                     '--out', '1000G_Phase3'])
            # Merge should be successful this time, but the user should double check.
            print("Successfully merged 1000G, though you should double check the log file to be sure.")

        elif os.path.exists('1000G_Phase3.bim'):
            # Merge was successful.
            print("Successfully merged 1000G.")
            # Remove the per chromosome 1000G files, since they're just taking up space now.
            for i in range(0, len(chr_1000G_Phase3_names)):
                subprocess.call(rm + chr_1000G_Phase3_names[i] + '.*', shell=True)

        else: # If merge didn't work for reasons other than merge warnings and missnps.
            print(Fore.RED + Style.BRIGHT)
            sys.exit("Unable to merge 1000G chromosome files. You should try to merge them on your own.")
            print(Style.RESET_ALL)

        ### Merge of house data and 1000G ####
        # Copying harmonized to 1000G files to this folder.
        shutil.copy2(os.path.join(orig_wd, geno_name + '.bed'), os.getcwd())
        shutil.copy2(os.path.join(orig_wd, geno_name + '.bim'), os.getcwd())
        shutil.copy2(os.path.join(orig_wd, geno_name + '.fam'), os.getcwd())

        # Perform initial merge
        subprocess.check_output([plink, '--bfile', geno_name, '--bmerge', '1000G_Phase3', '--geno', '0.01',
                                 '--make-bed', '--out', geno_name + '_1000G'])

        # Read in log file to see if anything went wrong
        logfile = pd.DataFrame()
        # Concatenate log file to pandas dataframe.
        with open(geno_name + '_1000G.log', 'r') as f:
            for line in f:
                logfile = pd.concat([logfile, pd.DataFrame([tuple(line.strip().split(" "))])], ignore_index=True)

        # Check to see if the logfile has warnings
        if logfile[0].str.contains('Warning:').any():
            # Identify the snps that made the warning.
            rsid_warnings = logfile.loc[logfile[0] == 'Warning:', 6].str.split("'", expand=True)
            # Put those SNPs in a file so we can remove them.
            rsid_warnings[1].dropna(how='any').to_csv(geno_name + '_1000G_MergeWarnings.txt', sep='\t', header=False,
                                                      index=False)
        # If there are merge warnings
        if os.path.exists(geno_name + '_1000G_MergeWarnings.txt'):
            if os.path.exists(geno_name + '_1000G-merge.missnp'): #and if there are triallelic snps that need to be
                # flipped (missnp)
                # If merge warnings and missnps exist, exclude warning snps from both 1000G and house dataset
                subprocess.check_output([plink, '--bfile', '1000G_Phase3', '--exclude',
                                         geno_name + '_1000G_MergeWarnings.txt', '--geno', '0.01', '--make-bed',
                                         '--out', '1000G_Phase3'])
                # Flip missnps in house dataset.
                subprocess.check_output([plink, '--bfile', geno_name, '--exclude',
                                         geno_name + '_1000G_MergeWarnings.txt', '--flip',
                                         geno_name + '_1000G-merge.missnp', '--geno', '0.01', '--make-bed', '--out',
                                         geno_name])
                # Remove old files.
                subprocess.call(rm + '*~', shell=True)
                # Retry merge.
                subprocess.check_output([plink, '--bfile', geno_name, '--bmerge', '1000G_Phase3', '--geno', '0.01',
                                         '--make-bed', '--out', geno_name + '_1000G_merge2'])

            else:  # If only mergewarnings exists, exclude warning snps from both 1000G and house dataset.
                subprocess.check_output([plink, '--bfile', '1000G_Phase3', '--exclude',
                                         geno_name + '_1000G_MergeWarnings.txt', '--geno', '0.01', '--make-bed',
                                         '--out', '1000G_Phase3'])
                # Exclude warning snps from house dataset.
                subprocess.check_output([plink, '--bfile', geno_name, '--exclude',
                                         geno_name + '_1000G_MergeWarnings.txt', '--geno', '0.01', '--make-bed',
                                         '--out', geno_name])
                # Remove old plink files.
                subprocess.call(rm + '*~', shell=True)
                # Retry merge
                subprocess.check_output([plink, '--bfile', geno_name, '--bmerge', '1000G_Phase3', '--geno', '0.01',
                                         '--make-bed', '--out', geno_name + '_1000G_merge2'])
        # If only the missnps exist, flip them in the house dataset.
        elif os.path.exists(geno_name + '_1000G-merge.missnp') \
                and not os.path.exists(geno_name + '_1000G_MergeWarnings.txt'):
            # Flip in house dataset.
            subprocess.check_output([plink, '--bfile', geno_name, '--flip', geno_name + '_1000G-merge.missnp',
                                     '--geno', '0.01', '--make-bed', '--out', geno_name])
            # Remove old plink files.
            subprocess.call(rm + '*~', shell=True)
            # Retry the merge.
            subprocess.check_output([plink, '--bfile', geno_name, '--bmerge', '1000G_Phase3', '--geno', '0.01',
                                     '--make-bed', '--out', geno_name + '_1000G_merge2'])

        elif os.path.exists(geno_name + '_1000G.bim'):
            # If mergewarnings and missnips don't exist, then hopefully the merge happened successfully on the first
            # try.
            print("Successfully merged house dataset with 1000G on the first try, though you should double check. I "
                  "can't predict every error.")
            # Copy successfully merged files to original working directory.
            shutil.copy2(geno_name + '_1000G.bed', orig_wd)
            shutil.copy2(geno_name + '_1000G.bim', orig_wd)
            shutil.copy2(geno_name + '_1000G.fam', orig_wd)
            shutil.copy2(geno_name + '_1000G.log', orig_wd)

            # Change back to original working directory.
            os.chdir(orig_wd)

        else:
            print(Fore.RED + Style.BRIGHT)
            sys.exit("The house dataset did not merge properly with 1000G, but not because of SNP merge "
                     "warnings or SNPs that needed to be flipped. I'm sorry, you'll have to perform the merge on your "
                     "own.")
            print(Style.RESET_ALL)

        # If we had to perform a second merge because of one of the reasons above.
        if os.path.exists(geno_name + '_1000G_merge2.log'):
            # Read in log file
            logfile = pd.DataFrame()
            # Convert log file to pandas dataframe
            with open(geno_name + '_1000G_merge2.log', 'r') as f:
                for line in f:
                    logfile = pd.concat([logfile, pd.DataFrame([tuple(line.strip().split(" "))])], ignore_index=True)

            # Check to see if the logfile has warnings
            if logfile[0].str.contains('Warning:').any():
                # Figure out what snps caused the warnings.
                rsid_warnings = logfile.loc[logfile[0] == 'Warning:', 6].str.split("'", expand=True)
                # Write these warnings to a text file for plink to use.
                rsid_warnings[1].dropna(how='any').to_csv(geno_name + '_1000G_merge2_warnings.txt', sep='\t',
                                                          header=False,
                                                          index=False)
            # If warnings exist
            if os.path.exists(geno_name + '_1000G_merge2_warnings.txt'):
                if os.path.exists(geno_name + '_1000G_merge2-merge.missnp'): # And if there are still triallelic snps
                    # Import missnp file to pandas
                    missnp = pd.read_csv(geno_name + '_1000G_merge2-merge.missnp', sep='\t', header=None)
                    # Merge this with merge warnings snps
                    warnings_missnp = pd.concat([rsid_warnings, missnp], axis=0)
                    # Write to text file for plink to use.
                    warnings_missnp[1].dropna(how='any').to_csv(geno_name + '_1000G_merge2_warnings_missnp.txt',
                                                                sep='\t', header=False, index=False)
                    # Remove all of these snps from 1000G dataset.
                    subprocess.check_output([plink, '--bfile', '1000G_Phase3', '--exclude',
                                             geno_name + '_1000G_merge2_warnings_missnp.txt', '--geno', '0.01',
                                             '--make-bed', '--out', '1000G_Phase3'])
                    # Remove all of these snps from house dataset.
                    subprocess.check_output([plink, '--bfile', geno_name, '--exclude',
                                             geno_name + '_1000G_merge2_warnings_missnp.txt', '--geno', '0.01',
                                             '--make-bed', '--out', geno_name])
                    # Remove old plink files.
                    subprocess.call(rm + '*~', shell=True)
                    # Try merge a third time.
                    subprocess.check_output([plink, '--bfile', geno_name, '--bmerge', '1000G_Phase3', '--geno',
                                             '0.01', '--make-bed','--out', geno_name + '_1000G_merge3'])

                else:  # If only merge warnings still exist, exclude from both 1000G and house dataset.
                    # Exclude from 1000G dataset.
                    subprocess.check_output([plink, '--bfile', '1000G_Phase3', '--exclude',
                                             geno_name + '_1000G_merge2_warnings.txt', '--geno', '0.01', '--make-bed',
                                             '--out', '1000G_Phase3'])
                    # Exclude from house dataset
                    subprocess.check_output([plink, '--bfile', geno_name, '--exclude',
                                             geno_name + '_1000G_merge2_warnings.txt', '--geno', '0.01', '--make-bed',
                                             '--out', geno_name])
                    # Remove old plink files.
                    subprocess.call(rm + '*~', shell=True)
                    # Retry merge a third time.
                    subprocess.check_output([plink, '--bfile', geno_name, '--bmerge', '1000G_Phase3', '--geno',
                                             '0.01', '--make-bed', '--out', geno_name + '_1000G_merge3'])
            # If only the missnps still exist, remove them in both datasets.
            elif os.path.exists(geno_name + '_1000G_merge2-merge.missnp') \
                    and not os.path.exists(geno_name + '_1000G_merge2_warnings.txt'):
                # Exclude from house dataset.
                subprocess.check_output([plink, '--bfile', geno_name, '--exclude',
                                         geno_name + '_1000G_merge2-merge.missnp', '--geno', '0.01', '--make-bed',
                                         '--out', geno_name])
                # Exclude from 1000G dataset.
                subprocess.check_output([plink, '--bfile', '1000G_Phase3', '--exclude',
                                         geno_name + '_1000G_merge2-merge.missnp', '--geno', '0.01', '--make-bed',
                                         '--out', '1000G_Phase3'])
                # Remove old files.
                subprocess.call(rm + '*~', shell=True)
                # Retry merge a third time.
                subprocess.check_output([plink, '--bfile', geno_name, '--bmerge', '1000G_Phase3', '--geno', '0.01',
                                         '--make-bed', '--out', geno_name + '_1000G_merge3'])
            # If we don't find warnings or missnps
            elif os.path.exists(geno_name + '_1000G_merge2.bim'):
                # Merge should have happened, but user should check.
                print("Successfully merged house dataset with 1000G on 2nd try, though you should double check. I "
                      "can't predict every error.")
                # Copy merged files to original working directory
                shutil.copy2(geno_name + '_1000G_merge2.bed', orig_wd)
                shutil.copy2(geno_name + '_1000G_merge2.bim', orig_wd)
                shutil.copy2(geno_name + '_1000G_merge2.fam', orig_wd)
                shutil.copy2(geno_name + '_1000G_merge2.log', orig_wd)

                # Change back to original working directory.
                os.chdir(orig_wd)
            # If the merge didn't happen, but not because of merge warnings or missnps.
            else:
                print(Fore.RED + Style.BRIGHT)
                sys.exit("House dataset and 1000G did not successfully merge 2nd time, though not because of merge "
                         "warnings or SNPs that needed to be flipped. You'll have to perform the merge on your own. "
                         "I'm sorry!")
                print(Style.RESET_ALL)

        if os.path.exists(geno_name + '_1000G_merge3.log'):
            # If merge 3 log exists, read it into pandas dataframe.
            logfile = pd.DataFrame()
            with open(geno_name + '_1000G_merge3.log', 'r') as f:
                for line in f:
                    logfile = pd.concat([logfile, pd.DataFrame([tuple(line.strip().split(" "))])], ignore_index=True)

            # If the bim file doesn't exist, try to figure out why.
            if not os.path.exists(geno_name + '_1000G_merge3.bim'):
                # If logfile still has warnings, or there are still missnps, the user will need to identify them and
                # take care of them manually.
                if logfile[0].str.contains('Warning:').any() \
                        and os.path.exists(geno_name + '_1000G_merge3-merge.missnp'):
                    print(Fore.RED + Style.BRIGHT)
                    sys.exit("I'm sorry, the logfile still has warnings, even after removing snps that "
                             "threw errors in the first two tries. There are also still snps with 3+ variants present, "
                             "even after flipping some and removing the ones that the flip didn't solve. You'll have "
                             "to manually deal with these using the merge3 log and the merge3-merge.missnp file.")
                    print(Style.RESET_ALL)
                # Check to see if the logfile still has warnings
                elif logfile[0].str.contains('Warning:').any():
                    print(Fore.RED + Style.BRIGHT)
                    sys.exit("I'm sorry, the logfile still has warnings, even after removing snps that threw errors in "
                             "the first two tries. You'll have to manually deal with these using the merge3 log file.")
                    print(Style.RESET_ALL)
                # If triallelic snps still exist
                elif os.path.exists(geno_name + '_1000G_merge3-merge.missnp'):
                    print(Fore.RED + Style.BRIGHT)
                    sys.exit("I'm sorry, there are still snps with 3+ variants present, even after flipping some and "
                             "removing the ones that the flip didn't solve. You'll have to deal with these manually "
                             "using the merge3-merge.missnp file")
                    print(Style.RESET_ALL)

            # If the bim file exists, then the merge should have happened correctly.
            if os.path.exists(geno_name + '_1000G_merge3.bim'):
                print("House dataset and 1000G merged correctly on 3rd try. Though you should double-check the log "
                      "file, I can't predict every error.")
                # Copy files to original working directory.
                shutil.copy2(geno_name + '_1000G_merge3.bed', orig_wd)
                shutil.copy2(geno_name + '_1000G_merge3.bim', orig_wd)
                shutil.copy2(geno_name + '_1000G_merge3.fam', orig_wd)
                shutil.copy2(geno_name + '_1000G_merge3.log', orig_wd)

                # Change back to original working directory.
                os.chdir(orig_wd)

    # End the program if the user did not harmonize first.
    elif merge_proceed in ('n', 'no'):
        sys.exit("Please harmonize your genotypes with 1000G first.")

    # End the program if the user did not give a correctly formatted answer.
    else:
        sys.exit('Please answer yes or no.')