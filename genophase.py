def phase(geno_name, allocation_name):
    import platform
    import urllib.request
    import os
    import sys
    import pandas as pd
    import numpy as np

    # If it doesn't exist,
    if not os.path.exists('Phasing'):
        os.makedirs('Phasing')

    # Since shapeit only works on linux or mac, we need to first check what system they are on.
    system_check = platform.system()

    if system_check == "Linux":
        # Now we need to know where they have shapeit, if they have it.
        plink = './plink'
        shapeit_exists = input("\u001b[35;1m Do you already have the linux shapeit program unpacked? (y/n):  "
                               "\u001b[0m").lower()
        # If yes, then ask for path of program
        if shapeit_exists in ('yes', 'y'):
            shapeit_path = input("\u001b[36;1m Please tell me the path where you have the shapeit program. "
                                 "i.e. C:\\Users\\Julie White\\Box Sync\\Software\\shapeit\\bin\\ \u001b[0m")
        # If no, download and unpack shapeit.
        elif shapeit_exists in ('no', 'n'):
            print('\u001b[36;1m Downloading shapeit to this directory now. \u001b[0m')
            urllib.request.urlretrieve(
                'https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz',
                'shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz')
            # Making directory to store program
            os.makedirs('Shapeit_v2.12_Linux_Static')
            # Unpacking
            os.system('tar -zxvf shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz -C Shapeit_v2.12_Linux_Static/')
            shapeit_path = os.path.join(os.getcwd(),'Shapeit_v2.12_Linux_Static/bin/')
        else:
            sys.exit('\u001b[36;1m You did not answer "y" or "no" when asked where shapeit was. Exiting now. \u001b[0m')

    # If the user is on a mac
    elif system_check == "Darwin":
        # Ask if they already have shapeit
        plink = './plink'
        shapeit_exists = input("\u001b[34;1m Do you already have the mac shapeit program unpacked? (y/n):  "
                               "\u001b[0m").lower()
        if shapeit_exists in ('yes', 'y'):
            # Ask where shapeit is located.
            shapeit_path = input("\u001b[35;1m Please tell me the path where you have the shapeit program. "
                                 "i.e. C:\\Users\\Julie White\\Box Sync\\Software\\shapeit\\bin\\ \u001b[0m")
        elif shapeit_exists in ('no', 'n'):
            print('\u001b[35;1m Downloading shapeit now. \u001b[0m')
            # Download shapeit
            urllib.request.urlretrieve(
                'https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.MacOSX.tgz',
                'shapeit.v2.r837.MacOSX.tgz')
            # Create directory for shapeit.
            os.makedirs('Shapeit_v2.20_Mac')
            # Untar shapeit to that directory.
            os.system('tar -zxvf shapeit.v2.r837.MacOSX.tgz -C Shapeit_v2.20_Mac/')
            shapeit_path = os.path.join(os.getcwd(), 'Shapeit_v2.20_Mac/bin/')
        else:
            sys.exit('\u001b[35;1m You did not answer "y" or "no" when asked where shapeit was. Exiting now. \u001b[0m')

    # If they are running this on a windows machine, they cannot proceed because shapeit is *nix only.
    elif system_check == ("Windows"):
        sys.exit("\u001b[35;1m I'm sorry, you need access to a linux or mac system to make this part work. If you have "
                 "access to the Penn State clusters, you should run this script from there (they are linux). \u001b[0m")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("\u001b[35;1m I cannot detect the system you are working on. Exiting now. \u001b[0m")

    # Ask the user if they already have the 1000G Phase 3 Hap/Legend/Sample files.
    ref_exists = input('\u001b[35;1m Have you already downloaded the 1000G Phase 3 Hap/Legend/Sample files? (y/n): '
                       '\u001b[0m').lower()
    # If yes
    if ref_exists in ('y', 'yes'):
        # Ask the user where the ref files are.
        ref_path = input('\u001b[35;1m Please enter the pathname of where your 1000G hap/legend/sample files are '
                         '(i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\ etc.): \u001b[0m')
    # If no
    elif ref_exists in ('n', 'no'):
        # Get genodownload module
        import genodownload
        # Call download HapLegendSample command
        genodownload.hls_1000g_phase3()
        # Saving legend path
        ref_path = os.path.join(os.getcwd(), '1000G_Phase3_HapLegendSample')
    # If user gives non-recognized answer.
    else:
        sys.exit('Please answer yes or no. Quitting now because no hap/legend/sample files.')

    # Ask if the user is on the cluster right now to determine if we should submit the files for them
    on_cluster = input('Are you currently running this from the Penn State ACI-B cluster? If yes, I can submit the jobs'
                       ' for you. If not, you will need to submit the files yourself').lower()

    # Names of files needed.
    genetic_map_names = ['genetic_map_chr%d_combined_b37.txt' % x for x in range(1, 23)]
    genetic_map_names.extend(['genetic_map_chrX_nonPAR_combined_b37.txt'])
    hap_names = ['1000GP_Phase3_chr%d.hap.gz' % x for x in range(1,23)]
    hap_names.extend(['1000GP_Phase3_chrX_NONPAR.hap.gz'])
    legend_names = ['1000GP_Phase3_chr%d.legend.gz' % x for x in range(1,23)]
    legend_names.extend(['1000GP_Phase3_chrX_NONPAR.legend.gz'])
    geno_bed_names = ['Phasing/' + geno_name + '.chr%d.bed' % x for x in range(1,24)]
    geno_bim_names = ['Phasing/' + geno_name + '.chr%d.bim' % x for x in range(1, 24)]
    geno_fam_names = ['Phasing/' + geno_name + '.chr%d.fam' % x for x in range(1, 24)]
    check_log_names = ['Phasing/' + geno_name + '_PhaseCheck.chr%d' % x for x in range(1, 24)]
    strand_exclude_names = ['Phasing/' + geno_name + '_PhaseCheck.chr%d.snp.strand.exclude' % x for x in
                                range(1,24)]
    ind_me_names = ['Phasing/' + geno_name + '_PhaseCheck.chr%d.ind.me' % x for x in range(1, 24)]
    snp_me_names = ['Phasing/' + geno_name + '_PhaseCheck.chr%d.snp.me' % x for x in range(1, 24)]
    output_max_names = ['Phasing/' + geno_name + '_PhasedTo1000G.chr%d' % x for x in range(1,24)]
    output_log_names = ['Phasing/' + geno_name + '_PhasedTo1000G.chr%d' % x for x in range(1,24)]

    # Use plink to set mendel errors to missing.
    os.system(plink + ' --bfile ' + geno_name + ' --me 1 1 --set-me-missing --make-bed --out ' + geno_name)
    # Remove old files
    os.system('rm *~')

    # Perform phasing check per chromosome.
    for i in range(0,23):
        # Split the geno file into separate chromosomes and put it in the Phasing folder.
        os.system(plink + ' --bfile ' + geno_name + ' --chr ' + str(i + 1) + ' --make-bed --out Phasing/'
                  + geno_name + '.chr' + str(i + 1))
        # If on autosomes
        if i < 22:
            # Perform phasing check. We are telling it to ignore the pedigree information for now because we don't want
            # it to error out when it encounters a pedigree and I'm going to use the --duohmm flag later.
            os.system(os.path.join(shapeit_path,'shapeit') + ' -check --input-bed ' + geno_bed_names[i] + ' '
                      + geno_bim_names[i] + ' ' + geno_fam_names[i] + ' --noped --input-map '
                      + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                      + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i]) + ' '
                      + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --output-log ' + check_log_names[i])
        # If working on the X chromosome
        if i == 22:
            # In the 1000G sample file, need to change male to 1 and female to 2.
            ref_sample = pd.read_csv(os.path.join(ref_path, '1000GP_Phase3.sample'), sep=" ", header = 0)
            # Replace female with 2
            ref_sample.iloc[:,3].replace('female', '2', inplace=True)
            # Replace male with 1
            ref_sample.iloc[:,3].replace('male', '1', inplace=True)
            # Write file
            ref_sample.to_csv(os.path.join(ref_path, '1000GP_Phase3.sample'))

            # Perform phasing check. Phasing chrX specifically considers all people unrelated. They might change this
            # later.
            os.system(os.path.join(shapeit_path, 'shapeit') + ' -check --chrX --input-bed ' + geno_bed_names[i] + ' '
                      + geno_bim_names[i] + ' ' + geno_fam_names[i] + ' --input-map '
                      + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                      + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i]) + ' '
                      + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --output-log ' + check_log_names[i])

    # Now that check is performed, we're on the lookout for these files:
        #  myLogFile.snp.strand.exclude that  gives a list of the physical positions of all Case1 and Case3 problems
            # found for easy exclusion using --exclude-snp option.
        #  gwas.checks.ind.me contains the Mendel errors at the individual level (N lines).
        #  gwas.checks.snp.me contains the Mendel errors at the SNP level (L lines).
    # Additional files just for chrX:
        # gwas.phased.snp.hh reports the number of haploid heterozygous per variant site.
        # gwas.phased.ind.hh reports the number of haploid heterozygous per sample.
    # The mendel errors should not happen because we set them as missing in plink, but putting them in here just in case.

    # Prepare files for phasing and submit to cluster.
    for i in range(0,23):
        if i < 22:
            # If a log file with snp.strand exclude is produced, then these SNPs need to be removed.
            if os.path.exists(strand_exclude_names[i]):
                # Read in strand_exclude file to just get positions to remove.
                strand_exclude = pd.read_csv(strand_exclude_names[i], sep=" ", header=None)
                log = ['a']
            else:
                log = []

            # If a *.ind.me file is produced and the error rates are nonzero, tell the user because they should
            # investigate these people. Could indicate that their family assignments are incorrect.
            if os.path.exists(ind_me_names[i]):
                # Read in ind file.
                snp_ind_file = pd.read_csv(ind_me_names[i], sep='\t', header = 1, dtype = {0:str,1:str,2:float,3:str,
                                                                                           4:float,5:float,6:int})
                father_me_errors = snp_ind_file[snp_ind_file[2] > 0]
                mother_me_errors = snp_ind_file[snp_ind_file[4] > 0]
                if (len(father_me_errors) > 0) or (len(mother_me_errors) > 0):
                    print("\u001b[31;1m Your files have people with non-zero mendel errors. You should investigate "
                          "the " + ind_me_names[i] + ' file and take a careful look at the people with high values in '
                                                     'the father_mendel & mother_mendel column. This result suggests '
                                                     'that your paternity/maternity assignment could be incorrect. '
                                                     '\u001b[0m')
            else:
                pass

            # If a *.snp.me file is produced, see if there are any snps with high mendel error rates.
            if os.path.exists(snp_me_names[i]):
                # Read in snp.me file
                snp_me_file = pd.read_csv(snp_me_names[i], sep='\t', header=None, skiprows=1,
                                          dtype={0: str, 1: int, 2: float, 3: float})
                # Create mendel error column where we see which snps have a error rate of > 0.05
                snp_me_file['MendelError'] = np.where((snp_me_file[2] / snp_me_file[3]) > 0.05, 'Yes', 'No')
                # Create new dataframe with only the positions where Mendel Error was yes
                me_exclude = snp_me_file[snp_me_file[1]['MendelError'] == 'Yes']
                # If this is non-empty, make note.
                if len(me_exclude) > 0:
                    log.extend(['b'])
                # If it is empty (i.e. no Mendel Errors exist) do nothing.
                else:
                    pass
            else:
                pass

            # If both snp.strand.exclude exists and there are mendel error snps:
            if log == ['a', 'b']:
                # Concatenate the two lists of snps to exclude
                snp_exclude = pd.concat([strand_exclude, me_exclude], axis=0)
                # Drop duplicates
                snp_exclude.drop_duplicates(keep='first', inplace=True)
                # Write file
                snp_exclude.to_csv(check_log_names[i] + '.ExcludeSnps',
                                   sep="\t", header=False, index=False)
                # Save name of file we just created
                snp_exclude_name = check_log_names[i] + '.ExcludeSnps'

            # If just the snp.strand.exclude exists.
            elif log == ['a']:
                # Since it was just the strand problem log that was written, the exclude snp name is the same as that
                # log file name
                snp_exclude_name = strand_exclude_names[i]

            # If just the snp.me errors exist.
            elif log == ['b']:
                # Write file
                me_exclude.to_csv(check_log_names[i] + '.ExcludeSnps', sep="\t", header=False, index=False)
                # Save name of file we just created.
                snp_exclude_name = check_log_names[i] + '.ExcludeSnps'
            # If these files don't exist, do nothing for now.
            elif log == []:
                pass
            else:
                sys.exit("I shouldn't have gotten here. There is something wrong with the way the script is checking "
                         "for the phase check files.")

            # If I just filled the snp_exclude_name variable with a name, then we have snps to remove.
            if 'snp_exclude_name' in locals():
                # Write pbs script
                with open('Phasing/' + geno_name + '_PhaseScript.chr' + str(i + 1) + '.pbs', 'w') as file:
                    file.write('#!/bin/bash\n'
                               '#PBS -l walltime=30:00:00\n'
                               '#PBS -l nodes=1:ppn=8\n'
                               '#PBS -l pmem=7gb\n'
                               '#PBS -A ' + allocation_name +
                               '\n'
                               '#PBS -o Phasing/\n'
                               '#PBS -j oe\n'
                               'cd $PBS_O_WORKDIR\n'
                               '\n' +
                               os.path.join(shapeit_path, 'shapeit') + ' --input-bed ' + geno_bed_names[i] + ' '
                               + geno_bim_names[i] + ' ' + geno_fam_names[i] + ' --duohmm --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --exclude-snp '
                               + snp_exclude_name + ' --output-max ' + output_max_names[i] + ' --output-log '
                               + output_log_names[i])
            # If snp_exclude_names isn't filled, then phase with all SNPs.
            else:
                # Write pbs file.
                with open ('Phasing/' + geno_name + '_PhaseScript.chr' + str(i+1) + '.pbs', 'w') as file:
                    file.write('#!/bin/bash\n'
                               '#PBS -l walltime=30:00:00\n'
                               '#PBS -l nodes=1:ppn=8\n'
                               '#PBS -l pmem=7gb\n'
                               '#PBS -A ' + allocation_name +
                               '\n'
                               '#PBS -o Phasing/\n'
                               '#PBS -j oe\n'
                               'cd $PBS_O_WORKDIR\n'
                               '\n' +
                               os.path.join(shapeit_path, 'shapeit') + ' --input-bed ' + geno_bed_names[i] + ' '
                               + geno_bim_names[i] + ' ' + geno_fam_names[i] + ' --duohmm --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --output-max '
                               + output_max_names[i] + ' --output-log ' + output_log_names[i])

        # Need to phase chrX specially.
        if i == 22:
            # If snp.strand.exclude exists:
            if os.path.exists(strand_exclude_names[i]):
                # Read in strand_exclude file to just get positions to remove.
                strand_exclude = pd.read_csv(strand_exclude_names[i], sep=" ", header=None)
                log = ['a']
            else:
                log = []

            # If a *.ind.me file is produced and the error rates are nonzero, tell the user because they should
            # investigate these people. Could indicate that their family assignments are incorrect.
            if os.path.exists(ind_me_names[i]):
                # Read in ind file.
                snp_ind_file = pd.read_csv(ind_me_names[i], sep='\t', header=1,
                                           dtype={0: str, 1: str, 2: float, 3: str,
                                                  4: float, 5: float, 6: int})
                father_me_errors = snp_ind_file[snp_ind_file[2] > 0]
                mother_me_errors = snp_ind_file[snp_ind_file[4] > 0]
                if (len(father_me_errors) > 0) or (len(mother_me_errors) > 0):
                    print("\u001b[31;1m Your files have people with non-zero mendel errors. You should investigate "
                          "the " + ind_me_names[
                              i] + ' file and take a careful look at the people with high values in '
                                   'the father_mendel & mother_mendel column. This result suggests '
                                   'that your paternity/maternity assignment could be incorrect. '
                                   '\u001b[0m')
                else:
                    pass

            # If a *.snp.me file is produced, see if there are any snps with high mendel error rates.
            if os.path.exists(snp_me_names[i]):
                # Read in snp.me file
                snp_me_file = pd.read_csv(snp_me_names[i], sep='\t', header=None, skiprows=1,
                                          dtype={0: str, 1: int, 2: float, 3: float})
                # Create mendel error column where we see which snps have a error rate of > 0.05
                snp_me_file['MendelError'] = np.where((snp_me_file[2] / snp_me_file[3]) > 0.05, 'Yes', 'No')
                # Create new dataframe with only the positions where Mendel Error was yes
                me_exclude = snp_me_file[snp_me_file[1]['MendelError'] == 'Yes']
                # If this is non-empty, make note.
                if len(me_exclude) > 0:
                    log.extend(['b'])
                # If it is empty (i.e. no Mendel Errors exist) do nothing.
                else:
                    pass
            else:
                pass

            # If a *.snp.hh file is produced, see if there are snps with high haploid heterozygosity rate.
            if os.path.exists(check_log_names[i] + '.snp.hh'):
                # Read in haploid het file
                snp_hh_file = pd.read_csv(check_log_names[i] + '.snp.hh', sep=" ", header=None, skiprows=1,
                                          dtype={0: str, 1: int, 2: float, 3: float})
                # Create a new column with 'Yes' for snps with mendel error > 1%
                snp_hh_file['HHError'] = np.where((snp_hh_file[2] / snp_hh_file[3]) > 0.01, 'Yes', 'No')
                # Create new dataframe with only the positions with high HH rates
                snp_hh_exclude = snp_hh_file[snp_hh_file[1]['HHError'] == 'Yes']
                # If this is non-empty, make note:
                if len(snp_hh_exclude) > 0:
                    log.extend(['c'])
                # If it is empty (i.e. no HH errors exist) do nothing
                else:
                    pass
            else:
                pass

            # If a *.ind.hh file is produced, see if there are males with high haploid heterozygosity rate.
            if os.path.exists(check_log_names[i] + '.ind.hh')
                # Read in file of people to check.
                ind_hh_file = pd.read_csv(check_log_names[i] + '.ind.hh', sep = " ", header = None, skiprows=1,
                                          dtype={0: str, 1: int, 2: float, 3: float})
                # Create new column with info on whether the ind has mendel error > 1%
                ind_hh_file['HHError'] = np.where((ind_hh_file[2] / ind_hh_file[3]) > 0.01, 'Yes', 'No')
                # Create new dataframe of people with high mendel errors.
                ind_hh_exclude = ind_hh_file[ind_hh_file[1]['HHError'] == 'Yes']
                # If this list is non-zero, write file of individuals to exclude.
                if len(ind_hh_exclude) > 0:
                    ind_hh_exclude.to_csv(check_log_names[i] + '.ind.hh.exclude')
                # If not, then do nothing
                else:
                    ind_hh_exclude = []
            else:
                ind_hh_exclude = []

            # If snp.strand.exclude, me.snp, and snp.hh exists:
            if log == ['a', 'b', 'c']:
                # Concatenate the three lists of snps to exclude
                snp_exclude = pd.concat([strand_exclude, me_exclude, snp_hh_exclude], axis=0)
                # Drop duplicates
                snp_exclude.drop_duplicates(keep='first', inplace=True)
                # Write file
                snp_exclude.to_csv(check_log_names[i] + '.ExcludeSnps', sep="\t", header=False, index=False)
                # Save name of file we just created
                snp_exclude_name = check_log_names[i] + '.ExcludeSnps'

            # If snp.strand.exclude and me.snp exist:
            elif log == ['a', 'b']:
                # Concatenate the two lists of snps to exclude
                snp_exclude = pd.concat([strand_exclude, me_exclude], axis=0)
                # Drop duplicates
                snp_exclude.drop_duplicates(keep='first', inplace=True)
                # Write file
                snp_exclude.to_csv(check_log_names[i] + '.ExcludeSnps', sep="\t", header=False, index=False)
                # Save name of file we just created
                snp_exclude_name = check_log_names[i] + '.ExcludeSnps'


            # If snp.strand.exclude and .snp.hh exist:
            elif log == ['a', 'c']:
                # Concatenate the two lists of snps to exclude
                snp_exclude = pd.concat([strand_exclude, snp_hh_exclude], axis=0)
                # Drop duplicates
                snp_exclude.drop_duplicates(keep='first', inplace=True)
                # Write file
                snp_exclude.to_csv(check_log_names[i] + '.ExcludeSnps',sep="\t", header=False, index=False)
                # Save name of file we just created
                snp_exclude_name = check_log_names[i] + '.ExcludeSnps'

            # If me.snp and snp.hh exist:
            elif log == ['b', 'c']:
                # Concatenate the two lists of snps to exclude
                snp_exclude = pd.concat([strand_exclude, snp_hh_exclude], axis=0)
                # Drop duplicates
                snp_exclude.drop_duplicates(keep='first', inplace=True)
                # Write file
                snp_exclude.to_csv(check_log_names[i] + '.ExcludeSnps', sep="\t", header=False, index=False)
                # Save name of file we just created
                snp_exclude_name = check_log_names[i] + '.ExcludeSnps'

            # If just the snp.strand.exclude exists.
            elif log == ['a']:
                # Since it was just the strand problem log that was written, the exclude snp name is the same as that
                # log file name
                snp_exclude_name = strand_exclude_names[i]

            # If just the snp.me errors exist.
            elif log == ['b']:
                # Write file
                me_exclude.to_csv(check_log_names[i] + '.ExcludeSnps', sep="\t", header=False, index=False)
                # Save name of file we just created.
                snp_exclude_name = check_log_names[i] + '.ExcludeSnps'

            # If just the snp.hh errors exit.
            elif log == ['c']:
                # Write file
                snp_hh_exclude.to_csv(check_log_names[i] + '.ExcludeSnps', sep="\t", header=False, index=False)
                # Save name of file we just created.
                snp_exclude_name = check_log_names[i] + '.ExcludeSnps'
            # If these files don't exist, do nothing for now.
            elif log == []:
                pass
            else:
                sys.exit(
                    "I shouldn't have gotten here. There is something wrong with the way the script is checking "
                    "for the phase check files.")

            # If I just filled the snp_exclude_name variable, then we have snps to remove.
            # If there are people in ind_hh_exclude, then we have people to remove.
            if ('snp_exclude_name' in locals()) and (len(ind_hh_exclude) > 0):
                # Write pbs script
                with open('Phasing/' + geno_name + '_PhaseScript.chr' + str(i + 1) + '.pbs', 'w') as file:
                    file.write('#!/bin/bash\n'
                               '#PBS -l walltime=30:00:00\n'
                               '#PBS -l nodes=1:ppn=8\n'
                               '#PBS -l pmem=7gb\n'
                               '#PBS -A ' + allocation_name +
                               '\n'
                               '#PBS -o Phasing/\n'
                               '#PBS -j oe\n'
                               'cd $PBS_O_WORKDIR\n'
                               '\n' +
                               os.path.join(shapeit_path, 'shapeit') + ' --input-bed ' + geno_bed_names[i] + ' '
                               + geno_bim_names[i] + ' ' + geno_fam_names[i] + ' --chrX --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --exclude-snp '
                               + snp_exclude_name + ' --exclude-ind ' + check_log_names[i] + '.ind.hh.exclude'
                               + ' --output-max ' + output_max_names[i] + ' --output-log ' + output_log_names[i])
            # If only snp_exclude_name variable is filled, then we have snps to remove.
            elif 'snp_exclude_name' in locals():
                # Write pbs script
                with open('Phasing/' + geno_name + '_PhaseScript.chr' + str(i + 1) + '.pbs', 'w') as file:
                    file.write('#!/bin/bash\n'
                               '#PBS -l walltime=30:00:00\n'
                               '#PBS -l nodes=1:ppn=8\n'
                               '#PBS -l pmem=7gb\n'
                               '#PBS -A ' + allocation_name +
                               '\n'
                               '#PBS -o Phasing/\n'
                               '#PBS -j oe\n'
                               'cd $PBS_O_WORKDIR\n'
                               '\n' +
                               os.path.join(shapeit_path, 'shapeit') + ' --input-bed ' + geno_bed_names[i] + ' '
                               + geno_bim_names[i] + ' ' + geno_fam_names[i] + ' --chrX --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --exclude-snp '
                               + snp_exclude_name + ' --output-max ' + output_max_names[i] + ' --output-log '
                               + output_log_names[i])
            # If there are only people in ind_hh_exclude, then we have people to remove.
            elif len(ind_hh_exclude) > 0:
                # Write pbs script
                with open('Phasing/' + geno_name + '_PhaseScript.chr' + str(i + 1) + '.pbs', 'w') as file:
                    file.write('#!/bin/bash\n'
                               '#PBS -l walltime=30:00:00\n'
                               '#PBS -l nodes=1:ppn=8\n'
                               '#PBS -l pmem=7gb\n'
                               '#PBS -A ' + allocation_name +
                               '\n'
                               '#PBS -o Phasing/\n'
                               '#PBS -j oe\n'
                               'cd $PBS_O_WORKDIR\n'
                               '\n' +
                               os.path.join(shapeit_path, 'shapeit') + ' --input-bed ' + geno_bed_names[i] + ' '
                               + geno_bim_names[i] + ' ' + geno_fam_names[i] + ' --chrX --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --exclude-ind '
                               + check_log_names[i] + '.ind.hh.exclude' + ' --output-max ' + output_max_names[i]
                               + ' --output-log ' + output_log_names[i])
            # If none of these are filled, then phase with all SNPs.
            else:
                # Write pbs file.
                with open('Phasing/' + geno_name + '_PhaseScript.chr' + str(i + 1) + '.pbs', 'w') as file:
                    file.write('#!/bin/bash\n'
                               '#PBS -l walltime=30:00:00\n'
                               '#PBS -l nodes=1:ppn=8\n'
                               '#PBS -l pmem=7gb\n'
                               '#PBS -A ' + allocation_name +
                               '\n'
                               '#PBS -o Phasing/\n'
                               '#PBS -j oe\n'
                               'cd $PBS_O_WORKDIR\n'
                               '\n' +
                               os.path.join(shapeit_path, 'shapeit') + ' --input-bed ' + geno_bed_names[i] + ' '
                               + geno_bim_names[i] + ' ' + geno_fam_names[i] + ' --chrX --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --output-max '
                               + output_max_names[i] + ' --output-log ' + output_log_names[i])

    #I If the user is currently on the cluster, then submit the pbs files to start running.
    if on_cluster in ('y', 'yes'):
        for i in range(0,23):
            # Submit to cluster
            os.system('qsub Phasing/' + geno_name + '_PhaseScript.chr' + str(i + 1) + '.pbs')
    elif on_cluster in ('no', 'no'):
        sys.exit('Since you are not on the Penn State cluster right now, you should transfer the genotype files, 1000G '
                 'genetic maps, 1000G legend files, 1000G hap files, and 1000G sample file, shapeit, and the phasing '
                 'pbs files to the cluster. You can submit the files by using qsub name_of_file.pbs')
    else:
        sys.exit(" You didn't answer 'yes' or 'no' when I asked whether you were on the cluster or not, so I'm assuming "
                 "you're not. You should transfer the genotype files, 1000G genetic maps, 1000G legend files, 1000G hap "
                 "files, and 1000G sample file, shapeit, and the phasing pbs files to the cluster. You can submit the "
                 "files by using qsub name_of_file.pbs")