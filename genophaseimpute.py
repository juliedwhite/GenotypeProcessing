import platform
import os
import sys
import subprocess
from os.path import expanduser
try:
    import colorama
except ImportError:
    import genodownload
    genodownload.getcolorama()

from colorama import init, Fore, Style
init()

home = expanduser("~")
bindir = os.path.join(home, 'software', 'bin')
system_check = platform.system()

# If they are running this on a windows machine, they cannot proceed because shapeit, bcftools, vcftools are *nix only.
if system_check == "Windows":
    print(Fore.RED + Style.BRIGHT)
    sys.exit("I'm sorry, you need access to a linux or mac system to make this part work. If you have "
             "access to the Penn State clusters, you should run this script from there (they are linux).")
    print(Style.RESET_ALL)

# Since we use plink a lot, I'm going to go ahead and set a plink variable with the system-specific plink name.
plink = "plink"
rm = "rm "

# Determine if they have plink, if not download it.
if os.path.exists(os.path.join(bindir, plink)):
    pass
else:
    import genodownload
    genodownload.plink()


def phase(geno_name, allocation_name):
    try:
        import pandas as pd
    except (ImportError, ModuleNotFoundError):
        import genodownload
        genodownload.getpandas()
        import pandas as pd

    try:
        import numpy as np
    except (ImportError, ModuleNotFoundError):
        import genodownload
        genodownload.getnumpy()
        import numpy as np

    # Ask if the user is on the cluster right now to determine if we should submit the files for them
    print(Fore.BLUE + Style.BRIGHT)
    on_cluster = input('Are you currently running this from the Penn State ACI-B cluster? If yes, I can submit the jobs'
                       ' for you. If not, you will need to submit the files yourself. (y/n): ').lower()
    print(Style.RESET_ALL)

    # If it doesn't exist, create Phasing folder
    if not os.path.exists('Phasing'):
        os.makedirs('Phasing')

    # Determine if they have shapeit, if not download it.
    if os.path.exists(os.path.join(bindir, 'shapeit')):
        pass
    else:
        import genodownload
        genodownload.shapeit()

    # Ask the user if they already have the 1000G Phase 3 Hap/Legend/Sample files.
    print(Fore.CYAN)
    ref_exists = input('Have you already downloaded the 1000G Phase 3 Hap/Legend/Sample files? (y/n): ').lower()
    print(Style.RESET_ALL)
    # If yes
    if ref_exists in ('y', 'yes'):
        # Ask the user where the ref files are.
        print(Fore.BLUE + Style.BRIGHT)
        ref_path = input('Please enter the pathname of where your 1000G hap/legend/sample files are '
                         '(i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\ etc.): ')
        print(Style.RESET_ALL)
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
    output_names = ['Phasing/' + geno_name + '_PhasedTo1000G.chr%d' % x for x in range(1,24)]
    output_vcf_names = ['Phasing/' + geno_name + '_PhasedTo1000G.chr%d.vcf' % x for x in range(1,24)]
    output_vcf_log = ['Phasing/' + geno_name + '_PhasedTo1000G.chr%d.vcf.log' % x for x in range(1,24)]

    # Use plink to set mendel errors to missing.
    subprocess.check_output([plink, '--bfile', geno_name, '--me', '1', '1', '--set-me-missing', '--make-bed', '--out',
                             geno_name])
    # Remove old files
    subprocess.call(rm + '*~', shell=True)

    # Perform phasing check per chromosome.
    for i in range(0,23):
        # Split the geno file into separate chromosomes and put it in the Phasing folder.
        subprocess.check_output([plink, '--bfile', geno_name, '--chr', str(i + 1), '--make-bed', '--out',
                                 'Phasing/' + geno_name + '.chr' + str(i+1)])
        # If on autosomes
        if i < 22:
            # Perform phasing check.
            subprocess.call(['shapeit', '-check', '--input-bed', geno_bed_names[i], geno_bim_names[i],
                             geno_fam_names[i], '--input-map', os.path.join(ref_path, genetic_map_names[i]),
                             '--input-ref', os.path.join(ref_path, hap_names[i]),
                             os.path.join(ref_path, legend_names[i]), os.path.join(ref_path, '1000GP_Phase3.sample'),
                             '--output-log', check_log_names[i]])
            print("Pre-phasing check done on chr" + str(i + 1))

        # If working on the X chromosome
        if i == 22:
            # In the 1000G sample file, need to change male to 1 and female to 2.
            ref_sample = pd.read_csv(os.path.join(ref_path, '1000GP_Phase3.sample'), sep=" ", header = 0)
            # Replace female with 2
            ref_sample.iloc[:,3].replace('female', '2', inplace=True)
            # Replace male with 1
            ref_sample.iloc[:,3].replace('male', '1', inplace=True)
            # Write file
            ref_sample.to_csv(os.path.join(ref_path, '1000GP_Phase3.sample'), sep = " ", header = True, index = False)

            # Perform phasing check. Phasing chrX specifically considers all people unrelated. They might change this
            # later.
            subprocess.check_output(['shapeit', '-check', '--chrX', '--input-bed', geno_bed_names[i], geno_bim_names[i],
                                     geno_fam_names[i], '--input-map', os.path.join(ref_path, genetic_map_names[i]),
                                     '--input-ref', os.path.join(ref_path, hap_names[i]),
                                     os.path.join(ref_path, legend_names[i]),
                                     os.path.join(ref_path, '1000GP_Phase3.sample'), '--output-log',
                                     check_log_names[i]])
            print("Pre-phasing check done on chr" + str(i + 1))

    # Now that check is performed, we're on the lookout for these files:
        #  myLogFile.snp.strand.exclude that  gives a list of the physical positions of all Case1 and Case3 problems
            # found for easy exclusion using --exclude-snp option.
        #  gwas.checks.ind.me contains the Mendel errors at the individual level (N lines).
        #  gwas.checks.snp.me contains the Mendel errors at the SNP level (L lines).
    # Additional files just for chrX:
        # gwas.phased.snp.hh reports the number of haploid heterozygous per variant site.
        # gwas.phased.ind.hh reports the number of haploid heterozygous per sample.
    # The mendel errors should not happen because we set them as missing in plink, but putting them in here just in
    # case.

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
                ind_me_file = pd.read_csv(ind_me_names[i], sep='\t', header = 0, dtype = {0:str,1:str,2:float,3:str,
                                                                                           4:float,5:float,6:int})
                father_me_errors = ind_me_file[ind_me_file['father_mendel'] > 0]
                mother_me_errors = ind_me_file[ind_me_file['mother_mendel'] > 0]
                if (len(father_me_errors) > 0) or (len(mother_me_errors) > 0):
                    print("Your files have people with non-zero mendel errors. You should investigate the "
                          + ind_me_names[i] + ' file and take a careful look at the people with high values in the '
                                              'father_mendel & mother_mendel column. This result suggests that your '
                                              'paternity/maternity assignment could be incorrect.')
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
                me_exclude = snp_me_file[snp_me_file['MendelError'] == 'Yes'][1]
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
                               'shapeit --input-bed ' + geno_bed_names[i] + ' ' + geno_bim_names[i] + ' '
                               + geno_fam_names[i] + ' --duohmm --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --exclude-snp '
                               + snp_exclude_name + ' --output-max ' + output_names[i] + ' --output-log '
                               + output_names[i]
                               + '\n'
                               + 'shapeit -convert --input-haps ' + output_names[i] + '.haps '
                               + output_names[i] + '.sample --output-vcf ' + output_vcf_names[i] + ' --output-log '
                               + output_vcf_log[i])
            # If snp_exclude_names isn't filled, then phase with all SNPs.
            else:
                # Write pbs file.
                with open('Phasing/' + geno_name + '_PhaseScript.chr' + str(i+1) + '.pbs', 'w') as file:
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
                               'shapeit --input-bed ' + geno_bed_names[i] + ' ' + geno_bim_names[i] + ' '
                               + geno_fam_names[i] + ' --duohmm --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --output-max '
                               + output_names[i] + ' --output-log ' + output_names[i]
                               + '\n'
                               + 'shapeit -convert --input-haps ' + output_names[i] + '.haps '
                               + output_names[i] + '.sample --output-vcf ' + output_vcf_names[i] + ' --output-log '
                               + output_vcf_log[i])
            print("Done preparing chr" + str(i + 1) + " for phasing")

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
                ind_me_file = pd.read_csv(ind_me_names[i], sep='\t', header=0,
                                           dtype={0: str, 1: str, 2: float, 3: str, 4: float, 5: float, 6: int})
                father_me_errors = ind_me_file[ind_me_file['father_mendel'] > 0]
                mother_me_errors = ind_me_file[ind_me_file['mother_mendel'] > 0]
                if (len(father_me_errors) > 0) or (len(mother_me_errors) > 0):
                    print("Your files have people with non-zero mendel errors. You should investigate the "
                          + ind_me_names[i] + ' file and take a careful look at the people with high values in the '
                                              'father_mendel & mother_mendel column. This result suggests that your '
                                              'paternity/maternity assignment could be incorrect.')
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
                me_exclude = snp_me_file[snp_me_file['MendelError'] == 'Yes'][1]
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
                snp_hh_file = pd.read_csv(check_log_names[i] + '.snp.hh', sep='\t', header=None, skiprows=1,
                                          dtype={0: str, 1: int, 2: float, 3: float})
                # Create a new column with 'Yes' for snps with mendel error > 1%
                snp_hh_file['HHError'] = np.where((snp_hh_file[2] / snp_hh_file[3]) > 0.01, 'Yes', 'No')
                # Create new dataframe with only the positions with high HH rates
                snp_hh_exclude = snp_hh_file[snp_hh_file['HHError'] == 'Yes'][1]
                # If this is non-empty, make note:
                if len(snp_hh_exclude) > 0:
                    log.extend(['c'])
                # If it is empty (i.e. no HH errors exist) do nothing
                else:
                    pass
            else:
                pass

            # If a *.ind.hh file is produced, see if there are males with high haploid heterozygosity rate.
            if os.path.exists(check_log_names[i] + '.ind.hh'):
                # Read in file of people to check.
                ind_hh_file = pd.read_csv(check_log_names[i] + '.ind.hh', sep = '\t', header = None, skiprows=1,
                                          dtype={0: str, 1: int, 2: float, 3: float})
                # Create new column with info on whether the ind has mendel error > 1%
                ind_hh_file['HHError'] = np.where((ind_hh_file[2] / ind_hh_file[3]) > 0.01, 'Yes', 'No')
                # Create new dataframe of people with high mendel errors.
                ind_hh_exclude = ind_hh_file[ind_hh_file['HHError'] == 'Yes'][1]
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
                sys.exit("I shouldn't have gotten here. There is something wrong with the way the script is checking "
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
                               'shapeit --input-bed ' + geno_bed_names[i] + ' ' + geno_bim_names[i] + ' '
                               + geno_fam_names[i] + ' --chrX --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --exclude-snp '
                               + snp_exclude_name + ' --exclude-ind ' + check_log_names[i] + '.ind.hh.exclude'
                               + ' --output-max ' + output_names[i] + ' --output-log ' + output_names[i]
                               + '\n'
                               + 'shapeit -convert --input-haps ' + output_names[i] + '.haps '
                               + output_names[i] + '.sample --output-vcf ' + output_vcf_names[i] + ' --output-log '
                               + output_vcf_log[i])

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
                               'shapeit --input-bed ' + geno_bed_names[i] + ' ' + geno_bim_names[i] + ' '
                               + geno_fam_names[i] + ' --chrX --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --exclude-snp '
                               + snp_exclude_name + ' --output-max ' + output_names[i] + ' --output-log '
                               + output_names[i]
                               + '\n'
                               + 'shapeit -convert --input-haps ' + output_names[i] + '.haps '
                               + output_names[i] + '.sample --output-vcf ' + output_vcf_names[i] + ' --output-log '
                               + output_vcf_log[i])
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
                               'shapeit --input-bed ' + geno_bed_names[i] + ' ' + geno_bim_names[i] + ' '
                               + geno_fam_names[i] + ' --chrX --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --exclude-ind '
                               + check_log_names[i] + '.ind.hh.exclude' + ' --output-max ' + output_names[i]
                               + ' --output-log ' + output_names[i]
                               + '\n'
                               + 'shapeit -convert --input-haps ' + output_names[i] + '.haps '
                               + output_names[i] + '.sample --output-vcf ' + output_vcf_names[i] + ' --output-log '
                               + output_vcf_log[i])
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
                               'shapeit --input-bed ' + geno_bed_names[i] + ' ' + geno_bim_names[i] + ' '
                               + geno_fam_names[i] + ' --chrX --input-map '
                               + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                               + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i])
                               + ' ' + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --output-max '
                               + output_names[i] + ' --output-log ' + output_names[i]
                               + '\n'
                               + 'shapeit -convert --input-haps ' + output_names[i] + '.haps '
                               + output_names[i] + '.sample --output-vcf ' + output_vcf_names[i] + ' --output-log '
                               + output_vcf_log[i])

            print("Done preparing chr" + str(i+1) + " for phasing")

    # If the user is currently on the cluster, then submit the pbs files to start running.
    if on_cluster in ('y', 'yes'):
        for i in range(0, 23):
            # Submit to cluster
            subprocess.check_output(['qsub', 'Phasing/' + geno_name + '_PhaseScript.chr' + str(i+1) + '.pbs'])
    elif on_cluster in ('n', 'no'):
        sys.exit('Since you are not on the Penn State cluster right now, you should transfer the genotype files, 1000G '
                 'genetic maps, 1000G legend files, 1000G hap files, and 1000G sample file, shapeit, and the phasing '
                 'pbs files to the cluster. You can submit the files by using qsub name_of_file.pbs')
    else:
        sys.exit("You didn't answer 'yes' or 'no' when I asked whether you were on the cluster or not, so I'm assuming "
                 "you're not. You should transfer the genotype files, 1000G genetic maps, 1000G legend files, 1000G "
                 "hap files, and 1000G sample file, shapeit, and the phasing pbs files to the cluster. You can submit "
                 "the files by using qsub name_of_file.pbs")

    print("Done")


def impute():
    # Required by the sanger imputation server.
    # If not requesting pre-phasing, then all sites and samples should be phased with no missing data. - genophase
    # All alleles on the forward strand - genoharmonize
    # A single VCF file, not one file per-chromosome - this script
    # Coordinates are on GRCh37 - user, before genoharmonize
    # Chromosome names should be 1, 2, 3, etc. Not chr1, chr2, chr3, etc - this script
    # Records are sorted by genomic position (chromosomal order is not important) - this script
    # REF allele matches GRCh37. See the resources for help checking and fixing the REF allele - genoharmonize should do
    # this, but in this script we will double check.
    # Valid VCF - this script.

    import glob
    import urllib.request
    from subprocess import Popen, PIPE
    import shutil
    import re

    # Needed for sorting lists by numeric instead of string.
    def sorted_nicely(list):
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(list, key=alphanum_key)

    # Ask user what they want the eventual file to be named.
    print(Fore.BLUE + Style.BRIGHT)
    vcf_name = input("Please tell me what you would like your phased VCF file to be called. I'm going to concatenate "
                     "all of the phased per chromosomes into one whole genome file and give it this name.: ")
    print(Style.RESET_ALL)

    # Use glob to find files ending with the endings set in geno phase.
    raw_vcf_list = glob.glob('Phasing/*_PhasedTo1000G.chr*.vcf')

    # If glob wasn't able to find files ending in the above, ask the user where the phased files are.
    if raw_vcf_list == []:
        # Get path to phased files
        print(Fore.MAGENTA + Style.BRIGHT)
        vcf_path = input("Please tell me where your phased vcf and haps/sample files are. Make sure they are in their "
                         "own folder, as this part of the script looks for anything with the ending '.vcf' "
                         "(i.e. C:\\Users\\Julie White\\Box Sync\\Genotypes\\Phasing\\ etc.): ")
        print(Style.RESET_ALL)
        # Make list of anything with .vcf ending
        raw_vcf_list = glob.glob(os.path.join(vcf_path,'*.vcf'))
    else:
        vcf_path = 'Phasing'

    # Sort the vcf_list by numeric instead of string.
    vcf_list = sorted_nicely(raw_vcf_list)

    # If it doesn't exist, create SangerImputationFolder
    if not os.path.exists('SangerImputation'):
        os.makedirs('SangerImputation')

    # Determine if they have bcftools, vcftools, samtools, if not download them.
    if os.path.exists(os.path.join(bindir, 'bcftools')):
        pass
    else:
        import genodownload
        genodownload.bcftools()

    if os.path.exists(os.path.join(bindir, 'vcftools')):
        pass
    else:
        import genodownload
        genodownload.vcftools()

    if os.path.exists(os.path.join(bindir, 'samtools')):
        pass
    else:
        import genodownload
        genodownload.samtools()

    # hg19 fasta files.
    print(Fore.BLUE + Style.BRIGHT)
    fasta_exists = input('Have you already downloaded the human_g1k_v37.fasta.gz file? (y/n): ').lower()
    print(Style.RESET_ALL)

    if fasta_exists in ('y', 'yes'):
        # Ask the user where the fasta file is.
        print(Fore.MAGENTA + Style.BRIGHT)
        fasta_path = input('Please enter the pathname of where the your 1000G hg19 fasta file is '
                           '(i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\Fasta\\ etc.): ')
        print(Style.RESET_ALL)

    elif fasta_exists in ('n', 'no'):
        # Get geno download module
        import genodownload
        # Call download fasta command
        genodownload.fasta_1000G_hg19()
        # Saving fasta path
        fasta_path = os.path.join(os.getcwd(), '1000G_hg19_fasta')
    else:
        sys.exit('Please answer yes or no. Quitting now because no fasta file.')

    # See if the user already has htslib installed, if not, install it.
    p = Popen('bgzip', stdout = PIPE, stderr = PIPE)
    if p.returncode == None:
        pass
    if p.returncode != None:
        import genodownload
        genodownload.htslib()

    # Need to add a contig tag to all of the vcf files for bcftools to use them.
    # Must be zipped using bgzip, then indexed using samtools.
    awkstatement = str('{printf("##contig=<ID=%s,length=%d\\n",$1,$2);}')
    if os.path.exists(os.path.join(fasta_path, 'human_g1k_v37.fasta')):
        # zip using bgzip
        print("Zipping with bgzip")
        subprocess.check_output(['bgzip', os.path.join(fasta_path, 'human_g1k_v37.fasta')])
        # See if command will run on our now re-zipped file.
        p = Popen(['samtools', 'faidx', os.path.join(fasta_path, 'human_g1k_v37.fasta.gz')], stdout=PIPE,
                  stderr=PIPE)
        if p.returncode == None:
            subprocess.check_output(['samtools', 'faidx', os.path.join(fasta_path, 'human_g1k_v37.fasta.gz')])
            # Write the contig file.
            f = open(os.path.join(vcf_path, "Contigs.txt"), "w")
            subprocess.call(['awk', awkstatement, os.path.join(fasta_path, 'human_g1k_v37.fasta.gz.fai')], stdout=f)
            # Replace X with 23 in Contigs.txt file.
            sedstatement = str('s/ID=X/ID=23/g')
            subprocess.call(['sed', '-i', sedstatement, os.path.join(vcf_path, "Contigs.txt")])
            # Use bcftools to add contigs to each vcf.
            for vcf in vcf_list:
                subprocess.check_output(['bcftools', 'annotate', '-h', os.path.join(vcf_path, 'Contigs.txt'), vcf,
                                         '-Ov', '-o', vcf + '.contigs'])
        elif p.returncode != None:
            sys.exit(p.returncode)
    elif os.path.exists(os.path.join(fasta_path, 'human_g1k_v37.fasta.gz')):
        print("Unzipping with gunzip")
        # Unzip human_g1k_v37 file because samtools cannot index things that were indexed using gzip
        subprocess.check_output(['gunzip', os.path.join(fasta_path, 'human_g1k_v37.fasta.gz')])
        # Rezip using bgzip
        print("Rezipping with bgzip")
        subprocess.check_output(['bgzip', os.path.join(fasta_path, 'human_g1k_v37.fasta')])
        # Create fasta.fai file.
        print("Create fasta.fai file")
        subprocess.check_output(['samtools', 'faidx', os.path.join(fasta_path, 'human_g1k_v37.fasta.gz')])
        # Write the contig file.
        f = open(os.path.join(vcf_path, "Contigs.txt"), "w")
        subprocess.call(['awk', awkstatement, os.path.join(fasta_path, 'human_g1k_v37.fasta.gz.fai')], stdout=f)
        # Replace X with 23 in Contigs.txt file.
        sedstatement = str('s/ID=X/ID=23/g')
        subprocess.call(['sed', '-i', sedstatement, os.path.join(vcf_path, "Contigs.txt")])
        # Use bcftools to add contigs to each vcf.
        for vcf in vcf_list:
            subprocess.check_output(['bcftools', 'annotate', '-h', os.path.join(vcf_path, 'Contigs.txt'), vcf,
                                     '-Ov', '-o', vcf + '.contigs'])

    # Rename the .contigs vcf to be normal vcf names.
    for vcf in vcf_list:
        subprocess.call(['mv', vcf + '.contigs', vcf])

    # Use bcftools to concatenate the per chromosome phased vcf files. Save this in SangerImputation
    concat_vcf_list = str(" ".join(str(x) for x in vcf_list))
    subprocess.call('bcftools concat -Ov ' + concat_vcf_list + ' > SangerImputation/' + vcf_name + '.vcf', shell=True)

    # Bgzip this file
    subprocess.check_output('bgzip ' + os.path.join('SangerImputation', vcf_name + '.vcf'), shell=True)

    # Check VCF is sorted:
    subprocess.check_output(['bcftools', 'index', os.path.join('SangerImputation', vcf_name + '.vcf.gz')])

    # Save in phasing too.
    shutil.copy2('SangerImputation/' + vcf_name + '.vcf.gz', 'Phasing/' + vcf_name + '.vcf.gz')

    # Use bcftools and bgzip to zip the vcf files so that they stop taking up space.
    for file in vcf_list:
        subprocess.check_output('bgzip --index ' + file, shell=True)

    # Use bcftools to change the chromosome names, if needed.
    urllib.request.urlretrieve('https://imputation.sanger.ac.uk/www/plink2ensembl.txt',
                               'SangerImputation/plink2ensembl.txt')
    subprocess.check_output(['bcftools', 'annotate', '--rename-chrs', 'SangerImputation/plink2ensembl.txt', '-Oz', '-o',
                             os.path.join('SangerImputation', vcf_name + '_ChrUpdated.vcf.gz'),
                             os.path.join('SangerImputation', vcf_name + '.vcf.gz')])
    print("Done checking chromosome names")
    subprocess.call(['mv', os.path.join('SangerImputation', vcf_name + '_ChrUpdated.vcf.gz'),
                     os.path.join('SangerImputation', vcf_name + '.vcf.gz')])
    # Reindex VCF file cause we just changed it:
    subprocess.check_output(['bcftools', 'index', os.path.join('SangerImputation', vcf_name + '.vcf.gz')])

    # Check that our data matches the reference allele of 1000G Phase 3. Fix those that do not match.
    subprocess.check_output(['bcftools', 'norm', '--check-ref', 'ws', '--rm-dup', 'both', '--fasta-ref',
                             os.path.join(fasta_path, 'human_g1k_v37.fasta.gz'), '-Oz', '-o',
                             os.path.join('SangerImputation', vcf_name + '_RefChecked.vcf.gz'),
                             os.path.join('SangerImputation', vcf_name + '.vcf.gz')])
    print("Done checking if reference allele matches 1000G Phase3 and fixing those that don't match")
    subprocess.call(['mv', os.path.join('SangerImputation', vcf_name + '_RefChecked.vcf.gz'),
                     os.path.join('SangerImputation', vcf_name + '.vcf.gz')])
    # Reindex VCF file cause we just changed it:
    subprocess.check_output(['bcftools', 'index', os.path.join('SangerImputation', vcf_name + '.vcf.gz')])

    # Double check sex and fix ploidy
    subprocess.call('bcftools +guess-ploidy --genome b37 --tag GT '
                    + os.path.join('SangerImputation', vcf_name + '.vcf.gz') + ' > '
                    + os.path.join('SangerImputation', vcf_name + '_SexEst.txt'), shell=True)

    subprocess.check_output(['bcftools', '+fixploidy', '-Oz', '-o',
                             os.path.join('SangerImputation', vcf_name + '_PloidyChecked.vcf.gz'),
                             os.path.join('SangerImputation', vcf_name + '.vcf.gz'),
                             '--', '--sex', os.path.join('SangerImputation', vcf_name + '_SexEst.txt')])
    print("Done checking for sex and fixing ploidy")
    subprocess.call(['mv', os.path.join('SangerImputation', vcf_name + '_PloidyChecked.vcf.gz'),
                     os.path.join('SangerImputation', vcf_name + '.vcf.gz')])
    # Reindex VCF file cause we just changed it:
    subprocess.check_output(['bcftools', 'index', os.path.join('SangerImputation', vcf_name + '.vcf.gz')])

    # Check that you have a valid vcf
    print("Checking to see if you have a valid vcf now.")
    subprocess.check_output(['vcf-validator', os.path.join('SangerImputation', vcf_name + '.vcf.gz')])

    sys.exit("Your vcf has been checked and is ready for imputation on the Sanger Imputation Server. Please follow the "
             "directions here (https://imputation.sanger.ac.uk/?instructions=1) for uploading your data. The file is "
             "called " + vcf_name + ".vcf.gz")


# Extract the info quality scores from your imputed VCF files.
def getinfo(imputed_path, geno_name):
    # Determine if they have vcftools, if not download.
    if os.path.exists(os.path.join(bindir, 'vcftools')):
        pass
    else:
        import genodownload
        genodownload.vcftools()

    # Ask the user if they are running this from the PSU cluster. This is the preferred way to go, since this takes a
    # while and will clog their computer's RAM if they try to do it from their personal computer.
    print(Fore.CYAN)
    on_cluster = input('Are you currently running this from the Penn State ACI-B cluster? If yes, I can submit the jobs'
                       ' for you. If not, I will try to run this using your local RAM, but this might make your '
                       'computer very slow and will probably take several hours. (y/n): ').lower()
    print(Style.RESET_ALL)

    if on_cluster in ("yes", "y"):
        # I based this formatting off of PSU cluster users, so they need to have a PSU cluster allocation.
        print(Fore.BLUE + Style.BRIGHT)
        allocation_name = input('Please enter the name of your cluster allocation: ')
        print(Style.RESET_ALL)

        # Write pbs file
        with open(os.path.join(imputed_path, 'GetImputationInfo.pbs'), 'w') as file:
            file.write('#!/bin/bash\n'
                       '#PBS -l walltime=24:00:00\n'
                       '#PBS -l nodes=1:ppn=8\n'
                       '#PBS -l pmem=7gb\n'
                       '#PBS -A ' + allocation_name +
                       '\n'
                       '#PBS -j oe\n'
                       'cd $PBS_O_WORKDIR\n'
                       '\n' +
                       'for i in {1..23}; do vcftools --gzvcf ' + geno_name
                       + '_chr$i.vcf.gz --get-INFO INFO --get-INFO RefPanelAF --out ' + geno_name
                       + '_chr$[i]_InfoScoreAF; done')
        # Change directory into imputed path
        os.chdir(imputed_path)
        # Submit to cluster
        subprocess.check_output(['qsub', 'GetImputationInfo.pbs'])
    # If they aren't on the cluster, then try to run this from their computer.
    elif on_cluster in ('no', 'n'):
        subprocess.check_output('for i in {1..23}; do vcftools --gzvcf '
                                + os.path.join(imputed_path, geno_name + '_chr$i.vcf.gz')
                                + ' --get-INFO INFO --get-INFO RefPanelAF --out '
                                + os.path.join(imputed_path, geno_name + '_chr$[i]_InfoScoreAF') + '; done', shell=True)
        print("Done. Your files have the ending .INFO and are in " + imputed_path)

    else:
        sys.exit("You didn't answer 'yes' or 'no' when I asked whether you were on the cluster or not, so I don't know "
                 "how to proceed. Exiting now.")


# Make plots of imputation quality score.
def qualscoreplot(info_path):
    import glob
    import re

    try:
        import pandas as pd
    except (ImportError, ModuleNotFoundError):
        import genodownload
        genodownload.getpandas()
        import pandas as pd

    try:
        import numpy as np
    except (ImportError, ModuleNotFoundError):
        import genodownload
        genodownload.getnumpy()
        import numpy as np

    try:
        import matplotlib
        matplotlib.use('Agg')
    except ImportError:
        import genodownload
        genodownload.getmatplotlib()
        import matplotlib
        matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    # Needed for sorting lists by numeric instead of string.
    def sorted_nicely(list):
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(list, key=alphanum_key)

    # Use glob to find files ending with the ending .INFO
    info_list = glob.glob(os.path.join(info_path, '*.INFO'))

    # If list is empty, exit.
    if not info_list:
        sys.exit("I can't find any .INFO files in the path that you specified. Please double check that they are there "
                 "and have the ending .INFO and try again.")

    sorted_info_list = sorted_nicely(info_list)

    # First need to combine all of the .INFO files.
    pandas_info = ['INFO_%d' % x for x in range(1, 24)]
    for i in range(0, len(sorted_info_list)):
        if i < 22:
            pandas_info[i] = pd.read_csv(sorted_info_list[i], sep='\t', header=0,
                                         dtype={'CHROM': int, 'POS': int, 'REF': str, 'ALT': str, 'INFO': float,
                                                'RefPanelAF': float})
        elif i == 22:
            pandas_info[i] = pd.read_csv(sorted_info_list[i], sep='\t', header=0,
                                         dtype={'CHROM': str, 'POS': int, 'REF': str, 'ALT': str, 'INFO': float,
                                                'RefPanelAF': float})
            pandas_info[i].iloc[:, 0].replace('X', 23, inplace=True)
            #pandas_info[i]['CHROM'] = pandas_info[i]['CHROM'].astype(int)

    # Concatenate all of the id updates into one file.
    combined_info = pd.concat([pandas_info[0], pandas_info[1], pandas_info[2], pandas_info[3], pandas_info[4],
                               pandas_info[5], pandas_info[6], pandas_info[7], pandas_info[8], pandas_info[9],
                               pandas_info[10], pandas_info[11], pandas_info[12], pandas_info[13], pandas_info[14],
                               pandas_info[15], pandas_info[16], pandas_info[17], pandas_info[18], pandas_info[19],
                               pandas_info[20], pandas_info[21], pandas_info[22]])

    INFO = np.asarray(combined_info)[:,4]

    f = plt.figure()
    hist, bins = np.histogram(INFO, bins = 100)
    width = np.diff(bins)
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, width=width, color="green", edgecolor="black")
    plt.xlabel('Imputation Quality Score')
    plt.ylabel('Number of SNPs')
    plt.title('Histogram of Genome-Wide Imputation Quality Scores', y=1.08)
    plt.show()
    f.savefig("ImputationQualScores_Hist.pdf")
