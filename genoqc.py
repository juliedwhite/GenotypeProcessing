import platform
import os
from os.path import expanduser
home = expanduser("~")
bindir = os.path.join(home, 'software', 'bin')

# Since we use plink a lot, I'm going to go ahead and set a plink variable with the system-specific plink name.
system_check = platform.system()
if system_check in ("Linux", "Darwin"):
    plink = "plink"
    rm = "rm "
elif system_check == "Windows":
    plink = 'plink.exe'
    rm = "del "

# Determine if they have plink, if not download it.
if os.path.exists(os.path.join(bindir, plink)):
    pass
else:
    import genodownload
    genodownload.plink()


def estimate_sex(geno_name):
    import subprocess

    print("This can take a while, so sit back for a bit and don't worry.")

    # Checks the sex listed in the file against the chromosomal sex.
    subprocess.check_output([plink, '--bfile', geno_name, '--split-x', 'hg19', 'no-fail', '--out', geno_name])
    subprocess.check_output([plink, '--bfile', geno_name, '--indep-pairphase', '20000', '2000', '0.5', '--out',
                             geno_name])
    subprocess.check_output([plink, '--bfile', geno_name, '--exclude', geno_name + '.prune.out', '--check-sex',
                             'ycount', '0.3', '0.8', '0', '0', '--out', geno_name])

    print("Finished. Your sex estimation files will have the name " + geno_name
          + ".sexcheck. You should check this file for any problems and decide if you want to update your file.")


def update_sex(geno_name, update_sex_filename):
    # File for updating sex should have:
    #   1) FID
    #   2) IID
    #   3) Sex (1 = M, 2 = F, 0 = missing)
    import subprocess

    subprocess.check_output([plink, '--bfile', geno_name, '--update-sex', update_sex_filename, '--make-bed', '--out',
                             geno_name + '_SexUpdated'])

    print("Finished. Your genotype files with sex updated will have the name " + geno_name + "_SexUpdated.")


def missing_call_rate(geno_name):
    # Exclude SNPs (geno) and people (mind) with missing call rates > 10%
    import subprocess
    subprocess.check_output([plink, '--bfile', geno_name, '--geno', '0.1', '--make-bed', '--out',
                             geno_name + '_geno0.1'])
    subprocess.check_output([plink, '--bfile', geno_name + '_geno0.1', '--mind', '0.1', '--make-bed', '--out',
                             geno_name + '_geno0.1_mind0.1'])
    subprocess.call(rm + geno_name + '_geno0.1.*', shell=True)

    print("Finished. Your pruned genotype files will have the name " + geno_name + "_geno0.01_mind0.01")


def het(geno_name):
    # Identifies individuals with extreme heterozygosity values (more than +- 3 SD)
    # Getting extra required modules
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

    import subprocess

    # Use plink to calculate the heterozygosity, paying attention to geno and mind.
    subprocess.check_output([plink, '--bfile', geno_name, '--geno', '0.1', '--mind', '0.1', '--het', '--out',
                             geno_name])

    # Read het file into pandas
    het_file = pd.read_csv(geno_name + '.het', sep='\s+', header=0)

    # Create new column with formula: (N(NM)-O(HOM))/N(NM)
    het_file['HET'] = (het_file['N(NM)'] - het_file['O(HOM)']) / het_file['N(NM)']

    # Getting standard deviation and average of HET column
    het_sd = np.std(het_file['HET'])
    het_avg = np.mean(het_file['HET'])

    # Add label 'keep' to people within 3*SD of the average het value, give 'remove' to everyone else.
    het_file['HET_Filter'] = np.where(
        (het_file['HET'] > het_avg - 3 * het_sd) & (het_file['HET'] < het_avg + 3 * het_sd), 'Keep', 'Remove')
    # Write this file so the user has it.
    het_file.to_csv(geno_name + '.het', sep='\t', header=True, index=False)
    # Make a list of the people who pass the filter.
    het_keep = het_file[het_file['HET_Filter'] == 'Keep']
    # Write this file so that we can use it later to filter people out.
    het_keep[['FID', 'IID']].to_csv(geno_name + '_KeptAfterHetCheck.txt', sep='\t', header=False, index=False)
    # Make a list of the people who fail the filter.
    het_rem = het_file[het_file['HET_Filter'] == 'Remove']
    # Write this to file so we have the record.
    het_rem.to_csv(geno_name + '_RemAfterHetCheck.txt', sep='\t', header=True, index=False)

    # Make new plink file with people passing het check.
    subprocess.check_output([plink, '--bfile', geno_name, '--keep', geno_name + '_KeptAfterHetCheck.txt', '--geno',
                             '0.1', '--make-bed', '--out', geno_name + '_HetChecked'])

    print("Done. Your new file of people with non-extreme heterozygosity values will be called "
          + geno_name + "_HetChecked")
