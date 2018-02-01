import platform
import os
from os.path import expanduser

try:
    import colorama
    from colorama import init, Fore, Style
    init()
except ImportError:
    import genodownload
    genodownload.getcolorama()
    import colorama
    from colorama import init, Fore, Style
    init()

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


def prep(admix_name):
    import os
    import subprocess
    import sys
    import shutil
    import pandas as pd

    # Check if folder called 'Admixture' exists, if not, create it.
    if not os.path.exists('Admixture'):
        os.makedirs('Admixture')

    # Ask if the user is on the cluster right now to determine if we should submit the files for them
    print(Fore.BLUE + Style.BRIGHT)
    on_cluster = input('Are you currently running this from the Penn State ACI-B cluster? If yes, I make this '
                       'process as a job and submit it to run. If you are not on the cluster, then you will have to '
                       'submit the jobs yourself. (y/n): ').lower()
    print(Style.RESET_ALL)

    # If they are on the cluster, then see if they have already downloaded the admixture program.
    if on_cluster in ("yes", "y"):
        # Determine if they have admixture, if not download it.
        if os.path.exists(os.path.join(bindir, 'admixture')):
            pass
        else:
            import genodownload
            genodownload.admixture()
    # If they are not on the cluster. Tell them to get it before they submit the jobs.
    elif on_cluster in ("no", "n"):
        print(Fore.RED + Style.BRIGHT + "Make sure to download the admixture program on the cluster and have it in the "
                                        "same directory when you submit your jobs.")
        print(Style.RESET_ALL)

    # I based this formatting off of PSU cluster users, so they need to have a PSU cluster allocation.
    print(Fore.BLUE + Style.BRIGHT)
    allocation_name = input('Please enter the name of your cluster allocation: ')

    if allocation_name == "open":
        walltime = "48:00:00"
        print("Since you are on the open queue, you have a walltime limit of 48hrs. I'm not sure if the admixture will"
              "finish by then, so after the job is done or gets killed you should check your log files to see which "
              "finished. If any didn't finish, then modify the .pbs files to change the 'for K in {n..n}' part to "
              "reflect the K values that did not complete. Then use qsub filename.pbs to resubmit the jobs.")
    else:
        walltime = "150:00:00"
    print(Style.RESET_ALL)

    # Ask the user what k-values they would like to run
    print(Fore.GREEN + "Admixture will partition your genetic variation into K different groups. As K is user-defined, "
                       "I'm going to ask you what values of K you would like me to run. Usually I start with K=2 and go"
                       " to K=12.")
    k_start = input("What K value would you like to start at (i.e. 2)?: ")
    k_end = input("What K value would you like to end at (i.e. 12)?: ")

    if k_start.isdigit():
        pass
    else:
        sys.exit("Please enter an integer for the starting K value. Exiting now.")

    if k_end.isdigit():
        pass
    else:
        sys.exit("Please enter an integer for the ending K value. Exiting now.")

    k_values = list(range(int(k_start), int(k_end)+1))

    # Move files for running admixture to Admixture folder
    shutil.copy2(admix_name + '.bed', 'Admixture')
    shutil.copy2(admix_name + '.bim', 'Admixture')
    shutil.copy2(admix_name + '.fam', 'Admixture')

    # Move to Admixture directory
    os.chdir('Admixture')

    # We want the LD correction to be the same for all sets, so do this on the full genotype file and put it in the
    # Admixture file.
    subprocess.check_output([plink, '--bfile', admix_name, '--indep', '50', '10', '2', '--out', admix_name])
    subprocess.check_output([plink, '--bfile', admix_name, '--extract', admix_name + '.prune.in', '--make-bed',
                             '--out', admix_name + '_LDPruned'])

    # Create pbs files for admixture
    for i in range(0, len(k_values)):
        with open(admix_name + '_Admixture_k' + str(k_values[i]) + '.pbs', 'w') as file:
            file.write('#!/bin/bash\n'
                       '#PBS -l walltime=' + walltime + '\n'
                       '#PBS -l nodes=1:ppn=8\n'
                       '#PBS -l pmem=8gb\n'
                       '#PBS -A ' + allocation_name + '\n'
                       '#PBS -j oe\n'
                       'cd $PBS_O_WORKDIR\n'
                       '\n'
                       'admixture -j5 --cv ' + admix_name + '_LDPruned.bed ' + str(k_values[i]) + ' | tee '
                       + admix_name + '_LDPruned.log' + str(k_values[i]) + '.out')

    if on_cluster in ("yes", "y"):
        # Submit the jobs
        for i in range(0, len(k_values)):
            os.system('qsub ' + admix_name + '_Admixture_k' + str(k_values[i]) + '.pbs')
        print("Your admixture jobs have been submitted. You can check their status using qstat -u usrname. When "
              "you get your results, you should evaluate them to see which makes sense given your study population"
              "and which has the lowest CV value (located in the log)")

    if on_cluster in ("no", "n"):
        print("Transfer " + admix_name + "_LDPruned bed/bim/fam files and all of the .pbs files to the cluster. "
                                         "Also make sure you have the admixture program in the same location as "
                                         "these files or in your PATH.\n"
              + "Submit them using qsub " + admix_name + "Admixture_k{N}.pbs\n"
              + "When you get your results, you should evaluate them to see which makes sense given your study "
                "population and which has the lowest CV value (located in the log).")
