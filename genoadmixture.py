from colorama import init, Fore, Style
init()

import platform
# Since we use plink a lot, I'm going to go ahead and set a plink variable with the system-specific plink name.
system_check = platform.system()
if system_check in ("Linux", "Darwin"):
    plink = "./plink"
elif system_check == "Windows":
    plink = 'plink.exe'


def prep(admix_name):
    import os
    import subprocess
    
    # I based this formatting off of PSU cluster users, so they need to have a PSU cluster allocation.
    print(Fore.BLUE + Style.BRIGHT)
    allocation_name = input('Please enter the name of your cluster allocation: ')
    print(Style.RESET_ALL)

    # Check if folder called 'Admixture' exists, if not, create it.
    if not os.path.exists('Admixture'):
        os.makedirs('Admixture')

    # Ask if they have relatives in their sample.
    print(Fore.MAGENTA + Style.BRIGHT)
    relative_check = input('Do you have relatives in your sample? Perhaps those identified in an IBD analysis '
                           '(y/n): ').lower()
    print(Style.RESET_ALL)

    # We want the LD correction to be the same for all sets, so do this on the full genotype file and put it in the
    # Admixture file.
    subprocess.check_output([plink, '--bfile', admix_name, '--indep-pairwise', '50', '10', '2', '--out',
                             'Admixture/' + admix_name])

    if relative_check in ('y', 'yes'):
        # If they have relatives in their sample, get a list of the filenames for each set of people.
        print(Fore.GREEN)
        user_sets = input('Please give me a comma separated list of your set list filenames (with file extenstion). '
                          'I.e. dataset_setA.txt, dataset_setB.txt, etc. To do this, break up your entire dataset (not '
                          'just related individuals) across sets, making sure that there are not related individuals '
                          'within each set. These lists should be space or tab delimited with FID then IID: ')
        print(Style.RESET_ALL)

        # Convert the user given list to a python list.
        set_list = user_sets.split(', ')
        print(set_list)

        # Perform the admixture prep separately on each set.
        for i in range(0,len(set_list)):
            # Tell the user what they gave as file names and what I'm going to output as file names (SetA, SetB, SetC,
            # etc.)
            set_name = chr(ord('a') + i).upper()
            print(set_list[i] + ' = Set' + set_name)

            # For each set, extract those people from the working genotype file and remove SNPs in LD. These files are
            # what the user should put on the cluster.
            subprocess.check_output([plink, '--bfile', admix_name, '--keep', set_list[i], '--extract',
                                     'Admixture/' + admix_name + '.prune.in', '--make-bed', '--out',
                                     'Admixture/' + admix_name + '_Set' + set_name + '_LDPruned'])

            # For each set, write a pbs script for admixture k = 3..6
            with open('Admixture/' + admix_name + '_Set' + set_name + '_Admixture_k3to6.pbs', 'w') as file:
                file.write('#!/bin/bash\n'
                           '#PBS -l walltime=150:00:00\n'
                           '#PBS -l nodes=1:ppn=8\n'
                           '#PBS -l pmem=8gb\n'
                           '#PBS -A ' + allocation_name + '\n'
                           '#PBS -j oe\n'
                           'cd $PBS_O_WORKDIR\n'
                           '\n'
                           'for K in {3..6}; do ./admixture --cv '
                           + admix_name + '_Set' + set_name + '_LDPruned.bed $K | tee '
                           + admix_name + '_Set' + set_name + '_LDPruned.log${K}.out; done')

            # For each set, write a pbs script for admixture k = 7..9
            with open('Admixture/' + admix_name + '_Set' + set_name + '_Admixture_k7to9.pbs', 'w') as file:
                file.write('#!/bin/bash\n'
                           '#PBS -l walltime=150:00:00\n'
                           '#PBS -l nodes=1:ppn=8\n'
                           '#PBS -l pmem=8gb\n'
                           '#PBS -A ' + allocation_name + '\n'
                           '#PBS -j oe\n'
                           'cd $PBS_O_WORKDIR\n'
                           '\n'
                           'for K in {7..9}; do ./admixture --cv '
                           + admix_name + '_Set' + set_name + '_LDPruned.bed $K | tee '
                           + admix_name + '_Set' + set_name + '_LDPruned.log${K}.out; done')

        # Tell the user it's finished and give them directions.
        print("For each set, transfer " + admix_name + "_LDPruned bed/bim/fam files, " + admix_name
              + "_Admixture_k3to6.pbs, and " + admix_name + "_Admixture_7to9.pbs files to the cluster.\n"
              "For each set, submit them using qsub " + admix_name + "Admixture_k3to6.pbs and qsub " + admix_name
              + "Admixture_k7to9.pbs\n"
              "When you get your results, you should evaluate them to see which makes sense given your study "
              "population and which has the lowest CV value.")

    # If there's no relatives, we can do all of this on the full genotype file.
    elif relative_check in ('n', 'no'):
        # For all people, create file that is LD pruned.
        subprocess.check_output([plink, '--bfile', admix_name, '--extract',
                                 'Admixture/' + admix_name + '.prune.in', '--make-bed', '--out',
                                 'Admixture/' + admix_name + '_LDPruned'])

        # For all people, create a pbs file for admixture k = 3..6
        with open ('Admixture/' + admix_name + '_Admixture_k3to6.pbs', 'w') as file:
            file.write('#PBS -l walltime=150:00:00\n'
                       '#PBS -l nodes=1:ppn=8\n'
                       '#PBS -l pmem=8gb\n'
                       '#PBS -A ' + allocation_name + '\n'
                       '#PBS -j oe\n'
                       'cd $PBS_O_WORKDIR\n'
                       '\n'
                       'for K in {3..6}; do ./admixture --cv ' + admix_name + '_LDPruned.bed $K | tee '
                       + admix_name + '_LDPruned.log${K}.out; done')

        # For all people, create a pbs file for admixture k = 7..9
        with open('Admixture/' + admix_name + '_Admixture_k7to9.pbs', 'w') as file:
            file.write('#PBS -l walltime=150:00:00\n'
                       '#PBS -l nodes=1:ppn=8\n'
                       '#PBS -l pmem=8gb\n'
                       '#PBS -A ' + allocation_name + '\n'
                       '#PBS -j oe\n'
                       'cd $PBS_O_WORKDIR\n'
                       '\n'
                       'for K in {7..9}; do ./admixture --cv ' + admix_name + '_LDPruned.bed $K | tee '
                       + admix_name + '_LDPruned.log${K}.out; done')

        # End and give directions.
        print("Transfer " + admix_name + "_LDPruned bed/bim/fam files, " + admix_name + "_Admixture_k3to6.pbs, and "
              + admix_name + "_Admixture_7to9.pbs files to the cluster.\n"
              "Submit them using qsub " + admix_name + "Admixture_k3to6.pbs and qsub " + admix_name
              + "Admixture_k7to9.pbs\n"
                "When you get your results, you should evaluate them to see which makes sense given your study "
                "population and which has the lowest CV value.")
