import platform

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

# Since we use plink a lot, I'm going to go ahead and set a plink variable with the system-specific plink name.
system_check = platform.system()
if system_check in ("Linux", "Darwin"):
    plink = "./plink"
elif system_check == "Windows":
    plink = 'plink.exe'


def prep(admix_name):
    import os
    import subprocess
    import sys
    import shutil

    # Check if folder called 'Admixture' exists, if not, create it.
    if not os.path.exists('Admixture'):
        os.makedirs('Admixture')

    # Ask if the user is on the cluster right now to determine if we should submit the files for them
    print(Fore.BLUE + Style.BRIGHT)
    on_cluster = input('Are you currently running this from the Penn State ACI-B cluster? If yes, I make this '
                       'process as a job and submit it to run. If you are not on the cluster, then you will have to '
                       'submit the jobs yourself. (y/n): ').lower()
    print(Style.RESET_ALL)

    # If they are on the cluster, then ask if they have already downloaded the admixture program.
    if on_cluster in ("yes", "y"):
        print(Fore.GREEN)
        admixture_exists = input('Have you already downloaded the admixture program? (y/n): ').lower()
        print(Style.RESET_ALL)
        # If yes, then get path to admixture
        if admixture_exists in ('y', 'yes'):
            print(Fore.CYAN)
            admixture_path = input('Please enter the pathname of where the admixture program is '
                                   '(i.e. /storage/home/jdw345/software/Admixture_1.3.0_Linux/): ')
            print(Style.RESET_ALL)
        # If no, then download admixture
        elif admixture_exists in ('n', 'no'):
            import genodownload
            genodownload.admixture()
            admixture_path = os.path.join(os.getcwd(), 'Admixture_1.3.0_Linux/')
        else:
            sys.exit('You did not give a yes or no answer when I asked if you had admixture. Quitting now.')
        # Copy admixture program to the admixture folder.
        shutil.copy2(os.path.join(admixture_path, 'admixture'), 'Admixture')
    # If they are not on the cluster. Tell them to get it before they submit the jobs.
    elif on_cluster in ("no", "n"):
        print(Fore.RED + Style.BRIGHT + "Make sure to download the admixture program on the cluster and have it in the "
                                        "same directory when you submit your jobs.")
        print(Style.RESET_ALL)

    # I based this formatting off of PSU cluster users, so they need to have a PSU cluster allocation.
    print(Fore.BLUE + Style.BRIGHT)
    allocation_name = input('Please enter the name of your cluster allocation: ')

    if allocation_name == "open":
        walltime = "24:00:00"
        print("Since you are on the open queue, you have a walltime limit of 24hrs. I'm not sure if the admixture will"
              "finish by then, so after the job is done or gets killed you should check your log files to see which "
              "finished. If any didn't finish, then modify the .pbs files to change the 'for K in {n..n}' part to "
              "reflect the K values that did not complete. Then use qsub filename.pbs to resubmit the jobs.")
    else:
        walltime = "150:00:00"
    print(Style.RESET_ALL)

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
        user_sets = input('Please give me a comma separated list of your set list filenames (with file extension). '
                          'I.e. dataset_setA.txt, dataset_setB.txt, etc. To do this, break up your entire dataset (not '
                          'just related individuals) across sets, making sure that there are not related individuals '
                          'within each set. These lists should be space or tab delimited with FID then IID. Please '
                          'enter the comma separated list here: ')
        print(Style.RESET_ALL)

        # Convert the user given list to a python list.
        set_list = user_sets.split(', ')
        print(set_list)

        # Perform the admixture prep separately on each set.
        for i in range(0, len(set_list)):
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
                           '#PBS -l walltime=' + walltime + '\n'
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
                           '#PBS -l walltime=' + walltime + '\n'
                           '#PBS -l nodes=1:ppn=8\n'
                           '#PBS -l pmem=8gb\n'
                           '#PBS -A ' + allocation_name + '\n'
                           '#PBS -j oe\n'
                           'cd $PBS_O_WORKDIR\n'
                           '\n'
                           'for K in {7..9}; do ./admixture --cv '
                           + admix_name + '_Set' + set_name + '_LDPruned.bed $K | tee '
                           + admix_name + '_Set' + set_name + '_LDPruned.log${K}.out; done')

        if on_cluster in ("yes", "y"):
            os.chdir('Admixture')
            # Submit the jobs
            subprocess.call(['qsub', '*.pbs'], shell=True)
            print("Your admixture jobs have been submitted. You can check their status using qstat -u usrname. When "
                  "you get your results, you should evaluate them to see which makes sense given your study population"
                  "and which has the lowest CV value (located in the log)")

        if on_cluster in ("no", "n"):
            # Tell the user it's finished and give them directions.
            print("For each set, transfer the " + admix_name + "_LDPruned bed/bim/fam files, " + admix_name
                  + "_Admixture_k3to6.pbs, and " + admix_name + "_Admixture_7to9.pbs files to the cluster Also make "
                                                                "sure you have the admixture program in the same "
                                                                "folder as these files on the cluster.\n"
                  "For each set, submit them using qsub " + admix_name + "Admixture_k3to6.pbs and qsub " + admix_name
                  + "Admixture_k7to9.pbs\n"
                  "When you get your results, you should evaluate them to see which makes sense given your study "
                  "population and which has the lowest CV value (located in the log).")

    # If there's no relatives, we can do all of this on the full genotype file.
    elif relative_check in ('n', 'no'):
        # For all people, create file that is LD pruned.
        subprocess.check_output([plink, '--bfile', admix_name, '--extract',
                                 'Admixture/' + admix_name + '.prune.in', '--make-bed', '--out',
                                 'Admixture/' + admix_name + '_LDPruned'])

        # For all people, create a pbs file for admixture k = 3..6
        with open('Admixture/' + admix_name + '_Admixture_k3to6.pbs', 'w') as file:
            file.write('#!/bin/bash\n'
                       '#PBS -l walltime=' + walltime + '\n'
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
            file.write('#!/bin/bash\n'
                       '#PBS -l walltime=' + walltime + '\n'
                       '#PBS -l nodes=1:ppn=8\n'
                       '#PBS -l pmem=8gb\n'
                       '#PBS -A ' + allocation_name + '\n'
                       '#PBS -j oe\n'
                       'cd $PBS_O_WORKDIR\n'
                       '\n'
                       'for K in {7..9}; do ./admixture --cv ' + admix_name + '_LDPruned.bed $K | tee '
                       + admix_name + '_LDPruned.log${K}.out; done')

        if on_cluster in ("yes", "y"):
            os.chdir('Admixture')
            # Submit the jobs
            subprocess.call(['qsub', admix_name + '_Admixture_k3to6.pbs'], shell=True)
            subprocess.call(['qsub', admix_name + '_Admixture_k7to9.pbs'], shell=True)
            print("Your admixture jobs have been submitted. You can check their status using qstat -u usrname. When "
                  "you get your results, you should evaluate them to see which makes sense given your study population"
                  "and which has the lowest CV value (located in the log)")

        if on_cluster in ("no", "n"):
            print("Transfer " + admix_name + "_LDPruned bed/bim/fam files, " + admix_name + "_Admixture_k3to6.pbs, and "
                  + admix_name + "_Admixture_7to9.pbs files to the cluster. Also make sure you have the admixture "
                                 "program in the same location as these files.\n"
                  "Submit them using qsub " + admix_name + "Admixture_k3to6.pbs and qsub " + admix_name
                  + "Admixture_k7to9.pbs\n"
                    "When you get your results, you should evaluate them to see which makes sense given your study "
                    "population and which has the lowest CV value (located in the log).")
