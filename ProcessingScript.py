# Julie's Processing script for Shriver Lab genotype data
# Date 9.5.17

# WHAT THIS DOES #
# 1) Download Plink
# 2) Update sex
# 3) Produce new dataset with people and SNPs with missing call rate < 10%
# 4) Run IBD matrix through plink
# 5) Update family (FID) and individual IDs (IID)
# 6) Update maternal and paternal IDs
# 7) Prepare for ADMIXTURE with k = 3...9. Relatedness matters


# 7) Merges your data with 1000 Genomes

# 9) Prepares files for phasing using SHAPEIT2. Relatedness matters
# 10) Prepares files for imputation using the Sanger Imputation Server. Relatedness matters.

# REQUIREMENTS #
# You must have plink 1.9 (https://www.cog-genomics.org/plink2) in the same directory as your genotype files and this
#   script, if you do not have plink already, run step #1 and we will download it.
# You must have downloaded the 1000G Phase 3 legend files from http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/.
#   These can be anywhere you want, you'll tell me the path later
# You must have downloaded the 1000G Phase 3 VCF files from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/.
#   These can be anywhere you want, you'll tell me the path later

#GenoQC
# Update sex
# Missing call rate

#Relatives
# Run IBD to identify relatives
# Update FID IID information
# Update parental IDs

#Admixture Steps:
# Missing call rate
# Run IBD to identify relatives
# Update FID IID information
# Update parental IDs
# Harmonize with 1000G Phase 3
# Merge with 1000G
# Prepare for ADMIXTURE with k = 3..9
# If on PSU cluster, can submit.

#Imputation steps:
# Must be on hg19, user should check this.
# Missing call rate
# Heterozygosity check (MAF threshold, HWE p-value)
# Set haploid genotypes (male chr X) as missing
# Check for gender mismatches
# Harmonize with 1000G Phase3
# Prephasing check
# Phasing file preparation
# If on PSU cluster, can submit phasing to shapeit.
# Imputation file preparation.
#   -Convert to single VCF, not one file per chromosome
#   -Check validity of VCF
#   -Records sorted by genomic position
#   -Chromosome names should be 1, 2, 3, etc… not chr1, chr2, chr3, etc… Some programs will represent X as 23, Y as 24, etc…. Please remove or replace these names.
#   -Give sample file (samples.txt) with Sample Name then M or F.
# User submits to Sanger Imputation Server (https://imputation.sanger.ac.uk/?instructions=1)

# Getting the needed modules.
import os
import shutil
import glob
import sys

to_do = input('\u001b[31;1m What would you like to do?\n'
              '1) Download Plink\n'
              '2) Download 1000G Phase 3 VCF files\n'
              '3) Download 1000G Phase 3 Hap/Legend/Sample files \n'
              '4) Download Genotype Harmonizer\n'
              '5) Update sex. You need a file with FID, IID, Sex (M=1, F=2, Unknown=0) (in that order, no column headings)\n'
              '6) Produce a new dataset of people and SNPs with missing call rate < 10%\n'
              '7) Run an IBD analysis to identify relatives. All you need are plink bed/bim/fam files.\n'
              '8) Update FID or IID information. You need a file with the following information Old FID, Old IID, '
              'New FID, New IID.\n'
              '9) Update parental IDs. You need a file with FID, IID, Paternal IID, and Maternal IID.\n'
              '10) Prepare for ADMIXTURE with 1000G Phase 3 files'
              ') Nothing. \n'
              'Please enter a number (i.e. 2): \u001b[0m')

#GenoDownload: Download Plink
if to_do == '1':
    #Get the module for downloading stuff.
    import genodownload
    #Call download plink command
    genodownload.plink()

#GenoDownload: Download 1000G VCF files.
elif to_do == '2':
    #Get the module for downloading stuff.
    import genodownload
    #Call the download 1000G phase 3 VCF command.
    genodownload.vcf_1000g_phase3()

#GenoDownload: Download 1000G HapLegendSample files.
elif to_do == '3':
    #Get the module for downloading stuff.
    import genodownload
    #Call the download 1000G Phase 3 HapLegendSample command
    genodownload.hls_1000g_phase3()

#GenoDownload: Download Genotype Harmonizer.
elif to_do == '4':
    #Get the module
    import genodownload
    #Call the download Genotype Harmonizer command
    genodownload.genotype_harmonizer()

#GenoQC: Update sex
elif to_do == '5':
    #Get name of genotype file
    geno_name = input("\u001b[32;1m Please enter the name of the plink genotype files you'd like to update sex in "
                      "(without bed/bim/fam extension: \u001b[0m")

    #Get name of file to be used for updating sex
    update_sex_filename = input('\u001b[34;1m Please enter the name of your text file for updating sex (with file extension): \u001b[0m')

    #Import module where this command is.
    import genoqc

    #Call UpdateSex command using geno name and update sex filename as input
    genoqc.update_sex(geno_name, update_parents_filename)

#GenoQC: Clean dataset by missing call rate > 10%
elif to_do == '6':
    #Get name of genotype file
    geno_name = input('\u001b[32;1m Please enter the name of the genotype files (without bed/bim/fam extension: \u001b[0m')

    #Import module and call command
    import genoqc
    genoqc.missing_call_rate(geno_name)

#GenoRelatives: Run IBD
elif to_do == '7':
    # Identity-by-descent in Plink
    # This part of the script will prune for LD, calculate IBD, and exclude individuals who have IBD < 0.2
    # The IBD results will have .genome appended to your file name. I have also included a line to convert the IBD results
    #   from whitespace to tab delimited. This will have .tab.genome appended to your filename.

    # Important values of Pi-hat
    #   -First-degree relative = 0.5 (full sibs, parent-offspring)
    #   -Second-degree relative = 0.25 (half-sibs, uncle/aunt-nephew/niece, grandparent-grandchild)
    #   -Third-degree relative = 0.125 (cousins, etc.)
    #   -Fourth-degree relative = 0.0625
    #   -Fifth-degree relative = 0.03125
    # A good cutoff to use for Pi_Hat is 0.1875. This represents the halfway point between 2nd and 3rd degree relatives.

    #Get name of genotype file
    geno_name = input('\u001b[32;1m Please enter the name of the genotype files to run an IBD on (without bed/bim/fam extension: \u001b[0m')

    #Import module and call command.
    import genorelatives
    genorelatives.ibd(geno_name)

#GenoRelatives: Update FID or IID
elif to_do == '8':
    #Just making sure the user knows what is needed.
    print("The tab delimited text file for updating FID or IID should have four fields: \n"
          "1) Old FID\n"
          "2) Old IID\n"
          "3) New FID\n"
          "4) New IID")
    #Getting name of working file.
    geno_name = input('\u001b[32;1m Please enter the name of your genotype files that you would like to update FID/IID '
                      'in (without bed/bim/fam extension): \u001b[0m')
    #Name of file to be used to update the genotype files.
    update_id_filename = input('\u001b[34;1m Please enter the name of your text file for updating FID or IID '
                                '(with file extension): \u001b[0m')
    #Import module and call command.
    import genorelatives
    genorelatives.update_id(geno_name, update_id_filename)

#GenoRelatives: Update parental IDs
elif to_do == '9':
    #Just making sure the user knows what is needed.
    print("The tab delimited text file for updating parents should have four fields: \n"
          "1) FID\n"
          "2) IID\n"
          "3) Paternal IID\n"
          "4) Maternal IID")
    #Getting name of working file.
    geno_name = input('\u001b[32;1m Please enter the name of your genotype files that you would like to update parents '
                      'in (without bed/bim/fam extension): \u001b[0m')
    #Getting name of file to be used for update
    update_parents_filename = input('\u001b[34;1m Please enter the name of your text file for updating parents '
                                    '(with file extension): \u001b[0m')

    #Import module and call command.
    import genorelatives
    genorelatives.update_parental(geno_name, update_parents_filename)

#GenoHarmonize: Harmonize with 1000G
elif to_do == '10':
    #Get GenoName
    geno_name = input('\u001b[33;1m Please enter the name of the genotype file you would like to harmonize with 1000G Phase 3 '
                      '(without bed/bim/fam extension: \u001b[0m')

    # Harmonize with 1000G Phase 3
    import genoharmonize
    genoharmonize.harmonize_with_1000g(geno_name)

#GenoMerge: Merge with 1000G
elif to_do == '11':
    geno_name = input('\u001b[34;1m Please enter the name of the genotype files you would like to merge with 1000G '
                      '(without bed/bim/fam extension: \u001b[0m')
    import genomerge
    genomerge.merge(geno_name)

#PrepAdmixture: Prepares files for running ADMIXTURE, using 1000G as reference.
#Steps:
# Harmonize with 1000G Phase 3
# Merge with 1000G
# Prepare for ADMIXTURE with k = 3..9
# If on PSU cluster, can submit.
elif to_do == '10':
    # Make sure the reader knows what they're getting into.
    admixture_proceed_check = input("\u001b[32;1m This will merge your data with the 1000G data to and prepare files for"
                                    " an unsupervised ADMIXTURE analysis. Some cautions/notes before you perform this step:\n"
                                    "1) You should perform the steps 5-9 BEFORE this one (in roughly that order).\n"
                                    "2) IT WILL TAKE A LONG TIME (~10 hrs) TO MERGE YOUR DATA WITH 1000G\n"
                                    "3) There should not be related individuals when you perform admixture. If you have "
                                    "related individuals in your sample, you should create set lists so that the people "
                                    "in each set are unrelated (using information from the IBD analysis\n"
                                    "4) This will prepare files to run ADMIXTURE from k = 3 - 9. If you'd like other "
                                    "admixture runs performed, then you should change the PrepAdmixture.py code to "
                                    "reflect that.\n"
                                    "5) You must have a Penn State ACI cluster allocation to perofrm this step. We are "
                                    "using the cluster because ADMIXTURE takes a long time to run. I will ask you for "
                                    "your cluster name.\n"
                                    "6) This will write the files that you need, but you are responsible for the memory, "
                                    "node, and time usage (walltime = 150 hrs, nodes 1, ppn = 8, pmem = 8gb) and for "
                                    "putting them on the cluster and submitting them\n"
                                    "7) On the cluster, You will need the admixture program either on your path or in "
                                    "the same folder where you will submit this job.\n"
                                    "8) You will need to transfer the pbs files and genotype bed/bim/fam files to your "
                                    "cluster before running. I'll make a folder called 'Admixture' with all the files "
                                    "for you to transfer.\n"
                                    "Are you sure you want to proceed? (y/n): \u001b[0m").lower()
    if admixture_proceed_check in ('y', 'yes'):
        # Ask the user if they've already harmonized their data.
        harmonize_check = input('\u001b[33;1m Have you already harmonized your data with 1000G Phase 3? (y/n): \u001b[0m').lower()
        #If yes
        if harmonize_check in ('y', 'yes'):
           #Ask the user if they've already merged their data.
           merge_check = input('\u001b[34;1m Have you already merged your data with 1000G Phase 3 (y/n): \u001b[0m').lower()
           #If yes, proceed
           if merge_check in ('y', 'yes'):
               #This will kick it out of the if statement to proceed with the admixture prep.
           #If no, merge the data.
           elif merge_check in ('n','no'):
               #Ask for name of harmonized genotype files, which we will use to merge.
               geno_name = input('\u001b[33;1m Please enter the name of your harmonized genotype files that you would '
                                  'like to merge with 1000G (without bed/bim/fam extension): \u001b[0m')
               #Ask the user where their harmonized files are.
               harmonize_path = input('\u001b[34;1m Please enter the path name where your harmonized genotype files are '
                                      '(i.e. C:\\Users\\Julie White\\Box Sync\\Harmonized\\ etc.): \u001b[0m')
               #Get module for merging
               import genomerge
               genomerge.merge(geno_name, harmonized_path)
               #After, Should come back to this script and continue below whith admixture prep.

           #If user gives non yes or no response:
           else:
               sys.exit("Please answer yes or no to merge question. Quitting now.")

        #If they haven't harmonized, then harmonize and merge.
        elif harmonize_check in ('n', 'no'):
            #Ask for name of genotype file, which we will use to harmonize and then merge.
            geno_name = input('\u001b[33;1m Please enter the name of the genotype file you would like to harmonize, merge,'
                              ' then prepare for ADMIXTURE (without bed/bim/fam extension): \u001b[0m')

            #Harmonize with 1000G Phase 3
            import genoharmonize
            genoharmonize.harmonize_with_1000g(geno_name)

            # Since we've just harmonized, I know what the path is.
            harmonize_path = os.path.join(os.getcwd(), 'Harmonized_To_1000G')

            #Merge with 1000G Phase 3
            import genomerge
            genomerge.merge(geno_name, harmonize_path)

        else:
            sys.exit("Please give a yes or no answer. Quitting now.")

    #If they do not want to perform admixture at this time.
    elif admixture_proceed_check in ('n', 'no'):
        sys.exit("Okay we will not perform admixture at this time.")
    #If they give a non yes or no answer.
    else:
        sys.exit('Please give a yes or no answer. Quitting now.')

'''
#prepare for ADMIXTURE
elif to_do == '8':
    #Prepares files for an admixture run k = 3...9
    if admixture_proceed_check in ('y', 'yes'):
        #Get filename to run admixture on
        geno_name = input('\u001b[34;1m Please enter the name of the genotype files that you would like to perform '
                          'admixture on (without bed/bim/fam extension: \u001b[0m')

        #I based this formatting off of PSU cluster users, so they need to have a PSU cluster allocation.
        allocation_name = input('\u001b[35;1m Please enter the name of your cluster allocation: \u001b[0m')

        #Check if folder called 'Admixture' exists, if not, create it.
        if not os.path.exists('Admixture'):
            os.makedirs('Admixture')

        #Ask if they have relatives in their sample.
        relative_check = input('\u001b[32;1m Do you have relatives in your sample? Perhaps those identified in step 2. (y/n): \u001b[0m').lower()

        # We want the LD correction to be the same for all sets, so do this on the full genotype file and put it in the Admixture file.
        os.system('plink --bfile ' + geno_name + ' --indep-pairwise 50 10 0.1 --out Admixture/' + geno_name)

        if relative_check in ('y', 'yes'):
            # If they have relatives in their sample, get a list of the filenames for each set of people.
            user_set_input = input('\u001b[34;1m Please give me a comma separated list of your set list filenames (with '
                                   'file extenstion). I.e. dataset_setA.txt, dataset_setB.txt, etc. To do this,'
                                   'break up your entire dataset (not just related individuals) across sets, making sure '
                                   'that there are not related individuals within each set. These lists should be space '
                                   'or tab delimited with FID then IID: \u001b[0m')
            #Convert the user given list to a python list.
            set_list = user_set_input.split(', ')
            print(set_list)

            #Perform the admixture prep separately on each set.
            for i in range(0,len(set_list)):
                #Tell the user what they gave as file names and what I'm going to output as file names (SetA, SetB, SetC, etc.)
                set_name = chr(ord('a') + i).upper()
                print(set_list[i] + ' = Set' + set_name)

                #For each set, extract those people from the working genotype file and remove SNPs in LD. These files are what the user should put on the cluster.
                os.system('plink --bfile ' + geno_name + ' --keep ' + set_list[i] + ' --extract Admixture/'
                          + geno_name + '.prune.in --make-bed --out Admixture/' + geno_name + '_Set' + set_name + '_LDPruned')

                #For each set, write a pbs script for admixture k = 3..6
                with open('Admixture/' + geno_name + '_Set' + set_name + '_Admixture_k3to6.pbs', 'w') as file:
                    file.write('#PBS -l walltime=150:00:00\n'
                               '#PBS -l nodes=1:ppn=8\n'
                               '#PBS -l pmem=8gb\n'
                               '#PBS -A ' + allocation_name + '\n'
                               '#PBS -j oe\n'
                               'cd $PBS_O_WORKDIR\n'
                               'for K in {3..6}; do ./admixture --cv '
                               + geno_name + '_Set' + set_name + '_LDPruned.bed $K | tee '
                               + geno_name + '_Set' + set_name + '_LDPruned.log${K}.out; done')

                #For each set, write a pbs script for admixture k = 7..9
                with open('Admixture/' + geno_name + '_Set' + set_name + '_Admixture_k7to9.pbs', 'w') as file:
                    file.write('#PBS -l walltime=150:00:00\n'
                               '#PBS -l nodes=1:ppn=8\n'
                               '#PBS -l pmem=8gb\n'
                               '#PBS -A ' + allocation_name + '\n'
                               '#PBS -j oe\n'
                               'cd $PBS_O_WORKDIR\n'
                               'for K in {7..9}; do ./admixture --cv '
                               + geno_name + '_Set' + set_name + '_LDPruned.bed $K | tee '
                               + geno_name + '_Set' + set_name + '_LDPruned.log${K}.out; done')

            #Tell the user it's finished and give them directions.
            print("\u001b[36;1m Transfer " + geno_name + "_LDPruned bed/bim/fam files for each set, "
                  + geno_name + "_Admixture_k3to6.pbs, and " + geno_name + "_Admixture_7to9.pbs files for each set to the cluster.\n"
                  "Submit them using qsub " + geno_name + "Admixture_k3to6.pbs and qsub " + geno_name + "Admixture_k7to9.pbs\n"
                  "When you get your results, you should evaluate them to see which makes sense given your study "
                  "population and which has the lowest CV value\u001b[0m")

        #If there's no relatives, we can do all of this on the full genotype file.
        elif relative_check in ('n', 'no'):
            #For all people, create file that is LD pruned.
            os.system('plink --bfile ' + geno_name + ' --extract Admixture/' + geno_name + '.prune.in --make-bed --out Admixture/' + geno_name + '_LDPruned')

            #For all people, create a pbs file for admixture k = 3..6
            with open ('Admixture/' + geno_name + '_Admixture_k3to6.pbs', 'w') as file:
                file.write('#PBS -l walltime=150:00:00\n'
                           '#PBS -l nodes=1:ppn=8\n'
                           '#PBS -l pmem=8gb\n'
                           '#PBS -A ' + allocation_name + '\n'
                           '#PBS -j oe\n'
                           'cd $PBS_O_WORKDIR\n'
                           'for K in {3..6}; do ./admixture --cv ' + geno_name + '_LDPruned.bed $K | tee ' + geno_name + '_LDPruned.log${K}.out; done')

            #For all people, create a pbs file for admixture k = 7..9
            with open('Admixture/' + geno_name + '_Admixture_k7to9.pbs', 'w') as file:
                file.write('#PBS -l walltime=150:00:00\n'
                           '#PBS -l nodes=1:ppn=8\n'
                           '#PBS -l pmem=8gb\n'
                           '#PBS -A ' + allocation_name + '\n'
                           '#PBS -j oe\n'
                           'cd $PBS_O_WORKDIR\n'
                           'for K in {7..9}; do ./admixture --cv ' + geno_name + '_LDPruned.bed $K | tee ' + geno_name + '_LDPruned.log${K}.out; done')

            #End and give directions.
            print("\u001b[36;1m Transfer " + geno_name + "_LDPruned bed/bim/fam files, " + geno_name + "_Admixture_k3to6.pbs, and "
                  + geno_name + "_Admixture_7to9.pbs files to the cluster.\n"
                  "Submit them using qsub " + geno_name + "Admixture_k3to6.pbs and qsub " + geno_name + "Admixture_k7to9.pbs\n"
                  "When you get your results, you should evaluate them to see which makes sense given your study population and which has the lowest CV value\u001b[0m")

    #If the user says they do not want to perform admixture at this time
    elif admixture_proceed_check in ('n', 'no'):
        print("\u001b[36;1m Okay, we will not perform admixture at this time.\u001b[0m")

    #If the user does not return a recognizable answer.
    else:
        print("\u001b[36;1m Please answer yes or no\u001b[0m")

#Heterozygosity
elif to_do == '9':
    #Identifies individuals with extreme heterozygosity values (more than +- 3 SD)
    #Getting extra required modules
    import pandas as pd
    import numpy as np

    #Asking user what genotype files we're using
    geno_name = input('\u001b[34;1m Please enter the name of the genotype files that you would like to run a '
                      'heterozygosity check on (without bed/bim/fam extension: \u001b[0m')

    #Uses plink to calculate the heterozygosity, paying attention to geno and mind.
    os.system('plink --bfile ' + geno_name + ' --geno 0.1 --mind 0.1 --het --out ' + geno_name)

    #Read het file into pandas
    het_file = pd.read_csv(geno_name + '.het', sep='\s+', header=0)

    #Create new column with formula: (N(NM)-O(HOM))/N(NM)
    het_file['HET'] = (het_file['N(NM)'] - het_file['O(HOM)']) / het_file['N(NM)']

    #Getting standard deviation and average of HET column
    het_sd = np.std(het_file['HET'])
    het_avg = np.mean(het_file['HET'])

    #Add label 'keep' to people within 3*SD of the average het value, give 'remove' to everyone else.
    het_file['HET_Filter'] = np.where((het_file['HET'] > het_avg - 3 * het_sd) & (het_file['HET'] < het_avg + 3 * het_sd), 'Keep', 'Remove')
    #Write this file so the user has it.
    het_file.to_csv(geno_name + '.het', sep='\t', header=True, index=False)
    #Make a list of the people who pass the filter.
    het_keep = het_file[het_file['HET_Filter'] == 'Keep']
    #Write this file so that we can use it later to filter people out.
    het_keep[['FID', 'IID']].to_csv(geno_name + '_KeptAfterHetCheck.txt', sep='\t', header=False, index=False)

#Prephasing check
elif to_do == '10':
    #Ask the user what to run the phasing check on.
    geno_name = input('\u001b[32;1m Please enter the name of the genotype files that you would like to run a '
                      'phasing check on (without bed/bim/fam extension: \u001b[0m')

    #Since shapeit only works on linux or mac, we need to first check what system they are on.
    system_check = platform.system()

    if system_check == "Linux":
        #Now we need to know where they have shapeit, if they have it.
        shapeit_exists = input("\u001b[35;1m Do you already have the linux shapeit program unpacked? (y/n):  \u001b[0m").lower()
        #If yes, then ask for path of program
        if shapeit_exists in ('yes', 'y'):
            shapeit_path = input("\u001b[36;1m Please tell me the path where you have the shapeit program. "
                                 "i.e. C:\\Users\\Julie White\\Box Sync\\Software\\shapeit\\bin\\ \u001b[0m")
        #If no, download and unpack shapeit.
        elif shapeit_exits in ('no', 'n'):
            print('\u001b[36;1m Downloading shapeit to this directory now. \u001b[0m')
            urllib.request.urlretrieve('https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.20.Linux.static.tgz, '
                                       'shapeit.v2.r837.GLIBCv2.20.Linux.static.tgz')
            # Making directory to store program
            os.makedirs('Shapeit_v2.20_Linux_Static')
            #Unpacking
            os.system('tar -zxvf shapeit.v2.r837.GLIBCv2.20.Linux.static.tgz -C /Shapeit_v2.20_Linux_Static/')
        else:
            sys.exit('\u001b[36;1m You did not answer "y" or "no" when asked where shapeit was. Exiting now. \u001b[0m')

        #This part is unfinished.
        hap_legend_sample_path = input('\u001b[35;1m Please enter the pathname of where your 1000G Phase3 '
                                       'hap/legend/sample files are (i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\Hap_Legend_Sample etc.): \u001b[0m')
        
        
        os.system('./shapeit -check -B ' + geno_name + ' -M genetic_map_chr'
                  + [i] + '_combined_b37.txt --input-ref 1000GP_Phase3_chr'
                  + [i] + '.hap.gz 1000GP_Phase3_chr' + [i] + '.legend.gz 1000GP_Phase3.sample --output-log '
                  + geno_name + '_PhaseCheck')
    
    #If the user is on a mac
    elif system_check == "Darwin":
        #Ask if they already have shapeit
        shapeit_exists = input("\u001b[34;1m Great, do you already have the mac shapeit program unpacked? (y/n):  \u001b[0m").lower()
        if shapeit_exists in ('yes', 'y'):
            #Ask where shapeit is located.
            shapeit_path = input("\u001b[35;1m Please tell me the path where you have the shapeit program. "
                                 "i.e. C:\\Users\\Julie White\\Box Sync\\Software\\shapeit\\bin\\ \u001b[0m")
        elif shapeit_exits in ('no', 'n'):
            print('\u001b[35;1m Downloading shapeit to this directory now. \u001b[0m')
            #Download shapeit
            urllib.request.urlretrieve('https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.MacOSX.tgz, '
                                       'shapeit.v2.r837.MacOSX.tgz')
            #Create directory for shapeit.
            os.makedirs('Shapeit_v2.20_Mac')
            os.system('tar -zxvf shapeit.v2.r837.MacOSX.tgz -C /Shapeit_v2.20_Mac/')
        else:
            sys.exit('\u001b[35;1m You did not answer "y" or "no" when asked where shapeit was. Exiting now. \u001b[0m')

    elif system_check == ("Windows"):
        print("\u001b[35;1m I'm sorry, you need access to a linux or mac system to make this part work. If you have "
              "access to the Penn State clusters, you should run this script from there (they are linux). \u001b[0m")

    else:
        sys.exit("\u001b[35;1m I have detected that you are not running Linux, Mac, or Windows. Exiting now. \u001b[0m")

#Phasing ##### Unfinished.
elif to_do == '11':

    #Prepares files for phasing using shapeit
    phasing_proceed_check = input("\u001b[32;1m Some cautions/notes before you perform this step:\n"
                                    "1) You must perform step 1-6 before this step.\n"
                                    "2) You should have an ACI-B cluster allocation at Penn State to perform this step.\n"
                                    "3) This will write the files that you need, but you are responsible for the memory, node, and "
                                    "time usage (walltime = 150 hrs, nodes 1, ppn = 8, pmem = 8gb) and for putting them "
                                    "on the cluster and submitting them to SHAPEIT \n"
                                    "5) On the cluster, You will need the SHAPEIT program either on your path or in the same folder where "
                                    "you will submit this job.\n"
                                    "6) You will need to transfer the pbs file and genotype bed/bim/fam files to your cluster before running.\n"
                                    "7) Are you okay with all of this? (y/n): \u001b[0m").lower()
    if phasing_proceed_check in ('y', 'yes'):

        geno_name = input('\u001b[34;1m Please enter the name of the genotype files that you would like to phase on '
                          '(aka the name of the _HarmonizedTo1000G file produced from #5 (without bed/bim/fam extension: \u001b[0m')

        if not os.path.exists('Phasing'):
            os.makedirs('Phasing')

        shutil.copy2(geno_name + '.bed', 'Phasing')
        shutil.copy2(geno_name + '.bim', 'Phasing')
        shutil.copy2(geno_name + '.fam', 'Phasing')

        for file in glob.glob(r'plink*'):
            print(file)
            shutil.copy2(file, 'Phasing')

        os.chdir('Phasing')

#Nothing
elif to_do == '12':
    sys.exit("\u001b[36;1m You go, couch potato\u001b[0m")

else:
    print("\u001b[36;1m Please enter a number 1-9.\u001b[0m")
'''

