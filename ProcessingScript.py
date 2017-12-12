# Julie's Processing script for Shriver Lab genotype data
# Date 9.5.17

# WHAT THIS DOES #
# Sees what version of Python you're running. Downloads Python 3.6.3 if you're on linux and exits telling you to
    # download it if on Mac or Windows.
# Asks the user if they want to do the following:
    # 1) Download reference files and programs
    # 2) Update sex
    # 3) Produce new dataset with people and SNPs with missing call rate < 10%
    # 4) Run IBD matrix through plink
    # 5) Update family (FID) and individual IDs (IID)
    # 6) Update maternal and paternal IDs
    # 7) Harmonize with 1000G
    # 8) Filter for extreme (+-3SD) heterozygosity values
    # 9) Merge with 1000G
    # 10) Prepare for ADMIXTURE with k = 3...9.
    # 11) Run a phasing check and prepares files for phasing using SHAPEIT2.
    # 12) Prepares files for imputation using the Sanger Imputation Server.
    # 13) Nothing.

# MODULE FUNCTIONS #
# genodownload
    # Python 3.6.3
    # Plink 1.9
    # 1000G Phase 3 VCF
    # 1000G Phase 3 Hap/Legend/Sample'
    # GRCh37/hg19 1000G FASTA file'
    # Genotype Harmonizer'
    # pip'
    # snpflip'
    # shapeit'
    # vcftools'
    # bcftools'
    # htslib'

# genoqc
    # Update sex
    # Missing call rate
    # Heterozygosity check

# genorelatives
    # Run IBD to identify relatives
    # Update FID IID information
    # Update parental IDs

# genoharmonize
    # Harmonize with 1000G

# genomerge
    # Merge with 1000G

# genoadmixture
    # Prepare for admixture, submit job if on Penn Sate ACI-B cluster
    # Things that should be done before running this:
        # Missing call rate - genoqc.missing_call_rate()
        # Run IBD to identify relatives - genorelatives.ibd
        # Update FID IID information - genorelatives.update_id
        # Update parental IDs - genorelatives.update_parental
        # Harmonize with 1000G Phase 3 - genoharmonize.harmonize_with_1000G
        # Merge with 1000G - genomerge.merge1000g

# genophaseimpute
    # phase - Checks data and prepares for phasing using shapeit.
    # impute - checks phased data for imputation

        # Things that should be done before running this:
            # Must be on hg19, user should check this.
            # Missing call rate - genoqc.missing_call_rate
            # Heterozygosity check - genoqc.het
            # MAF threshold - genoharmonize.harmonize_with_1000G
            # HWE p-value - genoharmonize.harmonize_with_1000G
            # Set haploid genotypes (male chr X) as missing - DO I NEED TO DO THIS?
            # Check for gender mismatches - genoqc.update_sex, genophaseimpute.impute
            # Harmonize with 1000G Phase3 - genoharmonize.harmonize_with_1000G

# Getting the needed modules.
import os
import sys

if (sys.version_info > (3, 0)):
   pass
else:
    print("I've detected that you are not running Python3. This script was written with that in mind, so I'm going to "
          "download it (or ask you to download it) and exit. Then you should re-run this script.")
    import genodownload
    genodownload.python3()
    sys.exit("Exiting now, please re-run the script now that we've downloaded python3.")

# Ask the user what they'd like to do.
to_do = input('\u001b[31;1m What would you like to do?\n'
              '1) Download reference files or programs.\n'
              '2) Update sex. You need a file with FID, IID, Sex (M=1, F=2, Unknown=0) (in that order, no column '
              'headings)\n'
              '3) Produce a new dataset of people and SNPs with missing call rate < 10%\n'
              '4) Run an IBD analysis to identify relatives. All you need are plink bed/bim/fam files.\n'
              '5) Update FID or IID information. You need a file with the following information Old FID, Old IID, '
              'New FID, New IID.\n'
              '6) Update parental IDs. You need a file with FID, IID, Paternal IID, and Maternal IID.\n'
              '7) Harmonize with 1000G\n'
              '8) Filter for extreme (+-3SD) heterozygosity values\n'
              '9) Merge with 1000G\n'
              '10) Prepare for ADMIXTURE with 1000G Phase 3 files\n'
              '11) Run a phasing check and prepare files for phasing using SHAPEIT\n'
              '12) Prepare for imputation on the sanger imputation server\n'
              '13) Nothing. \n'
              'Please enter a number (i.e. 2): \u001b[0m')

# Use genodownload to figure out what to download.
if to_do == '1':
    # Get the module for downloading stuff.
    import genodownload
    # Call command
    genodownload.todownload()

# GenoQC: Update sex
elif to_do == '2':
    # Get name of genotype file
    geno_name = input("\u001b[32;1m Please enter the name of the plink genotype files you'd like to update sex in "
                      "(without bed/bim/fam extension: \u001b[0m")
    # Get name of file to be used for updating sex
    update_sex_filename = input('\u001b[34;1m Please enter the name of your text file for updating sex (with file '
                                'extension): \u001b[0m')
    # Import module where this command is.
    import genoqc
    # Call UpdateSex command using geno name and update sex filename as input
    genoqc.update_sex(geno_name, update_sex_filename)

# GenoQC: Clean dataset by missing call rate > 10%
elif to_do == '3':
    # Get name of genotype file
    geno_name = input('\u001b[32;1m Please enter the name of the genotype files (without bed/bim/fam extension:'
                      ' \u001b[0m')
    # Import module and call command
    import genoqc
    genoqc.missing_call_rate(geno_name)

# GenoRelatives: Run IBD
elif to_do == '4':
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

    # Get name of genotype file
    geno_name = input('\u001b[32;1m Please enter the name of the genotype files to run an IBD on (without bed/bim/fam extension: \u001b[0m')
    # Import module and call command.
    import genorelatives
    genorelatives.ibd(geno_name)

# GenoRelatives: Update FID or IID
elif to_do == '5':
    # Just making sure the user knows what is needed.
    print("The tab delimited text file for updating FID or IID should have four fields: \n"
          "1) Old FID\n"
          "2) Old IID\n"
          "3) New FID\n"
          "4) New IID")
    # Getting name of working file.
    geno_name = input('\u001b[32;1m Please enter the name of your genotype files that you would like to update FID/IID '
                      'in (without bed/bim/fam extension): \u001b[0m')
    # Name of file to be used to update the genotype files.
    update_id_filename = input('\u001b[34;1m Please enter the name of your text file for updating FID or IID '
                                '(with file extension): \u001b[0m')
    # Import module and call command.
    import genorelatives
    genorelatives.update_id(geno_name, update_id_filename)

# GenoRelatives: Update parental IDs
elif to_do == '6':
    # Just making sure the user knows what is needed.
    print("The tab delimited text file for updating parents should have four fields: \n"
          "1) FID\n"
          "2) IID\n"
          "3) Paternal IID\n"
          "4) Maternal IID")
    # Getting name of working file.
    geno_name = input('\u001b[32;1m Please enter the name of your genotype files that you would like to update parents '
                      'in (without bed/bim/fam extension): \u001b[0m')
    # Getting name of file to be used for update
    update_parents_filename = input('\u001b[34;1m Please enter the name of your text file for updating parents '
                                    '(with file extension): \u001b[0m')
    # Import module and call command.
    import genorelatives
    genorelatives.update_parental(geno_name, update_parents_filename)

# GenoHarmonize: Harmonize with 1000G
elif to_do == '7':
    print("Before we harmonize your data, please make sure your genotype data are on GRCh37/hg19 by comparing some of "
          "the positions in your bim file to the position for that rsid on the UCSC Genome Browser: "
          "https://genome.ucsc.edu/cgi-bin/hgGateway You should search your risd after selecting Human and "
          "Feb. 2009 (GRCh37/hg19). Your snp position should be directly in the middle of the snp ranges it gives when "
          "you search (i.e. rs1042522 at chr17:7579472-7579472). If you find that your snps don't match, then you "
          "should figure out what GRCh/hg version your snps are on and use the LiftOver tool "
          "(https://genome.ucsc.edu/cgi-bin/hgLiftOver) to change the coordinates to GRCh37/hg19.")

    coord_check = input("Have you checked that your data are on GRCh37/hg19? (y/n): ").lower()

    if coord_check in ('yes', 'y'):
        # Get name of genotypes.
        geno_name = input('\u001b[33;1m Please enter the name of the genotype file you would like to harmonize with '
                          '1000G Phase 3 (without bed/bim/fam extension: \u001b[0m')
        # Harmonize with 1000G Phase 3
        import genoharmonize
        genoharmonize.harmonize_with_1000g(geno_name)

    elif coord_check in ('no', 'n'):
        sys.exit("Please make sure your data are on GRCh37/hg19 then re-run this script. Exiting now.")

    else:
        sys.exit("Please answer 'yes' or 'no'. Exiting now.")

# GenoQC: Remove individuals with  extreme heterozygosity values (more than +- 3 SD)
elif to_do == '8':
    geno_name = input('\u001b[34;1m Please enter the name of the genotype files that you would like to run a '
                      'heterozygosity check on (without bed/bim/fam extension: \u001b[0m')
    # Call module and function.
    import genoqc
    genoqc.het(geno_name)

# GenoMerge: Merge with 1000G
elif to_do == '9':
    # Ask user genotype names.
    geno_name = input('\u001b[34;1m Please enter the name of the genotype files you would like to merge with 1000G '
                      '(without bed/bim/fam extension: \u001b[0m')
    # If there are genotype files with _HarmonizedTo1000G as ending in this working directory, then I know the path
    if os.path.exists(geno_name + '_HarmonizedTo1000G.bed'):
        harmonize_path = os.getcwd()
    else: # If I can't find the files in this working directory, ask the user where their harmonized files are.
        harmonize_path = input('\u001b[34;1m Please enter the path name where your harmonized genotype files are '
                               '(i.e. C:\\Users\\Julie White\\Box Sync\\Harmonized\\ etc.): \u001b[0m')
    #Import module and run.
    import genomerge
    genomerge.merge1000g(geno_name, harmonize_path)

# PrepAdmixture: Prepares files for running ADMIXTURE, using 1000G as reference.
# Steps:
#   Harmonize with 1000G Phase 3
#   Merge with 1000G
#   Prepare for ADMIXTURE with k = 3..9
elif to_do == '10':
    # Make sure the reader knows what they're getting into.
    admixture_proceed_check = input("\u001b[32;1m This will merge your data with the 1000G data to and prepare files "
                                    "for an unsupervised ADMIXTURE analysis. Some cautions/notes before you perform "
                                    "this step:\n"
                                    "1) You should perform the steps 2-9 BEFORE this one (in roughly that order).\n"
                                    "2) IT WILL TAKE A LONG TIME (~10 hrs) TO MERGE YOUR DATA WITH 1000G\n"
                                    "3) There should not be related individuals when you perform admixture. If you have"
                                    " related individuals in your sample, you should create set lists so that the "
                                    "people in each set are unrelated (using information from the IBD analysis\n"
                                    "4) This will prepare files to run ADMIXTURE from k = 3 - 9. If you'd like other "
                                    "admixture runs performed, then you should change the PrepAdmixture.py code to "
                                    "reflect that.\n"
                                    "5) You must have a Penn State ACI cluster allocation to perofrm this step. We are "
                                    "using the cluster because ADMIXTURE takes a long time to run. I will ask you for "
                                    "your cluster name.\n"
                                    "6) This will write the files that you need, but you are responsible for the "
                                    "memory, node, and time usage (walltime = 150 hrs, nodes 1, ppn = 8, pmem = 8gb) "
                                    "and for putting them on the cluster and submitting them\n"
                                    "7) On the cluster, You will need the admixture program either on your path or in "
                                    "the same folder where you will submit this job.\n"
                                    "8) You will need to transfer the pbs files and genotype bed/bim/fam files to your "
                                    "cluster before running. I'll make a folder called 'Admixture' with all the files "
                                    "for you to transfer.\n"
                                    "Are you sure you want to proceed? (y/n): \u001b[0m").lower()

    if admixture_proceed_check in ('y', 'yes'):
        # Ask the user if they've already harmonized their data.
        harmonize_check = input('\u001b[33;1m Have you already harmonized your data with 1000G Phase 3? (y/n): '
                                '\u001b[0m').lower()
        # If yes
        if harmonize_check in ('y', 'yes'):
            # Ask the user if they've already merged their data.
            merge_check = input('\u001b[34;1m Have you already merged your data with 1000G Phase 3 (y/n): '
                                '\u001b[0m').lower()
            # If yes, proceed, but first ask what they called the files.
            if merge_check in ('y', 'yes'):
                admix_name = input('\u001b[34;1m Please enter the name of the genotype files that you would like to '
                                   'perform admixture on (without bed/bim/fam extension: \u001b[0m')
            # If no, merge the data.
            elif merge_check in ('n','no'):
                # Ask for name of harmonized genotype files, which we will use to merge.
                geno_name = input('\u001b[33;1m Please enter the name of your harmonized genotype files that you would '
                                  'like to merge with 1000G (without bed/bim/fam extension): \u001b[0m')
                # Ask the user where their harmonized files are.
                harmonized_path = input('\u001b[34;1m Please enter the path name where your harmonized genotype files '
                                        'are (i.e. C:\\Users\\Julie White\\Box Sync\\Harmonized\\ etc.): \u001b[0m')
                # Get module for merging
                import genomerge
                genomerge.merge(geno_name, harmonized_path)
                # After, Should come back to this script and continue below with admixture prep.

                # Figure out what the final name of the merged file was.
                if os.path.exists(geno_name + '1000G.bed'):
                    admix_name = geno_name + '1000G.bed'
                elif os.path.exists(geno_name + '1000G_merge2.bed'):
                    admix_name = geno_name + '1000G_merge2.bed'
                elif os.path.exists(geno_name + '1000G_merge3.bed'):
                    admix_name = geno_name + '1000G_merge3.bed'
                else:
                    admix_name = input('\u001b[34;1m Please enter the name of the genotype files that you would like '
                                       'to perform admixture on (without bed/bim/fam extension: \u001b[0m')

            # If user gives non yes or no response:
            else:
                sys.exit("Please answer yes or no to merge question. Quitting now.")

        # If they haven't harmonized, then harmonize and merge.
        elif harmonize_check in ('n', 'no'):
            # Ask for name of genotype file, which we will use to harmonize and then merge.
            geno_name = input('\u001b[33;1m Please enter the name of the genotype file you would like to harmonize, '
                              'merge, then prepare for ADMIXTURE (without bed/bim/fam extension): \u001b[0m')

            # Harmonize with 1000G Phase 3
            import genoharmonize
            genoharmonize.harmonize_with_1000g(geno_name)

            # Since we've just harmonized, I know what the path is.
            harmonize_path = os.path.join(os.getcwd(), 'Harmonized_To_1000G')

            #Create new name because they've just been harmonized and I know what the ending should be.
            harmonized_name = geno_name + '_HarmonizedTo1000G_StrandChecked'
            # Merge with 1000G Phase 3
            import genomerge
            genomerge.merge(harmonized_name, harmonize_path)

            # Figure out what the final name of the merged file was.
            if os.path.exists(harmonized_name + '1000G.bed'):
                admix_name = (harmonized_name + '1000G.bed')
            elif os.path.exists(harmonized_name + '1000G_merge2.bed'):
                admix_name = harmonized_name + '1000G_merge2.bed'
            elif os.path.exists(harmonized_name + '1000G_merge3.bed'):
                admix_name = harmonized_name + '1000G_merge3.bed'
            else:
                admix_name = input('\u001b[34;1m Please enter the name of the genotype files that you would like to '
                                   'perform admixture on (without bed/bim/fam extension: \u001b[0m')
        # If user gives non-recognized answer.
        else:
            sys.exit("Please give a yes or no answer. Quitting now.")

        # Prep for admixture and done.
        import genoadmixture
        genoadmixture.prep(admix_name)

    # If they do not want to perform admixture at this time.
    elif admixture_proceed_check in ('n', 'no'):
        sys.exit("Okay we will not perform admixture at this time.")

    # If they give a non yes or no answer.
    else:
        sys.exit('Please give a yes or no answer. Quitting now.')

# GenoPhase: Run pre-phasing check; prepare and submit files for phasing.
elif to_do == '11':
    # Ask the user what to run the phasing check on.
    geno_name = input('\u001b[32;1m Please enter the name of the genotype files that you would like to run a '
                      'phasing check on (without bed/bim/fam extension). You should only do this after running steps '
                      '2-8: \u001b[0m')

    # I based this formatting off of PSU cluster users, so they need to have a PSU cluster allocation.
    allocation_name = input('\u001b[35;1m Please enter the name of your cluster allocation: \u001b[0m')
    # Import module
    import genophase
    # Call function
    genophase.phase(geno_name, allocation_name)

elif to_do == '12':
    # The user should only do this after phasing.
    print("I will prepare files for imputation now. It is very important that you do this AFTER phasing.")
    # Import module
    import genoimpute
    # Call function
    genoimpute.prep()

# Nothing
elif to_do == '13':
    sys.exit("\u001b[36;1m You go, couch potato\u001b[0m")

else:
    print("\u001b[36;1m Please enter a number 1-9.\u001b[0m")


