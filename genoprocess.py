# Getting the needed modules.
import os
import sys

from os.path import expanduser
home = expanduser("~")

if sys.version_info > (3, 0):
    pass
else:
    print("I've detected that you are not running Python3. This script was written with that in mind, so I'm going to "
          "download it (or ask you to download it) and exit. Then you should re-run this script.")
    import getpython
    getpython.python3()
    sys.exit("Exiting now, please move to " + os.path.join(home, 'software')
             + " and type in 'source InstallPython'. Then re-run this script.")

try:
    from colorama import init, Fore, Style
    init()
except ImportError:
    import genodownload
    genodownload.getcolorama()
    from colorama import init, Fore, Style
    init()


# Ask the user what they'd like to do.
print(Fore.RED + Style.BRIGHT +
      'What would you like to do?\n'
      '1) Download reference files or programs.\n'
      '2) Produce a new dataset of people and SNPs with missing call rate < 10%\n'
      '3) Estimate the sex of your samples. All you need are plink bed/bim/fam files.\n'
      '4) Update sex. You need a file with FID, IID, Sex (M=1, F=2, Unknown=0) (in that order, no column headings)\n'
      '5) Run an IBD analysis to identify relatives. All you need are plink bed/bim/fam files.\n'
      '6) Update FID or IID information. You need a file with the following information Old FID, Old IID, New FID, '
      'New IID.\n'
      '7) Update parental IDs. You need a file with FID, IID, Paternal IID, and Maternal IID.\n'
      '8) Harmonize with 1000G\n'
      '9) Filter for extreme (+-3SD) heterozygosity values\n'
      '10) Merge with 1000G\n'
      '11) Prepare for ADMIXTURE with 1000G Phase 3 files\n'
      '12) Run a phasing check and prepare files for phasing using SHAPEIT\n'
      '13) Prepare for imputation on the Sanger Imputation Server\n'
      '14) Extract imputation quality score and reference allele frequency from Sanger Imputation Server VCF files\n'
      '15) Plot imputation quality scores\n'
      '16) Nothing.')
to_do = input("Please enter the number that references the task you would like to perform above (i.e. 2): ")
print(Style.RESET_ALL)

# Use genodownload to figure out what to download.
if to_do == '1':
    # Get the module for downloading stuff.
    import genodownload
    # Call command
    genodownload.todownload()

# GenoQC: Clean dataset by missing call rate > 10%
elif to_do == '2':
    # Get name of genotype file
    print(Fore.BLUE + Style.BRIGHT)
    geno_name = input('Please enter the name of the genotype files (without bed/bim/fam extension: ')
    print(Style.RESET_ALL)

    # Import module and call command
    import genoqc
    genoqc.missing_call_rate(geno_name)

# GenoQC: Estimate sex
elif to_do == '3':
    # Get name of genotype file
    print(Fore.BLUE + Style.BRIGHT)
    geno_name = input("Please enter the name of the genotype files you'd like to estimate sex in "
                      "(without bed/bim/fam extension: ")
    print(Style.RESET_ALL)

    # Import module where this command is.
    import genoqc

    # Call UpdateSex command using geno name and update sex filename as input
    genoqc.estimate_sex(geno_name)

# GenoQC: Update sex
elif to_do == '4':
    # Get name of genotype file
    print(Fore.BLUE + Style.BRIGHT)
    geno_name = input("Please enter the name of the genotype files you'd like to update sex in "
                      "(without bed/bim/fam extension: ")
    print(Style.RESET_ALL)

    # Get name of file to be used for updating sex
    print(Fore.MAGENTA + Style.BRIGHT)
    update_sex_filename = input('Please enter the name of your text file for updating sex (with file extension): ')
    print(Style.RESET_ALL)
    # Import module where this command is.
    import genoqc
    # Call UpdateSex command using geno name and update sex filename as input
    genoqc.update_sex(geno_name, update_sex_filename)

# GenoRelatives: Run IBD
elif to_do == '5':
    # Identity-by-descent in Plink
    # This part of the script will prune for LD, calculate IBD, and exclude individuals who have IBD < 0.2
    # The IBD results will have .genome appended to your file name. I have also included a line to convert the IBD
    #   results from whitespace to tab delimited. This will have .tab.genome appended to your filename.

    # Important values of Pi-hat
    #   -First-degree relative = 0.5 (full sibs, parent-offspring)
    #   -Second-degree relative = 0.25 (half-sibs, uncle/aunt-nephew/niece, grandparent-grandchild)
    #   -Third-degree relative = 0.125 (cousins, etc.)
    #   -Fourth-degree relative = 0.0625
    #   -Fifth-degree relative = 0.03125
    # A good cutoff to use for Pi_Hat is 0.1875. This represents the halfway point between 2nd and 3rd degree relatives.

    # Get name of genotype file
    print(Fore.BLUE + Style.BRIGHT)
    geno_name = input('Please enter the name of the genotype files to run an IBD on (without bed/bim/fam extension: ')
    print(Style.RESET_ALL)
    # Import module and call command.
    import genorelatives
    genorelatives.ibd(geno_name)

# GenoRelatives: Update FID or IID
elif to_do == '6':
    # Just making sure the user knows what is needed.
    print("The tab delimited text file for updating FID or IID should have four fields: \n"
          "1) Old FID\n"
          "2) Old IID\n"
          "3) New FID\n"
          "4) New IID")
    # Getting name of working file.
    print(Fore.BLUE + Style.BRIGHT)
    geno_name = input('Please enter the name of your genotype files that you would like to update FID/IID in '
                      '(without bed/bim/fam extension): ')
    print(Style.RESET_ALL)

    # Name of file to be used to update the genotype files.
    print(Fore.MAGENTA + Style.BRIGHT)
    update_id_filename = input('Please enter the name of your text file for updating FID or IID '
                               '(with file extension): ')
    print(Style.RESET_ALL)

    # Import module and call command.
    import genorelatives
    genorelatives.update_id(geno_name, update_id_filename)

# GenoRelatives: Update parental IDs
elif to_do == '7':
    # Just making sure the user knows what is needed.
    print("The tab delimited text file for updating parents should have four fields: \n"
          "1) FID\n"
          "2) IID\n"
          "3) Paternal IID\n"
          "4) Maternal IID")
    # Getting name of working file.
    print(Fore.BLUE + Style.BRIGHT)
    geno_name = input('Please enter the name of your genotype files that you would like to update parents in '
                      '(without bed/bim/fam extension): ')
    print(Style.RESET_ALL)

    # Getting name of file to be used for update
    print(Fore.MAGENTA + Style.BRIGHT)
    update_parents_filename = input('Please enter the name of your text file for updating parents '
                                    '(with file extension): ')
    print(Style.RESET_ALL)

    # Import module and call command.
    import genorelatives
    genorelatives.update_parental(geno_name, update_parents_filename)

# GenoHarmonize: Harmonize with 1000G
elif to_do == '8':
    print(Fore.BLUE + Style.BRIGHT
          + "Before we harmonize your data, please make sure your genotype data are on GRCh37/hg19 by comparing some "
            "of the positions in your bim file to the position for that rsid on the UCSC Genome Browser: "
            "https://genome.ucsc.edu/cgi-bin/hgGateway You should search your risd after selecting Human and Feb. 2009 "
            "(GRCh37/hg19). Your snp position should be directly in the middle of the snp ranges it gives when you "
            "search (i.e. rs1042522 at chr17:7579472-7579472). If you find that your snps don't match, then you should "
            "figure out what GRCh/hg version your snps are on and use the LiftOver tool "
            "(https://genome.ucsc.edu/cgi-bin/hgLiftOver) to change the coordinates to GRCh37/hg19.")
    print(Fore.CYAN + Style.RESET_ALL)
    coord_check = input("Have you checked that your data are on GRCh37/hg19? (y/n): ").lower()
    print(Style.RESET_ALL)

    if coord_check in ('yes', 'y'):
        # Get name of genotypes.
        print(Fore.GREEN)
        geno_name = input('Please enter the name of the genotype file you would like to harmonize with 1000G Phase 3 '
                          '(without bed/bim/fam extension): ')
        print(Style.RESET_ALL)

        # Getting required reference files and genotype harmonizer program
        # Ask the user if they already have the 1000G Phase 3 vcf files.
        print(Fore.BLUE + Style.BRIGHT)
        vcf_exists = input('Have you already downloaded the 1000G Phase3 VCF files? (y/n): ').lower()
        print(Style.RESET_ALL)
        # If yes, then get user's path to VCF files
        if vcf_exists in ('y', 'yes'):
            print(Fore.MAGENTA + Style.BRIGHT)
            vcf_path = input('Please enter the pathname of where your 1000G VCF files are '
                             '(i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\VCF\\ etc.): ')
            print(Style.RESET_ALL)
        # If no, then download vcf files.
        elif vcf_exists in ('n', 'no'):
            import genodownload
            genodownload.vcf_1000g_phase3()
            # Saving VCF path
            vcf_path = os.path.join(os.getcwd(), '1000G_Phase3_VCF')
        else:
            sys.exit("Please answer yes or no. Quitting now because no VCF files.")

        # Ask the user if they already have the 1000G Phase 3 Hap/Legend/Sample files.
        print(Fore.GREEN)
        legend_exists = input('Have you already downloaded the 1000G Phase 3 Hap/Legend/Sample files? (y/n): ').lower()
        print(Style.RESET_ALL)
        # If yes, then get legend_path
        if legend_exists in ('y', 'yes'):
            print(Fore.CYAN)
            legend_path = input('Please enter the pathname of where your 1000G legend files are '
                                '(i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\ etc.): ')
            print(Style.RESET_ALL)
        # If no, then download legend files
        elif legend_exists in ('n', 'no'):
            import genodownload
            genodownload.hls_1000g_phase3()
            # Saving legend path
            legend_path = os.path.join(os.getcwd(), '1000G_Phase3_HapLegendSample')
        else:
            sys.exit('Please answer yes or no. Quitting now because no legend files.')

        # Ask if they have the hg19 fasta files.
        print(Fore.BLUE + Style.BRIGHT)
        fasta_exists = input('Have you already downloaded the 1000G hg19 fasta file? (y/n): ').lower()
        print(Style.RESET_ALL)
        # If yes, then get fasta path
        if fasta_exists in ('y', 'yes'):
            # Ask the user where the fasta file is.
            print(Fore.MAGENTA + Style.BRIGHT)
            fasta_path = input('Please enter the pathname of where the your 1000G hg19 fasta file is '
                               '(i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\Fasta\\ etc.): ')
            print(Style.RESET_ALL)
        # If no, then download fasta files.
        elif fasta_exists in ('n', 'no'):
            import genodownload
            genodownload.fasta_1000G_hg19()
            # Saving fasta path
            fasta_path = os.path.join(os.getcwd(), '1000G_hg19_fasta')
        else:
            sys.exit('Please answer yes or no. Quitting now because no fasta file.')

        # Ask if they have genotype harmonizer.
        print(Fore.GREEN)
        harmonizer_exists = input('Have you already downloaded Genotype Harmonizer? (y/n): ').lower()
        print(Style.RESET_ALL)
        # If yes, then get path to genotype harmonizer.
        if harmonizer_exists in ('y', 'yes'):
            print(Fore.CYAN)
            harmonizer_path = input('Please enter the pathname of where the GenotypeHarmonizer.jar file is '
                                    '(i.e. C:\\Users\\Julie White\\Box Sync\\Software\\GenotypeHarmonizer-1.4.20\\): ')
            print(Style.RESET_ALL)
        # If no, then download genotype harmonizer
        elif harmonizer_exists in ('n', 'no'):
            import genodownload
            genodownload.genotype_harmonizer()
            # Harmonize path now that we've downloaded it.
            harmonizer_path = os.path.join(os.getcwd(), 'GenotypeHarmonizer-1.4.20/GenotypeHarmonizer-1.4.20-SNAPSHOT/')
        else:
            sys.exit('Please write yes or no. Quitting now because no Genotype Harmonizer.')

        # Ask if the user is on the cluster right now to determine if we should submit the files for them
        print(Fore.BLUE + Style.BRIGHT)
        on_cluster = input('Are you currently running this from the Penn State ACI-B cluster? If yes, I make this '
                           'process as a job and submit it to run. If you are not on the cluster, then this will run '
                           'locally and will take approximately 15 hours. (y/n): ').lower()
        print(Style.RESET_ALL)
        # If they are on the cluster, then run as a job.
        if on_cluster in ('yes', 'y'):
            # I based this formatting off of PSU cluster users, so they need to have a PSU cluster allocation.
            print(Fore.MAGENTA + Style.BRIGHT)
            allocation_name = input('Please enter the name of your cluster allocation: ')
            print(Style.RESET_ALL)
            # Run
            import genoharmonize
            genoharmonize.cluster(geno_name, allocation_name, harmonizer_path, vcf_path, legend_path, fasta_path)

        # If they are not on the cluster, then run on their local machine.
        elif on_cluster in ('no', 'n'):
            # Run harmonization script on their local machine.
            import genoharmonize
            genoharmonize.local(geno_name, harmonizer_path, vcf_path, legend_path, fasta_path)

    # If they have not checked that they are on hg19, quit
    elif coord_check in ('no', 'n'):
        sys.exit("Please make sure your data are on GRCh37/hg19 then re-run this script. Exiting now.")
    # If they give a bad answer.
    else:
        sys.exit("Please answer 'yes' or 'no'. Exiting now.")


# GenoQC: Remove individuals with  extreme heterozygosity values (more than +- 3 SD)
elif to_do == '9':
    print(Fore.BLUE + Style.BRIGHT)
    geno_name = input('Please enter the name of the genotype files that you would like to run a heterozygosity check '
                      'on (without bed/bim/fam extension: ')
    print(Style.RESET_ALL)

    # Call module and function.
    import genoqc
    genoqc.het(geno_name)

# GenoMerge: Merge with 1000G
elif to_do == '10':
    # Ask user genotype names.
    print(Fore.BLUE + Style.BRIGHT)
    harmonized_name = input('Please enter the name of the harmonized genotype files you would like to merge with 1000G '
                            '(without bed/bim/fam extension: ')
    print(Style.RESET_ALL)

    # If there are genotype files with _HarmonizedTo1000G as ending in this working directory, then I know the path
    if os.path.exists(harmonized_name + '.bed'):
        harmonize_path = os.getcwd()
    else:  # If I can't find the files in this working directory, ask the user where their harmonized files are.
        print(Fore.MAGENTA + Style.BRIGHT)
        harmonize_path = input('Please enter the path name where your harmonized genotype files are '
                               '(i.e. C:\\Users\\Julie White\\Box Sync\\Harmonized\\ etc.): ')
        print(Style.RESET_ALL)

    # Import module and run.
    import genomerge
    genomerge.merge1000g(harmonized_name, harmonize_path)

# PrepAdmixture: Prepares files for running ADMIXTURE, using 1000G as reference.
# Steps:
#   Harmonize with 1000G Phase 3
#   Merge with 1000G
#   Prepare for ADMIXTURE with k = 3..9
elif to_do == '11':
    # Make sure the reader knows what they're getting into.
    print("This will merge your data with the 1000G data to and prepare files for an unsupervised ADMIXTURE analysis. "
          "Some cautions/notes before you perform this step:\n"
          "1) You should perform the steps 2-10 BEFORE this one (in roughly that order).\n"
          "2) IT WILL TAKE A LONG TIME (~10 hrs) TO MERGE YOUR DATA WITH 1000G\n"
          "3) There should not be related individuals when you perform admixture. If you have related individuals in "
          "your sample, you should create set lists so that the people in each set are unrelated (using information "
          "from the IBD analysis\n"
          "4) This will prepare files to run ADMIXTURE from k = 3 - 9. If you'd like other admixture runs performed, "
          "then you should (carefully) change the genoadmixture.py code, or the produced .pbs files, to reflect that.\n"
          "5) You must have a Penn State ACI cluster allocation to perform this step. We are using the cluster because "
          "ADMIXTURE takes a long time to run. I will ask you for your cluster name.\n")
    print(Fore.BLUE + Style.BRIGHT)
    admixture_proceed_check = input("Are you sure you want to proceed? (y/n): ").lower()
    print(Style.RESET_ALL)

    if admixture_proceed_check in ('y', 'yes'):
        # Ask the user if they've already harmonized their data.
        print(Fore.MAGENTA + Style.BRIGHT)
        harmonize_check = input('Have you already harmonized your data with 1000G Phase 3? (y/n): ').lower()
        print(Style.RESET_ALL)

        # If yes
        if harmonize_check in ('y', 'yes'):
            # Ask the user if they've already merged their data.
            print(Fore.GREEN)
            merge_check = input('Have you already merged your data with 1000G Phase 3 (y/n): ').lower()
            print(Style.RESET_ALL)

            # If yes, proceed, but first ask what they called the files.
            if merge_check in ('y', 'yes'):
                print(Fore.CYAN)
                admix_name = input('Please enter the name of the genotype files that you would like to perform '
                                   'admixture on (without bed/bim/fam extension: ')
                print(Style.RESET_ALL)

                # Prep for admixture and done.
                import genoadmixture
                genoadmixture.prep(admix_name)

            # If no, they must merge the data
            elif merge_check in ('n', 'no'):
                print(Fore.BLUE + Style.BRIGHT)
                sys.exit("Please merge your data with 1000G and then re-run this part of the script.")

            # If user gives non yes or no response:
            else:
                sys.exit("Please answer yes or no to merge question. Quitting now.")

        # If they haven't harmonized, the user must harmonize and merge.
        elif harmonize_check in ('n', 'no'):
            print(Fore.BLUE + Style.BRIGHT)
            sys.exit("Please harmonize then merge your data with 1000G and then re-run this part of the script.")

        # If user gives non-recognized answer.
        else:
            sys.exit("Please give a yes or no answer. Quitting now.")

    # If they do not want to perform admixture at this time.
    elif admixture_proceed_check in ('n', 'no'):
        sys.exit("Okay we will not perform admixture at this time.")

    # If they give a non yes or no answer.
    else:
        sys.exit('Please give a yes or no answer. Quitting now.')


# GenoPhase: Run pre-phasing check; prepare and submit files for phasing.
elif to_do == '12':
    # Ask the user what to run the phasing check on.
    print(Fore.BLUE + Style.BRIGHT)
    geno_name = input('Please enter the name of the genotype files that you would like to run a phasing check on '
                      '(without bed/bim/fam extension). You should only do this after running steps 2-8: ')
    print(Style.RESET_ALL)

    # I based this formatting off of PSU cluster users, so they need to have a PSU cluster allocation.
    print(Fore.MAGENTA + Style.BRIGHT)
    allocation_name = input('Please enter the name of your cluster allocation: ')
    print(Style.RESET_ALL)

    # Import module
    import genophaseimpute
    # Call function
    genophaseimpute.phase(geno_name, allocation_name)

# Prepare for imputation.
elif to_do == '13':
    # The user should only do this after phasing.
    print("I will prepare files for imputation now. It is very important that you do this AFTER phasing.")
    # Import module
    import genophaseimpute
    # Call function
    genophaseimpute.impute()

# Get info from imputed files.
elif to_do == '14':
    print(Fore.BLUE + Style.BRIGHT)
    imputed_path = input('Please enter the path for your Sanger Imputed VCF files '
                         '(e.g C:\\Users\\Julie White\\Box Sync\\SangerImputation\\) : ')
    print(Style.RESET_ALL)

    print(Fore.MAGENTA + Style.BRIGHT)
    geno_name = input('Please enter the name of the imputed VCF files that you would like to extract imputation '
                      'quality scores and reference allele frequency from. The VCF file names should be in the format '
                      'YourStudy_chr1, YourStudy_chr2...YourStudy_chr23" and I would like to know what "YourStudy" is, '
                      'I will fill in the _chr# in the script.: ')
    print(Style.RESET_ALL)

    # Import module
    import genophaseimpute

    # Call function
    genophaseimpute.getinfo(imputed_path, geno_name)

# Make plots of imputation quality scores.
elif to_do == '15':
    print(Fore.BLUE + Style.BRIGHT)
    info_path = input("Please tell me where your .INFO files produced by step #13 are. Make sure they are in a folder "
                      "with no other .INFO files, as this part of the script looks for anything with the ending '.INFO'"
                      " (e.g C:\\Users\\Julie White\\Box Sync\\SangerImputation\\ etc.): ")
    print(Style.RESET_ALL)

    # Import module
    import genophaseimpute

    # Call function
    genophaseimpute.qualscoreplot(info_path)

# Nothing
elif to_do == '16':
    sys.exit("You go, couch potato")

else:
    print("Please enter a number 1-16.")
