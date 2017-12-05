def phase(geno_name):
    import platform
    import urllib.request
    import os
    import sys

    # If it doesn't exist,
    if not os.path.exists('Phasing'):
        os.makedirs('Phasing')

    # Get original working directory in case we change it later.
    orig_wd = os.getcwd()

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

    # Remove SNPs and families with high Mendel error rates.
    os.system(plink + ' --bfile ' + geno_name + ' --me 0.05 0.1 --make-bed --out ' + geno_name + '_MeFilter')

    # Names of reference files needed.
    genetic_map_names = ['genetic_map_chr%d_combined_b37.txt' % x for x in range(1, 24)]
    hap_names = ['1000GP_Phase3_chr%d.hap.gz' % x for x in range(1,24)]
    legend_names = ['1000GP_Phase3_chr%d.legend.gz' % x for x in range(1,24)]
    geno_bed_names = ['Phasing/' + geno_name + '_MeFilter.chr%d.bed' % x for x in range(1,24)]
    geno_bim_names = ['Phasing/' + geno_name + '_MeFilter.chr%d.bim' % x for x in range(1, 24)]
    geno_fam_names = ['Phasing/' + geno_name + '_MeFilter.chr%d.fam' % x for x in range(1, 24)]

    # Perform phasing check per chromosome.
    for i in range(0,23):
        # Split the geno file into separate chromosomes and put it in the Phasing folder.
        os.system(plink + ' --bfile ' + geno_name + '_MeFilter --chr ' + str(i+1) + ' --make-bed --out Phasing/'
                  + geno_name + '_MeFilter.chr' + str(i+1))
        # Perform phasing check.
        os.system(os.path.join(shapeit_path,'shapeit') + ' -check --input-bed ' + geno_bed_names[i] + ' '
                  + geno_bim_names[i] + ' ' + geno_fam_names[i] + ' --input-map '
                  + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                  + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i]) + ' '
                  + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --output-log Phasing/'
                  + geno_name + '_MeFilter_PhaseCheck.chr' + str(i+1))

    for i in range(0,23):
        # Identify if there were Mendel errors (shouldn't be, but maybe)
        if (os.path.exists('Phasing/' + geno_name + '_MeFilter_PhaseCheck.chr' + str(i+1) + '.snp.me')) or \
                (os.path.exists('Phasing/' + geno_name + '_MeFilter_PhaseCheck.chr' + str(i+1) + '.ind.me')):
            sys.exit("There are SNPs or people in your sample with high rates of Mendel error. Please remove these and "
                     "rerun.")

    for i in range(0,23):
        # If a log file with snp.strand exclude is produced, then these SNPs need to be removed.
        if os.path.exists('Phasing/' + geno_name + '_MeFilter_PhaseCheck.chr' + str(i+1) + 'snp.strand.exclude'):
            print("You have SNPs in your sample that are not on the same strand as the reference or are not in the "
                  "reference. I'll remove these when phasing")

            # Run the phase after removing these SNPs
            os.system(os.path.join(shapeit_path, 'shapeit') + ' --input-bed ' + geno_bed_names[i] + ' '
                      + geno_bim_names[i] + ' ' + geno_fam_names[i] + ' --duohmm --input-map '
                      + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                      + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i]) + ' '
                      + os.path.join(ref_path, '1000GP_Phase3.sample') + '--exclude-snp Phasing/' + geno_name
                      + '_MeFilter_PhaseCheck.chr' + str(i+1) + '.snp.strand.exclude --output-log Phasing/'
                      + geno_name + '_MeFilter_PhasedTo1000G.chr' + str(i + 1))

        # If this file doesn't exist, then run the phase on all snps
        else:
            print("Phasing now.")
            # Run the phase
            os.system(os.path.join(shapeit_path, 'shapeit') + ' --input-bed ' + geno_bed_names[i] + ' '
                      + geno_bim_names[i] + ' ' + geno_fam_names[i] + ' --duohmm --input-map '
                      + os.path.join(ref_path, genetic_map_names[i]) + ' --input-ref '
                      + os.path.join(ref_path, hap_names[i]) + ' ' + os.path.join(ref_path, legend_names[i]) + ' '
                      + os.path.join(ref_path, '1000GP_Phase3.sample') + ' --output-log Phasing/'
                      + geno_name + 'MeFilter_PhasedTo1000G.chr' + str(i + 1))
