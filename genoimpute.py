def prep(haps_list, sample_list):
'''
Required by the sanger imputation server.
    If not requesting pre-phasing, then all sites and samples should be phased with no missing data. - genophase
    All alleles on the forward strand - in genoharmonize
    A single VCF file, not one file per-chromosome - this script
    Coordinates are on GRCh37 - user, before genoharmonize
    Chromosome names should be 1, 2, 3, etc… not chr1, chr2, chr3, etc… - Start HERE.

    Records are sorted by genomic position (chromosomal order is not important) -


    Valid VCF - here.
    REF allele matches GRCh37. See the resources for help checking and fixing the REF allele.
'''
    import platform
    import os
    import sys
    import glob
    try:
        import pip
    except ImportError:
        import genodownload
        genodownload.pip()

    try:
        import PyVCF
    except ImportError:
        import genodownload
        genodownload.PyVCF()

    vcf_name == input("Please tell me what you would like your phased VCF file to be called. I'm going to concatenate "
                      "all of the phased per chromosomes into one whole genome file and give it this name.: ")

    # Use glob to find files ending with the endings set in geno phase.
    vcf_list = glob.glob('Phasing/*_PhasedTo1000G.chr*.vcf')

    # If glob wasn't able to find files ending in the above, ask the user where the phased files are.
    if vcf_list == []:
        # Get path to phased files
        vcf_path = input("Please tell me where your phased vcf and haps/sample files are. Make sure they are in their "
                         "own folder, as this part of the script looks for anything with the ending '.vcf' "
                          "(i.e. C:\\Users\\Julie White\\Box Sync\\Genotypes\\Phasing\\ etc.): ")
        # Make list of anything with .vcf ending
        vcf_list = glob.glob(os.path.join(vcf_path,'*.vcf'))
    else:
        pass

    # If it doesn't exist, create SangerImputationFolder
    if not os.path.exists('SangerImputation'):
        os.makedirs('SangerImputation')

    # Since bcftools only works on linux or mac, we need to first check what system they are on.
    system_check = platform.system()

    if system_check == "Linux":
        # Now we need to know where they have bcftools, if they have it.
        bcftools_exists = input("\u001b[35;1m Do you already have the program bcftools unpacked and installed? (y/n): "
                                "\u001b[0m").lower()
        # If yes, then ask for path of program
        if bcftools_exists in ('yes', 'y'):
            bcftools_path = input("\u001b[36;1m Please tell me the path where you have the bcftools program. "
                                 "i.e. C:\\Users\\Julie White\\Box Sync\\Software\\bcftools\\ \u001b[0m")
            # Setting path of where bcftools is since it's hard to use it without this.
            os.system('export PATH=$PATH:' + bcftools_path)
        # If no, download and unpack bcftools.
        elif bcftools_exists in ('no', 'n'):
            import genodownload
            genodownload.bcftools()
        else:
            sys.exit('\u001b[36;1m You did not answer "y" or "no" when asked if you had bcftools. Exiting now. '
                     '\u001b[0m')

    # If the user is on a mac
    elif system_check == "Darwin":
        # Ask if they already have bcftools
        bcftools_exists = input("\u001b[34;1m Do you already have the bcftools program unpacked & installed? (y/n):  "
                               "\u001b[0m").lower()
        if bcftools_exists in ('yes', 'y'):
            # Ask where bcftools is located.
            bcftools_path = input("\u001b[35;1m Please tell me the path where you have the bcftools program. "
                                 "i.e. C:\\Users\\Julie White\\Box Sync\\Software\\bcftools\\ \u001b[0m")
            os.system('export PATH=$PATH:' + bcftools_path)
        elif bcftools_exists in ('no', 'n'):
            import genodownload
            genodownload.bcftools()
        else:
            sys.exit('\u001b[35;1m You did not answer "y" or "no" when asked if you had bcftools. Exiting now. '
                     '\u001b[0m')

    # If they are running this on a windows machine, they cannot proceed because bcftools is *nix only.
    elif system_check == ("Windows"):
        sys.exit("\u001b[35;1m I'm sorry, you need access to a linux or mac system to make this part work. If you have "
                 "access to the Penn State clusters, you should run this script from there (they are linux). \u001b[0m")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("\u001b[35;1m I cannot detect the system you are working on. Exiting now. \u001b[0m")

    # Use bcftools to concatenate the per chromosome phased vcf files. Save this in SangerImputation and Phasing
    os.system('bcftools concat -Oz ' + print(" ".join(str(x) for x in vcf_list)) + ' > SangerImputation/' + vcf_name
              + '.vcf')
    # Save in phasing too.
    os.system('cp SangerImputation/' + vcf_name + '.vcf.gz Phasing/' + vcf_name + '.vcf.gz')

    # Use bgzip and bcftools to zip and order file.
    # bgzip -c input.vcf > input.vcf.gz
    # bcftools index input.vcf.gz




