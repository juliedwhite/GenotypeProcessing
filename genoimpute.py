def prep(haps_list, sample_list):
    '''
    Required by the sanger imputation server.
    If not requesting pre-phasing, then all sites and samples should be phased with no missing data. - genophase
    All alleles on the forward strand - genoharmonize
    A single VCF file, not one file per-chromosome - this script
    Coordinates are on GRCh37 - user, before genoharmonize
    Chromosome names should be 1, 2, 3, etc… not chr1, chr2, chr3, etc… - this script
    Records are sorted by genomic position (chromosomal order is not important) - this script
    REF allele matches GRCh37. See the resources for help checking and fixing the REF allele - genoharmonize should do
    this, but in this script we will double check.
    Valid VCF - this script.
    '''

    import platform
    import os
    import sys
    import glob
    import urllib.request
    import subprocess
    try:
        import pip
    except ImportError:
        import genodownload
        genodownload.pip()

    vcf_name = input("Please tell me what you would like your phased VCF file to be called. I'm going to concatenate "
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

    if system_check in ("Linux", "Darwin"):
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

    # If they are running this on a windows machine, they cannot proceed because bcftools is *nix only.
    elif system_check == ("Windows"):
        sys.exit("\u001b[35;1m I'm sorry, you need access to a linux or mac system to make this part work. If you have "
                 "access to the Penn State clusters, you should run this script from there (they are linux). \u001b[0m")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("\u001b[35;1m I cannot detect the system you are working on. Exiting now. \u001b[0m")

    if system_check in ("Linux", "Darwin"):
        # Now we need to know where they have vcftools, if they have it.
        vcftools_exists = input("\u001b[35;1m Do you already have the program vcftools unpacked and installed? (y/n): "
                                "\u001b[0m").lower()
        # If yes, then ask for path of program
        if vcftools_exists in ('yes', 'y'):
            vcftools_path = input("\u001b[36;1m Please tell me the path where you have the vcftools program. "
                                  "i.e. C:\\Users\\Julie White\\Box Sync\\Software\\vcftools\\ \u001b[0m")
            # Setting path of where vcftools is since it's hard to use it without this.
            os.system('export PATH=$PATH:' + vcftools_path)
        # If no, download and unpack vcftools.
        elif vcftools_exists in ('no', 'n'):
            import genodownload
            # Download vcftools
            genodownload.vcftools()
        else:
            sys.exit('\u001b[36;1m You did not answer "y" or "no" when asked if you had vcftools. Exiting now. '
                     '\u001b[0m')

    # If they are running this on a windows machine, they cannot proceed because vcftools is *nix only.
    elif system_check == ("Windows"):
        sys.exit("\u001b[35;1m I'm sorry, you need access to a linux or mac system to make this part work. If you have "
                 "access to the Penn State clusters, you should run this script from there (they are linux). \u001b[0m")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("\u001b[35;1m I cannot detect the system you are working on. Exiting now. \u001b[0m")

    # Ask if they have the hg19 fasta files.
    fasta_exists = input('\u001b[34;1m Have you already downloaded the 1000G hg19 fasta file? (y/n): \u001b[0m').lower()

    if fasta_exists in ('y', 'yes'):
        # Ask the user where the fasta file is.
        fasta_path = input('\u001b[34;1m Please enter the pathname of where the your 1000G hg19 fasta file is '
                           '(i.e. C:\\Users\\Julie White\\Box Sync\\1000GP\\Fasta\\ etc.): '
                           '\u001b[0m')
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
    try:
        subprocess.call(['bgzip'])
    except OSError:
        import genodownload
        genodownload.htslib()

    # Use bcftools and bgzip to zip the vcf files so that they stop taking up space.
    for file in vcf_list:
        os.system('bgzip -c ' + file + ' > ' + file + '.vcf.gz')

    # Appending gz to the vcf_list names because we just gzipped them all.
    vcf_gz_list = [name + '.gz' for name in vcf_list]

    # Use bcftools to concatenate the per chromosome phased vcf files. Save this in SangerImputation and Phasing
    os.system('bcftools concat -Oz ' + print(" ".join(str(x) for x in vcf_gz_list)) + ' > SangerImputation/' + vcf_name
              + '.vcf.gz')
    # Save in phasing too.
    os.system('cp SangerImputation/' + vcf_name + '.vcf.gz Phasing/' + vcf_name + '.vcf.gz')

    # Use bcftools to change the chromosome names, if needed. Save this in SangerImputation.
    urllib.request.urlretrieve('https://imputation.sanger.ac.uk/www/plink2ensembl.txt',
                               'SangerImputation/plink2ensembl.txt')
    os.system('bcftools annotate -Oz --rename-chrs plink2ensembl.txt ' + vcf_name + '.vcf.gz > ' + vcf_name + '.vcf.gz')

    # Check VCF is sorted:
    os.system('bcftools index ' + vcf_name + '.vcf.gz')

    # Check that our data matches the reference allele of 1000G Phase 3. Fix those that do not match.
    os.system('bcftools norm -check-ref ws --fasta-ref ' + os.path.join(fasta_path,'human_g1k_v37.fasta.gz') + vcf_name
              + '.vcf.gz')

    # Double check sex and fix ploidy
    os.system('bcftools +guess-ploidy -g hg19 ' + vcf_name + '.vcf.gz > ' + vcf_name + '_SexEst.txt')
    os.system('bcftools +fixploidy ' + vcf_name + '.vcf.gz -Oz -o ' + vcf_name + '_SexUpdated.vcf.gz -- -s '
              + vcf_name + '_SexEst.txt')

    # Check that you have a valid vcf
    os.system('vcf-validator ' + vcf_name + '_SexUpdated.vcf.gz')

    sys.exit("Your vcf has been checked and is ready for imputation on the Sanger Imputation Server. Please follow the "
             "directions here (https://imputation.sanger.ac.uk/?instructions=1) for uploading your data. The file is "
             "called " + vcf_name + "_SexUpdated.vcf.gz")
