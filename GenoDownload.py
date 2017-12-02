def Plink():
    import os
    import platform
    import zipfile
    import shutil
    import urllib.request
    import sys

    #Get what system the user is using
    system_check = platform.system()
    #Get the version of that system
    architecture_check = platform.architecture()[0]

    if system_check == "Linux":
        if architecture_check == '64bit':
            # Download 64 bit Linux Plink 1.9 https://www.cog-genomics.org/static/bin/plink171114/plink_linux_x86_64.zip
            print("Downloading Linux Plink 1.9 to this directory now.")
            urllib.request.urlretrieve('https://www.cog-genomics.org/static/bin/plink171114/plink_linux_x86_64.zip',
                                       'Plink_1.9_Linux64.zip')
            # Making directory to store program
            os.makedirs('Plink_1.9_Linux64')
            # Unpacking to this directory
            with zipfile.ZipFile("Plink_1.9_Linux64.zip", "r") as zip_ref:
                zip_ref.extractall("Plink_1.9_Linux64")
            # Copy plink from archive folder to current working directory.
            shutil.copy2('Plink_1.9_Linux64/plink', os.getcwd())

        elif architecture_check == '32bit':
            # Download 32 bit Linux Plink 1.9 https://www.cog-genomics.org/static/bin/plink171114/plink_linux_i686.zip
            print("Downloading Linux Plink 1.9 to this directory now.")
            urllib.request.urlretrieve('https://www.cog-genomics.org/static/bin/plink171114/plink_linux_i686.zip',
                                       'Plink_1.9_Linux32.zip')
            # Making directory to store program
            os.makedirs('Plink_1.9_Linux32')
            # Unpacking to this directory
            with zipfile.ZipFile("Plink_1.9_Linux32.zip", "r") as zip_ref:
                zip_ref.extractall("Plink_1.9_Linux32")
            # Copy plink from archive folder to current working directory.
            shutil.copy2('Plink_1.9_Linux32/plink', os.getcwd())

        else:
            sys.exit("I'm sorry, I could not determine what Linux Plink version to download.")

    elif system_check == "Darwin":
        # Download Mac Plink 1.9 https://www.cog-genomics.org/static/bin/plink171114/plink_mac.zip
        print("Downloading Mac Plink 1.9 to this directory now.")
        urllib.request.urlretrieve('https://www.cog-genomics.org/static/bin/plink171114/plink_mac.zip',
                                   'Plink_1.9_Mac.zip')
        # Making directory to store program
        os.makedirs('Plink_1.9_Mac')
        # Unpacking to this directory
        with zipfile.ZipFile("Plink_1.9_Mac.zip", "r") as zip_ref:
            zip_ref.extractall("Plink_1.9_Mac")
        # Copy plink from archive folder to current working directory.
        shutil.copy2('Plink_1.9_Mac/plink', os.getcwd())

    elif system_check == "Windows":
        if architecture_check == '64bit':
            # Download 64 bit Windows Plink 1.9 https://www.cog-genomics.org/static/bin/plink171114/plink_win64.zip
            print("Downloading Windows Plink 1.9 to this directory now.")
            urllib.request.urlretrieve('https://www.cog-genomics.org/static/bin/plink171114/plink_win64.zip',
                                       'Plink_1.9_Win64.zip')
            # Making directory to store program
            os.makedirs('Plink_1.9_Win64')
            # Unpacking to this directory
            with zipfile.ZipFile("Plink_1.9_Win64.zip", "r") as zip_ref:
                zip_ref.extractall("Plink_1.9_Win64")
            # Copy plink from archive folder to current working directory.
            shutil.copy2('Plink_1.9_Win64/plink.exe', os.getcwd())

        elif architecture_check == '32bit':
            # Download 32 bit Windows Plink 1.9 https://www.cog-genomics.org/static/bin/plink171114/plink_win32.zip
            print("Downloading Windows Plink 1.9 to this directory now.")
            urllib.request.urlretrieve('https://www.cog-genomics.org/static/bin/plink171114/plink_win32.zip',
                                       'Plink_1.9_Win32.zip')
            # Making directory to store program
            os.makedirs('Plink_1.9_Win32')
            # Unpacking to this directory
            with zipfile.ZipFile("Plink_1.9_Win32.zip", "r") as zip_ref:
                zip_ref.extractall("Plink_1.9_Win32")
            # Copy plink from archive folder to current working directory.
            shutil.copy2('Plink_1.9_Win32/plink.exe', os.getcwd())

        else:
            sys.exit("I'm sorry, I could not determine what Windows Plink version to download.")

    else:
        sys.exit("I'm sorry, I could not determine what operating system you are running. Please download the latest "
              "version of Plink at https://www.cog-genomics.org/plink2")

def VCF1000GPhase3:
    print('Downloading 1000G Phase 3 VCF files now, putting them in "1000G_Phase3_VCF" folder. '
          'This will create a ~16G folder on your computer and WILL take a while.')

    # Create folder:
    if not os.path.exists('1000G_Phase3_VCF'):
        os.makedirs('1000G_Phase3_VCF')

    # List of VCF files that we're going to download.
    vcf_file_names = ['ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' % x for x in
                      range(1, 23)]
    vcf_file_names.extend(['ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'])

    # Open ftp connection
    ftp = ftplib.FTP('ftp.1000genomes.ebi.ac.uk')
    ftp.login()
    ftp.cwd('/vol1/ftp/release/20130502/')

    # Download files and put them in 1000G_Phase3_VCF folder.
    for filename in vcf_file_names:
        local_filename = os.path.join(os.getcwd(), '1000G_Phase3_VCF', filename)
        file = open(local_filename, 'wb')
        ftp.retrbinary('RETR ' + filename, file.write)
        file.close()
    # Once we are done downloading, close the connection to the ftp server.
    ftp.quit()

def HapLegendSample1000GPhase3():
    print('Downloading 1000G Phase 3 files now, putting them in "1000G_Phase3_HapLegendSample" folder. '
          'This will create a ~12G folder on your computer and WILL take a while.')

    # Create folder:
    if not os.path.exists('1000G_Phase3_HapLegendSample'):
        os.makedirs('1000G_Phase3_HapLegendSample')

        # List of file names that we're going to need.

        file_list = ['1000GP_Phase3.sample']
        file_list.extend(['1000GP_Phase3_chr%d.legend.gz' % x for x in range(1, 23)])
        file_list.extend(['1000GP_Phase3_chr%d.hap.gz' % x for x in range(1,23)])
        file_list.extend(['1000GP_Phase3_chrX_NONPAR.hap.gz', '1000GP_Phase3_chrX_NONPAR.legend.gz',
                          '1000GP_Phase3_chrX_PAR1.hap.gz', '1000GP_Phase3_chrX_PAR1.legend.gz',
                          '1000GP_Phase3_chrX_PAR2.hap.gz', '1000GP_Phase3_chrX_PAR2.legend.gz'])
        file_list.extend(['genetic_map_chr%d_combined_b37.txt' % x for x in range(1,23)])
        file_list.extend(['genetic_map_chrX_PAR1_combined.b37.txt', 'genetic_map_chrX_PAR2_combined.b37.txt',
                          'genetic_map_chrX_nonPAR_combined_b37.txt'])

        legend_server = "http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/"

        # Download files and put them in 1000G_Phase3_VCF folder.
        for filename in file_list:
            urllib.request.urlretrieve(os.path.join(legend_server, filename),
                                       os.path.join(os.getcwd(), '1000G_Phase3_HapLegendSample', filename))

def GenotypeHarmonizer:
    print('\u001b[36;1m Downloading genotype harmonizer now. \u001b[0m')
    # Download genotype harmonizer zip file
    urllib.request.urlretrieve(
        'http://www.molgenis.org/downloads/GenotypeHarmonizer/GenotypeHarmonizer-1.4.20-dist.zip',
        'GenotypeHarmonizer-1.4.20.zip')
    # Unzip Genotype Harmonizer
    zip_ref = zipfile.ZipFile('GenotypeHarmonizer-1.4.20.zip', 'r')
    zip_ref.extractall('GenotypeHarmonizer-1.4.20')
    zip_ref.close()