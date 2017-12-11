def todownload():
    import sys

    item = input('\u001b[31;1m What would you like to download?\n'
                 '1) Python 3.6.3\n'
                 '2) Plink\n'
                 '3) 1000G Phase 3 VCF\n'
                 '4) 1000G Phase 3 Hap/Legend/Sample\n'
                 '5) GRCh37/hg19 1000G FASTA file\n'
                 '6) Genotype Harmonizer\n'
                 '7) pip\n'
                 '8) snpflip\n'
                 '9) shapeit\n'
                 '10) vcftools\n'
                 '11) bcftools\n'
                 '12) htslib\n'
                 '13) Nothing\n'
                 'Please enter a number (i.e. 2): \u001b[0m')

    if item == '1':
        python3()
    elif item == '2':
        plink()
    elif item == '3':
        vcf_1000g_phase3()
    elif item == '4':
        hls_1000g_phase3()
    elif item == '5':
        fasta_1000G_hg19()
    elif item == '6':
        genotype_harmonizer()
    elif item == '7':
        pip()
    elif item == '8':
        snpflip()
    elif item == '9':
        shapeit()
    elif item == '10':
        vcftools()
    elif item == '11':
        bcftools()
    elif item == '12':
        htslib()
    elif item == '13':
        sys.exit("Exiting now")
    else:
        sys.exit("Quitting because you did not give a recognizable number when asked what to download.")


def python3():
    import platform
    import urllib.request
    import os
    import sys

    from os.path import expanduser
    home = expanduser("~")

    # Get what system the user is using
    system_check = platform.system()

    if system_check == "Linux":
        if not os.path.exists(os.path.join(home,'software')):
            os.makedirs(os.path.join(home,'software'))

        print("\u001b[36;1m Downloading python3 to " + os.path.join(home,'software') + " now. \u001b[0m")
        urllib.request.urlretrieve('https://www.python.org/ftp/python/3.6.3/Python-3.6.3.tgz',
                                   os.path.join(home,'software/Python-3.6.3.tgz'))
        # Unpacking
        os.system('tar -zxvf Python-3.6.3.tgz -C ' + os.path.join(home,'software'))
        # Moving into the bcftools folder
        os.chdir(os.path.join(home,'software/Python-3.6.3'))
        # Running configuration and installation steps
        os.system('./configure --prefix=' + os.path.join(home,'software/'))
        os.system('make')
        os.system('make install')
        os.system('export PATH=$PATH:' + os.path.join(home, 'software/bin'))
        os.system('export PYTHONPATH=' + os.path.join(home,'software/Python-3.6.3/'))

    # If the user is on a mac
    elif system_check == "Darwin":
        sys.exit("I've detected that you are working on a Mac. Your best bet is to follow these directions "
                 "(http://docs.python-guide.org/en/latest/starting/install3/osx/) and download"
                 "python3 using xcode and homebrew. You will probably need homebrew later in this script too.")

    # If they are running this on a windows machine
    elif system_check == "Windows":
        sys.exit("\u001b[35;1m I've detected that you are working on a Windows computer. Your best bet is to follow "
                 "these directions (http://docs.python-guide.org/en/latest/starting/install3/win/) and download "
                 "python3 mnaually. As a word of warning, you can do many functions of this script on a Windows "
                 "computer, but some (admixture, phasing, imputation) requires access to a linux based computational "
                 "cluster (like the Penn State ACI-B system) \u001b[0m")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("\u001b[35;1m I cannot detect the system you are working on. Exiting now. \u001b[0m")


def plink():
    import os
    import platform
    import zipfile
    import shutil
    import urllib.request
    import sys

    # Get what system the user is using
    system_check = platform.system()
    # Get the version of that system
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


def vcf_1000g_phase3():
    import os
    import ftplib

    print('Downloading 1000G Phase 3 VCF files now, putting them in "1000G_Phase3_VCF" folder. '
          'This will create a ~16G folder on your computer and WILL take a while.')

    # Create folder:
    if not os.path.exists('1000G_Phase3_VCF'):
        os.makedirs('1000G_Phase3_VCF')

    # List of VCF files that we're going to download.
    vcf_file_names = ['ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' % x for x in
                      range(1, 23)]
    vcf_file_names.extend(['ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz',
                           'ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz'])

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


def hls_1000g_phase3():
    import os
    import urllib.request
    import tarfile

    print('Downloading 1000G Phase 3 files now, putting them in "1000G_Phase3_HapLegendSample" folder. '
          'This will create a ~12G folder on your computer and WILL take a while.')

    orig_wd = os.getcwd()

    # Create folder:
    if not os.path.exists('1000G_Phase3_HapLegendSample'):
        os.makedirs('1000G_Phase3_HapLegendSample')

        # Where the files are
        legend_server = "http://mathgen.stats.ox.ac.uk/impute/"

        # List of file names that we're going to need.
        file_list = ['1000GP_Phase3.tgz','1000GP_Phase3_chrX.tgz']

        # Download files and put them in 1000G_Phase3_VCF folder.
        for filename in file_list:
            urllib.request.urlretrieve(os.path.join(legend_server, filename),
                                       os.path.join(os.getcwd(), '1000G_Phase3_HapLegendSample', filename))

        # Change directory where tgz files are
        os.chdir('1000G_Phase3_HapLegendSample')
        # Unpack tar files
        tar = tarfile.open('1000GP_Phase3.tgz', 'r:gz')
        for item in tar:
            tar.extract(item)
            print('Done extracting' + str(item))
        # Unpack tar files
        tar = tarfile.open('1000GP_Phase3_chrX.tgz', 'r:gz')
        for item in tar:
            tar.extract(item)
            print('Done extracting' + str(item))

        # Change back to original working directory.
        os.chdir(orig_wd)


def fasta_1000G_hg19():
    print('Downloading 1000G hg19 fasta file now, putting it in the "1000G_hg19_fasta" folder.')

    # Create folder:
    if not os.path.exists('1000G_hg19_fasta'):
        os.makedirs('1000G_hg19_fasta')

    # Open ftp connection
    ftp = ftplib.FTP('ftp.1000genomes.ebi.ac.uk')
    ftp.login()
    ftp.cwd('/vol1/ftp/technical/reference/')

    # Download files and put them in 1000G_hg19_fasta.
    file = open('1000G_hg19_fasta/human_g1k_v37.fasta.gz', 'wb')
    ftp.retrbinary('RETR human_g1k_v37.fasta.gz', file.write)
    ftp.quit()
    file.close()


def genotype_harmonizer():
    import urllib.request
    import zipfile

    print('\u001b[36;1m Downloading genotype harmonizer now. \u001b[0m')
    # Download genotype harmonizer zip file
    urllib.request.urlretrieve(
        'http://www.molgenis.org/downloads/GenotypeHarmonizer/GenotypeHarmonizer-1.4.20-dist.zip',
        'GenotypeHarmonizer-1.4.20.zip')
    # Unzip Genotype Harmonizer
    zip_ref = zipfile.ZipFile('GenotypeHarmonizer-1.4.20.zip', 'r')
    zip_ref.extractall('GenotypeHarmonizer-1.4.20')
    zip_ref.close()


def pip():
    import os
    import urllib.request
    import warnings

    print('\u001b[36;1m Downloading pip now. \u001b[0m')
    # Getting get-pip module
    urllib.request.urlretrieve('https://bootstrap.pypa.io/get-pip.py', 'get-pip.py')
    with warnings.simplefilter("ignore"):
        try:
            os.system('python get-pip.py')
        except:
            # Run it and install pip into user directory.
            os.system('python get-pip.py --user')
            os.system('PATH=$PATH:~/.local/bin')


def snpflip():
    try:
        import pip
        # Use pip to install snpflip and it's dependencies.
        pip.main(['install', 'snpflip'])
    except ImportError:
        pip()
        import pip
        pip.main(['install', 'snpflip'])


def shapeit():
    import platform
    import urllib.request
    import os
    import sys

    # Since shapeit only works on linux or mac, we need to first check what system they are on.
    system_check = platform.system()

    if system_check == "Linux":
        print('\u001b[36;1m Downloading shapeit now. \u001b[0m')
        urllib.request.urlretrieve(
            'https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz',
            'shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz')
        # Making directory to store program
        os.makedirs('Shapeit_v2.12_Linux_Static')
        # Unpacking
        os.system('tar -zxvf shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz -C Shapeit_v2.12_Linux_Static/')

    # If the user is on a mac
    elif system_check == "Darwin":
        print('\u001b[35;1m Downloading shapeit now. \u001b[0m')
        # Download shapeit
        urllib.request.urlretrieve(
            'https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.MacOSX.tgz',
            'shapeit.v2.r837.MacOSX.tgz')
        # Create directory for shapeit.
        os.makedirs('Shapeit_v2.20_Mac')
        # Untar shapeit to that directory.
        os.system('tar -zxvf shapeit.v2.r837.MacOSX.tgz -C Shapeit_v2.20_Mac/')

    # If they are running this on a windows machine, they cannot proceed because shapeit is *nix only.
    elif system_check == "Windows":
        sys.exit("\u001b[35;1m I'm sorry, I've detected that you're working on a Windows computer and shapeit is a "
                 "linux or unix program only. . If you have access to the Penn State clusters, you should run this "
                 "script from there (they are linux). \u001b[0m")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("\u001b[35;1m I cannot detect the system you are working on. Exiting now. \u001b[0m")


def vcftools():
    import platform
    import urllib.request
    import os
    import sys

    from os.path import expanduser
    home = expanduser("~")

    # vcftools is easily installed on a linux system, more difficult to install on a mac, and does not have a Windows
    # distribution, so we need to check the platform.
    system_check = platform.system()

    if system_check == "Linux":
        print('\u001b[36;1m Downloading vcftools now. \u001b[0m')
        if not os.path.exists(os.path.join(home,'software')):
            os.makedirs(os.path.join(home,"software"))
        urllib.request.urlretrieve('https://github.com/vcftools/vcftools/tarball/master',
                                   os.path.join(home,'software/vcftools.tgz'))
        # Unpacking
        os.system('tar -xvf ' + os.path.join(home, '/software/vcftools.tgz') + ' -C ' + os.path.join(home, 'software'))
        # Moving into the vcftools folder
        os.system('cd ' + os.path.join(home, 'software/vcftools'))
        # Running configuration and installation steps
        os.system('./autogen.sh')
        os.system('./configure --prefix=' + os.path.join(home,'software'))
        os.system('make')
        os.system('make install')
        os.system('export PATH=$PATH:' + os.path.join(home,'software/bin'))

    # If the user is on a mac
    elif system_check == "Darwin":
        sys.exit('\u001b[35;1m To download vcftools on a Mac you should get homebrew. On your own, download homebrew from '
              'https://brew.sh/ then download vcftools by typing "brew install homebrew/science/vcftools" in your '
              'command terminal. Then, follow the directions in the README file for configuring. \u001b[0m')

    # If they are running this on a windows machine, they cannot proceed because shapeit is *nix only.
    elif system_check == "Windows":
        sys.exit("\u001b[35;1m I'm sorry, I've detected that you're working on a Windows computer and vcftools is a "
                 "linux or unix program only. . If you have access to the Penn State clusters, you should run this "
                 "script from there (they are linux). \u001b[0m")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("\u001b[35;1m I cannot detect the system you are working on. Exiting now. \u001b[0m")


def bcftools():
    import platform
    import urllib.request
    import os
    import sys

    from os.path import expanduser
    home = expanduser("~")

    # bcftools is easily installed on a linux system, more difficult to install on a mac, and does not have a Windows
    # distribution, so we need to check the platform.
    system_check = platform.system()

    if system_check == "Linux":
        print('\u001b[36;1m Downloading bcftools now. \u001b[0m')
        if not os.path.exists(os.path.join(home,'software')):
            os.mkdir(os.path.join(home,'software'))
        urllib.request.urlretrieve('https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2',
                                   os.path.join(home,'software/bcftools-1.6.tar.bz2'))
        # Unpacking
        os.system('tar -xvjf ' + os.path.join(home,'software/bcftools-1.6.tar.bz2') + ' -C '
                  + os.path.join(home,'software'))
        # Moving into the bcftools folder
        os.chdir(os.path.join(home,'software/bcftools-1.6/'))
        # Running configuration and installation steps
        os.system('./configure --prefix=' + os.path.join(home,'software'))
        os.system('make')
        os.system('make install')
        os.system('export PATH=$PATH:' + os.path.join(home, 'software/bin'))

    # If the user is on a mac
    elif system_check == "Darwin":
        sys.exit("\u001b[35;1m It's possible to install bcftools on a Mac, but it might require that you install other "
            "libraries too. You can try following this tutorial "
            "http://www.danielecook.com/installing-tabix-and-samtools-on-mac/ and download homebrew (if you haven't "
            "already) and xcode. Then run brew install homebrew/science/bcftools. \u001b[0m")

    # If they are running this on a windows machine, they cannot proceed because bcftools is *nix only.
    elif system_check == "Windows":
        sys.exit("\u001b[35;1m I'm sorry, I've detected that you're working on a Windows computer and bcftools is a "
                 "linux or unix program only. If you have access to the Penn State clusters, you should run this "
                 "script from there (they are linux). \u001b[0m")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("\u001b[35;1m I cannot detect the system you are working on. Exiting now. \u001b[0m")


def htslib():
    import platform
    import urllib.request
    import os
    import sys

    from os.path import expanduser
    home = expanduser("~")

    # htslib is easily installed on a linux system, more difficult to install on a mac, and does not have a Windows
    # distribution, so we need to check the platform.
    system_check = platform.system()

    if system_check == "Linux":
        print('\u001b[36;1m Downloading htslib now. \u001b[0m')
        if not os.path.exists(os.path.join(home,'software')):
            os.mkdir(os.path.join(home,'software'))
        urllib.request.urlretrieve('https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2',
                                   os.path.join(home,'software,htslib-1.6.tar.bz2'))
        # Making directory to store program
        # os.makedirs(os.path.join(home,'software/htslib_1.6'))
        # Unpacking
        os.system('tar -xvjf ' + os.path.join(home,'software/htslib-1.6.tar.bz2') + ' -C '
                  + os.path.join(home,'software'))
        # Moving into the htslib folder
        os.system('cd ' + os.path.join(home, 'software/htslib-1.6/'))
        # Running configuration and installation steps
        os.system('./configure --prefix=' + os.path.join(home, 'software'))
        os.system('make')
        os.system('make install')
        os.system('export PATH=$PATH:' + os.path.join(home,'software/bin'))

    # If the user is on a mac
    elif system_check == "Darwin":
        sys.exit("\u001b[35;1m It's possible to install htslib on a Mac, but it might require that you install other "
            "libraries too. You can try following this tutorial "
            "http://www.danielecook.com/installing-tabix-and-samtools-on-mac/ and download homebrew (if you haven't "
            "already) and xcode. Then run brew install homebrew/science/htslib. \u001b[0m")

    # If they are running this on a windows machine, they cannot proceed because bcftools is *nix only.
    elif system_check == "Windows":
        sys.exit(
            "\u001b[35;1m I'm sorry, I've detected that you're working on a Windows computer and htslib is a "
            "linux or unix program only. If you have access to the Penn State clusters, you should run this "
            "script from there (they are linux). \u001b[0m")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("\u001b[35;1m I cannot detect the system you are working on. Exiting now. \u001b[0m")
