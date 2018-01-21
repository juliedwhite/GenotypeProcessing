import platform
import subprocess
import urllib.request
import os
import sys
from os.path import expanduser

system_check = platform.system()
home = expanduser("~")


# Define get pip
def pip():
    print('Downloading pip now. You need this for downloading other things.')
    # Getting get-pip module
    urllib.request.urlretrieve('https://bootstrap.pypa.io/get-pip.py', 'get-pip.py')
    try:
        subprocess.check_output(['python', 'get-pip.py'])
        print("Done downloading pip")
    except:
        # Run it and install pip into user directory.
        subprocess.check_output(['python', 'get-pip.py', '--user'])
        os.remove('get-pip.py')
        # Tell the user that they should add the following folders to their PATH before rerunning the script:
        if system_check in ("Linux", "Darwin"):
            sys.exit("Because I couldn't use root access, I installed pip into a local directory. You should type "
                     "'echo $PATH' to check that " + os.path.join(home, '.local/bin')
                     + " is in your $PATH variable. If it isn't, then you should add it to your .bash_profile. If you "
                       "don't know how to do this, ask the internet or myself. Then re-run this script.")

        elif system_check == "Windows":
            sys.exit("Because I couldn't use root access, I installed pip into a local directory. You should type "
                     "'PATH' to check that " + os.path.join(home, '.local/bin')
                     + " is in your PATH variable. If it isn't, then you should add it to your environment variables. "
                       "If you won't know how to do this, ask the internet or myself. Then re-run this script.")


# Define module to get colorama
def getcolorama():
    try:
        import pip
    except ImportError:
        pip()
        import pip

    # Try to download
    try:
        pip.main(['install', 'colorama'])
    except:
        pip.main(['install', 'colorama', '--user'])

# Try to import colorama, if that doesn't work, download it using pip.
try:
    import colorama
    from colorama import init, Fore, Style
except ImportError:
    getcolorama()
    import colorama
    from colorama import init, Fore, Style
    init()


# What to download?
def todownload():
    print(Fore.BLUE + Style.BRIGHT)
    item = input('What would you like to download?\n'
                 '1) Everything - this will take about an hour and you should be careful to read the messages that '
                 'appear on the screen. Some installations require your help to complete.\n'
                 '2) Plink 1.9\n'
                 '3) 1000G Phase 3 VCF\n'
                 '4) 1000G Phase 3 Hap/Legend/Sample\n'
                 '5) GRCh37/hg19 1000G FASTA file\n'
                 '6) Genotype Harmonizer\n'
                 '7) snpflip\n'
                 '8) shapeit\n'
                 '9) htslib\n'
                 '10) vcftools\n'
                 '11) bcftools\n'
                 '12) samtools\n'
                 '13) Nothing\n'
                 'Please enter the number of the reference file or program that you would like to download (i.e. 2): ')
    print(Style.RESET_ALL)

    if item == '1':
        plink()
        vcf_1000g_phase3()
        hls_1000g_phase3()
        fasta_1000G_hg19()
        genotype_harmonizer()
        snpflip()
        shapeit()
        htslib()
        vcftools()
        bcftools()
        samtools()
        if system_check == "Linux":
            print(Fore.RED + Style.BRIGHT + "I installed several programs into " + os.path.join(home, 'software/bin')
                  + ". You should type 'echo $PATH' to check that " + os.path.join(home, 'software/bin')
                  + " is in your $PATH variable. If it isn't, I highly recommend you add this path to your "
                    ".bash_profile. If you don't know how to do this, ask the internet or myself.")
            print(Style.RESET_ALL)
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
        snpflip()
    elif item == '8':
        shapeit()
    elif item == '9':
        htslib()
        if system_check == "Linux":
            print(Fore.RED + Style.BRIGHT + "I installed htslib into " + os.path.join(home, 'software/bin')
                  + ". You should type 'echo $PATH' to check that " + os.path.join(home, 'software/bin')
                  + " is in your $PATH variable. If it isn't, I highly recommend you add this path to your "
                    ".bash_profile. If you don't know how to do this, ask the internet or myself.")
            print(Style.RESET_ALL)
    elif item == '10':
        vcftools()
        if system_check == "Linux":
            print(Fore.RED + Style.BRIGHT + "I installed vcftools into " + os.path.join(home, 'software/bin')
                  + ". You should type 'echo $PATH' to check that " + os.path.join(home, 'software/bin')
                  + " is in your $PATH variable. If it isn't, I highly recommend you add this path to your "
                    ".bash_profile. If you don't know how to do this, ask the internet or myself.")
            print(Style.RESET_ALL)
    elif item == '11':
        bcftools()
        if system_check == "Linux":
            print(Fore.RED + Style.BRIGHT + "I installed bcftools into " + os.path.join(home, 'software/bin')
                  + ". You should type 'echo $PATH' to check that " + os.path.join(home, 'software/bin')
                  + " is in your $PATH variable. If it isn't, I highly recommend you add this path to your "
                    ".bash_profile. If you don't know how to do this, ask the internet or myself.")
            print(Style.RESET_ALL)
    elif item == '12':
        samtools()
        if system_check == "Linux":
            print(Fore.RED + Style.BRIGHT + "I installed samtools into " + os.path.join(home, 'software/bin')
                  + ". You should type 'echo $PATH' to check that " + os.path.join(home, 'software/bin')
                  + " is in your $PATH variable. If it isn't, I highly recommend you add this path to your "
                    ".bash_profile. If you don't know how to do this, ask the internet or myself.")
            print(Style.RESET_ALL)
    elif item == '13':
        sys.exit("Exiting now")
    else:
        sys.exit("Quitting because you did not give a recognizable number when asked what to download.")


def plink():
    import zipfile
    import shutil

    try:
        import pip
    except ImportError:
        pip()
        import pip

    try:
        import requests
    except (ImportError, ModuleNotFoundError):
        try:
            pip.main(['install', 'requests'])
        except:
            pip.main(['install', 'requests', '--user'])
        import requests

    try:
        import lxml.html
    except (ImportError, ModuleNotFoundError):
        try:
            pip.main(['install', 'lxml'])
        except:
            pip.main(['install', 'lxml', '--user'])
        import lxml.html

    try:
        import cssselect
    except (ImportError, ModuleNotFoundError):
        try:
            pip.main(['install', 'cssselect'])
        except:
            pip.main(['install', 'cssselect', '--user'])
        import cssselect

    # Get the version of the user's system
    architecture_check = platform.architecture()[0]

    # Get plink urls
    start_url = 'https://www.cog-genomics.org/plink/1.9/'
    response = requests.get(start_url)
    tree = lxml.html.fromstring(response.text)
    links = tree.cssselect('a')  # or tree.xpath('//a')
    out = []
    for link in links:
        # we use this if just in case some <a> tags lack an href attribute
        if 'href' in link.attrib:
            out.append(link.attrib['href'])

    linux_x86 = [x for x in out if '/plink_linux_x86_64.zip' in x]
    linux_x32 = [x for x in out if '/plink_linux_i686.zip' in x]
    mac = [x for x in out if '/plink_mac.zip' in x]
    win64 = [x for x in out if '/plink_win64.zip' in x]
    win32 = [x for x in out if '/plink_win32.zip' in x]

    base_url = 'https://www.cog-genomics.org'

    # Download based on system
    if system_check == "Linux":
        if architecture_check == '64bit':
            # Download 64 bit Linux Plink 1.9
            print("Downloading Linux Plink 1.9 to this directory now.")
            urllib.request.urlretrieve(requests.compat.urljoin(base_url, linux_x86[0]), 'Plink_1.9_Linux64.zip')
            # Making directory to store program
            os.makedirs('Plink_1.9_Linux64')
            # Unpacking to this directory
            with zipfile.ZipFile("Plink_1.9_Linux64.zip", "r") as zip_ref:
                zip_ref.extractall("Plink_1.9_Linux64")
            # Copy plink from archive folder to current working directory.
            shutil.copy2('Plink_1.9_Linux64/plink', os.getcwd())

        elif architecture_check == '32bit':
            # Download 32 bit Linux Plink 1.9
            print("Downloading Linux Plink 1.9 to this directory now.")
            urllib.request.urlretrieve(requests.compat.urljoin(base_url, linux_x32[0]), 'Plink_1.9_Linux32.zip')
            # Making directory to store program
            os.makedirs('Plink_1.9_Linux32')
            # Unpacking to this directory
            with zipfile.ZipFile("Plink_1.9_Linux32.zip", "r") as zip_ref:
                zip_ref.extractall("Plink_1.9_Linux32")
            # Copy plink from archive folder to current working directory.
            shutil.copy2('Plink_1.9_Linux32/plink', os.getcwd())

        else:
            print("I'm sorry, I could not determine what Linux Plink version to download.")

        if os.path.exists('plink'):
            subprocess.check_output(['chmod', '777', 'plink'])

    elif system_check == "Darwin":
        # Download Mac Plink 1.9
        print("Downloading Mac Plink 1.9 to this directory now.")
        urllib.request.urlretrieve(requests.compat.urljoin(base_url, mac[0]), 'Plink_1.9_Mac.zip')
        # Making directory to store program
        os.makedirs('Plink_1.9_Mac')
        # Unpacking to this directory
        with zipfile.ZipFile("Plink_1.9_Mac.zip", "r") as zip_ref:
            zip_ref.extractall("Plink_1.9_Mac")
        # Copy plink from archive folder to current working directory.
        shutil.copy2('Plink_1.9_Mac/plink', os.getcwd())

    elif system_check == "Windows":
        if architecture_check == '64bit':
            # Download 64 bit Windows Plink 1.9
            print("Downloading Windows Plink 1.9 to this directory now.")
            urllib.request.urlretrieve(requests.compat.urljoin(base_url, win64[0]), 'Plink_1.9_Win64.zip')
            # Making directory to store program
            os.makedirs('Plink_1.9_Win64')
            # Unpacking to this directory
            with zipfile.ZipFile("Plink_1.9_Win64.zip", "r") as zip_ref:
                zip_ref.extractall("Plink_1.9_Win64")
            # Copy plink from archive folder to current working directory.
            shutil.copy2('Plink_1.9_Win64/plink.exe', os.getcwd())

        elif architecture_check == '32bit':
            # Download 32 bit Windows Plink 1.9
            print("Downloading Windows Plink 1.9 to this directory now.")
            urllib.request.urlretrieve(requests.compat.urljoin(base_url, win32[0]), 'Plink_1.9_Win32.zip')
            # Making directory to store program
            os.makedirs('Plink_1.9_Win32')
            # Unpacking to this directory
            with zipfile.ZipFile("Plink_1.9_Win32.zip", "r") as zip_ref:
                zip_ref.extractall("Plink_1.9_Win32")
            # Copy plink from archive folder to current working directory.
            shutil.copy2('Plink_1.9_Win32/plink.exe', os.getcwd())

        else:
            print(Fore.RED + Style.BRIGHT + "I'm sorry, I could not determine what Windows Plink version to download.")
            print(Style.RESET_ALL)
    else:
        print(Fore.RED + Style.BRIGHT + "I'm sorry, I could not determine what operating system you are running. "
                                        "Please download the latest version of Plink at "
                                        "https://www.cog-genomics.org/plink2")
        print(Style.RESET_ALL)

    print("Done downloading plink")


def vcf_1000g_phase3():
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
    tbi_file_names = [s + '.tbi' for s in vcf_file_names]

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
    for filename in tbi_file_names:
        local_filename = os.path.join(os.getcwd(), '1000G_Phase3_VCF', filename)
        file = open(local_filename, 'wb')
        ftp.retrbinary('RETR ' + filename, file.write)
        file.close()
    # Once we are done downloading, close the connection to the ftp server.
    ftp.quit()

    print("Done downloading 1000G Phase3 VCF files")


def hls_1000g_phase3():
    print('Downloading 1000G Phase 3 files now, putting them in "1000G_Phase3_HapLegendSample" folder. '
          'This will create a ~12G folder on your computer and will take a little while.')

    orig_wd = os.getcwd()

    # Create folder:
    if not os.path.exists('1000G_Phase3_HapLegendSample'):
        os.makedirs('1000G_Phase3_HapLegendSample')

        # Where the files are
        legend_server = "http://mathgen.stats.ox.ac.uk/impute/"

        # List of file names that we're going to need.
        file_list = ['1000GP_Phase3.tgz', '1000GP_Phase3_chrX.tgz']

        # Download files and put them in 1000G_Phase3_VCF folder.
        for filename in file_list:
            urllib.request.urlretrieve(os.path.join(legend_server, filename),
                                       os.path.join(os.getcwd(), '1000G_Phase3_HapLegendSample', filename))

        # Change directory where tgz files are
        os.chdir('1000G_Phase3_HapLegendSample')
        # Unpack tar files
        if system_check in ("Linux", "Darwin"):
            subprocess.check_output(['tar', '-xvf', '1000GP_Phase3.tgz', '--strip-components', '1'])
            subprocess.check_output(['tar', '-xvf', '1000GP_Phase3_chrX.tgz'])
        elif system_check == "Windows":
            print(Fore.RED + Style.BRIGHT + "Please use 7-Zip to unzip your .tgz file in the "
                                            "1000G_Phase3_HapLegendSample folder, then use 7-Zip again to untar the "
                                            ".tar file. Do this for both 1000GP_Phase3.tgz and 1000GP_Phase3_chrX.tgz")
            print(Style.RESET_ALL)
        # Print when done.
        print('Done extracting 1000G_Phase3_HapLegendSample')

        # Change back to original working directory.
        os.chdir(orig_wd)

    print("Done downloading 1000G Phase3 Hap/Legend/Sample files")


def fasta_1000G_hg19():
    import ftplib
    import gzip
    import shutil

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

    # Fasta file needs to be unzipped for snpflip to work.
    try:
        with gzip.open(os.path.join('1000G_hg19_fasta', 'human_g1k_v37.fasta.gz'), 'rb') as f_in, \
                open(os.path.join('1000G_hg19_fasta', 'human_g1k_v37.fasta'), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    except:
        if system_check in ("Linux", "Darwin"):
            subprocess.call(['gunzip', '-c', os.path.join('1000G_hg19_fasta', 'human_g1k_v37.fasta.gz'), '>',
                             os.path.join('1000G_hg19_fasta', 'human_g1k_v37.fasta')])
        elif system_check == "Windows":
            zip_path = []
            for r, d, f in os.walk(os.path.join('C:\\', 'Program Files')):
                for files in f:
                    if files == "7zG.exe":
                        zip_path = os.path.join(r, files)
            subprocess.check_output([zip_path, 'e', os.path.join('1000G_hg19_fasta', 'human_gik_v37.fasta.gz')])

    print("Done downloading hg19 fasta file")


def genotype_harmonizer():
    import zipfile
    from distutils.dir_util import copy_tree

    print('Downloading genotype harmonizer now.')
    # Download genotype harmonizer zip file
    urllib.request.urlretrieve(
        'http://www.molgenis.org/downloads/GenotypeHarmonizer/GenotypeHarmonizer-1.4.20-dist.zip',
        'GenotypeHarmonizer-1.4.20.zip')
    # Unzip Genotype Harmonizer
    zip_ref = zipfile.ZipFile('GenotypeHarmonizer-1.4.20.zip', 'r')
    zip_ref.extractall('GenotypeHarmonizer-1.4.20')
    zip_ref.close()

    # copy subdirectory up
    fromdirectory = os.path.join('GenotypeHarmonizer-1.4.20', 'GenotypeHarmonizer_1.4.20-SNAPSHOT')
    todirectory = 'GenotypeHarmonizer-1.4.20'
    copy_tree(fromdirectory, todirectory)

    print("Done downloading genotype harmonizer")


def snpflip():
    try:
        import pip
    except ImportError:
        pip()
        import pip
    # Use pip to install snpflip and it's dependencies
    try:
        pip.main(['install', 'snpflip'])
    except:
        pip.main(['install', 'snpflip', '--user'])
    print("Done installing snpflip")


def shapeit():
    if system_check == "Linux":
        print('Downloading shapeit now.')
        urllib.request.urlretrieve(
            'https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r900.glibcv2.12.linux.tar.gz',
            'shapeit.v2.r900.glibcv2.12.linux.tar.gz')
        # Making directory to store program
        os.makedirs('Shapeit_v2.r900_Linux_Static')
        # Unpacking
        subprocess.check_output(['tar', '-xzvf', 'shapeit.v2.r900.glibcv2.12.linux.tar.gz', '-C',
                                 'Shapeit_v2.r900_Linux_Static/'])

        print("Done downloading shapeit")

    # If the user is on a mac
    elif system_check == "Darwin":
        print('Downloading shapeit now.')
        # Download shapeit
        urllib.request.urlretrieve(
            'https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.MacOSX.tgz',
            'shapeit.v2.r837.MacOSX.tgz')
        # Create directory for shapeit.
        os.makedirs('Shapeit_v2.20_Mac')
        # Untar shapeit to that directory.
        subprocess.check_output(['tar', '-zxvf', 'shapeit.v2.r837.MacOSX.tgz', '-C', 'Shapeit_v2.20_Mac/'])

        print("Done downloading shapeit")

    # If they are running this on a windows machine, they cannot proceed because shapeit is *nix only.
    elif system_check == "Windows":
        print(Fore.RED + Style.BRIGHT + "I'm sorry, I've detected that you're working on a Windows computer and "
                                        "shapeit is a linux or unix program only. If you have access to the Penn State "
                                        "clusters, you should run this script from there (they are linux).")

    # If I cannot detect what system they're on, force exit.
    else:
        print(Fore.RED + Style.BRIGHT + "I cannot detect the system you are working on. Please download Shapeit on "
                                        "your own.")
        print(Style.RESET_ALL)


def vcftools():
    if system_check == "Linux":
        print('Downloading vcftools now.')
        if not os.path.exists(os.path.join(home, 'software')):
            os.makedirs(os.path.join(home, "software"))
        urllib.request.urlretrieve('https://github.com/vcftools/vcftools/tarball/master',
                                   os.path.join(home, 'software/vcftools.tgz'))
        # Unpacking
        subprocess.check_output(['tar', '-xvf', os.path.join(home, 'software/vcftools.tgz'), '-C',
                                 os.path.join(home, 'software')])
        # Moving into the vcftools folder
        os.chdir(os.path.join(home, 'software/vcftools-vcftools-ea875e2'))
        # Running configuration and installation steps
        subprocess.check_output('./autogen.sh')
        subprocess.check_output(['./configure', '--prefix=' + os.path.join(home, 'software')])
        subprocess.check_output('make')
        subprocess.check_output(['make', 'install'])
        # Tell the user that they should add the following folders to their PATH:
        print(Fore.RED + Style.BRIGHT + "For vcftools to work properly, you should set your PERL5LIB. I highly "
                                        "recommend that you add the following line to your .bash_profile: "
                                        "'export PERL5LIB="
              + os.path.join(home, 'software/vcftools-vcftools-ea875e2/src/perl'))
        print(Style.RESET_ALL)

        print("Done downloading vcftools")

    # If the user is on a mac
    elif system_check == "Darwin":
        print('To download vcftools on a Mac you should get homebrew. On your own, download homebrew from '
              'https://brew.sh/ then download vcftools by typing "brew install homebrew/science/vcftools" in your '
              'command terminal. Then, follow the directions in the README file for configuring. Then make sure the'
              'install folder is in your PATH variable.')

    # If they are running this on a windows machine, they cannot proceed because shapeit is *nix only.
    elif system_check == "Windows":
        print("I'm sorry, I've detected that you're working on a Windows computer and vcftools is a linux or unix "
              "program only. If you have access to the Penn State clusters, you should run this script from there "
              "(they are linux).")

    # If I cannot detect what system they're on, force exit.
    else:
        print(Fore.RED + Style.BRIGHT + "I cannot detect the system you are working on. Please download vcftools on "
                                        "your own.")
        print(Style.RESET_ALL)


def bcftools():
    if system_check == "Linux":
        print('Downloading bcftools now.')
        if not os.path.exists(os.path.join(home, 'software')):
            os.mkdir(os.path.join(home, 'software'))
        urllib.request.urlretrieve('https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2',
                                   os.path.join(home, 'software/bcftools-1.6.tar.bz2'))
        # Unpacking
        subprocess.check_output(['tar', '-xvjf', os.path.join(home, 'software/bcftools-1.6.tar.bz2'), '-C',
                                 os.path.join(home, 'software')])
        # Moving into the bcftools folder
        os.chdir(os.path.join(home, 'software/bcftools-1.6/'))
        # Running configuration and installation steps
        subprocess.check_output(['./configure', '--prefix=' + os.path.join(home, 'software')])
        subprocess.check_output('make')
        subprocess.check_output(['make', 'install'])

        print("Done downloading bcftools")

    # If the user is on a mac
    elif system_check == "Darwin":
        print("It's possible to install bcftools on a Mac, but it might require that you install other libraries too. "
              "You can try following this tutorial (not for bcftools specifically, but related programs) "
              "http://www.danielecook.com/installing-tabix-and-samtools-on-mac/ and download homebrew (if you haven't "
              "already) and xcode. Then run brew install homebrew/science/bcftools and install using the directions in "
              "the README file. You should also add the install location to your PATH if it is not already there.")

    # If they are running this on a windows machine, they cannot proceed because bcftools is *nix only.
    elif system_check == "Windows":
        print("I'm sorry, I've detected that you're working on a Windows computer and bcftools is a linux or unix "
              "program only. If you have access to the Penn State clusters, you should run this script from there "
              "(they are linux).")

    # If I cannot detect what system they're on, force exit.
    else:
        print(Fore.RED + Style.BRIGHT + "I cannot detect the system you are working on. Please download bcftools on "
                                        "your own.")
        print(Style.RESET_ALL)


def htslib():
    if system_check == "Linux":
        print('Downloading htslib now.')
        if not os.path.exists(os.path.join(home, 'software')):
            os.mkdir(os.path.join(home, 'software'))
        urllib.request.urlretrieve('https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2',
                                   os.path.join(home, 'software/htslib-1.6.tar.bz2'))
        # Making directory to store program
        # os.makedirs(os.path.join(home,'software/htslib_1.6'))
        # Unpacking
        subprocess.check_output(['tar', '-xvjf', os.path.join(home, 'software/htslib-1.6.tar.bz2'), '-C',
                                 os.path.join(home, 'software')])
        # Moving into the htslib folder
        os.chdir(os.path.join(home, 'software/htslib-1.6/'))
        # Running configuration and installation steps
        subprocess.check_output(['./configure', '--prefix=' + os.path.join(home, 'software')])
        subprocess.check_output('make')
        subprocess.check_output(['make', 'install'])
        print("Done downloading htslib")
    # If the user is on a mac
    elif system_check == "Darwin":
        print("It's possible to install htslib on a Mac, but it might require that you install other libraries too. "
              "You can try following this tutorial (not for htslib specifically, but related programs) "
              "http://www.danielecook.com/installing-tabix-and-samtools-on-mac/ and download homebrew (if you haven't "
              "already) and xcode. Then run brew install homebrew/science/htslib and install using the directions in "
              "the README file. You should also add the install location to your PATH if it is not already there.")
    # If they are running this on a windows machine, they cannot proceed because bcftools is *nix only.
    elif system_check == "Windows":
        print("I'm sorry, I've detected that you're working on a Windows computer and htslib is a linux or unix "
              "program only. If you have access to the Penn State clusters, you should run this script from there "
              "(they are linux).")

    # If I cannot detect what system they're on, force exit.
    else:
        print(Fore.RED + Style.BRIGHT + "I cannot detect the system you are working on. Please download htslib on "
                                        "your own.")
        print(Style.RESET_ALL)


def samtools():
    if system_check == "Linux":
        print('Downloading samtools now.')
        if not os.path.exists(os.path.join(home, 'software')):
            os.mkdir(os.path.join(home, 'software'))
        urllib.request.urlretrieve('https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2',
                                   os.path.join(home, 'software/samtools-1.6.tar.bz2'))

        # Unpacking
        subprocess.check_output(['tar', '-xvjf', os.path.join(home, 'software/samtools-1.6.tar.bz2'), '-C',
                                 os.path.join(home, 'software')])
        # Moving into the samtools folder
        os.chdir(os.path.join(home, 'software/samtools-1.6/'))
        # Running configuration and installation steps
        subprocess.check_output(['./configure', '--prefix=' + os.path.join(home, 'software')])
        subprocess.check_output('make')
        subprocess.check_output(['make', 'install'])

        print("Done downloading samtools")

    # If the user is on a mac
    elif system_check == "Darwin":
        print("It's possible to install samtools on a Mac, but it might require that you install other libraries too. "
              "You can try following this tutorial http://www.danielecook.com/installing-tabix-and-samtools-on-mac/ "
              "and download homebrew (if you haven't already) and xcode. You should also add the install location to "
              "your PATH if it is not already there.")

    # If they are running this on a windows machine, they cannot proceed because samtools is *nix only.
    elif system_check == "Windows":
        print("I'm sorry, I've detected that you're working on a Windows computer and samtools is a linux or unix "
              "program only. If you have access to the Penn State clusters, you should run this script from there "
              "(they are linux).")

    # If I cannot detect what system they're on, force exit.
    else:
        print(Fore.RED + Style.BRIGHT + "I cannot detect the system you are working on. Please download samtools on "
                                        "your own.")
        print(Style.RESET_ALL)


def getmatplotlib():
    try:
        import pip
    except ImportError:
        pip()
        import pip

    try:
        pip.main(['install', 'matplotlib'])
    except ImportError:
        pip.main(['install', 'matplotlib', '--user'])

    print("Done installing matplotlib")


def getargparse():
    try:
        import pip
    except ImportError:
        pip()
        import pip

    # Try to download
    try:
        pip.main(['install', 'argparse'])
    except:
        pip.main(['install', 'argparse', '--user'])


def getpandas():
    try:
        import pip
    except ImportError:
        pip()
        import pip

    # Try to download
    try:
        pip.main(['install', 'pandas'])
    except:
        pip.main(['install', 'pandas', '--user'])


def getnumpy():
    try:
        import pip
    except ImportError:
        pip()
        import pip

    # Try to download
    try:
        pip.main(['install', 'numpy'])
    except:
        pip.main(['install', 'numpy', '--user'])