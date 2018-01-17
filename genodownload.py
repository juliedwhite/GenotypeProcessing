def pip():
    import urllib.request
    import platform
    import subprocess

    print('Downloading pip now')
    # Getting get-pip module
    urllib.request.urlretrieve('https://bootstrap.pypa.io/get-pip.py', 'get-pip.py')
    try:
        subprocess.check_output(['python','get-pip.py'])
    except:
        # Run it and install pip into user directory.
        subprocess.check_output(['python','get-pip.py','--user'])

        # Tell the user that they should add the following folders to their PATH:
        system_check = platform.system()
        if system_check in ("Linux", "Darwin"):
            print(Fore.RED + Style.BRIGHT)
            print("Because I couldn't use root access, I installed pip into a local directory. You should type "
                  "'echo $PATH' to check that '~/.local/bin' is in your $PATH variable. IF it isn't, then type you "
                  "should add it to your .bash_profile. If you don't know how to do this, ask the internet or myself.")
            print(Style.RESET_ALL)

        elif system_check == "Windows":
            print(Fore.RED + Style.BRIGHT)
            print("Because I couldn't use root access, I installed pip into a local directory. You should type 'PATH' "
                  "to check that '~/.local/bin' is in your PATH variable. IF it isn't, then you should add it to your "
                  "environment variables. If you won't know how to do this, ask the internet or myself.")
            print(Style.RESET_ALL)

    print("Done downloading pip")


def getcolorama():
    try:
        import pip
    except ImportError:
        pip()
        import pip

    try:
        pip.main(['install', 'colorama'])
    except:
        pip.main(['install', 'colorama', '--user'])


try:
    import colorama
except ImportError:
    getcolorama()


from colorama import init, Fore, Style
init()

def todownload():
    import sys
    print(Fore.BLUE + Style.BRIGHT)
    item = input('What would you like to download?\n'
                 '1) Plink 1.9\n'
                 '2) 1000G Phase 3 VCF\n'
                 '3) 1000G Phase 3 Hap/Legend/Sample\n'
                 '4) GRCh37/hg19 1000G FASTA file\n'
                 '5) Genotype Harmonizer\n'
                 '6) pip\n'
                 '7) snpflip\n'
                 '8) shapeit\n'
                 '9) vcftools\n'
                 '10) bcftools\n'
                 '11) htslib\n'
                 '12) samtools\n'
                 '13) Nothing\n'
                 'Please enter a number (i.e. 2): ')
    print(Style.RESET_ALL)

    if item == '1':
        plink()
    elif item == '2':
        vcf_1000g_phase3()
    elif item == '3':
        hls_1000g_phase3()
    elif item == '4':
        fasta_1000G_hg19()
    elif item == '5':
        genotype_harmonizer()
    elif item == '6':
        pip()
    elif item == '7':
        snpflip()
    elif item == '8':
        shapeit()
    elif item == '9':
        vcftools()
    elif item == '10':
        bcftools()
    elif item == '11':
        htslib()
    elif item == '12':
        samtools()
    elif item == '13':
        sys.exit("Exiting now")
    else:
        sys.exit("Quitting because you did not give a recognizable number when asked what to download.")


def plink():
    import os
    import platform
    import zipfile
    import shutil
    import urllib.request
    import sys
    import subprocess

    try:
        import pip
    except ImportError:
        pip()
        import pip

    try:
        import requests
    except ImportError:
        pip.main(['install', 'requests'])
        import requests

    try:
        import lxml.html
    except ImportError:
        pip.main(['install', 'lxml'])
        import lxml.html

    try:
        import cssselect
    except ImportError:
        pip.main(['install', 'cssselect'])
        import cssselect

    # Get what system the user is using
    system_check = platform.system()
    # Get the version of that system
    architecture_check = platform.architecture()[0]

    #Get plink urls
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
            sys.exit("I'm sorry, I could not determine what Linux Plink version to download.")

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
            sys.exit("I'm sorry, I could not determine what Windows Plink version to download.")

    else:
        sys.exit("I'm sorry, I could not determine what operating system you are running. Please download the latest "
                 "version of Plink at https://www.cog-genomics.org/plink2")

    print("Done downloading plink")


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
    import os
    import urllib.request
    import subprocess

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
        subprocess.check_output(['tar', '-xvf', '1000GP_Phase3.tgz', '--strip-components', '1'])
        subprocess.check_output(['tar', '-xvf', '1000GP_Phase3_chrX.tgz'])

        # Print when done.
        print('Done extracting 1000G_Phase3_HapLegendSample')

        # Change back to original working directory.
        os.chdir(orig_wd)

    print("Done downloading 1000G Phase3 Hap/Legend/Sample files")


def fasta_1000G_hg19():
    import os
    import ftplib

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

    print("Done downloading hg19 fasta file")


def genotype_harmonizer():
    import urllib.request
    import zipfile

    print('Downloading genotype harmonizer now.')
    # Download genotype harmonizer zip file
    urllib.request.urlretrieve(
        'http://www.molgenis.org/downloads/GenotypeHarmonizer/GenotypeHarmonizer-1.4.20-dist.zip',
        'GenotypeHarmonizer-1.4.20.zip')
    # Unzip Genotype Harmonizer
    zip_ref = zipfile.ZipFile('GenotypeHarmonizer-1.4.20.zip', 'r')
    zip_ref.extractall('GenotypeHarmonizer-1.4.20')
    zip_ref.close()

    print("Done downloading genotype harmonizer")



def snpflip():
    try:
        import pip
        # Use pip to install snpflip and it's dependencies.
        pip.main(['install', 'snpflip'])
    except ImportError:
        pip()
        import pip
        pip.main(['install', 'snpflip'])

    print("Done installing snpflip")


def shapeit():
    import platform
    import urllib.request
    import os
    import sys
    import subprocess

    # Since shapeit only works on linux or mac, we need to first check what system they are on.
    system_check = platform.system()

    if system_check == "Linux":
        print('Downloading shapeit now.')
        urllib.request.urlretrieve(
            'https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r900.glibcv2.12.linux.tar.gz',
            'shapeit.v2.r900.glibcv2.12.linux.tar.gz')
        # Making directory to store program
        os.makedirs('Shapeit_v2.r900_Linux_Static')
        # Unpacking
        subprocess.check_output(['tar','-xzvf','shapeit.v2.r900.glibcv2.12.linux.tar.gz','-C',
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
        subprocess.check_output(['tar','-zxvf','shapeit.v2.r837.MacOSX.tgz','-C','Shapeit_v2.20_Mac/'])

        print("Done downloading shapeit")

    # If they are running this on a windows machine, they cannot proceed because shapeit is *nix only.
    elif system_check == "Windows":
        sys.exit("I'm sorry, I've detected that you're working on a Windows computer and shapeit is a "
                 "linux or unix program only. If you have access to the Penn State clusters, you should run this "
                 "script from there (they are linux).")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("I cannot detect the system you are working on. Exiting now.")


def vcftools():
    import platform
    import urllib.request
    import os
    import sys
    import subprocess

    from os.path import expanduser
    home = expanduser("~")

    # vcftools is easily installed on a linux system, more difficult to install on a mac, and does not have a Windows
    # distribution, so we need to check the platform.
    system_check = platform.system()

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
        subprocess.check_output(['./configure','--prefix=' + os.path.join(home, 'software')])
        subprocess.check_output('make')
        subprocess.check_output(['make','install'])
        # Tell the user that they should add the following folders to their PATH:
        print("I installed vcftools into " + os.path.join(home,'software/bin')
              + ". You should type 'echo $PATH' to check that " + os.path.join(home, 'software/bin')
              + " is in your $PATH variable. If it isn't, then type 'export PATH=$PATH:"
              + os.path.join(home, 'software/bin')
              + " to add it. You should also type in the following to set your PERL5LIB: 'export PERL5LIB="
              + os.path.join(home, 'software/vcftools-vcftools-ea875e2/src/perl')
              + ". I HIGHLY recommend you add both of these lines to your .bash_profile file, or else you'll have to "
                "set this path every time you open a new terminal window.")

        print("Done downloading vcftools")

    # If the user is on a mac
    elif system_check == "Darwin":
        sys.exit('To download vcftools on a Mac you should get homebrew. On your own, download homebrew from '
                 'https://brew.sh/ then download vcftools by typing "brew install homebrew/science/vcftools" in your '
                 'command terminal. Then, follow the directions in the README file for configuring. Then make sure the'
                 'install folder is in your PATH variable.')

    # If they are running this on a windows machine, they cannot proceed because shapeit is *nix only.
    elif system_check == "Windows":
        sys.exit("I'm sorry, I've detected that you're working on a Windows computer and vcftools is a "
                 "linux or unix program only. If you have access to the Penn State clusters, you should run this "
                 "script from there (they are linux).")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("I cannot detect the system you are working on. Exiting now.")


def bcftools():
    import platform
    import urllib.request
    import os
    import sys
    import subprocess

    from os.path import expanduser
    home = expanduser("~")

    # bcftools is easily installed on a linux system, more difficult to install on a mac, and does not have a Windows
    # distribution, so we need to check the platform.
    system_check = platform.system()

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
        subprocess.check_output(['make','install'])
        # Tell the user that they should add the following folders to their PATH:
        print("I installed bcftools into " + os.path.join(home, 'software/bin')
              + " . You should type '$PATH' to check that " + os.path.join(home, 'software/bin')
              + " is in your $PATH variable. IF it isn't, then type 'export PATH=$PATH:"
              + os.path.join(home, 'software/bin') + " to add it. I HIGHLY recommend that you add the export command to"
                                                     "your .bash_profile.")

        print("Done downloading bcftools")

    # If the user is on a mac
    elif system_check == "Darwin":
        sys.exit("It's possible to install bcftools on a Mac, but it might require that you install other "
                 "libraries too. You can try following this tutorial (not for bcftools specifically, but related "
                 "programs) http://www.danielecook.com/installing-tabix-and-samtools-on-mac/ and download homebrew "
                 "(if you haven't already) and xcode. Then run brew install homebrew/science/bcftools and install using"
                 " the directions in the README file. You should also add the install location to your PATH if it is "
                 "not already there.")

    # If they are running this on a windows machine, they cannot proceed because bcftools is *nix only.
    elif system_check == "Windows":
        sys.exit("I'm sorry, I've detected that you're working on a Windows computer and bcftools is a "
                 "linux or unix program only. If you have access to the Penn State clusters, you should run this "
                 "script from there (they are linux).")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("I cannot detect the system you are working on. Exiting now.")


def htslib():
    import platform
    import urllib.request
    import os
    import sys
    import subprocess

    from os.path import expanduser
    home = expanduser("~")

    # htslib is easily installed on a linux system, more difficult to install on a mac, and does not have a Windows
    # distribution, so we need to check the platform.
    system_check = platform.system()

    if system_check == "Linux":
        print('Downloading htslib now.')
        if not os.path.exists(os.path.join(home, 'software')):
            os.mkdir(os.path.join(home, 'software'))
        urllib.request.urlretrieve('https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2',
                                   os.path.join(home, 'software/htslib-1.6.tar.bz2'))
        # Making directory to store program
        # os.makedirs(os.path.join(home,'software/htslib_1.6'))
        # Unpacking
        subprocess.check_output(['tar','-xvjf', os.path.join(home, 'software/htslib-1.6.tar.bz2'),'-C',
                                 os.path.join(home, 'software')])
        # Moving into the htslib folder
        os.chdir(os.path.join(home, 'software/htslib-1.6/'))
        # Running configuration and installation steps
        subprocess.check_output(['./configure','--prefix=' + os.path.join(home, 'software')])
        subprocess.check_output('make')
        subprocess.check_output(['make','install'])
        print("I installed htslib into " + os.path.join(home, 'software/bin')
              + " . You should type '$PATH' to check that " + os.path.join(home, 'software/bin')
              + " is in your $PATH variable. IF it isn't, then type 'export PATH=$PATH:"
              + os.path.join(home, 'software/bin') + " to add it. I highly recommend that you add the export command to"
                                                     "your .bash_profile.")

    # If the user is on a mac
    elif system_check == "Darwin":
        sys.exit("It's possible to install htslib on a Mac, but it might require that you install other "
                 "libraries too. You can try following this tutorial (not for htslib specifically, but related "
                 "programs) http://www.danielecook.com/installing-tabix-and-samtools-on-mac/ and download homebrew "
                 "(if you haven't already) and xcode. Then run brew install homebrew/science/htslib and install using"
                 " the directions in the README file. You should also add the install location to your PATH if it is "
                 "not already there.")
    # If they are running this on a windows machine, they cannot proceed because bcftools is *nix only.
    elif system_check == "Windows":
        sys.exit(
            "I'm sorry, I've detected that you're working on a Windows computer and htslib is a "
            "linux or unix program only. If you have access to the Penn State clusters, you should run this "
            "script from there (they are linux).")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("I cannot detect the system you are working on. Exiting now.")

    print("Done installing htslib")


def samtools():
    import platform
    import urllib.request
    import os
    import sys
    import subprocess

    from os.path import expanduser
    home = expanduser("~")

    # samtools is easily installed on a linux system, more difficult to install on a mac, and does not have a Windows
    # distribution, so we need to check the platform.
    system_check = platform.system()

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
        subprocess.check_output(['make','install'])
        print("I installed samtools into " + os.path.join(home, 'software/bin')
              + " . You should type '$PATH' to check that " + os.path.join(home, 'software/bin')
              + " is in your $PATH variable. IF it isn't, then type 'export PATH=$PATH:"
              + os.path.join(home, 'software/bin') + " to add it. I highly recommend that you add the export command to"
                                                     " your .bash_profile.")

        print("Done downloading samtools")

    # If the user is on a mac
    elif system_check == "Darwin":
        sys.exit("It's possible to install samtools on a Mac, but it might require that you install other "
                 "libraries too. You can try following this tutorial "
                 "http://www.danielecook.com/installing-tabix-and-samtools-on-mac/ and download homebrew "
                 "(if you haven't already) and xcode. You should also add the install location to your PATH if it is "
                 "not already there.")

    # If they are running this on a windows machine, they cannot proceed because samtools is *nix only.
    elif system_check == "Windows":
        sys.exit("I'm sorry, I've detected that you're working on a Windows computer and samtools is a "
                 "linux or unix program only. If you have access to the Penn State clusters, you should run this "
                 "script from there (they are linux).")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("I cannot detect the system you are working on. Exiting now.")


def getmatplotlib():
    try:
        import pip
        # Use pip to install snpflip and it's dependencies.
        pip.main(['install', 'matplotlib'])
    except ImportError:
        pip()
        import pip
        pip.main(['install', 'matplotlib'])

    print("Done installing matplotlib")



