def python3():
    import platform
    import urllib2
    import os
    import sys
    import subprocess

    from os.path import expanduser
    home = expanduser("~")

    # Get what system the user is using
    system_check = platform.system()

    if system_check == "Linux":
        if not os.path.exists(os.path.join(home, 'software')):
            os.makedirs(os.path.join(home, 'software'))
        print("Downloading python3 to " + os.path.join(home, 'software') + " now.")

        f = urllib2.urlopen("https://www.python.org/ftp/python/3.6.3/Python-3.6.3.tgz")
        with open(os.path.join(home, 'software/Python-3.6.3.tgz'), "wb") as code:
            code.write(f.read())

        # Unpacking
        subprocess.Popen(['tar','-xzvf',os.path.join(home, 'software/Python-3.6.3.tgz'),'-C',
                                 os.path.join(home, 'software')])

        # Moving into the Python folder
        os.chdir(os.path.join(home, 'software/Python-3.6.3'))
        with open('InstallPython', "w") as file:
            file.write('./configure --prefix=' + os.path.join(home, 'software')
                       + '\n'
                         'make\n'
                         'make install\n'
                         'cd ' + home
                       + '\n'
                         "echo 'export PATH=" + os.path.join(home, 'software/bin/') + ":$PATH' >> .bashrc"
                       + '\n'
                         "echo 'export PATH=" + os.path.join(home, 'software/Python-3.6.3/') + ":$PATH' >> .bashrc"
                       + '\n'
                         'echo "export PYTHONPATH=' + os.path.join(home, 'software/Python-3.6.3/" >> .bashrc')
                       + '\n'
                         'source .bashrc')

        sys.exit("Please move to the " + os.path.join(home, 'software/Python-3.6.3')
                 + " folder and type in 'source InstallPython'. Then re-run this script.")

    # If the user is on a mac
    elif system_check == "Darwin":
        sys.exit("I've detected that you are working on a Mac. Your best bet is to follow these directions "
                 "(http://docs.python-guide.org/en/latest/starting/install3/osx/) and download"
                 "python3 using xcode and homebrew. You will probably need homebrew later in this script too.")

    # If they are running this on a windows machine
    elif system_check == "Windows":
        sys.exit("I've detected that you are working on a Windows computer. Your best bet is to follow "
                 "these directions (http://docs.python-guide.org/en/latest/starting/install3/win/) and download "
                 "python3 mnaually. As a word of warning, you can do many functions of this script on a Windows "
                 "computer, but some (admixture, phasing, imputation) requires access to a linux based computational "
                 "cluster (like the Penn State ACI-B system)")

    # If I cannot detect what system they're on, force exit.
    else:
        sys.exit("I cannot detect the system you are working on. Exiting now.")
