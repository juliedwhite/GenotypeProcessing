def python3():
    import platform
    import urllib2
    import os
    import sys
    import tarfile
    from os.path import expanduser
    home = expanduser("~")

    # Get what system the user is using
    system_check = platform.system()

    if system_check == "Linux":
        if not os.path.exists(os.path.join(home, 'software')):
            os.makedirs(os.path.join(home, 'software'))
        softwaredir = os.path.join(home, 'software')
        print("Downloading python3 to " + softwaredir)

        f = urllib2.urlopen("https://www.python.org/ftp/python/3.6.3/Python-3.6.3.tgz")
        with open(os.path.join(softwaredir, 'Python-3.6.3.tgz'), "wb") as code:
            code.write(f.read())

        with open(os.path.join(home, 'software/InstallPython'), "w") as file:
            file.write('cd ' + os.path.join(softwaredir, 'Python-3.6.3') + '\n'
                       + './configure --prefix=' + softwaredir + '\n'
                       + 'make\n'
                         'make install\n'
                         'cd ' + home + '\n'
                       + "echo 'export PATH=" + os.path.join(softwaredir, 'bin/') + ":$PATH' >> .bash_profile\n"
                       + "echo 'export PATH=" + os.path.join(softwaredir, 'Python-3.6.3/') + ":$PATH' >> "
                                                                                             ".bash_profile\n"
                       + "echo 'export PYTHONPATH=" + os.path.join(softwaredir, 'Python-3.6.3/') + "' >> "
                                                                                                   ".bash_profile\n"
                       + 'source ' + os.path.join(home, '.bash_profile'))

        # Unpacking
        tar = tarfile.open(os.path.join(softwaredir, 'Python-3.6.3.tgz'))
        tar.extractall(softwaredir)
        tar.close()

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
