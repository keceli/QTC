#!/usr/bin/env python
"""
Contains IO and OS related tools.
Since python have different modules for various
IO related functionalities, it is good to have
a single module simplifying their usage.
TODO: Add unit tests. Seperate IO vs OS
"""
import time
import os
from os.path import isfile

__updated__ = "2017-05-19"


def get_date():
    """
    Returns the current date and time
    """
    import datetime
    return datetime.datetime.now()

def touch(fname, times=None):
    """
    Creates a file with the given fname, aka unix touch
    See http://stackoverflow.com/questions/1158076/implement-touch-using-python
    """
    import os
    with open(fname, 'a'):
        os.utime(fname, times)
    return


def rm(fname):
    """
    Deletes a file with the given fname.
    """
    import os
    os.remove(fname)
    return


def mkdir(path):
    """
    Creates directory if it doesn't exist.
    Check http://stackoverflow.com/questions/273192/how-to-check-if-a-directory-exists-and-create-it-if-necessary
    """
    import os
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    return


def cd(path):
    """
    Change working directory
    """
    import os
    os.chdir(path)
    return


def pwd():
    """
    Return current directory
    """
    import os
    return os.getcwd()


def join_path(*paths):
    """
    Concatenes strings into a portable path using correct seperators.
    """
    import os
    return os.path.join(*paths)

def write_file(s, filename='newfile'):
    """
    Writes s string to a file with the given'filename'.
    """
    with open(filename, 'w') as f:
        f.write(s)
    return


def read_file(filename, aslines=False):
    """
    Reads a file and return either a list of lines or a string.
    """
    with open(filename, 'r') as f:
        if aslines:
            tmp = f.readlines()
        else:
            tmp = f.read()
    return tmp


def check_file(filename, timeout=0):
    """
    Returns True (False) if a file exists (doesn't exist).
    If timeout>0 is given, then checks file in a loop until
    timeout seconds pass.
    """
    import time
    from os.path import isfile
    exists = isfile(filename)
    if timeout > 0 and not exists:
        t0 = time.time()
        while time.time() < t0 + timeout:
            if isfile(filename):
                exists = True
                break
    return exists


def check_dir(dirname, timeout=0):
    """
    Returns True (False) if a file exists (doesn't exist).
    If timeout>0 is given, then checks file in a loop until
    timeout seconds pass.
    """
    import time
    from os.path import isdir
    exists = isdir(dirname)
    if timeout > 0 and not exists:
        t0 = time.time()
        while time.time() < t0 + timeout:
            if isdir(dirname):
                exists = True
                break
    return exists


def check_exe(exename):
    """
    Check if an executable is available.
    TODO:
    """
    import distutils.spawn as ds

    return ds.find_executable(exename)


def cp(source, target):
    """
    Copies a file, source and target are strings for paths.
    Target can be a directory or a file.
    """
    from shutil import copy
    copy(source, target)
    return


def get_path(f,executable=False):
    """
    Returns absolute path for a file or folder.
    """
    import os
    import distutils.spawn as ds
    if executable:
        return ds.find_executable(f)
    else:
        return os.path.abspath(f)


def get_line_number(keyword, lines=None, filename=None):
    """
    Returns the line number of a keyword found in given lines of string.
    Returns -1 if keyword is not found
    """
    num = -1
    if lines is None and filename is None:
        print 'List of lines or a filename to be read is required for get_line_number'
    elif filename:
        lines = read_file(filename, aslines=True)

    for n, line in enumerate(lines):
        if keyword in line:
            num = n
    return num

def get_git_version():
    """
    Return git version.
    """
    import subprocess
    return subprocess.check_output(["git", "describe"])


def read_list(listfile):
    """
    Return a list of strings from all lines in a text file.
    Skips blank lines.
    """
    with open(listfile, 'r') as f:
        lines = filter(None, (line.rstrip() for line in f))
    return lines


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

