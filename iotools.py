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

__updated__ = "2017-06-23"


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


def mv(oldname,newname):
    """
    Renames or moves a file
    """
    import os
    os.rename(oldname,newname)
    return 

def pwd():
    """
    Return current directory
    """
    import os
    return os.getcwd()


def find_files(directory, pattern):
    """
    Yields files in a directory (including subdirectories) with a given pattern
    https://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
    """
    import os, fnmatch
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename

                
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


def check_file(filename, timeout=0, verbose=False):
    """
    Returns True (False) if a file exists (doesn't exist).
    If timeout>0 is given, then checks file in a loop until
    timeout seconds pass.
    If verbose is True and file does not exist, prints an error message
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
    if not exists and verbose:
        print('"{0}" file not found.'.format(filename))
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


def get_line_number(keyword, lines=None, filename=None,getlastone=False):
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
            if getlastone:
                num = n
            else:
                return n
    return num


def get_line_numbers(keyword, lines=None, filename=None):
    """
    Returns the line numbers as a list for a keyword in given lines of string.
    Returns -1 if keyword is not found
    """
    if lines is None and filename is None:
        print 'List of lines or a filename to be read is required for get_line_numbers'
    elif filename:
        lines = read_file(filename, aslines=True)
    nums = []
    for n, line in enumerate(lines):
        if keyword in line:
            nums.append(n)
    if len(nums) == 0:
        nums = -1
    return nums


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


def execute_old(exe,inp=None,out=None):
    """
    Executes a calculation for a given input file inp and executable exe.
    """
    from subprocess import Popen, PIPE
    if not check_exe(exe):
        return 'Executable "{0}" not found.\n'.format(exe)
    if inp:
        if not check_file(inp,1):
            return 'Input file "{0}" not found.\n'.format(get_path(inp))
    if inp and out:
        process = Popen([exe, inp, out], stdout=PIPE, stderr=PIPE)
    elif inp:
        process = Popen([exe, inp], stdout=PIPE, stderr=PIPE)
    else:
        process = Popen([exe], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    if stderr is None or stderr == '':
        msg = 'Run {0} {1}: Success.\n'.format(exe, inp)
    else:
        errstr = """ERROR in "{0}"\n
        STDOUT:\n{1}\n
        STDERR:\n{2}""".format(inp, stdout, stderr)
        errfile = inp + '.err'
        write_file(errstr, errfile)
        msg = 'Run {0} {1}: Failed, see "{2}"\n'.format(exe, inp, get_path(errfile))
    return msg


def execute(command, stdoutfile=None, stderrfile=None, merge=False):
    """
    Executes a given command, and optionally write stderr and/or stdout.
    Parameters
    ----------
    command: List of strings, where a command line is seperated into words.
    stderrfile: None or a string for a file name to write stderr
    stdoutfile: None or a string for a file name to write stdout
    
    Returns
    ---------
    If stdoutfile:
        A string describing the success or failure of the calculation
    Else:
        stdout
        
    Doctest
    ---------    
    >>> io.execute(['echo','this works'])
    'this works\n'     
    >>> io.execute('echo this also works')
    'this also works\n'     
    """
    from subprocess import Popen, PIPE
    if type(command) == str:
        commandstr = command
        command = command.split()
    else:
        commandstr = ' '.join(command)
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    msg = ''
    if merge:
        if type(out) == str and type(err) == str:
            out += err
    if out is None or out == '':
        msg = 'No STDOUT\n'
    else:
        if stdoutfile:
            write_file(out,stdoutfile)
        else:
            msg += 'STDOUT\n'  
    if err is None or err == '':
        msg += 'Run {0}: Success.\n'.format(commandstr)
    else:
        msg += '\nRun {0}: Error: \n "{1}"\n'.format(commandstr, get_path(err))
        if stderrfile:
            write_file(msg + err, stderrfile)      
    return msg

def db_head_path(db_location=None):

    if db_location == None:
        return '/home/elliott/testdirectory/'
    else:
        return db_location


def db_smiles_path(smiles, db_location = None):

    import obtools as ob
    mol = ob.get_mol(smiles)
    directory = ob.get_smiles_path(mol)
    if db_location == None:
        db_location = db_head_path()

    return join_path(db_location, directory)


def prog_meth_bas_path(prog, method, basis):
    return prog + '__' + method + '__' + basis + '/'


def db_opt_path(prog, method, basis, db_location=None, smiles=None):

    if smiles == None:
        directory = db_head_path(db_location)
    else:
        directory = db_smiles_path(smiles, db_location)
    if method == None:
        return join_path(directory, 'sp/')
    return join_path(directory, prog_meth_bas_path(prog, method, basis))


def db_sp_path(prog, method, basis, db_location=None, smiles=None, optprog=None, optmethod=None, optbasis=None):
    
    directory = db_opt_path(optprog, optmethod, optbasis, db_location, smiles)
    if prog == None:
        if optprog == None:
            return db_head_path(db_location)
        return directory
    return join_path(directory, prog_meth_bas_path(prog,method,basis))


def db_store_opt_prop(s, smiles, typ='zmat', db_location=None, prog=None, method=None, basis=None):

    if prog == None:
        directory =  db_head_path(db_location)
    else:
        directory =  db_opt_path(prog, method, basis, db_location, smiles)

    mkdir(directory)
    write_file(s, join_path(directory, smiles + '.' + typ))
    return

def db_get_opt_prop(smiles, typ='zmat', db_location=None, prog=None, method=None, basis=None):

    if prog == None:
        directory =  db_head_path(db_location)
    else:
        directory =  db_opt_path(prog, method, basis, db_location, smiles)

    if check_file(join_path(directory, smiles + '.' + typ)):
        return read_file(join_path(directory, smiles + '.' + typ))
    return 

def db_store_sp_prop(s, smiles, typ='ene', db_location=None, prog=None, method=None, basis=None, optprog=None, optmethod=None, optbasis=None):
       
    if prog == None:
        directory =  db_head_path(db_location)
    elif optprog == None:
        directory = db_sp_path(prog, method, basis, db_location, smiles, prog, method, basis)
    else:
        directory =  db_sp_path(prog, method, basis, db_location, smiles, optprog, optmethod, optbasis)

    mkdir(directory) 
    write_file(s, join_path(directory, smiles + '.' + typ))
    return 
    
def db_get_sp_prop(smiles, typ='ene', db_location=None, prog=None, method=None, basis=None, optprog=None, optmethod=None, optbasis=None):

    if prog == None:
        directory =  db_head_path(db_location)
    elif optprog == None:
        directory = db_sp_path(prog, method, basis, db_location, smiles, prog, method, basis)
    else:
        directory =  db_sp_path(prog, method, basis, db_location, smiles, optprog, optmethod, optbasis)

    if check_file(join_path(directory, smiles + '.' + typ)):
        return read_file(join_path(directory, smiles + '.' + typ))
    return 

def db_logfile_dir(db_location, smiles=None, prog=None, method=None, basis=None, optprog=None, optmethod=None, optbasis=None):
    directory =  db_sp_path(prog, method, basis, db_location, smiles, optprog, optmethod, optbasis)
    mkdir(join_path(directory, 'logfiles/')) 
    return  join_path(directory, 'logfiles/')

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

