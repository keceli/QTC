#!/usr/bin/env python
"""
Quantum chemistry tools.

"""
import subprocess
import obtools as ob
import argparse
import iotools as io
import numpy as np
from iotools import write_file


__updated__ = "2017-05-03"


def check_mopac():
        return True


def run_mopac(s, mopacexe='mopac', method='pm7', mopackeys='precise nosym threads=1 opt', mult=0, overwrite=False):
    """
    Runs mopac calculation.
    If overwrite=False, no calculation if output file already exists.
    """
    mol = ob.get_mol(s, make3D=True)
    if mult == 0:
        mult = ob.get_multiplicity(mol)
    inptext = get_mopac_input(mol, method=method, keys=mopackeys, mult=mult, dothermo=True)
    prefix = ob.get_unique_key(mol, mult)
    inpfile = prefix + '.mop'
    outfile = prefix + '.out'
    errcode = 0
#     formula = ob.get_formula(mol)
#     formula_noH = ob.get_formula(mol, stoichemetry=True, hydrogens=False)
#     elements_noH = ob.get_formula(mol, stoichemetry=False, hydrogens=False)
#     uniquekey = ob.get_unique_key(mol, mult)
#     dirs = 'database', elements_noH, formula_noH, formula, uniquekey, method
#     dirpath = io.join_path(*dirs)
#     io.mkdir(dirpath)
#     cwd = io.pwd()
#     io.cd(dirpath)
    if io.check_file(outfile, timeout=1):
        if overwrite:
            print "Overwriting previous calculation {0}".format(outfile)
            run = True
        else:
            print 'Skipping calculation, found {0}'.format(outfile)
            run = False
    else:
        run = True
    if run:
        if not io.check_file(inpfile, timeout=1):
            io.write_file(inptext, inpfile)
        if io.check_file(inpfile, timeout=1):
            errcode = execute_mopac(inpfile, mopacexe)
        else:
            errcode = 3
        if not io.check_file(outfile, timeout=1):
            errcode += 10
    return '{0} : {1}'.format(io.get_path(outfile), errcode)


def execute_mopac(inp, mopacexe='mopac'):
    """
    Runs mopac calculation.
    Mopac is a fortran code that does not return an error code
    but writes error to stderr.
    If there is no error stderr = None
    For some keyword errors, no stderr output is provided,
    but still
    it does not also write to stdout.
    """
    from subprocess import Popen, PIPE
    import iotools as io
#    subprocess.call([mopacexe, inp])
    process = Popen([mopacexe, inp], stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    if err is None:
        errcode = 0
    elif err == '':
        errcode = 1
    else:
        errcode = 2
        errstr = """ERROR in {0}\n
        STDOUT:\n{1}\n
        STDERR:\n{2}""".format(inp, out, err)
        errfile = inp + '.err'
        io.write_file(errstr, errfile)
    return errcode


def run_qc(qcexe, inputfile, stdout=False):
    import subprocess
    #          opt={'k': 'pm3  precise nosym THREADS=1 opt tctools'})
    outputfile = inputfile + '.out'
    if stdout:
        cmd = [qcexe, inputfile, outputfile]
    else:
        cmd = [qcexe, inputfile]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out, err = proc.communicate()
    if err is None:
        io.write_file(out, outputfile)
    else:
        print err
    return


def run_nwchem(s, nwchem='nwchem'):
    import iotools as io
    import obtools as ob

    mol = ob.get_mol(s)
    input = get_input_text(mol=mol)
    inchikey = mol.write(format='inchikey').replace('\n', '')
    inputfile = inchikey + 'nw'
    io.write_file(input, filename=inputfile)
    run_qc('nwchem', inputfile, stdout=True)
    return s


def get_input_text(mol=None, s=None, template='qc_template.txt'):
    """
    Returns the text for quantum chemistry input file based on a template.
    """
    if mol is None:
        mol = ob.get_mol(s)
    inchikey = mol.write(format='inchikey').replace('\n', '')
    xyzlines = mol.write(format='xyz').splitlines(True)
    xyzstr = ''.join(xyzlines[2:])
    with open(template, 'r') as f:
        tmp = f.read()
    tmp = tmp.replace('QCPUT(GEO)\n', xyzstr)
    tmp = tmp.replace('QCPUT(UNIQUEKEY)', inchikey)
    tmp = tmp.replace('QCPUT(MULTIPLICITY)', str(mol.spin))
    tmp = tmp.replace('QCPUT(BASIS)', 'cc-pvdz')
    if mol.spin == 1:
        tmp = tmp.replace('QCPUT(MULTIPLICITY_WORD)', 'singlet')
        tmp = tmp.replace('QCPUT(SCFTYPE)', 'rhf')
    elif mol.spin == 2:
        tmp = tmp.replace('QCPUT(MULTIPLICITY_WORD)', 'doublet')
        tmp = tmp.replace('QCPUT(SCFTYPE)', 'uhf')
    elif mol.spin == 3:
        tmp = tmp.replace('QCPUT(MULTIPLICITY_WORD)', 'triplet')
        tmp = tmp.replace('QCPUT(SCFTYPE)', 'uhf')
    return tmp


def get_mopac_input(x, method='pm3', keys='precise nosym threads=1 opt', mult=1, dothermo=False):
    """
    Returns mopac input as a string.
    Note: For doctest I had to escape newline characters \n as \\n
    Since it gives EOL error.
    Note2: Doctest is also sensitive to whitespace at the end of lines.
    Hence, I used .strip() to awoid unnecessary whitespace.
    >>> xyz = "2\\n \\n H 0. 0. 0.\\n H 0. 0. 0.9\\n  \\n"
    >>> print get_mopac_input(xyz,method='pm7',dothermo=True)
    pm7 precise nosym threads=1 opt
    <BLANKLINE>
    <BLANKLINE>
    H   0.00000 1  0.00000 1  0.00000 1
    H   0.00000 1  0.00000 1  0.90000 1
    <BLANKLINE>
    pm7 precise nosym threads=1 oldgeo thermo
    """
    if type(x) is str:
        mol = ob.get_mol(x)
    else:
        mol = x
    multDictionary = {
                        1: '',
                        2: 'uhf doublet',
                        3: 'uhf triplet',
                        4: 'uhf quartet',
                        5: 'uhf quintet',
                        6: 'uhf sextet',
                        7: 'uhf septet',
                        8: 'uhf octet',
                        9: 'uhf nonet',
                    }
    keys = method + ' ' + keys + ' ' + multDictionary[mult]
    inp = ob.get_mop(mol, keys=keys.strip())
    if dothermo:
        inp += '\n' + keys.replace('opt', 'oldgeo thermo').strip()
    return inp


def get_mopac_natom(lines):
    """
    Return the number of atoms from mopac output
    >>> s = io.read_file('test/input.out')
    >>> print get_mopac_natom(s)
    5
    """
    import iotools as io
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'Empirical Formula'
    n = io.get_line_number(keyword, lines=lines)
    natom = int(lines[n].split()[-2])
    return natom


def get_mopac_xyz(lines):
    """
    Returns xyz string from mopac output lines.
    >>> s = io.read_file('test/input.out')
    >>> print get_mopac_xyz(s)
    5
    <BLANKLINE>
    C 0.0000 -0.0000 0.0000
    H 1.0870 -0.0000 0.0000
    H -0.3623 1.0248 0.0000
    H -0.3624 -0.5124 0.8875
    H -0.3624 -0.5124 -0.8875
    <BLANKLINE>
    """
    import iotools as io
    if type(lines) == str:
        lines = lines.splitlines()
    natom = get_mopac_natom(lines)
    keyword = "ORIENTATION OF MOLECULE IN FORCE CALCULATION"
    xyzline = io.get_line_number(keyword, lines=lines) + 4
    xyz = '{0}\n'.format(natom)
    comment = '\n'
    xyz += comment
    for i in range(natom):
        xyz += ' '.join(lines[xyzline + i].split()[1:]) + '\n'
    return xyz


def get_mopac_freq(lines):
    """
    Returns a float list of vibrational frequencies in cm-1.
    >>> s = io.read_file('test/input.out')
    >>> print get_mopac_freq(s)
    [ 1362.21  1362.44  1362.56  1451.04  1451.06  3207.4   3207.46  3207.59
      3310.99]
    """
    import numpy as np
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'FREQ.'
    natom = get_mopac_natom(lines)
    freqs = np.zeros(3 * natom - 5)
    i = 0
    for line in lines:
        if keyword in line:
            freqs[i] = float(line.split()[1])
            i += 1
    return freqs[:i]


def get_mopac_zpe(lines):
    """
    Return zero point energy in kcal/mol from mopac output.
    >>> s = io.read_file('test/input.out')
    >>> print get_mopac_zpe(s)
    28.481
    """
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'ZERO POINT ENERGY'
    n = io.get_line_number(keyword, lines=lines)
    return float(lines[n].split()[3])


def get_mopac_deltaH(lines):
    """
    Return delta H in kcal/mol from mopac output.
    >>> s = io.read_file('test/input.out')
    >>> print get_mopac_deltaH(s)
    -13.02534
    """
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'FINAL HEAT OF FORMATION'
    n = io.get_line_number(keyword, lines=lines)
    return float(lines[n].split()[5])


def print_list(s):
    return s

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

