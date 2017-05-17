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


__updated__ = "2017-05-17"


def check_mopac():
        return True


def get_symbol(atomno):
    """
    Returns the element symbol for a given atomic number.
    Returns 'X' for atomno=0
    >>>print get_symbol(1)
    >>>H
    """
    syms = ['X',
            'H', 'He'
            'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
    return syms[atomno]


def run_qcscript(qcscriptpath, inputpath, geopath, multiplicity):
    """
    Submit jobs using Ahren's script.
     Usage: qcscript.pl input.qcscript initialgeometry.geo multiplicity
     Sample geo file:
     C          0.96737       -0.07578        0.02761
     O          2.33437       -0.07578        0.02761
    """
    from subprocess import Popen, PIPE
    process = Popen(['perl',qcscriptpath, inputpath, geopath, str(multiplicity)], stdout=PIPE, stderr=PIPE)
    msg, err = process.communicate()
    if err:
        msg = 'Failed {0}'.format(err)
    return msg


def execute(inp, exe):
    """
    Executes a calculation for a given input file inp and executable exe.
    """
    from subprocess import Popen, PIPE
    process = Popen([exe, inp], stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    
    if err is None or err == '':
        msg = 'Run {0} {1}: Success.\n'.format(exe, inp)
    else:
        errstr = """ERROR in {0}\n
        STDOUT:\n{1}\n
        STDERR:\n{2}""".format(inp, out, err)
        errfile = inp + '.err'
        io.write_file(errstr, errfile)
        msg = 'Run {0} {1}: Failed, see {2}.\n'.format(exe, inp, io.get_path(errfile))
    return msg



def execute_gaussian(inp, exe='g09'):
    """
    Runs gaussian calculation.
    """
    from subprocess import Popen, PIPE
    process = Popen([exe, inp], stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    
    if err is None or err == '':
        msg = 'Run {0} {1}: Success.'.format(exe, inp)
    else:
        errstr = """ERROR in {0}\n
        STDOUT:\n{1}\n
        STDERR:\n{2}""".format(inp, out, err)
        errfile = inp + '.err'
        io.write_file(errstr, errfile)
        msg = 'Run {0} {1}: Failed, see {2}.'.format(exe, inp, io.get_path(errfile))
    return msg


def run_gaussian(s, exe='g09', template='qctemplate.txt',mult=0,overwrite=False):
    """
    Runs gaussian calculation
    """
    mol = ob.get_mol(s, make3D=True)
    if mult == 0:
        mult = ob.get_multiplicity(mol)
    tmp = io.read_file(template)    
    inptext = get_gaussian_input(mol, tmp, mult)
    prefix = ob.get_unique_key(mol, mult)
    inpfile = prefix + '.g09'  
    outfile = prefix + '.log'
    if io.check_file(outfile, timeout=1):
        if overwrite:
            msg = "Overwriting previous calculation {0}\n".format(io.get_path(outfile))
            run = True
        else:
            msg = 'Skipping calculation, found {0}\n'.format(io.get_path(outfile))
            run = False
    else:
        run = True
    if run:
        if not io.check_file(inpfile, timeout=1):
            io.write_file(inptext, inpfile)
        if io.check_file(inpfile, timeout=1):
            msg = execute(inpfile, exe)
            if io.check_file(outfile, timeout=1):
                msg += ' Output file: {0}\n'.format(io.get_path(outfile))
        else:
            msg = 'Failed, cannot find input file {0}\n'.format(io.get_path(inpfile))
    return msg


def run_mopac(s, exe='mopac', method='pm7', mopackeys='precise nosym threads=1 opt', mult=0, overwrite=False):
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
    if io.check_file(outfile, timeout=1):
        if overwrite:
            msg = "Overwriting previous calculation {0}\n".format(io.get_path(outfile))
            run = True
        else:
            msg = 'Skipping calculation, found {0}\n'.format(io.get_path(outfile))
            run = False
    else:
        run = True
    if run:
        if not io.check_file(inpfile, timeout=1):
            io.write_file(inptext, inpfile)
        if io.check_file(inpfile, timeout=1):
            msg = execute(inpfile, exe)
            if io.check_file(outfile, timeout=1):
                msg += ' Output file: {0}.\n'.format(io.get_path(outfile))
        else:
            msg = ''
        if not io.check_file(outfile, timeout=1):
            msg = 'Failed, can not find {0}.\n'.format(io.get_path(outfile))
    return msg


def execute_mopac(inp, exe='mopac'):
    """
    Runs mopac calculation.
    Mopac is a fortran code that does not return an error code
    but writes error to stderr.
    If there is no error stderr = None or ''
    For some keyword errors, still no stderr or stdout is provided,
    so additional checks are required.
    """
    from subprocess import Popen, PIPE
    import iotools as io
#    subprocess.call([mopacexe, inp])
    process = Popen([exe, inp], stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    if err is None or err == '':
        msg = 'Run {0} {1}: Success.'.format(exe, inp)
    else:
        errstr = """ERROR in {0}\n
        STDOUT:\n{1}\n
        STDERR:\n{2}""".format(inp, out, err)
        errfile = inp + '.err'
        io.write_file(errstr, errfile)
        msg = 'Run {0} {1}: Failed, see {2}\n'.format(exe, inp, errfile)
    return msg


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


def get_gaussian_input(x, template, mult=0):
    """
    Returns Gaussian input file based on a given template.
    """
    if type(x) is str:
        mol = ob.get_mol(x)
    else:
        mol = x
    zmat = ob.get_zmat(mol)
    if mult == 0:
        mult = ob.get_multiplicity(mol)
    charge = ob.get_charge(mol)
    uniquekey = ob.get_unique_key(mol, mult)
    inp = template.replace("QTC(CHARGE)", str(charge))
    inp = inp.replace("QTC(MULTIPLICITY)", str(mult))
    inp = inp.replace("QTC(UNIQUEKEY)", uniquekey)
    inp = inp.replace("QTC(ZMAT)", zmat)
    if "QTC(" in inp:
        print("Error in template file:\n" + inp)
        return
    return inp


def get_gaussian_natom(lines):
    """
    NAtoms=     30 NQM=       30 NQMF=       0 NMMI=      0 NMMIF=      0
    """
    import iotools as io
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'NAtoms='
    n = io.get_line_number(keyword, lines=lines)
    return int(lines[n].split()[1])


def get_gaussian_xyz(lines):
    """
    TODO
                             Ref.Geom.
 ---------------------------------------------------------------------
 Center     Atomic                     Coordinates (Angstroms)
 Number     Number                        X           Y           Z
 ---------------------------------------------------------------------
      1          6                    0.639026    0.420560    0.000000
      2          1                    0.639026    1.074418    0.870857
      3          1                    0.639026    1.074418   -0.870857
      4          6                   -0.639026   -0.420560    0.000000
      5          1                   -0.639026   -1.074418    0.870857
      6          1                   -0.639026   -1.074418   -0.870857
    """
    import iotools as io
    if type(lines) == str:
        lines = lines.splitlines()
    natom = get_gaussian_natom(lines)
    keyword = 'Input orientation:'
    n = io.get_line_number(keyword, lines=lines)
    
    return 


def get_gaussian_xmatrix(lines):
    """
    TODO
    Return anharmonic X matrix from Gaussian log file.
     X matrix of Anharmonic Constants (cm-1)
                1             2             3             4             5
      1       -16.516
      2       -66.957       -16.697
      3       -58.351       -61.032       -14.183
      4       -59.910       -58.652       -56.823       -14.169
      5         0.070         7.271       -12.551        -5.852        -3.033
      6        -6.314       -10.945        -3.626       -10.897        -1.176
      7        -8.142        -6.253        -4.086        -4.940        -9.433
      8        -6.350        -0.155        -3.440        -5.016       -15.112
      9        -5.865        -5.516        -6.812        -9.987        -3.641
     10       -11.927        -8.007        -4.714        -4.513        -5.103
     11        -9.270       -13.121        -6.251        -5.246        -6.107
     12        -2.600        -5.099        -2.873         1.830        -0.766
                6             7             8             9            10
      6        -2.352
      7        -4.938        -0.630
      8        -5.973        -2.335        -0.776
      9         0.893        -2.619         0.133        -1.539
     10        -2.750        -1.897         1.474        -8.260         0.565
     11        -1.733        -1.237         3.059        -6.832        -2.110
     12        -2.883        -0.846        -5.684         1.507         5.788
               11            12
     11         0.933
     12         1.571         1.650
     
     Alternative:
     
     No anharmonic for linear molecules
     NO Anharmonic anal. for Linear Molecules
    """
    return


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

