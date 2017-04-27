#!/usr/bin/env python

"""
Quantum chemistry tools.
"""
import subprocess
import obtools as ob    
import argparse
import iotools as io




def get_heavyatomlist(formula):
    """
    Return a string that only contains nonhydrogen elements of the formula
    without stoichemetry information
    """
    s = ''.join(i for i in formula if not i.isdigit() and not i =='H' )
    if s == '':
        s = 'H'
    return s


def get_heavyatomformula(formula):
    """
    Return a string that exludes Hydrogens in the molecular formula
    """
    n = len(formula)
    s=''
    hdigit = False
    for i in range(n):
        if formula[i].isdigit():
            if not hdigit:
                s += formula[i]
        elif formula[i] == 'H':
            hdigit = True
        else:
            hdigit = False
            s += formula[i]
    if s == '':
        s = formula                        
    return s


def run_mopac(s=None, mol=None,mopackeys=None):
    import subprocess
    import obtools as ob
    import iotools as io
    if mol is None and s is not None:    
        mol = ob.get_mol(s)
    key = ob.get_uniquekey(mol)    
    if mopackeys is None:
        mopackeys = 'pm3 precise nosym THREADS=1 opt'
    inputfile = key + '.mop'
    outfile = key + '.out'
    if not io.check_file(inputfile):
        inputstr  = mol.write(format='mop',
                  opt={'k': mopackeys})
        inputstr += '\n' + mopackeys.replace('opt','oldgeo tctools')
        io.write_file(inputstr,filename=inputfile)
        io.check_file(inputfile,timeout=3)
    if not io.check_file(outfile):    
        mopacexe='mopac'
        subprocess.call([mopacexe,inputfile])
        io.check_file(outfile,timeout = 1)
    lines = io.read_file(outfile)
    return lines

def run_qc(qcexe,inputfile,stdout=False):
    import subprocess
    #          opt={'k': 'pm3  precise nosym THREADS=1 opt tctools'})
    outputfile = inputfile + '.out'
    if stdout:
        cmd = [qcexe,inputfile,outputfile]
    else:
        cmd  = [qcexe,inputfile]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out, err = proc.communicate()
    if err is None:
        io.write_file(out,outputfile)
    else:
        print err
    return

def run_nwchem(s,nwchem='nwchem'):
    import iotools as io
    import obtools as ob
    
    mol = ob.get_mol(s)
    input = get_input_text(mol=mol)
    inchikey = mol.write(format='inchikey').replace('\n','')
    inputfile = inchikey+'nw'
    io.write_file(input,filename=inputfile)
    run_qc('nwchem',inputfile,stdout=True)
    return s


def get_input_text(mol=None,s=None,template='qc_template.txt'):
    """
    Returns the text for quantum chemistry input file based on a template.
    """
    if mol is None:
        mol = ob.get_mol(s)
    inchikey = mol.write(format='inchikey').replace('\n','')
    xyzlines = mol.write(format='xyz').splitlines(True)
    xyzstr   = ''.join(xyzlines[2:])
    with open(template,'r') as f:
        tmp = f.read()
    tmp = tmp.replace('QCPUT(GEO)\n',xyzstr)
    tmp = tmp.replace('QCPUT(UNIQUEKEY)',inchikey)
    tmp = tmp.replace('QCPUT(MULTIPLICITY)',str(mol.spin))
    tmp = tmp.replace('QCPUT(BASIS)','cc-pvdz')
    if mol.spin == 1:
        tmp = tmp.replace('QCPUT(MULTIPLICITY_WORD)','singlet')
        tmp = tmp.replace('QCPUT(SCFTYPE)','rhf')
    elif mol.spin == 2:
        tmp = tmp.replace('QCPUT(MULTIPLICITY_WORD)','doublet')
        tmp = tmp.replace('QCPUT(SCFTYPE)','uhf')
    elif mol.spin == 3:
        tmp = tmp.replace('QCPUT(MULTIPLICITY_WORD)','triplet')
        tmp = tmp.replace('QCPUT(SCFTYPE)','uhf')
    return tmp


def get_natom_mopac(lines):
    """
    Return the number of atoms from mopac output
    """
    import iotools as io
    keyword = 'Empirical Formula'
    n = io.get_line_number(keyword,lines=lines)
    natom = int(lines[n].split()[-2])
    return natom


def get_xyz_mopac(lines):
    """
    Returns xyz as a list of strings from mopac output lines.
    mopac output:
    ORIENTATION OF MOLECULE IN FORCE CALCULATION

    NO.       ATOM         X         Y         Z

     1          O      -0.5848    0.0000    0.0000
     2          O       0.5848    0.0000   -0.0000
     xyz is a list of natom+2 string lines
     natom
     comment
     atom1 x1 y1 z1
     .
     .
     .
     atomn xn yn zn
    """
    import iotools as io
    natom = get_natom_mopac(lines)
    keyword = "ORIENTATION OF MOLECULE IN FORCE CALCULATION"
    xyzline = io.get_line_number(keyword,lines=lines) +4
    xyz = '{0}\n'.format(natom)
    comment = '\n'
    xyz += comment
    for i in range(natom):
        xyz += ' '.join(lines[xyzline+i].split()[1:]) + '\n'
    return xyz


def get_freq_mopac(lines):
    """
    Returns a float list of vibrational frequencies in cm-1.
    """
    keyword = 'FREQ.' 
    natom = get_natom_mopac(lines)
    freqs = ['']*(3*natom - 5)
    i = 0
    for line in lines:
        if keyword in line:
            freqs[i] = line.split()[1]
            i += 1
    return freqs[:i]  


def get_zpe_mopac(lines):
    """
    Return zero point energy in kcal/mol from mopac output.
    
    """
    keyword = 'ZERO POINT ENERGY'
    n = io.get_line_number(keyword,lines=lines)
    return float(lines[n].split()[3])


def get_deltaH_mopac(lines):
    """
    Return delta H in kcal/mol from mopac output.
    
    """
    keyword = 'FINAL HEAT OF FORMATION'
    n = io.get_line_number(keyword,lines=lines)
    return float(lines[n].split()[5])


def print_list(s):
    return s

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

