#!/usr/bin/env python
import argparse
import numpy as np
import datetime
import time
import subprocess
import os
from os.path import isfile
import iotools as io
import obtools as ob
import qctools as qc

"""
Thermochemistry tools.
Requires:
Quantum chemistry code, NWChem, MOPAC, ...
MESS partition function code
PAC99 
thermp
"""    
__updated__ = "2017-04-27"
  
def get_stoichometry(formula,element):
    """
    Returns the stoichometry (count) of an element in a given formula
    Note: Case insensitive
    >>> [get_stoichometry(f,'H') for f in ['C2H4', 'O2', 'CH', 'H2O', 'HO', 'CH3OH']]
    [4, 0, 1, 2, 1, 4]
    >>> [get_stoichometry('H2SO4',e) for e in ['N', 'O', 'S', 'H']]
    [0, 4, 1, 2]
    """
    formula = formula.upper()
    element = element.upper()
    length  = len(element)
    idx     = formula.find(element)
    n = 0
    while idx > -1:
        try:
            n  += int(formula[idx+length])
            idx = formula.find(element,idx+length)
        except:
            n += 1
            idx = formula.find(element,idx+1) 
    return n


def parse_line16(s):
    """
    Return a numpy array of numbers parsed from a
    string of numbers located in every 16 chars.
    Note: Numbers may not have a space in between.
    >>> parse_line16(' 2.807326142D-08-7.923286750D-12 0.000000000D+00 3.329428940D+04 3.816278870D+01')
    array([  2.80732614e-08,  -7.92328675e-12,   0.00000000e+00,
             3.32942894e+04,   3.81627887e+01])
    """
    assert len(s) % 16 == 0, 'Given string for parse_line should have 16n chararacters, n={1,2,...}' 
    assert len(s) > 0, 'Given string for parse_line should have 16n chararacters, n={1,2,...}'
    
    n = len(s) / 16
    #replace fortran exponent D to E
    tmp = s.replace('D','E')
    nums = [0.]*n
    nums = np.zeros(n)
    for i in range(n):
        nums[i] = float(tmp[i*16:(i+1)*16])
    return nums


def get_comment_lines(tag,deltaH):
    """
    Returns 3 line string that includes the comment based on tag, deltaH and date.
    Based on Franklin Goldsmith's NASA_CKIN.py
    e.g.:
    line 1:!
    line 2:!DHf(0K) =    25.00 [kcal/mol], taken from SJK ANL0
    line 3:!Q(T) from CI+QC/cc-pVTZ by SJK on  18Apr2017
    
    TODO: To simplify, instead of a tag, comments could be given directly or
    it could be read from a database of tags, which can be accessed outside
    the code.
    """
    import datetime
    date = datetime.datetime.now().strftime("%d%b%Y")
    line1 = '!\n'
    if tag=='SJKB0':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL0\n"%(deltaH)
        line3 = '!Q(T) from B3LYP/6-311++G(d,p) by CFG on  ' + date + '\n'
    elif tag=='SJKB20':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL0\n"%(deltaH)
        line3 = '!Q(T) from B2PLYPD3/cc-pVTZ by SJK on  ' + date + '\n'
    elif tag=='SJKT0':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL0\n"%(deltaH)
        line3 = '!Q(T) from CCSD(T)/cc-pVTZ by SJK on  ' + date + '\n'
    elif tag=='SJKQ0':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL0\n"%(deltaH)
        line3 = '!Q(T) from CCSD(T)/cc-pVQZ by SJK on  ' + date + '\n'
    elif tag=='SJKCBS0':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL0\n"%(deltaH)
        line3 = '!Q(T) from CCSD(T)/CBSby SJK on  ' + date + '\n'
    elif tag=='SJKCIT0':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL0\n"%(deltaH)
        line3 = '!Q(T) from CI+QC/cc-pVTZ by SJK on  ' + date + '\n'
    elif tag=='SJKCIQ0':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL0\n"%(deltaH)
        line3 = '!Q(T) from CI+QC/cc-pVQZ by SJK on  ' + date + '\n'
    elif tag=='SJKCICBS0':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL0\n"%(deltaH)
        line3 = '!Q(T) from CI+QC/cc-pVQZ by SJK on  ' + date + '\n'
    elif tag=='SJKPT2T0':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL0\n"%(deltaH)
        line3 = '!Q(T) from CASPT2/cc-pVTZ by SJK on  ' + date + '\n'
    elif tag=='SJKPT2Q0':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL0\n"%(deltaH)
        line3 = '!Q(T) from CASPT2/cc-pVQZ by SJK on  ' + date + '\n'
    elif tag=='SJKB1':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL1\n"%(deltaH)
        line3 = '!Q(T) from B3LYP/6-311++G(d,p) by CFG on  ' + date + '\n'
    elif tag=='SJKT1':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL1\n"%(deltaH)
        line3 = '!Q(T) from CCSD(T)/cc-pVTZ by SJK on  ' + date + '\n'
    elif tag=='SJKQ1':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL1\n"%(deltaH)
        line3 = '!Q(T) from CCSD(T)/cc-pVQZ by SJK on  ' + date + '\n'
    elif tag=='SJKCBS1':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL1\n"%(deltaH)
        line3 = '!Q(T) from CCSD(T)/CBS by SJK on  ' + date + '\n'
    elif tag=='SJKCBSA1':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL1\n"%(deltaH)
        line3 = '!Q(T) from CCSD(T)/CBS + anh by SJK on  ' + date + '\n'
    elif tag=='SJKCIT1':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL1\n"%(deltaH)
        line3 = '!Q(T) from CI+QC/cc-pVTZ by SJK on  ' + date + '\n'
    elif tag=='SJKCIQ1':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL1\n"%(deltaH)
        line3 = '!Q(T) from CI+QC/cc-pVQZ by SJK on  ' + date + '\n'
    elif tag=='SJKPT2T1':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL1\n"%(deltaH)
        line3 = '!Q(T) from CASPT2/cc-pVTZ by SJK on  ' + date + '\n'
    elif tag=='SJKPT2Q1':
        line2 ="!DHf(0K) = %8.2F [kcal/mol], taken from SJK ANL1\n"%(deltaH)
        line3 = '!Q(T) from CASPT2/cc-pVQZ by SJK on  ' + date + '\n'
    else:
        line2 = '!{0}.\n'.format(date)
        line3 = '!{0}.n'.format(tag)
    return line1 + line2 + line3


def get_coefficients(c97filename):
    """
    Returns a string of 3 lines containing NASA polynomial
    coefficients in chemkin format
    *.c97 file:
    C2H3                                                                      
    3 201704 C   2.00H   3.00    0.00    0.00    0.00 0   27.0452200     296391.000
    100.000   200.000 2  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0        10698.000
    3.587567650D+00 3.894470300D-03 0.000000000D+00 0.000000000D+00 0.000000000D+00
    0.000000000D+00 0.000000000D+00 0.000000000D+00 3.437879620D+04 6.439970150D+00
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10698.000
    2.881119522D+00 4.825191250D-03 1.818030931D-05-2.828286454D-08 1.209654946D-11
    0.000000000D+00 0.000000000D+00 0.000000000D+00 3.446352960D+04 9.703788500D+00
    1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10698.000
    2.984627599D+00 1.078391826D-02-5.158601830D-06 1.200731137D-09-1.103701620D-13
    0.000000000D+00 0.000000000D+00 0.000000000D+00 3.423078010D+04 7.923373280D+00
    
        coefficients in chemkin format:
     2.98462760E+00 1.07839183E-02-5.15860183E-06 1.20073114E-09-1.10370162E-13    2
     3.42307801E+04 7.92337328E+00 2.88111952E+00 4.82519125E-03 1.81803093E-05    3
    -2.82828645E-08 1.20965495E-11 3.44635296E+04 9.70378850E+00                   4
    """
    with open(c97filename,'r') as f:
        lines = f.readlines()
    las  = np.zeros(7)
    has  = np.zeros(7)
    las[0:5] = parse_line16(lines[6][0:80])
    las[5:7] = parse_line16(lines[7][48:80])
    has[0:5] = parse_line16(lines[9][0:80])
    has[5:7] = parse_line16(lines[10][48:80])

    line2 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    2\n"%(has[0], has[1], has[2], has[3], has[4])
    line3 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    3\n"%(has[5], has[6], las[0], las[1], las[2])
    line4 = "% 15.8E% 15.8E% 15.8E% 15.8E                   4\n"%(las[3], las[4], las[5], las[6])
    return line2+line3+line4


def get_name_from_messpf(inputfile='pf.inp'):
    """
    Returns species formula by parsing messpf input file.
    input file:
    Temperature(step[K],size)        100.   30
    RelativeTemperatureIncrement            0.001
    Species   C2H3
    !      RRHO
    Geometry[angstrom]     5   ! CCSD(T)/cc-pVQZ
    C          0.0000000000        0.0229607853       -0.6194779279
    C          0.0000000000       -0.0846208019        0.6897727967
    H          0.0000000000        0.9977802761       -1.1068370147
    H          0.0000000000       -0.8465433753       -1.2673164024
    H          0.0000000000        0.5835275291        1.5364927734
    Core     RigidRotor
    SymmetryFactor       1
    End
    Frequencies[1/cm]      9   ! CCSD(T)/CBS; anh
    692    796    894    1019    1353    1580    2895    3023    3114
    ZeroEnergy[kcal/mol]            -34.41  ! ANL1
    ElectronicLevels[1/cm]        1
    0   2
    End
    """
    with open(inputfile,'r') as f:
        lines = f.readlines()
    for line in lines:
        if 'Species' in line:
            name = line.split()[-1]   
    return name


def write_chemkin_file(deltaH,tag,formula,filename):
    """
    Given formula string, tag string, deltaH float and a filename string,
    writes a file containing NASA polynomials in chemkin format:
    !
    !DHf(0K) =    25.00 [kcal/mol], taken from SJK ANL0
    !Q(T) from CI+QC/cc-pVTZ by SJK on  18Apr2017
    C2H3                    H   3C   2O   0N   0G   200.00   3000.00  1000.00      1
     2.98462760E+00 1.07839183E-02-5.15860183E-06 1.20073114E-09-1.10370162E-13    2
     3.42307801E+04 7.92337328E+00 2.88111952E+00 4.82519125E-03 1.81803093E-05    3
    -2.82828645E-08 1.20965495E-11 3.44635296E+04 9.70378850E+00                   4    
    """
    lines1to3 = get_comment_lines(tag, deltaH)
    nH = get_stoichometry(formula, 'H')
    nC = get_stoichometry(formula, 'C')
    nN = get_stoichometry(formula, 'N')
    nO = get_stoichometry(formula, 'O')
    line4 = "%s        H%4dC%4dO%4dN%4dG%9.2F%10.2F%9.2F      1\n"%(formula.ljust(16)[0:16], nH, nC, nO, nN, 200.0, 3000.0, 1000.0)
    lines5to7 = get_coefficients(formula+'.c97')
    with open(filename,'w') as f:
        f.write(lines1to3 + line4 +lines5to7)   
    return 0   


def write_thermp_input(formula,deltaH,enthalpyT=0.,breakT=1000.,filename='thermp.dat'):
    """
    Write thermp input file with given formula string, deltaH float
    e.g.
    1, 0
    Nwell, Nprod
    30
    nt
    71.82    0.                          //enthalpy at specified T (in kcal/mol), with T = 0 or 298 K
    C2H3                                    //name of species
    C  2                                    //composition in terms of
    H  3                                    //element_name  element_count
    **
    1000.                                   //temperature to break the fit
    """
    nH = get_stoichometry(formula, 'H')
    nC = get_stoichometry(formula, 'C')
    nN = get_stoichometry(formula, 'N')
    nO = get_stoichometry(formula, 'O')
    with open(filename,'w') as f:
        f.write('1, 0\n')
        f.write('Nwell, Nprod\n')
        f.write('30\n')
        f.write('nt\n')
        f.write('{0} {1}\n'.format(deltaH, enthalpyT))
        f.write('{0}\n'.format(formula))
        if nC > 0:
            f.write('C {0}\n'.format(nC))
        if nH > 0:
            f.write('H {0}\n'.format(nH))
        if nO > 0:
            f.write('O {0}\n'.format(nO))
        if nN > 0:
            f.write('N {0}\n'.format(nN))
        f.write('**\n')
        f.write('{0}'.format(breakT))
    return


def get_pf_input(mol,method,zpe,xyz,freqs):
    """
    Write input file for mess partition function program
    """
    import iotools as io
    optmethod  = method
    freqmethod = method
    tagmethod  = method
    sym = 1
    natom = len(mol.atoms)
    formula = mol.formula
    multiplicity = mol.spin
    inp  = 'Temperature(step[K],size)        100.   30\n'
    inp += 'RelativeTemperatureIncrement            0.001\n'
    inp += 'Species {0}\n'.format(formula)
    inp += 'RRHO\n'
    inp += 'Geometry[angstrom] {0} !{1}\n'.format(natom,optmethod)
    inp += ''.join(xyz.splitlines(True)[2:])
    inp += 'Core RigidRotor\n'
    inp += 'SymmetryFactor {0}\n'.format(sym)
    inp += 'End\n'
    inp += 'Frequencies[1/cm] {0} !{1}\n'.format(len(freqs),freqmethod)
    inp += ' '.join(freqs) + '\n'
    inp += 'ZeroEnergy[kcal/mol] {0} ! {1}\n'.format(zpe,tagmethod)
    inp += 'ElectronicLevels[1/cm]  1\n'
    inp += '0 {0}\n'.format(multiplicity) 
    inp += 'End\n'
    return inp
    
def run_pf(messpf='messpf',inputfile='pf.inp'):
    """
    Runs mess to generate partition function
    Requires an input file,i.e. pf.inp.
    '/tcghome/ygeorgi/fock/crossrate/bin/partition_function'
    Output is pf.log
    e.g.
    Temperature(step[K],size)        100.   30
    RelativeTemperatureIncrement            0.001
    Species   C2H3
    RRHO
    Geometry[angstrom]     5   ! CCSD(T)/cc-pVQZ
    C          0.0000000000        0.0229607853       -0.6194779279
    C          0.0000000000       -0.0846208019        0.6897727967
    H          0.0000000000        0.9977802761       -1.1068370147
    H          0.0000000000       -0.8465433753       -1.2673164024
    H          0.0000000000        0.5835275291        1.5364927734
    Core     RigidRotor
    SymmetryFactor       1
    End
    Frequencies[1/cm]      9   ! CCSD(T)/CBS; anh
    692    796    894    1019    1353    1580    2895    3023    3114
    ZeroEnergy[kcal/mol]            -34.41  ! ANL1
    ElectronicLevels[1/cm]        1
    0   2
    End    
    """
    import subprocess
    import iotools as io

    if io.check_exe(messpf):
        if io.check_file(inputfile,1):
            subprocess.call([messpf,inputfile])
        else:
            print "{0} input file does not exist".format(inputfile)

    else:
        print "{0} mess partitition function executable does not exist".format(messpf)
    return io.check_file('pf.log',1)

def run_thermp(thermp='thermp',thermpfile='thermp.dat',pffile='pf.dat'):
    """
    Runs thermp.exe
    Requires pf.dat and thermp.dat files as input
    linus
    /tcghome/sjk/gen/aux_me/therm/thermp.exe
    """
    import os
    import subprocess
    import iotools as io

    if not io.check_exe(thermp):
        print "{0} thermp executable does not exist".format(thermp)
        return
    pflog = pffile.replace('dat','log')
    if io.check_file(pflog) and not io.check_file(pffile):
        os.rename(pflog,pffile)
    if io.check_file(thermpfile,3):
        if io.check_file(pffile,3):
            subprocess.call([thermp])
        else:
            print "{0} file does not exist".format(pffile)
    return


def run_pac99(formula,pac99='pac99'):
    """
    Run pac99 for a given species name (formula)
    https://www.grc.nasa.gov/WWW/CEAWeb/readme_pac99.htm
    requires formula+'i97' and new.groups files
    TODO: Maybe add delete empty files
    linus
    pac99='/tcghome/sjk/gen/aux_me/therm/pac99.x'
    """
    from subprocess import Popen, PIPE
    import iotools as io

    c97file = formula +'.c97'
    i97file = formula +'.i97'
    o97file = formula +'.o97'
    if io.check_exe(pac99):
        if io.check_file(i97file):
            if io.check_file('new.groups'):
                p = Popen(pac99, stdin=PIPE)
                p.communicate(formula)
            else:
                print 'new.groups file is required to run pac99'
        else:        
            print '{0} file not found'.format(i97file)        

    else:
        print '{0} file not found'.format(pac99)        
    if io.check_file(c97file) and io.check_file(o97file):
        print "{0} {1} files are written".format(c97file,o97file)
    return


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
