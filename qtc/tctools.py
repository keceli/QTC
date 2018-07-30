#!/usr/bin/env python
import argparse
import datetime
import time
import subprocess
import os
from os.path import isfile
import iotools as io
import obtools as ob
import qctools as qc
import patools as pa
import logging
import unittools as ut
import math
import numpy as np
"""
Thermochemistry tools.
Requires:
Quantum chemistry code, NWChem, MOPAC, ...
MESS partition function code
PAC99
thermp
"""
__updated__ = "2017-12-15"


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
    n = 0
    if len(formula.split(element)) > 1:
          formula = formula.split(element)[1]
          formula = formula.split('C')[0]
          formula = formula.split('H')[0]
          formula = formula.split('N')[0]
          formula = formula.split('O')[0]
          if len(formula) > 0:
              n =  int(formula)
          else:
              n =  1
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
    nums = [0.] * n
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
        line3 = '!{0}.\n'.format(tag)
    return line1 + line2 + line3


def get_coefficients(c97text):
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
    lines = c97text.splitlines()
    las  = [0.] * 7
    has  = [0.] * 7
    msg  = None
    if len (lines) > 1:
        las[0:5] = parse_line16(lines[6][0:80])
        las[5:7] = parse_line16(lines[7][48:80])
        has[0:5] = parse_line16(lines[9][0:80])
        has[5:7] = parse_line16(lines[10][48:80])
    else:
        msg = 'pacc has failed'
    return las, has, msg


def get_coefficients_str(las,has):
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
    line2 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    2\n"%(has[0], has[1], has[2], has[3], has[4])
    line3 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    3\n"%(has[5], has[6], las[0], las[1], las[2])
    line4 = "% 15.8E% 15.8E% 15.8E% 15.8E                   4\n"%(las[3], las[4], las[5], las[6])
    return line2+line3+line4


def convert_chemkin2rmg(ckin):
    """
! OOCCC(O[O])C_m2      torsscan/m062x/cc-pvtz/gaussian
! deltaH(0) -32.7035898546 kcal/mol
! deltaH(298) -39.3632755502 kcal/mol
C4H9O4                  H   9C   4O   4N   0G   200.00   3000.00  1000.00      1
 1.57361381E+01 3.29113185E-02-1.59602886E-05 3.77476959E-09-3.52930905E-13    2
-2.63226479E+04-4.77325306E+01 4.60561015E-01 9.96060317E-02-1.28949962E-04    3
 9.08860280E-08-2.58936516E-11-2.34008733E+04 2.49356355E+01                   4
 
     line2 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    2\n"%(has[0], has[1], has[2], has[3], has[4])
    line3 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    3\n"%(has[5], has[6], las[0], las[1], las[2])
    line4 = "% 15.8E% 15.8E% 15.8E% 15.8E                   4\n"%(las[3], las[4], las[5], las[6])
    """
    lines = ckin.splitlines()
    nline = len(lines)
    has = [0]*7
    las = [0]*7
    if nline < 4:
        logging.error('Bad format for chemkin file at convert_chemkin2rmg')
    else:
        lines = lines[-4:]
        tlow, thigh, tmed = lines[0].split()[6:9]
        tlow = float(tlow)
        tmed = float(tmed)
        thigh = float(thigh)
        has[0] = float(lines[1][0:15])
        has[1] = float(lines[1][15:30])
        has[2] = float(lines[1][30:45])
        has[3] = float(lines[1][45:60])
        has[4] = float(lines[1][60:75])
        has[5] = float(lines[2][0:15])
        has[6] = float(lines[2][15:30])
        las[0] = float(lines[2][30:45])
        las[1] = float(lines[2][45:60])
        las[2] = float(lines[2][60:75])
        las[3] = float(lines[2][60:75])
        las[3] = float(lines[3][0:15])
        las[4] = float(lines[3][15:30])
        las[5] = float(lines[3][30:45])
        las[6] = float(lines[3][45:60])
    return get_rmg_polynomial(las, has, temps=[tlow,tmed,tmed,thigh])

    
def get_rmg_polynomial(las, has,temps=[200.,1000.,1000.,3000.]):
    """
    Return NASA polynomial as a dictionary in RMG format:
    
    NASA Polynomial, seven or nine coefficients, Tmax and Tmin = valid temperature range
    polynomials = [{'coeffs':[2.3443,0.00798042,-1.94779e-05,2.0157e-08,-7.37603e-12,-917.924,0.683002], 'Tmin':(200,'K'), 'Tmax':(1000,'K')},
                      {'coeffs':[2.93283,0.000826598,-1.46401e-07,1.54099e-11,-6.88796e-16,-813.056,-1.02432], 'Tmin':(1000,'K'), 'Tmax':(6000,'K')}]
                      Tmin = (200,'K')
                      Tmax = (6000,'K')
                      NASAPolynomial = {'polynomials':polynomials,'Tmin':Tmin,'Tmax':Tmax}
    """
    p = [{'coeffs': las, 'Tmin':(temps[0],'K'), 'Tmax':(temps[1],'K')},
         {'coeffs': has, 'Tmin':(temps[2],'K'), 'Tmax':(temps[3],'K')}]
    
    return {'polynomials': p, 'Tmin' : (min(temps), 'K'), 'Tmax' : (max(temps),'K')}


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


def get_chemkin_str(deltaH,tag,formula,filename):
    """
    Given formula string, tag string, deltaH float and a filename string,
    returns a string for NASA polynomials in chemkin format:
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
    lines5to7 = get_coefficients_str(formula+'.c97')

    return lines1to3 + line4 +lines5to7


def write_chemkin_file(slabel, qlabel, hof, hof298, formula, mid, las, has, filename):
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
    comments   = '! {} \t {}\n'.format(slabel, qlabel)
    comments += '! deltaH(0) {:.2f} kcal/mol\n'.format(hof)
    comments += '! deltaH(298) {:.2f} kcal/mol\n'.format(hof298)
    nH = get_stoichometry(formula, 'H')
    nC = get_stoichometry(formula, 'C')
    nN = get_stoichometry(formula, 'N')
    nO = get_stoichometry(formula, 'O')
    cformula = '{}_{}'.format(formula,str(mid))
    line4 = "%s        H%4dC%4dO%4dN%4dG%9.2F%10.2F%9.2F      1\n"%(cformula.ljust(16)[0:16], nH, nC, nO, nN, 200.0, 3000.0, 1000.0)
    lines5to7 = get_coefficients_str(las,has)
    s = comments + line4 +lines5to7
    io.write_file(s, filename)
    return s


def get_thermp_input(formula,deltaH,enthalpyT=0.,breakT=1000.):
    """
    Returns thermp input text as a string for a given formula string and deltaH (float)
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
    tmp = '1, 0\n'
    tmp +='Nwell, Nprod\n'
    tmp +='30\n'
    tmp +='nt\n'
    tmp +='{0} {1}\n'.format(deltaH, enthalpyT)
    tmp +='{0}\n'.format(formula)
    if nC > 0:
        tmp +='C {0}\n'.format(nC)
    if nH > 0:
        tmp +='H {0}\n'.format(nH)
    if nO > 0:
        tmp +='O {0}\n'.format(nO)
    if nN > 0:
        tmp +='N {0}\n'.format(nN)
    tmp +='**\n'
    tmp +='{0}'.format(breakT)
    return tmp


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
    if io.check_file(filename, 1):
        msg = 'Thermp input file "{0}" written.\n'.format(filename)
    else:
        msg = 'Failed writing thermp input file "{0}".\n'.format(filename)
    return msg


# def get_pf_input(mol,method,xyz,freqs,zpe=0., xmat=[], hindered=None):
#     """
#     Write input file for mess partition function program
#     Temperature(step[K],size)        100.   30
#     RelativeTemperatureIncrement            0.001
#     Species CH4
#     RRHO
#     Geometry[angstrom] 5 !pm3
#     C -0.0000 0.0000 0.0000
#     H 1.0870 0.0000 0.0000
#     H -0.3623 1.0249 0.0000
#     H -0.3623 -0.5124 0.8876
#     H -0.3623 -0.5124 -0.8876
#     Core RigidRotor
#     SymmetryFactor 1
#     End
#     Frequencies[1/cm] 9 !pm3
#     1362.16 1362.39 1362.48 1451.03 1451.06 3207.47 3207.48 3207.50 3311.01
#     ZeroEnergy[kcal/mol] 28.481 ! pm3
#     ElectronicLevels[1/cm]  1
#     0 1
#     End
#     """
#     optmethod  = method
#     freqmethod = method
#     tagmethod  = method
#     sym = 1
#     natom = len(mol.atoms)
#     formula = mol.formula
#     multiplicity = mol.spin
#     inp  = 'Temperature(step[K],size)        100.   30\n'
#     inp += 'RelativeTemperatureIncrement            0.001\n'
#     inp += 'Species {0}\n'.format(formula)
#     inp += 'RRHO\n'
#     inp += 'Geometry[angstrom] {0} !{1}\n'.format(natom,optmethod)
#     inp += ''.join(xyz.splitlines(True)[2:])
#     inp += '\nCore RigidRotor\n'
#     inp += 'SymmetryFactor {0}\n'.format(sym)
#     inp += 'End\n'
#     if hindered:
#         inp += hindered
#     inp += 'Frequencies[1/cm] {0} !{1}\n'.format(len(freqs),freqmethod)
#     inp += ' '.join(freqs) + '\n'
#     if len(xmat) > 0:
#         inp += ' Anharmonicities[1/cm]\n'
#         for i in range( len(xmat)):
#             for j in range(i+1):
#                 inp += '  ' + str(i) + ' ' + str(j) + ' ' + str(xmat[i,j]) + '\n'
#         inp += ' End\n'
#     inp += 'ZeroEnergy[kcal/mol] {0} ! {1}\n'.format(zpe,tagmethod)
#     inp += 'ElectronicLevels[1/cm]  1\n'
#     inp += '0 {0}\n'.format(multiplicity)
#     inp += 'End\n'
#     return inp


def get_messpf_input(mol,parameters):
    """
    TODO: Anharmonic frequencies 
    Write input file for mess partition function program
    AtomDistanceMin[angstrom] 0.6
    Temperature(step[K],size)        100.   30
    RelativeTemperatureIncrement            0.001
    Species CH4
    RRHO
    Geometry[angstrom] 5 !pm3
    C -0.0000 0.0000 0.0000
    H 1.0870 0.0000 0.0000
    H -0.3623 1.0249 0.0000
    H -0.3623 -0.5124 0.8876
    H -0.3623 -0.5124 -0.8876
    Core RigidRotor
    SymmetryFactor 1
    End
    Frequencies[1/cm] 9 !pm3
    1362.16 1362.39 1362.48 1451.03 1451.06 3207.47 3207.48 3207.50 3311.01
    ZeroEnergy[kcal/mol] 28.481 ! pm3
    ElectronicLevels[1/cm]  1
    0 1
    End
    
    For a single atom:    

    AtomDistanceMin[angstrom] 0.6
    Temperature(step[K],size)        100.   30
    RelativeTemperatureIncrement            0.001
    Species H
      Atom
         Mass[amu]     1

         ElectronicLevels[1/cm]     1
            0    2
       End
    """
    import unittools as ut
    import anharm
    natom = parameters['natom']
    label = parameters['qlabel']
    results = parameters['results']
    multiplicity = parameters['mult']
    xyz = results['xyz']
    sym = results['sym']
    natom = len(mol.atoms)
    formula = mol.formula
    #multiplicity = mol.spin
    freqs = []
    xmat = []
    rotconsts = []
    posfreqs = []
    zpe = 0
    rotdists = ''
    scale = 0
    scaletype = None
    vibrots = None
    emax = 500 #kcal/mol, not sure
    if 'azpve' in results:
        zpve = results['azpve']
    elif 'zpve' in results:
        zpve = results['zpve']
    else:
        zpve = 0.
    if 'pfreqs' in results:
        freqs = results['pfreqs']
    elif 'freqs' in results:
        freqs = results['freqs']    
        if len(freqs) > 0:
            for freq in freqs:
                if float(freq) > 0:
                    posfreqs.append(freq)
            if len(posfreqs) < len(freqs):
                logging.warning('Imaginary frequencies are ignored')
            results['freqs'] = posfreqs
            freqs = posfreqs
    if 'rotconsts' in results:
        rotconsts = results['rotconsts']  
    if 'vibrots' in results:
        vibrots = results['vibrots']    
    else:
        vibrots = None
    if 'rotdists' in results:
        rotdists = results['rotdists']  
    if 'xmat' in results:
        xmat = np.asarray(results['xmat'])
        if 'pfreqs' in results:
            freqs, fill, xmat, fill2, fill3, vibrots = anharm.main(results, vibrots)
        elif 'afreqs' in results:
            freqs = results['afreqs']
        #xmat = anharm.mess_x(xmat)
    if 'scale' in parameters:
        scale = parameters['scale']
    if 'scaletype' in parameters:
        scaletype = parameters['scaletype']

    coreIsMd = False
    if  'hindered potential' in results:
        if  'Core' in results['hindered potential']:
            coreIsMd = True

    #########HEADER
    inp  = 'AtomDistanceMin[angstrom] 0.6\n'
    inp += 'Temperature(step[K],size)        100.   30\n'
    inp += 'RelativeTemperatureIncrement            0.001\n'
    ###  BEGIN INPUT
    inp += 'Species {0}\n'.format(formula)
    if natom == 1:
        inp += 'Atom\n'
        inp += 'Mass[amu] {}\n'.format(ut.atommasses[formula])
        inp += 'End\n'
    else:
        ###  BEGIN RRHO
        inp += 'RRHO\n'
        inp += '  Geometry[angstrom] {0} !{1}\n\t  '.format(natom,label)
        inp += '\t  '.join(xyz.splitlines(True)[2:])
        inp += '\n  ZeroEnergy[kcal/mol] {0} ! {1}\n'.format(zpve,label)
        inp += '  ElectronicLevels[1/cm]  1\n'
        inp += '     0 {0}\n'.format(multiplicity)

        ###  BEGIN CORE
        coreline  = '   Core RigidRotor\n'
        coreline += '      ZeroPointEnergy[1/cm] {}\n'.format(zpe)
        hindlines = ''
        if 'hindered potential' in results: 
            if  coreIsMd:
                coreline  = '  Core MultiRotor\n'
                hindlines = '     {}'.format('     '.join(results['hindered potential' ].splitlines(True)[3:]))
            else:
                hindlines  = '   End\n'
                hindpot = results['hindered potential']
                if scale and scaletype:
                     if scaletype.startswith('h'):
                         hindpot = hindpot.split('Potential[kcal/mol]')
                         if len(hindpot) > 1:
                             for h, pot in enumerate(hindpot[1:]):
                                 pot, end = pot.split('End')
                                 num = pot.split()[0] 
                                 pot = pot.split()[1:]
                                 newpot = ' {}\n    '.format(num)
                                 for val in pot:
                                     newpot += '   {:.3f}'.format(float(scale) * float(val))
                                 hindpot[h+1] = newpot + '\n End' + end
                         hindpot = 'Potential[kcal/mol]'.join(hindpot)
                hindlines += '  {}'.format('  '.join(hindpot.splitlines(True)))
        inp += coreline
        inp += '      InterpolationEnergyMax[kcal/mol] {}\n'.format(emax)
        inp += '      SymmetryFactor {0}\n'.format(sym)
        if coreIsMd:
            inp += hindlines  ###   END RRHO
        #freqs
        if len(freqs) > 0:
            inp += '      Frequencies[1/cm] {0} !{1}\n'.format(len(freqs),label)
            inp += '      ' + ' '.join([str(x) for x in freqs]) + '\n'
        if scaletype and scale:
            if 'f' in scaletype:
                inp += '      FrequencyScalingFactor {:.4f}\n'.format(scale)
        #anharmonics
        if len(xmat) > 0:
            inp += '      Anharmonicities[1/cm]\n'
            for i in range( len(xmat)):
                inp += '\t\t' + ' '.join([str(xmat[i][j]) for j in range(i+1)]) + '\n'
        if not coreIsMd:
            if 'norot' in parameters:
                if not parameters['norot']:
                    #if len(rotconsts) > 0:
                    #    inp += '      RotationalConstants[1/cm] '
                    #    inp += ' '.join(rotconsts) + '\n'
                    if vibrots:
                        vibrots = vibrots.splitlines(True)
                        if len(freqs) == len(vibrots):
                            inp += '      RovibrationalCouplings[1/cm]\n'
                            inp += '\t   ' + '\t   '.join(vibrots) + '\n'
                        else:
                            logging.warning("Rotational Couplings length does not match freqs -- removed from pf.inp")
                    if len(rotdists) > 0:
                        inp += '      RotationalDistortion[1/cm]\n'
                        inp += '\t   ' + '\t   '.join(rotdists.splitlines(True)) + '\n'
                        inp += '      End\n' ###   END CORE
            inp += hindlines  ###   END RRHO
        if not 'hindered potential' in results:
            inp += '   End\n'
        inp += 'End\n' ###   END CORE
    return inp

def run_pf(messpf='messpf',inputfile='pf.inp'):
    """
    Runs mess to generate partition function
    Requires an input file,i.e. pf.inp.
    '/tcghome/ygeorgi/fock/crossrate/bin/partition_function'
    Output is input_prefix + ".log"
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
    msg = ''
    if io.check_exe(messpf):
        if io.check_file(inputfile,1):
            if io.check_exe(messpf):
                subprocess.call([messpf,inputfile])
                if io.check_file('pf.log',1):
                    msg += '{0} generated by mess.\n'.format(inputfile)
            else:
                msg += 'Mess executable not found {0}\n'.format(messpf)
                return msg
        else:
            msg += "{0} input file does not exist.\n".format(inputfile)

    else:
        msg += "{0} mess partitition function executable does not exist.\n".format(messpf)
    return msg

def run_thermp(thermpinput,thermpfile='thermp.dat',pffile='pf.out', thermpexe='thermp'):
    """
    Runs thermp.exe
    Requires pffile and thermpfile to be present
    linus
    /tcghome/sjk/gen/aux_me/therm/thermp.exe
    """
    import iotools as io
    msg = ''
    io.write_file(thermpinput, thermpfile)
    if not io.check_file(thermpfile,1):
        return "{0} file not found.\n".format(thermpfile)
    pfdat = pffile.replace('out','dat')
    if io.check_file(pffile):
        io.mv(pffile,pfdat)
    if io.check_file(pfdat,1):
        msg += io.execute(thermpexe)
    else:
        msg += "{0} file not found.\n".format(pffile)
    return msg


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
    msg = ''
    c97file = formula +'.c97'
    i97file = formula +'.i97'
    o97file = formula +'.o97'
    if io.check_exe(pac99):
        if io.check_file(i97file):
            if io.check_file('new.groups'):
                if io.check_exe(pac99):
                    p = Popen(pac99, stdin=PIPE)
                    p.communicate(formula)
                else:
                    msg += 'pac99 not found.\n'
                    return msg
            else:
                msg += 'new.groups file is required to run pac99.\n'
        else:
            msg += '{0} file not found.\n'.format(i97file)

    else:
        msg += '{0} file not found.\n'.format(pac99)
    if io.check_file(c97file) and io.check_file(o97file):
        msg += "{0} {1} files are written.\n".format(c97file,o97file)
    return msg


# def write_chemkin_polynomial(mol, xyz, freqs, deltaH,parameters, xmat=[], zpe=0.):
#     """
#     A driver to perform all operations to write NASA polynomial in
#     chemkin format. Assumes quantum chemistry calculation is performed.
#     """
#     messpfinput = 'pf.inp'
#     messpfoutput = 'pf.log'
#     name = mol.formula
#     tag = parameters['qcmethod']
#     inp = get_pf_input(mol, tag, xyz, freqs, xmat=xmat, zpe=0.)
#     io.write_file(inp, messpfinput)
#     msg = 'Running {0} to generate partition function.\n'.format(parameters['messpf'])
#     msg += io.execute([parameters['messpf'],messpfinput])
#     msg += 'Running thermp .\n'
#     inp = get_thermp_input(mol.formula, deltaH)
#     msg = run_thermp(inp, 'thermp.dat', messpfoutput, parameters['thermp']) 
#     msg += 'Running pac99.\n'
#     msg += run_pac99(name)
#     msg += 'Converting to chemkin format.\n'
#     chemkinfile = name + '.ckin'
#     msg += 'Writing chemkin file {0}.\n'.format(chemkinfile)
#     try:
#         msg += write_chemkin_file(deltaH, tag, name, chemkinfile)
#     except:
#         "Failed to write polynomials"
#     return msg


def write_chemkin_polynomial(mol, parameters):
    """
    A driver to perform all operations to write NASA polynomial in
    chemkin format. Assumes quantum chemistry calculation is performed.
    """
    messpfinput = 'pf.inp'
#   messpfoutput = 'pf.log'
#    messpfoutput = 'pf.out'
    messpfoutput = 'pf.out'
    formula = mol.formula
    qlabel = parameters['qlabel']
    slabel = parameters['slabel']
    mid    = parameters['mol_index']
    hof = parameters['results']['deltaH0']
    if parameters['skippf']:
        logging.debug('Skipping pf generation...')
    else:
        inp = get_messpf_input(mol, parameters)
        io.write_file(inp, messpfinput)
        logging.debug('Running {0} to generate partition function...'.format(parameters['messpf']))
        msg = io.execute([parameters['messpf'],messpfinput])
        logging.debug(msg)
    logging.debug('Running thermp...')
    inp = get_thermp_input(mol.formula, hof)
    msg = run_thermp(inp, 'thermp.dat', messpfoutput, parameters['thermp']) 
    logging.debug(msg)
    logging.debug('Running pac99...')
    msg = run_pac99(formula)
    hof298 = 0
    chemkininput = ''
    rmgpoly = {}
    logging.debug(msg)
    if io.check_file('thermp.out'):
        lines = io.read_file('thermp.out')
        hof298 = pa.get_298(lines)
        logging.info('delHf(298) = {0} kcal/mol'.format(hof298))
    else:
        logging.error('Failed to create thermp.out')
    c97file = formula + '.c97'
    if io.check_file(c97file):
        c97text  = io.read_file(c97file)
        las, has, msg = get_coefficients(c97text)
        if msg:
            logging.info(msg)
        logging.debug('Converting to chemkin format.')
        chemkinfile = formula + '.ckin'
        logging.debug('Writing chemkin file {0}.\n'.format(chemkinfile))    
        try:
            chemkininput = write_chemkin_file(slabel,qlabel, hof, hof298, formula, mid, las, has, chemkinfile)
            rmgpoly = get_rmg_polynomial(las,has)
        except:
            logging.error("Failed to write chemkin polynomials")
    else:
        logging.error('Cannot find {}.'.format(c97file))
    return hof298, chemkininput, rmgpoly


def get_heat_capacity(rmgpoly,T):
    """
    rmgpoly is a dictionary in the following format
    {u'Tmax': [3000.0, u'K'],
     u'Tmin': [200.0, u'K'],
     u'polynomials': [{u'Tmax': [1000.0, u'K'],
                       u'Tmin': [200.0, u'K'],
                       u'coeffs': [3.47200416,
                               0.0002877246707,
                               -1.014581759e-06,
                               1.344086535e-09,
                               -4.50222436e-13,
                               2548.182033,
                               1.600701434]},
                  {u'Tmax': [3000.0, u'K'],
                   u'Tmin': [1000.0, u'K'],
                   u'coeffs': [3.25257139,
                               0.000111709217,
                               4.54124982e-07,
                               -2.061316072e-10,
                               2.673718745e-14,
                               2658.216226,
                               2.95565086]}]}
                               
    Formulas for calculation:
    Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
    H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
    S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
    where a1, a2, a3, a4, a5, a6, and a7 are the numerical coefficients 
    supplied in NASA thermodynamic files. 
    The first 7 numbers starting on the second line of each species entry 
    (five of the second line and the first two of the third line) are the 
    seven coefficients (a1 through a7, respectively) for the high-temperature 
    range (above 1000 K, the upper boundary is specified on the first line of 
    the species entry). The following seven numbers are the coefficients 
    (a1 through a7, respectively) for the low-temperature range 
    (below 1000 K, the lower boundary is specified on the first line of the species entry).
    H in the above equation is defined as
    H(T) = Delta Hf(298) + [ H(T) - H(298) ]
    so that, in general, H(T) is not equal to 
    Delta Hf(T) and one needs to have the data for the reference elements to calculate Delta Hf(T).
    """
    alist = []
    cp = 0.
    for poly in rmgpoly['polynomials']:
        Tmax = poly['Tmax'][0]
        Tmin = poly['Tmin'][0]
        if T >= Tmin and T <= Tmax:
            alist = poly['coeffs']
    if len(alist) > 4:
        cp = alist[0] + alist[1]*T + alist[2]*T**2 + alist[3]*T**3 + alist[4]*T**4
        cp = cp * ut.Rinkcal
    else:
        logging.error['{} K is outside the temperature range of the given NASA polynomials [{},{}]'.format
                      (T,rmgpoly['Tmin'][0],rmgpoly['Tmax'][0])]
    return cp


def get_entropy(rmgpoly,T):
    """
    rmgpoly is a dictionary in the following format
    {u'Tmax': [3000.0, u'K'],
     u'Tmin': [200.0, u'K'],
     u'polynomials': [{u'Tmax': [1000.0, u'K'],
                       u'Tmin': [200.0, u'K'],
                       u'coeffs': [3.47200416,
                               0.0002877246707,
                               -1.014581759e-06,
                               1.344086535e-09,
                               -4.50222436e-13,
                               2548.182033,
                               1.600701434]},
                  {u'Tmax': [3000.0, u'K'],
                   u'Tmin': [1000.0, u'K'],
                   u'coeffs': [3.25257139,
                               0.000111709217,
                               4.54124982e-07,
                               -2.061316072e-10,
                               2.673718745e-14,
                               2658.216226,
                               2.95565086]}]}
                               
    Formulas for calculation:
    Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
    H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
    S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
    where a1, a2, a3, a4, a5, a6, and a7 are the numerical coefficients 
    supplied in NASA thermodynamic files. 
    The first 7 numbers starting on the second line of each species entry 
    (five of the second line and the first two of the third line) are the 
    seven coefficients (a1 through a7, respectively) for the high-temperature 
    range (above 1000 K, the upper boundary is specified on the first line of 
    the species entry). The following seven numbers are the coefficients 
    (a1 through a7, respectively) for the low-temperature range 
    (below 1000 K, the lower boundary is specified on the first line of the species entry).
    H in the above equation is defined as
    H(T) = Delta Hf(298) + [ H(T) - H(298) ]
    so that, in general, H(T) is not equal to 
    Delta Hf(T) and one needs to have the data for the reference elements to calculate Delta Hf(T).
    """
    alist = []
    S = 0.
    for poly in rmgpoly['polynomials']:
        Tmax = poly['Tmax'][0]
        Tmin = poly['Tmin'][0]
        if T >= Tmin and T <= Tmax:
            alist = poly['coeffs']
    if len(alist) > 6:
        S = alist[0] * math.log(T) + alist[1]*T + alist[2]*T**2/2. + alist[3]*T**3/3. + alist[4]*T**4/4 + alist[6]
        S = S * ut.Rinkcal
    else:
        logging.error['{} K is outside the temperature range of the given NASA polynomials [{},{}]'.format
                      (T,rmgpoly['Tmin'][0],rmgpoly['Tmax'][0])]
    return S


def get_enthalpy(rmgpoly,T):
    """
    rmgpoly is a dictionary in the following format
    {u'Tmax': [3000.0, u'K'],
     u'Tmin': [200.0, u'K'],
     u'polynomials': [{u'Tmax': [1000.0, u'K'],
                       u'Tmin': [200.0, u'K'],
                       u'coeffs': [3.47200416,
                               0.0002877246707,
                               -1.014581759e-06,
                               1.344086535e-09,
                               -4.50222436e-13,
                               2548.182033,
                               1.600701434]},
                  {u'Tmax': [3000.0, u'K'],
                   u'Tmin': [1000.0, u'K'],
                   u'coeffs': [3.25257139,
                               0.000111709217,
                               4.54124982e-07,
                               -2.061316072e-10,
                               2.673718745e-14,
                               2658.216226,
                               2.95565086]}]}
                               
    Formulas for calculation:
    H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
    where a1, a2, a3, a4, a5, a6, and a7 are the numerical coefficients 
    supplied in NASA thermodynamic files. 
    H in the above equation is defined as
    H(T) = Delta Hf(298) + [ H(T) - H(298) ]
    so that, in general, H(T) is not equal to 
    Delta Hf(T) and one needs to have the data for the reference elements to calculate Delta Hf(T).
    """
    alist = []
    H = 0.
    for poly in rmgpoly['polynomials']:
        Tmax = poly['Tmax'][0]
        Tmin = poly['Tmin'][0]
        if T >= Tmin and T <= Tmax:
            alist = poly['coeffs']
    if len(alist) > 6:
        H = alist[0] + alist[1]*T/2 + alist[2]*T**2/3. + alist[3]*T**3/4. + alist[4]*T**4/5 + alist[5]/T
        H = H * ut.Rinkcal * T / 1000.
    else:
        logging.error['{} K is outside the temperature range of the given NASA polynomials [{},{}]'.format
                      (T,rmgpoly['Tmin'][0],rmgpoly['Tmax'][0])]
    return H


def get_hindered_potential(s,report=False):
    """
 Rotor                             Hindered
 Group    4   5   6   7   8   9  10  11  12  13  14  15  16  17
 Axis             2           1
 Symmetry            1
 Potential[kcal/mol]           12
    0.00    1.52    2.20    2.14    1.83    1.58    2.03    3.78    6.35    7.46    5.81    1.99
 End
 Rotor                             Hindered
 Group    5   6   7   8   9  10  11  12  13  14  15  16  17
 Axis             4           2
 Symmetry            1
 Potential[kcal/mol]           12
    0.00    3.72    8.34    7.05    3.27    2.22    4.05    3.24    1.90    3.35    5.10    3.42
    """
    s = s.lower()
    lines = s.splitlines()
    pot = []
    for line in lines:
        if line.islower(): #Check if line has any letter
            pass
        elif line.strip():
            items = line.split()
            newpot = [float(x) for x in items]
            pot.append(newpot)
            if report:
                logging.info(line)
            if any(p < 0. for p in newpot):
                logging.error('Negative hindered potential detected')
                if not report:
                    break
    return pot


def get_new_groups():
    s = """
H2SO4
 2 J 9/77 H   2.00S   1.00O   4.00    0.00    0.00 0     98.07948    -735148.376
    300.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  5.0  0.0  0.0
  1.07256800D+00  4.37692260D-02 -5.53332430D-05  3.55182530D-08 -9.06773580D-12
                                                 -9.02597580D+04  1.89395820D+01
   1000.000  5000.000 5  0.0  1.0  2.0  3.0  4.0  5.0  0.0  0.0
  1.08895320D+01  7.50041780D-03 -2.92104780D-06  5.25955130D-10 -3.57894150D-14
                                                 -9.24713640D+04 -2.94047820D+01
H2S
 2 J 6/77 H   2.00S   1.00    0.00    0.00    0.00 0     34.08188     -20502.254
    300.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  5.0  0.0  0.0
  3.93234760D+00 -5.02609050D-04  4.59284730D-06 -3.18072140D-09  6.64975610D-13
                                                 -3.65053590D+03  2.31579050D+00
   1000.000  5000.000 5  0.0  1.0  2.0  3.0  4.0  5.0  0.0  0.0
  2.74521990D+00  4.04346070D-03 -1.53845100D-06  2.75202490D-10 -1.85920950D-14
                                                 -3.41994440D+03  8.05467450D+00
SH
 2 J 6/77 S   1.00H   1.00    0.00    0.00    0.00 0     33.07394     139332.329
    300.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  5.0  0.0  0.0
  4.44203220D+00 -2.43591970D-03  1.90645760D-06  9.91666300D-10 -9.57407620D-13
                                                  1.55232580D+04 -1.14449035D+00
   1000.000  5000.000 5  0.0  1.0  2.0  3.0  4.0  5.0  0.0  0.0
  3.00145370D+00  1.33949570D-03 -4.67896630D-07  7.88040150D-11 -5.02804530D-15
                                                  1.59053200D+04  6.28462715D+00
SO3
 2 J 9/65 S   1.00O   3.00    0.00    0.00    0.00 0     80.06420    -395752.673
    300.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  5.0  0.0  0.0
  2.57803850D+00  1.45563350D-02 -9.17641730D-06 -7.92030220D-10  1.97094730D-12
                                                 -4.89317530D+04  1.22651384D+01
   1000.000  5000.000 5  0.0  1.0  2.0  3.0  4.0  5.0  0.0  0.0
  7.07573760D+00  3.17633870D-03 -1.35357600D-06  2.56309120D-10 -1.79360440D-14
                                                 -5.02113760D+04 -1.11875176D+01
SO2
 2 J 6/61 S   1.00O   2.00    0.00    0.00    0.00 0     64.06480    -296834.548
    300.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  5.0  0.0  0.0
  3.26653380D+00  5.32379020D-03  6.84375520D-07 -5.28100470D-09  2.55904540D-12
                                                 -3.69081480D+04  9.66465108D+00
   1000.000  5000.000 5  0.0  1.0  2.0  3.0  4.0  5.0  0.0  0.0
  5.24513640D+00  1.97042040D-03 -8.03757690D-07  1.51499690D-10 -1.05580040D-14
                                                 -3.75582270D+04 -1.07404892D+00
AR                Argon. NSRDS-NBS 35, vl, 1971. Temperature cutoff.
 2 L 6/88 AR  1.00    0.00    0.00    0.00    0.00 0     39.94800          0.000
    200.000  1000.000 1  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0         6197.428
 2.500000000D+00 0.000000000D+00 0.000000000D+00 0.000000000D+00 0.000000000D+00
 0.000000000D+00 0.000000000D+00 0.000000000D+00-7.453750000D+02 4.379674910D+00
   1000.000  6000.000 1  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0         6197.428
 2.500000000D+00 0.000000000D+00 0.000000000D+00 0.000000000D+00 0.000000000D+00
 0.000000000D+00 0.000000000D+00 0.000000000D+00-7.453750000D+02 4.379674910D+00
N2                Nitrogen. GLUSHKO ET.AL. v1, pt2, p207, 1978.
 2 TPIS78 N   2.00    0.00    0.00    0.00    0.00 0     28.01348          0.000
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0         8670.104
 3.531005280D+00-1.236609870D-04-5.029994360D-07 2.435306118D-09-1.408812347D-12
 0.000000000D+00 0.000000000D+00 0.000000000D+00-1.046976280D+03 2.967474680D+00
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0         8670.104
 2.952576373D+00 1.396900385D-03-4.926315980D-07 7.860101870D-11-4.607551990D-15
 0.000000000D+00 0.000000000D+00 0.000000000D+00-9.239486900D+02 5.871891890D+00
O2                Oxygen. Gurvich et al. v1, pt 2, p9, 1989.
 2 TPIS89 O   2.00    0.00    0.00    0.00    0.00 0     31.99880          0.000
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0         8680.104
 3.782456360D+00-2.996734156D-03 9.847302010D-06-9.681295090D-09 3.243728370D-12
 0.000000000D+00 0.000000000D+00 0.000000000D+00-1.063943564D+03 3.657675730D+00
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0         8680.104
 3.660960650D+00 6.563658110D-04-1.411496268D-07 2.057979356D-11-1.299134362D-15
 0.000000000D+00 0.000000000D+00 0.000000000D+00-1.215977179D+03 3.415362790D+00
CO2               Props & Hf298: TPIS v2,pt1,1991,p27.
 2 L 7/88 C   1.00O   2.00    0.00    0.00    0.00 0     44.00980    -393510.000
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0         9365.469
 2.356773524D+00 8.984596770D-03-7.123562690D-06 2.459190224D-09-1.436995477D-13
 0.000000000D+00 0.000000000D+00 0.000000000D+00-4.837196970D+04 9.901052220D+00
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0         9365.469
 4.636594930D+00 2.741319907D-03-9.958285310D-07 1.603730114D-10-9.161034680D-15
 0.000000000D+00 0.000000000D+00 0.000000000D+00-4.902493410D+04-1.935348550D+00
CA
 2 BEN76  C   1.00    0.00    0.00    0.00    0.00 0     12.01100      17210.010
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  6.20131064d-01  6.55685438d-03 -8.40971939d-06  5.34237173d-09 -1.34221334d-12
                                                  1.67980617d+04 -2.13964424d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  1.70014134d+00  2.13036557d-03 -1.46491180d-06  4.57867726d-10 -5.30959813d-14
                                                  1.65796625d+04 -7.34011722d+00
CBC
 2 BEN76  C   1.00    0.00    0.00    0.00    0.00 0     12.01100       2772.724
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  1.10379898d+00 -1.40207754d-03  1.01721889d-05 -1.04915786d-08  3.35530569d-12
                                                  2.43522265d+03 -1.01067689d+01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  1.18859776d+00  2.62679100d-03 -1.40397450d-06  3.62315999d-10 -3.62315999d-14
                                                  2.15954456d+03 -1.17034055d+01
CBCB
 2 S&F85  C   1.00H   0.00    0.00    0.00    0.00 0     12.01100          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -9.07225320d-01  1.21274570d-02 -1.59921140d-05  1.06772290d-08 -2.88771260d-12
                                                  2.34891290d+03 -2.17231890d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  1.64937250d+00  2.43505820d-03 -1.40110120d-06  3.74624350d-10 -3.82283850d-14
                                                  1.78059830d+03 -1.47139470d+01
CBCD
 2 S&F85  C   1.00H   0.00    0.00    0.00    0.00 0     12.01100          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  1.89361290d+00 -3.63320180d-03  1.33476670d-05 -1.32062840d-08  4.38875210d-12
                                                  2.36121080d+03 -1.41160640d+01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  6.57951850d-01  3.79153180d-03 -2.18052300d-06  5.80353650d-10 -5.87829120d-14
                                                  2.50341600d+03 -8.72472110d+00
CBCT
 2 S&F85  C   1.00H   0.00    0.00    0.00    0.00 0     12.01100          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  1.89361290d+00 -3.63320180d-03  1.33476670d-05 -1.32062840d-08  4.38875210d-12
                                                  2.36121080d+03 -1.58018030d+01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  6.57951850d-01  3.79153180d-03 -2.18052300d-06  5.80353650d-10 -5.87829120d-14
                                                  2.50341600d+03 -1.04104590d+01
CBH
 2 S&F85  C   1.00H   1.00    0.00    0.00    0.00 0     13.01894          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -8.59079180d-01  1.01575080d-02 -6.05790130d-06  1.11817290d-10  8.76799660d-13
                                                  1.51812920d+03  7.93471810d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  8.12920500d-01  5.70326930d-03 -2.94688640d-06  7.32625720d-10 -7.11722860d-14
                                                  1.07063610d+03 -6.86258470d-01
CDC2
 2 BEN76  C   1.00    0.00    0.00    0.00    0.00 0     12.01100       5203.260
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  1.57567141d-02  1.20990725d-02 -2.39398157d-05  2.39017773d-08 -9.01711547d-12
                                                  4.82932602d+03 -9.21726687d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  2.25642353d+00  1.36203978d-03 -7.27986777d-07  1.87867555d-10 -1.87867555d-14
                                                  4.34871092d+03 -1.99090723d+01
CDCBC
 2 BEN76  C   1.00    0.00    0.00    0.00    0.00 0     12.01100       4347.792
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -2.15050677d+00  2.56218701d-02 -4.82075840d-05  4.17875044d-08 -1.37198166d-11
                                                  4.19996194d+03 -9.32804085d-01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  3.04446083d+00  4.86442777d-04 -2.59995277d-07  6.70955554d-11 -6.70955554d-15
                                                  3.27765919d+03 -2.51782099d+01
CDCDC
 2 BEN76  C   1.00    0.00    0.00    0.00    0.00 0     12.01100       4468.564
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -2.15050677d+00  2.56218701d-02 -4.82075840d-05  4.17875044d-08 -1.37198166d-11
                                                  4.32073394d+03 -9.32804085d-01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  3.04446083d+00  4.86442777d-04 -2.59995277d-07  6.70955554d-11 -6.70955554d-15
                                                  3.39843119d+03 -2.51782099d+01
CDHC
 2 BEN76  C   1.00H   1.00    0.00    0.00    0.00 0     13.01894       4322.631
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  3.87658435d-01  6.91462635d-03 -4.92473068d-06  3.06564147d-09 -1.19103098d-12
                                                  3.93773106d+03 -6.55229555d-02
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  1.76953318d+00  3.63735034d-03 -1.32392850d-06  1.71399542d-10 -2.18997389d-15
                                                  3.48001920d+03 -7.46676354d+00
CDHCB
 2 S&F85  C   1.00H   1.00    0.00    0.00    0.00 0     13.01894          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -1.77479380d+00  2.03672900d-02 -2.92063140d-05  2.13900470d-08 -6.19476560d-12
                                                  3.25431500d+03  8.37139370d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  2.16078460d+00  3.89973620d-03 -1.87505280d-06  4.40575310d-10 -4.10819930d-14
                                                  2.44872450d+03 -1.05679580d+01
CDHCD
 2 S&F85  C   1.00H   1.00    0.00    0.00    0.00 0     13.01894          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -1.77479380d+00  2.03672900d-02 -2.92063140d-05  2.13900470d-08 -6.19476560d-12
                                                  3.25431500d+03  8.37139370d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  2.16078460d+00  3.89973620d-03 -1.87505280d-06  4.40575310d-10 -4.10819930d-14
                                                  2.44872440d+03 -1.05679580d+01
CDHCT
 2 S&F85  C   1.00H   1.00    0.00    0.00    0.00 0     13.01894          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -1.77479380d+00  2.03672900d-02 -2.92063140d-05  2.13900470d-08 -6.19476560d-12
                                                  3.25431500d+03  9.20168260d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  2.16078460d+00  3.89973620d-03 -1.87505280d-06  4.40575310d-10 -4.10819930d-14
                                                  2.44872450d+03 -9.73766950d+00
CDH2
 2 S&F85  C   1.00H   2.00    0.00    0.00    0.00 0     14.02688          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  7.08636360d-01  5.71738370d-03  3.97432860d-06 -8.14882140d-09  3.39759220d-12
                                                  2.66405290d+03  8.03997270d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  7.62035270d-01  7.90072810d-03 -3.83366760d-06  9.04921890d-10 -8.42780610d-14
                                                  2.55458540d+03  7.24431290d+00
CHC3
 2 BEN76  C   1.00H   1.00    0.00    0.00    0.00 0     13.01894       -956.112
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -1.21199942d+00  1.57135839d-02 -1.58129928d-05  8.16308832d-09 -1.79427286d-12
                                                 -1.16875169d+03 -3.21908314d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  4.21180248d+00  8.18242540d-05  1.33950798d-06 -6.69862359d-10  9.41347774d-14
                                                 -2.66361789d+03 -3.11576551d+01
CHCBC2
 2 BEN76  C   1.00H   1.00    0.00    0.00    0.00 0     13.01894       -493.152
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -2.70626696d+00  2.55442809d-02 -3.35585641d-05  2.09055818d-08 -5.02705367d-12
                                                 -5.64094599d+02  3.00591586d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  3.92458678d+00  2.09170394d-03 -1.11797969d-06  2.88510888d-10 -2.88510888d-14
                                                 -2.12756082d+03 -2.99433080d+01
CHCDC2
 2 BEN76  C   1.00H   1.00    0.00    0.00    0.00 0     13.01894       -744.761
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -2.69178809d+00  2.27174202d-02 -2.67597650d-05  1.50516772d-08 -3.18914467d-12
                                                 -7.43740880d+02  3.74363064d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  2.00129268d+00  5.30222627d-03 -2.83394852d-06  7.31341554d-10 -7.31341554d-14
                                                 -1.74761503d+03 -1.92282945d+01
CHCTC2
 2 BEN76  C   1.00H   1.00    0.00    0.00    0.00 0     13.01894       -865.533
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -2.29967144d+00  2.07315703d-02 -2.60278088d-05  1.75696487d-08 -4.99182528d-12
                                                 -9.03749304d+02  2.30198068d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  1.48247630d+00  5.93460188d-03 -3.17194238d-06  8.18565776d-10 -8.18565777d-14
                                                 -1.70025797d+03 -1.60989330d+01
CH2C2
 2 BEN76  C   1.00H   2.00    0.00    0.00    0.00 0     14.02688      -2480.858
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  2.95576901d-01  8.26548665d-03  2.02730929d-06 -8.22251499d-09  3.84394234d-12
                                                 -2.93983603d+03  5.66809199d-01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  8.95359335d-01  8.82630650d-03 -4.51610338d-06  1.11218836d-09 -1.07950628d-13
                                                 -3.18218807d+03 -2.98904913d+00
CH2CBC
 2 BEN76  C   1.00H   2.00    0.00    0.00    0.00 0     14.02688      -2445.633
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -1.61079486d+00  2.13020824d-02 -2.44229173d-05  1.42423100d-08 -3.20026830d-12
                                                 -2.72304502d+03  8.49250374d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  2.81096829d+00  5.93460188d-03 -3.17194238d-06  8.18565777d-10 -8.18565777d-14
                                                 -3.81247249d+03 -1.36149824d+01
CH2CBD
 2 BEN76  C   1.00H   2.00    0.00    0.00    0.00 0     14.02688      -2158.799
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -2.99489148d+00  2.49586869d-02 -2.77613624d-05  1.47398773d-08 -2.60172081d-12
                                                 -2.15783856d+03  1.58638825d+01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  1.17752700d+00  8.75596999d-03 -4.67991499d-06  1.20772000d-09 -1.20772000d-13
                                                 -3.03586486d+03 -4.40568818d+00
CH2CD2
 2 BEN76  C   1.00H   2.00    0.00    0.00    0.00 0     14.02688      -2158.799
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -2.99489148d+00  2.49586869d-02 -2.77613624d-05  1.47398773d-08 -2.60172081d-12
                                                 -2.15783856d+03  1.58638825d+01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  1.17752700d+00  8.75596999d-03 -4.67991499d-06  1.20772000d-09 -1.20772000d-13
                                                 -3.03586486d+03 -4.40568818d+00
CH2CDC
 2 BEN76  C   1.00H   2.00    0.00    0.00    0.00 0     14.02688      -2395.311
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -9.84166104d-01  1.41467063d-02 -7.15413971d-06 -2.15080256d-09  2.42261065d-12
                                                 -2.66434595d+03  6.65325780d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  8.87674198d-01  9.14512421d-03 -4.88791121d-06  1.26139644d-09 -1.26139644d-13
                                                 -3.13410437d+03 -2.90870113d+00
CH2CTC
 2 BEN76  C   1.00H   2.00    0.00    0.00    0.00 0     14.02688      -2380.215
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -7.61851269d-01  1.30343429d-02 -7.51763250d-06  4.44629802d-10  9.34733367d-13
                                                 -2.66730666d+03  5.96602298d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  3.68857816d-01  9.77749982d-03 -5.22590507d-06  1.34862066d-09 -1.34862066d-13
                                                 -2.94558197d+03  2.32409280d-01
CH3C
 2 BEN76  C   1.00H   3.00    0.00    0.00    0.00 0     15.03482      -5132.810
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  9.67091211d-01  4.54272496d-03  1.40931220d-05 -2.03529587d-08  8.18255263d-12
                                                 -5.71121159d+03  7.97556073d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
 -6.27511183d-01  1.36690420d-02 -7.30586729d-06  1.88538511d-09 -1.88538511d-13
                                                 -5.43213903d+03  1.52438529d+01
CTC
 2 BEN76  C   1.00    0.00    0.00    0.00    0.00 0     12.01100      13863.619
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  8.02437417d-01  3.49837285d-03 -4.15397270d-06  4.07988802d-09 -1.75083632d-12
                                                  1.34983448d+04 -2.26753348d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  1.24244195d+00  2.09170394d-03 -1.11797969d-06  2.88510888d-10 -2.88510888d-14
                                                  1.33531243d+04 -4.58500863d+00
CTCB
 2 S&F85  C   1.00H   0.00    0.00    0.00    0.00 0     12.01100          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
 -3.49384520d+00  2.72321320d-02 -4.76891040d-05  3.86559630d-08 -1.19225380d-11
                                                  1.33155360d+04  1.68245410d+01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  2.28560130d+00  8.22574900d-04 -4.14361390d-07  1.03588260d-10 -1.03983680d-14
                                                  1.22382860d+04 -1.04535180d+01
CTCD
 2 S&F85  C   1.00H   0.00    0.00    0.00    0.00 0     12.01100          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  7.72587840d-01  1.44159610d-03  2.24535910d-06 -3.27337130d-09  1.18176800d-12
                                                  1.38820450d+04 -1.66930970d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  4.49136330d-01  3.36057940d-03 -1.88820530d-06  4.95862840d-10 -4.98295890d-14
                                                  1.39278700d+04 -2.35699140d-01
CTCT
 2 S&F85  C   1.00H   0.00    0.00    0.00    0.00 0     12.01100          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  9.09050340d-03  9.68282400d-03 -1.61930660d-05  1.34477210d-08 -4.31250880d-12
                                                  1.25675030d+04  6.29563560d-01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  1.61375010d+00  1.75286410d-03 -9.50741270d-07  2.44108940d-10 -2.42282270d-14
                                                  1.22902950d+04 -6.81710090d+00
CTH
 2 S&F85  C   1.00H   1.00    0.00    0.00    0.00 0     13.01894          0.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0
  3.22062990d-01  1.23444940d-02 -1.94881630d-05  1.58382730d-08 -4.95581450d-12
                                                  1.30498420d+04  7.64972930d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0
  2.08408360d+00  3.00946820d-03 -1.27530500d-06  2.64569910d-10 -2.18420530d-14
                                                  1.27910130d+04 -3.35539540d-01
HVIN              C2H4 - C2H3
 2 L 2/91 H   1.00    0.00    0.00    0.00    0.00 0      1.00794    -625071.953
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0            0.000
  7.46735360d-01 -9.08531026d-03  3.11780593d-05 -3.33930685d-08  1.22733472d-11
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -7.52125326d+04 -2.33376993d+00
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0            0.000
 -3.59257940d-01  2.99014152d-03 -1.07409461d-06  1.73348003d-10 -1.03738147d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -7.53284072d+04  1.23891111d+00
HVINS             C2H4 - C2H3 + 8 kcal correction on H.
 2 L 2/91 H   1.00    0.00    0.00    0.00    0.00 0      1.00794    -591599.953
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0            0.000
  7.46735360d-01 -9.08531026d-03  3.11780593d-05 -3.33930685d-08  1.22733472d-11
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -7.11867993d+04 -2.33376993d+00
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0            0.000
 -3.59257940d-01  2.99014152d-03 -1.07409461d-06  1.73348003d-10 -1.03738147d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -7.13026738d+04  1.23891111d+00
HPHEN             C6H6 - C6H5
 2 L 1/91 H   1.00    0.00    0.00    0.00    0.00 0      1.00794    -254320.000
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0            0.000
 -2.06302352d-01 -8.15243200d-04  1.43769315d-05 -1.95955059d-08  8.17474280d-12
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -3.05819013d+04  1.10333770d+00
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0            0.000
  3.06920100d-01  2.32195120d-03 -8.16396860d-07  1.29838420d-10 -7.68981960d-15
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -3.08941277d+04 -2.51191430d+00
HC2H              C2H2 - C2H1
 2 L 3/91 H   1.00    0.00    0.00    0.00    0.00 0      1.00794    -331613.472
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0            0.000
 -3.22796852d+00  1.87214852d-02 -3.02272137d-05  2.39352994d-08 -7.26880375d-12
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -3.95302119d+04  1.38493147d+01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0            0.000
 -5.51834000d-02  3.11952519d-03 -1.41953706d-06  3.15909714d-10 -2.79158072d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -4.00479060d+04 -8.06478801d-01
C6H5              PHENYL RADICAL. NASA TM 83800, 1985. TRC 10/89.
 2 L 1/91 C   6.00H   5.00    0.00    0.00    0.00 0     77.10570     337200.000
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        14004.789
  7.09733118d-01  1.93298588d-02  5.94082169d-05 -9.85091151d-08  4.25428989d-11
  0.00000000d+00  0.00000000d+00  0.00000000d+00  3.91345677d+04  2.30298910d+01
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        14004.789
  1.07700130d+01  1.83851527d-02 -6.70001899d-06  1.09228975d-09 -6.58438624d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00  3.52041168d+04 -3.50135190d+01
C6H6              Benzene. NASA TM 83800, 1985. TRC 10/86.
 2 L 1/91 C   6.00H   6.00    0.00    0.00    0.00 0     78.11364      82880.000
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        14194.792
  5.03430766d-01  1.85146156d-02  7.37851484d-05 -1.18104621d-07  5.07176417d-11
  0.00000000d+00  0.00000000d+00  0.00000000d+00  8.55266640d+03  2.41332287d+01
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        14194.792
  1.10769331d+01  2.07071039d-02 -7.51641585d-06  1.22212817d-09 -7.35336820d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00  4.30998914d+03 -3.75254333d+01
C2H3              VINYL RAD.  Ervin, JACS 1990, v112, p5750. Taylor;Ames.
 2 L 2/91 C   2.00H   3.00    0.00    0.00    0.00 0     27.04582     677571.953
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10575.049
  3.21246076d+00  1.51485596d-03  2.59207153d-05 -3.57655130d-08  1.47149755d-11
  0.00000000d+00  0.00000000d+00  0.00000000d+00  8.03023088d+04  7.81741050d+00
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10575.049
  4.35099402d+00  7.49338358d-03 -2.64318675d-06  4.21293814d-10 -2.49901527d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00  7.95971039d+04 -1.21152748d-01
C2H4              ETHYLENE.  VARIOUS REFS.
 2 L 1/07 C   2.00H   4.00    0.00    0.00    0.00 0     28.05376      52500.000
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10518.688
  3.95919612d+00 -7.57045430d-03  5.70987746d-05 -6.91585815d-08  2.69883227d-11
  0.00000000d+00  0.00000000d+00  0.00000000d+00  5.08977621d+03  5.48364057d+00
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10518.688
  3.99173608d+00  1.04835251d-02 -3.71728136d-06  5.94641817d-10 -3.53639674d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00  4.26869674d+03  1.11775836d+00
C2H2              GLUSHKO CONSTANTS,1978. DelH TRC tables, 10/31/88.
 2 L 1/91 C   2.00H   2.00    0.00    0.00    0.00 0     26.03788     228200.000
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10012.261
  5.72658274d-01  2.51337936d-02 -4.00165415d-05  3.27263088d-08 -1.02964640d-11
  0.00000000d+00  0.00000000d+00  0.00000000d+00  2.64518701d+04  1.56116943d+01
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10012.261
  4.12608035d+00  6.10531164d-03 -2.61558836d-06  5.50068670d-10 -4.61171092d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00  2.59396953d+04 -4.43643704d-01
C2H               ETHYNYL. Delh: Ervin,JACS v112,1990. Jacox, 1988.
 2 L 1/91 C   2.00H   1.00    0.00    0.00    0.00 0     25.02994     559813.472
    298.150  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10454.472
  3.80062679d+00  6.41230844d-03 -9.78932779d-06  8.79100936d-09 -3.02766025d-12
  0.00000000d+00  0.00000000d+00  0.00000000d+00  6.59820820d+04  1.76237958d+00
   1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10454.472
  4.18126375d+00  2.98578645d-03 -1.19605130d-06  2.34158956d-10 -1.82013020d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00  6.59876013d+04  3.62835097d-01
CRESOL            EQL MIXTURE. KUDCHADKER ET AL JPCRD,V7,N2,P417,1978.
 2 L 6/87 C   7.00H   8.00O   1.00    0.00    0.00 0    108.13992    -132298.000
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        21109.369
  4.22633160d-01  4.55511152d-02  3.20141892d-05 -8.11240449d-08  3.76665388d-11
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -1.82026239d+04  2.60327123d+01
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        21109.369
  1.59327193d+01  2.70116915d-02 -9.94520997d-06  1.62975034d-09 -9.85197125d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -2.35919689d+04 -5.97313776d+01
C6H5OH            PHENOL. NASA TM 83800, JAN 1985.
 2 L 6/90 C   6.00H   6.00O   1.00    0.00    0.00 0     94.11304     -96399.000
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        17496.655
 -2.91043722d-01  4.08566827d-02  2.42825579d-05 -7.14481301d-08  3.46006031d-11
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -1.34129213d+04  2.68748641d+01
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        17496.655
  1.41556495d+01  1.99344395d-02 -7.18195552d-06  1.16224855d-09 -6.97121653d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -1.81288534d+04 -5.18007382d+01
C6H5O             PHENOXY RADICAL. NASA TM-83800, 1985.
 2 L 6/90 C   6.00H   5.00O   1.00    0.00    0.00 0     93.10510      47697.600
    200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        16208.083
  7.75754295d-02  3.30579420d-02  3.60341330d-05 -7.93148251d-08  3.64321927d-11
  0.00000000d+00  0.00000000d+00  0.00000000d+00  4.06540008d+03  2.57601045d+01
   1000.000  6000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        16208.083
  1.31516170d+01  1.90163347d-02 -6.94688180d-06  1.13441103d-09 -6.84628699d-14
  0.00000000d+00  0.00000000d+00  0.00000000d+00 -4.73010834d+02 -4.67113087d+01
"""
    return s

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
