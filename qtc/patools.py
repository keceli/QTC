#!usr/bin/python

import re
import iotools as io
import unittools as ut
import numpy as np
import logging
"""
Module for parsing logfiles

How to use parser on lines from logfile:
   import patools as pa
   pa.program_key(lines)
   or
   pa.key(lines) (to automatically determine program)

Program: 
   gaussian
   molpro
Key:
   basisset
   method
   energy i   (can specify method or just will use last method in logfile)
   freqs
   zmat
   xyz      (gaussian only thus far)
"""


def get_prog(lines):
   
   if 'Gaussian' in lines:
       return 'gaussian'
   elif 'molpro' in lines or 'MOLPRO' in lines:
       return 'molpro'
   elif 'Q-Chem' in lines:
       return 'qchem'
   return

################################################
############     Gaussian PARSER    ############
################################################
def gaussian_islinear(s):
    """
    Returns true if the molecule is linear for the given log.
    """
    if "Linear Molecule" in s or gaussian_natom(s) == 2:
        return True
    else:
        return False
    
def gaussian_natom(s):
    """
    NAtoms=     30 NQM=       30 NQMF=       0 NMMI=      0 NMMIF=      0
    """
    import iotools as io
    if type(s) == str:
        lines = s.splitlines()
    keyword = 'NAtoms='
    n = io.get_line_number(keyword, lines=lines)
    return int(lines[n].split()[1])

def gaussian_nfreq(s):
    """
    Return the number of vibrational degrees of freedom for
    a given log.
    """
    natom = gaussian_natom(s)
    if gaussian_islinear(s):
        nvdof = 3*natom - 5
    else:
        nvdof = 3*natom - 6
    return nvdof

  
def gaussian_basisset(lines):

    bas = 'Standard basis:\s*(\S*)'
    bas = re.findall(bas,lines)
    bas[-1] = bas[-1].replace('(d)','*')
    return bas[-1]

def gaussian_method(lines):

    method = 'Done:\s*E\((\w+)'
    method = re.findall(method,lines)
    #if method[-1].strip().upper() == 'CORR' or method[-1].strip().upper() == 'Z':
    if 'CCSD(T)' in lines:
        return 'CCSD(T)'
    elif 'CCSD' in lines:
        return 'CCSD'
    method = method[-1].lstrip('r').lstrip('u').lstrip('R').lstrip('U')
    if 'HF' in method:
        if 'MP2' in lines:
            return 'MP2'
        if 'MP4' in lines:
            return 'MP4'
    return method

def gaussian_energy(lines,method=''):

    if method == '':
        method = gaussian_method(lines)
    if 'CCSD' in method or 'MP' in method:
        method = method.replace('(','\(').replace(')','\)')
    #    energ  = method + '=([u,U,r,R]*[\w,\.,\s,-,D,\+]*)'
        energ  = method + '=([u,U,r,R]*[\w,\.,\s,-]*)'
        energ  = re.findall(energ,lines.replace('\n','').replace(' ',''))
    #    return (method,float(energ[-1].replace('D','E').replace('\n','').replace(' ','')))
        return (method,float(energ[-1].replace('\n','').replace(' ','')))
    else:
        lines = lines.strip().replace('\n','').replace(' ','')
        if 'anharm' in lines:
            energ = 'MP2=\s*([\d,\-,\.,D,\+]*)'
            energ = re.findall(energ,lines)
            if energ:
                return (method, float(energ[-1].replace('D','E')))
            else:
                energ = 'HF=\s*([\d,\-,\.,D,\+]*)'
                energ = re.findall(energ,lines)
                return (method, float(energ[-1].replace('D','E')))
       # energ = '(\S+)\s*A\.U\.'
        energ = 'E\([u,U,r,R]*' + method + '\)\s*=\s*([\d,\-,\.,D,\+]*)'
        energ = re.findall(energ,lines)
        return (method, float(energ[-1].replace('D','E')))
    return 

def gaussian_opt_zmat_params(lines):
    
    params = ''
    if not 'Optimized Parameters' in lines:
        return None
    varis  = lines.split('Optimized Parameters')[1].split('---------------------------------------')
    if 'Definition' in varis[0]:
        varis = varis[2].split('\n')
        for var in varis[1:-1]:
            var = var.split()
            params += ' '+  var[1] + '\t' + var[3] + '\n'
    else:
        varis = varis[1].split('\n')
        for var in varis[1:-1]:
            var = var.split()
            params += ' '+  var[1] + '\t' + var[2] + '\n'
    return params

def gaussian_zmat(lines):

    lines  = lines.split('Z-matrix:\n')
    zmat   = lines[1].split('       Variables:')[0]
    zmat  += 'Variables:\n'
    zmat   = zmat.replace('Charge = ','')
    zmat   = zmat.replace('Multiplicity =','')
    optzmat = gaussian_opt_zmat_params(lines[1])
    if optzmat == None:
        return None
    return zmat + optzmat

def gaussian_freqs(lines):
    """
    Return harmonic frequencies.
    """
    nfreq = gaussian_nfreq(lines)

    freqs = []
    lines = lines.splitlines()
    key = 'Fundamental Bands (DE w.r.t. Ground State)'
    iline = io.get_line_number(key,lines=lines)
    if iline > 0:
        for i in range(nfreq):
            iline += 1
            line = lines[iline]
            cols = line.split()
            freqs.append(cols[-5])
    else:
        lines = '\n'.join(lines)
        kw = 'Frequencies --  (.+)'
        freqlines = re.findall(kw, lines)
        freqs = []
        k = 0
        if len(freqlines) > 0:
            freqs = ['']*nfreq
            for i in range(len(freqlines)):
                tokens = freqlines[i].split()
                for j in range(len(tokens)):
                    freqs[k] = tokens[j]
                    k += 1
                if k == nfreq:
                    break
    return freqs

def gaussian_hessian(lines):
    startkey = 'Force constants in Cartesian coordinates:'
    endkey   = 'Force constants in internal coordinates:'
    lines= lines.split('Harmonic vibro-rotational analysis')[-1]
    lines = lines.splitlines()
    sline = io.get_line_number(startkey,lines=lines)
    eline = io.get_line_number(endkey,lines=lines)
    if sline < 0:
        return ''
    hess   = '\n'.join(lines[sline+1:eline]).replace('D','E')
    return hess

def gaussian_zpve(lines):

    zpve = 'Zero\-point\s*correction=\s*([\d,\.,\-]*)'
    zpve = re.findall(zpve, lines)
    if len(zpve) > 0:
        return float(zpve[-1])
    return 0.0 

def gaussian_anzpve(lines):
    """
    Returns anharmonic zpve in au.
    """
    from unittools import rcm2au,kj2au
    zpve = 'ZPE\(anh\)=\s*([\d,\w,\+,\.,\-]*)'
    zpve = re.findall(zpve, lines)
    if len(zpve) > 0:
        return float(zpve[-1].replace('D','E')) * kj2au
    zpve = 'Total Anharm\s*:\s*cm\-1\s*=\s*([\d,\.,\-]*)'
    zpve = re.findall(zpve, lines)
    if len(zpve) > 0:
        return float(zpve[-1].replace('D','E')) * rcm2au
    return
 
def gaussian_calc(lines):
    if 'Optimization complete' in lines:
        return 'geometry optimization'
    else:
        return ''

def gaussian_xyz_foresk(lines):
    atom = lines.split('Distance matrix')[-1].split('Symm')[0]
    if len(atom.split('\n')) > 8:
        atom = atom.split(' 6 ')[0] + ' 6 ' + atom.split(' 6 ')[1]
        atom = atom.split('\n')[2:-2]
    length = len(atom)
    atoms  = []
    for at in atom:
        atoms.extend(at.split()[1])
    xyz = 'Geometry ' + str(length) + ' Angstrom\n'
    if 'Eckart' in lines:
        lines = lines.split('Gaussian Orientation')[-1].split('Eckart')[0]
        lines = lines.split('\n')[5:-2]
        for i,line in enumerate(lines):
            line = line.split()
            xyz += atoms[i] + '  ' + line[2] + '  ' + line[3] + '  ' + line[4] + '\n'
    else:
        lines = lines.split('Coordinates (Angstroms)')[-1].split(' Distance matrix')[0]
        lines = lines.split('\n')[3:-2]
        for i,line in enumerate(lines):
            line = line.split()
            xyz +=  atoms[i] + '  ' + line[3] + '  ' + line[4] + '  ' + line[5] + '\n'
    return xyz 
    
def gaussian_geo(lines):
    atomnum = {'1':'H','6':'C','7':'N','8':'O'}
    xyz = ''
    try:
        if 'Eckart' in lines:
            lines = lines.split('Gaussian Orientation')[-1].split('Eckart')[0]
            lines = lines.split('\n')[5:-2]
            for i,line in enumerate(lines):
                line = line.split()
                xyz += atomnum[line[1]]+ '  ' + line[2] + '  ' + line[3] + '  ' + line[4] + '\n'
        else:
            lines = lines.split('Coordinates (Angstroms)')[-1].split(' Distance matrix')[0].split(' Rotation')[0].split('Symm')[0]
            lines = lines.split('\n')[3:-2]
            for i,line in enumerate(lines):
                line = line.split()
                xyz += ' ' + atomnum[line[1]] + '  ' + line[3] + '  ' + line[4] + '  ' + line[5] + '\n'
    except:
        logging.error('Cannot parse xyz')
    return xyz
   
def gaussian_xyz(lines):
    geo = gaussian_geo(lines) 
    if geo:
        n   = str(len(geo.splitlines()))
        xyz = n + '\n\n' +  geo
    else:
        xyz = ''
    return xyz

def gaussian_rotconstscent(lines):
    startkey = 'Effective Rotational Constants'
    lines = lines.splitlines()
    sline = io.get_line_number(startkey,lines=lines)
    if sline < 0:
        return ''
    rotlines   =  lines[sline+4:sline+7]
    constants = []
    for line in rotlines:
        constants.append(line.split()[1])
    return constants

def gaussian_rotconsts(lines):
    rot = 'Rotational constants\s*\(GHZ\):\s*([\s,\d,\.,\-]*)'     
    rot = re.findall(rot,lines)
    if len(rot) > 0: 
        rot = rot[-1].split()
    ndof  = gaussian_nfreq(lines)
    if ndof < 2:
        rot = rot[1:]
    if len(rot) > 0:
         if abs(float(rot[0])) < 0.000001:
             rot = rot[1:]
    return rot
 
def gaussian_rotdists (lines):
    startkey = 'Quartic Centrifugal Distortion Constants Tau Prime'
    endkey   = 'Asymmetric Top Reduction'
    lines = lines.splitlines()
    sline = io.get_line_number(startkey,lines=lines)
    if sline < 0:
        return ''
    lines  = lines[sline+3:sline+9]
    distlines = []
    for line in lines:
        splitline = line.split()
        if splitline[0] == 'TauP': 
           distlines.append('\t'.join(splitline[1:3]))
        else:
           break
    constants   = '\n'.join(distlines).replace('D','e')
    return constants

def gaussian_vibrot(lines):
    startkey = 'Vibro-Rot alpha Matrix (in cm^-1)'
    ndof  = gaussian_nfreq(lines)
    lines = lines.splitlines()
    sline = io.get_line_number(startkey,lines=lines)
    if sline < 0:
        return ''
    lines =  lines[sline+3:sline+3+ndof]
    for i in range(len(lines)):
       if ')' in lines[i]:
           lines[i] = lines[i].split(')')[1]
       if ndof < 2:
          lines[i] = '\t'.join(lines[i].split()[:-1])
    mat   = '\n'.join(lines).split('---------------')[0]
    return mat
    
##############################################
############      MOLPRO PARSER    ###########
##############################################
def molpro_energy(lines,method=''):
    if method == '':
        method = molpro_method(lines)
    method = method + '[a,b]?'
    method = method.replace('(','\(').replace(')','\)')
    if 'OPTG' in lines:
        energ = 'E\([U,u,R,r]*' + method +'\) \/ Hartree\s*[\d,\-,\.]*\s*([\d,\-,\.]*)'
        energ  = re.findall(energ,lines)
        if len(energ) != 0:
            return (method.rstrip('[a,b]?'), float(energ[-1].replace('\n','').replace(' ','')))

    if 'CCSD' in method:
        energ  = method + ' total energy\s*([\d,\-,\.]+)'
        energ  = re.findall(energ,lines)
        if len(energ) == 0:
            energ  = '!\w*\-\s*[\U,\R]' + method + '\s*energy\s*([\d,\-,\.]+)'
            energ  = re.findall(energ,lines)
            if len(energ) > 0:
                return (method.rstrip('[a,b]?'), float(energ[-1].replace('\n','').replace(' ','')))
        else:
            return (method.rstrip('[a,b]?'), float(energ[-1].replace('\n','').replace(' ','')))

    elif 'HF' in method:
        energ = method + ' STATE\s*\d\.\d\s*Energy\s*([\w,\-,\.]+)'
        energ = re.findall(energ,lines)
        if len(energ) > 0:
            return (method.rstrip('[a,b]?'), float(energ[-1].replace('\n','').replace(' ','')))
    
    elif 'MP' in method:
        energ = ' ' + method + ' total energy:\s*([\w,\-,\.]+)'
        energ = re.findall(energ,lines)
        if  len(energ) > 0:
            return (method.rstrip('[a,b]?'), float(energ[-1].replace('\n','').replace(' ','')))

    energ  = 'SETTING ENERGY\s*=\s*([\w,\.,-]+)'
    energ  = re.findall(energ,lines)
    if len(energ) == 0:
        energ  = 'SETTING CBSEN\s*=\s*([\w,\.,-]+)'
        energ  = re.findall(energ,lines)
        if len(energ) > 0:
            return (method.rstrip('[a,b]?'),float(energ[-1].replace('\n','').replace(' ','')))
    else:
        return (method.rstrip('[a,b]?'),float(energ[-1].replace('\n','').replace(' ','')))
    if len(energ) == 0:
        print 'energy not found'
    return 0 
   
def  molpro_freqs(lines):

    freqs = 'Wavenumbers \[cm-1\]   (.+)'
    freqlines = re.findall(freqs, lines)
    freqs = []
    if freqlines == []:
        return []
    for line in freqlines:
        if line.split()[0].strip() != '0.00': 
            freqs.extend(line.split())
    return freqs

def molpro_hessian(lines):
    startkey = 'Force Constants (Second Derivatives of the Energy)'
    endkey   = 'Atomic Masses'
    lines = lines.splitlines()
    sline = io.get_line_number(startkey,lines=lines)
    eline = io.get_line_number(endkey,lines=lines)
    if sline < 0:
        return ''
    hess = ''
    for line in lines[sline+1:eline-2]:
       hessline = ''
       for val in line.split():
           if 'G'  in val:
               if 'GX' in val:
                   add = 1
                   val = val.replace('GX','')
               elif 'GY' in val:
                   add = 2
                   val = val.replace('GY','')
               else:
                   add = 3
                   val = val.replace('GZ','')
               val =  str( (int(val) - 1 ) * 3 + add)
           hessline += '\t' +  val
       hess +=  hessline + '\n'
    return hess


def molpro_zpve(lines):
    zpve = 'Zero point energy:\s*([\d,\-,\.]*)'
    zpve = re.findall(zpve, lines)
    if len(zpve) > 0:
        return float(zpve[-1])
    return 0.0

def molpro_method(lines):
    method  = '1PROGRAM\s*\*\s*(\S*)'
    method  = re.findall(method,lines)
    if 'CCSD(T)' in lines:
        method =  ['CCSD(T)']
    elif 'MP2' in lines:
        method = ['MP2']
    elif 'CCSD' in lines:
        method = ['CCSD']
    if method[-1] == 'DFT':
        method = 'dft=\[([\d,\w]*)\]'
        method  = re.findall(method,lines)
    method = method[-1].lstrip('r').lstrip('u').lstrip('R').lstrip('U')
    if 'F12' in lines:
        method += '-F12'
    return method

def molpro_calc(lines):
    if 'optg' in lines:
       return 'geometry optimization'
    else:
       return ''

def molpro_basisset(lines):

    basis  = 'basis=(\S*)' 
    basis  = re.findall(basis,lines)
    basis[-1] = basis[-1].replace('(d)','*')
    return basis[-1]

def molpro_rotconsts(lines):

    rot  = 'Rotational constants:\s*([\s,\d,\.,\-]*)' 
    rot = re.findall(rot,lines)
    rot = rot[-1].split()
    return rot
 
def molpro_zmat(lines):
    geolines =  lines.split('geometry={')[1].split('}')[0].split('\n')[1:-1]
    zmat = 'geometry={angstrom \n' + '\n'.join(geolines) + '\n}\n'
    optzmat = False
    if 'OPTG' in lines:
        optzmat = True
        lines = lines.split('END OF GEOMETRY OPTIMIZATION')[0].split('Variable')[-1].split('\n')
        lines = lines[3:-3]
        for line in lines:
            zmat += line.split()[0].lower() + ' =  ' + line.split()[4] + '\n'
    if optzmat:
        return zmat
    return None

def molpro_geo(lines):
    lines =  lines.split('Current geometry (xyz format')
    if len(lines) > 1:
        lines = lines[-1].split('************')[0]
        lines =  lines.split('\n')[4:]
        return '\n'.join(lines)
    else:
        lines =  lines[0].split('Dump information in style XYZ')
        if len(lines) > 1:
            lines = lines[-1].split('************')[0]
            lines =  lines.split('\n')[4:]
            return '\n'.join(lines)
    return

def molpro_xyz(lines):
    lines =  lines.split('Current geometry (xyz format')
    if len(lines) > 1:
        lines = lines[-1].split('************')[0]
        lines =  lines.split('\n')[2:]
        return '\n'.join(lines)
    else:
        lines =  lines[0].split('Dump information in style XYZ')
        if len(lines) > 1:
            lines = lines[-1].split('************')[0]
            lines =  lines.split('\n')[2:]
            return '\n'.join(lines)
    return


##############################################
############     QCHEM PARSER    ###########
##############################################
def qchem_geo(lines):
    xyz = ''
    try:
        if 'OPTIMIZATION CONVERGED' in lines:
            lines = lines.split('OPTIMIZATION CONVERGED')[1]
            lines = lines.split('\n\n')[1]
            lines = lines.splitlines(True)[2:]
            for line in lines:
                line = line.split()
                xyz += ' ' + line[1] + '  ' + line[2] + '  ' + line[3] + '  ' + line[4] +'\n' 
    except:
        logging.error('Cannot parse xyz')
    return xyz

def qchem_xyz(lines):
    geo = qchem_geo(lines) 
    if geo:
        n   = str(len(geo.splitlines()))
        xyz = n + '\n\n' +  geo
    else:
        xyz = ''
    return xyz
   
def qchem_energy(lines):
    energy  = 'Final energy is\s*([\d,\.,-]*)' 
    energy  = re.findall(energy,lines)
    if len(energy) < 1:
        energy  = 'energy in the final basis set =\s*([\d,\.,-]*)'
        energy  = re.findall(energy,lines)
    return float(energy[-1])

def qchem_method(lines):
    method  = 'method\s*(\S*)' 
    method  = re.findall(method,lines)
    return  method[-1]

def qchem_basisset(lines):
    method  = 'Requested basis set is\s*(\S*)' 
    method  = re.findall(method,lines)
    return  method[-1].lower()

def qchem_freqs(lines):
    """
    Return harmonic frequencies.
    """
    kw = 'Frequency: (.+)'
    freqlines = re.findall(kw, lines)
    freqs = []
    for line in freqlines:
        freqs.extend(line.split())
    return freqs

def qchem_zpve(lines):

    zpve = 'Zero point vibrational energy:\s*([\d,\.,\-]*)'
    zpve = re.findall(zpve, lines)
    if len(zpve) > 0:
        return float(zpve[-1]) / ut.au2kcal
    return 0.0 

def qchem_calc(lines):
    if 'OPTIMIZATION COMPLETE' in lines:
        return 'geometry optimization'
    else:
        return ''
##############################################
############     EStokTP PARSER    ###########
##############################################

def EStokTP_freqs(lines):
    """
    Pulls the frequencies out from EStokTP me output file 
    INPUT:
    lines -    lines from EStokTP output file (reac1_fr.me or reac1_unpfr.me)
    OUTPUT:
    freqs    - frequencies obtained from output file
    """
    import numpy as np
 
    lines  = lines.strip('\n')
    lines  = lines.split('[1/cm]')[1].split('Zero')[0] 
    lines  = lines.split()
    nfreqs = lines[0]
    freqs  = lines[1:]
    freqs  = np.array(map(float, freqs))
    freqs  = np.sort(freqs)[::-1]
    return freqs.tolist()

###########################
#####  FOR GENERAL ########
############################

def method(lines):
    prog = get_prog(lines)
    if prog == 'gaussian':
        return gaussian_method(lines)
    if prog == 'molpro':
        return molpro_method(lines)
    print 'program not recognized as gaussian or molpro'
    return

def basisset(lines):
    prog = get_prog(lines)
    if prog == 'gaussian':
        return gaussian_basisset(lines)
    if prog == 'molpro':
        return molpro_basisset(lines)
    print 'program not recognized as gaussian or molpro'
    return

def energy(lines):
    prog = get_prog(lines)
    if prog == 'gaussian':
        return gaussian_energy(lines)
    if prog == 'molpro':
        return molpro_energy(lines)
    print 'program not recognized as gaussian or molpro'
    return

def zmat(lines):
    prog = get_prog(lines)
    if prog == 'gaussian':
        return gaussian_zmat(lines)
    if prog == 'molpro':
        return molpro_zmat(lines)
    print 'program not recognized as gaussian or molpro'
    return

def freqs(lines):
    prog = get_prog(lines)
    if prog == 'gaussian':
        return gaussian_freqs(lines)
    if prog == 'molpro':
        return molpro_freqs(lines)
    print 'program not recognized as gaussian or molpro'
    return

def zpve(lines):
    prog = get_prog(lines)
    zpve = 0.
    if prog == 'gaussian':
        zpve = gaussian_zpve(lines)
    elif prog == 'molpro':
        zpve = molpro_zpve(lines)
    else:
        'program not recognized as gaussian or molpro'
    return 0.0

def anzpve(lines):
    prog = get_prog(lines)
    if prog == 'gaussian':
        return gaussian_anzpve(lines)
    if prog == 'molpro':
       return #molpro_anzpve(lines)
    print 'program not recognized as gaussian or molpro'
    return 

def xyz(lines):
    prog = get_prog(lines)
    if prog == 'gaussian':
        return gaussian_xyz(lines)
    if prog == 'molpro':
        return molpro_xyz(lines)
    print 'program not recognized as gaussian or molpro'
    return

def geo(lines):
    prog = get_prog(lines)
    if prog == 'gaussian':
        return gaussian_geo(lines)
    if prog == 'molpro':
        return molpro_geo(lines)
    print 'program not recognized as gaussian or molpro'
    return

def rotconsts(lines):
    prog = get_prog(lines)
    if prog == 'gaussian':
        return gaussian_rotconsts(lines)
    if prog == 'molpro':
        return molpro_rotconsts(lines)
    print 'program not recognized as gaussian or molpro'
    return

def get_298(lines):
    deltaH298 = ' h298 final\s*([\d,\-,\.]*)'
    lines = lines.splitlines()
    tmp = ''
    for line in lines:
        if 'h298 final' in line:
            tmp = float(line.split()[-1])
    return tmp
    
 
