#!usr/bin/python

import re

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
       return 'g09'
   elif 'molpro' in lines:
       return 'molpro'
   return

################################################
############     Gaussian PARSER    ############
################################################

def gaussian_basisset(lines):

    bas = 'Standard basis:\s*(\S*)'
    bas = re.findall(bas,lines)

    return bas[-1]

def gaussian_method(lines):

    method = 'E\((\w+)'
    method = re.findall(method,lines)
    if method[-1] == 'CORR':
       if 'CCSD(T)' in lines:
           return 'CCSD(T)'
       elif 'CCSD' in lines:
           return 'CCSD'
    return method[-1]

def gaussian_energy(lines,method=''):

    if method == '':
        method = gaussian_method(lines)
    
    if 'CCSD' in method:
        method = method.replace('(','\(').replace(')','\)')
        energ  = method + '=([\w,\.,\s,-]*)'
        energ  = re.findall(energ,lines)
        return float(energ[-1].replace('\n','').replace(' ',''))

    else:
        energ = '(\S+)\s*A\.U\.'
        energ = re.findall(energ,lines)
        return (method, float(energ[-1]))

def gaussian_opt_zmat_params(lines):
    
    params = ''
    varis  = lines.split('Optimized Parameters')[1].split('--------------------------------------')
    if 'Definition' in varis[0]:
        varis = varis[1].split('\n')
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

    lines  = lines.split('Z-matrix:')
    zmat   = lines[1].split('       Variables:')[0]
    zmat  += 'Variables:\n'
    zmat   = zmat.replace('Charge = ','')
    zmat   = zmat.replace('Multiplicity =','')
    zmat  += gaussian_opt_zmat_params(lines[1])
    return zmat

def  gaussian_freqs(lines):

    freqs = 'Frequencies --  (.+)'
    freqlines = re.findall(freqs, lines)
    freqs = []
    for line in freqlines:
        freqs.extend(line.split())
    return freqs

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
    
def gaussian_xyz(lines):
    atom = lines.split('Distance matrix')[-1].split('Symm')[0]
    if len(atom.split('\n')) > 8:
        atom = atom.split(' 6 ')[0] + ' 6 ' + atom.split(' 6 ')[1]
        atom = atom.split('\n')[2:-2]
    else:
        atom = atom.split('\n')[2:-1]
    length = len(atom)
    atoms  = []
    for at in atom:
        atoms.extend(at.split()[1])
    xyz = ''
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
            xyz += ' ' + atoms[i] + '  ' + line[3] + '  ' + line[4] + '  ' + line[5] + '\n'
    return xyz 
    
     
##############################################
############      MOLPRO PARSER    ###########
##############################################
def molpro_energy(lines,method=''):
    
    if method == '':
        method = molpro_method(lines)
    method = method.replace('(','\(').replace(')','\)')
    if 'OPTG' in lines:
        energ = 'E\(' + method +'\) \/ Hartree\s*[\d,\-,\.]*\s*([\d,\-,\.]*)'
        energ  = re.findall(energ,lines)
        if len(energ) != 0:
            return float(energ[-1].replace('\n','').replace(' ',''))

    if 'CCSD' in method:
        energ  = method + ' total energy\s*([\d,\-,\.]+)'
        energ  = re.findall(energ,lines)
        if len(energ) == 0:
            energ  = '!\w*\-\s*[\U,\R]' + method + '\s*energy\s*([\d,\-,\.]+)'
            energ  = re.findall(energ,lines)
            if len(energ) == 0:
                print 'energy not found'
            else:
                return float(energ[-1].replace('\n','').replace(' ',''))
        else:
            return float(energ[-1].replace('\n','').replace(' ',''))

    elif 'HF' in method:
        energ = method + ' STATE\s*\d\.\d\s*Energy\s*([\w,\-,\.]+)'
        energ = re.findall(energ,lines)
        print energ
        if len(energ) == 0:
            print 'energy not found'
        else:
            return float(energ[-1].replace('\n','').replace(' ',''))
    
    elif 'MP' in method:
        energ = ' ' + method + ' total energy:\s*([\w,\-,\.]+)'
        energ = re.findall(energ,lines)
        if len(energ) == 0:
            print 'energy not found'
        else:
            return float(energ[-1].replace('\n','').replace(' ',''))

    energ  = 'SETTING ENERGY\s*=\s*([\w,\.,-]+)'
    energ  = re.findall(energ,lines)
    if len(energ) == 0:
        print 'energy not found'
    else:
        return (method,float(energ[-1].replace('\n','').replace(' ','')))

    return 0 
   
def  molpro_freqs(lines):

    freqs = 'Wavenumbers \[cm-1\]   (.+)'
    freqlines = re.findall(freqs, lines)
    freqs = []
    for line in freqlines:
        if line.split()[0].strip() != '0.00': 
            freqs.extend(line.split())
    return freqs

def molpro_method(lines):
    
    method  = '1PROGRAM \* (\S*)'
    method  = re.findall(method,lines)
    if 'CCSD(T)' in lines:
        return 'CCSD(T)'
    return method[-1]

def molpro_basisset(lines):

    basis = 'basis=(\S*)' 
    basis  = re.findall(basis,lines)
    return basis[-1]
      
def molpro_zmat(lines):
    geolines = lines.split('geometry={')[1].split('}')[0].split('\n')[1:-1]
    zmat = '\n'.join(geolines) + '\n'
    if 'OPTG' in lines:
        lines = lines.split('END OF GEOMETRY OPTIMIZATION')[0].split('Variable')[-1].split('\n')
        lines = lines[3:-3]
        for line in lines:
            zmat += line.split()[0] + '   ' + line.split()[4] + '\n'

    return zmat
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
    if prog == 'g09':
        return gaussian_method(lines)
    if prog == 'molpro':
        return molpro_method(lines)

def basisset(lines):
    prog = get_prog(lines)
    if prog == 'g09':
        return gaussian_basisset(lines)
    if prog == 'molpro':
        return molpro_basisset(lines)

def energy(lines):
    prog = get_prog(lines)
    if prog == 'g09':
        return gaussian_energy(lines)
    if prog == 'molpro':
        return molpro_energy(lines)

def zmat(lines):
    prog = get_prog(lines)
    if prog == 'g09':
        return gaussian_zmat(lines)
    if prog == 'molpro':
        return molpro_zmat(lines)

def freqs(lines):
    prog = get_prog(lines)
    if prog == 'g09':
        return gaussian_freqs(lines)
    if prog == 'molpro':
        return molpro_freqs(lines)

def xyz(lines):
    prog = get_prog(lines)
    if prog == 'g09':
        return gaussian_xyz(lines)
    if prog == 'molpro':
        return ''
