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
import tctools as tc
try:
    runserial = False
    from scoop import futures
    from scoop import utils
except:
    runserial = True
    print "scoop not available, no concurency"

def get_args():
    """
    Returns args object that contains command line options.
    """
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=
    """
    April 18, 2017
    Murat Keceli
    
    Performs quantum chemistry calculations to calculate thermochemical parameters.
    Writes NASA polynomials in different formats.
    Uses different codes for these purposes
    """)
    parser.add_argument('-i','--input',type=argparse.FileType('r'),nargs=1,
                        default='qc_list.txt',
                        help='List of inchi or smiles for species to be calculated')
    parser.add_argument('-m','--messpf',type=str,nargs=1,
                        default='messpf',
                        help='Path for mess partition function executable')
    parser.add_argument('-t','--thermp',type=str,nargs=1,
                        default='thermp',
                        help='Path for thermp executable')
    parser.add_argument('-p','--pac99',type=str,nargs=1,
                        default='pac99',
                        help='Path for pac99 executable')
    parser.add_argument('-r','--runserial',action='store_true',
                        help='Run serial')
    return parser.parse_args()

def get_chemkin_polynomial(mol,method,zpe,xyz,freqs,deltaH):
    """
    A driver to perform all operations to write NASA polynomial in
    chemkin format. Assumes quantum chemistry calculation is performed.
    """
    inputfile = 'pf.inp'
    name = mol.formula
    tag = method
    inp = tc.get_pf_input(mol, method, zpe, xyz,freqs)
    
 #   print 'Running mess partition function'
    run_pf()
 #   print 'Generate thermp input'
    write_thermp_input(mol.formula, deltaH)
 #   print 'Running thermp'
    run_thermp()
#    print 'Running pac99'
    run_pac99(name)
#    print 'Converting to chemkin format'
    chemkinfile= name 
    write_chemkin_file(deltaH, tag, name, name+'.ckin')
    return 

def run_pm3(s):
    import os
    import qctools as qc
    import obtools as ob
    method = 'pm3' 
    mol = ob.get_mol(s)
    natom = ob.get_natom(mol)
    formula = mol.formula
    hformula = qc.get_heavyatomformula(formula)
    hlist    = qc.get_heavyatomlist(hformula)
    key = ob.get_uniquekey(mol)
    dir = 'database/' + hlist + '/' + hformula + '/' + formula + '/' + key + '/' + method
    groupsfile='new.groups'
    io.mkdir(dir)
    if io.check_file(groupsfile):
        if io.check_file(dir,1):
            io.cp(groupsfile,dir)
            io.cd(dir)
        else:
            print 'dir error', dir
    else:
        print 'groups file error'
        return s    
    lines = qc.run_mopac(mol=mol)
    xyz = qc.get_xyz_mopac(lines)
    freqs = qc.get_freq_mopac(lines)
    zpe = qc.get_zpe_mopac(lines)
    deltaH = qc.get_deltaH_mopac(lines)
    return mol,method,zpe,xyz,freqs,deltaH


def run_all(s):
    import os
    import qctools as qc
    import obtools as ob
    method = 'pm3' 
    mol = ob.get_mol(s)
    natom = ob.get_natom(mol)
    formula = mol.formula
    hformula = qc.get_heavyatomformula(formula)
    hlist    = qc.get_heavyatomlist(hformula)
    key = ob.get_uniquekey(mol)
    dir = 'database/' + hlist + '/' + hformula + '/' + formula + '/' + key + '/' + method
    groupsfile='new.groups'
    io.mkdir(dir)
    if io.check_file(groupsfile):
        io.cp(groupsfile,dir)
        io.cd(dir)
    else:
        print 'new.groups file required in working directory'
        return -1
    if not io.check_file(groupsfile,1):
        print 'new.groups file required in target directory'
        return -1    
    lines = qc.run_mopac(mol=mol)
    xyz = qc.get_xyz_mopac(lines)
    freqs = qc.get_freq_mopac(lines)
    zpe = qc.get_zpe_mopac(lines)
    deltaH = qc.get_deltaH_mopac(lines)
    get_chemkin_polynomial(mol,method,zpe,xyz,freqs,deltaH)
    return 0    
    
if __name__ == "__main__":
    args   = get_args()
    mylist = io.read_list('qc_list.txt')
    if runserial:
        returnValues = list(map(run_all,mylist))
    else:
        returnValues = list(futures.map(run_all,mylist))

    