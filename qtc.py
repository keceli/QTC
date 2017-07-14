#!/usr/bin/env python
import argparse
import multiprocessing
import iotools as io
import obtools as ob
import qctools as qc
import tctools as tc
import dbtools as db
import heatform as hf
from pprint import pprint

__updated__ = "2017-07-13"
__authors__ = 'Murat Keceli, Sarah Elliott'
__logo__ = """
***************************************

     <===>   <=============>   <=>    
  <=>     <=>      <=>      <=>   <=> 
<=>         <=>    <=>     <=>        
<=>         <=>    <=>     <=>        
<=>         <=>    <=>     <=>        
  <=>     <=>      <=>      <=>   <=> 
     <===><>       <=>         <=>   
         <<>>                           
                                      
***************************************
For computation of accurate thermochemistry
By ECP-PACC team                                  
"""
def get_args():
    """
    Returns args object that contains command line options.
    """
    import argparse
    parser = argparse.ArgumentParser(#formatter_class=argparse.RawDescriptionHelpFormatter,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=
    """
    April 18, 2017

    Performs quantum chemistry calculations to calculate thermochemical parameters.
    Writes NASA polynomials in different formats.
    Provides a unified interface for different codes.
    """)
    parser.add_argument('-i', '--input', type=str,
                        default='qclist.txt',
                        help='List of inchi or smiles for species to be calculated')
    parser.add_argument('-f', '--first', type=int,
                        default=0,
                        help='Beginning index of the species list')
    parser.add_argument('-l', '--last', type=int,
                        help='Ending index of the species list')
    parser.add_argument('-n', '--nproc', type=int,
                        default=1,
                        help='Number of processors, default is one processor')
    parser.add_argument('-c', '--qccalculation', type=str,
                        default='',
                        help='Type of quantum chemistry calculation: "optimization", "frequency", "energy", "anharmonic"')
    parser.add_argument('-k', '--qckeyword', type=str,
                        default='',
                        help='Keyword string that can define package, method, basis, task i.e.: "opt/ccsd/cc-pvdz/gaussian, nwchem/ccsdt/cc-pvdz/energy,')
    parser.add_argument('-m', '--qcmethod', type=str,
                        default='',
                        help='Quantum chemistry method to be used')
    parser.add_argument('-b', '--qcbasis', type=str,
                        default='',
                        help='Basis-set in quantum chemistry calculations')
    parser.add_argument('-a', '--qctask', type=str,
                        default='',
                        help='Task for quantum chemistry calculations (opt,freq,energy)')
    parser.add_argument('-z', '--qcpackage', type=str,
                        default='',
                        help='Quantum chemistry package ("gausian","mopac","nwchem","qcscript") to be used')
    parser.add_argument('-p', '--qcnproc', type=int,
                        default=1,
                        help='Number of processors for quantum chemistry calculations')
    parser.add_argument('-t', '--qctemplate', type=str,
                        default='',
                        help='Template for quantum chemistry input file')
    parser.add_argument('-d', '--qcdirectory', type=str,
                        default='',
                        help='Path for the directory for running qc jobs')
    parser.add_argument('-e', '--qcexe', type=str,
                        default='',
                        help='Path for the executable for qc calculations')
    parser.add_argument('-x', '--xyzpath', type=str,
                        default='',
                        help='Path for the directory for xyz file')
    parser.add_argument('-o', '--qcoutput', type=str,
                        default='',
                        help='Path for the qc output file')
    parser.add_argument('-B', '--hfbasis', type=str,
                        default='auto',
                        help='Heat of formation basis molecules')
    parser.add_argument('-D', '--database', type=str,
                        default= io.pwd(),
                        help='Heat of formation basis molecules')
    parser.add_argument('-Q', '--runqc', action='store_true',
                        help='Run quantum chemistry calculation')
    parser.add_argument('-P', '--parseqc', action='store_true',
                        help='Parses quantum chemistry output')
    parser.add_argument('-S', '--subqc', action='store_true',
                        help='Submit quantum chemistry calculations')
    parser.add_argument('-T', '--runthermo', action='store_true',
                        help='Run thermochemistry calculations')
    parser.add_argument('-W', '--writefiles', action='store_true',
                        help='Write .xyz, .ene files')
    parser.add_argument('-w', '--storefiles', action='store_true',
                        help='Store .xyz, .ene files')
    parser.add_argument('-I', '--runinteractive', action='store_true',
                        help='Interactive mode for QTC')
    parser.add_argument('-O', '--overwrite', action='store_true',
                        help='Overwrite existing calculations. Be careful, data will be lost.')
    parser.add_argument('-A', '--anharmonic', action='store_true',
                        help='Anharmonic corrections')
    parser.add_argument('--mopac', type=str,
                        default='mopac',
                        help='Path for mopac executable')
    parser.add_argument('--nwchem', type=str,
                        default='nwchem',
                        help='Path for nwchem executable')
    parser.add_argument('--molpro', type=str,
                        default='molpro',
                        help='Path for molpro executable')
    parser.add_argument('--gaussian', type=str,
                        default='g09',
                        help='Path for gaussian executable')
    parser.add_argument('--torsscan', type=str,
                        default='/home/elliott/Packages/TorsScan/torsional_scan.py',
                        help='Path for torsscan executable')
    parser.add_argument('--messpf', type=str,
                        default='messpf',
                        help='Path for MESS partition function executable')
    parser.add_argument('--thermp', type=str,
                        default='thermp',
                        help='Path for thermp executable')
    parser.add_argument('--pac99', type=str,
                        default='pac99',
                        help='Path for pac99 executable')
    parser.add_argument('--qcscript', type=str,
#                        default='/lcrc/project/PACC/test-awj/builddb/bin/qcscript.pl',
                        default='',
                        help='Path for qcscript perl script')
    return parser.parse_args()



def run(s):
    """
    A driver function to run quantum chemistry and thermochemistry calculations.
    It updates parameters dictionary.
    """
    print("***************************************\n")
    global parameters
    runqc = parameters['runqc']
    parseqc = parameters['parseqc']
    package = parameters['qcpackage'].lower()
    runthermo = parameters['runthermo']
    qcnproc = parameters['qcnproc']
    if package in ['nwchem', 'molpro', 'mopac', 'gaussian', 'torsscan' ]:
        print('Quantum chemistry package: {0}'.format(package))
        if qcnproc > 1:
            if package.startswith('nwc'):
                parameters['qcexe'] = 'mpirun -n {0} nwchem'.format(qcnproc)
            elif package.startswith('mol'):
                parameters['qcexe'] = '{0} -n {1}'.format(parameters['molpro'],qcnproc)
            else:
                parameters['qcexe'] = parameters[package]
        else:
            parameters['qcexe'] = parameters[package]
    if not package:
        if parameters['qctemplate']:
            parameters['qcpackage'] = qc.get_package(parameters['qctemplate'])
            print('Quantum chemistry package (based on template): {0}'.format(package))
    elif not package.startswith('extrap'):
        templatename = package + '_template' + '.txt'
        parameters['qctemplate'] = io.join_path(*[parameters['qtcdirectory'],'templates',templatename])
    if parameters['writefiles']:
        parameters['parseqc'] = True
    msg = "Smiles: {0}\n".format(s)
    mult = ob.get_mult(s)
    mol = ob.get_mol(s)
    smilesname = ob.get_smiles_filename(s)
    smilesdir =  ob.get_smiles_path(mol, mult, method='', basis='')
    parameters['smilesdir'] = smilesdir
    workdirectory  = io.join_path(*[smilesdir,parameters['qcdirectory']])
    parameters['qctemplate'] = io.get_path(parameters['qctemplate'])
    qcpackage = parameters['qcpackage']
    qcscript = io.get_path(parameters['qcscript'])
    qcoutput = smilesname + '_' + qcpackage + '.out'
    xyzpath = parameters['xyzpath']
    xyzfile = qc.find_xyzfile(xyzpath, smilesdir)
    if xyzfile:
        msg += "Using xyz file in '{0}'\n".format(xyzfile)
        xyz = io.read_file(xyzfile)
        coords = ob.get_coordinates_array(xyz)
        mol = ob.set_xyz(mol, coords)
    print(msg)
    msg = ''
    cwd = io.pwd()
    io.mkdir(workdirectory)
    if io.check_dir(workdirectory, 1):
        io.cd(workdirectory)
        msg += "cd '{0}'\n".format(workdirectory)
    else:
        msg += ('I/O error, {0} directory not found.\n'.format(workdirectory))
        return -1
    print(msg)
    msg = ''
    available_packages=['nwchem', 'molpro', 'mopac', 'gaussian', 'extrapolation', 'torsscan']          
    if runqc:
        if qcpackage in available_packages:
            print('Running qcpackage: {0}'.format(qcpackage))
            msg += qc.run(mol, parameters, mult)
        elif qcpackage == 'qcscript':
            msg += "Running qcscript...\n"
            geofile = smilesname + '.geo'
            geo = ob.get_geo(mol)
            io.write_file(geo, geofile)
            if io.check_file(geofile, 1):
                msg += qc.run_qcscript(qcscript, parameters['qctemplate'], geofile, mult)
        else:
            msg = '{0} package not implemented\n'.format(qcpackage)
            msg += 'Available packages are {0}'.format(available_packages)
            exit(msg)   
        print(msg)
        msg = ''
    if parseqc:
        if io.check_file(qcoutput, timeout=1,verbose=False):
            out = io.read_file(qcoutput,aslines=False)
            d = qc.parse_output(out,smilesname, parameters['writefiles'], parameters['storefiles'], parameters['optlevel'])
            pprint(d)
                                   
    if runthermo:
        groupstext = tc.get_new_groups()
        io.write_file(groupstext, 'new.groups')
        msg += "Parsing qc logfile '{0}'\n".format(io.get_path(qcoutput))
        if parameters['qcpackage'] == 'torsscan':
            msg += io.read_file(io.get_path(qcoutput))
        else: 
            if parameters['qckeyword'].split(',')[parameters['calcindex']].startswith('extrap'):
                d = {}
            else:
                d = qc.parse_output(out,smilesname, parameters['writefiles'], parameters['storefiles'])
                xyz = next(db.gen_dict_extract('xyz',d))
                freqs = next(db.gen_dict_extract('harmonic frequencies',d))
                energy = float(next(db.gen_dict_extract('energy',d)))
                xmat   = next(db.gen_dict_extract('xmat',d))
            package   = parameters['qcpackage']
            optlevel  = parameters['optlevel']
            enlevel   = optlevel
            freqlevel = None
            if 'enlevel' in parameters:
                enlevel   = parameters['enlevel']
            if 'freqlevel' in parameters:
                freqlevel   = parameters['freqlevel']
            hfbasis = parameters['hfbasis']
            #hf0, hfset = hf.main(s,E=energy,optlevel=optlevel,enlevel=enlevel,
             #                   freqlevel=freqlevel, basis=hfbasis,anharm=parameters['anharmonic'])
            hf0, hfset = hf.main_keyword(s,parameters['qckeyword'], parameters['calcindex'],basis=hfbasis,anharm=parameters['anharmonic'],dbdir=parameters['database'])
            hftxt  = 'Energy (kcal)\tBasis\n----------------------------------'
            hftxt += '\n' + str(hf0) + '\t' + '  '.join(hfset) 
            io.write_file(hftxt,s + '.hf0k')
            if len(freqs) > 0:
                msg += tc.write_chemkin_polynomial(mol, xyz, freqs, hf0, parameters,xmat=xmat)
                if io.check_file('thermp.out'):
                    import patools as pa
                    lines = io.read_file('thermp.out')
                    msg += 'delHf(298) = {} kcal'.format(pa.get_298(lines))
    io.cd(cwd)
    print(msg)
    return


def main():
    from socket import gethostname
    from timeit import default_timer as timer
    import os
    global parameters
    start  = timer()
    args = get_args()
    parameters = vars(args)
    if parameters['qckeyword']:
        ncalc = len(parameters['qckeyword'].split(','))
    else:
        ncalc = 1
    parameters['qtcdirectory'] = os.path.dirname(os.path.realpath(__file__))
    parameters['number_of_calculations'] = ncalc
    parameters['optlevel'] = 'sp' #TODO
    beginindex = args.first
    endindex = args.last
    inp = args.input
    nproc = args.nproc
    if io.check_file(inp):
        mylist = io.read_list(inp)
    else:
        mylist = inp.split(',')
    if endindex:
        mylist = mylist[beginindex:endindex]
    else:
        mylist = mylist[beginindex:]        
    print(__logo__)
    print("QTC: Date and time           = {0}".format(io.get_date()))
    print("QTC: Last update             = {0}".format(__updated__))
    print("QTC: Number of processes     = {0}".format(nproc))
    print("QTC: Hostname                = {0}".format(gethostname()))
    print('QTC: Given arguments         =')
    for param in parameters:
        print('                             --{0:20s}\t{1}'.format(param, getattr(args, param)))
    print("QTC: Number of species       = {0}".format(len(mylist)))
    init = timer()
    print("QTC: Initialization time (s) = {0:.2f}".format(init-start))
    for i in range(ncalc):
        parameters['calcindex'] = i
        if parameters['qckeyword']:
            qc.parse_qckeyword(parameters, calcindex=i)
        if nproc > 1:
            pool = multiprocessing.Pool(nproc)
            pool.map(run, mylist)
        else:
            map(run, mylist)
    end = timer()
    print("QTC: Calculations time (s)   = {0:.2f}".format(end - init))
    print("QTC: Total time (s)          = {0:.2f}".format(end-start))
    print("QTC: Date and time           = {0}".format(io.get_date()))

if __name__ == "__main__":
    main()

