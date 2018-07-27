#!/usr/bin/env python
import argparse
import collections
import iotools as io
import obtools as ob
import qctools as qc
import tctools as tc
import unittools as ut
import dbtools as db
import heatform as hf
import pprint
import sys
import os
import logging
import numpy as np
from patools import energy
from timeit import default_timer as timer
__updated__ = "2018-03-03"
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
"""
TODO: 
When optimization fails, restart from the last geometry
Delete NWChem scratch files
Add torsX task/template
Add auto stamp for logfile
Check for negative energies in hindered potential
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
    September 11, 2017

    Performs quantum chemistry calculations to calculate thermochemical parameters.
    Writes NASA polynomials in different formats.
    Provides a unified interface for different codes.
    """)
    parser.add_argument('-i', '--input', type=str,
                        default='qclist.txt',
                        help='INPUT can be a text file containing a list of inchi or smiles strings, or it can be a single string containing inchi or smiles strings separated by commas')
    parser.add_argument('-j', '--jsoninput', type=str,
                        default='queue.json',
                        help='Input file in  json format, i.e. RMG queue.json files')
    parser.add_argument('-p', '--qcnproc', type=int,
                        default=io.get_ppn(),
                        help='Number of processors for quantum chemistry calculations. Default is total number of cores on a node.')
    parser.add_argument('-n', '--nproc', type=int,
                        default=1,
                        help='Number of processors for qtc calculations, to run different species in parallel')
    parser.add_argument('-m', '--machinefile', type=str,
                        default='',
                        help='Machinefile for mpirun')
    parser.add_argument('-k', '--qckeyword', type=str,
                        default='',
                        help='Keyword string that defines quantum chemistry calculations i.e.: "opt/ccsd/cc-pvdz/gaussian,energy/ccsd/cc-pvtz/nwchem,extrapolation/cbs/energy=0.3*E0+0.7*E1" Note that each calculation is separated by a comma (,) and calculations are defined by TASK/METHOD/BASIS/PACKAGE. TASK can be opt, freq, anharm,extrapolation.METHOD and BASIS are simply copied into quantum chemistry input file as defined in the templates folder. PACKAGE can be gaussian, molpro or nwchem')
    parser.add_argument('-t', '--qctemplate', type=str,
                        default='/home/elliott/Packages/QTC/templates',
                        help='Path for the templates directory. Templates have a specific format for filenames. See qtc/templates.')
    parser.add_argument('-l', '--loglevel', type=int,
                        default=-1,
                        help='Verbosity level of logging, -1: qtc decides,  0: errors, 1: 0 + warnings, 2: 1 + info, 3: 2 + debug')
    parser.add_argument('-f', '--logfile', type=str,
                        default= 'none',
                        help='Log file prefix, use none for logging to STDOUT, include DATE if you want a date stamp')
    parser.add_argument('-s', '--scratch', type=str,
                        default=io.get_env('TMPDIR',default='/tmp'),
                        help='Scratch directory. If not given checks TMPDIR env. variable, if not defined uses /tmp.')
    parser.add_argument('-x', '--xyzdir', type=str,
                        default= 'xyz',
                        help='Path for the xyz directory, where initial geometries will be read and written.')
    parser.add_argument('-g', '--logindex', type=str,
                        default= '',
                        help='Log file index')
    parser.add_argument('-d', '--database', type=str,
                        default= 'database',
                        help='Path for database directory')
    parser.add_argument('-o', '--qcoutput', type=str,
                        default='',
                        help='Path for the qc output file')
    parser.add_argument('-b', '--first', type=int,
                        default=1,
                        help='Beginning index of the species list')
    parser.add_argument('-e', '--last', type=int,
                        help='Ending index of the species list')
    parser.add_argument('-a', '--abcd', type=str,
                        default= '3,1,3,100',
                        help='a,b,c,d parameters for the number of MC sampling points based on min(a + b * c**nrotors, d)')
    parser.add_argument('-r', '--reference', type=str,
                        default='auto',
                        help='List of SMILES (seperated by commas) for reference thermo species')
    parser.add_argument('-z', '--task_seperator', type=str,
                        default=',',
                        help='The characther used for seperating different tasks in a qckeyword.')
    parser.add_argument('-B', '--bac', action='store_true',
                        help='Bond additivity correction using the bac.txt template file')
    parser.add_argument('-G', '--generate', action='store_true',
                        help='Generates a sorted list of species')
    parser.add_argument('-M', '--sortbymass', action='store_true',
                        help='Generates a sorted list of species')
    parser.add_argument('-Q', '--runqc', action='store_true',
                        help='Run quantum chemistry calculation')
    parser.add_argument('-P', '--parseqc', action='store_true',
                        help='Parses quantum chemistry output')
    parser.add_argument('-S', '--swift', action='store_true',
                        help='Required to run with swift')
    parser.add_argument('-T', '--runthermo', action='store_true',
                        help='Run thermochemistry calculations')
    parser.add_argument('-W', '--writefiles', action='store_true',
                        help='Write .xyz, .ene files')
    parser.add_argument('-w', '--storefiles', action='store_true',
                        help='Store .xyz, .ene files')
    parser.add_argument('-I', '--ignorerunningjobs', action='store_true',
                        help='Ignores RUNNING.tmp files and run the calculations')
    parser.add_argument('-R', '--recover', action='store_true',
                        help='Attempt to recover failed calculations')
    parser.add_argument('-O', '--overwrite', type=str,
                        default = 'none',
                        help='Overwrite existing calculations. Be careful, data will be lost.')
    parser.add_argument('-X', '--excel', action='store_true',
                        help='Generate excel file')
    parser.add_argument('-J', '--dumpjsonfile', action='store_true',
                        help='Writes a json file containing all the results')
    parser.add_argument('-A', '--anharmonic', action='store_true',
                        help='Anharmonic corrections')
    parser.add_argument('-D', '--debug', action='store_true',
                        help='Run in debug mode')
    parser.add_argument('--skippf', action='store_true',
                        help='use to skip pf input and output file generation when collecting thermo')
    parser.add_argument('--getBAC', action='store_true',
                        help='calculate the parameters for a BAC')
    parser.add_argument('-N', '--norot', action='store_true',
                        help='Turns off rotational pf input')
    parser.add_argument('--uncertainty', type=str, default='',
                        help='use to generate uncertainty analysis')
    parser.add_argument('--fix', type=int,
                        default=0,
                        help='If FIX > 0, interpolate negative energies in hindered potential input for mess')
    parser.add_argument('--mopac', type=str,
                        default='mopac',
                        help='Path for mopac executable')
    parser.add_argument('--nwchem', type=str,
                        default='nwchem',
                        help='Path for nwchem executable')
    parser.add_argument('--molpro', type=str,
                        default='molpro',
                        help='Path for molpro executable')
    parser.add_argument('--qchem', type=str,
                        default='qchem',
                        help='Path for molpro executable')
    parser.add_argument('--gaussian', type=str,
                        default='g09',
                        help='Path for gaussian executable')
    parser.add_argument('--torsscan', type=str,
                        #default='/home/elliott/Packages/TorsScan/torsional_scan.py',
                        default='/home/keceli/qtc/TorsScan/torsional_scan.py',
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
    parser.add_argument('--x2z', type=str,
                        default='x2z',
                        help='Path for x2z executable')
    parser.add_argument('--suppress_printing', action='store_true')
    parser.add_argument('--qcscript', type=str,
#                        default='/lcrc/project/PACC/test-awj/builddb/bin/qcscript.pl',
                        default='',
                        help='Path for qcscript perl script')
    return parser.parse_args()

def printp(string):
    if not parameters['suppress_printing']:
        print (string)
    return

def run(s):
    """
    A driver function to run quantum chemistry and thermochemistry calculations for a given molecule object identified by a SMILES or InChI string.
    It uses and modifies 'parameters', which is defined as a global variable.
    Not pure.
    """
    global parameters
    runqc = parameters['runqc']
    parseqc = parameters['parseqc']
    task    = parameters['qctask']
    runthermo = parameters['runthermo']
    qckeyword = parameters['qckeyword']
    calcindex = parameters['calcindex']
    ignore = parameters['ignorerunningjobs']
    overwrite=parameters['overwrite']
    over = False
    if overwrite == 'all':
        over = True
    elif len(overwrite.split('.')) > 1:
        overwrite = overwrite.split('.')
        if overwrite[0] == 'a':
            if overwrite[1] == 'a':
                over = True
            elif overwrite[1] == str(calcindex):
                over = True
        elif overwrite[0] == str(parameters['mol_index']):
            if overwrite[1] == 'a':
                over = True
            elif overwrite[1] == str(calcindex):
                over = True
        else:
            over = False
        overwrite = '.'.join(overwrite)
    elif overwrite == str(parameters['mol_index']):
        over = True
    parameters['overwrite'] = over
    scratch = parameters['scratch']
    xyz = parameters['xyz']
    mol = ob.get_mol(xyz)
    inchi = ob.get_inchi(mol)
    mult = parameters['mult']
    formula = parameters['formula']
    nrotor = parameters['nrotor']
    natom = parameters['natom']
    qlabel = qc.get_qlabel(qckeyword, calcindex)
    slabel = qc.get_slabel(s)
    parameters['all results'][slabel].update({qlabel:{}})
    parameters['qlabel'] = qlabel
    parameters['slabel'] = slabel
    parameters['tmpdir']  = io.join_path(*[scratch,ob.get_smiles_filename(s)])
    results = parameters['results']
    smilesname = ob.get_smiles_filename(s)
    smilesdir = io.join_path(parameters['database'], formula, smilesname)
    rundir  = io.join_path(*[smilesdir,parameters['qcdirectory']])
    qcpackage = parameters['qcpackage']
    qcscript = io.get_path(parameters['qcscript'])
    qcoutput = formula + '.out'
    parameters['smilesname' ] = smilesname
    parameters['smilesdir'] = smilesdir
    parameters['qcoutput'] = qcoutput
    if parameters['qctemplate']:
        parameters['qctemplate'] = io.get_path(parameters['qctemplate'])
        parameters['bacdirectory'] = parameters['qctemplate']
    if io.check_dir(parameters['qctemplate']):
        pass
    else:
        parameters['qctemplate'] = io.join_path(*[parameters['qtcdirectory'],'templates'])
    if not parameters['qckeyword']:
        runqc = False
        parseqc = False
        runthermo = False      
    if task=='composite':
        parameters['qctemplate'] = ''
    if parameters['writefiles'] or parameters['runthermo']:
        parameters['parseqc'] = True

    msg  = 'Mol. index   = {0}\n'.format(parameters['mol_index'])
    msg += 'Calc index   = {0}\n'.format(calcindex)
    msg += "Formula      = {0}\n".format(formula)
    msg += "SMILES       = {0}\n".format(s)
    msg += "InChI        = {0}\n".format(inchi)
    msg += "Multiplicity = {0}\n".format(mult)
    msg += "N_atoms      = {0}\n".format(natom)
    msg += "N_rotor      = {0}\n".format(nrotor)
    msg += "N_allrotor   = {0}\n".format(parameters['nallrotor'])
    msg += "N_methyl     = {0}\n".format(parameters['nmethyl'])
    msg += 'Task         = {0}\n'.format(task)
    msg += 'Method       = {0}\n'.format(parameters['qcmethod'])
    msg += 'Basis        = {0}\n'.format(parameters['qcbasis'])
    msg += 'Package      = {0}\n'.format(qcpackage)
    msg += 'QLabel       = {0}\n'.format(qlabel)
    msg += 'TemplateDir  = {0}\n'.format(parameters['qctemplate'])
    msg += 'RunDir       = {0}\n'.format(rundir)
    msg += 'XYZ          = \n{0}\n'.format(xyz)
    logging.info(msg)
    cwd = io.pwd()
    io.mkdir(rundir)
    if io.check_dir(rundir, 1):
        io.cd(rundir)
    else:
        logging.error('I/O error, {0} directory not found.\n'.format(rundir))
        return -1
    available_packages=['nwchem', 'molpro', 'mopac', 'gaussian','qchem']
    runfile = 'RUNNING.tmp'
    if io.check_file(runfile):
        if over:
            runqc = True
            logging.info('Overwriting the calculation...')
        elif task is 'composite':
            runqc = True
            logging.info('Composite calculation...')            
        else:
            runqc = False
            parseqc = False
            runthermo = False
            if 'opt' in task:
                parameters['break'] = True
            logging.warning('Skipping calculation since it is already running. Use -O {1}.{2} to overwrite or delete "{0}" file'.format(io.get_path(runfile),parameters['mol_index'],calcindex))
    elif ignore and task is not 'composite':
        runqc = False
    if task is 'composite':
        runqc = True
    runtime = 0
    xyzfilename = formula+'_initial.xyz'
    io.write_file(xyz, xyzfilename)
    parameters['xyzpath']=io.get_path(xyzfilename)
    if runqc:
        io.touch(runfile)
        try:
            if qcpackage in available_packages:
                runstart = timer()
                qc.run(s, parameters, mult=mult)
                runtime = timer() - runstart
                logging.debug("Runtime = {:15.3f} s".format(runtime))
                if runtime > 1.0:
                    logging.info("Runtime = {:15.3f} s".format(runtime))
                io.rm(runfile)
            elif task == 'composite':
                qc.run_composite(parameters)
                io.rm(runfile)
            elif qcpackage == 'qcscript':
                geofile = formula + '.geo'
                geo = ''.join(parameters['xyz'][2:natom+2]) 
                io.write_file(geo, geofile)
                if io.check_file(geofile, 1):
                    qc.run_qcscript(qcscript, parameters['qctemplate'], geofile, mult)
            elif parameters['qcmethod'] == 'given':
                qc.use_given_hof(parameters)
                io.rm(runfile)
            else:
                logging.error('{0} package not implemented.\nAvailable packages are {1}'.format(qcpackage,available_packages))
        except KeyboardInterrupt:
            logging.error('CTRL+C command...')
            logging.info('Deleting lock file {}'.format(io.get_path(runfile)))
            io.rm(runfile)
            sys.exit()
        except Exception as e:
            if parameters['debug']:
                raise
            else:
                logging.error('Error in running quantum chemistry calculations')
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                logging.error('Exception {}: {} {} {}'.format( e, exc_type, fname, exc_tb.tb_lineno))         
    parameters['overwrite'] = overwrite
    if parseqc:
        logging.info('Parsing output...')
        if parameters['qctask'] == 'composite':
            enefile = formula + '.ene'
            if io.check_file(enefile,0,True):
                energy = float(io.read_file(enefile).strip())
                parameters['results']['energy'] = energy
            else:
                if runthermo:
                    logging.error('Cannot run thermo, no energy file "{0}".\n'.format(enefile))
                    runthermo = False
                else:
                    logging.debug('No energy file "{0}".\n'.format(enefile))
        elif parameters['qcmethod'] == 'given':
            enefile = formula + '.ene'
            if io.check_file(enefile,0,True):
                energy = float(io.read_file(enefile).strip())
                parameters['results']['energy'] = energy
        elif io.check_file(qcoutput, timeout=1,verbose=False):
            out = io.read_file(qcoutput, aslines=False)
            if qc.check_output(out):
                try:
                    results = qc.parse_output(out,formula,parameters['writefiles'])
                except Exception as e:
                    if 'opt' in task:
                        parameters['break'] = True
                    if parameters['debug']:
                        raise
                    else:
                        logging.error('Error in parsing {}'.format(io.get_path(qcoutput)))
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        logging.error('Exception {}: {} {} {}'.format( e, exc_type, fname, exc_tb.tb_lineno))        
                for key in results.keys():
                    val = results[key]
                    if hasattr(val, '__iter__'):
                        if len(list(val))>0:
                            parameters['results'][key] = results[key]
                    else:
                        if val:
                            parameters['results'].update({key: results[key]})
            else:
                if 'opt' in task:
                    parameters['break'] = True
                if runthermo:
                    logging.error('Cannot run thermo failed calculation in "{0}"'.format(qcoutput))
                    runthermo = False
                else:
                    logging.error('Found failed calculation "{0}"'.format(qcoutput))
            if 'xyz' in results and natom > 1:
                final_xyz  = parameters['results']['xyz']
                parameters['xyzpath']=io.get_path(xyzfilename)
                try:
                    inchifinal = ob.get_inchi(final_xyz)
                except Exception as e:
                    logging.error('Exception {}. Final xyz can not be converted to INCHI -- Using original inchi'.format(e))       
                    inchifinal = inchi 
                logging.info('Final xyz = \n {}'.format(final_xyz))
                if inchi.strip() == inchifinal.strip():
                    parameters['xyz'] = final_xyz
                    mol = ob.get_mol(final_xyz)
                else:
                    logging.error('InChI mismatch: \n{} --> \n{}'.format(inchi,inchifinal))
            else:
                if 'opt' in task:
                    logging.error('Cannot find xyz, optimization failed')
                    parameters['break'] = True
                    runthermo = False
############################ BLUES specific
            if 'Hessian' in parameters['results'] and 'RPHt input' in parameters['results']:
                RPHtexe = '/lcrc/project/PACC/codes/EStokTP/exe/RPHt.exe'
                if not task.startswith('tors') and not task.startswith('md'):
                    RPHt, geolines = parameters['results']['RPHt input'].split('geometry')
                    geolines, gradlines = geolines.split('gradient')
                    xyz = io.read_file(formula + '_initial.xyz').splitlines()[2:]
                    RPHt += 'geometry\n'
                    for i, line in enumerate( geolines.splitlines()[1:]):
                        RPHt += '\t' + '\t'.join(line.split()[0:])+ '\n'
                        #RPHt += '\t' + '\t'.join(line.split()[0:3]) + '\t' + '\t'.join(xyz[i].split()[1:]) + '\n'
                    RPHt += 'gradient' +  gradlines.split('Hessian')[0] + 'Hessian\n' + parameters['results']['Hessian']
                    parameters['results']['RPHt input'] = RPHt
                    RPHtfile = 'RPHt_input_data.dat'
                if io.check_file(RPHtexe):
                    if not task.startswith('tors') and not task.startswith('md'):
                        io.write_file(RPHt,RPHtfile)
                        io.execute(RPHtexe  + ' ' + RPHtfile)
                    #if io.check_file( 'hrproj_freq.dat'):
                    if io.check_file( 'me_files/reac1_fr.me'):
                        out = io.read_file('me_files/reac1_fr.me', aslines=False)
                        pfreqs = qc.get_mess_frequencies(out)
                        parameters['results']['pfreqs'] = pfreqs
                    elif io.check_file( 'hrproj_freq.dat'):
                        pfreqs = io.read_file('hrproj_freq.dat').split(' 0.00')[0].split('\n')[:-1]
                        parameters['results']['pfreqs'] = pfreqs
                    else:
                        logging.warning('hrproj_freq.dat file not found')
                else:
                    logging.warning('{} not found.'.format(RPHtexe))
#########################
        else:
            if 'opt' in task:
                parameters['break'] = True
            logging.error('Output file "{0}" not found in {1}.'.format(qcoutput,io.pwd()))
            if runthermo:
                logging.error('Cannot run thermo')
                runthermo = False
        parameters['all results'][slabel][qlabel]['energy'] = float('nan')
        if parameters['bac']:
            bonds = {}
            x2zinp = ''
            if io.check_file(formula + '.xyz') :
                x2zinp = (formula + '.xyz')
            elif 'xyz' in parameters['results'] and natom > 1:
                x2zinp = parameters['results']['xyz']
            logging.info('Running x2z for BAC bonds')
            if x2zinp and natom > 1:
                try:
                    x2zout = qc.run_x2z(x2zinp, parameters['x2z'])
                    bonds =  qc.get_x2z_bonds(x2zout)
                    logging.info('Bonds found = {}'.format(bonds))
                except:
                    logging.error('x2z run failed, no bonds will be corrected')
            else:
               bonds = {}
            bac  = hf.calc_bac(parameters, bonds, parameters['qlabel'])
            bondstr = ''
            for key in bonds:
                bondstr += '{} x{}\n'.format(key, bonds[key])
            io.write_file('{:.3f}\n{}'.format(bac, bondstr), formula + '.bac')
            parameters['all results'][slabel][qlabel]['bac'] = bac                   
        parameters['all results'][slabel][qlabel]['energy'] = float('nan')
    if 'zpve' in parameters['results']:
        parameters['all results'][slabel][qlabel]['zpve'] = parameters['results']['zpve']   
    else:
        parameters['all results'][slabel][qlabel]['zpve'] = 0.0
    parameters['all results'][slabel][qlabel]['path'] = rundir   
    if 'hindered potential' in parameters['results'] and (task.startswith('tors') or task.startswith('md')):
        tc.get_hindered_potential(parameters['results']['hindered potential'],report=parameters['debug'])
        
    #parameters['all results'][slabel]['mol_index'] = parameters['mol_index']  
    for key in results.keys():
        val = results[key]
        if hasattr(val, '__iter__'):
            if len(list(val))>0:
                parameters['all results'][slabel][qlabel][key] = results[key]
                if 'freqs' in key:
                    floatfreqs = sorted([float(freq) for freq in results[key]])
                    parameters['freqdir'] = parameters['qcdirectory']
                    logging.info('{:10s} = {}'.format(key,['{:6.1f}'.format(freq) for freq in floatfreqs]))
                    if any(freq < 0 for freq in floatfreqs) and runthermo:
                        runthermo = False
                        logging.error('Cannot run thermo')
                if 'pfreqs' in key:
                    floatfreqs = sorted([float(freq) for freq in results[key]])
                    parameters['freqdir'] = parameters['qcdirectory']
                    logging.info('{:10s} = {}'.format(key,['{:6.1f}'.format(freq) for freq in floatfreqs]))
                    if any(freq < 0 for freq in floatfreqs) and runthermo:
                        runthermo = False
                        logging.error('Cannot run thermo')

        else:
            if val:
                parameters['all results'][slabel][qlabel].update({key: results[key]})
                if key == 'energy' or 'zpve' in key:
                    logging.info('{:10s} = {:10.8f}'.format(key,results[key]))
    parameters['results']['deltaH0'] = float('nan')
    parameters['all results'][slabel][qlabel]['deltaH0'] = float('nan')   
    parameters['results']['deltaH298'] = float('nan')
    parameters['all results'][slabel][qlabel]['deltaH298'] = float('nan')                   
    parameters['all results'][slabel][qlabel]['chemkin'] = ''
    if runtime > 1:
        parameters['all results'][slabel][qlabel]['runtime'] = runtime
    bonds = True
    if runthermo:
        sym = parameters['symm']
        hfcoeff = []
        hfset = []
        if sym:
            pass
        else:
            sym = 1
            if natom == 1:
                logging.info('Single atom, sym set to 1.')
                pass
            else:
                x2zinp = ''
                if io.check_file(formula + '.xyz') :
                    x2zinp = (formula + '.xyz')
                elif 'xyz' in parameters['results'] and natom > 1:
                    x2zinp = parameters['results']['xyz']
                logging.info('Running x2z for symmetry number')
                if x2zinp and natom > 1:
                    try:
                        x2zout = qc.run_x2z(x2zinp, parameters['x2z'])
                        sym = qc.get_x2z_sym(x2zout) 
                        sym = sym / ob.get_ent(s) 
                        bonds =  qc.get_x2z_bonds(x2zout)
                        logging.info('Symmetry number = {}'.format(sym))
                    except:
                        logging.error('x2z run failed, sym. number is set to 1. Probably a failed xyz')
                else:
                    bonds = {}
                    logging.error('xyz file cannot be found')
        parameters['results']['sym'] = sym
        hof = parameters['hof']
        if hof and parameters['qcmethod'] == 'given':
            hfset  = ['N/A (Precomputed)']
            hftxt  = 'Energy (kcal/mol)\tBasis\n----------------------------------'
            hftxt += '\n' + str(hof) + '\t' + '  '.join(hfset) 
        else:
            if formula in ['XH2','XO2','XN2']: #Remove X to use the definided values
                hof = 0.
                hfset = 'Definition'
                logging.info('Heat of formation of {} is set to 0 by definition.'.format(formula))
            else:
                hof, hfset, hfcoeff = hf.main_keyword(s,parameters)
            hftxt  = 'Energy (kcal/mol)\tBasis\n----------------------------------'
            hftxt += '\n' + str(hof) + '\t' + '  '.join(hfset) 
        io.write_file(hftxt,formula + '.hofk')
        parameters['results']['deltaH0'] = hof
        parameters['results']['heat of formation basis'] = hfset
        parameters['all results'][slabel][qlabel]['deltaH0'] = hof
        parameters['all results'][slabel][qlabel]['heat of formation basis'] = hfset
        if not type(hfcoeff) == list:
            hfcoeff = hfcoeff.tolist()
        parameters['all results'][slabel][qlabel]['heat of formation coeff'] = hfcoeff
        hof298 = 0.
        chemkintext = ''
        rmgpoly = {}
        pfout   = ''
        if formula in ['XH2','XO2','XN2']: #Remove X to use the definided values
            logging.info('Heat of formation of {} is set to 0 by definition.'.format(formula))
        else:
            if not io.check_file('new.groups'):
                groupstext = tc.get_new_groups()
                io.write_file(groupstext, 'new.groups')
            try:
                hof298, chemkintext, rmgpoly = tc.write_chemkin_polynomial(mol, parameters)
                if io.check_file('pf.dat'):
                    pfout = io.read_file('pf.dat')
                else:
                    pfout = ''
                ####STILL WORKING ON NEXT SECTION (SENSITIVITY/UNCERTAINTY DATA)
                ##(1) Scale HoF
                ##(2) Fix DoF to rovibrational from vibrational for QMCPSI, and set c etc
                ##(3) Scale Vibs should affect HoF
                if parameters['uncertainty']:
                   io.mkdir('uncertainty')
                   io.cd('uncertainty')
                   uncertainty = parameters['uncertainty'].split(',')
                   for uncert in uncertainty:
                       uncert = uncert.split('-')
                       if len(uncert) > 2:
                           typ,lscale, hscale = uncert
                       elif(len(uncert) > 1):
                           typ,lscale = uncert
                       else:
                           typ = uncert[0]
                       types = {'f':'freqs','h':'hindered potential'}
                       if typ[0].lower() in types:
                           if types[typ[0].lower()] in parameters['results']:
                               lscale = float(lscale)
                               hscale = float(hscale)
                               scales = np.arange(lscale,hscale,(hscale-lscale) / 9)
                               scales = np.append(scales,scales[-1] +(hscale-lscale) / 9)
                               for scale in scales:
                                   parameters['scale'] = scale
                                   parameters['scaletype'] = typ
                                   io.mkdir('scale{}_{:.3f}'.format(typ,scale))
                                   io.cd('scale{}_{:.3f}'.format(typ,scale))
                                   if not io.check_file('new.groups'):
                                       groupstext = tc.get_new_groups()
                                       io.write_file(groupstext, 'new.groups')
                                   _, tempchemkintext,_ = tc.write_chemkin_polynomial(mol, parameters)
                                   parameters['all results'][slabel][qlabel]['scalechemkin_{}-{:.3f}'.format(typ,scale)] = tempchemkintext
                                   parameters['scaletype'] = None
                                   io.cd('..')
                       elif typ.lower().startswith('q'):
                           if pfout:
                               dof = 3*natom - 6
                               if 'moltype' in parameters:
                                  if parameters['moltype'].lower() == 'linear':
                                      dof += 1
                               c = float(lscale)
                               io.mkdir('scale{}_{:.3f}'.format(typ,c))
                               io.cd('scale{}_{:.3f}'.format(typ,c))
                               hdr = '\n'.join(pfout.splitlines()[:2])
                               pfnew = pfout.splitlines()[2:]
                               for l, line in enumerate(pfnew):
                                   T, zeroth, first, second = line.split()
                                   T, zeroth, first, second = float(T), float(zeroth), float(first), float(second)
                                   newzeroth = zeroth + dof * np.log(1 + c * T / 1000)
                                   newfirst = first + dof / (1000 / c + T)
                                   newsecond = second + dof / (1000 /c + T )**2
                                   pfnew[l] = '{:3.1f}\t{:.4f}\t{:.8f}\t{:.5e}'.format(T, newzeroth, newfirst, newsecond)
                               pfnew = hdr + '\n' + '\n'.join(pfnew)
                               io.write_file(pfnew, 'pf.dat')
                               if not io.check_file('new.groups'):
                                   groupstext = tc.get_new_groups()
                                   io.write_file(groupstext, 'new.groups')
                               skippf = parameters['skippf']
                               parameters['skippf'] = True
                               _, tempchemkintext,_ = tc.write_chemkin_polynomial(mol, parameters)
                               parameters['skippf'] = skippf
                               parameters['all results'][slabel][qlabel]['scalechemkin_{}-{:.3f}'.format(typ,c)] = tempchemkintext
                               parameters['scaletype'] = None
                               io.cd('..')
                           else:
                               logging.info('Cannot scale Q, no pf.dat in {}'.format(os.path.sep.join(io.pwd().split(os.path.sep)[:-1])))
                       else:
                           logging.info('Invalid uncertainty analysis type {}'.format(typ))
                       io.cd('..')
                   parameters['scale'] = 0
                   #######END UNCERTAINTY SECTION
            except Exception as e:
                if parameters['debug']:
                    raise
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                logging.error('Failed in chemkin polynomial generation')
                logging.error('Exception {}: {} {} {}'.format(e, exc_type, fname, exc_tb.tb_lineno))         
        parameters['results']['deltaH298'] = hof298
        parameters['all results'][slabel][qlabel]['deltaH298'] = hof298   
        parameters['all results'][slabel][qlabel]['chemkin'] = chemkintext
        parameters['all results'][slabel][qlabel]['NASAPolynomial'] = rmgpoly
        parameters['all results'][slabel][qlabel]['partition function'] = pfout
     
    io.cd(cwd)
    return


def main(arg_update={}):
    from socket import gethostname
    from timeit import default_timer as timer
    import os
    from time import strftime as get_date_time
    global parameters
    mpirank = io.get_mpi_rank()
    mpisize = io.get_mpi_size(default=1)
    start  = timer()
    args = get_args()
    parameters = collections.OrderedDict()
    parameters.update(vars(args))
    qcnproc = parameters['qcnproc']
    beginindex = parameters['first']
    if parameters['swift']:
        qtcruninfo = 'Running QTC with Swift...'
        if mpirank:
            sys.exit()
    else:
        if mpisize > 1:
            qtcruninfo = 'Running QTC in parallel...'
            if (mpirank % qcnproc) == 0:
                qcrank = mpirank // qcnproc
                qtcruninfo += 'Using {} MPI ranks in total'.format(mpisize)
                qtcruninfo += 'Using {} cores for a single quantum chemistry calculation'.format(qcnproc)
                beginindex = parameters['first'] + qcrank
                qtcruninfo += 'Begin index shifted to {} for qcrank {} mpirank {}'.format(beginindex,qcrank, mpirank)
            else:
                sys.exit()
        else:
            qtcruninfo = 'Running QTC...'
    
    endindex = parameters['last']   
    parameters['all results'] = collections.OrderedDict()
    logfile = parameters['logfile']
    logindex = parameters['logindex']
    hostname = gethostname()
    loglevel = logging.INFO
    if parameters['loglevel'] == -1:
        if mpisize > 1  and mpirank > 0:
            loglevel = logging.ERROR
        else:
            loglevel = logging.INFO
    elif parameters['loglevel'] == 0:
        loglevel = logging.ERROR
    elif parameters['loglevel'] == 1:
        loglevel = logging.WARNING
    elif parameters['loglevel'] == 2:
        loglevel = logging.INFO
    elif parameters['loglevel'] > 2:
        loglevel = logging.DEBUG
    if parameters['debug']:
        loglevel = logging.DEBUG
    logging.addLevelName(logging.INFO, '')
    logging.addLevelName(logging.DEBUG, 'Debug:')
    logging.addLevelName(logging.ERROR, 'ERROR:')
    logging.addLevelName(logging.WARNING, 'WARNING:')
    if not logindex and mpisize > 1:
        logindex = str(mpirank)
    if logfile is 'none':
        if logindex:
            logfile = 'qtc_' + logindex + '_' + hostname + '.log'
        else:
            logging.basicConfig(format='%(levelname)s%(message)s', level=loglevel)
    else:
        logfile = logfile + logindex + '_' + hostname + '.log'
        logfile = logfile.replace('DATE', get_date_time("%y%m%d_%H%M%S"))
        logfile = io.get_unique_filename(logfile)
        logging.basicConfig(format='%(levelname)s%(message)s', filename=logfile, level=loglevel)
    for key in arg_update:
        parameters[key] = arg_update[key]
    logging.info(__logo__)
    logging.info("QTC: Date and time           = {0}".format(io.get_date()))
    logging.info("QTC: Last update             = {0}".format(__updated__))
    logging.info("QTC: Hostname                = {0}".format(hostname))
    logging.info(qtcruninfo)
    logging.info('QTC: Given arguments         =')
    for param in parameters:
        if not hasattr(args, param):
            setattr(args, param, {})
        logging.info('                             --{0:20s}\t{1}'.format(param, getattr(args, param)))
    if parameters['qckeyword']:
        parameters['qckeyword'] = qc.fix_qckeyword(parameters['qckeyword'])
        ncalc = len(parameters['qckeyword'].split(parameters['task_seperator']))
    else:
        ncalc = 1
    parameters['qtcdirectory'] = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    parameters['number_of_calculations'] = ncalc
    parameters['optlevel'] = 'sp' #TODO
    templatedir = parameters['qctemplate']
    if not os.path.sep in parameters['xyzdir']:
        parameters['xyzdir'] = io.join_path(*[io.pwd(),parameters['xyzdir']])
    inp = args.input
    jsonfile = args.jsoninput
    jlist = []
    if io.check_file(inp):
        if inp.split('.')[-1] == 'json':
            jlist = db.load_json(inp)
            mylist = qc.get_slabels_from_json(jlist)
        else:
            mylist = io.read_list2(inp)
    else:
        mylist = inp.split(',')
    if endindex:
        mylist = mylist[beginindex-1:endindex]
    else:
        mylist = mylist[beginindex-1:]
    # Convert to open-babel canonical smiles and add  _mX for multiplicity if not specified.
    mylist = qc.update_smiles_list(mylist)
    init = timer()
    logging.info("QTC: Initialization time (s) = {0:.2f}".format(init-start))
    runthermo = parameters['runthermo']
    if runthermo: 
        logging.info("QTC: Number of species       = {0}".format(len(mylist)))
        molbasis = []
        for s in mylist:
            s = qc.get_slabel(s)
            formula = ob.get_formula(s)
            ref = parameters['reference']
            if  'cbh' in ref.lower():
                _, basismolecules, _ = hf.cbh_coefficients(s.split('_')[0], ref, parameters)
            else: 
                _, basismolecules, _ = hf.comp_coefficients(formula, basis=parameters['reference'].split(','))
            for basismol in basismolecules:
                smi = qc.get_slabel(basismol)
                molbasis.append(basismol)
                if smi not in mylist:
                    msg = '{0} added to input list for heat of formation calculation of {1}'.format(basismol,s)
                    mylist = [smi] + mylist
                    logging.info(msg)
        molbasis = set(molbasis)
        molbasis = '\n{}'.format('\n'.join(molbasis))
        logging.info('QTC: Species required for ref{}{}'.format(ref, molbasis))
        logging.info("QTC: Number of species required for thermo= {0}".format(len(mylist)))
    if parameters['generate']:
        mylist = qc.sort_species_list(mylist, printinfo=True)
        myliststr = '\n'.join(mylist)
        sortedfile = 'sorted.txt'
        io.write_file(myliststr, sortedfile)
        if io.check_file(sortedfile,1):
            logging.info('Sorted SMILES file = {}'.format(sortedfile))
            logging.info('You can use qtc -b 1 -e 5, to compute species with indices 1,2,3,4,5.')
        else:
            logging.error('Problem in writing sorted SMILES file {}'.format(sortedfile))
    if parameters['sortbymass']:
        mylist = qc.sort_species_list(mylist, printinfo=True,byMass=True)
        myliststr = '\n'.join(mylist)
        sortedfile = 'sortedbymass.txt'
        io.write_file(myliststr, sortedfile)
        if io.check_file(sortedfile,1):
            logging.info('Sorted SMILES file = {}'.format(sortedfile))
        else:
            logging.error('Problem in writing sorted SMILES file {}'.format(sortedfile))
    if parameters['qckeyword']:
        logging.info('List of species')
        logging.info(pprint.pformat(mylist))
        for mid,s in enumerate(mylist):
            parameters['runthermo'] = False
            parameters['optdir'] = ''
            parameters['freqdir'] = ''
            parameters['anharmdir'] = ''
            parameters['qcdirectory'] = ''
            parameters['optlevel'] = ''
            parameters['freqlevel'] = ''
            parameters['mol_index'] = mid + 1
            parameters['break'] = False
            parameters['results'] = {}
            parameters['all results'].update({qc.get_slabel(s):{}}) 
            parameters = qc.add_species_info(s,parameters)
            for i in range(ncalc):
                if parameters['break']:
                    logging.info('Skipping next calculations for {}'.format(s))
                    break
                else:
                    parameters['calcindex'] = i
                    parameters['qctemplate'] = templatedir
                    logging.info('\n' + 100*'*' + '\n')
                    parameters = qc.parse_qckeyword(parameters, calcindex=i)
                    run(s)
        if runthermo:
            logging.info('\n' + 120*'#' + '\n')
            logging.info("Starting thermo calculations")
            for mid,s in enumerate(mylist):
                parameters['runthermo'] = runthermo
                parameters['runqc'] = False
                parameters['optdir'] = ''
                parameters['freqdir'] = ''
                parameters['anharmdir'] = ''
                parameters['qcdirectory'] = ''
                parameters['optlevel'] = ''
                parameters['freqlevel'] = ''
                parameters['mol_index'] = mid + 1
                parameters['break'] = False
                parameters['results'] = {}
                parameters['heat'] = None
                parameters = qc.add_species_info(s,parameters)
                for i in range(ncalc):
                    if parameters['break']:
                        logging.info('Skipping next calculations for {}'.format(s))
                        break
                    parameters['calcindex'] = i
                    parameters['qctemplate'] = templatedir
                    logging.info('\n' + 50*'-' + '\n')
                    parameters = qc.parse_qckeyword(parameters, calcindex=i)
                    run(s)            
    else:
        logging.info("You need to specify qckeyword with -k to run calculations")
    end = timer()
    if parameters['all results']:
        logging.info('\n' + 100*'-' + '\n')
        pathtitle = 'Path in {}'.format(parameters['database'])
        out   = '{0:5s} {1:30s} {2:>15s} {3:>15s}\t   {4}\n'.format('IDX','SMILES', 'Energy', 'ZPVE', pathtitle)
        out  += '{0:5s} {1:30s} {2:>15s} {3:>15s}\t   {4}\n'.format('   ','      ', '[Hartree]', '[Hartree]', '  ')
        for i,s in enumerate(mylist):
            sresults = parameters['all results'][qc.get_slabel(s)]
            for qcresultkey, qcresultval in sorted(sresults.iteritems(),key= lambda x: x[0]):
                runpath = qcresultval['path'].split('/database/')[-1]
                zpve = 0.
                if 'azpve' in qcresultval:
                    zpve = qcresultval['azpve']
                elif 'zpve' in qcresultval:
                    zpve = qcresultval['zpve']
                if not 'energy' in qcresultval:
                    qcresultval['energy'] = float('NaN') 
                out += '{0:5s} {1:30s} {2:15.5f} {3:15.5f}\t   {4}\n'.format(
                        str(i+1), qc.get_slabel(s), qcresultval['energy'],zpve,runpath)
        logging.info(out)
        logging.info('\n' + 100*'-' + '\n')
        if runthermo:
            ckin = ''
            out   = '{0:5s} {1:30s} {2:>15s} {3:>15s}\t   {4}\n'.format('IDX','SMILES', 'DeltaH(0)', 'DeltaH(298)', 'Key')
            out  += '{0:5s} {1:30s} {2:>15s} {3:>15s}\t   {4}\n'.format('   ','      ', '[kj/mol]', '[kj/mol]', '  ')
            for i,s in enumerate(mylist):
                sresults = parameters['all results'][qc.get_slabel(s)]
                for qcresultkey, qcresultval in sorted(sresults.iteritems(),key= lambda x: x[0]):
                    if qcresultval['deltaH298']:
                        out += '{0:5s} {1:30s} {2:15.5f} {3:15.5f}\t   {4}\n'.format(
                            str(i+1), qc.get_slabel(s),qcresultval['deltaH0']*ut.kcal2kj,qcresultval['deltaH298']*ut.kcal2kj,qcresultkey)
                        ckin += qcresultval['chemkin']
                    else:
                        out += s + '  not included in ckin because there is no pf output for ' + qcresultkey  + '\n'
            logging.info(out)
            out = ''
            if parameters['getBAC']:
                out += hf.get_bac(parameters, mylist, 0, 10)
                logging.info(out)
            ckinfile = 'chemkin_' + parameters['logfile'] + get_date_time("_%y%m%d_%H%M%S") + '.txt'
            ckinfile = io.get_unique_filename(ckinfile)
            io.write_file(ckin,ckinfile)
            logging.info('Written all chemkin polynomials in {}'.format(io.get_path(ckinfile)))
        if parameters['dumpjsonfile']:
            jsonfile = 'qtc_parameters_' +  get_date_time("%y%m%d_%H%M%S") + '.json'
            jsonfile = io.get_unique_filename(jsonfile)
            logging.info('Writing parameters json file {}'.format(jsonfile))
            db.dump_json(parameters, jsonfile)
            jsonfile = 'qtc_thermo_' +  get_date_time("%y%m%d_%H%M%S") + '.json'
            jsonfile = io.get_unique_filename(jsonfile)
            logging.info('Writing thermo json file {}'.format(jsonfile))
            db.dump_json(parameters['all results'], jsonfile)
            if jlist:
                prefix = inp.split('.')[0]
                for i in range(ncalc):
                    csvfile = prefix + '_method_' + str(i) +  get_date_time("_%y%m%d_%H%M%S") + '.csv'
                    csvfile = io.get_unique_filename(csvfile)
                    csvtext = '{},{},{},{},{},{},{},{},{},{}\n'.format(
                        'Slabel', 'RMGlabel', 'deltaH(0)', 'deltaH(298)', 'H298', 'S298', 'Cp(300)', 'Cp(500)','Cp(1000)', 'Cp(1500)')
                    for d in jlist:
                        name = str(d['name'])
                        smi  = str(d['SMILES'])
                        mult = int(d['multiplicity'])
                        qlabel = qc.get_qlabel(parameters['qckeyword'], i) 
                        s    = qc.get_slabel(smi,mult)
                        try:
                            thermoresults = parameters['all results'][slabel][qlabel]
                            deltaH0   = thermoresults['deltaH0']
                            deltaH298 = thermoresults['deltaH298']
                            poly      = thermoresults['NASAPolynomial']
                            Cplist    = [tc.get_heat_capacity(poly,T) for T in [300,500,1000,1500]]#cal/mol*K
                            S298      = tc.get_entropy(poly,298.15) #cal/mol*K
                            H298      = tc.get_enthalpy(poly,298.15) #kcal/mol
                            csvtext += '{},{},{},{},{},{},{},{},{},{}\n'.format(
                                s, name, deltaH0, deltaH298, H298, S298, Cplist[0], Cplist[1],Cplist[2], Cplist[3])
                        except:
                            csvtext += '{},{},{},{},{},{},{},{},{},{}\n'.format(
                                s, name,'NA',     'NA',      'NA', 'NA', 'NA',        'NA',    'NA', '   NA')
                    if csvtext:
                        logging.info('Writing csv file {}'.format(csvfile))
                        io.write_file(csvtext,csvfile)
                csvfile = prefix + '_rmg_' +  get_date_time("_%y%m%d_%H%M%S") + '.csv'
                csvfile = io.get_unique_filename(csvfile)
                csvtext = '{},{},{},{},{},{},{},{},{},{},{}\n'.format(
                    'Slabel', 'RMGlabel', 'Sensitivity', 'Uncertainty', 'Value', 'H298', 'S298', 'Cp(300)', 'Cp(500)','Cp(1000)', 'Cp(1500)')
                for d in jlist:
                    name = str(d['name'])
                    smi  = str(d['SMILES'])
                    mult = int(d['multiplicity'])
                    s    = qc.get_slabel(smi,mult)
                    try:
                        sensitivity = float(d['Sensitivity'])
                        uncertainty = float(d['Uncertainty'])
                        value       = float(d['Value'])
                        Cplist      = [float(cp) / ut.kcal2kj for cp in [d['Cp300'],d['Cp500'],d['Cp1000'],d['Cp1500']]]
                        S298        = float(d['S298']) / ut.kcal2kj
                        H298        = float(d['H298']) / ut.kcal2kj / 1000.
                        csvtext += '{},{},{},{},{},{},{},{},{},{},{}\n'.format(s, name, sensitivity, uncertainty, value, H298, S298, Cplist[0], Cplist[1],Cplist[2], Cplist[3])
                    except:
                        csvtext += '{},{},{},{},{},{},{},{},{},{},{}\n'.format(s, name,        'NA',        'NA',  'NA', 'NA', 'NA',      'NA',      'NA',     'NA',      'NA')
                if csvtext:
                    logging.info('Writing csv file {}'.format(csvfile))
                    io.write_file(csvtext,csvfile)
    logging.info("QTC: Calculations time (s)   = {0:.2f}".format(end - init))
    logging.info("QTC: Total time (s)          = {0:.2f}".format(end - start))
    logging.info("QTC: Date and time           = {0}".format(io.get_date()))

if __name__ == "__main__":
    main()

