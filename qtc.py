#!/usr/bin/env python
import argparse
import multiprocessing
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
from patools import energy
from timeit import default_timer as timer
__updated__ = "2017-10-27"
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
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     #formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
                        default='',
                        help='Path for the templates directory. Templates have a specific format for filenames. See qtc/templates.')
    parser.add_argument('-l', '--loglevel', type=int,
                        default=2,
                        help='Verbosity level of logging, 0 for only errors, 3 for debugging')
    parser.add_argument('-f', '--logfile', type=str,
                        default= 'none',
                        help='Log file prefix, use none for logging to STDOUT, include DATE if you want a date stamp')
    parser.add_argument('-s', '--scratch', type=str,
                        default=io.get_env('TMPDIR'),
                        help='Scratch directory. Default is obtained from TMPDIR env. variable.')
    parser.add_argument('-g', '--logindex', type=str,
                        default= '',
                        help='Log file index')
    parser.add_argument('-d', '--database', type=str,
                        default= io.pwd(),
                        help='Path for database directory')
    parser.add_argument('-o', '--qcoutput', type=str,
                        default='',
                        help='Path for the qc output file')
    parser.add_argument('-b', '--first', type=int,
                        default=1,
                        help='Beginning index of the species list')
    parser.add_argument('-e', '--last', type=int,
                        help='Ending index of the species list')
    parser.add_argument('-B', '--hfbasis', type=str,
                        default='auto',
                        help='List of SMILES (seperated by commas) for reference thermo species')
    parser.add_argument('-G', '--generate', action='store_true',
                        help='Generates a sorted list of species')
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
    parser.add_argument('-I', '--ignorerunningjobs', action='store_true',
                        help='Ignores RUNNING.tmp files and run the calculations')
    parser.add_argument('-R', '--recover', action='store_true',
                        help='Attempt to recover failed calculations')
    parser.add_argument('-O', '--overwrite', action='store_true',
                        help='Overwrite existing calculations. Be careful, data will be lost.')
    parser.add_argument('-X', '--excel', action='store_true',
                        help='Generate excel file')
    parser.add_argument('-J', '--dumpjsonfile', action='store_true',
                        help='Writes a json file containing all the results')
    parser.add_argument('-A', '--anharmonic', action='store_true',
                        help='Anharmonic corrections')
    parser.add_argument('-D', '--debug', action='store_true',
                        help='Run in debug mode')
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
    package = parameters['qcpackage']
    task    = parameters['qctask']
    runthermo = parameters['runthermo']
    qcnproc = parameters['qcnproc']
    qckeyword = parameters['qckeyword']
    calcindex = parameters['calcindex']
    optdir = parameters['optdir']
    ignore = parameters['ignorerunningjobs']
    overwrite=parameters['overwrite']
    scratch = parameters['scratch']
    machinefile=io.get_path(parameters['machinefile'])
    mol = ob.get_mol(s,make3D=True)
    mult = ob.get_mult(mol)
    formula = ob.get_formula(mol)
    nrotor = ob.get_nrotor(mol)
    natom = ob.get_natom(mol)
    label = qc.get_qc_label(natom, qckeyword, calcindex)
    parameters['all results'][s].update({label:{}})
    parameters['natom'] = natom
    parameters['nrotor'] = nrotor
    parameters['label'] = label
    parameters['slabel'] = s
    parameters['tmpdir']  = io.join_path(*[scratch,ob.get_smiles_filename(s)])
    results = parameters['results']
    if parameters['qctemplate']:
        parameters['qctemplate'] = io.get_path(parameters['qctemplate'])
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
    msg += "Formula      = {0}\n".format(formula)
    msg += "SMILES       = {0}\n".format(s)
    msg += "Multiplicity = {0}\n".format(mult)
    msg += "N_atoms      = {0}\n".format(natom)
    msg += "N_obrotors   = {0}\n".format(nrotor)
    msg += 'Task         = {0}\n'.format(parameters['qctask'])
    msg += 'Method       = {0}\n'.format(parameters['qcmethod'])
    msg += 'Basis        = {0}\n'.format(parameters['qcbasis'])
    msg += 'Package      = {0}\n'.format(parameters['qcpackage'])
    msg += 'Label        = {0}\n'.format(parameters['label'])
    msg += 'TemplateDir  = {0}\n'.format(parameters['qctemplate'])
    logging.info(msg)
    smilesname = ob.get_smiles_filename(s)
    parameters['smilesname' ] = smilesname
    smilesdir =  ob.get_smiles_path(s, mult)
    smilesdir = io.join_path(parameters['database'], smilesdir)
    parameters['smilesdir'] = smilesdir
    workdirectory  = io.join_path(*[smilesdir,parameters['qcdirectory']])
    qcpackage = parameters['qcpackage']
    qcscript = io.get_path(parameters['qcscript'])
    if task.startswith('tors'):
        qcoutput = smilesname + '_' + task + '.out'
    else:
        qcoutput = smilesname + '_' + qcpackage + '.out'
    parameters['qcoutput'] = qcoutput
    xyzfile = ''
    if optdir:
        xyzfilename = smilesname + '.xyz'
        if io.check_file(io.join_path(*[smilesdir,optdir,xyzfilename])):
            xyzfile = io.join_path(*[smilesdir,optdir,xyzfilename])
    if xyzfile and natom > 1:
        logging.info("XYZ file     = '{0}'".format(xyzfile))
        try:
            mol = ob.get_mol(xyzfile)
        except:
            logging.error('Not a valid xyz file {0}. Skipping following calculations.'.format(xyzfile))
            runqc = False
            parseqc = False
            runthermo = False
    elif 'results' in parameters.keys():
        results = parameters['results']
        if 'xyz' in results.keys():
            logging.info('Using a previously calculated xyz:\n {}'.format(results['xyz']))
            mol = ob.get_mol(results['xyz'])
    else:
        if natom == 1:
            logging.info("XYZ not required, single atom calculation")
        elif task.startswith('tors') or task.startswith('opt'):
            logging.info("XYZ file not found in optdir '{0}'".format(optdir))
        else:
            logging.error('No optimized geometry found, skipping subsequent calculations')
            runqc = False
            parseqc = False
            runthermo = False
    cwd = io.pwd()
    io.mkdir(workdirectory)
    if io.check_dir(workdirectory, 1):
        io.cd(workdirectory)
        logging.info("Rundir       = '{0}'".format(workdirectory))
    else:
        logging.error('I/O error, {0} directory not found.\n'.format(workdirectory))
        return -1
    available_packages=['nwchem', 'molpro', 'mopac', 'gaussian']
    runfile = 'RUNNING.tmp'
    if io.check_file(runfile):
        if overwrite or task is 'composite' :
            runqc = True
            logging.info('Running composite calculation...')
        else:
            runqc = False
            parseqc = False
            runthermo = False
            logging.warning('Skipping calculation since it is already running. Use -O to overwrite or delete "{}" file'.format(io.get_path(runfile)))
    elif ignore and task is not 'composite':
        runqc = False
    if task is 'composite':
        runqc = True
    runtime = 0
    if runqc:
        io.touch(runfile)
        try:
           # if natom > 1: #uncomment when x2z can exlude methyl rotors and can be installed on all plaforms
           #     if io.check_exe(parameters('x2z')):
           #         test_out = qc.run_x2z(ob.get_xyz(mol), parameters['x2z'])
           #         nrotor = qc.get_x2z_nrotor(test_out)
           #         #parameters['nrotor'] = nrotor # We may want to uncomment in the future
           #         logging.info("Number of rotors (x2z) = {0}\n".format(nrotor))
           #     else:
           #         logging.warning("x2z not found")
            if qcpackage in available_packages:
                runstart = timer()
                qc.run(mol, parameters, mult=mult)
                runtime = timer() - runstart
                logging.debug("Runtime = {:15.3f} s".format(runtime))
                if runtime > 1.0:
                    logging.info("Runtime = {:15.3f} s".format(runtime))
            elif task == 'composite':
                qc.run_composite(parameters)
            elif qcpackage == 'qcscript':
                geofile = smilesname + '.geo'
                geo = ob.get_geo(mol)
                io.write_file(geo, geofile)
                if io.check_file(geofile, 1):
                    qc.run_qcscript(qcscript, parameters['qctemplate'], geofile, mult)
            else:
                logging.error('{0} package not implemented.\nAvailable packages are {1}'.format(qcpackage,available_packages))
            io.rm(runfile)
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
    if parseqc:
        logging.info('Parsing output...')
        if io.check_file('geom1.xyz'):
            parameters['optdir'] = io.pwd()
        if parameters['qctask'] == 'composite':
            enefile = smilesname + '.ene'
            if io.check_file(enefile):
                energy = float(io.read_file(enefile).strip())
                parameters['results']['energy'] = energy
            else:
                parameters['break'] = True
                if runthermo:
                    logging.error('Cannot run thermo, no energy file "{0}".\n'.format(enefile))
                    runthermo = False
                else:
                    logging.debug('No energy file "{0}".\n'.format(enefile))
                    
        elif io.check_file(qcoutput, timeout=1,verbose=False):
            out = io.read_file(qcoutput, aslines=False)
            if qc.check_output(out):
                try:
                    results = qc.parse_output(out,smilesname,parameters['writefiles'],parameters['storefiles'],parameters['optlevel'])
                except Exception as e:
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
                parameters['break'] = True
                if runthermo:
                    logging.error('Cannot run thermo failed calculation in "{0}"'.format(qcoutput))
                    runthermo = False
                else:
                    logging.error('Found failed calculation "{0}"'.format(qcoutput))
            if 'xyz' in results or natom == 1:
                pass
            else:
                parameters['break'] = True
                if runthermo:
                    logging.error('Cannot run thermo, no xyz')
                    runthermo = False
############################ BLUES specific
            if 'Hessian' in parameters['results'] and 'RPHt input' in parameters['results']:
                RPHt, geolines = parameters['results']['RPHt input'].split('geometry')
                geolines, gradlines = geolines.split('gradient')
                xyz = parameters['results']['xyz'].splitlines()[2:]
                RPHt += 'geometry\n'
                for i, line in enumerate( geolines.splitlines()[1:]):
                    RPHt += '\t' + '\t'.join(line.split()[0:3]) + '\t' + '\t'.join(xyz[i].split()[1:]) + '\n'
                RPHt += 'gradient' +  gradlines.split('Hessian')[0] + 'Hessian\n' + parameters['results']['Hessian']
                parameters['results']['RPHtinput'] = RPHt
                RPHtexe = '/lcrc/project/PACC/codes/EStokTP/exe/RPHt.exe'
                RPHtfile = 'RPHt_input_data.dat'
                if io.check_file(RPHtexe):
                    io.write_file(RPHt,RPHtfile)
                    io.execute(RPHtexe  + ' ' + RPHtfile)
                    if io.check_file( 'hrproj_freq.dat'):
                        pfreqs = io.read_file('hrproj_freq.dat').split('0.0')[0].split('\n')[:-1]
                        parameters['results']['pfreqs'] = pfreqs
                    else:
                        logging.warning('hrproj_freq.dat file not found')
                else:
                    logging.warning('{} not found.'.format(RPHtexe))
#########################
        else:
            parameters['break'] = True
            logging.error('Output file "{0}" not found in {1}.'.format(qcoutput,io.pwd()))
            if runthermo:
                logging.error('Cannot run thermo')
                runthermo = False
    parameters['all results'][s][label]['energy'] = 0.
    if 'zpve' in parameters['results']:
        parameters['all results'][s][label]['zpve'] = parameters['results']['zpve']   
    else:
        parameters['all results'][s][label]['zpve'] = 0.
    parameters['all results'][s][label]['path'] = workdirectory   
    #parameters['all results'][s]['mol_index'] = parameters['mol_index']  
    for key in results.keys():
        val = results[key]
        if hasattr(val, '__iter__'):
            if len(list(val))>0:
                parameters['all results'][s][label][key] = results[key]
                if 'freqs' in key:
                    floatfreqs = sorted([float(freq) for freq in results[key]])
                    logging.info('{:10s} = {}'.format(key,['{:6.1f}'.format(freq) for freq in floatfreqs]))
                    if any(freq < 0 for freq in floatfreqs) and runthermo:
                        runthermo = False
                        logging.error('Cannot run thermo')

        else:
            if val:
                parameters['all results'][s][label].update({key: results[key]})
                if key == 'energy' or 'zpve' in key:
                    logging.info('{:10s} = {:10.8f}'.format(key,results[key]))
    parameters['results']['deltaH0'] = 0
    parameters['all results'][s][label]['deltaH0'] = 0   
    parameters['results']['deltaH298'] = 0
    parameters['all results'][s][label]['deltaH298'] = 0                                    
    parameters['all results'][s][label]['chemkin'] = ''
    if runtime > 1:
        parameters['all results'][s][label]['runtime'] = runtime
    if runthermo:
        sym = 1
        if natom == 1:
            logging.info('Single atom, sym set to 1.')
            pass
        else:
            test_inp = ''
            if io.check_file(smilesname + '.xyz') :
                test_inp = (smilesname + '.xyz')
            elif 'xyz' in parameters['results'] and natom > 1:
                test_inp = parameters['results']['xyz']
            logging.info('Running x2z for symmetry number')
            if test_inp:
                try:
                    out_x2z = qc.run_x2z(test_inp, parameters['x2z'])
                    sym = qc.get_x2z_sym(out_x2z) 
                    logging.info('Symmetry number = {}'.format(sym))
                except:
                    logging.error('x2z run failed, sym. number is set to 1. Probably a failed xyz')
            else:
                logging.error('xyz file cannot be found')
        parameters['results']['sym'] = sym
        if formula in ['H2','O2','N2']:
            hof = 0.
            hfset = 'Definition'
            logging.info('Heat of formation of {} is set to 0 by definition.'.format(formula))
        else:
            hof, hfset = hf.main_keyword(s,parameters)
        hftxt  = 'Energy (kcal/mol)\tBasis\n----------------------------------'
        hftxt += '\n' + str(hof) + '\t' + '  '.join(hfset) 
        parameters['results']['deltaH0'] = hof
        parameters['results']['heat of formation basis'] = hfset
        parameters['all results'][s][label]['deltaH0'] = hof
        parameters['all results'][s][label]['heat of formation basis'] = hfset
        io.write_file(hftxt,smilesname + '.hofk')
        hof298 = 0.
        chemkintext = ''
        rmgpoly = {}
        if formula in ['H2','O2','N2']:
            logging.info('Heat of formation of {} is set to 0 by definition.'.format(formula))
        else:
            if not io.check_file('new.groups'):
                groupstext = tc.get_new_groups()
                io.write_file(groupstext, 'new.groups')
            try:
                hof298, chemkintext, rmgpoly = tc.write_chemkin_polynomial(mol, parameters)
            except Exception as e:
                if parameters['debug']:
                    raise
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                logging.error('Failed in chemkin polynomial generation')
                logging.error('Exception {}: {} {} {}'.format(e, exc_type, fname, exc_tb.tb_lineno))         
        parameters['results']['deltaH298'] = hof298
        parameters['all results'][s][label]['deltaH298'] = hof298   
        parameters['all results'][s][label]['chemkin'] = chemkintext
        parameters['all results'][s][label]['NASAPolynomial'] = rmgpoly
    io.cd(cwd)
    return


def main(arg_update={}):
    from socket import gethostname
    from timeit import default_timer as timer
    import os
    from time import strftime as get_date_time
    global parameters
    mpirank = io.get_mpi_rank()
    if mpirank:
        sys.exit()
    start  = timer()
    args = get_args()
    parameters = vars(args)
    parameters['all results'] = {}
    logfile = parameters['logfile']
    logindex = parameters['logindex']
    hostname = gethostname()
    if parameters['loglevel'] == 0:
        loglevel = logging.ERROR
    elif parameters['loglevel'] == 1:
        loglevel = logging.WARNING
    elif parameters['loglevel'] == 2:
        loglevel = logging.INFO
    elif parameters['loglevel'] == 3:
        loglevel = logging.DEBUG
    if parameters['debug']:
        loglevel = logging.DEBUG
    logging.addLevelName(logging.INFO, '')
    logging.addLevelName(logging.DEBUG, 'Debug:')
    logging.addLevelName(logging.ERROR, 'ERROR:')
    logging.addLevelName(logging.WARNING, 'WARNING:')
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
    logging.info('QTC: Given arguments         =')
    for param in parameters:
        logging.info('                             --{0:20s}\t{1}'.format(param, getattr(args, param)))
    if parameters['qckeyword']:
        parameters['qckeyword'] = qc.update_qckeyword(parameters['qckeyword'])
        ncalc = len(parameters['qckeyword'].split(','))
    else:
        ncalc = 1
    parameters['qtcdirectory'] = os.path.dirname(os.path.realpath(__file__))
    parameters['number_of_calculations'] = ncalc
    parameters['optlevel'] = 'sp' #TODO
    templatedir = parameters['qctemplate']
    beginindex = args.first
    endindex = args.last
    inp = args.input
    jsonfile = args.jsoninput
    nproc = args.nproc
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
    mylist = qc.update_smiles_list(mylist)

    init = timer()
    logging.info("QTC: Initialization time (s) = {0:.2f}".format(init-start))
    runthermo = parameters['runthermo']
    if runthermo:
        logging.info("QTC: Number of species       = {0}".format(len(mylist)))
        for s in mylist:
            formula = ob.get_formula(s)
            _, basismolecules, _ = hf.comp_coefficients(formula, basis=parameters['hfbasis'].split(','))
            for basismol in basismolecules:
                smi = qc.get_slabel(basismol)
                if smi not in mylist:
                    msg = '{0} added to input list for heat of formation calculation of {1}'.format(basismol,s)
                    mylist = [smi] + mylist
                    logging.info(msg)
        logging.info("QTC: Number of species required for thermo= {0}".format(len(mylist)))
    if parameters['generate']:
        mylist = qc.sort_species_list(mylist, printinfo=True)
        myliststr = '\n'.join(mylist)
        sortedfile = io.join_path(io.get_path(inp, directory=True),'sorted.txt')
        io.write_file(myliststr, sortedfile)
        if io.check_file(sortedfile,1):
            logging.info('Sorted SMILES file = {}'.format(sortedfile))
            logging.info('You can use qtc -b 1 -e 5, to compute species with indices 1,2,3,4,5.')
        else:
            logging.error('Problem in writing sorted SMILES file {}'.format(sortedfile))
    elif parameters['qckeyword']:
        logging.info('List of species')
        logging.info(pprint.pformat(mylist))
        if nproc == 1:
            for mid,s in enumerate(mylist):
                parameters['runthermo'] = False
                parameters['optdir'] = ''
                parameters['freqdir'] = ''
                parameters['anharmdir'] = ''
                parameters['qcdirectory'] = parameters['database']
                parameters['optlevel'] = ''
                parameters['freqlevel'] = ''
                parameters['mol_index'] = mid + 1
                parameters['break'] = False
                parameters['results'] = {}
                parameters['all results'].update({s:{}}) 
                mol = ob.get_mol(s,make3D=True)
                parameters['natom'] = ob.get_natom(mol)
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
                parameters['qcdirectory'] = parameters['database']
                parameters['optlevel'] = ''
                parameters['freqlevel'] = ''
                parameters['mol_index'] = mid + 1
                parameters['break'] = False
                parameters['results'] = {}
                parameters['heat'] = None
                for i in range(ncalc):
                    if parameters['break']:
                        logging.info('Skipping next calculations for {}'.format(s))
                        break
                    mol = ob.get_mol(s,make3D=True)
                    parameters['natom'] = ob.get_natom(mol)
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
            sresults = parameters['all results'][s]
            for qcresultkey, qcresultval in sorted(sresults.iteritems(),key= lambda x: x[0]):
                runpath = 'database/' + qcresultval['path'].split('/database/')[-1]
                out += '{0:5s} {1:30s} {2:15.5f} {3:15.5f}\t   {4}\n'.format(
                        str(i+1), s, qcresultval['energy'],qcresultval['zpve'],runpath)
        logging.info(out)
        logging.info('\n' + 100*'-' + '\n')
        if runthermo:
            ckin = ''
            out   = '{0:30s} {1:>15s} {2:>15s}\t   {3:50s}\n'.format('SMILES', 'deltaH(0)', 'deltaH(298)', 'Key')
            out += '{0:30s} {1:>15s} {2:>15s}\t   {3:50s}\n'.format('      ', ' [kj/mol]', '[kj/mol]', ' ')
            for resultkey,resultval in parameters['all results'].iteritems():
                for qcresultkey, qcresultval in sorted(resultval.iteritems(),key= lambda x: x[0]):
                    out += '{0:30s} {1:15.5f} {2:15.5f}\t   {3:50s}\n'.format(
                        resultkey, qcresultval['deltaH0']*ut.kcal2kj,qcresultval['deltaH298']*ut.kcal2kj,qcresultkey)
                    ckin += qcresultval['chemkin']
            logging.info(out)
            ckinfile = 'chemkin_' + get_date_time("%y%m%d_%H%M%S") + '.txt'
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
                        qlabel = qc.get_qc_label(ob.get_natom(smi), parameters['qckeyword'], i)
                        s    = qc.get_slabel(smi,mult)
                        try:
                            thermoresults = parameters['all results'][s][qlabel]
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

