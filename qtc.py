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
import sys
from patools import energy
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
"""
TODO: 
When optimization fails, restart from the last geometry
hof=0 at 0 K should also be at room T
Delete NWChem scratch files
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
                        default=1,
                        help='Number of processors for each quantum chemistry calculation')
    parser.add_argument('-n', '--nproc', type=int,
                        default=1,
                        help='Number of processors for qtc calculations, to run different species in parallel')
    parser.add_argument('-k', '--qckeyword', type=str,
                        default='',
                        help='Keyword string that defines quantum chemistry calculations i.e.: "opt/ccsd/cc-pvdz/gaussian,energy/ccsd/cc-pvtz/nwchem,extrapolation/cbs/energy=0.3*E0+0.7*E1" Note that each calculation is separated by a comma (,) and calculations are defined by TASK/METHOD/BASIS/PACKAGE. TASK can be opt, freq, anharm,extrapolation.METHOD and BASIS are simply copied into quantum chemistry input file as defined in the templates folder. PACKAGE can be gaussian, molpro or nwchem')
    parser.add_argument('-t', '--qctemplate', type=str,
                        default='',
                        help='Path for the templates directory. Templates have a specific format for filenames. See qtc/templates.')
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
    parser.add_argument('-f', '--first', type=int,
                        default=0,
                        help='Beginning index of the species list')
    parser.add_argument('-l', '--last', type=int,
                        help='Ending index of the species list')
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
    mol = ob.get_mol(s,make3D=True)
    mult = ob.get_mult(mol)
    formula = ob.get_formula(mol)
    nrotor = ob.get_nrotor(mol)
    natom = ob.get_natom(mol)
    label = qc.get_qc_label(natom, qckeyword, calcindex)
    parameters['nrotor'] = nrotor
    parameters['label'] = label
    if package in ['nwchem', 'molpro', 'mopac', 'gaussian', 'torsscan','torsopt' ]:
        if package.startswith('nwc'):
            parameters['qcexe'] = 'mpirun -n {0} nwchem'.format(qcnproc)
        elif package.startswith('mol'):
            parameters['qcexe'] = '{0} -n {1}'.format(parameters['molpro'],qcnproc)
        else:
            parameters['qcexe'] = parameters[package]
    if task.startswith('tors'):
        parameters['qcexe'] = parameters['torsscan']
    if io.check_dir(parameters['qctemplate']):
        tempdir = parameters['qctemplate']
    else:
        tempdir = io.join_path(*[parameters['qtcdirectory'],'templates'])
    if not parameters['qckeyword']:
        runqc = False
        parseqc = False
        runthermo = False      
    if task=='composite':
        parameters['qctemplate'] = ''
    elif task.startswith('tors'):
        templatename = task + '_template' + '.txt'
        parameters['qctemplate'] = io.join_path(*[tempdir,templatename])
    else:
        templatename = '{0}_{1}_template.txt'.format(task,package)
        templatename =  io.join_path(*[tempdir,templatename])
        if not io.check_file(templatename):
            templatename = '{0}_template.txt'.format(package)
            templatename =  io.join_path(*[tempdir,templatename])            
        parameters['qctemplate'] = templatename
    if parameters['writefiles']:
        parameters['parseqc'] = True
    if parameters['qctemplate']:
        parameters['qctemplate'] = io.get_path(parameters['qctemplate'])

    msg  = "Formula = {0}\n".format(formula)
    msg += "SMILES = {0}\n".format(s)
    msg += "Multiplicity = {0}\n".format(mult)
    msg += "Number of rotors = {0}\n".format(nrotor)
#    if nrotor == 0 and task.startswith('tors'):
#        msg += 'Passing {0} for {1}.'.format(parameters['qctask'], s)
#        return
    msg += 'Task = {0}\n'.format(parameters['qctask'])
    msg += 'Method = {0}\n'.format(parameters['qcmethod'])
    msg += 'Basis = {0}\n'.format(parameters['qcbasis'])
    msg += 'Package = {0}\n'.format(parameters['qcpackage'])
    msg += 'Label = {0}\n'.format(parameters['label'])
    msg += 'Template = {0}\n'.format(parameters['qctemplate'])
    smilesname = ob.get_smiles_filename(s)
    smilesdir =  ob.get_smiles_path(mol, mult, method='', basis='')
    smilesdir = io.join_path(parameters['database'], smilesdir)
    parameters['smilesdir'] = smilesdir
    workdirectory  = io.join_path(*[smilesdir,parameters['qcdirectory']])
    qcpackage = parameters['qcpackage']
    qcscript = io.get_path(parameters['qcscript'])
    qcoutput = smilesname + '_' + qcpackage + '.out'
    if task.startswith('tors'):
        qcoutput = smilesname + '_' + task + '.out'
    else:
        qcoutput = smilesname + '_' + qcpackage + '.out'
    xyzpath = parameters['xyzpath']
    xyzfile = ''
    if xyzpath:
        xyzfile = qc.find_xyzfile(xyzpath, smilesdir)
        msg += 'xyzpath = {0}\n'.format(xyzpath)
    elif optdir:
        xyzfilename = smilesname + '.xyz'
        msg += 'Optdir = {0}\n'.format(optdir)
        if io.check_file(io.join_path(*[smilesdir,optdir,xyzfilename])):
            xyzfile = io.join_path(*[smilesdir,optdir,xyzfilename])
    if xyzfile:
        msg += "XYZ file = '{0}'\n".format(xyzfile)
        mol = ob.get_mol(xyzfile)
    else:
        msg += "XYZ file not found in optdir '{0}' or xyzpath '{1}' \n".format(optdir,xyzpath)
    printp(msg)
    msg = ''
    cwd = io.pwd()
    io.mkdir(workdirectory)
    if io.check_dir(workdirectory, 1):
        io.cd(workdirectory)
        msg += "Run directory = '{0}'\n".format(workdirectory)
    else:
        msg += ('I/O error, {0} directory not found.\n'.format(workdirectory))
        return -1
    printp(msg)
    msg = ''
    available_packages=['nwchem', 'molpro', 'mopac', 'gaussian']          
    if runqc:
        if qcpackage in available_packages:
            msg += qc.run(mol, parameters, mult)
        elif task == 'composite':
            msg = qc.run_extrapolation(mol, parameters)
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
        printp(msg)
        msg = ''
    if parseqc:
        printp('Parsing output...')
        if io.check_file('geom1.xyz'):
            parameters['optdir'] = io.pwd()
        if parameters['qctask'] == 'composite':
            enefile = smilesname + '.ene'
            if io.check_file(enefile):
                energy = float(io.read_file(enefile).strip())
                parameters['results']['energy'] = energy
            else:
                msg += 'Can not find "{0}".\n'.format(enefile)
                msg += 'Can not run thermo\n'
                runthermo = False
        elif io.check_file(qcoutput, timeout=1,verbose=False):
            out = io.read_file(qcoutput,aslines=False)
            if qc.check_output(out):
                results = qc.parse_output(out,smilesname, parameters['writefiles'], parameters['storefiles'], parameters['optlevel'])
                for key in results.keys():
                    val = results[key]
                    if hasattr(val, "any"):
                        if val.any():
                            parameters['results'][key] = results[key]
                    else:
                        if val:
                            parameters['results'][key] = results[key]
            else:
                msg += 'Failed calculation in "{0}".\n'.format(qcoutput)
                msg += 'Can not run thermo\n'
                runthermo = False                
        else:
            msg += 'Output file "{0}" not found.\n'.format(qcoutput)
            msg += 'Can not run thermo\n'
            runthermo = False
    printp(msg)
    pprint(parameters['results'])
    parameters['all results'][s].update({label:parameters['results']})
    msg = ''                                   
    if runthermo:
        hof, hfset = hf.main_keyword(mol,parameters)
        hftxt  = 'Energy (kcal/mol)\tBasis\n----------------------------------'
        hftxt += '\n' + str(hof) + '\t' + '  '.join(hfset) 
        parameters['results']['heat of formation at 0 K'] = hof
        parameters['results']['heat of formation basis'] = hfset
        io.write_file(hftxt,s + '.hofk')
        if not io.check_file('new.groups'):
            groupstext = tc.get_new_groups()
            io.write_file(groupstext, 'new.groups')
#whoops i'll remove this after everything has merged fine
#        if task.startswith('comp'):
#            d = {}
#            freqs = []
#        else:
#            d = qc.parse_output(out,smilesname, parameters['writefiles'], parameters['storefiles'])
#            xyz = next(db.gen_dict_extract('xyz',d))
#            freqs = next(db.gen_dict_extract('harmonic frequencies',d))
#            hessian = next(db.gen_dict_extract('Hessian',d))
#            RPHt    = next(db.gen_dict_extract('RPHt input',d))
#            parameters['prjfreqs'] = next(db.gen_dict_extract('projected frequencies',d))
#            parameters['messhindered'] = next(db.gen_dict_extract('mess hindered input',d))
#            if RPHt:
#                parameters['RPHtinput'] = next(db.gen_dict_extract('RPHt input',d))
#            if hessian and 'RPHtinput' in parameters:
#                RPHt, geolines = parameters['RPHtinput'].split('geometry')
#                geolines, gradlines = geolines.split('gradient')
#                xyz = xyz.splitlines()[2:]
#                RPHt += 'geometry\n'
#                for i, line in enumerate( geolines.splitlines()[1:]):
#                    RPHt += '\t' + '\t'.join(line.split()[0:3]) + '\t' + '\t'.join(xyz[i].split()[1:]) + '\n'
#                RPHt += 'gradient' +  gradlines.split('Hessian')[0] + 'Hessian\n' + hessian
#                parameters['RPHtinput'] = RPHt
#                RPHtexe = '/lcrc/project/PACC/codes/EStokTP/exe/RPHt.exe'
#                RPHtfile = 'RPHt_input_data.dat'
#                io.write_file(RPHt,RPHtfile)
#                io.execute(RPHtexe  + ' ' + RPHtfile)
#            xmat   = next(db.gen_dict_extract('xmat',d))
#        hf0, hfset = hf.main_keyword(mol,parameters)
#        hftxt  = 'Energy (kcal/mol)\tBasis\n----------------------------------'
#        hftxt += '\n' + str(hf0) + '\t' + '  '.join(hfset) 
#        parameters['heat'] = hf0
#        io.write_file(hftxt,s + '.hf0k')
#        if len(freqs) > 0:
#            msg += tc.write_chemkin_polynomial(mol, xyz, freqs, hf0, parameters,xmat=xmat)
        if 'harmonic frequencies' in parameters['results']:
            msg += tc.write_chemkin_polynomial2(mol, parameters)
            if io.check_file('thermp.out'):
                import patools as pa
                lines = io.read_file('thermp.out')
                hof298 = pa.get_298(lines)
                parameters['results']['heat of formation at 298 K'] = hof298
                msg += 'delHf(298) = {0} kcal/mol'.format(hof298)
        parameters['all results'][s].update({label:parameters['results']})
    io.cd(cwd)
    printp(msg)
    return


def main(arg_update={}):
    from socket import gethostname
    from timeit import default_timer as timer
    import os
    global parameters
    start  = timer()
    args = get_args()
    parameters = vars(args)
    parameters['all results'] = {}
    for key in arg_update:
        parameters[key] = arg_update[key]
    if parameters['qckeyword']:
        parameters['qckeyword'] = qc.update_qckeyword(parameters['qckeyword'])
        ncalc = len(parameters['qckeyword'].split(','))
    else:
        ncalc = 1
    parameters['qtcdirectory'] = os.path.dirname(os.path.realpath(__file__))
    parameters['number_of_calculations'] = ncalc
    parameters['optlevel'] = 'sp' #TODO
    beginindex = args.first
    endindex = args.last
    inp = args.input
    jsonfile = args.jsoninput
    nproc = args.nproc
    if io.check_file(inp):
        mylist = io.read_list(inp)
    elif io.check_file(jsonfile):
        mylist = db.get_smiles_from_json(jsonfile)  
    else:
        mylist = inp.split(',')
    if endindex:
        mylist = mylist[beginindex:endindex]
    else:
        mylist = mylist[beginindex:]        
    printp(__logo__)
    printp("QTC: Date and time           = {0}".format(io.get_date()))
    printp("QTC: Last update             = {0}".format(__updated__))
    printp("QTC: Number of processes     = {0}".format(nproc))
    printp("QTC: Hostname                = {0}".format(gethostname()))
    printp('QTC: Given arguments         =')
    for param in parameters:
        printp('                             --{0:20s}\t{1}'.format(param, getattr(args, param)))
    init = timer()
    printp("QTC: Initialization time (s) = {0:.2f}".format(init-start))
    runthermo = parameters['runthermo']
    if runthermo:
        printp("QTC: Initial number of species       = {0}".format(len(mylist)))
        for s in mylist:
            formula = ob.get_formula(s)
            _, basismolecules, _ = hf.comp_coefficients(formula, basis=parameters['hfbasis'].split())
            for basismol in basismolecules:
                if basismol not in mylist:
                    msg = '{0} added to input list for heat of formation calculation of {1}'.format(basismol,s)
                    mylist = [basismol] + mylist
                    printp(msg)
    printp("QTC: Number of species       = {0}".format(len(mylist)))            
    if parameters['qckeyword']:
        if nproc == 1:
            for s in mylist:
                parameters['runthermo'] = False
                parameters['optdir'] = ''
                parameters['freqdir'] = ''
                parameters['anharmdir'] = ''
                parameters['qcdirectory'] = ''
                parameters['optlevel'] = ''
                parameters['freqlevel'] = ''
                parameters['results'] = {}
                parameters['all results'].update({s:{}}) 
                mol = ob.get_mol(s,make3D=True)
                parameters['natom'] = ob.get_natom(mol)
                for i in range(ncalc):
                    parameters['calcindex'] = i                        
                    printp('\n' + 100*'*' + '\n')
                    printp('Running QTC')
                    parameters = qc.parse_qckeyword(parameters, calcindex=i)
                    run(s)
        else:
            for i in range(ncalc):
                parameters['calcindex'] = i
                if parameters['qckeyword']:
                    qc.parse_qckeyword(parameters, calcindex=i)
                pool = multiprocessing.Pool(nproc)
        if runthermo:
            printp('\n' + 100*'#' + '\n')
            printp("Starting thermo calculations")
            for s in mylist:
                parameters['runthermo'] = runthermo
                parameters['optdir'] = ''
                parameters['freqdir'] = ''
                parameters['anharmdir'] = ''
                parameters['qcdirectory'] = ''
                parameters['optlevel'] = ''
                parameters['freqlevel'] = ''
                parameters['heat'] = None
                mol = ob.get_mol(s,make3D=True)
                parameters['natom'] = ob.get_natom(mol)
                for i in range(ncalc):
                    parameters['calcindex'] = i
                    printp('\n' + 100*'*' + '\n')
                    printp('Running QTC')
                    parameters = qc.parse_qckeyword(parameters, calcindex=i)
                    run(s)            
    else:
        pprint(mylist)
        printp("You need to specify qckeyword with -k to run calculations")
    end = timer()
    if 'all results' in parameters:
        printp('\n' + 100*'-' + '\n')
        pprint('All results')
        pprint(parameters['all results'])
    printp("QTC: Calculations time (s)   = {0:.2f}".format(end - init))
    printp("QTC: Total time (s)          = {0:.2f}".format(end-start))
    printp("QTC: Date and time           = {0}".format(io.get_date()))

if __name__ == "__main__":
    main()

