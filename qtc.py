#!/usr/bin/env python
import argparse
import subprocess
import multiprocessing

import iotools as io
import obtools as ob
import qctools as qc
import tctools as tc

__updated__ = "2017-06-13"
__author__ = "Murat Keceli"
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

# _mopacexe = 'mopac'
# _nwchemexe = 'nwchem'
# _gaussianexe = 'g09'
# _messpfexe = 'messpf'
# _thermpexe = 'thermp'
# _pac99exe = 'pac99'
# _qcmethod = 'pm3'
# _qccode = 'mopac'
# _runqc = False
# _runthermo = False


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
                        help='Number of processors, default is all processors')
    parser.add_argument('-m', '--qcmethod', type=str,
                        default='',
                        help='Quantum chemistry method to be used (obsolete, use --rundir)')
    parser.add_argument('-b', '--qcbasis', type=str,
                        default='',
                        help='Basis-set level in quantum chemistry calculations (obsolete, use --rundir)')
    parser.add_argument('-c', '--qccode', type=str,
                        default='mopac',
                        help='Quantum chemistry code ("gausian","mopac","qcscript") to be used')
    parser.add_argument('-t', '--qctemplate', type=str,
                        default='qctemplate.txt',
                        help='Template for quantum chemistry input file')
    parser.add_argument('-p', '--qcpath', type=str,
                        default='',
                        help='Path for the directory for running qc jobs')
    parser.add_argument('-x', '--xyzpath', type=str,
                        default='',
                        help='Path for the directory for xyz file')
    parser.add_argument('-Q', '--runqc', action='store_true',
                        help='Run quantum chemistry calculation')
    parser.add_argument('-S', '--subqc', action='store_true',
                        help='Submit quantum chemistry calculations')
    parser.add_argument('-T', '--runthermo', action='store_true',
                        help='Run thermochemistry calculations')
    parser.add_argument('-W', '--writefiles', action='store_true',
                        help='Write .xyz, .ene files')
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
                        default='/lcrc/project/PACC/test-awj/builddb/bin/qcscript.pl',
                        help='Path for qcscript perl script')
    return parser.parse_args()


def write_chemkin_polynomial(mol, method, zpe, xyz, freqs, deltaH):
    """
    A driver to perform all operations to write NASA polynomial in
    chemkin format. Assumes quantum chemistry calculation is performed.
    """
    messpfinput = 'pf.inp'
    thermpinput = 'thermp.dat'
    messpfoutput = 'pf.log'
    name = mol.formula
    tag = method
    inp = tc.get_pf_input(mol, method, zpe, xyz, freqs)
    io.write_file(inp, messpfinput)
    msg = 'Running {0} to generate partition function.\n'.format(_messpf)
    msg += io.execute([_messpf,messpfinput])
    msg += 'Running thermp .\n'
    inp = tc.get_thermp_input(mol.formula, deltaH)
    msg = tc.run_thermp(inp, 'thermp.dat', messpfoutput, _thermp) 
    msg += 'Running pac99.\n'
    msg += tc.run_pac99(name)
    msg += 'Converting to chemkin format.\n'
    chemkinfile = name + '.ckin'
    msg += 'Writing chemking file {0}.\n'.format(chemkinfile)
    try:
        msg += tc.write_chemkin_file(deltaH, tag, name, chemkinfile)
    except:
        "Failed to write polynomials"
    return msg


def run(s):
    """
    A driver function to run quantum chemistry and thermochemistry calculations based
    on command line options:
    --qcmethod
    --qccode
    """
    import qctools as qc
    import obtools as ob
    import iotools as io
    mult = 0
    msg = "***************************************\n"
    msg += "{0}\n".format(s)
    mult = ob.get_mult(s)
    mol = ob.get_mol(s)
    smilesname = ob.get_smiles_filename(s)
    smilesdir =  ob.get_smiles_path(mol, mult, method='', basis='', geopath='')
    qcpath  = io.join_path(*(smilesdir,_qcpath))
    xyzpath = None
    if _xyzpath:
        if io.check_dir(_xyzpath):
            xyzpath = _xyzpath
        elif io.check_dir(io.join_path(*(smilesdir,_xyzpath))):
            xyzpath = io.join_path(*(smilesdir,_xyzpath))
        else:
            msg += "xyz path not found {0}".format(_xyzpath)
    if xyzpath:                
        xyzfile = io.find_files(xyzpath, '*.xyz')
        try:
            xyzfile = next(xyzfile)
            xyz = io.read_file(xyzfile)
            coords = ob.get_coordinates_array(xyz)
            mol = ob.set_xyz(mol, coords)
        except StopIteration:
            msg += "xyz file not found in {0}".format(xyzpath)
    runfile = smilesname + '.run'
    runqc = _runqc
    runthermo = _runthermo
    runanharmonic = _anharmonic
    io.mkdir(qcpath)
    cwd = io.pwd()
    if io.check_dir(qcpath, 1):
        io.cd(qcpath)
        msg += "cd '{0}'\n".format(qcpath)
    else:
        msg += ('I/O error, {0} directory not found.\n'.format(qcpath))
        return -1
    if _qccode == 'mopac':
        qclog = smilesname + '.out'
    elif _qccode == 'gaussian'  :
        qclog = smilesname + '_gaussian.log'
    elif _qccode == 'nwchem'  :
        qclog = smilesname + '_nwchem.log'
    elif _qccode == 'molpro'  :
        qclog = smilesname + '_molpro.log'
    else:
        qclog = smilesname + '_qc.log'
    print(msg)
    msg = ''          
    if runqc:
        if io.check_file(runfile):
            if _overwrite:
                msg += 'Overwriting previous attempt.\n'
            else:
                msg += ('Skipping {0}\n'.format(smilesname))
        else:
            io.touch(runfile)
        if _qccode == 'mopac':
            msg += "Running mopac...\n"
            msg += qc.run_mopac(s, exe=_mopac, template=_qctemplate, mult=mult,overwrite=_overwrite)
        elif _qccode == 'gaussian':
            msg += "Running gaussian...\n"
            msg += qc.run_gaussian(s, exe=_gaussian, template=_qctemplate, mult=mult,overwrite=_overwrite)
        elif _qccode == 'nwchem':
            msg += "Running nwchem...\n"
            msg += qc.run_nwchem(s, exe=_nwchem, template=_qctemplate, mult=mult,overwrite=_overwrite)
        elif _qccode == 'qcscript':
            msg += "Running qcscript...\n"
            geofile = smilesname + '.geo'
            xyzlines = ob.get_xyz(mol).splitlines()
            natom = int(xyzlines[0].strip())
            geo = ob.get_geo(mol)
            io.write_file(geo, geofile)
            if io.check_file(geofile, 1):
                msg += qc.run_qcscript(_qcscript, _qctemplate, geofile, mult)
        elif _qccode == 'submit':
            print(s)        
        print(msg)
        msg = ''
    if _writefiles:
        if io.check_file(qclog):
            newmsg, xyz,freqs,zpe,deltaH,afreqs,xmat = qc.parse_qclog_cclib(qclog, anharmonic=runanharmonic)
            msg += newmsg
            print(msg)
            io.write_file(xyz, smilesname + '.xyz')
            if freqs:
                io.write_file(qc.get_string(freqs), smilesname + '.hrm')
                if afreqs:
                    io.write_file(qc.get_string(afreqs), smilesname + '.anhrm')
                if zpe:    
                    io.write_file(str(zpe), smilesname + '.zpe')

        else:
            msg += "Log file '{0}' not found".format(qclog)
        print(msg)
        msg = ''                            
    if runthermo:
        groupstext = tc.get_new_groups()
        io.write_file(groupstext, 'new.groups')
        msg += "Parsing qc logfile '{0}'\n".format(io.get_path(qclog))
        newmsg, xyz,freqs,zpe,deltaH,afreqs,xmat = qc.parse_qclog(qclog, _qccode, anharmonic=runanharmonic)
        msg += newmsg
        if xyz is not None:
            msg += "Optimized xyz in Angstroms:\n{0} \n".format(xyz)
        else:
            runthermo = False
        if freqs is not None:
            msg += "Harmonic frequencies in cm-1:\n {0} \n".format(freqs)
        else:
            runthermo = False        
        if afreqs is not None:
            msg += "Anharmonic frequencies in cm-1:\n {0}\n".format(afreqs)
        else:
            runanharmonic = False        
        if zpe is not None:
            msg += 'ZPE = {0} kcal/mol\n'.format(zpe)
        else:
            runthermo = False        
        if deltaH is not None:
            msg += 'deltaH = {0} kcal/mol\n'.format(deltaH)
        else:
            runthermo = False        
        if xmat is not None:
            msg += 'Xmat = {0} kcal/mol\n'.format(xmat)   
        else:
            runanharmonic = False        
        if runthermo:    
            msg += write_chemkin_polynomial(mol, _qcmethod, zpe, xyz, freqs, deltaH)
    io.cd(cwd)
    print(msg)
    return


if __name__ == "__main__":
    import socket
    import iotools as io
    from timeit import default_timer as timer
    start  = timer()
    args = get_args()
    global _runqc, _runthermo, _qcmethod, _qccode,_qctemplate,_qcscript
    global _writefiles, _anharmonic,_overwrite
    global _mopac, _nwchem, _molpro, _gaussian, _messpf, _thermp, _pac99
    _runinteractive = args.runinteractive
    _runqc = args.runqc
    _runthermo = args.runthermo
    _qcmethod = args.qcmethod
    _qcbasis = args.qcbasis
    _qccode = args.qccode
    _qcpath = args.qcpath
    _anharmonic = args.anharmonic
    _overwrite = args.overwrite
    _xyzpath = args.xyzpath
    _writefiles = args.writefiles
    _qctemplate = io.get_path(args.qctemplate)
    _qcscript = io.get_path(args.qcscript,executable=True)
    _mopac = io.get_path(args.mopac,executable=True)
    _nwchem = io.get_path(args.nwchem,executable=True)
    _gaussian = io.get_path(args.gaussian,executable=True)
    _molpro = io.get_path(args.molpro,executable=True)
    _messpf = io.get_path(args.messpf,executable=True)
    _thermp = io.get_path(args.thermp,executable=True)
    _pac99 = io.get_path(args.pac99,executable=True)
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
    print("QTC: Hostname                = {0}".format(socket.gethostname()))
    print('QTC: Given arguments         =')
    for arg in vars(args):
        print('                             --{0:20s}\t{1}'.format(arg, getattr(args, arg)))
    print("QTC: Number of species       = {0}".format(len(mylist)))
    init = timer()
    print("QTC: Initialization time (s) = {0:.2f}".format(init-start))
#    print("QTC: Calculations started...")
    if nproc > 1:
        pool = multiprocessing.Pool(nproc)
        pool.map(run, mylist)
    else:
        results = map(run, mylist)
    end = timer()
    print("QTC: Calculations time (s)   = {0:.2f}".format(end - init))
    #print("QTC: Printing logs for the calculations...")
    #for result in results:
    #    print(result)
    print("QTC: Total time (s)          = {0:.2f}".format(end-start))
    print("QTC: Date and time           = {0}".format(io.get_date()))

