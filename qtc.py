#!/usr/bin/env python
import argparse
import subprocess
import multiprocessing

import iotools as io
import obtools as ob
import qctools as qc
import tctools as tc

__updated__ = "2017-05-15"
__author__ = "Murat Keceli"
# _mopacexe = 'mopac'
# _nwchemexe = 'nwchem'
# _gaussianexe = 'g09'
# _messpfexe = 'messpf'
# _thermpexe = 'thermp'
# _pac99exe = 'pac99'
# _qcmethod = 'pm3'
# _qctype = 'mopac'
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
    Integrates different codes for these purposes.
    """)
    parser.add_argument('-i', '--input', type=str,
                        default='qclist.txt',
                        help='List of inchi or smiles for species to be calculated')
    parser.add_argument('-n', '--nproc', type=int,
                        default=multiprocessing.cpu_count(),
                        help='Number of processors, default is all processors')
    parser.add_argument('-m', '--qcmethod', type=str,
                        default='pm3',
                        help='Quantum chemistry method to be used')
    parser.add_argument('-t', '--qctype', type=str,
                        default='mopac',
                        help='Quantum chemistry calculation type ("gausian","mopac","qcscript") to be used')
    parser.add_argument('-p', '--qctemplate', type=str,
                        default='qctemplate.txt',
                        help='Template for quantum chemistry input file')
    parser.add_argument('-Q', '--runqc', action='store_true',
                        help='Run quantum chemistry calculation')
    parser.add_argument('-S', '--subqc', action='store_true',
                        help='Submit quantum chemistry calculations')
    parser.add_argument('-J', '--qcscript', action='store_true',
                        help='Submit quantum chemistry calculations using qcscript.pl')
    parser.add_argument('-T', '--runthermo', action='store_true',
                        help='Run thermochemistry calculations')
    parser.add_argument('-I', '--runinteractive', action='store_true',
                        help='Interactive mode for QTC')
    parser.add_argument('--qcscript', type=str,
                        default='qcscript.pl',
                        help='Path for qcscript')
    parser.add_argument('--mopacexe', type=str,
                        default='mopac',
                        help='Path for mopac executable')
    parser.add_argument('--nwchemexe', type=str,
                        default='nwchem',
                        help='Path for nwchem executable')
    parser.add_argument('--gaussianexe', type=str,
                        default='g09',
                        help='Path for gaussian executable')
    parser.add_argument('--messpfexe', type=str,
                        default='messpf',
                        help='Path for MESS partition function executable')
    parser.add_argument('--thermpexe', type=str,
                        default='thermp',
                        help='Path for thermp executable')
    parser.add_argument('--pac99exe', type=str,
                        default='pac99',
                        help='Path for pac99 executable')
    parser.add_argument('--qcscriptexe', type=str,
                        default='/lcrc/project/PACC/test-awj/builddb/bin/qcscript.pl',
                        help='Path for qcscript perl script')
    return parser.parse_args()


def write_chemkin_polynomial(mol, method, zpe, xyz, freqs, deltaH):
    """
    A driver to perform all operations to write NASA polynomial in
    chemkin format. Assumes quantum chemistry calculation is performed.
    """
    inputfile = 'pf.inp'
    name = mol.formula
    tag = method
    inp = tc.get_pf_input(mol, method, zpe, xyz, freqs)

 #   print 'Running mess partition function'
    tc.run_pf()
 #   print 'Generate thermp input'
    tc.write_thermp_input(mol.formula, deltaH)
 #   print 'Running thermp'
    tc.run_thermp()
#    print 'Running pac99'
    tc.run_pac99(name)
#    print 'Converting to chemkin format'
    chemkinfile = name + '.ckin'
    tc.write_chemkin_file(deltaH, tag, name, chemkinfile)
    return


def subqc(s,methodfile, mult=0):
    """
    TODO
    """
    pass
    return

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
    mol = ob.get_mol(s)
    mult = ob.get_multiplicity(mol)
    uniquekey = ob.get_unique_key(mol,mult,extra=_qcmethod)
    dirpath = ob.get_unique_path(mol, method=_qcmethod, mult=mult)
    groupsfile = 'new.groups'
    runfile = uniquekey + '.run'
    io.mkdir(dirpath)
    cwd = io.pwd()
    if _runthermo:
        if io.check_file(groupsfile):
            io.cp(groupsfile, dirpath)
            if not io.check_file(groupsfile, 1):
                print 'Could not copy new.groups file to target directory {0}'.format(dirpath)
                return -1
        else:
            print 'new.groups file required in working directory'
            return -1
    if io.check_dir(dirpath, 1):
        io.cd(dirpath)
    else:
        print 'I/O error, {0} directory not found'.format(dirpath)
        return -1
    if _runqc:
        if io.check_file(runfile):
            print('Skipping {0}, it is already running'.format(uniquekey))
        else:
            io.touch(runfile)
        if _qctype == 'mopac':
            outstr = qc.run_mopac(s, exe=_mopacexe, method=_qcmethod, mult=mult)
            outfile = outstr.split(' : ')[0]
        elif _qctype == 'gaussian':
            outstr = qc.run_gaussian(s, exe=_gaussianexe, template=_qctemplate, mult=mult,overwrite=False)
            outfile = outstr.split(' : ')[0]                
        elif _qctype == 'qcscript':
            geofile = uniquekey + '.geo'
            xyzlines = ob.get_xyz(mol).splitlines()
            geo = ''.join(xyzlines[2:])
            io.write_file(geo, geofile)
            if io.check_file(geofile, 1):
                qc.run_qcscript(_qcscript, _qctemplate, geofile, mult)
        elif _qctype == 'submit':
            print(s)        
                
    if _runthermo:
        lines = io.read_file(outfile, aslines=True)
        xyz = qc.get_mopac_xyz(lines)
        freqs = qc.get_mopac_freq(lines)
        zpe = qc.get_mopac_zpe(lines)
        deltaH = qc.get_mopac_deltaH(lines)
        write_chemkin_polynomial(mol, _qcmethod, zpe, xyz, freqs, deltaH)
    io.cd(cwd)
    return outstr


if __name__ == "__main__":
    args = get_args()
    global _runqc, _runthermo, _qcmethod, _qctype,_qctemplate,_qcscript
    global _mopacexe, _nwchemexe, _gaussianexe, _messpfexe, _thermpexe, _pac99exe
    _runqc = args.runqc
    _runthermo = args.runthermo
    _qcmethod = args.qcmethod
    _qctype = args.qccode
    _qctemplate = args.qctemplate
    _qcscript = args.qcscript
    _mopacexe = args.mopacexe
    _nwchemexe = args.nwchemexe
    _gaussianexe = args.gaussianexe
    _messpfexe = args.messpfexe
    _thermpexe = args.thermpexe
    _pac99exe = args.pac99exe
    _runinteractive = args.runinteractive
    inp = args.input
    nproc = args.nproc
    mylist = io.read_list(inp)
    print("Input file {0} has {1} species".format(inp,len(mylist)))
    pool = multiprocessing.Pool(nproc)
    results = pool.map(run, mylist)
    print('Given arguments')
    for arg in args:
        print(arg)
    print 'Output file : Error code'
    for result in results:
        print result

