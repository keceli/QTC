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
    parser.add_argument('-T', '--runthermo', action='store_true',
                        help='Run thermochemistry calculations')
    parser.add_argument('-I', '--runinteractive', action='store_true',
                        help='Interactive mode for QTC')
    parser.add_argument('--mopacpath', type=str,
                        default='mopac',
                        help='Path for mopac executable')
    parser.add_argument('--nwchempath', type=str,
                        default='nwchem',
                        help='Path for nwchem executable')
    parser.add_argument('--gaussianpath', type=str,
                        default='g09',
                        help='Path for gaussian executable')
    parser.add_argument('--messpfpath', type=str,
                        default='messpf',
                        help='Path for MESS partition function executable')
    parser.add_argument('--thermppath', type=str,
                        default='thermp',
                        help='Path for thermp executable')
    parser.add_argument('--pac99path', type=str,
                        default='pac99',
                        help='Path for pac99 executable')
    parser.add_argument('--qcscriptpath', type=str,
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

    msg = 'Running mess partition function.\n'
    tc.run_pf()
    msg += 'Generating thermp input.\n'
    tc.write_thermp_input(mol.formula, deltaH)
    msg += 'Running thermp.\n'
    tc.run_thermp()
    msg += 'Running pac99.\n'
    tc.run_pac99(name)
    msg += 'Converting to chemkin format.\n'
    chemkinfile = name + '.ckin'
    msg += 'Writing chemking file {0}.\n'.format(chemkinfile)
    tc.write_chemkin_file(deltaH, tag, name, chemkinfile)
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
    msg = "{0}\n".format(s)
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
                msg += ('Could not copy new.groups file to target directory {0}.\n'.format(dirpath))
                return -1
        else:
            msg += ('new.groups file required in working directory.\n')
            return -1
    if io.check_dir(dirpath, 1):
        io.cd(dirpath)
    else:
        msg += ('I/O error, {0} directory not found.\n'.format(dirpath))
        return -1
    if _runqc:
        if io.check_file(runfile):
            msg += ('Skipping {0}.\n'.format(uniquekey))
        else:
            io.touch(runfile)
        if _qctype == 'mopac':
            msg += "Running mopac...\n"
            msg += qc.run_mopac(s, exe=_mopacpath, method=_qcmethod, mult=mult)
            outfile = msg.split(' : ')[0]
        elif _qctype == 'gaussian':
            msg += "Running gaussian...\n"
            msg += qc.run_gaussian(s, exe=_gaussianpath, template=_qctemplate, mult=mult,overwrite=False)
            outfile = msg.split(' : ')[0]                
        elif _qctype == 'qcscript':
            msg += "Running qcscript...\n"
            geofile = uniquekey + '.geo'
            xyzlines = ob.get_xyz(mol).splitlines()
            geo = ''.join(xyzlines[2:])
            io.write_file(geo, geofile)
            if io.check_file(geofile, 1):
                msg += qc.run_qcscript(_qcscriptpath, _qctemplate, geofile, mult)
        elif _qctype == 'submit':
            print(s)        
                
    if _runthermo:
        msg += "Parsing qc output...\n"
        lines = io.read_file(outfile, aslines=True)
        xyz = qc.get_mopac_xyz(lines)
        freqs = qc.get_mopac_freq(lines)
        zpe = qc.get_mopac_zpe(lines)
        deltaH = qc.get_mopac_deltaH(lines)
        msg += write_chemkin_polynomial(mol, _qcmethod, zpe, xyz, freqs, deltaH)
    io.cd(cwd)
    return msg


if __name__ == "__main__":
    import socket
    args = get_args()
    global _runqc, _runthermo, _qcmethod, _qctype,_qctemplate,_qcscriptpath
    global _mopacpath, _nwchempath, _gaussianpath, _messpfpath, _thermppath, _pac99path
    _runqc = args.runqc
    _runthermo = args.runthermo
    _qcmethod = args.qcmethod
    _qctype = args.qctype
    _qctemplate = args.qctemplate
    _qcscriptpath = args.qcscriptpath
    _mopacpath = args.mopacpath
    _nwchempath = args.nwchempath
    _gaussianpath = args.gaussianpath
    _messpfpath = args.messpfpath
    _thermppath = args.thermppath
    _pac99path = args.pac99path
    _runinteractive = args.runinteractive
    inp = args.input
    nproc = args.nproc
    mylist = io.read_list(inp)
    print('Given arguments')
    print(args)
    print("Input file {0} has {1} species".format(inp,len(mylist)))
    print("QTC running with {0} processes on host {1}.\n".format(nproc, socket.gethostname()))
    pool = multiprocessing.Pool(nproc)
    results = pool.map(run, mylist)

    for result in results:
        print result

