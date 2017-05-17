#!/usr/bin/env python
import argparse
import subprocess

import iotools as io
import obtools as ob
import qctools as qc
import tctools as tc
try:
    _runserial = False
    from scoop import futures
    from scoop import utils
except:
    _runserial = True
    print "No scoop, no concurency \n Running in serial mode..."

__updated__ = "2017-05-03"
_mopacexe = 'mopac'
_nwchemexe = 'nwchem'
_gaussianexe = 'mopac'
_messpexe = 'messpf'
_thermpexe = 'thermp'
_pac99exe = 'pac99'
_qcmethod = 'pm3'
_qccode = 'mopac'
_runqc = False
_runthermo = False


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
    parser.add_argument('-n', '--nproc', type=int,
                        default=multiprocessing.cpu_count(),
                        help='Number of processors, default is all processors')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), nargs=1,
                        default='qc_list.txt',
                        help='List of inchi or smiles for species to be calculated')
    parser.add_argument('-m', '--qcmethod', type=str, nargs=1,
                        default='pm3',
                        help='Quantum chemistry method to be used')
    parser.add_argument('-c', '--qccode', type=str, nargs=1,
                        default='mopac',
                        help='Quantum chemistry code to be used')
    parser.add_argument('-q', '--runqc', action='store_true',
                        help='Run quantum chemistry calculation')
    parser.add_argument('-t', '--runthermo', action='store_true',
                        help='Run thermochemistry calculations')
    parser.add_argument('--mopacexe', type=str, nargs=1,
                        default='mopac',
                        help='Path for mopac executable')
    parser.add_argument('--messpf', type=str, nargs=1,
                        default='messpf',
                        help='Path for MESS partition function executable')
    parser.add_argument('--thermp', type=str, nargs=1,
                        default='thermp',
                        help='Path for thermp executable')
    parser.add_argument('--pac99', type=str, nargs=1,
                        default='pac99',
                        help='Path for pac99 executable')
    return parser.parse_args()


def get_chemkin_polynomial(mol, method, zpe, xyz, freqs, deltaH):
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


def run(s):
    """
    A driver function to run quantum chemistry and thermochemistry calculations based
    on command line options:
    --qcmethod
    --qccode
    """
    import qctools as qc
    import obtools as ob
    import tctools as tc
    import iotools as io
    mol = ob.get_mol(s)
    mult = ob.get_multiplicity(mol)
    dirpath = ob.get_unique_path(mol, method=_qcmethod, mult=mult)
    groupsfile = 'new.groups'
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
        if _qccode == 'mopac':
            outstr = qc.run_mopac(s, mopacexe=_mopacexe, method=_qcmethod, mult=mult)
            outfile = outstr.split(' : ')[0]
            if _runthermo:
                lines = io.read_file(outfile, aslines=True)
                xyz = qc.get_mopac_xyz(lines)
                freqs = qc.get_mopac_freq(lines)
                zpe = qc.get_mopac_zpe(lines)
                deltaH = qc.get_mopac_deltaH(lines)
                get_chemkin_polynomial(mol, _qcmethod, zpe, xyz, freqs, deltaH)
    io.cd(cwd)
    return outstr


if __name__ == "__main__":
    args = get_args()
    print args
    global _runqc
    _runqc = args.runqc
    _runthermo = args.runthermo
    _qcmethod = args.qcmethod
    _qccode = args.qccode
    nproc = args.nproc
    mylist = io.read_list('qc_list.txt')
    results = pool.map(run, mylist)
    print 'Output file : Error code'
    for result in results:
        print result

