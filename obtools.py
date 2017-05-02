#!/usr/bin/env python
import pybel
import openbabel
import os

"""
Module for simplifying and enhancing the usage of Open Babel.
Open Babel is a tool-box mainly used for cheminformatics.
It allows us to make conversions among different chemical formats,
such as inchi, smiles, xyz, zmat, etc
and extract molecular information such as force field optimized geometry,
spin, bonding information, conformers, etc.
More info on:  openbabel.org
Documentation on: http://openbabel.org/docs/current/index.html
This module is useful for a new user of Open Babel since it
provides information on the functionalities and how to use them
in python.
"""
__updated__ = "2017-05-02"

def get_format(s):
    """
    Returns the Open Babel format of the given string.
    Note: It is primitive, distinguishes only xyz, smiles, inchi formats.
    >>> print get_format('C')
    smi
    >>> print get_format('InChI=1S/H2O/h1H2')
    inchi
    >>> print get_format(get_xyz('C'))
    xyz
    """
    frm = 'unknown'
    lines = s.splitlines()
    n = len(lines)  # number of lines
    if n == 1:
        if s.startswith('InChI'):
            frm = 'inchi'
        else:
            frm = 'smi'
    else:
        try:
            natom = int(lines[0].strip())
            if n >= natom + 2:
                frm = 'xyz'
        except:
            pass

    return frm


def get_mol(s, make3D=True):
    """
    Returns open-babel mol object from a given inchi or smiles string
    >>> mol = get_mol('[O][O]')
    >>> print mol.formula
    O2
    >>> print mol.spin
    3
    >>> mol = get_mol('InChI=1S/H2O/h1H2')
    >>> print mol.formula
    H2O
    >>> print mol.spin
    1
    """
    import pybel
    frm = get_format(s)
    mol = pybel.readstring(frm, s)
    if make3D and not frm == 'xyz':
        mol.make3D()
    return mol

def get_multiplicity(x):
    """
    Returns the spin multiplicity (2S+1) of the molecule, where S is the
    total spin angular momentum.
    Usage:
    >>> mol = get_mol('C')
    >>> get_multiplicity(mol)
    1
    >>> mols = [get_mol(s) for s in ['[CH3]', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    >>> [get_multiplicity(mol) for mol in mols]
    [2, 1, 3, 1]
    """
    if type(x) is str:
        mol = get_mol(x)
    else:
        mol = x
    return mol.spin


def get_formula(x, hydrogens=True, stoichemetry=True):
    """
    Returns the molecular formula.
    I think it is based on Hill system.
    Print first carbon atom stoichemetry then hydrogens,
    and all the other elements in an alphabetical order.
    Usage:
    >>> get_formula('CCC')
    'C3H8'
    >>> [get_formula(s) for s in ['C', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    ['CH4', 'C2H6O', 'O2', 'H2O']
    >>> get_formula('CCC',hydrogens=False)
    'C3'
    >>> get_formula('CCC',hydrogens=False,stoichemetry=False)
    'C'
    >>> get_formula('CCC',hydrogens=True,stoichemetry=False)
    'CH'
    """
    if type(x) is str:
        mol = get_mol(x)
    else:
        mol = x

    formula = mol.formula
    if stoichemetry:
        if hydrogens:
            s = formula
        else:
            n = len(formula)
            s = ''
            hdigit = False
            for i in range(n):
                if formula[i].isdigit():
                    if not hdigit:
                        s += formula[i]
                elif formula[i] == 'H':
                    hdigit = True
                else:
                    hdigit = False
                    s += formula[i]
            if s == '':
                s = formula
    else:
        if hydrogens:
            s = ''.join(i for i in formula if not i.isdigit())
        else:
            s = ''.join(i for i in formula if not i.isdigit() and not i == 'H')
            if s == '':
                s = 'H'
    return s


def get_natom(x):
    """
    Return number of atoms in mol.
    >>> mol = get_mol('CC')
    >>> get_natom(mol)
    8
    >>> mols = [get_mol(s) for s in ['[CH3]', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    >>> [get_natom(mol) for mol in mols]
    [4, 9, 2, 3]
    """
    if type(x) is str:
        mol = get_mol(x)
    else:
        mol = x
    return len(mol.atoms)


def get_xyz(x):
    """
    Returns coordinates as a string in xyz format.
    Note: xyz coordinates are not deterministic.
    Each run gives a different set of coordinates.
    >>> mol = get_mol('CCCC')
    >>> print get_xyz(mol).splitlines()[0]
    14
    """
    if type(x) is str:
        mol = get_mol(x)
    else:
        mol = x
    return mol.write(format='xyz')


def get_zmat(x):
    """
    Returns internal coordinates as as string suitable for Gaussian zmat input.
    Note: zmat coordinates are deterministic.
    >>> print get_zmat('C')
    C
    H  1  r2
    H  1  r3  2  a3
    H  1  r4  2  a4  3  d4
    H  1  r5  2  a5  3  d5
    Variables:
    r2= 1.0922
    r3= 1.0922
    a3= 109.47
    r4= 1.0922
    a4= 109.47
    d4= 240.00
    r5= 1.0922
    """
    if type(x) is str:
        mol = get_mol(x)
    else:
        mol = x
    return '\n'.join(mol.write('gzmat').splitlines()[5:-3])


def get_mop(x, keys='pm3 precise nosym threads=1 opt'):
    """
    Returns mopac input as a string.
    Note: For doctest I had to escape newline characters \n as \\n
    Since it gives EOL error.
    >>> xyz = "2\\n \\n H 0. 0. 0.\\n H 0. 0. 0.9\\n  \\n"
    >>> print get_mop(xyz)
    pm3 precise nosym threads=1 opt
    <BLANKLINE>
    <BLANKLINE>
    H   0.00000 1  0.00000 1  0.00000 1
    H   0.00000 1  0.00000 1  0.90000 1
    <BLANKLINE>
    """
    if type(x) is str:
        mol = get_mol(x)
    else:
        mol = x
    return mol.write(format='mop', opt={'k': keys})


def get_unique_key(x, mult=0, extra=''):
    """
    Returns a unique key composed of inchikey and multiplicity
    >>> mol = get_mol('[O][O]')
    >>> get_unique_key(mol)
    'MYMOFIZGZYHOMD-UHFFFAOYSA-N3'
    """
    if type(x) is str:
        mol = get_mol(x)
    else:
        mol = x
    if mult == 0:
        mult = mol.spin
    return mol.write("inchikey").strip() + str(mult) + extra


def get_unique_path(x, method='', mult=0):
    """
    Returns a portable unique path for database directory.
    >>> import os
    >>> if os.path.sep == '/': print get_unique_path('C',method='pm6')
    database/C/C/CH4/VNWKTOKETHGBQD-UHFFFAOYSA-N1/pm6
    """
    import iotools as io
    if type(x) is str:
        mol = get_mol(x)
    else:
        mol = x
    if mult == 0:
        mult = mol.spin
    formula = get_formula(mol)
    formula_noH = get_formula(mol, stoichemetry=True, hydrogens=False)
    elements_noH = get_formula(mol, stoichemetry=False, hydrogens=False)
    uniquekey = get_unique_key(mol, mult)
    dirs = 'database', elements_noH, formula_noH, formula, uniquekey, method
    return io.join_path(*dirs)


def get_formats():
    """
    Return available write formats in Open Babel.
    >>> for k, v in get_formats().items():
    ...     if 'z-matrix' in v.lower():
    ...         print v, k
    ...
    ...
    Gaussian Z-Matrix Input gzmat
    Fenske-Hall Z-Matrix format fh
    """
    return pybel.outformats


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)


