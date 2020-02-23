"""
Module for simplifying and enhancing the usage of Open Babel.
Open Babel is a tool-box mainly used for cheminformatics.
It enables conversions among different chemical formats,
such as inchi, smiles, xyz, zmat, etc
and extract molecular information such as force field optimized geometry,
spin, bonding information, conformers, etc.
More info:  openbabel.org
Documentation: http://openbabel.org/docs/current/index.html
"""
import pybel
import openbabel
import os
import logging
__updated__ = "2018-07-23"
__authors__ = "Murat Keceli, Sarah Elliott, Stephen Klippenstein"
__license__ = """Copyright 2017-2018 Murat Keceli, Sarah Elliott, Stephen Klippenstein"

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""


def get_periodic_table():
    """
    Return the periodic table as a list.
    Includes elements with atomic number less than 55.
    >>> pt = get_periodic_table()
    >>> print(len(pt))
    54
    """
    pt = ['X',
          'H', 'He',
          'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
          'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'
          'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
          'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe']
    return pt


def get_symbol(atomno):
    """
    Returns the element symbol for a given atomic number.
    Returns 'X' for atomno=0
    >>> print(get_symbol(1))
    H
    """
    pt = get_periodic_table()
    return pt[atomno]


def get_atomno(symbol):
    """
    Return the atomic number for a given element symbol.
    >>> print(get_atomno('H'))
    1
    """
    pt = get_periodic_table()
    symbol = symbol.capitalize()
    return pt.index(symbol)


def get_format(s):
    """
    Returns the Open Babel format of the given string.
    Note: It is primitive, distinguishes only xyz, smiles, inchi formats.
    >>> print(get_format('C'))
    smi
    >>> print(get_format('InChI=1S/H2O/h1H2'))
    inchi
    >>> print(get_format(get_xyz('C')))
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


def get_mol(s, make3D=False, mult=None):
    """
    Returns open-babel mol object from a given slabel, smiles string
    >>> mol = get_mol('[O][O]')
    >>> print(mol.formula)
    O2
    >>> print(mol.spin)
    3
    >>> mol = get_mol('InChI=1S/H2O/h1H2')
    >>> print(mol.formula)
    H2O
    >>> print(mol.spin)
    1
    """
    import pybel
    if isinstance(s, pybel.Molecule):
        mol = s
    elif isinstance(s, str) or 'unicode' in str(type(s)):
        if '_e' in s:
            s, ene = s.split('_e')
        if '_s' in s:
            s, symm = s.split('_s')
        if s.endswith('.xyz'):
            mol = next(pybel.readfile('xyz', s))
        elif '_m' in s and len(s.splitlines()) == 1:
            s, mult = s.split('_m')
            mol = set_mult(s, int(mult))
        else:
            frm = get_format(s)
            mol = pybel.readstring(frm, s)
    else:
        logging.error('Incompetible type {} for ob.get_mol'.format(type(s)))
        return None
    if make3D and mol.dim < 3:
        mol.make3D()
    if mult:
        mol.OBMol.SetTotalSpinMultiplicity(int(mult))
    return mol


def get_multiplicity(s):
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
    return get_mult(s)


def get_mult(s):
    """
    Returns spin multiplicity as an integer for a given smiles or inchi string
    """
    mult = None
    if isinstance(s, str):
        if '_m' in s:
            try:
                mult = s.split('_m')[-1]
                if '_e' in mult:
                    mult, ene = mult.split('_e')
                if '_s' in mult:
                    mult = int(mult.split('_s')[0])
            except:
                logging.debug('Multiplicity format problem, get_mult failed in s.split for s= {}'.format(s))
    if mult is None:
        mol = get_mol(s, make3D=False)
        mult = mol.spin
    return int(mult)


def get_symm(s):
    """
    Returns the symmetry factor specified in the smiles input string
    """
    sym = None
    if isinstance(s, str):
        if '_s' in s:
            try:
                sym = s.split('_s')[-1]
                if '_e' in sym:
                    sym, ene = sym.split('_e')
                sym = float(sym)
            except:
                logging.debug('Symmetry format problem, get_symm failed in s.split for s= {}'.format(s))
    return sym


def get_smileshof(s):
    """
    Returns the heat of formation at 0k specified in the smiles input string
    """
    ene = None
    if isinstance(s, str):
        if '_e' in s:
            try:
                ene = s.split('_e')[-1]
                ene = float(ene)
            except:
                logging.debug('Energy format problem, get_smileshof failed in s.split for s= {}'.format(s))
    return ene


def get_slabel(s,mult=None):
    """
    slabel is a unique smiles string for labeling species in QTC.
    Composed of two parts 'canonical smiles' and 'multiplicity'
    Canical smiles strings are unique only for a given code.
    QTC uses open babel.
    slabel = s + '_m' + str(mult)
    """
    symm = 0.
    if '_e' in s:
        s, ene = s.split('_e')
    if '_s' in s:
        s, symm = s.split('_s')
    if '_m' in s:
        s, mult = s.split('_m')
    s = get_smiles(s)
    if not mult:
        mult = get_multiplicity(s)
    return s + '_m' + str(mult)


def get_ent(s):
    """
    Return 1 if there is no enantiomer and 2 if there is based on
    @ sign appearing in the smiles
    """
    sites = 0

    xyz = get_xyz(s)
    s2  = get_smiles(xyz)
    frags = s2.split('[')
    for frag in frags:
        if '@' in frag:
            sites += 1
    ent = 2.**sites
    return ent

def get_isomers_old(s):
    """
    Return a list of stereoisomers for a given smiles, slabel, xyz or mol object.
    Currently only works for cis/trans isomerism with one double bond.
    Generated isomers include chirality but only one of them.
    Depends on xyz generated by open_babel.
    Note:
    For multiple double bonds algorithm is complicated when there is branching.
    """
    import logging
    xyz = get_xyz(s)
    s2  = get_smiles(xyz)
    ndouble  = s2.count('=')
    nslash   = s2.count('/') + s2.count('\\')
    nchiral  = s2.count('@')
    mult = None
    if '_m' in s:
        mult = get_mult(s)
        s2 = get_slabel(s2, mult)
    if s.strip() == s2.strip():
        pass
    else:
        logging.debug("SMILES changed after open babel xyz is used (get_isomers) {} --> {}".format(s, s2))
    isomers = [s2]
    if ndouble > 1:
        logging.debug('{0} double bonds in {1}'.format(ndouble, s))
        logging.debug('Can only work with one double bond')
    elif ndouble == 1 and nslash > 1:
        logging.debug('One double bond and a stereocenter found in {0}'.format(s))
        if nslash > 3:
            logging.debug('More than 3 slashes {} --> {}'.format(s, s2))
        left, right = s2.split('=')
        newright = ''
        for char in right:
            if char == '/':
                newchar = '\\'
            elif char == '\\':
                newchar = '/'
            else:
                newchar = char
            newright += newchar
        news = left + '=' + newright
        if mult:
            slabel = get_slabel(news, mult)
            news = get_slabel(news, mult)
        isomers.append(news)
    if nchiral > 0:
        logging.debug('{0} chiral centers in {1}'.format(nchiral, s))
    if len(isomers)==1:
        if '_m' in s:
            isomers = [s]
        else:
            mult = get_mult(s)
            slabel = get_slabel(s, mult)
            logging.debug("Multiplicity {} assigned by open babel for {}".format(mult, s))
            isomers = [slabel]
    return isomers


def get_isomers(s):
    """
    Return a list of cis/trans isomers.
    Could be generalized to include stereoisomers - but this has been denigrated for now.
    Depends on xyz generated by open_babel.
    For multiple double bonds algorithm only works for up to 6 slashes; i.e., 2 or 3 double bonds
    Algorithm first expands all slashes to make all combinations
    Then checks to see whether they are different
    Reduce back to original
    Coding could be made much cleaner, but Yuri is working on an x2z version that should be better anyways.
    """
    import logging
    import re
    xyz = get_xyz(s)
    s2  = get_smiles(xyz)
    nslash   = s2.count('/') + s2.count('\\')
    nchiral  = s2.count('@')
    mult = None
    if nslash == 0 and nchiral == 0:
        s2 = s
    if '_m' in s:
        mult = get_mult(s)
        s3 = get_slabel(s2, mult)
#   isomers = [s2]
    if s.strip() == s3.strip():
        pass
    else:
        logging.debug("SMILES changed after open babel xyz is used (get_isomers) {} --> {}".format(s, s2))
    segments = re.split('\\\\|/', s2)
    if nslash == 1:
       sp0 = ''
       sp0 += segments[0]
       sp0 += '/'
       sp0 += segments[1]

       sp1 = ''
       sp1 += segments[0]
       sp1 += '\\'
       sp1 += segments[1]

    elif nslash == 2:
       sp0 = ''
       sp0 += segments[0]
       sp0 += '/'
       sp0 += segments[1]
       sp0 += '/'
       sp0 += segments[2]

       sp1 = ''
       sp1 += segments[0]
       sp1 += '/'
       sp1 += segments[1]
       sp1 += '\\'
       sp1 += segments[2]

    elif nslash == 3:
       sp0 = ''
       sp0 += segments[0]
       sp0 += '/'
       sp0 += segments[1]
       sp0 += '/'
       sp0 += segments[2]
       sp0 += '/'
       sp0 += segments[3]

       sp1 = ''
       sp1 += segments[0]
       sp1 += '/'
       sp1 += segments[1]
       sp1 += '/'
       sp1 += segments[2]
       sp1 += '\\'
       sp1 += segments[3]

       sp2 = ''
       sp2 += segments[0]
       sp2 += '/'
       sp2 += segments[1]
       sp2 += '\\'
       sp2 += segments[2]
       sp2 += '/'
       sp2 += segments[3]

       sp3 = ''
       sp3 += segments[0]
       sp3 += '/'
       sp3 += segments[1]
       sp3 += '\\'
       sp3 += segments[2]
       sp3 += '\\'
       sp3 += segments[3]

    elif nslash == 4:
       sp0 = ''
       sp0 += segments[0]
       sp0 += '/'
       sp0 += segments[1]
       sp0 += '/'
       sp0 += segments[2]
       sp0 += '/'
       sp0 += segments[3]
       sp0 += '/'
       sp0 += segments[4]

       sp1 = ''
       sp1 += segments[0]
       sp1 += '/'
       sp1 += segments[1]
       sp1 += '/'
       sp1 += segments[2]
       sp1 += '/'
       sp1 += segments[3]
       sp1 += '/'
       sp1 += segments[4]

       sp2 = ''
       sp2 += segments[0]
       sp2 += '/'
       sp2 += segments[1]
       sp2 += '/'
       sp2 += segments[2]
       sp2 += '\\'
       sp2 += segments[3]
       sp2 += '/'
       sp2 += segments[4]

       sp3 = ''
       sp3 += segments[0]
       sp3 += '/'
       sp3 += segments[1]
       sp3 += '\\'
       sp3 += segments[2]
       sp3 += '/'
       sp3 += segments[3]
       sp3 += '/'
       sp3 += segments[4]

       sp4 = ''
       sp4 += segments[0]
       sp4 += '/'
       sp4 += segments[1]
       sp4 += '/'
       sp4 += segments[2]
       sp4 += '\\'
       sp4 += segments[3]
       sp4 += '\\'
       sp4 += segments[4]

       sp5 = ''
       sp5 += segments[0]
       sp5 += '/'
       sp5 += segments[1]
       sp5 += '\\'
       sp5 += segments[2]
       sp5 += '/'
       sp5 += segments[3]
       sp5 += '\\'
       sp5 += segments[4]

       sp6 = ''
       sp6 += segments[0]
       sp6 += '/'
       sp6 += segments[1]
       sp6 += '\\'
       sp6 += segments[2]
       sp6 += '\\'
       sp6 += segments[3]
       sp6 += '/'
       sp6 += segments[4]

       sp7 = ''
       sp7 += segments[0]
       sp7 += '/'
       sp7 += segments[1]
       sp7 += '\\'
       sp7 += segments[2]
       sp7 += '\\'
       sp7 += segments[3]
       sp7 += '\\'
       sp7 += segments[4]

    elif nslash == 5:
       sp0 = ''
       sp0 += segments[0]
       sp0 += '/'
       sp0 += segments[1]
       sp0 += '/'
       sp0 += segments[2]
       sp0 += '/'
       sp0 += segments[3]
       sp0 += '/'
       sp0 += segments[4]
       sp0 += '/'
       sp0 += segments[5]

       sp1 = ''
       sp1 += segments[0]
       sp1 += '/'
       sp1 += segments[1]
       sp1 += '/'
       sp1 += segments[2]
       sp1 += '/'
       sp1 += segments[3]
       sp1 += '/'
       sp1 += segments[4]
       sp1 += '\\'
       sp1 += segments[5]

       sp2 = ''
       sp2 += segments[0]
       sp2 += '/'
       sp2 += segments[1]
       sp2 += '/'
       sp2 += segments[2]
       sp2 += '/'
       sp2 += segments[3]
       sp2 += '\\'
       sp2 += segments[4]
       sp2 += '/'
       sp2 += segments[5]

       sp3 = ''
       sp3 += segments[0]
       sp3 += '/'
       sp3 += segments[1]
       sp3 += '/'
       sp3 += segments[2]
       sp3 += '\\'
       sp3 += segments[3]
       sp3 += '/'
       sp3 += segments[4]
       sp3 += '/'
       sp3 += segments[5]

       sp4 = ''
       sp4 += segments[0]
       sp4 += '/'
       sp4 += segments[1]
       sp4 += '\\'
       sp4 += segments[2]
       sp4 += '/'
       sp4 += segments[3]
       sp4 += '/'
       sp4 += segments[4]
       sp4 += '/'
       sp4 += segments[5]

       sp5 = ''
       sp5 += segments[0]
       sp5 += '/'
       sp5 += segments[1]
       sp5 += '/'
       sp5 += segments[2]
       sp5 += '/'
       sp5 += segments[3]
       sp5 += '\\'
       sp5 += segments[4]
       sp5 += '\\'
       sp5 += segments[5]

       sp6 = ''
       sp6 += segments[0]
       sp6 += '/'
       sp6 += segments[1]
       sp6 += '/'
       sp6 += segments[2]
       sp6 += '\\'
       sp6 += segments[3]
       sp6 += '/'
       sp6 += segments[4]
       sp6 += '\\'
       sp6 += segments[5]

       sp7 = ''
       sp7 += segments[0]
       sp7 += '/'
       sp7 += segments[1]
       sp7 += '\\'
       sp7 += segments[2]
       sp7 += '/'
       sp7 += segments[3]
       sp7 += '/'
       sp7 += segments[4]
       sp7 += '\\'
       sp7 += segments[5]

       sp8 = ''
       sp8 += segments[0]
       sp8 += '/'
       sp8 += segments[1]
       sp8 += '/'
       sp8 += segments[2]
       sp8 += '\\'
       sp8 += segments[3]
       sp8 += '\\'
       sp8 += segments[4]
       sp8 += '/'
       sp8 += segments[5]

       sp9 = ''
       sp9 += segments[0]
       sp9 += '/'
       sp9 += segments[1]
       sp9 += '\\'
       sp9 += segments[2]
       sp9 += '/'
       sp9 += segments[3]
       sp9 += '\\'
       sp9 += segments[4]
       sp9 += '/'
       sp9 += segments[5]

       sp10 = ''
       sp10 += segments[0]
       sp10 += '/'
       sp10 += segments[1]
       sp10 += '\\'
       sp10 += segments[2]
       sp10 += '\\'
       sp10 += segments[3]
       sp10 += '/'
       sp10 += segments[4]
       sp10 += '/'
       sp10 += segments[5]

       sp11 = ''
       sp11 += segments[0]
       sp11 += '/'
       sp11 += segments[1]
       sp11 += '/'
       sp11 += segments[2]
       sp11 += '\\'
       sp11 += segments[3]
       sp11 += '\\'
       sp11 += segments[4]
       sp11 += '\\'
       sp11 += segments[5]

       sp12 = ''
       sp12 += segments[0]
       sp12 += '/'
       sp12 += segments[1]
       sp12 += '\\'
       sp12 += segments[2]
       sp12 += '/'
       sp12 += segments[3]
       sp12 += '\\'
       sp12 += segments[4]
       sp12 += '\\'
       sp12 += segments[5]

       sp13 = ''
       sp13 += segments[0]
       sp13 += '/'
       sp13 += segments[1]
       sp13 += '\\'
       sp13 += segments[2]
       sp13 += '\\'
       sp13 += segments[3]
       sp13 += '/'
       sp13 += segments[4]
       sp13 += '\\'
       sp13 += segments[5]

       sp14 = ''
       sp14 += segments[0]
       sp14 += '/'
       sp14 += segments[1]
       sp14 += '\\'
       sp14 += segments[2]
       sp14 += '\\'
       sp14 += segments[3]
       sp14 += '\\'
       sp14 += segments[4]
       sp14 += '/'
       sp14 += segments[5]

       sp15 = ''
       sp15 += segments[0]
       sp15 += '/'
       sp15 += segments[1]
       sp15 += '\\'
       sp15 += segments[2]
       sp15 += '\\'
       sp15 += segments[3]
       sp15 += '\\'
       sp15 += segments[4]
       sp15 += '\\'
       sp15 += segments[5]

    elif nslash == 6:
       sp0 = ''
       sp0 += segments[0]
       sp0 += '/'
       sp0 += segments[1]
       sp0 += '/'
       sp0 += segments[2]
       sp0 += '/'
       sp0 += segments[3]
       sp0 += '/'
       sp0 += segments[4]
       sp0 += '/'
       sp0 += segments[5]
       sp0 += '/'
       sp0 += segments[6]

       sp1 = ''
       sp1 += segments[0]
       sp1 += '/'
       sp1 += segments[1]
       sp1 += '/'
       sp1 += segments[2]
       sp1 += '/'
       sp1 += segments[3]
       sp1 += '/'
       sp1 += segments[4]
       sp1 += '/'
       sp1 += segments[5]
       sp1 += '\\'
       sp1 += segments[6]

       sp2 = ''
       sp2 += segments[0]
       sp2 += '/'
       sp2 += segments[1]
       sp2 += '/'
       sp2 += segments[2]
       sp2 += '/'
       sp2 += segments[3]
       sp2 += '/'
       sp2 += segments[4]
       sp2 += '\\'
       sp2 += segments[5]
       sp2 += '/'
       sp2 += segments[6]

       sp3 = ''
       sp3 += segments[0]
       sp3 += '/'
       sp3 += segments[1]
       sp3 += '/'
       sp3 += segments[2]
       sp3 += '/'
       sp3 += segments[3]
       sp3 += '\\'
       sp3 += segments[4]
       sp3 += '/'
       sp3 += segments[5]
       sp3 += '/'
       sp3 += segments[6]

       sp4 = ''
       sp4 += segments[0]
       sp4 += '/'
       sp4 += segments[1]
       sp4 += '/'
       sp4 += segments[2]
       sp4 += '\\'
       sp4 += segments[3]
       sp4 += '/'
       sp4 += segments[4]
       sp4 += '/'
       sp4 += segments[5]
       sp4 += '/'
       sp4 += segments[6]

       sp5 = ''
       sp5 += segments[0]
       sp5 += '/'
       sp5 += segments[1]
       sp5 += '\\'
       sp5 += segments[2]
       sp5 += '/'
       sp5 += segments[3]
       sp5 += '/'
       sp5 += segments[4]
       sp5 += '/'
       sp5 += segments[5]
       sp5 += '/'
       sp5 += segments[6]

       sp6 = ''
       sp6 += segments[0]
       sp6 += '/'
       sp6 += segments[1]
       sp6 += '/'
       sp6 += segments[2]
       sp6 += '/'
       sp6 += segments[3]
       sp6 += '/'
       sp6 += segments[4]
       sp6 += '\\'
       sp6 += segments[5]
       sp6 += '\\'
       sp6 += segments[6]

       sp7 = ''
       sp7 += segments[0]
       sp7 += '/'
       sp7 += segments[1]
       sp7 += '/'
       sp7 += segments[2]
       sp7 += '/'
       sp7 += segments[3]
       sp7 += '\\'
       sp7 += segments[4]
       sp7 += '/'
       sp7 += segments[5]
       sp7 += '\\'
       sp7 += segments[6]

       sp8 = ''
       sp8 += segments[0]
       sp8 += '/'
       sp8 += segments[1]
       sp8 += '/'
       sp8 += segments[2]
       sp8 += '\\'
       sp8 += segments[3]
       sp8 += '/'
       sp8 += segments[4]
       sp8 += '/'
       sp8 += segments[5]
       sp8 += '\\'
       sp8 += segments[6]

       sp9 = ''
       sp9 += segments[0]
       sp9 += '/'
       sp9 += segments[1]
       sp9 += '\\'
       sp9 += segments[2]
       sp9 += '/'
       sp9 += segments[3]
       sp9 += '/'
       sp9 += segments[4]
       sp9 += '/'
       sp9 += segments[5]
       sp9 += '\\'
       sp9 += segments[6]

       sp10 = ''
       sp10 += segments[0]
       sp10 += '/'
       sp10 += segments[1]
       sp10 += '/'
       sp10 += segments[2]
       sp10 += '/'
       sp10 += segments[3]
       sp10 += '\\'
       sp10 += segments[4]
       sp10 += '\\'
       sp10 += segments[5]
       sp10 += '/'
       sp10 += segments[6]

       sp11 = ''
       sp11 += segments[0]
       sp11 += '/'
       sp11 += segments[1]
       sp11 += '/'
       sp11 += segments[2]
       sp11 += '\\'
       sp11 += segments[3]
       sp11 += '/'
       sp11 += segments[4]
       sp11 += '\\'
       sp11 += segments[5]
       sp11 += '/'
       sp11 += segments[6]

       sp12 = ''
       sp12 += segments[0]
       sp12 += '/'
       sp12 += segments[1]
       sp12 += '\\'
       sp12 += segments[2]
       sp12 += '/'
       sp12 += segments[3]
       sp12 += '/'
       sp12 += segments[4]
       sp12 += '\\'
       sp12 += segments[5]
       sp12 += '/'
       sp12 += segments[6]

       sp13 = ''
       sp13 += segments[0]
       sp13 += '/'
       sp13 += segments[1]
       sp13 += '/'
       sp13 += segments[2]
       sp13 += '\\'
       sp13 += segments[3]
       sp13 += '\\'
       sp13 += segments[4]
       sp13 += '/'
       sp13 += segments[5]
       sp13 += '/'
       sp13 += segments[6]

       sp14 = ''
       sp14 += segments[0]
       sp14 += '/'
       sp14 += segments[1]
       sp14 += '\\'
       sp14 += segments[2]
       sp14 += '/'
       sp14 += segments[3]
       sp14 += '\\'
       sp14 += segments[4]
       sp14 += '/'
       sp14 += segments[5]
       sp14 += '/'
       sp14 += segments[6]

       sp15 = ''
       sp15 += segments[0]
       sp15 += '/'
       sp15 += segments[1]
       sp15 += '\\'
       sp15 += segments[2]
       sp15 += '\\'
       sp15 += segments[3]
       sp15 += '/'
       sp15 += segments[4]
       sp15 += '/'
       sp15 += segments[5]
       sp15 += '/'
       sp15 += segments[6]

       sp16 = ''
       sp16 += segments[0]
       sp16 += '/'
       sp16 += segments[1]
       sp16 += '/'
       sp16 += segments[2]
       sp16 += '/'
       sp16 += segments[3]
       sp16 += '\\'
       sp16 += segments[4]
       sp16 += '\\'
       sp16 += segments[5]
       sp16 += '\\'
       sp16 += segments[6]

       sp17 = ''
       sp17 += segments[0]
       sp17 += '/'
       sp17 += segments[1]
       sp17 += '/'
       sp17 += segments[2]
       sp17 += '\\'
       sp17 += segments[3]
       sp17 += '/'
       sp17 += segments[4]
       sp17 += '\\'
       sp17 += segments[5]
       sp17 += '\\'
       sp17 += segments[6]

       sp18 = ''
       sp18 += segments[0]
       sp18 += '/'
       sp18 += segments[1]
       sp18 += '\\'
       sp18 += segments[2]
       sp18 += '/'
       sp18 += segments[3]
       sp18 += '/'
       sp18 += segments[4]
       sp18 += '\\'
       sp18 += segments[5]
       sp18 += '\\'
       sp18 += segments[6]

       sp19 = ''
       sp19 += segments[0]
       sp19 += '/'
       sp19 += segments[1]
       sp19 += '/'
       sp19 += segments[2]
       sp19 += '\\'
       sp19 += segments[3]
       sp19 += '\\'
       sp19 += segments[4]
       sp19 += '/'
       sp19 += segments[5]
       sp19 += '\\'
       sp19 += segments[6]

       sp20 = ''
       sp20 += segments[0]
       sp20 += '/'
       sp20 += segments[1]
       sp20 += '\\'
       sp20 += segments[2]
       sp20 += '/'
       sp20 += segments[3]
       sp20 += '\\'
       sp20 += segments[4]
       sp20 += '/'
       sp20 += segments[5]
       sp20 += '\\'
       sp20 += segments[6]

       sp21 = ''
       sp21 += segments[0]
       sp21 += '/'
       sp21 += segments[1]
       sp21 += '\\'
       sp21 += segments[2]
       sp21 += '\\'
       sp21 += segments[3]
       sp21 += '/'
       sp21 += segments[4]
       sp21 += '/'
       sp21 += segments[5]
       sp21 += '\\'
       sp21 += segments[6]

       sp22 = ''
       sp22 += segments[0]
       sp22 += '/'
       sp22 += segments[1]
       sp22 += '/'
       sp22 += segments[2]
       sp22 += '\\'
       sp22 += segments[3]
       sp22 += '\\'
       sp22 += segments[4]
       sp22 += '\\'
       sp22 += segments[5]
       sp22 += '/'
       sp22 += segments[6]

       sp23 = ''
       sp23 += segments[0]
       sp23 += '/'
       sp23 += segments[1]
       sp23 += '\\'
       sp23 += segments[2]
       sp23 += '/'
       sp23 += segments[3]
       sp23 += '\\'
       sp23 += segments[4]
       sp23 += '\\'
       sp23 += segments[5]
       sp23 += '/'
       sp23 += segments[6]

       sp24 = ''
       sp24 += segments[0]
       sp24 += '/'
       sp24 += segments[1]
       sp24 += '\\'
       sp24 += segments[2]
       sp24 += '\\'
       sp24 += segments[3]
       sp24 += '/'
       sp24 += segments[4]
       sp24 += '\\'
       sp24 += segments[5]
       sp24 += '/'
       sp24 += segments[6]

       sp25 = ''
       sp25 += segments[0]
       sp25 += '/'
       sp25 += segments[1]
       sp25 += '\\'
       sp25 += segments[2]
       sp25 += '\\'
       sp25 += segments[3]
       sp25 += '\\'
       sp25 += segments[4]
       sp25 += '/'
       sp25 += segments[5]
       sp25 += '/'
       sp25 += segments[6]

       sp26 = ''
       sp26 += segments[0]
       sp26 += '/'
       sp26 += segments[1]
       sp26 += '/'
       sp26 += segments[2]
       sp26 += '\\'
       sp26 += segments[3]
       sp26 += '\\'
       sp26 += segments[4]
       sp26 += '\\'
       sp26 += segments[5]
       sp26 += '\\'
       sp26 += segments[6]

       sp27 = ''
       sp27 += segments[0]
       sp27 += '/'
       sp27 += segments[1]
       sp27 += '\\'
       sp27 += segments[2]
       sp27 += '/'
       sp27 += segments[3]
       sp27 += '\\'
       sp27 += segments[4]
       sp27 += '\\'
       sp27 += segments[5]
       sp27 += '\\'
       sp27 += segments[6]

       sp28 = ''
       sp28 += segments[0]
       sp28 += '/'
       sp28 += segments[1]
       sp28 += '\\'
       sp28 += segments[2]
       sp28 += '\\'
       sp28 += segments[3]
       sp28 += '/'
       sp28 += segments[4]
       sp28 += '\\'
       sp28 += segments[5]
       sp28 += '\\'
       sp28 += segments[6]

       sp29 = ''
       sp29 += segments[0]
       sp29 += '/'
       sp29 += segments[1]
       sp29 += '\\'
       sp29 += segments[2]
       sp29 += '\\'
       sp29 += segments[3]
       sp29 += '\\'
       sp29 += segments[4]
       sp29 += '/'
       sp29 += segments[5]
       sp29 += '\\'
       sp29 += segments[6]

       sp30 = ''
       sp30 += segments[0]
       sp30 += '\\'
       sp30 += segments[1]
       sp30 += '\\'
       sp30 += segments[2]
       sp30 += '\\'
       sp30 += segments[3]
       sp30 += '\\'
       sp30 += segments[4]
       sp30 += '\\'
       sp30 += segments[5]
       sp30 += '/'
       sp30 += segments[6]

       sp31 = ''
       sp31 += segments[0]
       sp31 += '/'
       sp31 += segments[1]
       sp31 += '\\'
       sp31 += segments[2]
       sp31 += '\\'
       sp31 += segments[3]
       sp31 += '\\'
       sp31 += segments[4]
       sp31 += '\\'
       sp31 += segments[5]
       sp31 += '\\'
       sp31 += segments[6]
    if nslash > 6:
        logging.debug('{0} slashes in {1}'.format(nslash, s))
        logging.debug('Can only work with 6 or fewer slashes')
    if '_m' in s:
        mult = get_mult(s)
        s3 = get_slabel(s2, mult)
    if nslash  == 0:
        isomers = [s3]
    if nslash  > 0:
        xyz0 = get_xyz(sp0)
        spp0 = get_smiles(xyz0)
        if mult:
            spp0m = get_slabel(spp0, mult)
        isomers = [spp0m]
        logging.debug("spp0 {}".format(spp0m))

        xyz1 = get_xyz(sp1)
        spp1 = get_smiles(xyz1)
        if spp1 != spp0:
            if mult:
                spp1m = get_slabel(spp1, mult)
            isomers.append(spp1m)
            logging.debug("spp1 {}".format(spp1m))

    if nslash > 2:
        xyz2 = get_xyz(sp2)
        spp2 = get_smiles(xyz2)
        if spp2 != spp0 and spp2 != spp1:
            if mult:
                spp2m = get_slabel(spp2, mult)
            isomers.append(spp2m)
            logging.debug("spp2 {}".format(spp2m))

        xyz3 = get_xyz(sp3)
        spp3 = get_smiles(xyz3)
        if spp3 != spp0 and spp3 != spp1 and spp3 != spp2:
            if mult:
                spp3m = get_slabel(spp3, mult)
            isomers.append(spp3m)
            logging.debug("spp3 {}".format(spp3m))

    if nslash > 3:
        xyz4 = get_xyz(sp4)
        spp4 = get_smiles(xyz4)
        if spp4 != spp0 and spp4 != spp1 and spp4 != spp2 and spp4 != spp3:
            if mult:
                spp4m = get_slabel(spp4, mult)
            isomers.append(spp4m)
            logging.debug("spp4 {}".format(spp4m))

        xyz5 = get_xyz(sp5)
        spp5 = get_smiles(xyz5)
        if spp5 != spp0 and spp5 != spp1 and spp5 != spp2 and spp5 != spp3 and spp5 != spp4:
            if mult:
                spp5m = get_slabel(spp5, mult)
            isomers.append(spp5m)
            logging.debug("spp5 {}".format(spp5m))

        xyz6 = get_xyz(sp6)
        spp6 = get_smiles(xyz6)
        if spp6 != spp0 and spp6 != spp1 and spp6 != spp2 and spp6 != spp3 and spp6 != spp4 and spp6 != spp5:
            if mult:
                spp6m = get_slabel(spp6, mult)
            isomers.append(spp6m)
            logging.debug("spp6 {}".format(spp6m))

        xyz7 = get_xyz(sp7)
        spp7 = get_smiles(xyz7)
        if spp7 != spp0 and spp7 != spp1 and spp7 != spp2 and spp7 != spp3 and spp7 != spp4 and spp7 != spp5 and spp7 != spp6:
            if mult:
                spp7m = get_slabel(spp7, mult)
            isomers.append(spp7m)
            logging.debug("spp7 {}".format(spp7m))

    if nslash > 4:
        xyz8 = get_xyz(sp8)
        spp8 = get_smiles(xyz8)
        if spp8 != spp0 and spp8 != spp1 and spp8 != spp2 and spp8 != spp3 and spp8 != spp4 and spp8 != spp5 and spp8 != spp6 and spp8 != spp7:
            if mult:
                spp8m = get_slabel(spp8, mult)
            isomers.append(spp8m)
            logging.debug("spp8 {}".format(spp8m))

        xyz = get_xyz(sp9)
        spp9 = get_smiles(xyz9)
        if spp9 != spp0 and spp9 != spp1 and spp9 != spp2 and spp9 != spp3 and spp9 != spp4 and spp9 != spp5 and spp9 != spp6 and spp9 != spp7:
            if spp9 != spp8:
                if mult:
                    spp9m = get_slabel(spp9, mult)
                isomers.append(spp9m)
                logging.debug("spp9 {}".format(spp9m))

        xyz10 = get_xyz(sp10)
        spp10 = get_smiles(xyz10)
        if spp10 != spp0 and spp10 != spp1 and spp10 != spp2 and spp10 != spp3 and spp10 != spp4 and spp10 != spp5 and spp10 != spp6 and spp10 != spp7:
            if spp10 != spp8 and spp10 != spp9:
                if mult:
                    spp10m = get_slabel(spp10, mult)
                isomers.append(spp10m)
                logging.debug("spp10 {}".format(spp10m))

        xyz11 = get_xyz(sp11)
        spp11 = get_smiles(xyz11)
        if spp11 != spp0 and spp11 != spp1 and spp11 != spp2 and spp11 != spp3 and spp11 != spp4 and spp11 != spp5 and spp11 != spp6 and spp11 != spp7:
            if spp11 != spp8 and spp11 != spp9 and spp11 != spp10:
                if mult:
                    spp11m = get_slabel(spp11, mult)
                isomers.append(spp11m)
                logging.debug("spp11 {}".format(spp11m))

        xyz12 = get_xyz(sp12)
        spp12 = get_smiles(xyz12)
        if spp12 != spp0 and spp12 != spp1 and spp12 != spp2 and spp12 != spp3 and spp12 != spp4 and spp12 != spp5 and spp12 != spp6 and spp12 != spp7:
            if spp12 != spp8 and spp12 != spp9 and spp12 != spp10 and spp12 != spp11:
                if mult:
                    spp12m = get_slabel(spp12, mult)
                isomers.append(spp12m)
                logging.debug("spp12 {}".format(spp12m))

        xyz13 = get_xyz(sp13)
        spp13 = get_smiles(xyz13)
        if spp13 != spp0 and spp13 != spp1 and spp13 != spp2 and spp13 != spp3 and spp13 != spp4 and spp13 != spp5 and spp13 != spp6 and spp13 != spp7:
            if spp13 != spp8 and spp13 != spp9 and spp13 != spp10 and spp13 != spp11 and spp13 != spp12:
                if mult:
                    spp13m = get_slabel(spp13, mult)
                isomers.append(spp13m)
                logging.debug("spp13 {}".format(spp13m))

        xyz14 = get_xyz(sp14)
        spp14 = get_smiles(xyz14)
        if spp14 != spp0 and spp14 != spp1 and spp14 != spp2 and spp14 != spp3 and spp14 != spp4 and spp14 != spp5 and spp14 != spp6 and spp14 != spp7:
            if spp14 != spp8 and spp14 != spp9 and spp14 != spp10 and spp14 != spp11 and spp14 != spp12 and spp14 != spp13:
                if mult:
                    spp14m = get_slabel(spp14, mult)
                isomers.append(spp14m)
                logging.debug("spp14 {}".format(spp14m))

        xyz15 = get_xyz(sp15)
        spp15 = get_smiles(xyz15)
        if spp15 != spp0 and spp15 != spp1 and spp15 != spp2 and spp15 != spp3 and spp15 != spp4 and spp15 != spp5 and spp15 != spp6 and spp15 != spp7:
            if spp15 != spp8 and spp15 != spp9 and spp15 != spp10 and spp15 != spp11 and spp15 != spp12 and spp15 != spp13 and spp15 != spp14:
                if mult:
                    spp15m = get_slabel(spp15, mult)
                isomers.append(spp15m)
                logging.debug("spp15 {}".format(spp15m))

    if nslash > 5:
        xyz16 = get_xyz(sp16)
        spp16 = get_smiles(xyz16)
        if spp16 != spp0 and spp16 != spp1 and spp16 != spp2 and spp16 != spp3 and spp16 != spp4 and spp16 != spp5 and spp16 != spp6 and spp16 != spp7:
            if spp16 != spp8 and spp16 != spp9 and spp16 != spp10 and spp16 != spp11 and spp16 != spp12 and spp16 != spp13 and spp16 != spp14 and spp16 != spp15:
                if mult:
                    spp16m = get_slabel(spp16, mult)
                isomers.append(spp16m)
                logging.debug("spp16 {}".format(spp16m))

        xyz17 = get_xyz(sp17)
        spp17 = get_smiles(xyz17)
        if spp17 != spp0 and spp17 != spp1 and spp17 != spp2 and spp17 != spp3 and spp17 != spp4 and spp17 != spp5 and spp17 != spp6 and spp17 != spp7:
            if spp17 != spp8 and spp17 != spp9 and spp17 != spp10 and spp17 != spp11 and spp17 != spp12 and spp17 != spp13 and spp17 != spp14 and spp17 != spp15:
                if spp17 != spp16:
                    if mult:
                        spp17m = get_slabel(spp17, mult)
                    isomers.append(spp17m)
                    logging.debug("spp17 {}".format(spp17m))

        xyz18 = get_xyz(sp18)
        spp18 = get_smiles(xyz18)
        if spp18 != spp0 and spp18 != spp1 and spp18 != spp2 and spp18 != spp3 and spp18 != spp4 and spp18 != spp5 and spp18 != spp6 and spp18 != spp7:
            if spp18 != spp8 and spp18 != spp9 and spp18 != spp10 and spp18 != spp11 and spp18 != spp12 and spp18 != spp13 and spp18 != spp14 and spp18 != spp15:
                if spp18 != spp16 and spp18 != spp17:
                    if mult:
                        spp18m = get_slabel(spp18, mult)
                    isomers.append(spp18m)
                    logging.debug("spp18 {}".format(spp18m))

        xyz19 = get_xyz(sp19)
        spp19 = get_smiles(xyz19)
        if spp19 != spp0 and spp19 != spp1 and spp19 != spp2 and spp19 != spp3 and spp19 != spp4 and spp19 != spp5 and spp19 != spp6 and spp19 != spp7:
            if spp19 != spp8 and spp19 != spp9 and spp19 != spp10 and spp19 != spp11 and spp19 != spp12 and spp19 != spp13 and spp19 != spp14 and spp19 != spp15:
                if spp19 != spp16 and spp19 != spp17 and spp19 != spp18:
                    if mult:
                        spp19m = get_slabel(spp19, mult)
                    isomers.append(spp19m)
                    logging.debug("spp19 {}".format(spp19m))

        xyz20 = get_xyz(sp20)
        spp20 = get_smiles(xyz20)
        if spp20 != spp0 and spp20 != spp1 and spp20 != spp2 and spp20 != spp3 and spp20 != spp4 and spp20 != spp5 and spp20 != spp6 and spp20 != spp7:
            if spp20 != spp8 and spp20 != spp9 and spp20 != spp10 and spp20 != spp11 and spp20 != spp12 and spp20 != spp13 and spp20 != spp14 and spp20 != spp15:
                if spp20 != spp16 and spp20 != spp17 and spp20 != spp18 and spp20 != spp19:
                    if mult:
                        spp20m = get_slabel(spp20, mult)
                    isomers.append(spp20m)
                    logging.debug("spp20 {}".format(spp20m))

        xyz21 = get_xyz(sp21)
        spp21 = get_smiles(xyz21)
        if spp21 != spp0 and spp21 != spp1 and spp21 != spp2 and spp21 != spp3 and spp21 != spp4 and spp21 != spp5 and spp21 != spp6 and spp21 != spp7:
            if spp21 != spp8 and spp21 != spp9 and spp21 != spp10 and spp21 != spp11 and spp21 != spp12 and spp21 != spp13 and spp21 != spp14 and spp21 != spp15:
                if spp21 != spp16 and spp21 != spp17 and spp21 != spp18 and spp21 != spp19 and spp21 != spp20:
                    if mult:
                        spp21m = get_slabel(spp21, mult)
                    isomers.append(spp21m)
                    logging.debug("spp21 {}".format(spp21m))

        xyz22 = get_xyz(sp22)
        spp22 = get_smiles(xyz22)
        if spp22 != spp0 and spp22 != spp1 and spp22 != spp2 and spp22 != spp3 and spp22 != spp4 and spp22 != spp5 and spp22 != spp6 and spp22 != spp7:
            if spp22 != spp8 and spp22 != spp9 and spp22 != spp10 and spp22 != spp11 and spp22 != spp12 and spp22 != spp13 and spp22 != spp14 and spp22 != spp15:
                if spp22 != spp16 and spp22 != spp17 and spp22 != spp18 and spp22 != spp19 and spp22 != spp20 and spp22 != spp21:
                    if mult:
                        spp22m = get_slabel(spp22, mult)
                    isomers.append(spp22m)
                    logging.debug("spp22 {}".format(spp22m))

        xyz23 = get_xyz(sp23)
        spp23 = get_smiles(xyz23)
        if spp23 != spp0 and spp23 != spp1 and spp23 != spp2 and spp23 != spp3 and spp23 != spp4 and spp23 != spp5 and spp23 != spp6 and spp23 != spp7:
            if spp23 != spp8 and spp23 != spp9 and spp23 != spp10 and spp23 != spp11 and spp23 != spp12 and spp23 != spp13 and spp23 != spp14 and spp23 != spp15:
                if spp23 != spp16 and spp23 != spp17 and spp23 != spp18 and spp23 != spp19 and spp23 != spp20 and spp23 != spp21 and spp23 != spp22:
                    if mult:
                        spp23m = get_slabel(spp23, mult)
                    isomers.append(spp23m)
                    logging.debug("spp23 {}".format(spp23m))

        xyz24 = get_xyz(sp24)
        spp24 = get_smiles(xyz24)
        if spp24 != spp0 and spp24 != spp1 and spp24 != spp2 and spp24 != spp3 and spp24 != spp4 and spp24 != spp5 and spp24 != spp6 and spp24 != spp7:
            if spp24 != spp8 and spp24 != spp9 and spp24 != spp10 and spp24 != spp11 and spp24 != spp12 and spp24 != spp13 and spp24 != spp14 and spp24 != spp15:
                if spp24 != spp16 and spp24 != spp17 and spp24 != spp18 and spp24 != spp19 and spp24 != spp20 and spp24 != spp21 and spp24 != spp22 and spp24 != spp23:
                    if mult:
                        spp24m = get_slabel(spp24, mult)
                    isomers.append(spp24m)
                    logging.debug("spp24 {}".format(spp24m))

        xyz25 = get_xyz(sp25)
        spp25 = get_smiles(xyz25)
        if spp25 != spp0 and spp25 != spp1 and spp25 != spp2 and spp25 != spp3 and spp25 != spp4 and spp25 != spp5 and spp25 != spp6 and spp25 != spp7:
            if spp25 != spp8 and spp25 != spp9 and spp25 != spp10 and spp25 != spp11 and spp25 != spp12 and spp25 != spp13 and spp25 != spp14 and spp25 != spp15:
                if spp25 != spp16 and spp25 != spp17 and spp25 != spp18 and spp25 != spp19 and spp25 != spp20 and spp25 != spp21 and spp25 != spp22 and spp25 != spp23:
                    if spp25 != spp24:
                        if mult:
                            spp25m = get_slabel(spp25, mult)
                        isomers.append(spp25m)
                        logging.debug("spp25 {}".format(spp25m))

        xyz26 = get_xyz(sp26)
        spp26 = get_smiles(xyz26)
        if spp26 != spp0 and spp26 != spp1 and spp26 != spp2 and spp26 != spp3 and spp26 != spp4 and spp26 != spp5 and spp26 != spp6 and spp26 != spp7:
            if spp26 != spp8 and spp26 != spp9 and spp26 != spp10 and spp26 != spp11 and spp26 != spp12 and spp26 != spp13 and spp26 != spp14 and spp26 != spp15:
                if spp26 != spp16 and spp26 != spp17 and spp26 != spp18 and spp26 != spp19 and spp26 != spp20 and spp26 != spp21 and spp26 != spp22 and spp26 != spp23:
                    if spp26 != spp24 and spp26 != spp25:
                        if mult:
                            spp26m = get_slabel(spp26, mult)
                        isomers.append(spp26m)
                        logging.debug("spp26 {}".format(spp26m))

        xyz27 = get_xyz(sp27)
        spp27 = get_smiles(xyz27)
        if spp27 != spp0 and spp27 != spp1 and spp27 != spp2 and spp27 != spp3 and spp27 != spp4 and spp27 != spp5 and spp27 != spp6 and spp27 != spp7:
            if spp27 != spp8 and spp27 != spp9 and spp27 != spp10 and spp27 != spp11 and spp27 != spp12 and spp27 != spp13 and spp27 != spp14 and spp27 != spp15:
                if spp27 != spp16 and spp27 != spp17 and spp27 != spp18 and spp27 != spp19 and spp27 != spp20 and spp27 != spp21 and spp27 != spp22 and spp27 != spp23:
                    if spp27 != spp24 and spp27 != spp25 and spp27 != spp26:
                        if mult:
                            spp27m = get_slabel(spp27, mult)
                        isomers.append(spp27m)
                        logging.debug("spp27 {}".format(spp27m))

        xyz28 = get_xyz(sp28)
        spp28 = get_smiles(xyz28)
        if spp28 != spp0 and spp28 != spp1 and spp28 != spp2 and spp28 != spp3 and spp28 != spp4 and spp28 != spp5 and spp28 != spp6 and spp28 != spp7:
            if spp28 != spp8 and spp28 != spp9 and spp28 != spp10 and spp28 != spp11 and spp28 != spp12 and spp28 != spp13 and spp28 != spp14 and spp28 != spp15:
                if spp28 != spp16 and spp28 != spp17 and spp28 != spp18 and spp28 != spp19 and spp28 != spp20 and spp28 != spp21 and spp28 != spp22 and spp28 != spp23:
                    if spp28 != spp24 and spp28 != spp25 and spp28 != spp26 and spp28 != spp27:
                        if mult:
                            spp28m = get_slabel(spp28, mult)
                        isomers.append(spp28m)
                        logging.debug("spp28 {}".format(spp28m))

        xyz29 = get_xyz(sp29)
        spp29 = get_smiles(xyz29)
        if spp29 != spp0 and spp29 != spp1 and spp29 != spp2 and spp29 != spp3 and spp29 != spp4 and spp29 != spp5 and spp29 != spp6 and spp29 != spp7:
            if spp29 != spp8 and spp29 != spp9 and spp29 != spp10 and spp29 != spp11 and spp29 != spp12 and spp29 != spp13 and spp29 != spp14 and spp29 != spp15:
                if spp29 != spp16 and spp29 != spp17 and spp29 != spp18 and spp29 != spp19 and spp29 != spp20 and spp29 != spp21 and spp29 != spp22 and spp29 != spp23:
                    if spp29 != spp24 and spp29 != spp25 and spp29 != spp26 and spp29 != spp27 and spp29 != spp28:
                        if mult:
                            spp29m = get_slabel(spp29, mult)
                        isomers.append(spp29m)
                        logging.debug("spp29 {}".format(spp29m))

        xyz30 = get_xyz(sp30)
        spp30 = get_smiles(xyz30)
        if spp30 != spp0 and spp30 != spp1 and spp30 != spp2 and spp30 != spp3 and spp30 != spp4 and spp30 != spp5 and spp30 != spp6 and spp30 != spp7:
            if spp30 != spp8 and spp30 != spp9 and spp30 != spp10 and spp30 != spp11 and spp30 != spp12 and spp30 != spp13 and spp30 != spp14 and spp30 != spp15:
                if spp30 != spp16 and spp30 != spp17 and spp30 != spp18 and spp30 != spp19 and spp30 != spp20 and spp30 != spp21 and spp30 != spp22 and spp30 != spp23:
                    if spp30 != spp24 and spp30 != spp25 and spp30 != spp26 and spp30 != spp27 and spp30 != spp28 and spp30 != spp29:
                        if mult:
                            spp30m = get_slabel(spp30, mult)
                        isomers.append(spp30m)
                        logging.debug("spp30 {}".format(spp30m))

        xyz31 = get_xyz(sp31)
        spp31 = get_smiles(xyz31)
        if spp31 != spp0 and spp31 != spp1 and spp31 != spp2 and spp31 != spp3 and spp31 != spp4 and spp31 != spp5 and spp31 != spp6 and spp31 != spp7:
            if spp31 != spp8 and spp31 != spp9 and spp31 != spp10 and spp31 != spp11 and spp31 != spp12 and spp31 != spp13 and spp31 != spp14 and spp31 != spp15:
                if spp31 != spp16 and spp31 != spp17 and spp31 != spp18 and spp31 != spp19 and spp31 != spp20 and spp31 != spp21 and spp31 != spp22 and spp31 != spp23:
                    if spp31 != spp24 and spp31 != spp25 and spp31 != spp26 and spp31 != spp27 and spp31 != spp28 and spp31 != spp29:
                        if mult:
                            spp31m = get_slabel(spp31, mult)
                        isomers.append(spp31m)
                        logging.debug("spp31 {}".format(spp31m))

    if nchiral > 0:
        logging.debug('{0} chiral centers in {1}'.format(nchiral, s))
    if len(isomers)==1:
        if '_m' in s:
            #isomers = [s]
            isomers = [get_slabel(isomers[0])]
        else:
            mult = get_mult(s)
            slabel = get_slabel(isomer[0], mult)
            #slabel = get_slabel(s,mult)
            logging.debug("Multiplicity {} assigned by open babel for {}".format(mult, s))
            isomers = [slabel]
    return isomers


def write_isomers_list(listfile):
    """
    Writes a file containing the isomers of the species in a given list.
    The filename for the list is required as the input.
    The filename of the new file is returned.
    """
    from . import iotools as io
    slist = io.read_list(listfile)
    newlist = ''
    for s in slist:
        isomers = get_isomers(s)
        for isomer in isomers:
            newlist += '{}\n'.format(isomer)
    newfilename = listfile.split('.')[0] + '_isomers.txt'
    io.write_file(newlist, filename=newfilename)
    return newfilename


def get_formula(x, hydrogens=True, stoichiometry=True):
    """
    Returns the molecular formula.
    I think it is based on Hill system.
    Print first carbon atom stoichiometry then hydrogens,
    and all the other elements in an alphabetical order.
    Usage:
    >>> get_formula('CCC')
    'C3H8'
    >>> [get_formula(s) for s in ['C', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    ['CH4', 'C2H6O', 'O2', 'H2O']
    >>> get_formula('CCC',hydrogens=False)
    'C3'
    >>> get_formula('CCC',hydrogens=False,stoichiometry=False)
    'C'
    >>> get_formula('CCC',hydrogens=True,stoichiometry=False)
    'CH'
    """
    mol = get_mol(x)
    formula = mol.formula
    if stoichiometry:
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
    mol = get_mol(x, make3D=True)
    return len(mol.atoms)


def get_natom_heavy(x):
    """
    Return number of heavy atoms (nonHydrogens) in mol.
    >>> mol = get_mol('CC')
    >>> get_natom(mol)
    8
    >>> mols = [get_mol(s) for s in ['[CH3]', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    >>> [get_natom_heavy(mol) for mol in mols]
    [1, 3, 2, 1]
    """
    mol = get_mol(x, make3D=True)
    return mol.OBMol.NumHvyAtoms()


def get_nrotor(x):
    """
    Return number of rotors.
    """
    if isinstance(x, str):
        mol = get_mol(x, make3D=True)
    else:
        mol = get_mol(x)
    return mol.OBMol.NumRotors()


def get_nelectron(x):
    """
    Return number of electrons.
    >>> get_nelectron('C')
    10
    """
    mol = get_mol(x)
    n = 0
    for i in range(get_natom(mol)):
        a = mol.OBMol.GetAtomById(i)
        n += a.GetAtomicNum()
    return n


def get_atomic_masses(x):
    """
    Return number of electrons.
    >>> get_atomic_masses('O')
    [15.9994, 1.00794, 1.00794]
    """
    mol = get_mol(x)
    natom = get_natom(mol)
    masses = [0.0] * natom
    for i in range(natom):
        a = mol.OBMol.GetAtomById(i)
        masses[i] = a.GetAtomicMass()
    return masses


def get_atomic_numbers(x):
    """
    Return number of electrons.
    >>> get_atomic_numbers('O')
    [8, 1, 1]
    """
    mol = get_mol(x)
    natom = get_natom(mol)
    numbers = [0] * natom
    for i in range(natom):
        a = mol.OBMol.GetAtomById(i)
        numbers[i] = a.GetAtomicNum()
    return numbers


def get_xyz_dictionary(x):
    """
    Given an xyz, return a dictinary with the following keys:
    number_of_atoms
    comment
    coordinates
    coordinates_unit
    atom_symbols
    atom_numbers
    atom_masses
    atom_masses_unit
    """
    d = {}
    mol = get_mol(x, make3D=True)
    xyz = get_xyz(mol)
    lines = xyz.splitlines()
    natom = get_natom(mol)
    symbols = ['X']*natom
    coordinates = [[0.0, 0.0, 0.0]]*natom
    i = 0
    for line in lines[2:]:
        s, x, y, z = line.split()
        symbols[i] = s
        coordinates[i] = [float(x), float(y), float(z)]
        i += 1
    d['number_of_atoms'] = natom
    d['coordinates'] = coordinates
    d['coordinates_unit'] = 'Angstrom'
    d['atom_symbols'] = symbols
    d['atom_numbers'] = get_atomic_numbers(xyz)
    d['atom_masses']  = get_atomic_masses(xyz)
    d['atom_masses_unit'] = 'amu'
    return d


def get_charge(x):
    """
    Return charge.
    TODO
    """
    mol = get_mol(x, make3D=True)
    return mol.OBMol.GetTotalCharge()


def get_mass(x):
    """
    Return exact mass of mol.
    >>> mol = get_mol('CC')
    >>> print(int(get_mass(mol)))
    30
    """
    mol = get_mol(x, make3D=True)
    return mol.exactmass


def get_weight(x):
    """
    Return molecular weight  of mol.
    >>> mol = get_mol('CC')
    >>> get_weight(mol)
    30.069040000000008
    """
    mol = get_mol(x, make3D=True)
    return mol.molwt


def get_xyz(x):
    """
    Returns coordinates as a string in xyz format.
    Note: xyz coordinates are not deterministic.
    Each run gives a different set of coordinates.
    >>> mol = get_mol('CCCC')
    >>> print(get_xyz(mol).splitlines()[0])
    14
    """
    mol = get_mol(x, make3D=True)
    return mol.write(format='xyz')


def get_molpro_mol(logfile):
    """
    Returns xyz file from molpro logfile.
    """
    import pybel
    return next(pybel.readfile('mpo', logfile))


def get_gaussian_mol(logfile):
    """
    Returns mol file from gaussian logfile.
    """
    import pybel
    return next(pybel.readfile('g09', logfile))


def get_geo(x):
    """
    Returns coordinates-only as a string.
    Note: coordinates are not deterministic.
    Each run gives a different set of coordinates.
    """
    mol = get_mol(x, make3D=True)
    xyz = mol.write(format='xyz').splitlines(True)
    natom = int(xyz[0].strip())
    return ''.join(xyz[2:natom+2])


def get_zmat(x, qchem=False):
    """
    Returns internal coordinates as as string suitable for Gaussian zmat input.
    Note: zmat coordinates are deterministic.
    >>> print(get_zmat('C'))
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
    a5= 109.47
    d5= 120.00
    <BLANKLINE>
    """
    mol = get_mol(x, make3D=True)
    if qchem:
        zmat = mol.write('fh').splitlines()[2:]
        zmat[0] = zmat[0].split()[0]
        zmat = '\n'.join(zmat)
    else:
        zmat =  '\n'.join(mol.write('gzmat').splitlines()[5:])
    return zmat


def get_mop(x, keys='pm3 precise nosym threads=1 opt'):
    """
    Returns mopac input as a string.
    Note: For doctest I had to escape newline characters \n as \\n
    Since it gives EOL error.
    >>> xyz = "2\\n \\n H 0. 0. 0.\\n H 0. 0. 0.9\\n  \\n"
    >>> print(get_mop(xyz))
    pm3 precise nosym threads=1 opt
    <BLANKLINE>
    <BLANKLINE>
    H   0.00000 1  0.00000 1  0.00000 1
    H   0.00000 1  0.00000 1  0.90000 1
    <BLANKLINE>
    """
    mol = get_mol(x)
    return mol.write(format='mop', opt={'k': keys})


def get_inchi(x):
    """
    Returns a unique key composed of inchikey and multiplicity
    >>> mol = get_mol('[O][O]')
    >>> get_inchi_key(mol)
    'MYMOFIZGZYHOMD-UHFFFAOYSA-N3'
    """
    mol = get_mol(x)
    return mol.write("inchi").strip()


def get_inchi_key(x, mult=0, extra=''):
    """
    Returns a unique key composed of inchikey and multiplicity
    >>> mol = get_mol('[O][O]')
    >>> get_inchi_key(mol)
    'MYMOFIZGZYHOMD-UHFFFAOYSA-N3'
    """
    mol = get_mol(x)
    if mult == 0:
        mult = mol.spin
    return mol.write("inchikey").strip() + str(mult) + extra


def get_unique_name(x, mult=0, extra=''):
    """
    Returns a unique key composed of inchikey and multiplicity
    >>> mol = get_mol('[O][O]')
    >>> get_unique_name(mol)
    'MYMOFIZGZYHOMD-UHFFFAOYSA-N3'
    """
    mol = get_mol(x, make3D=True)
    if mult == 0:
        mult = mol.spin
    return mol.write("inchikey").strip() + str(mult) + extra


def get_unique_path(x, mult=0, method=''):
    """
    Returns a portable unique path based on inchikey for database directory.
    >>> import os
    >>> if os.path.sep == '/': print(get_unique_path('C',method='pm6'))
    database/C/C/CH4/VNWKTOKETHGBQD-UHFFFAOYSA-N1/pm6
    """
    import os
    mol = get_mol(x, make3D=True)
    if mult == 0:
        mult = mol.spin
    formula = get_formula(mol)
    formula_noH = get_formula(mol, stoichiometry=True, hydrogens=False)
    elements_noH = get_formula(mol, stoichiometry=False, hydrogens=False)
    uniquekey = get_inchi_key(mol, mult)
    dirs = 'database', elements_noH, formula_noH, formula, uniquekey, method
    return os.path.join(*dirs)


def get_formats():
    """
    Return available write formats in Open Babel.
    >>> for k, v in get_formats().items():
    ...     if 'z-matrix' in v.lower():
    ...         print('found')
    ...         break
    ...
    found
    """
    return pybel.outformats


def get_smiles_mult(slabel):
    """
    Splits slabel into smiles and multiplicity and returns them as
    a string and an integer, respectively.
    """
    smi = slabel
    mult = 0
    symm = 0.
    if '_e' in smi:
        slabel, ene = slabel.split('_e')
    if '_s' in smi:
        slabel, symm = slabel.split('_s')
    if '_m' in smi:
        smi, mult = slabel.split('_m')
    if not mult:
        mult = get_multiplicity(smi)
    return smi, int(mult)


def get_smiles_path(x, mult=0, db= 'database'):
    """
    Returns a smiles based path for database directory.
    Note that smiles strings are not unique. Even the
    canonical smiles strings are unique only for the same
    code that generates the smiles string.
    """
    from . import iotools as io
    if isinstance(x, pybel.Molecule):
        if mult == 0:
            mult = x.spin
        s = x.write(format='can').strip().split()[0]
        s = s + '_m' + str(mult)
    elif isinstance(x, str):
        if '_e' in x:
            x = x.split('_e')[0]
        if '_s' in x:
            x = x.split('_s')[0]
        if '_m' in x:
            s = x
        else:
            s = get_smiles(x)
            mult = get_mult(s)
            s = s + '_m' + str(mult)
    formula = get_formula(s)
#    formula_noH = get_formula(s, stoichiometry=True, hydrogens=False)
#    elements_noH = get_formula(s, stoichiometry=False, hydrogens=False)
    s = get_smiles_filename(s)
#    dirs = db, elements_noH, formula_noH, formula, s
    dirs = db, formula, s
    return io.join_path(*dirs)


def get_smiles(x):
    """
    Returns open-babel canonical smiles.
    >>> print(get_smiles('O'))
    O
    >>> print(get_smiles('[H][O][H]'))
    O
    >>> print(get_smiles('O-O'))
    OO
    >>> print(get_smiles('[O]=[O]'))
    O=O
    """
    mol = get_mol(x)
    s = mol.write(format='can').strip().split()[0]
    return s


def get_smiles_filename(x):
    """
    Returns a suitable filename for a given smiles.
    Smiles strings may contain characters not suitable for file names,
    such as \/:*?"<>|(). Not sure if all these characters appear, but here
    they are replaced by an underscore, '_' followed by:
    """
    if isinstance(x, pybel.Molecule):
        s = x.write(format='can').strip().split()[0]
    elif isinstance(x, str):
        if '_e' in x:
           x, ene = x.split('_e')
        if '_s' in x:
           x, sym = x.split('_s')
        s = x
    else:
        s = ''
    s = s.replace('[', '_b_')
    s = s.replace(']', '_d_')
    #s = s.replace('=','_e_')
    s = s.replace(':', '_i_')
    s = s.replace('|', '_j_')
    s = s.replace('\\', '_k_')
    s = s.replace('/', '_l_')
    s = s.replace('(', '_p_')
    s = s.replace(')', '_q_')
    s = s.replace('*', '_s_')
    s = s.replace('$', '_S_')
    s = s.replace('#', '_t_')
    s = s.replace('<', '_v_')
    s = s.replace('>', '_y_')
    s = s.replace('@@', '_aa_')
    s = s.replace('@', '_a_')
    return s


def smiles2formula(filename):
    from . import iotools as io
    mols = io.read_list(filename)
    s = ''
    for mol in mols:
        formula = get_formula(mol)
        s += '{0} {1}\n'.format(mol,  formula)
    return s


def get_coordinates_array(xyz):
    """
    Given xyz string, return natom*3 array that contains the
    coordinates of the atoms.
    """
    lines = xyz.splitlines()
    n = int(lines[0])
    coords = [0.]*(n*3)
    i = 0
    for line in lines[2:2+n]:
        coords[i:i+3] = [float(c) for c in line.split()[1:4]]
        i += 3
    return coords


def set_mult(x, mult):
    """
    Sets the total spin multiplicity.
    """
    mol = get_mol(x)
    mult = int(mult)
    assert mult > 0, 'Multiplicity should be a positve integer'
    mol.OBMol.SetTotalSpinMultiplicity(mult)
    return mol


def set_xyz(x, coords):
    """
    Parameters:
    mol : Open babel mol object, or anything that can be converted to mol object with get_mol.
    coords : One-d array of natom floats
    """
    mol = get_mol(x)
    mol.OBMol.SetCoordinates(openbabel.double_array(coords))
    return mol

try:
    import cirpy
except:
    r = 'cirpy module not installed, see http://cirpy.readthedocs.io/'
    cirpy=None
    pass
if cirpy:
    def fetch_smiles(s):
        """
        Returns the smiles string for a given chemical name.
        Requires cirpy module and internet connection
        >>> fetch_smiles('methane')
        'C'
        """
        return cirpy.resolve(s, 'smiles')


    def fetch_inchi(s):
        """
        Returns the smiles string for a given chemical name.
        Requires cirpy module and internet connection
        >>> fetch_inchi('methane')
        'InChI=1/CH4/h1H4'
        """
        return  cirpy.resolve(s, 'inchi')


    def fetch_IUPAC_name(s):
        """
        Return IUPAC name for a given smiles or inchi string.
        Requires cirpy module and internet connection
        >>> print(fetch_IUPAC_name('C=O'))
        FORMALDEHYDE
        """
        frm = get_format(s)
        if frm == 'smi':
            name = cirpy.resolve(s, 'iupac_name', resolvers=['smiles'])
        elif frm == 'inchi':
            name = cirpy.resolve(s, 'iupac_name', resolvers=['inchi'])
        elif frm == 'xyz':
            mol = get_mol(s)
            name = cirpy.resolve(mol.write('inchi').strip(), 'iupac_name', resolvers=['inchi'])
        else:
            name = None
        return name

def get_ase_atoms_input(x):
    """
    Return atomic symbols and coordinates that can be used as inputs
    for ASE Atoms constructor.
    >>> symbols, coords = get_ase_atoms_input('C')
    >>> 'C' in symbols
    True
    >>> print(len(symbols))
    5
    >>> print(len(coords[0]))
    3
    >>> print(len(coords))
    5
    """
    xyz = get_xyz(x)
    lines = xyz.splitlines()
    n = int(lines[0])
    symbols = ['X']*n
    coords = [[0., 0., 0.]]*n
    for i, line in enumerate(lines[2:2+n]):
        items = line.split()
        symbols[i] = items[0]
        coords[i] = [float(c) for c in items[1:4]]
    return symbols, coords

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)


