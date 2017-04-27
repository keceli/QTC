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
__updated__ = "2017-04-27"

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
    if s.startswith('InChI'):
        mol = pybel.readstring("inchi", s)
    else:
        mol = pybel.readstring("smi", s)
    if make3D:
        mol.make3D()
    return mol

def get_spin(mol):
    """
    Returns the spin of the molecule.
    Usage:
    >>> mol = get_mol('C')
    >>> get_spin(mol)
    1
    >>> mols = [get_mol(s) for s in ['[CH3]', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    >>> [get_spin(mol) for mol in mols]
    [2, 1, 3, 1]
    """
    return mol.spin


def get_formula(mol):
    """
    Returns the molecular formula.
    I think it is based on Hill system.
    Print first carbon atom stoichemetry then hydrogens,
    and all the other elements in an alphabetical order.
    Usage:
    >>> mol = get_mol('CCC')
    >>> get_formula(mol)
    'C3H8'
    >>> mols = [get_mol(s) for s in ['C', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    >>> [get_formula(mol) for mol in mols]
    ['CH4', 'C2H6O', 'O2', 'H2O']
    """
    return mol.formula


def get_natom(mol):
    """
    Return number of atoms in mol.
    >>> mol = get_mol('CC')
    >>> get_natom(mol)
    8
    >>> mols = [get_mol(s) for s in ['[CH3]', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    >>> [get_natom(mol) for mol in mols]
    [4, 9, 2, 3]
    """
    return len(mol.atoms)


def get_xyz(mol):
    """
    Returns coordinates in xyz format.
    Note: xyz coordinates are not deterministic.
    Each run gives a different set of coordinates.
    >>> mol = get_mol('CCCC')
    >>> print get_xyz(mol).splitlines()[0]
    14
    """
    return mol.write(format='xyz')


def get_uniquekey(mol):
    """
    Returns a unique key composed of inchikey and multiplicity
    >>> mol = get_mol('[O][O]')
    >>> get_uniquekey(mol)
    'MYMOFIZGZYHOMD-UHFFFAOYSA-N3'
    """
    return mol.write("inchikey").strip() + str(mol.spin)


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)


