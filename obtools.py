#!/usr/bin/env python
import pybel
import os

def get_mol(s, make3D=True):
    """
    Returns open-babel mol object from a given inchi or smiles string
    """
    import pybel
    if s.startswith('InChI'):
        mol=pybel.readstring("inchi",s)
    else:
        mol=pybel.readstring("smi",s)
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
    >>> mol = get_mol('[O]')
    >>> get_xyz(mol)
    '1\n\nO          0.97366       -0.05852        0.07669\n'
    """
    return mol.write(format='xyz')


def get_uniquekey(mol):
    """
    Returns a unique key composed of inchikey and multiplicity
    """
    return mol.write("inchikey").strip()+str(mol.spin)


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)


