#!usr/bin/python

import re
import numpy as np
import os
import sys
#sys.path.insert(0, '/home/elliott/Packages/QTC/')
sys.path.insert(0, '/home/snelliott/projects/anl/QTC/')
import iotools as io
import qctools as qc
import patools as pa
import obtools as ob
import logging
import unittools as ut
"""
Heatform determines the heat of formation for a molecule by using a basis of 
molecules with well-determined heats of formation

To use this script from command line you can do the following:

(1) python heatform.py -s <smiles> -o prog/method/basis -f prog/method/basis -e prog/method/basis -r
    will run (-r) an optimization (-o) on the smiles molecule, will use that geometry for 
    an energy (-e) computation, and obtain a -f level zpve and do the same for the basis molecules.
    ONLY molpro and gaussian are available

(2) python heatform.py -l <logfile.log>
    this will automatically determine what stoichiometry, program, level of theory, basis set,
    and electronic energy heatform needs as parameters.  Will assume opt, freq, and en levels
    are all the same

(3) python heatform.py -s <smiles> -o prog/method/basis -f prog/method/basis -e prog/method/basis -E <energy>
    manually set all these parameters so that no logfile OR computation on main species is required
    if you specify a freqlevel, add zpve to your -E.  

(4) User can also use -B <R1 R2 R3> to manually select the basis of molecules or -b if they don't
    use smiles and -a will give anharmonic

To use this as a module in another script:

(1) main(stoich, logfile, electronic_energy, optlevel, freqlevel, enlevel, anharm, runE) will return 
    the heat of formation in kcal.  If logfile is set, the rest can be auto (see defaults), and if 
    logfile is gibberish, the rest must be set.
(2) main_keyword(parameters) for use from QTC

"""
#  TO DO:
#  Cleanup
#  Select basis based on type of bonds

def get_atomlist(mol):
    """
    Makes a list of all atoms in a molecule
    INPUT:
    mol      - stoichiometry of molecule
    OUTPUT:
    atomlist - list of distinct atoms in that molecule
    """
    atomlist = []
    elements = {'He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar',
      'Ca','Sc','Ti','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge', 
      'As','Se','Br','Kr','C','B','H','O','F','S','N','P','K','V'}
    for el in elements:
        if el in mol:
           atomlist.append(el)
           mol = mol.replace(el,'')
    return atomlist

def select_basis(atomlist,attempt=0):
    """
    Given a list of atoms, generates a list of molecules
    that is best suited to serve as a basis for those atoms
    INPUT:
    atomlist - list of atoms
    OUPUT:
    basis    - recommended basis as a list of stoichiometries
    """
    count = len(atomlist)-1
    basis = []
    i= 0
    if 'N' in atomlist and i <= count and attempt < 2:
        basis.append('N#N')
        i += 1
    if 'N' in atomlist and 'H' in atomlist and i <= count and attempt >1:
        basis.append('N')
        i += 1
    if 'S' in atomlist and i<= count:
        basis.append('O=S=O') 
        i += 1
    if 'H' in atomlist and i<= count and attempt < 2:
        basis.append('[H][H]')
        i += 1
    elif 'H' in atomlist and 'C' not in atomlist and i<= count and attempt < 3:
        basis.append('[H][H]')
        i += 1
    if 'O' in atomlist and i<= count and attempt < 3:
        basis.append('[O][O]')
        i += 1
    if 'C' in atomlist and i<= count and attempt < 4:
        basis.append('C')
        i += 1
    if 'O' in atomlist and 'H' in atomlist and  i<= count and attempt  < 4:
        basis.append('O')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count and attempt < 5:
        basis.append('C(=O)=O')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count and attempt < 5:
        basis.append('C=O')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count:
        basis.append('CO')
        i += 1
    if 'C' in atomlist and i<= count:
        basis.append('CC')
        i += 1
    if 'S' in atomlist and i<= count:
        basis.append('O=S=O') 
        i += 1
    if 'H' in atomlist  and i<= count and attempt < 1:
        basis.append('[H][H]')
        i += 1
    elif 'H' in atomlist and 'C' not in atomlist and i<= count and attempt < 3:
        basis.append('[H][H]')
        i += 1
    if 'O' in atomlist and i<= count and attempt < 2:
        basis.append('[O][O]')
        i += 1
    if 'C' in atomlist and i<= count and attempt < 3:
        basis.append('C')
        i += 1
    if 'O' in atomlist and 'H' in atomlist and  i<= count and attempt  < 3:
        basis.append('O')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count and attempt < 4:
        basis.append('C(=O)=O')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count and attempt < 4:
        basis.append('C=O')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count:
        basis.append('CO')
        i += 1
    if 'C' in atomlist and i<= count:
        basis.append('CC')
        i += 1
    return basis

def get_stoich(mol,atomlist):
    """
    Given a molecule's stoichiometry and a list of atoms, finds the
    number of each atom that the molecule contains
    INPUT:
    mol       - molecule stoichiometry
    atomlist  - list of atoms 

    OUTPUT:
    stoich    - list of numbers corresponding to the number of each atom 
                in the atomlist that the molecule contains
    """
    stoichlist = np.zeros(len(atomlist))

    for i, atom in enumerate(atomlist):
        val = 0
        if atom in mol:
            a = re.compile(atom + '(\d*)')
            a = a.findall(mol)
            for b in a:
                if b == '':
                    b = '1'
                val += float(b)
        stoichlist[i] = val

    return stoichlist


def form_mat(basis,atomlist):
    """
    Form a matrix for a given basis and atomlist
    INPUT:  
    basis     - basis of molecules
    atomlist  - list of atoms (all atoms that appear 
                in basis should be in atomlist)
    OUTPUT:
    mat       - matrix (length of basis by length of atomlist)
                (square if done right)
    """
    mat = np.zeros((len(atomlist),len(atomlist)))
    for i,mol in enumerate(basis):
        mat[i] = get_stoich(ob.get_formula(ob.get_mol(mol)),atomlist)
    mat = mat.T
     
    return mat

def comp_coeff(mat,stoich):
    """
    Finds the coefficients that solve C = M^-1 S.  For our purposes C are the coefficients [a,b,c,d..] for 
    the basis [R1, R2, R3, R4...] that give CxHyOzNw =  aR1 + bR2 + cR3 + dR4 where [x,y,z,w...] is S, our 
    stoichiometry vector for a molecule.  M is a nonsingular matrix that puts the basis in terms of 
    a list of atoms
    INPUT:
    mat    - nonsingular matrix 
    stoich - list of numbers that give a molecules stoichiometry in terms of a list of atoms

    OUTPUT:
    coeff  - coefficients [a,b,c,d...] as described above
    """
    mati = np.linalg.inv(mat)
    coeff = np.dot(mati,stoich)

    return coeff

#def select_better_basis(smiles, atomlist):    
##test script to pick basis based on bond type:
#
#    count = len(atomlist)-1
#    doublebond =[]
#    triplebond =[]
#    singlebond =[]
#    branch = '' 
#    
#    if '(' in smiles:
#        branches = smiles.split('(')
#        branch = branches[0][-1] + branches[1].split(')')[0]
#        smiles = branches[0] + branches[1].split(')')[1]
#    if '=' in branch:
#        frag = branch.split('=')
#        for i in range(len(frag)-1):
#            doublebond.append(frag[i][-1] + frag[i+1][0])
#    if '#' in branch:
#        frag = branch.split('#')
#        for j in range(len(frag)-1):
#             triplebond.append(frag[j][-1] + frag[j+1][0])
#    if '=' in smiles:
#        frag = smiles.split('=')
#        for i in range(len(frag)-1):
#            doublebond.append(frag[i][-1] + frag[i+1][0])
#    if '#' in smiles:
#        frag = smiles.split('#')
#        for j in range(len(frag)-1):
#             triplebond.append(frag[j][-1] + frag[j+1][0])
#    
#    if 'CC' in smiles:
#        singlebond.append('CC')
#    if 'CO' in smiles or 'OC' in smiles:
#        singlebond.append('CO')
#    if 'OO' in smiles:
#        singlebond.append('OO')
#    if 'CC' in branch:
#        singlebond.append('CC')
#    if 'CO' in branch or 'OC' in branch:
#        singlebond.append('CO')
#    if 'OO' in branch:
#        singlebond.append('OO')
#    if 'C' in smiles or 'C' in branch:
#        singlebond.append('CH')
#    if 'O' in smiles or 'O' in branch:
#        singlebond.append('OH')
#    if 'O' in smiles or 'O' in branch:
#        singlebond.append('OH')
#    singlebond.append('HH')
#    basis = []
#    for bond in doublebond:
#        if 'CO' == bond:
#            basis.append('C(=O)=O')
#        if 'CC' == bond:
#            basis.append('C=C')
#        if 'SO' == bond:
#            basis.append('O=S=O')
#    for bond in singlebond:
#        if 'CO' == bond:
#            basis.append('CO')
#        if 'CC' == bond:
#            basis.append('CC')
#        if 'OO' == bond:
#            basis.append('[O][O]')
#        if 'CH' == bond:
#            basis.append('C')
#        if 'OH' == bond:
#            basis.append('O')
#        if 'HH' == bond:
#            basis.append('[H][H]')
#    return basis[:count+1]
    
def is_auto(item):
    """
    Checks if a parameter should be automatically determined
    """
    if type(item) == float:
       if item == 9999.9:
           return True
    elif type(item) == str:
       if 'auto' in item.lower():
           return True
    return False

def getname_fromdirname():
    """
    Sets stoichiometry as the directory name
    if no -s is specified
    """
    cwd = os.getcwd()
    return cwd.split('/')[-1]

def nest_2_dic(bas,key1,key2=0):
    """
    Returns a dictionary value that requires two key values (is in singly nested loop )
    """
    #from heatform_db import db
    from newdb import db
 
    dicfound = False
    bas = ob.get_slabel(bas)
    for dic in db:
        if dic['_id'] == bas:
            dicfound = True
            break
    if dicfound:
        if key1 in dic:
            if key2 == 0:
                key2 = 'ANL1'
                if key2 in dic[key1]:
                    logging.info('{} {}: {:.2f}  pulled from {} in heatform_db'.format(bas, key1, dic[key1][key2], key2))
                    return dic[key1][key2]
                key2 = 'ATcTsig'
                if key2 in dic[key1]:
                    if dic[key1][key2] < 1:
                        if 'ATcT' in dic[key1]:
                            logging.info('{} {}: {:.2f}  pulled from {} in heatform_db (sigma = {:.2f})'.format(bas, key1, dic[key1]['ATcT'], 'ATcT', dic[key1][key2]))
                            return dic[key1]['ATcT']
                key2 = 'ANL0'
                if key2 in dic[key1]:
                    logging.info('{} {}: {:.2f}  pulled from {} in heatform_db'.format(bas, key1, dic[key1][key2], key2))
                    return dic[key1][key2]
            elif key2 == 'ANL0' or key2 == 'ANL1' or key2 == 'ATcT' or key2 == 'ATcTsig':
                if key2 in dic[key1]:
                    return dic[key1][key2]
                else:
                    return None
        logging.info('Value for ' + str(bas) + ',' + str(key1) + ' not found -- ommitting its contribution')
    return None

def get_gaussian_zmat(filename):
    """
    Forms a zmat from a gaussian logfile
    """
    lines  = io.read_file(filename)
    return pa.gaussian_zmat(lines)

def build_gauss(mol, theory, basisset, zmat='none', directory=None, opt=False, freq=False, anharm=False):
    """
    Builds a Guassian optimization inputfile for a molecule
    INPUT:
    dic     - dictionary for the molecule
    theory  - theory we want to optimize with
    basisset- basis set we want to optimize with
    OUPUT:
    None (but an inputfile now exists with name <stoich>.inp)
    """
    gauss  = '%Mem=25GB\n%nproc=8\n'
    gauss += '#P ' + theory.lstrip('R').lstrip('U') + '/' +  basisset +  ' int=ultrafine scf=verytight nosym '
    if opt:
        gauss += 'opt=internal '
    if freq:
        gauss += 'freq=VibRot '
    if anharm:
        gauss += 'freq=anharmonic '
    gauss += '\n'

    gauss += '\nEnergy for HeatForm\n\n'

    if type(mol) == dict:
        dic    = mol
        mol    = ob.get_mol(dic['_id'])
        stoich = dic['stoich']
        charge = dic['charge']
        mult   = dic['mult']
        #zmat   = dic['zmat']
    else:
        stoich = ob.get_formula(mol)
        mult   = ob.get_multiplicity(mol)
        charge = ob.get_charge(mol)

    if zmat == 'none':
        zmat = '\n' + ob.get_zmat(mol)
    elif zmat.startswith('geometry'):
        zmat = zmat.replace('geometry={angstrom','')
        zmat = zmat.replace('}','Variables:')
        zmat = zmat.replace('=','')
    zmat = str(charge) + '  ' + str(mult) + zmat
    gauss += zmat.lstrip('\n')

    directory = io.db_logfile_dir(directory)
    filename  = io.join_path(directory, stoich + '.inp')
    io.write_file(gauss, filename)

    return filename

def build_molpro(mol, theory, basisset, zmat='none', directory=None, opt=False, freq=False,anharm=False):
    """
    Builds a Guassian optimization inputfile for a molecule
    INPUT:
    dic     - dictionary for the molecule
    theory  - theory we want to optimize with
    basisset- basis set we want to optimize with
    OUPUT:
    None (but an inputfile now exists with name <stoich>.inp)
    """
    if type(mol) == dict:
        dic = mol
        mol = ob.get_mol(dic['_id'])
        stoich = dic['stoich']
        mult   = dic['mult']
    else:
        dic = {'g09':'nope'}
        stoich = ob.get_formula(mol)
        mult   = ob.get_mult(mol)

    molp  = 'memory,         200 ,m\nnosym\n'
    
    if zmat.startswith('geometry'):
        molp += zmat
    else:
        if zmat == 'none':
            zmat = '\n' + ob.get_zmat(mol).replace('=', '') 

        zmat1 = '\n'.join(zmat.split('Variables:')[0].split('\n')[1:] )
        zmat2 = '}\n'

        for line in zmat.split('Variables:')[1].split('\n')[1:-1]:
            zmat2 += line.split()[0] + '  =  ' + line.split()[1] + '\n'

        molp += 'geometry={angstrom\n' + zmat1 + zmat2
    spin  = 1/2.*(mult -1) * 2
    molp += '\nSET,SPIN=     ' + str(spin) + '\n\n'

    if spin == 0:
        molp += '!closed shell input\nbasis=' + basisset + '\nhf\n' 
        if 'ccsd' in theory.lower() or 'cisd' in theory.lower() or 'hf' in theory.lower() or 'mp' in theory.lower():
             pass
        else:
            molp += 'dft=[' + theory.lower() + ']\n'
            theory = 'dft'
        molp += theory.lower() 
    else:
        if 'ccsd' in theory.lower() or 'cisd' in theory.lower(): 
	    molp += '!open shell input\nbasis=' + basisset + '\nhf\nu' + theory.lower() 
        elif 'hf' in theory.lower():
	    molp += '!open shell input\nbasis=' + basisset + '\nuhf'
        elif 'mp' in theory.lower():
	    molp += '!open shell input\nbasis=' + basisset + '\nuhf\n' + 'u' + theory.lower() 
        else:
            molp += '!open shell input\nbasis=' + basisset + '\ndft=['+theory.lower() + ']\nuhf\ndft'
    if opt:
        molp += '\noptg'
    if freq:
        molp += '\nfrequencies'
    if anharm:
        molp += '\nlabel1\n{hf\nstart,atden}\n{mp2\ncphf,1}\n{surf,start1D=label1,sym=auto\n intensity,dipole=2}\nvscf\nvci'
        #molp += ',pmp=3,version=1\nvci,pmp=3,version=2\nvci,pmp=3,version=3'
    molp += '\nENERGY=energy'

    directory = io.db_logfile_dir(directory)
    filename  = io.join_path(directory, stoich + '.inp')
    io.write_file(molp, filename)
 
    return filename 

def run_gauss(filename):
    """
    Runs Gaussian on file: filename
    """
    os.system('soft add +g09; g09 ' + filename.replace('[','\[').replace(']','\]'))
    return

def run_molpro(filename):
    """
    Runs molpro on file: filename
    """
    os.system('/home/elliott/bin/molprop ' + filename.replace('[','\[').replace(']','\]'))
    
    return

def run_opt(mol, prog, meth, bas, mol_is_smiles=True):
    """
    Runs a g09 or molpro job for 
    INPUT:
    stoich   - stoichiometry of molecule
    theory   - theory energy should be listed or computed with
    basisset - basis set energy should be listed or computed with
    prog     - program an energy should be listed or computed with
    OUTPUT:
    E        - energy 
    """
    if mol_is_smiles:
        dic    = ob.get_mol(mol)
        stoich = ob.get_formula(mol)
        
    else:     #mol is dictionary
        stoich = mol['stoich']

    dic    = mol      
    
    directory = io.db_opt_path(prog, meth, bas, None, mol)
    
    if  prog == 'gaussian':
        print 'Running G09 optimization on ' + stoich + ' at ' + meth.lstrip('R').lstrip('U') + '/' + bas
        filename = build_gauss(dic, meth, bas, directory=directory, opt=True)
        run_gauss(filename)
        io.parse_all(mol, io.read_file(filename.replace('.inp','.log')))
        zmat = io.db_get_opt_prop(mol, 'zmat', db_location=directory)
        return zmat

    elif prog == 'molpro':
        print 'Running Molpro optimization on ' + stoich + ' at ' + meth.lstrip('R').lstrip('U') + '/' + bas
        filename = build_molpro(dic, meth, bas, directory=directory, opt=True)
        run_molpro(filename)
        io.parse_all(mol, io.read_file(filename.replace('.inp','.out')))
        zmat = io.db_get_opt_prop(mol, 'zmat', db_location=directory)
        return zmat

    return 

def run_energy(mol, optprog, optmeth, optbas, propprog, propmeth, propbas, entype, mol_is_smiles=True):
    """
    Runs a g09 or molpro job for 
    INPUT:
    stoich   - stoichiometry of molecule
    theory   - theory energy should be listed or computed with
    basisset - basis set energy should be listed or computed with
    prog     - program an energy should be listed or computed with
    OUTPUT:
    E        - energy 
    """
    if mol_is_smiles:
        dic    = ob.get_mol(mol)
        stoich = ob.get_formula(mol)
        
    else:     #mol is dictionary
        stoich = mol['stoich']

    dic    = mol      
    
    optdir    = io.db_opt_path(optprog, optmeth,  optbas, None, mol)
    coordfile = io.join_path(optdir, mol + '.zmat')
    if io.check_file(coordfile):
        zmat  = io.db_get_opt_prop(mol, 'zmat', db_location=optdir)
    else:
        zmat  = 'none'
    directory = io.db_sp_path(propprog, propmeth, propbas, None, mol, optprog, optmeth, optbas)
 
    anharm = False
    if 'zpve' in entype:
        freq = True
        if ob.get_natom(ob.get_mol(mol)) < 3:
            entype = 'zpve'
        if 'an' in entype and mol != 'C': 
            anharm = True
    else:
        freq = False

    if  propprog == 'gaussian':
        print 'Running G09 ' + entype + ' on ' + stoich + ' at ' + propmeth.lstrip('R').lstrip('U') + '/' + propbas
        filename = build_gauss(dic, propmeth, propbas, zmat=zmat, directory=directory, freq=freq, anharm=anharm)
        run_gauss(filename)
        io.parse_all(mol, io.read_file(filename.replace('.inp','.log')), optprog, optmeth, optbas)
        if io.check_file(directory + '/' + mol + '.' + entype):
            E = float(io.db_get_sp_prop(mol, entype, db_location=directory))
        elif entype == 'anzpve':
            E = float(io.db_get_sp_prop(mol, 'zpve', db_location=directory))
        return E

    elif propprog == 'molpro':
        print 'Running Molpro ' + entype + ' on ' + stoich + ' at ' + propmeth.lstrip('R').lstrip('U') + '/' + propbas
        filename = build_molpro(mol, propmeth, propbas, zmat=zmat, directory=directory, freq=freq,anharm=anharm)
        run_molpro(filename)
        io.parse_all(mol, io.read_file(filename.replace('.inp','.out')), optprog, optmeth, optbas)
        E = float(io.db_get_sp_prop(mol, entype, db_location = directory))
        return E

def find_E(bas, opt, en, freq, runE=True, anharm=False, dbdir='./'):
    """
    Checks a dictionary and directories for energy at a specified level of theory and basisset 
    Computes energy if it isn't in dictionary already
    INPUT:
    bas   - smiles string of molecule or basis moleule
    opt   - optlevel  (prog/method/basis)
    en    - enlevel   (prog/method/basis)
    freq  - freqlevel (prog/method/basis)
    runE  - should execute package if energies are not found
    anharm- anharmonic frequencies? True or False
    OUTPUT:
    E     - energy 
    """
    
    ### Check dictionary ###
    from heatform_db import db
    E, zpve = 0, 0
    zpvetype = 'zpve'
    
    #Check dictionary    
    #for dbdic in db:
    #    if dbdic['_id'] == bas:
    #      dic = dbdic  

    #if prog.lower() in dic: 
    #    if theory.lower().lstrip('r') in dic[prog.lower()]:
    #        if basisset.lower() in dic[prog.lower()][theory.lower().lstrip('r')]:
    #            E =  dic[prog.lower()][theory.lower().lstrip('r')][basisset.lower()][entype]
    #            print 'Energy ' + str(E) +  ' pulled from dictionary database' 
    #            return E

    ### Check directory ###
    if len(opt) > 1:
        optprog, optmethod, optbasis  = opt[0], opt[1], opt[2]
    else:
        optprog, optmethod, optbasis  = None, None, None
    enprog,   enmethod,  enbasis  =  en[0],  en[1],  en[2]
    cdire      = io.db_opt_path(optprog, optmethod, optbasis, None, bas)
    edire      = io.db_sp_path(enprog, enmethod, enbasis, None, bas, optprog, optmethod, optbasis)
    coordsfile = io.join_path(cdire, bas + '.zmat')
    enefile    = io.join_path(edire, bas + '.ene')

    if freq != None:
        freqprog,freqmethod, freqbasis  = freq[0], freq[1], freq[2]
        fdire      = io.db_sp_path(freqprog, freqmethod, freqbasis, None, bas, optprog, optmethod, optbasis)
        if anharm:
            zpvetype = 'anzpve'
        zpvefile = io.join_path(fdire, bas + '.' + zpvetype)
    if io.check_file(enefile):
        E = io.read_file(enefile).strip()
        logging.debug('{}    E:{:5} pulled from: {}'.format(bas, E, enefile))
        E = float(E)
    elif runE:
        if not io.check_file(coordsfile):
            run_opt(bas, optprog, optmethod, optbasis, True)
        E = run_energy(bas, optprog, optmethod, optbasis, enprog, enmethod, enbasis, 'ene')
        io.mkdir(edire)
        io.write_file(str(E), enefile)
        logging.debug('{}    E: {:5} saved to: {}'.format(bas, E,  enefile))

    if freq != None:
        if io.check_file(zpvefile):
            zpve= io.read_file(zpvefile).strip()
            zpve= float(zpve)
            logging.debug('{} ZPVE: {:5} pulled from: {}'.format(bas, zpve, zpvefile))
        elif runE:
            if not io.check_file(coordsfile):
                run_opt(bas, optprog, optmethod, optbasis, True)
            zpve = run_energy(bas, optprog, optmethod, optbasis, freqprog, freqmethod, freqbasis, zpvetype)
            io.mkdir(fdire)
            io.write_file(str(zpve), zpvefile)
            logging.debug('{} ZPVE: {:5} saved to: {}'.format(bas, zpve, zpvefile))
    else:
        logging.debug('Zero point vibrational energy NOT accounted for')
        zpve = 0

    #print 'no energy for ' +  dic['stoich'] + ' at '+ theory + '/' + basisset 
    #print 'No electronic energy found -- ommitting its contribution'
    return  float(E) + float(zpve)

def E_from_hfbasis(mol,basis,coefflist,E,opt, en, freq, anharm,dbdir='./'):
    """
    Uses the coefficients [a,b,c...] obtained from C = M^-1 S to find 
    delH(CxHyOz) = adelH(R1) + bdelH(R2) + cdelH(R3) + Eo(CxHyOz) - aEo(R1) - bEo(R2) -cEo(R3)
    where Rn are our basis molecules, delH(Rn) are their heats of formation, and Eo(Rn) are their
    electronic energies computed at the same level of theory as Eo(CxHyOz)
    INPUTS:
    mol       - molecule named stoichiometrically
    basis     - selected basis molecule list
    coefflist - coefficients [a,b,c,d...] described above
    E         - electronic energy of molecule
    OUTPUTS:
    E        - 0K heat of formation of molecule
   
    """
    for i,bas in enumerate(basis):
        h   = nest_2_dic(bas,'delHf',  0)
        if h is None:
            h = 0
        E  +=  coefflist[i] * h * ut.kj2au
        e    =  find_E(bas, opt, en, freq, anharm=anharm,dbdir=dbdir)
        E   -=  coefflist[i] * e

    return E

def E_hfbasis_QTC(mol,basis,coefflist,E,opt, en, freq, parameters):
    """
    Uses the coefficients [a,b,c...] obtained from C = M^-1 S to find 
    delH(CxHyOz) = adelH(R1) + bdelH(R2) + cdelH(R3) + Eo(CxHyOz) - aEo(R1) - bEo(R2) -cEo(R3)
    where Rn are our basis molecules, delH(Rn) are their heats of formation, and Eo(Rn) are their
    electronic energies computed at the same level of theory as Eo(CxHyOz)
    INPUTS:
    mol       - molecule named stoichiometrically
    basis     - selected basis molecule list
    coefflist - coefficients [a,b,c,d...] described above
    E         - electronic energy of molecule
    OUTPUTS:
    E        - 0K heat of formation of molecule
   
    """
    for i,bas in enumerate(basis):
        bas = ob.get_slabel(bas)
        smilesname = ob.get_smiles_filename(bas)
        formula   = ob.get_formula(bas)
        bas = bas.split('_')[0]
        h  =  nest_2_dic(bas,'delHf',  0)
        if h is None:
            smilesdir = io.join_path(parameters['database'], formula, smilesname)
            import os
            rundir =  io.join_path(*[os.getcwd().split(parameters['database'])[0],smilesdir,parameters['qcdirectory']])
            hoffile = io.join_path(*[rundir, formula + '.hofk'])
            if io.check_file(hoffile):
                h = float(io.read_file(hoffile).splitlines()[2].split()[0]) / ut.au2kcal / ut.kj2au
                logging.info('{} {}: {:5f}  pulled from {}'.format(bas, 'delHf', h, rundir))
            else:
                h = 0
        E  +=  coefflist[i] * h * ut.kj2au
        e    =  E_QTC(bas, opt, en, freq, parameters)
        if parameters['bac']:
            e += E_BAC(bas, parameters) / ut.au2kcal
        E   -=  coefflist[i] * e
    return E
    
def check(clist, basis,stoich,atomlist):
    """
    Makes sure nothing funky happened while computing coefficients
    """
    check = np.zeros(len(clist))
    statement = 'Coefficients produce correct stoichiometry\n'
    for i, c in enumerate(clist):
       check += c * get_stoich(ob.get_formula(ob.get_mol(basis[i])),atomlist)
    for i, sto in enumerate(stoich):
        if not check[i] == sto:
            statement = 'Coefficients do NOT produce correct stoichiometry'
            break
    return statement

def get_progmethbasis(level, loglines='', optlevel=''):

    level = update_level(level)
    if level == None:
        return None
    if is_auto(level):
        prog     = pa.get_prog(loglines)
        method   = pa.method(loglines)
        basisset = pa.basisset(loglines)
    elif level == 'optlevel':
        return get_progmethbasis(optlevel, loglines)
    else:
        prog, method, basisset = level.split('/')
    return [prog.replace('g09','gaussian'), method, basisset]

def update_level(level):
    if level != None:
        level = level.replace('\(','(').replace('\)',')')
    return level

def convert_to_smiles(basis):

    stoich_to_smiles = {'H2':'[H][H]','O2':'[O][O]','H2O':'O','NH3':'N','SO2':'O=S=O','CO2':'C(=O)=O',  
                        'CH3CH3': 'CC', 'C2H6':'CC','CH3OH':'CO','CH4O':'CO','H2CO':'C=O','CH2O':'C=O','CH4':'C'}
    for n, bas in enumerate(basis):
        if bas in stoich_to_smiles:
            basis[n] = stoich_to_smiles[bas]
    return basis

def AU_to_kcal(E, printout=True):
    if printout:
        lines =  '\n \t delHf(0K)'
        lines += '\n Hartree/atom \t'
        lines += str(E) + '\t'
        lines += '\n kJ/mol \t'
        lines += str(E/ .00038088) + '\t'
        lines += '\n kcal/mol \t'
        lines += str(E *  627.503) + '\t'
        logging.debug( lines)
    hf0k = E * ut.au2kcal #627.503
    return hf0k 

def comp_coefficients(molform, basis='auto'):
    ####BASIS SELECTION#### 
    basprint = 'manually select basis'
    atomlist = get_atomlist(molform)
    basisselection = 0
    if is_auto(basis[0]):
        basis = select_basis(atomlist)
        basisselection += 1
        basprint = 'automatically generate basis'
    elif basis[0] == 'basis.dat':
        basis = io.read_file('basis.dat').split()
        basprint = 'read basis from basis.dat'

    for bas in basis:
        bas = ob.get_formula(ob.get_mol(bas))
        atomlist.extend(get_atomlist(   bas))

    #COMPUTE Atomlist, stoichlist, matrix, and coefficients
    atomlist = list(set(atomlist))
    stoich = get_stoich(molform,atomlist)
    mat = form_mat(basis,atomlist)

    ##Pick a new basis if current one produces singular matrix
    for i in range(5):
        if np.linalg.det(mat) != 0:
            break

        basprint += '\nMatrix is singular -- select new basis'

        atomlist = get_atomlist(molform)
        if 'H' not in atomlist:
            atomlist.append('H')
        basis    = select_basis(atomlist,basisselection)
        basisselection += 1
       
        for bas in basis:
            bas = ob.get_formula(ob.get_mol(bas))
            atomlist.extend(get_atomlist(bas))

        atomlist = list(set(atomlist))
        stoich   = get_stoich(molform,atomlist)
        mat      = form_mat(basis,atomlist)
        basprint +='\n\nBasis is: ' + ', '.join(basis)
        logging.debug( basprint)
        #basprint +='\n'.join(['\t'.join([{}.format(el) for el in line] for line in mat])
    basprint += '\n  ' + molform + '\t\t' 
    basprint += '\t'.join([ob.get_formula(bas) for bas in basis])
    for i in range(len(mat)):
       basprint += '\n' + atomlist[i] + '  '
       basprint += str(stoich[i]) + '    \t'
       for el in mat[i]:
           basprint += str(el) + '\t'

    clist =  comp_coeff(mat,stoich)
    return clist, basis, basprint


def E_BAC(bas, parameters):
    ### Check dictionary ###
    from heatform_db import db
    import qctools as qc
    import iotools as io
    slabel = qc.get_slabel(bas)
    calcindex = parameters['calcindex']
    qckeyword = parameters['qckeyword']
    qlabel = qc.get_qlabel(qckeyword, calcindex)
    bac = 0.
    if 'bac' in parameters['all results'][slabel][qlabel]:
        bac = parameters['all results'][slabel][qlabel]['bac']
    if bac:
        logging.debug('BAC for {0} {1} = {2} kcal/mol'.format(slabel, qlabel,bac))
    else: 
        logging.error('BAC not found for {0} {1}'.format(slabel, qlabel))
    return  float(bac)



def E_QTC(bas, opt, en, freq, parameters):
    ### Check dictionary ###
    from heatform_db import db
    import qctools as qc
    import iotools as io
    natom = ob.get_natom(bas)
    slabel = qc.get_slabel(bas)
    parameters['natom'] = natom
    calcindex = parameters['calcindex']
    qckeyword = parameters['qckeyword']
    qlabel = qc.get_qlabel(qckeyword, calcindex)
    en, zpve, bac = 0., 0., 0.
    
    if qlabel in parameters['all results'][slabel]:
        if 'energy' in parameters['all results'][slabel][qlabel]:
            en = parameters['all results'][slabel][qlabel]['energy']
    if en:
        logging.debug('Energy for {0} {1} = {2} Hartree'.format(slabel, qlabel,en))
    else: 
        logging.error('Energy not found for {0} {1}'.format(slabel, qlabel))
    if qlabel in parameters['all results'][slabel]:
        if 'azpve' in parameters['all results'][slabel][qlabel]:
            zpve = parameters['all results'][slabel][qlabel]['azpve']
            zpvelabel = 'anharmonic ' + qlabel
        elif 'zpve' in parameters['all results'][slabel][qlabel]:
            zpve = parameters['all results'][slabel][qlabel]['zpve']
            zpvelabel = 'harmonic ' + qlabel
    else:
        for i in range(calcindex):
            qlabel = qc.get_qlabel(qckeyword, i)
            if qlabel in parameters['all results'][slabel]:
                if 'azpve' in parameters['all results'][slabel][qlabel]:
                    zpve = parameters['all results'][slabel][qlabel]['azpve']
                    zpvelabel = 'anharmonic ' + qlabel
                elif 'zpve' in parameters['all results'][slabel][qlabel]:
                    zpve = parameters['all results'][slabel][qlabel]['zpve']
                    zpvelabel = 'harmonic ' + qlabel
    if zpve:
        logging.debug('ZPVE (harmonic) for {0} {1} = {2} Hartree'.format(slabel,zpvelabel,zpve))
    else: 
        logging.warning('ZPVE not found for {0} {1}'.format(slabel, qlabel))
    return  float(en) + float(zpve)


def get_total_energy(mol, parameters):
    import qctools as qc
    natom = ob.get_natom(mol)
    parameters['natom'] = natom
    calcindex = parameters['calcindex']
#    parameters = qc.parse_qckeyword(parameters, calcindex)
    dbdir = parameters['database']
    anharm = parameters['anharmonic']
    if anharm:
        zpvetype = 'anzpve'    
    
    
def main_keyword(s,parameters):
    
   # mol    = ob.get_mol(s)
   # smi = ob.get_smiles(s)
    basis  = parameters['reference'].split(',')
    qckeys = parameters['qckeyword'].split(',')
    anharm = parameters['anharmonic']
    dbdir  = parameters['database']
    index  = parameters['calcindex']
    natom =  parameters['natom']
    optlevel  = 'sp'
    extrap    = False
    enlevel   = None
    freqlevel = None
    molform = ob.get_formula(s)
   
    if 'cbh' in basis[0]:
        clist, basis, basprint = cbh_coefficients(s.split('_')[0], basis[0], parameters)
    else:
        clist, basis, basprint = comp_coefficients(molform, basis)
 
#     lines =  ('\n___________________________________________________\n\n' +
#               'HEAT OF FORMATION FOR: ' + s + ' (' + molform + ')' +
#               '\nat ' + '/'.join(optlevel) + '//' + '/'.join(enlevel) + 
#               '\n\n___________________________________________________\n\n' +
#               '\nYou have chosen to ' + 
#               basprint) 
#     lines +=  '\n\nCoefficients are: '
#     lines += ', '.join(['{}'.format(co) for co in clist])
#     logging.debug(lines)
    params = parameters.copy()
    params['runthermo']=False
    params['xyzpath']=''
    params['suppress_printing']=True
    params['qckeyword'] = ','.join(qckeys[:index+1])
    E =  E_QTC(s, optlevel, enlevel, freqlevel, params)
    if parameters['bac']:
        E += E_BAC(s, params) / ut.au2kcal
    E =  E_hfbasis_QTC(molform, basis, clist, E, optlevel, enlevel, freqlevel, params)
    hf0k = AU_to_kcal(E)
    #parameters['runthermo']=True
    #parameters['suppress_printing']=False
    parameters['input']=s
    parameters['qckeyword'] = ','.join(qckeys)
    return hf0k, basis, clist

def main(mol,logfile='',basis='auto',E=9999.9,optlevel='auto/',freqlevel='optlevel',enlevel='optlevel',anharm=False,runE=False, bas_not_smiles=False):


    #Convert basis selection to smiles if it is mol. formula format
    if type(basis) is str:
        basis = basis.split()
    if bas_not_smiles:
        basis = convert_to_smiles(basis)

    ####AUTO SET NONUSERDEFINED PARAMETERS FROM LOGFILE###
    if logfile != '' and io.check_file(logfile):
        loglines = io.read_file(logfile)
    else: loglines = ''

    enlevel   = get_progmethbasis(enlevel,  optlevel = optlevel, loglines=loglines)
    freqlevel = get_progmethbasis(freqlevel,optlevel = optlevel, loglines=loglines)
    optlevel  = get_progmethbasis(optlevel, loglines = loglines)

    molform = ''
    if is_auto(mol):
        mol = getname_fromdirname() 
        molform = ob.get_formula(ob.get_mol(mol))

     
    #Search database for coords, energy, and zpve of mol and compute if its not there, unless told to parse it from logfile
    for spec in mol.split('_'):
        if not spec.startswith('m'):
            molform += ob.get_formula(ob.get_mol(spec))  #So transition states can be specified by two species
    clist, basis, basprint = comp_coefficients(molform, basis)

    ###PRINT STUFF OUT
    lines =  ('\n___________________________________________________\n\n' +
              'HEAT OF FORMATION FOR: ' + mol + ' (' + molform + ')' +
              '\nat ' + '/'.join(optlevel) + '//' + '/'.join(enlevel) + 
              '\n\n___________________________________________________\n\n' +
              '\nYou have chosen to ' + 
              basprint) 
    lines +=  '\n\nCoefficients are: '
    lines += ', '.join(['{}'.format(co) for co in clist])
    print lines
    logging.debug(lines + '\n')
    #print check(clist, basis, stoich,atomlist)

    #COMPUTE AND PRINT delH###
    if runE:
        E = find_E(mol, optlevel, enlevel, freqlevel, runE, anharm)
    elif is_auto(E):
        E = pa.energy(loglines)[1]
        logging.debug('{}    E: {:5} pulled from: {}'.format(mol, E, logfile))
        zpve = pa.zpve(loglines)
        if zpve == None:
            logging.debug( 'Zero point vibrational energy NOT accounted for')
            zpve = 0
        else:
            logging.debug( '{} ZPVE: {:5} pulled from: {}'.format(mol, zpve, logfile))
        E = E + zpve

    E =  E_from_hfbasis(molform,basis,clist,E, optlevel, enlevel, freqlevel,anharm)
    hf0k = AU_to_kcal(E)

    enprog, enmethod, enbasis = enlevel[0], enlevel[1], enlevel[2]
    optprog, optmethod, optbasis = optlevel[0], optlevel[1], optlevel[2]
    hfdire = io.db_sp_path(enprog, enmethod, enbasis, mol, None, optprog, optmethod, optbasis)
    if '_' not in mol:
        if not io.check_file(hfdire + '/' + mol + '.hf0k'):
            io.db_store_sp_prop('Energy (kcal/mol)\tBasis\n----------------------------------',mol,'hf0k',None,enprog,enmethod,enbasis, optprog, optmethod, optbasis)
        s = '\n' + str(hf0k) + '\t' + '  '.join(basis) 
        io.db_append_sp_prop(s,mol,'hf0k',None,enprog,enmethod,enbasis, optprog, optmethod, optbasis)
        
        lines  = '\n\nStored heats of formation:\n'
        lines += '----------------------------------\n' 
        lines += io.db_get_sp_prop(mol,'hf0k',None,enprog,enmethod,enbasis, optprog, optmethod, optbasis)
        lines += '\n\n_________________________________________________\n\n'
        logging.debug( lines)

    return hf0k, basis

###CBH
def get_balance(smiles, frags):
    stoichs = {}
    for frag in frags:
       stoich = qc.get_atom_stoich(ob.get_formula(frag))
       for atom in stoich:
           if atom in stoichs:
               stoichs[atom] += stoich[atom] * frags[frag]
           else:
               stoichs[atom]  = stoich[atom] * frags[frag]
    stoich = qc.get_atom_stoich(ob.get_formula(smiles))
    balance = {}
    for atom in stoich:
        if atom in stoichs:
            balance[atom] = stoich[atom] - stoichs[atom]
        else:
            balance[atom] = stoich[atom] 
    return balance

def make_balanced(smiles, frags):
    balance = get_balance(smiles, frags)
    if 'C' in balance:
        if balance['C'] != 0:
            if not 'C' in frags:
                frags['C'] = balance['C'] 
            else:
                frags['C'] += balance['C'] 
    if 'N' in balance:
        if balance['N'] != 0:
            if not 'N' in frags:
                frags['N'] = balance['N'] 
            else:
                frags['N'] += balance['N'] 
    if 'O' in balance:
        if balance['O'] != 0:
            if not 'O' in frags:
                frags['O'] = balance['O'] 
            else:
                frags['O'] = balance['O']

    balance = get_balance(smiles, frags)
    if 'H' in balance:
        if balance['H'] != 0:
            if not '[H][H]' in frags:
                frags['[H][H]'] = balance['H'] / 2.
            else:
                frags['[H][H]'] += balance['H'] / 2.
    return frags

def isUnbalanced(balance):
    bal = True
    for atom in balance:
       if abs(balance[atom]) > 0:
           bal = False
    if bal:
       return False
    return True

def get_adj_mat(lines):
    lines = lines.lower().replace('\n\n','\n')
    lines = lines.splitlines()
    rows = []
    rad  = []
    rads = []
    startkey = 'molecular structure:'
    endkey = 'z-matrix atom order:'
    endkey2 = 'free radical'
    save = False
    for line in lines:
        if line.startswith(endkey):
            save = False
        if line.startswith(endkey2):
            save = False
            if 'is' in line:
                rad = line.split('site is ')[1]
            elif 'are' in line:
                rad = line.split('sites are ')[1].split()
        if save:
            rows.append(line)
        if line.startswith(startkey):
            save = True
    
    if rad:
        for i, atom in enumerate(rows[0].split()[1:]):
            if atom.strip() in rad:
                rads.append(i)
    atoms = qc.get_atom_list(rows[0])
    mat = []
    for row in rows[1:]:
        if row:
            mat.append(row.split()[1:-1])
    return atoms, mat, rads

def bond_sum(mat, atoms, heavy, i):
    bonds = 0
    for j, col in enumerate(mat[i]):
        if atoms[j] in heavy:
            if float(col) > 0:
                bonds += int(round(float(col)-0.5))
    return bonds

def lhs_rhs(frags):
    rhs = {}
    lhs = {}
    for frag in frags:
        if frags[frag] > 0:
            rhs[frag] = frags[frag]
        elif frags[frag] < 0:
            lhs[frag] = - frags[frag]
    return lhs, rhs

def print_lhs_rhs(smiles, frags):
    lhs, rhs = lhs_rhs(frags)
    lhsprint = ob.get_formula(smiles)
    rhsprint = ''
    for frag in rhs:
        if rhsprint:
            rhsprint += ' +  {:.1f} {} '.format( rhs[frag], ob.get_formula(frag))
        else:
            rhsprint  = ' {:.1f} {} '.format( rhs[frag], ob.get_formula(frag))
    for frag in lhs:
            lhsprint += ' +  {:.1f} {} '.format( lhs[frag], ob.get_formula(frag))
    return '{} --> {}'.format(lhsprint, rhsprint)

def CBHzed(smiles, mat, atoms, heavy, bond, rads):
    frags = {}
    valence = {'C':4, 'O':2, 'N':3}
    for i, row in enumerate(mat):
       if atoms[i] in heavy:
          val = valence[atoms[i]]
          frag = ''
          frag = atoms[i]
          if i in rads:
              val -= 1
          saturate = ''
          if val == 1:
              saturate = 'H'
          elif val > 1:
              saturate = 'H{:g}'.format(val)
          frag = '[{}{}]'.format(frag, saturate)
          frag = ob.get_slabel(frag).split('_')[0]
          if frag in frags:
              frags[frag] += 1
          else: 
              frags[frag] = 1
    
    frags = make_balanced(smiles, frags)
    return frags

def CBHone(smiles, mat, atoms, heavy, bond, rads, parameters = {}):
    frags = {}
    valence = {'C':4, 'O':2, 'N':3}
    for i, row in enumerate(mat):
        for j, col in enumerate(row):
            if atoms[i] in heavy and atoms[j] in heavy:
               if float(col) > 0:
                   vali = valence[atoms[i]]
                   valj = valence[atoms[j]]
                   if i in rads:
                       vali -= 1
                   if j in rads:
                       valj -= 1
                   vali -= int(round(float(col)))
                   valj -= int(round(float(col)))
                   fragi = atoms[i]
                   fragj = atoms[j]
                   saturate = ''
                   if vali == 1:
                       saturate = 'H'
                   elif vali > 1:
                       saturate = 'H{:g}'.format(vali)
                   fragi = '[{}{}]'.format(fragi, saturate)
                   saturate = ''
                   if valj == 1:
                       saturate = 'H'
                   elif valj > 1:
                       saturate = 'H{:g}'.format(valj)
                   fragj = '[{}{}]'.format(fragj, saturate)
                   frag = fragi + bond[str(int(round(float(col))))] + fragj
                   frag = ob.get_slabel(frag).split('_')[0]
                   if frag in frags:
                       frags[frag] += 1
                   else: 
                       frags[frag] = 1
    #BALANCE
    newfrags = {}
    zedfrags = CBHzed(smiles, mat, atoms, heavy, bond, rads)
    new = {}
    frags =  {k: v for k, v in frags.items() if v}
    for frag in frags:
        newfrags[frag] = frags[frag]
        atoms, mat, rads = get_x2zparams(frag, parameters)
        new = CBHzed(frag, mat, atoms, heavy, bond, rads)
        for n in new:
           if n in newfrags:
               newfrags[n] -= new[n] * frags[frag]
           else:
               newfrags[n] =- new[n] * frags[frag]
    if not frags:
        frags = CBHzed(smiles, mat, atoms, heavy, bond, rads)
    for frag in zedfrags:
        if frag in newfrags:
            newfrags[frag] += zedfrags[frag]
        else:
            newfrags[frag]  = zedfrags[frag]
    frags = newfrags
    frags =  {k: v for k, v in frags.items() if v}
    frags = make_balanced(smiles, frags)
    return frags

def CBHtwo(smiles, mat, atoms, heavy, bond, rads):
    frags = {}
    fmat = square_mat(mat)
    valence = {'C':4, 'O':2, 'N':3}
    for i, row in enumerate(fmat):
        if atoms[i] in heavy:
            
            vali = valence[atoms[i]]
            fragi = atoms[i]
            if i in rads:
                vali -= 1
            vali -= bond_sum(fmat, atoms, heavy, i)
            saturate = ''
            if vali == 1:
                saturate = 'H'
            elif vali > 1:
                saturate = 'H{:g}'.format(vali)
            fragi = '[{}{}]'.format(fragi, saturate)
            saturate = ''
            count = 0
            frag = fragi
            for j, col in enumerate(row):
                if atoms[j] in heavy:
                    if float(col) > 0:
                        count += 1
                        valj = valence[atoms[j]]
                        if j in rads:
                            valj -= 1
                        valj -= int(round(float(col)-0.5))
                        fragj = atoms[j]
                        saturate = ''
                        if valj == 1:
                            saturate = 'H'
                        elif valj > 1:
                            saturate = 'H{:g}'.format(valj)
                        fragj = '[{}{}]'.format(fragj, saturate)
                        frag  += '(' + bond[col] + fragj + ')'
                        vali = valence[atoms[i]]
                        fragi = atoms[i]
                        if i in rads:
                            vali -= 1
                        vali -= int(round(float(col)-0.5))
                        saturate = ''
                        if vali == 1:
                            saturate = 'H'
                        elif vali > 1:
                            saturate = 'H{:g}'.format(vali)
                        fragi = '[{}{}]'.format(fragi, saturate)
                        subfrag = fragi +  bond[col] + fragj 
                        if i < j:
                            subfrag = ob.get_slabel(subfrag).split('_')[0]
                            if subfrag in frags:
                                frags[subfrag] -= 1
                            else: 
                                frags[subfrag] = -1
            frag = ob.get_slabel(frag).split('_')[0]
            if frag in frags:
                frags[frag] += 1
            else: 
                frags[frag] = 1
    frags = make_balanced(smiles, frags)
    return frags

def get_conns(i, axis, mat, atoms, heavy, bond):
    fmat = square_mat(mat)
    ret = atoms[i]
    row = fmat[i]
    for j, col in enumerate(row):
        if j != axis and atoms[j] in heavy and float(col) > 0:
             ret += '({0}{1})'.format(bond[col], atoms[j])
    return ret

def square_mat(mat):
    import copy
    fmat = copy.deepcopy(mat)
    for i, row in enumerate(fmat):
        for j, col in enumerate(fmat): 
            if i == j:
                fmat[i].append('0')
            elif i < j:
                fmat[i].append(mat[j][i])
    return fmat    

def get_recentxyz(s, results):
    xyz = ''
    slabel = ob.get_slabel(s)
    if slabel in results:
        for key in results:
            if 'xyz' in results:
                xyz = results[slabel]['xyz']
    return xyz

def choose_xyz(smiles, parameters):
    xyz = ''
    if 'all results' in parameters:
        xyz = get_recentxyz(smiles, parameters['all results'])
    if not xyz:
        if 'xyzdir' in parameters:
           xyz = io.read_xyzdir(smiles, parameters['xyzdir'])
    if not xyz:
            xyz = ob.get_xyz(smiles)
    return xyz

def get_x2zparams(smiles, parameters={}):
    """
    returns relevant x2z information for CBH calculation
    INPUT: 
    parameters -- dicitonary with qtc data (can be empty)
    OUTPUT:
    atoms      -- list of atoms in order of the x2z adjacency matrix
    mat        -- triangular connectivity matrix from x2z
    rads       -- list of radical sites zero indexed that correlate to the order of the atoms list
    EXAMPLE:
    >> get_x2zparams('C[C]=O',{})
    ['C', 'C', 'O', 'H', 'H', 'H'], [[], ['1'], ['0', '2'], ['1', '0', '0'], ['1', '0', '0', '0'], ['1', '0', '0', '0', '0']],  [1]
    """
    xyz = choose_xyz(smiles, parameters)
    io.write_file(xyz,'x2z.xyz')
    lines = qc.run_x2z('x2z.xyz', 'x2z')
    atoms, mat, rads = get_adj_mat(lines)
    return atoms, mat, rads

def cbh_coefficients(smiles, ref, parameters = {}):
   
    logging.info(smiles)
    if ob.get_natom(smiles) < 2:
        clist = [1]
        fraglist = [smiles]
        msg = 1
    else:
        atoms, mat, rads = get_x2zparams(smiles, parameters)
        heavy = ['C','O','N']
        bond  = {'1':'','2':'=','3':'#','4':'$','1.5':':','2.5':'=','1.7':'=','2.3':'=','1.3':'','2.7':'#'}
        
        msg = ''
        if ref.lower() == 'cbh0':
            frags = CBHzed(smiles, mat, atoms, heavy, bond, rads)
        elif ref.lower() == 'cbh1':
            frags = CBHone(smiles, mat, atoms, heavy, bond, rads, parameters)
        elif ref.lower() == 'cbh2':
            frags = CBHtwo(smiles, mat, atoms, heavy, bond, rads)

        output = print_lhs_rhs(smiles, frags)
        if isUnbalanced(get_balance(smiles, frags)):
            output += '\nERROR: unbalanced fragmentation'
        msg += output
        msg += '\n'

        fraglist = []
        clist = []
        for frag in frags:
            fraglist.append(frag)
            clist.append(frags[frag])
    return clist, fraglist, msg

def get_bac(parameters, mylist, samppercent = 0, errthresh = 10):
    """
    Computes BAC parameters (corrections to the HoF for each bond type for each level of theory)
    INPUT
    parameters  -- dictionary containing species info
    mylist      -- list of species
    samppercent -- percentage of species you want to form the set with (0 being all, 1 being none)
    errthresh   -- threshold for what species will be ignored when forming BAC based on kj difference from ANL0 
    OUTPUT
    out         -- msg to log
    bac.txt     -- BAC parameters file
    """
    out = 'begin least squares\n'
    bondlabel = []
    mols = []
    taskmat = {}
    testset = ''
    bacset = ''
    for i,s in enumerate(mylist):
        anl = nest_2_dic(qc.get_slabel(s), 'delHf','ANL0')
        bonds = None
        if anl:
            if np.random.rand() > samppercent:
                sresults = parameters['all results'][qc.get_slabel(s)]
                for qcresultkey, qcresultval in sorted(sresults.iteritems(),key= lambda x: x[0]):
                    logging.debug(s)
                    task = qcresultkey#qcresultval['chemkin'].splitlines()[0].split()[2]
                    if task in taskmat:
                        mat = taskmat[task]['mat']
                        err = taskmat[task]['errors']
                    else:
                        mat = []
                        err = []
                        taskmat[task] = {}
                        taskmat[task]['mat'] = mat
                        taskmat[task]['errors'] = err
                    energy = qcresultval['deltaH0']
                    error = anl / ut.kcal2kj - energy
                    #error = energy - anl/ ut.kcal2kj
                    if abs(error) < errthresh:
                        mols.append(qc.get_slabel(s))
                        natom = ob.get_natom(ob.get_mol(s))
                        if not bonds:
                            if natom > 1:
                                xyz = ''
                                if 'xyz' in qcresultval:
                                    xyz = qcresultval['xyz']
                                if not xyz:
                                    xyz = ob.get_xyz(ob.get_mol(s))
                                x2zout = qc.run_x2z(xyz, parameters['x2z'])
                                bacset += '{}\n'.format( qc.get_slabel(s)  )
                                x2zout = qc.run_x2z(xyz, parameters['x2z'])
                                bonds =  qc.get_x2z_bonds(x2zout)
                            else:
                                bonds = {}
                            hofbasis = qcresultval['heat of formation basis']
                            hofcoeff = qcresultval['heat of formation coeff']
                            hofbonds = {}
                            for c, bas in enumerate(hofbasis):
                                nnatom = ob.get_natom(ob.get_mol(bas))
                                if nnatom > 1:
                                    xyz = ''
                                    if 'xyz' in parameters['all results'][qc.get_slabel(bas)][task]:
                                        xyz = parameters['all results'][qc.get_slabel(bas)][task]['xyz']
                                    if not xyz:
                                        xyz = ob.get_xyz(ob.get_mol(s))
                                    x2zout = qc.run_x2z(xyz, parameters['x2z'])
                                    hofbonds =  qc.get_x2z_bonds(x2zout)
                                    for bond in hofbonds:
                                        if bond in bonds:
                                            bonds[bond] -= hofcoeff[c] * hofbonds[bond]
                                        else:
                                            bonds[bond] =- hofcoeff[c] * hofbonds[bond]
                        bondarray = np.zeros(100)
                        for bond in bonds:
                            new = True
                            j = 0
                            for label in bondlabel:
                                if bond == label:
                                    bondarray[j] = bonds[bond]
                                    new = False 
                                    break
                                j += 1
                            if new:
                                bondlabel.append(bond)
                                bondarray[j] = bonds[bond]
                        mat.append(bondarray.tolist())
                        err.append(error)
                        taskmat[task]['mat'] = mat
                    else:
                        out += 'ommitting {} from BAC, error is enormous\n'.format(s)
                    taskmat[task]['errors'] = err
            else:
                testset += '{}\n'.format(qc.get_slabel(s))
    io.write_file(bacset, 'formset.txt')
    io.write_file(testset, 'testset.txt')
    bacout = 'Formed from: {}\n'.format( ' ,'.join(mols)) 
    for task in taskmat:
        mat = taskmat[task]['mat']
        mat = [mat[k][:len(bondlabel)] for k in range(len(mat))]
        errors =  taskmat[task]['errors']
        findbac = True
        try:
            #print task
            #print 'AVG \t {:.4f} \nRMS \t {:.4f}'.format(np.average(errors),np.sqrt(np.average(np.square(errors))))
            coeff = np.linalg.lstsq(mat, errors)[0]
            bacout +=  '\n' + task + '\n'
            for k, label in enumerate(bondlabel):
                bacout += '{}  {:4.3f}\n'.format(label, coeff[k])
            out += 'BAC parameters successfully computed for {}, appended to bac.txt\n'.format(task)
        except:
           out += 'Cannot compute BAC for {}, nan found in errors\n'.format(task)
    io.write_file(bacout, 'strictestbac.txt')
    return out
      
def calc_bac(parameters, bonds, task):
    """
    Given a set of bonds, calculate the bond additivity correction
    using a  template specified in parameters
    """
    correction = 0.
    bacfile = parameters['bacdirectory'] +'/bac.txt'
    #bacfile = 'bac.txt'
    bacfile = io.read_file(bacfile)
    if task + '\n' in bacfile:
        bacfile = bacfile.split(task + '\n')[1].split('\n\n')[0]
        bacfile = bacfile.splitlines()
        bacdic  = {}
        for line in bacfile:
            line = line.split()
            if len(line) > 1:
                bacdic[line[0]] = float(line[1])
        for key in bonds:
            if key in bacdic:
                correction += bonds[key] * bacdic[key] 
            else:
                logging.info( '{} bond has no correction'.format(key))
    return correction

if __name__ == '__main__': 
    """
    Run heat of formation code
    """
    #SET PARAMETERS############
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                            description="SARAH!!! you haven't done this yet!!!")

    parser.add_argument('-m','--molecule',           type=str,   default='auto')

    parser.add_argument('-l','--logfile',            type=str,  default='geoms/reac1_l1.log')

    parser.add_argument('-b','--select_stoich_basis',type=str,   default='auto')
    parser.add_argument('-B','--select_smiles_basis',type=str,   default='auto')#Doesn't actually matter if you use b or B right now

    parser.add_argument('-E','--electronic_energy',  type=float, default=9999.9)
    parser.add_argument('-o','--optlevel',          type=str,   default='auto/')
    parser.add_argument('-f','--freqlevel',       type=str,  default='optlevel')
    parser.add_argument('-e','--energylevel',     type=str,  default='optlevel')

    parser.add_argument('-d','--database',           type=str, default='heatform_db')
    parser.add_argument('-s','--mol_not_smiles',            action='store_true')
    parser.add_argument('-r','--run_energy',                action='store_true')
    parser.add_argument('-a','--anharmonic',                action='store_true')

    ###########################
    args = parser.parse_args()

    mol    = args.molecule
    E      = args.electronic_energy
    runE   = args.run_energy
    basis  = args.select_stoich_basis
    bas_not_smiles = True
    if basis == 'auto':
        bas_not_smiles = False
        basis = args.select_smiles_basis
    optlevel  = args.optlevel
    freqlevel = args.freqlevel
    enlevel   = args.energylevel
    logfile= args.logfile
    anharm = args.anharmonic
    db     = args.database
    mol_not_smiles = args.mol_not_smiles
  
    main(mol,logfile=logfile, basis=basis, E=E, optlevel=optlevel, freqlevel=freqlevel, enlevel=enlevel, anharm=anharm, runE=runE,bas_not_smiles=bas_not_smiles) 
