#!usr/bin/python

import re
import numpy as np
import os
import sys
sys.path.insert(0, '/home/elliott/Packages/QTC/')
import iotools as io
import patools as pa
import obtools as ob

"""
Heatform determines the heat of formation for a molecule by using a basis of 
molecules with well-determined heats of formation

To use this script from command line you can do the following:

(1) python heatform.py -l <logfile.log>
    this will automatically determine what stoichiometry, program, level of theory, basis set,
    and electronic energy heatform needs as parameters 

(3) python heatform.py -s <smiles> -r -t <method>/<basisset> 
    will run an energy computation for the smiles molecule at method/basis level of theory 
    before it starts the basis molecules (-p will let you choose g09 or molpro)

(4) python heatform.py -s <stoichiometry> -p <program> -t <method>/<basisset> -e <energy>
    manually set all these parameters so that no logfile or first computation is required

(3) User can also use -b <R1 R2 R3> to manually select the basis of molecules

To use this as a module in another script:

(1) main(stoich, logfile, electronic_energy, basisset, theory, database, prog, runE) will return 
    the heat of formation in kcal.  If logfile is set, the rest can be auto (see defaults), and if 
    logfile is gibberish, the rest must be set.

"""
#  TO DO:
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
    if 'N' in atomlist and i <= count:
        basis.append('N')
        #basis.append('NH3')
        i += 1
    if 'S' in atomlist and i<= count:
        #basis.append('SO2') 
        basis.append('O=S=O') 
        i += 1
    if 'H' in atomlist  and i<= count and attempt < 1:
        #basis.append('H2')
        basis.append('[H][H]')
        i += 1
    elif 'H' in atomlist and 'C' not in atomlist and i<= count and attempt < 2:
        #basis.append('H2')
        basis.append('[H][H]')
        i += 1
    if 'O' in atomlist and i<= count and attempt < 2:
        #basis.append('O2')
        basis.append('[O][O]')
        i += 1
    if 'C' in atomlist and i<= count and attempt < 3:
        #basis.append('CH4')
        basis.append('C')
        i += 1
    if 'O' in atomlist and 'H' in atomlist and  i<= count and attempt  < 3:
        #basis.append('H2O')
        basis.append('O')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count and attempt < 4:
        #basis.append('CO2')
        basis.append('C(=O)=O')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count and attempt < 4:
        #basis.append('H2CO')
        basis.append('C=O')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count:
        #basis.append('CH3OH')
        basis.append('CO')
        i += 1
    if 'C' in atomlist and i<= count:
        #basis.append('CH3CH3')
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
    """
    cwd = os.getcwd()
    return cwd.split('/')[-1]

###FOR IF WE ARE USING A LOGFILE
def getprog_fromlogfile(lines):
    """
    returns name of program for the lines of a logfile
    """ 
    return pa.get_prog(lines)

def gettheory_fromlogfile(logfile):
    """
    return the last method run in the logfile
    """
    lines = io.read_file(logfile)
    return pa.method(lines).lstrip('R').lstrip('r') 

def getbasisset_fromlogfile(logfile):
    """
    returns basisset used in logfile
    """
    lines = io.read_file(logfile)
    return pa.basisset(lines)

def getenergy_fromlogfile(logfile,theory='',prog=''):
    """
    Returns by default final energy solved for in logfile, or 
    the energy corresponding to a given method in logfile
    """
    lines = io.read_file(logfile)
    print lines
    print prog
    if prog == '':
        prog  = getprog_fromlogfile(lines)
    if prog == 'g09':
       return pa.gaussian_energy(lines,theory)[1]
    elif prog == 'molpro':
       return pa.molpro_energy(lines,theory)[1]
    elif prog == None:
       return pa.molpro_energy(lines,theory)[1]
    print 'no energy found for this molecule, please use -e to manually set it'
    return

def nest_2_dic(bas,key1,key2):
    """
    Returns a dictionary value that requires two key values (is in singly nested loop )
    """
    from testdb import db

    for dic in db:
        if dic['_id'] == bas:
            break

    if key1 in dic:
        if key2 in dic[key1]:
            return dic[key1][key2]
    print 'Value for ' + str(key1) + ',' + str(key2) + ' not found -- ommitting its contribution'

    return 0

def get_gaussian_zmat(filename):
    """
    Forms a zmat from a gaussian logfile
    """
    lines  = io.read_file(filename)
    return pa.gaussian_zmat(lines)

def build_gauss(mol, theory, basisset,directory=None):
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
    gauss += '#P ' + theory.lstrip('R').lstrip('U') + '/' +  basisset +  ' opt=internal int=ultrafine scf=verytight nosym\n'

    gauss += '\nEnergy for HeatForm\n\n'

    meths = ['ccsdt','ccsd(t)','ccsd','m062x','b3lyp']
    bases = ['cc-pvqz','cc-pvtz','cc-pvdz','6-311+g(d,p)','6-31+g(d,p)']
    zmat  = 'none'
    
    if type(mol) == dict:
        dic    = mol
        mol    = ob.get_mol(dic['_id'])
        stoich = dic['stoich']
        charge = dic['charge']
        mult   = dic['mult']
    else:
        dic = {'g09':'nope'}
        stoich = ob.get_formula(mol)
        mult   = ob.get_multiplicity(mol)
        charge = ob.get_charge(mol)

    if theory.lower().lstrip('r') in dic['g09']:
        for j in range(len(bases)):
            if bases[j] in dic['g09'][theory.lower().lstrip('r')]:
                if zmat in dic['g09'][theory.lower().lstrip('r')][bases[j]]:
                    zmat = dic['g09'][theory.lower().lstrip('r')][bases[j]]['zmat']
    if zmat == 'none':
        for i in range(len(meths)):
            if meths[i] in dic['g09']:
                for j in range(len(bases)):
                    if bases[j]  in dic['g09'][meths[i]]:
                        if 'zmat'in dic['g09'][meths[i]][bases[j]]:
                            zmat =  dic['g09'][meths[i]][bases[j]]['zmat']

    if zmat == 'none':
        zmat = ob.get_zmat(mol)
        gauss += str(charge) + ' ' + str(mult) + '\n'
    gauss += zmat.lstrip('\n')

    directory = io.db_logfile_dir(directory)
    filename  = io.join_path(directory, stoich + '.inp')
    io.write_file(gauss, filename)

    return filename

def build_molpro(mol, theory, basisset, directory=None):
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

    molp  = 'memory,         200 ,m\nnosym\ngeometry={angstrom\n'
    meths = ['ccsdt','ccsd(t)','ccsd','m062x','b3lyp']
    bases = ['cc-pvqz','cc-pvtz','cc-pvdz','6-311+g(d,p)','6-31+g(d,p)']
    zmat  = 'none'
    
    if theory.lower().lstrip('r') in dic['g09']:
        for j in range(len(bases)):
            if bases[j] in dic['g09'][theory.lower().lstrip('r')]:
                if zmat in dic['g09'][theory.lower().lstrip('r')][bases[j]]:
                    zmat = dic['g09'][theory.lower().lstrip('r')][bases[j]]['zmat']
    if zmat == 'none':
        for i in range(len(meths)):
            if meths[i] in dic['g09']:
                for j in range(len(bases)):
                    if bases[j]  in dic['g09'][meths[i]]:
                        if 'zmat'in dic['g09'][meths[i]][bases[j]]:
                            zmat =  dic['g09'][meths[i]][bases[j]]['zmat']
    if zmat == 'none':
        zmat = '\n\n' + ob.get_zmat(mol).replace('=', '') 

    zmat1 = '\n'.join(zmat.split('Variables:')[0].split('\n')[2:] )
    zmat2 = '}\n'

    for line in zmat.split('Variables:')[1].split('\n')[1:-1]:
        zmat2 += line.split()[0] + '  =  ' + line.split()[1] + '\n'

    molp += zmat1 + zmat2
    spin  = 1/2.*(mult -1) * 2
    molp += '\nSET,SPIN=     ' + str(spin) + '\n\n'

    if spin == 0:
        molp += '!closed shell input\nbasis=' + basisset + '\nhf\n' 
        if 'ccsd' in theory.lower() or 'cisd' in theory.lower() or 'hf' in theory.lower() or 'mp' in theory.lower():
             pass
        else:
            molp += 'dft=[' + theory.lower() + ']\n'
            theory = 'dft'
        molp += theory.lower() + '\noptg\nENERGY=energy'
    else:
        if 'ccsd' in theory.lower() or 'cisd' in theory.lower(): 
	    molp += '!open shell input\nbasis=' + basisset + '\nhf\nu' + theory.lower() + '\noptg\nENERGY=energy'
        elif 'hf' in theory.lower():
	    molp += '!open shell input\nbasis=' + basisset + '\nuhf\noptg\nENERGY=energy'
        elif 'mp' in theory.lower():
	    molp += '!open shell input\nbasis=' + basisset + '\nrhf\n' + theory.lower() + '\noptg\nENERGY=energy'
        else:
            molp += '!open shell input\nbasis=' + basisset + '\ndft=['+theory.lower() + ']\nuhf\ndft\noptg\nENERGY=energy'

    directory = io.db_logfile_dir(directory)
    filename  = io.join_path(directory, stoich + '.inp')
    io.write_file(molp, filename)
 
    return filename 

def run_gauss(filename):
    """
    Runs Gaussian on file: filename
    """
    os.system('soft add +g09; g09 ' + filename.replace('(','\(').replace(')','\)').replace('[','\[').replace(']','\]'))
    return

def run_molpro(filename):
    """
    Runs molpro on file: filename
    """
    os.system('/home/elliott/bin/molprop ' + filename)
    
    return

def run_energy(mol, theory, basisset, prog, mol_is_smiles=True):
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
    
    directory = io.db_sp_path(prog, theory, basisset, None, mol, prog, theory, basisset)

    if  prog == 'g09':
        print 'Running G09 on ' + stoich + ' at ' + theory.lstrip('R').lstrip('U') + '/' + basisset
        filename = build_gauss(dic, theory, basisset, directory)
        run_gauss(filename)
        io.parse_all(mol, io.read_file(filename.replace('.inp','.log')))
        E = float(io.db_get_sp_prop(mol, 'ene', db_location=directory))
        print 'Energy found to be: ' + str(E)
        return E

    elif prog == 'molpro':
        print 'Running Molpro on ' + stoich + ' at ' + theory.lstrip('R').lstrip('U') + '/' + basisset
        filename = build_molpro(mol, theory, basisset, directory)
        run_molpro(filename)
        io.parse_all(mol, io.read_file(filename.replace('.inp','.out')))
        E = float(io.db_get_sp_prop(mol, 'ene', db_location=directory))
        print 'Energy found to be: ' + str(E)
        return E

def E_dic(bas,energORs,theory,basisset,prog):
    """
    Checks a dictionary and directories for energy at a specified level of theory and basisset 
    Computes energy if it isn't in dictionary already
    INPUT:
    bas      - smiles string of molecule or basis moleule
    energORs - are we checking for the energy or its uncertainty (uncertainty broken)
    theory   - theory energy should be listed or computed with
    basisset - basis set energy should be listed or computed with
    prog     - program an energy should be listed or computed with
    OUTPUT:
    E        - energy 
    """
    
    ### Check dictionary ###
    from testdb import db

    for dic in db:
        if dic['_id'] == bas:
            break

    if prog.lower() in dic: 
        if theory.lower().lstrip('r') in dic[prog.lower()]:
            if basisset.lower() in dic[prog.lower()][theory.lower().lstrip('r')]:
                E =  dic[prog.lower()][theory.lower().lstrip('r')][basisset.lower()][energORs]
                print 'Energy ' + str(E) +  ' pulled from dictionary database' 
                return E

    ### Check directory ###

    dire = io.db_sp_path(prog, theory, basisset, None, bas, prog, theory, basisset)
    enefile = io.join_path(dire, bas+ '.ene')
    if io.check_file(enefile):
        E = float(io.read_file(enefile).strip())
        print 'Energy ' + str(E) +  ' pulled from: ' + enefile
        return E

    if energORs == 'energy':
        E = run_energy(bas, theory, basisset, prog)
        io.db_store_sp_prop(str(E), bas, 'ene', dire)
        print 'Energy ' + str(E) +  ' saved to: ' + enefile
        return E

    elif energORs != 'energy':
        return .001
    #print 'no energy for ' +  dic['stoich'] + ' at '+ theory + '/' + basisset 
    print 'No electronic energy found -- ommitting its contribution'
    return 0

def comp_energy(mol,basis,coefflist,E,theory,basisset,prog):
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
    zE        - 0K heat of formation of molecule
    sigma     - uncertainty estimate (i.e., delH = lE +/- sigma)
   
    """
    zE        = E
    var       = 0
    for i,bas in enumerate(basis):
        zE  +=  coefflist[i] * nest_2_dic(bas,'HeatForm',  0) * 0.00038088

        var += (coefflist[i] * nest_2_dic(bas,'HeatForm','sigma') * 0.00038088   )**2
        E    =  E_dic(bas, 'energy',theory,basisset,prog)
        zE  -=  coefflist[i] * E
        #var += (coefflist[i] * E * E_dic(bas,'sigma',theory,basisset,prog) )**2

    sigma = np.sqrt(var)
    return zE, sigma
    
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

def main(mol,logfile='geoms/reac1_l1.log', basis='auto', E=9999.9, theory='auto/',db='tempdb',prog = 'auto', runE=False, mol_not_smiles=False):
   

    #Convert basis selection to smiles if it is mol. formula format
    basis = basis.split()
    stoich_to_smiles = {'H2':'[H][H]','O2':'[O][O]','H2O':'O','NH3':'N','SO2':'O=S=O','CO2':'C(=O)=O',  
                        'CH3CH3': 'CC', 'C2H6':'CC','CH3OH':'CO','CH4O':'CO','H2CO':'C=O','CH2O':'C=O','CH4':'C'}
    for n, bas in enumerate(basis):
        if bas in stoich_to_smiles:
            basis[n] = stoich_to_smiles[bas]

    ####AUTO SET NONUSERDEFINED PARAMETERS FROM LOGFILE###
    if is_auto(prog):
        prog = pa.get_prog(io.read_file(logfile))

    if is_auto(mol):
        mol = getname_fromdirname() 
        molform = ob.get_formula(ob.get_mol(mol))

    if is_auto(theory) and io.check_file(logfile):
        theory  = gettheory_fromlogfile(logfile)
        theory += '/'
        theory += getbasisset_fromlogfile(logfile)
    theory, basisset = theory.split('/')
    
    #Search database for energy of mol and compute if its not there, unless told to parse it from logfile
    if not mol_not_smiles:
        dire = io.db_sp_path(prog, theory, basisset, None, mol, prog, theory, basisset)
        enefile = io.join_path(dire, mol + '.ene')
        molform = ob.get_formula(ob.get_mol(mol))

        if io.check_file(enefile):
            E = float(io.read_file(enefile).strip())
            print 'Energy ' + str(E) +  ' pulled from: ' + enefile

        elif runE:
            dire = io.db_sp_path(prog, theory, basisset, None, mol, prog, theory, basisset)
            enefile = io.join_path(dire, mol + '.ene')
            E = run_energy(mol, theory, basisset, prog, True)
            io.write_file(str(E), enefile)
            print 'Energy ' + str(E) +  ' saved to: ' + enefile
    else:
        molform = mol
    if is_auto(E):
        E = getenergy_fromlogfile(logfile,theory)
   
    ####BASIS SELECTION#### 
    basprint = 'manually select basis'
    atomlist = get_atomlist(molform)
    basisselection = 0
    if is_auto(basis[0]):
        basis = select_basis(atomlist)
        #if mol_not_smiles:
        #    basis = select_basis(atomlist)
        #else:
        #    basis = select_better_basis(mol,atomlist)
        basisselection += 1
        basprint = 'automatically generate basis'
    elif basis[0] == 'basis.dat':
        basis = io.read_file('basis.dat').split()
        basprint = 'read basis from basis.dat'
    lines =  ('\n___________________________________________________\n\n' +
              'HEAT OF FORMATION FOR: ' + mol + ' (' + molform + ')' +
              '\n      at ' + theory + '/' +  basisset + 
              '\n\n___________________________________________________\n\n'
              + 'Electronic Energy is: ' +  str(E) + '\nYou have chosen to ' + 
              basprint + '\n\nBasis is: ' + ', '.join(basis))
    print lines 
 
    for bas in basis:
         bas = ob.get_formula(ob.get_mol(bas))
         atomlist.extend(get_atomlist(   bas))
    ####################################

    #COMPUTE Atomlist, stoichlist, matrix, and coefficients
    atomlist = list(set(atomlist))
    stoich = get_stoich(molform,atomlist)
    mat = form_mat(basis,atomlist)

    for i in range(5):
        if np.linalg.det(mat) != 0:
             break
        print 'Matrix is singular -- select new basis'
        atomlist = get_atomlist(molform)
        basis    = select_basis(atomlist,basisselection)
        basisselection += 1

        for bas in basis:
            bas = ob.get_formula(ob.get_mol(bas))
            atomlist.extend(get_atomlist(bas))

        atomlist = list(set(atomlist))
        stoich   = get_stoich(molform,atomlist)
        mat      = form_mat(basis,atomlist)
        print ('\n\nBasis is: ' + ', '.join(basis))
        print mat

    clist =  comp_coeff(mat,stoich)
    ######################################################
     
    ###PRINT STUFF OUT
    lines = '\n  ' + mol + '\t\t' +  '\t'.join(basis) 
    for i in range(len(mat)):
       lines += '\n' + atomlist[i] + '  '
       lines += str(stoich[i]) + '    \t'
       for el in mat[i]:
           lines += str(el) + '\t'
    lines +=  '\n\nCoefficients are: '
    for co in clist:  lines += str(co) + ' '
    print lines + '\n'
    print check(clist, basis, stoich,atomlist)
    ##################

    #COMPUTE AND PRINT delH###
    E =  comp_energy(molform,basis,clist,E,theory,basisset,prog)
    lines =  '\n        delHf(0K)'
    lines += '\nA.U. \t'
    for e in E[:1]:  lines += str(e) + '\t'
    lines += '\nkJ   \t'
    for e in E[:1]:  lines += str(e/ .00038088) + '\t'
    lines += '\nkcal   \t'
    for e in E[:1]:  lines += str(e *  627.503) + '\t'
    ##########################
    
    hf0k = E[0] * 627.503 
    if not mol_not_smiles:
        if not io.check_file(dire + '/' + mol + '.hf0k'):
            io.db_store_sp_prop('Energy (kcal)\tBasis\n----------------------------------',mol,'hf0k',None,prog,theory,basisset)
        s = '\n' + str(hf0k) + '\t' + '  '.join(basis) 
        io.db_append_sp_prop(s,mol,'hf0k',None,prog,theory,basisset)
        
        lines += '\n\nStored heats of formation:\n'
        lines += '----------------------------------\n' 
        lines += io.db_get_sp_prop(mol,'hf0k',None,prog,theory,basisset)
    lines += '\n\n_________________________________________________\n\n'
    print lines

    return hf0k, basis

       
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

    parser.add_argument('-t','--level_of_theory',    type=str,  default='auto/')
    parser.add_argument('-e','--electronic_energy',  type=float, default=9999.9)
    parser.add_argument('-p','--program',            type=str,   default='auto')

    parser.add_argument('-d','--database',           type=str, default='testdb')
    parser.add_argument('-s','--mol_not_smiles',            action='store_true')
    parser.add_argument('-r','--run_energy',                action='store_true')

    ###########################
    args = parser.parse_args()

    mol    = args.molecule
    E      = args.electronic_energy
    runE   = args.run_energy
    basis  = args.select_stoich_basis
    if basis == 'auto':
        basis  = args.select_smiles_basis
    theory = args.level_of_theory
    logfile= args.logfile
    db     = args.database
    prog   = args.program
    mol_not_smiles = args.mol_not_smiles

    main(mol,logfile, basis, E, theory, db, prog, runE, mol_not_smiles) 
