#!usr/bin/python

import re
import numpy as np
import os
import sys
sys.path.insert(0, '/home/elliott/Packages/QTC/')
import iotools as io
import patools as pa

"""
Heatform determines the heat of formation for a molecule by using a basis of 
molecules with well-determined heats of formation

To use this script from command line you can do the following:
(1) python heatform.py -l <logfile.log>
    this will automatically determine what stoichiometry, program, level of theory, basis set,
    and electronic energy heatform needs as parameters 
(2) python heatform.py -s <stoichiometry> -p <program> -t <method>/<basisset> -e <energy>
    manually set all these parameters so that no logfile is required
(3) after (1) or (2)'s arguements, user can also use -b <chosen basis> to manually select the
    basis of molecules
To use this as a module in another script:
(1) main(stoich, logfile, electronic_energy, basisset, theory, database, prog) will return 
    the heat of formation in kcal.  If logfile is set, the rest can be auto (see defaults), and if 
    logfile is gibberish, the rest must be set.
"""
#  TO DO:
#  Option to compute energy for molecule of interest and not just basis species
#  Search directory structures for energies instead of/ alongside using test database
#  Store energies when they are computed for basis setsOnce energy is computed for basis sets

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
        basis.append('NH3')
        i += 1
    if 'S' in atomlist and i<= count:
        basis.append('SO2') 
        i += 1
    if 'H' in atomlist  and i<= count and attempt < 1:
        basis.append('H2')
        i += 1
    elif 'H' in atomlist and 'C' not in atomlist and i<= count and attempt < 2:
        basis.append('H2')
        i += 1
    if 'O' in atomlist and i<= count and attempt < 2:
        basis.append('O2')
        i += 1
    if 'C' in atomlist and i<= count and attempt < 3:
        basis.append('CH4')
        i += 1
    if 'O' in atomlist and 'H' in atomlist and  i<= count and attempt  < 3:
        basis.append('H2O')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count and attempt < 4:
        basis.append('CO2')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count and attempt < 4:
        basis.append('H2CO')
        i += 1
    if 'C' in atomlist and 'O' in atomlist and  i<= count:
        basis.append('CH3OH')
        i += 1
    if 'C' in atomlist and i<= count:
        basis.append('CH3CH3')
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
        mat[i] = get_stoich(mol,atomlist)
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
    
#def update_dic(mol,dic):
#    """
#    Finds if we already have the energy for a molecule in our dictionary, and if not checks if we are in 
#    its EStokTP directory that we can get that information from and adds it to the dictionary
#    INPUT:
#    mol   - molecule
#    dic   - molecule/energy dictionary
#    OUTPUT:
#    dic   - updated molecule/energy dictionary
#    """
#    if mol in dic:
#        return dic
#    elif os.path.exists('geoms/reac1_l1.xyz'):
#        lines = open('geoms/reac1_l1.xyz','r').read()
#        E = float(lines.split('\n')[1])
#        dic[mol] = E
#        return dic
    
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
    return pa.method(lines) 

def getbasisset_fromlogfile(logfile):
    """
    returns basisset used in logfile
    """
    lines = io.read_file(logfile)
    return pa.basisset(lines)

def getenergy_fromlogfile(logfile,theory):
    """
    Returns by default final energy solved for in logfile, or 
    the energy corresponding to a given method in logfile
    """
    lines = io.read_file(logfile)
    prog  = getprog_fromlogfile(lines)
    if prog == 'g09':
       return pa.gaussian_energy(lines)[1]
    elif prog == 'molpro':
       return pa.molpro_energy(lines)[1]
    elif prog == None:
       return pa.molpro_energy(lines)[1]
    print 'no energy found for this molecule, please use -e to manually set it'
    return

def nest_2_dic(dic,key1,key2):
    """
    Returns a dictionary value that requires two key values (is in singly nested loop )
    """
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

def build_gauss(dic, theory, basisset):
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
        import obtools as ob
        mol = ob.get_mol(dic['_id'])
        zmat = ob.get_zmat(mol)
        gauss += str(dic['charge']) + ' ' + str(dic['mult']) + '\n'
    gauss += zmat.lstrip('\n')

    io.write_file(gauss, dic['stoich'] + '.inp')

    return

def build_molpro(dic, theory, basisset):
    """
    Builds a Guassian optimization inputfile for a molecule
    INPUT:
    dic     - dictionary for the molecule
    theory  - theory we want to optimize with
    basisset- basis set we want to optimize with
    OUPUT:
    None (but an inputfile now exists with name <stoich>.inp)
    """
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
    zmat1 = '\n'.join(zmat.split('Variables:')[0].split('\n')[2:] )
    zmat2 = '}\n'
    for line in zmat.split('Variables:')[1].split('\n')[1:-1]:
        zmat2 += line.split()[0] + '  =  ' + line.split()[1] + '\n'
    molp += zmat1 + zmat2
    spin  = 1/2.*(dic['mult']-1) * 2
    molp += '\nSET,SPIN=     ' + str(spin) + '\n\n'
    if spin == 0:
        molp += '!closed shell input\nbasis=' + basisset + '\nhf\n' + theory.lower() + '\noptg\nENERGY=energy'
    else:
        if 'ccsd' in theory.lower() or 'cisd' in theory.lower() or 'hf' in theory.lower() or 'mp' in theory.lower():
	    molp += '!open shell input\nbasis=' + basisset + '\nhf\nu' + theory.lower() + '\noptg\nENERGY=energy'
        else:
            molp += '!open shell input\nbasis=' + basisset + '\nuhf\n' + theory.lower() + '\noptg\nENERGY=energy'
    io.write_file(molp, dic['stoich'] + '.inp')
 
    return 

def run_gauss(filename):
    """
    Runs Gaussian on file: filename
    """
    os.system('soft add +g09; g09 ' + filename)
    
    return

def run_molpro(filename):
    """
    Runs molpro on file: filename
    """
    os.system('/home/elliott/bin/molprop ' + filename)
    
    return

def E_dic(dic,energORs,theory,basisset,prog):
    """
    Checks a dictionary for energy at a specified level of theory and basisset
    Computes energy if it isn't in dictionary already
    INPUT:
    dic      - dictionary we are checking
    energORs - are we checking for the energy or its uncertainty (uncertainty broken)
    theory   - theory energy should be listed or computed with
    basisset - basis set energy should be listed or computed with
    prog     - program an energy should be listed or computed with
    OUTPUT:
    E        - energy 
    """
    if prog.lower() in dic: 
        if theory.lower().lstrip('r') in dic[prog.lower()]:
            if basisset.lower() in dic[prog.lower()][theory.lower().lstrip('r')]:
                return dic[prog.lower()][theory.lower().lstrip('r')][basisset.lower()][energORs]

    if energORs == 'energy' and prog == 'g09':
        print 'Running G09 on ' + dic['stoich'] + ' at ' + theory.lstrip('R').lstrip('U') + '/' + basisset
        build_gauss(dic, theory, basisset)
        run_gauss(dic['stoich']+'.inp')
        E = getenergy_fromlogfile(dic['stoich']+'.log',theory)
        print 'Energy found to be: ' + str(E)
        return E

    elif energORs == 'energy' and prog == 'molpro':
        print 'Running Molpro on ' + dic['stoich'] + ' at ' + theory.lstrip('R').lstrip('U') + '/' + basisset
        build_molpro(dic, theory, basisset)
        run_molpro(dic['stoich']+'.inp')
        E = getenergy_fromlogfile(dic['stoich']+'.log',theory)
        print 'Energy found to be: ' + str(E)
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
    lE        - 0K heat of formation of molecule
    hE        - 298K heat of formation of molecule
    sigma     - uncertainty estimate (i.e., delH = lE +/- sigma)
   
    """
    from testdb import db

    if E == -9999.9:
        for dic in db:
            if dic['stoich'] == mol:
                bas = dic
        E = E_dic(bas, 'energy',theory,basisset,prog)
    lE        = E
    var       = 0
    for i,bas in enumerate(basis):
        for dic in db:
            if dic['stoich'] == bas:
                bas = dic
                break
        lE  +=  coefflist[i] * nest_2_dic(bas,'HeatForm',  0) * 0.00038088

        var += (coefflist[i] * nest_2_dic(bas,'HeatForm','sigma') * 0.00038088   )**2
        
        E    =  E_dic(bas, 'energy',theory,basisset,prog)
        lE  -=  coefflist[i] * E
        #var += (coefflist[i] * E * E_dic(bas,'sigma',theory,basisset,prog) )**2

    sigma = np.sqrt(var)
    return lE, sigma
    
def check(clist, basis,stoich,atomlist):
    """
    Makes sure nothing funky happened while computing coefficients
    """
    check = np.zeros(len(clist))
    statement = 'Coefficients produce correct stoichiometry'
    for i, c in enumerate(clist):
       check += c * get_stoich(basis[i],atomlist)
    for i, sto in enumerate(stoich):
        if not check[i] == sto:
            statement = 'Coefficients do NOT produce correct stoichiometry'
            break
    return statement

def main(mol,logfile='geoms/reac1_l1.log', E=9999.9, basis='auto', theory='auto/',db='tempdb',prog = 'auto'):
    
    basis = basis.split()
    #AUTO SET NONUSERDEFINED PARAMETRS##
    if is_auto(prog):
        prog = pa.get_prog(io.read_file(logfile))
    if is_auto(mol):
        mol = getname_fromdirname() 
    if is_auto(theory) and io.check_file(logfile):
        theory  = gettheory_fromlogfile(logfile)
        theory += '/'
        theory += getbasisset_fromlogfile(logfile)
    theory, basisset = theory.split('/')
    
    if is_auto(E):
         E = getenergy_fromlogfile(logfile,theory)

    basprint = 'manually select basis'
    atomlist = get_atomlist(mol)
    basisselection = 0
    if is_auto(basis[0]):
        basis = select_basis(atomlist)
        basisselection += 1
        basprint = 'automatically generate basis'
    elif basis[0] == 'basis.dat':
        basis = io.read_file('basis.dat').split()
        basprint = 'read basis from basis.dat'
    lines =  ('\n-------------------------------------------------\n\n' +
              'HEAT OF FORMATION FOR: ' + mol +
              '\n    at ' + theory + '/' +  basisset + 
              '\n\n-------------------------------------------------\n\nYou have chosen to ' + 
              basprint + '\n\nBasis is: ' + ', '.join(basis))
    print lines 
 
    for bas in basis:
         atomlist.extend(get_atomlist(bas))
    ####################################

    #COMPUTE Atomlist, stoichlist, matrix, and coefficients
    atomlist = list(set(atomlist))
    stoich = get_stoich(mol,atomlist)
    mat = form_mat(basis,atomlist)

    for i in range(5):
        if np.linalg.det(mat) != 0:
             break
        print 'Matrix is singular -- select new basis'
        atomlist = get_atomlist(mol)
        basis = select_basis(atomlist,basisselection)
        basisselection += 1
        print ('\n\nBasis is: ' + ', '.join(basis))
        for bas in basis:
            atomlist.extend(get_atomlist(bas))
        atomlist = list(set(atomlist))
        stoich = get_stoich(mol,atomlist)
        mat = form_mat(basis,atomlist)
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
    print check(clist, basis,stoich,atomlist)
    ##################

    #COMPUTE AND PRINT delH###
    E =  comp_energy(mol,basis,clist,E,theory,basisset,prog)
    lines =  '\n        delHf(0K)'
    lines += '\nA.U. \t'
    for e in E[:1]:  lines += str(e) + '\t'
    lines += '\nkJ   \t'
    for e in E[:1]:  lines += str(e/ .00038088) + '\t'
    lines += '\nkcal   \t'
    for e in E[:1]:  lines += str(e *  627.503) + '\t'
    lines += '\n\n-------------------------------------------------\n\n'
    print lines
    ##########################
    return E[0] * 627.503

if __name__ == '__main__': 
    """
    Run heat of formation code
    """
    #SET PARAMETERS############
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                            description="SARAH!!! you haven't done this yet!!!")

    parser.add_argument('-s','--stoichiometry',     type=str,  default='auto')
    parser.add_argument('-e','--electronic_energy', type=float,default=9999.9)
    parser.add_argument('-b','--select_basis',      type=str,  default='auto')
    parser.add_argument('-t','--level_of_theory',   type=str,  default='auto/')
    parser.add_argument('-l','--logfile',           type=str,  default='geoms/reac1_l1.log')
    parser.add_argument('-d','--database',          type=str,  default='testdb')
    parser.add_argument('-p','--program',           type=str,  default='g09')

    ###########################
    args = parser.parse_args()

    mol    = args.stoichiometry
    E      = args.electronic_energy
    basis  = args.select_basis
    theory = args.level_of_theory
    logfile= args.logfile
    db     = args.database
    prog   = args.program

    main(mol,logfile, E, basis, theory, db,prog) 
