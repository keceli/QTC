#!/usr/bin/env python
"""
Quantum chemistry tools.
"""
import obtools as ob
import iotools as io
import numpy as np
import patools as pa
import logging
try:
    import cclib
except:
    pass

__updated__ = "2018-03-03"
__authors__ = 'Murat Keceli, Sarah Elliott'


def sort_species_list(slist, printinfo=False, byMass=False):
    """
    Sorts a species list of smiles by number of rotors, electrons and atoms. 
    Optionally, prints info on the list
    """
    tmplist= []
    for s in slist:
        # s = get_slabel(s)
        isomers = ob.get_isomers(s)
        if len(isomers) > 1:
            logging.info('{} isomers found for {} : {}'.format(len(isomers),s,isomers))
        for isomer in isomers:     
            mol = ob.get_mol(isomer,make3D=True)
            nrotor = ob.get_nrotor(mol)
            nelec = ob.get_nelectron(mol)
            natom = ob.get_natom(mol)
            nheavy = ob.get_natom_heavy(mol)
            formula = ob.get_formula(mol)
            smult  = ob.get_multiplicity(isomer)
            obmult = ob.get_multiplicity(mol)
            mass   = ob.get_mass(mol)
            tmplist.append([isomer,formula,smult,obmult,nrotor,nelec,natom,nheavy,mass])
    if byMass:
        tmplist = sorted(tmplist,reverse=False,key=lambda x: (x[8],x[4],x[5],x[6]))
    else:
        tmplist = sorted(tmplist,reverse=True,key=lambda x: (x[4],x[5],x[6]))

    sortedlist = [x[0] for x in tmplist]
    if printinfo:
        logging.info('-'*100)
        if byMass:
            logging.info('{:>8s}\t{:30s} {:20s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}   {:>8s}'.format('Index', 'SMILES', 'Formula', 'Mult', 'OBMult', 'N_rot', 'N_elec', 'N_atom', 'N_heavy','Mass'))
        else:
            logging.info('{:>8s}\t{:30s} {:20s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}'.format('Index', 'SMILES', 'Formula', 'Mult', 'OBMult', 'N_rot', 'N_elec', 'N_atom', 'N_heavy'))
        i = 0
        for tmp in tmplist:
            i += 1
            if byMass:
                logging.info('{:8d}\t{:30s} {:20s} {:8d} {:8d} {:8d} {:8d} {:8d} {:8d}   {:6f}'.format(i,ob.get_smiles(tmp[0]),*tmp[1:]))
            else:
                logging.info('{:8d}\t{:30s} {:20s} {:8d} {:8d} {:8d} {:8d} {:8d} {:8d}'.format(i,ob.get_smiles(tmp[0]),*tmp[1:]))
        logging.info('-'*100)
    sortedlist = remove_dups(sortedlist)
    return sortedlist


def add_species_info(s, parameters):
    """
    Add info:
            nelec
            natom
            nheavy
            formula
            mult
            nrotor (from x2z excluding methyl rotors)
     for a species obtained by openbabel and x2z into a given dictionary that contain 'xyzdir' for
     the location of xyz file.
    Note:
    IO operations
    pwd
    mkdir
    cd
    file read/write
    """
    pwd = io.pwd()
    if parameters['xyzdir']:
        io.mkdir(parameters['xyzdir'])
    else:
        io.mkdir('xyz')
    io.cd(parameters['xyzdir'])
    parameters['symm'] = ob.get_symm(s) 
    parameters['hof'] = ob.get_smileshof(s)
    s = get_slabel(s)
    xyzfile = ob.get_smiles_filename(s) + '.xyz'
    xyzfile = io.fix_path(xyzfile)
    if io.check_file(xyzfile):
        xyz = io.read_file(xyzfile)
        mol = ob.get_mol(xyz)
    else:
        mol = ob.get_mol(s,make3D=True)
        xyz = ob.get_xyz(mol)       
        io.write_file(xyz, xyzfile)
    parameters['nelec']   = ob.get_nelectron(mol)
    parameters['natom']   = ob.get_natom(mol)
    parameters['nheavy']  = ob.get_natom_heavy(mol)
    parameters['formula'] = ob.get_formula(mol)
    parameters['mult']    = ob.get_multiplicity(s) 
    parameters['charge']  = ob.get_charge(mol)
    parameters['xyz'] = xyz
    parameters['nallrotor']  = 0
    parameters['nmethyl'] = 0
    parameters['nrotor'] = 0
    if parameters['natom'] > 1 and io.check_exe(parameters['x2z']):
        try:
            x2z_out = run_x2z(xyzfile, parameters['x2z'])
            nrotor = get_x2z_nrotor(x2z_out)
            nmethyl = get_x2z_nmethyl(x2z_out)
            x2zinfo = get_x2z_info(x2z_out)
            parameters['moltype'] = x2zinfo['moltype']
            parameters['nallrotor']  = nrotor 
            parameters['nrotor']  = nrotor - nmethyl
            parameters['nmethyl'] = nmethyl
        except:
            logging.error('x2z failed for {}'.format(s))
            parameters['nrotor'] = ob.get_nrotor(mol)
    
    else:
        logging.info('x2z not found, using number of rotors computed by open babel') 
        parameters['nrotor'] = ob.get_nrotor(mol)
    io.cd(pwd)
    return parameters  


def get_species_info_list(slist, parameters):
    """
    Return a list of dictionaries for a given species list with info obtained by
    openbabel and x2z.
    Note:
    IO operations
    pwd
    mkdir
    cd
    file read/write
    """
    infolist = []
    for i,s in enumerate(slist):
        d = {}
        pwd = io.pwd()
        io.mkdir(parameters['xyzdir'])
        io.cd(parameters['xyzdir'])
        xyzfile = s + '.xyz'
        xyzfile = io.fix_path(xyzfile)
        if io.check_file(xyzfile):
            xyz = io.read_file(xyzfile)
            mol = ob.get_mol(xyz)
        else:
            mol = ob.get_mol(s,make3D=True)
            xyz = ob.get_xyz(mol)       
            io.write_file(xyz, xyzfile)
        d['mol_id']  = i + 1
        d['nelec']   = ob.get_nelectron(mol)
        d['natom']   = ob.get_natom(mol)
        d['nheavy']  = ob.get_natom_heavy(mol)
        d['formula'] = ob.get_formula(mol)
        d['mult']    = ob.get_multiplicity(s) 
        if d['natom'] > 1 and io.check_exe(parameters['x2z']):
            try:
                x2z_out = run_x2z(xyzfile, parameters['x2z'])
                nrotor = get_x2z_nrotor(x2z_out)
                nmethyl = get_x2z_nmethyl(x2z_out)
                d['nrotor']  = nrotor - nmethyl
                d['nmethyl'] = nmethyl
            except:
                logging.error('x2z failed for {}'.format(s))
                d['nrotor'] = ob.get_nrotor(mol)
        else:
            d['nrotor'] = ob.get_nrotor(mol)
        infolist.append(d)
        io.cd(pwd)
    return infolist  
    
            
def get_input(x, template, parameters):
    """
    Returns input file text for a qc calculation based on a given template.
    """
    if 'xyz' in parameters['results']:
        xyz = parameters['results']['xyz']
    else:
        xyz = parameters['xyz']
    mol = ob.get_mol(xyz)
    mult = parameters['mult']
    nopen = int(mult) - 1
    charge = parameters['charge']
    formula = parameters['formula']
    natom = parameters['natom']
    geo = '\n'.join(xyz.splitlines()[2:natom+2]) 
    #zmat = ob.get_zmat(mol)
    x2zout = run_x2z(xyz, parameters['x2z'])
    if len(x2zout.splitlines()) > 6:
        zmat =  get_x2z_zmat(x2zout)
    else:
        logging.debug(x2zout)
        zmat = ob.get_zmat(mol)
    uniquename = ob.get_inchi_key(mol, mult)
    smilesname = ob.get_smiles_filename(mol)
    smiles = ob.get_smiles(mol)
    nelectron = ob.get_nelectron(mol)
    package = parameters['qcpackage'] 
    method  = parameters[ 'qcmethod'] 
    basis   = parameters[  'qcbasis']
    slabel  = parameters[   'slabel']
    nrotor  = parameters[   'nrotor']
    tmpdir  = parameters[   'tmpdir']
    xyzpath = parameters[  'xyzpath']
    abcd    = parameters[     'abcd'].split(',')
    extraline = ''
    if len(abcd) == 4:
        a,b,c,d = int(abcd[0]), int(abcd[1]), int(abcd[2]), int(abcd[3])
    else:
        a,b,c,d = 6,2,3,100
    heat    = 0
    if 'results' in parameters.keys():
        results = parameters['results']
        if 'deltaH0' in results.keys():
            heat = results['deltaH0']
        if 'xyz' in results.keys():
            xyz = results['xyz']
    task    = parameters['qctask']
    nproc   = parameters['qcnproc']
    totalmem  = int(io.get_total_memory() * 0.4) # in MB
    coremem   = int(float(totalmem)/nproc)  #in MB
    corememmw   = int(float(totalmem)/8./nproc)  # in MegaWords
    lines = template.splitlines()
    inp = ''
    for line in lines:
        if line.strip().startswith('##'):
            pass # Ignore comments
        else:
            inp += line + '\n'
    if task.startswith('tors') or task.startswith('md'):
        if io.check_file(xyzpath):
            inp = inp.replace("QTC(XYZPATH)",xyzpath)
        else:
            inp = inp.replace("QTC(XYZPATH)", 'false')
        if nrotor == 0:
            inp = inp.replace(      "QTC(NMC)", str(1))
        else:
            inp = inp.replace(      "QTC(NMC)", str(min(a+b*(c**nrotor),d)) )
        inp = inp.replace(  "QTC(REFERENCE)", parameters[  'reference'])
        inp = inp.replace(   "QTC(THERMO)", str(parameters['runthermo']))
        if heat:
            inp = inp.replace(   "QTC(HOF)", str(heat))
        else:
            inp = inp.replace(   "QTC(HOF)", 'false')

        if parameters['anharmonic'] == True:
            inp = inp.replace('QTC(ANHARMLOC)', parameters['optlevel'] + '/' + parameters['freqlevel'])
        else:
            inp = inp.replace('QTC(ANHARMLOC)', 'false')
    else:
        if package == 'nwchem':
            if task == 'opt':
                task = 'optimize'
            if ('cc') in method or method.startswith('tce') or 'mp' in method:
                task = 'tce ' + task
                method = method.replace('tce','')
            else:
                task = '{0} {1}'.format(method, task) 
        elif package == 'gaussian':
            if task == 'opt':
                task = 'opt=(maxcyc=50,internal)'
            elif task == 'optfreq':
                task = 'opt=(maxcyc=50,internal) freq'
            elif task == 'optanharm':
                task = 'opt=(maxcyc=50,internal) freq=(anharm,vibrot,readanharm)'
            elif task == 'freq':
                task = 'freq'
            elif task == 'energy':
                task = ''
            elif task == 'anharm':
                task = 'opt=(maxcyc=50,internal) freq=(anharm,vibrot,readanharm)'
        elif package == 'molpro':
            if "QTC(EXTRA)" in inp:
                if nopen == 0:
                    extraline = '{{multi;closed,{0};occ,{1};wf,{2},1,{3};canonical}}'.format(int(nelectron/2),int(nelectron/2),nelectron,nopen)
                elif nopen == 1:
                    extraline = '{{multi;closed,{0};occ,{1};wf,{2},1,{3};canonical}}'.format(int((nelectron-1)/2),int((nelectron+1)/2),nelectron,nopen)
                elif nopen == 2:
                    extraline = '{{multi;closed,{0};occ,{1};wf,{2},1,{3};canonical}}'.format(int((nelectron-2)/2),int((nelectron+2)/2),nelectron,nopen)
                elif nopen == 3:
                    extraline = '{{multi;closed,{0};occ,{1};wf,{2},1,{3};canonical}}'.format(int((nelectron-3)/2),int((nelectron+3)/2),nelectron,nopen)
                elif nopen == 4:
                    extraline = '{{multi;closed,{0};occ,{1};wf,{2},1,{3};canonical}}'.format(int((nelectron-4)/2),int((nelectron+4)/2),nelectron,nopen)
                else:
                    logging.warning('{} input not implemented for multiplicity {}'.format(method,mult))    
            if method.lower().startswith('ccsd'):
                if nopen > 0:
                    method = 'u'+method
            if method.lower().startswith('mp2'):
                if nopen > 0:
                    method = 'r'+method
            if task.lower().startswith('opt'):
                task = 'optg'
            elif task.lower().startswith('single'):
                task = ''
            elif task.lower().startswith('energy'):
                task = ''
            elif task.lower().startswith('freq'):
                task = '{frequencies;print,hessian}'
        elif package == 'qchem':
            zmat = ob.get_zmat(mol,True)
            if task == 'energy':
                task = 'sp'
    if nopen == 0:
        scftype = 'RHF'
        rhftype = 'RHF'
    else:
        scftype = 'UHF'
        rhftype = 'ROHF'
    inp = inp.replace("QTC(CHARGE)", str(charge))
    inp = inp.replace("QTC(MULTIPLICITY)", str(mult))
    inp = inp.replace("QTC(NOPEN)", str(nopen))
    inp = inp.replace("QTC(UNIQUENAME)", uniquename)
    inp = inp.replace("QTC(SMILESNAME)", smilesname)
    inp = inp.replace("QTC(TMPDIR)", tmpdir)
    if task.startswith('tors') or task.startswith('md'):
        inp = inp.replace("QTC(SMILES)", slabel)
    else:
        inp = inp.replace("QTC(SMILES)", smiles)
    inp = inp.replace("QTC(SLABEL)", slabel)
    inp = inp.replace("QTC(ZMAT)", zmat)
    inp = inp.replace("QTC(GEO)", geo)
    inp = inp.replace("QTC(XYZ)", xyz)
    inp = inp.replace("QTC(FORMULA)", formula)
    inp = inp.replace("QTC(METHOD)", method)
    inp = inp.replace("QTC(BASIS)", basis)
    inp = inp.replace("QTC(TASK)", task)
    inp = inp.replace("QTC(PACKAGE)", package)
    inp = inp.replace("QTC(RHF_OR_UHF)", scftype)
    inp = inp.replace("QTC(RHF_OR_ROHF)", rhftype)
    inp = inp.replace("QTC(NPROC)", str(nproc))   
    inp = inp.replace('QTC(ANHARMLOC)', 'false')
    inp = inp.replace('QTC(EXTRA)',extraline)
    inp = inp.replace('QTC(NODE_MEMORY_MB)', str(totalmem))
    inp = inp.replace('QTC(CORE_MEMORY_MB)', str(coremem))
    inp = inp.replace('QTC(CORE_MEMORY_MW)', str(corememmw))
    lines = inp.splitlines()

    if "QTC(" in inp:
        logging.info(66*'#')
        logging.info("Error in template file: \n" + inp)
        logging.info(66*'#')
    return inp

def fix_qckeyword(keyword):
    """
    Fix qckeyword  to minimize the differences
    based on user input.
    """
    keyword = keyword.lower()
    keyword = keyword.replace('g(d,p)','_gdp_')
    keyword = keyword.replace(' ','')
    keyword = keyword.replace('/adz','/aug-cc-pvdz')
    keyword = keyword.replace('/atz','/aug-cc-pvtz')
    keyword = keyword.replace('/aqz','/aug-cc-pvqz')
    keyword = keyword.replace('/dz','/cc-pvdz')
    keyword = keyword.replace('/tz','/cc-pvtz')
    keyword = keyword.replace('/qz','/cc-pvqz')
    keyword = keyword.replace('ccpvdz/','cc-pvdz/')
    keyword = keyword.replace('ccpvtz/','cc-pvtz/')
    keyword = keyword.replace('ccpvqz/','cc-pvqz/')
    keyword = keyword.replace('/augcc','/aug-cc')
    keyword = keyword.replace('/sto3g','/sto-3g')
    keyword = keyword.replace('631','6-31')
    keyword = keyword.replace('torscan/','torsscan/')
    if 'md' in keyword and not 'torsscan' in keyword:
        keywordstr = []
        keyword =  keyword.split(',')
        for key in keyword:
           if 'md' in key:
               keywordstr.append( key.replace('md','torsscan'))
           keywordstr.append(key)
        keyword = ','.join(keywordstr)
    return keyword

def get_slabels_from_json(j):
    """
    Builds a list of strings that contains slabels for 
    all species in json list.
    Needs to convert from unicode to string.
    """
    nitem = len(j)
    slabels = ['']*nitem
    i = 0
    for d in j :
        mult = d['multiplicity']
        slabels[i] = d['SMILES'].encode('ascii','ignore') + '_m' + str(mult)
        i += 1
    return slabels

def update_smiles_list(slist):
    """
    Replaces each smiles with open-babel canonical smiles
    and adds multiplicity with '_m' suffix.
    Removes all inert species.
    >>> sl = ['[Ne]','C','O-O','[O]O','O[O]']
    >>> print qc.update_smiles_list(sl)
    ['C', 'OO', 'O[O]', 'O[O]']
    """
    newlist = []
    for s in slist:
        symm =  None
        ene  =  None
       # if 'He' in s or 'Ne' in s or 'Ar' in s or 'Kr' in s or 'Xe' in s or 'Rn' in s:
        if 'He' in s or 'Ne' in s or 'Ar' in s or 'Kr' in s or 'Xe' in s or 'Rn' in s:
            logging.info('Inert species {0} is removed from the smiles list'.format(s))
        else:
            if '_e' in s:
                s, ene = s.split('_e')
            if '_s' in s:
                s, symm = s.split('_s')
            if '_m' in s:
                smi, mult = s.split('_m')
            else:
                smi = s
                mult = ob.get_mult(s)
                logging.debug("Multiplicity is set by open babel. {} : {}".format(smi,mult))
            canonical = ob.get_smiles(s)
            if canonical.strip() == smi:
                pass
            else:
                logging.debug('SMILES changed after open babel canonicalization {} --> {}'.format(smi,canonical)) 
            slabel = canonical + '_m' + str(mult)
            if symm:
                slabel += '_s' + symm
            if ene:
                slabel += '_e' + ene
            newlist.append(slabel)
    return newlist


def parse_qckeyword(parameters, calcindex=0):
    """
    Updates parameters based on qckeyword and given calcindex.
    Not pure!
    Specifically it updates the values of the following keys in parameters:
    optlevel
    qcpackage
    tspackage
    tsmethod
    tsbasis
    xyzpath
    qcdirectory
    qcpackage
    qcmethod
    qcbasis
    qctask
    parseqc
    writefiles
    anharmonic
    optlevel
    freqlevel
    enlevel
    """
    keyword = parameters['qckeyword']
    optdir  = parameters['optdir']
    package = 'nwchem'
    calcs   = keyword.split(parameters['task_seperator'])
    currentcalc = calcs[calcindex]
    tokens = currentcalc.split('/')
    task = tokens[0]
    if parameters['natom'] == 1:
        logging.info('Number of atoms = 1')
        logging.info('Task = {}'.format(task))
        optdir = ''
        parameters['optdir'] = ''
        if task.startswith('comp') or task.startswith('ene'):
            pass
        else: 
            logging.info('Replacing Task = {} with Task = energy, since there is only a single atom'.format(task))
            task = 'energy'
            parameters['task'] = task
    elif parameters['nallrotor'] == 0:
        if task.startswith('torsopt'):
            logging.info('Replacing Task = {} with Task = opt, since there are no torsions.'.format(task))
            task = 'opt'
        elif task.startswith('torsscan'):
            logging.info('Replacing Task = {} with Task = optfreq, since there are no torsions.'.format(task))
            task = 'optfreq'
        elif task.startswith('md'):
            logging.info('Replacing Task = {} with Task = optfreq, since there are no torsions.'.format(task))
            task = 'optfreq'
    elif parameters['nallrotor'] < 2:
        if task.startswith('md'):
            logging.info('Replacing Task = {} with Task = torsscan, since there is only one torsion.'.format(task))
            task = 'torsscan'
    if task.startswith('ext') or task.startswith('cbs') or task.startswith('comp'):
        task = 'composite'
        if len(tokens) > 2:
            method = tokens[1]
            parameters['composite formula'] = tokens[2:]
        elif len(tokens) == 2:
            method = 'generic'
            parameters['composite formula'] = tokens[1]
        else:
            logging.info('ERROR! Invalid qckeyword: {0}'.format(tokens))
            return
        qcdirectory = io.fix_path(io.join_path(*[optdir,task,method]))
        package = ''
        basis = ''
    else:
        if len(tokens) == 4:
            method = tokens[1]
            basis = tokens[2]
            package = tokens[3]
        elif len(tokens) == 3:
            method = tokens[1]
            basis = ''
            package = tokens[2]
        elif len(tokens) == 2:
            method = ''
            basis = ''
            package = tokens[1]
        else:
            logging.error('ERROR! Invalid qckeyword: {0}'.format(tokens))
        if task == 'torsopt':
            qcdirectory = io.fix_path(io.join_path(*[task,method,basis,package]))
            parameters['optdir'] = qcdirectory
        elif task == 'optfreq':
            qcdirectory = io.fix_path(io.join_path(*['opt',method,basis,package,'freq',method,basis,package]))
            parameters['optdir'] = io.fix_path(io.join_path(*['opt',method,basis,package]))
            parameters['freqdir'] = qcdirectory
        elif task == 'optanharm':
            qcdirectory = io.fix_path(io.join_path(*['opt',method,basis,package,'anharm',method,basis,package]))
            parameters['optdir'] = io.fix_path(io.join_path(*['opt',method,basis,package]))
            parameters['freqdir'] = qcdirectory
        elif 'opt' in task:
            qcdirectory = io.fix_path(io.join_path(*['opt',method,basis,package]))
            parameters['optdir'] = qcdirectory
        elif task == 'torsscan' or task == 'md':
            qcdirectory = io.fix_path(io.join_path(*[optdir,task,method,basis,package]))
            parameters['freqdir'] = qcdirectory
        elif task.startswith('freq') or task.startswith('harm') or task.startswith('hrm'):
            task = 'freq'
            qcdirectory = io.fix_path(io.join_path(*[optdir,task,method,basis,package]))
            parameters['freqdir'] = qcdirectory
        elif task.startswith('anh') or task.startswith('afre'):
            task = 'anharm'
            qcdirectory = io.fix_path(io.join_path(*[optdir,task,method,basis,package]))
            parameters['anharmdir'] = qcdirectory
            parameters['anharmonic'] = True
        elif task.startswith('sp') or task.startswith('single') or task.startswith('ene'):
            task = 'energy'
            if method == 'given':
                qcdirectory = io.fix_path(io.join_path(*[optdir,task,method,package]))
            else:
                qcdirectory = io.fix_path(io.join_path(*[optdir,task,method,basis,package]))
        else:
            logging.info('ERROR! Invalid qckeyword task: {0}'.format(task))
            return      
    parameters['qcdirectory'] = qcdirectory
    parameters['qcpackage'] = package
    parameters['qcmethod'] = method
    parameters['qcbasis'] = basis.replace('_gdp_','g(d,p)')
    parameters['qctask'] = task
    parameters['parseqc'] = True
    parameters['writefiles'] = True
    return parameters 


def get_slabel(smi,mult=None):
    """
    slabel is a unique smiles string for labeling species in QTC.
    Composed of two parts 'canonical smiles' and 'multiplicity'
    Canical smiles strings are unique only for a given code.
    QTC uses open babel.
    slabel = smi + '_m' + str(mult)
    """
    if '_e' in smi:
        smi, ene = smi.split('_e')
    if '_s' in smi:
        smi, symm = smi.split('_s')
    if '_m' in smi:
        smi, mult = smi.split('_m')
    smi = ob.get_smiles(smi)
    if not mult:
        mult = ob.get_multiplicity(smi)
    return smi + '_m' + str(mult)


def get_qlabel(qckeyword, calcindex):
    """
    Returns a string that can be used as a label for
    a given quantum chemistry calculation with qckeyword specified by the user
    and calculation index. The string is based on the dependencies of the
    corresponding calculation.
    """
    calcs = qckeyword.split(',')
    calc  = calcs[calcindex]
    if calc.startswith('compos'):
        label = ','.join(calcs[0:calcindex+1])
    else:
        labels = []
        for i in range(0,calcindex):
            if calcs[i].startswith('energy') or calcs[i].startswith('compos'):
                pass
            else:
                labels.append(calcs[i])
        labels.append(calcs[calcindex])
        label = ','.join(labels)
    return label


def get_qcdirectory(topdir,task,method,basis,package):
    if task.startswith('opt'):
        qcdir = io.join_path(*['opt',method,basis,package])
    else: 
        qcdir = io.join_path(*[topdir,method,basis,package])
    return qcdir


def get_xyz(out,package=None):
    """
    Return xyz as a string from a qc output.
    """
    xyz = ''
    if package is None:
        package = get_output_package(out)
    if package == 'nwchem':
        xyz = get_nwchem_xyz(out)
    elif package == 'gaussian':
        xyz = pa.gaussian_xyz(out)
    elif package == 'molpro':
        xyz = pa.molpro_xyz(out)
    elif package == 'qchem':
        xyz = pa.qchem_xyz(out)
    return xyz
        
        
def parse_output(s, formula, write=False):
    if type(s) is list:
        lines = s
        s = ''.join(lines)
    elif type(s) is str:
        lines = s.splitlines()
        if len(lines) == 1: # Check if s is a filename
            if io.check_file(s):
                s = io.read_file(s)
                lines = s.splitlines()
    else:
        logging.info("First parameter in parse_output should be a string or a list of strings")
    package = get_output_package(s)
    d = {}
    [method,calculation,xyz,basis] = ['na']*4
    nbasis = 0
    energy = 0.
    energies = {}
    freqs = []
    pfreqs = []
    afreqs = []
    xmat= []
    rotconsts = []
    vibrots  = ''
    rotdists = ''
    hessian= ''
    zpve= 0.0
    azpve = 0.0
    hof0 = 0.
    hof298 = 0.
    parsed = False
    messhindered = None
    RPHtinput = None
    if package == 'nwchem':
        method = get_nwchem_method(s)
        calculation = get_nwchem_calculation(s)
        xyz = get_nwchem_xyz(lines)
        nbasis = get_nwchem_nbasis(s)
        energies = get_nwchem_energies(lines)
        energy = energies[method]
        freqs = get_nwchem_frequencies(lines)
        zpve = get_zpve(freqs)
        if energy:
            parsed = True
    elif package == 'molpro':
        method, energy = pa.molpro_energy(s)
        method = method.replace('\(','(').replace('\)',')')  #will figureout source of this later
        zpve           = pa.molpro_zpve(s)
        xyz            = pa.molpro_xyz(s)
        geo            = pa.molpro_geo(s)
        calculation    = pa.molpro_calc(s)
        basis          = pa.molpro_basisset(s)
        zmat           = pa.molpro_zmat(s)
        hessian        = pa.molpro_hessian(s)
        freqs          = list(pa.molpro_freqs(s))
        if energy:
            parsed = True
    elif package == 'qchem':
        energy = pa.qchem_energy(s)
        method = pa.qchem_method(s)
        zpve           = pa.qchem_zpve(s)
        xyz            = pa.qchem_xyz(s)
        geo            = pa.qchem_geo(s)
        calculation    = pa.qchem_calc(s)
        basis          = pa.qchem_basisset(s)
        freqs          = list(pa.qchem_freqs(s))
        if energy:
            parsed = True
    elif package == 'gaussian':
	    method, energy = pa.gaussian_energy(s)
	    zpve           = pa.gaussian_zpve(s)
	    azpve          = pa.gaussian_anzpve(s)
	    calculation    = pa.gaussian_calc(s)
	    basis          = pa.gaussian_basisset(s)
	    zmat           = pa.gaussian_zmat(s)
	    xyz            = pa.gaussian_xyz(s)
	    geo            = pa.gaussian_geo(s)
	    hessian        = pa.gaussian_hessian(s)
	    rotconsts      = pa.gaussian_rotconsts(s) 
	    vibrots        = pa.gaussian_vibrot(s) 
	    rotdists       = pa.gaussian_rotdists(s) 
	    freqs          = list(pa.gaussian_freqs(s))
	    afreqs         = list(get_gaussian_fundamentals(s)[:,1])
	    if sum(afreqs) > 0:
	        xmat           = get_gaussian_xmatrix(s, get_gaussian_nfreq(s))
	        if type(xmat) == str:
		        xmat = []
	    else:
	        afreqs = []
	    if energy:
	        parsed = True
    elif package == 'mopac':
        method = 'SEMO'
        energy = get_mopac_energy(s)
        deltaH0 = get_mopac_deltaH0(s)
        deltaH298 = get_mopac_deltaH298(s)
        xyz = get_mopac_xyz(s)
        freqs = get_mopac_freq(s)
        zpve = get_mopac_zpe(s)
        if energy:
            parsed = True
    elif package.startswith('tors') or package.startswith('md'):
        #optlevel, method, energy = get_torsscan_info(s)
        outfile = 'geoms/reac1_l1.log'
        outfile2 = 'geom.log' #  For torsopt
        xyzfile = 'geoms/reac1_l1.xyz'
        geofile = 'geom.xyz'
        torsoptfile = 'torsopt.xyz'
        if io.check_file(outfile):
            try:
                out = io.read_file(outfile, aslines=False)
                xyz = pa.gaussian_xyz(out)
                method, energy = pa.gaussian_energy(out)
                freqs = pa.gaussian_freqs(out)
                parsed = True
            except:
                logging.error('parse_output: Cannot parse {}'.format(outfile))
                parsed = False
            if io.check_dir('me_files', 1):
                try:
                    xyznew, freqs, pfreqs, zpve, messhindered, RPHtinput = parse_me_files()
                    try:
                        ob.get_mol(xyznew)
                        xyz = xyznew
                    except:
                        logging.error('parse_output: torsscan optimization did not converge')
                except:
                    logging.error('parse_output: Cannot parse me_files {}'.format(io.get_path('me_files')))
        elif io.check_file(xyzfile):
            xyz = io.read_file(xyzfile)
            energy = float(xyz.splitlines()[1].strip())         
            parsed = True
        elif io.check_file(torsoptfile):
            xyz = io.read_file(torsoptfile).strip()
            energy = float(xyz.splitlines()[1].strip())
            parsed = True
        elif io.check_file(geofile):
            geo = io.read_file(geofile).strip()
            natom = len(geo.splitlines())
            energy = 0.
            xyz = '{}\n{}\n{}'.format(natom,energy,geo)
            parsed = True
        else:
            logging.debug('Error in parsing {}'.format(package))
    if parsed:
        if write:
            fname = formula + '.ene'
            io.write_file(str(energy), fname)
            if xyz:
                fname = formula + '.xyz'
                io.write_file(xyz, fname)
            if zpve:
                fname = formula + '.zpve'
                io.write_file(str(zpve), fname)
            if azpve:
                fname = formula + '.anzpve'
                io.write_file(str(azpve), fname)
            if len(freqs) > 0:
                if any(freq < 0 for freq in freqs):
                    logging.error('Imaginary frequency detected: {}'.format(['{:6.1f}'.format(freq) for freq in freqs]))
                fname = formula + '.hrm'
                io.write_file('\n'.join(str(x) for x in freqs), fname )
            if sum(afreqs) > 0:
                if any(freq < 0 for freq in afreqs):
                    logging.error('Imaginary frequency detected: {}'.format(['{:6.1f}'.format(freq) for freq in afreqs]))
                fname = formula + '.anhrm'
                io.write_file('\n'.join(str(x) for x in afreqs), fname)
        d = {'nbasis':nbasis,
               'energy':energy,
               'xyz':xyz,
               'freqs': [float(x) for x in freqs],
               'afreqs': afreqs,
               'pfreqs': pfreqs,
               'zpve': zpve,
               'azpve': azpve,
               'xmat': xmat,
               'rotconsts': rotconsts,
               'vibrots': vibrots,
               'rotdists': rotdists,
               'hindered potential': messhindered,
               'RPHt input': RPHtinput,
               'Hessian'   : hessian,
               'deltaH0': hof0,
               'deltaH298': hof298}

    return d


def get_output_data(out, package=None):
    """
    Parse the output text "out" and return a dictionary
    with parsed results that can be exported as a json file.
    """
    data = {}
    data['warnings'] = ['']
    if package is None:
        package = get_output_package(out)
    if package == 'nwchem':
        lines = out.splitlines()
        data['method'] = get_nwchem_method(out)
        data['calculation'] = get_nwchem_calculation(out)
        data['xyz'] = get_nwchem_xyz(lines)
        data['numberOfBasisFunctions'] = get_nwchem_nbasis(out)
        energies = get_nwchem_energies(lines)
        data['energy_hartree'] = energies[data['method']]
        freqs = get_nwchem_frequencies(lines)
        if len(freqs) > 0:
            data['harmonicFrequencies_cm-1'] = freqs
            data['zpve'] = get_zpve(data['harmonicFrequencies_cm-1'])
            if any(freq < 0 for freq in freqs) < 0:
                data['warnings'].append('Imaginary frequency dedected')
        if data['energy_hartree']:
            data['succesful'] = True
        else:
            data['succesful'] = False
    elif package == 'molpro':
        method, data['energy_hartree'] = pa.molpro_energy(out)
        data['method'] = method.replace('\(','(').replace('\)',')')  #will figureout source of this later
        data['zpve']   = pa.molpro_zpve(out)
        data['xyz']    = pa.molpro_xyz(out)
        data['calculation'] = pa.molpro_calc(out)
        data['basis']  = pa.molpro_basisset(out)
        data['zmat']   = pa.molpro_zmat(out)
        data['hessian'] = pa.molpro_hessian(out)
        data['harmonicFrequencies_cm-1'] = list(pa.molpro_freqs(out))
        if data['energy_hartree']:
            data['succesful'] = True
        else:
            data['succesful'] = False
    elif package == 'gaussian':
        data['method'], data['energy_hartree'] = pa.gaussian_energy(s)
        data['zpve']           = pa.gaussian_zpve(s)
        data['azpve']          = pa.gaussian_anzpve(s)
        data['calculation'] = pa.gaussian_calc(s)
        data['basis']  = pa.gaussian_basisset(s)
        data['zmat']   = pa.gaussian_zmat(s)
        data['xyz']     = pa.gaussian_xyz(s)
        hessian        = pa.gaussian_hessian(s)
        freqs          = list(pa.gaussian_freqs(s))
        afreqs         = list(get_gaussian_fundamentals(s)[:,1])
        if sum(afreqs) > 0:
            xmat           = get_gaussian_xmatrix(s, get_gaussian_nfreq(s))
            if type(xmat) == str:
                xmat = []
        else:
            afreqs = []
        if energy:
            parsed = True
    elif package == 'mopac':
        data['method'] = 'SEMO'
        energy = get_mopac_energy(s)
        deltaH0 = get_mopac_deltaH0(s)
        deltaH298 = get_mopac_deltaH298(s)
        xyz = get_mopac_xyz(s)
        freqs = get_mopac_freq(s)
        zpve = get_mopac_zpe(s)
        if energy:
            parsed = True
    else:
        data['']
    return data



def get_zpve(freqs):
    """
    Return zpve in au for given freq in 1/cm.
    """
    zpve = 0.0
    for freq in freqs:
        freq = float(freq)
        zpve += freq
    return zpve * 0.5 / 219474.63


def get_listofstrings(array):
    """
    Return a list of strings from a given array of strings
    """
    n = len(array)
    s = ['']*n
    for i in range(n):
        s[i] = '{0}\n'.format(array[i])
    return s


def parse_qclog_cclib(qclog,anharmonic=False):
    xyz = None
    freqs = None
    zpe = None
    deltaH = None
    xmat = None
    afreqs = None
    msg =''
    if io.check_file(qclog, 1):
        s = io.read_file(qclog, aslines=False)
    else:
        msg = 'File not found: "{0}"\n'.format(io.get_path(qclog))
        return msg,xyz,freqs,zpe,deltaH,afreqs,xmat
    if check_output(s):
        if cclib:
            ccdata = parse_cclib(qclog)
            xyz = ccdata.writexyz()
            try:
                freqs = ccdata.vibfreqs
                freqs = get_listofstrings(freqs)
                nfreq = len(freqs)
            except AttributeError:
                pass
            try:
                deltaH = ccdata.enthalpy
            except AttributeError:
                pass
            if anharmonic:
                xmat = ccdata.vibanharms
                afreqs = get_gaussian_fundamentals(s, nfreq)[:,1]
                afreqs = get_listofstrings(afreqs)
    else:
        msg = 'Failed job: "{0}"\n'.format(io.get_path(qclog))

    return msg,xyz,freqs,zpe,deltaH,afreqs,xmat


def parse_cclib(out):
    """
    Returns ccdata object that contains data extracted from
    the 'out' file.
    """
    import cclib
    return cclib.io.ccread(out)

def parse_me_files(path=None):
    """
    Parses files in me_files directory that TorsScan (EStoKTP) generates.
    """
    if path:
        pass
    else:
        path = 'me_files'
    xyz    = ''
    freqs  = []
    pfreqs = []
    zpve   = None
    messhindered = None
    RPHtinput = None
    fname = io.join_path(path,'reac1_ge.me')
    if io.check_file(fname):
        out = io.read_file(fname, aslines=False)
        xyz = get_mess_xyz(out)
    fname = io.join_path(path,'reac1_unpfr.me')
    if io.check_file(fname):
        out = io.read_file(fname, aslines=False)
        freqs = get_mess_frequencies(out)
        
    fname = io.join_path(path,'reac1_fr.me')
    if io.check_file(fname):
        out = io.read_file(fname, aslines=False)
        pfreqs = get_mess_frequencies(out)
        
    fname = io.join_path(path,'reac1_zpe.me')
    if io.check_file(fname):
        out = io.read_file(fname, aslines=False)
        try:
            zpve = float(out)  
        except:
            logging.error('cannot parse zpve for file {0}'.format(fname))             
    fname = io.join_path(path,'reac1_hr.me')
    if io.check_file(fname):
        messhindered = io.read_file(fname, aslines=False)   
    fname = 'RPHt_input_data.dat'
    if io.check_file(fname):
        RPHtinput = io.read_file(fname, aslines=False)

    if freqs == pfreqs:
        pfreqs = []
    return xyz, freqs, pfreqs, zpve, messhindered, RPHtinput
	    
def getcc_enthalpy(out):
    if type(out) is not cclib.parser.data.ccData_optdone_bool:
        if io.check_file(out, 1):
            ccdata = parse_cclib(out)
        else:
            return '{0} not found'.format(out)
    else:
        ccdata = out
    return ccdata.enthalpy

def getcc_entropy(out):
    if type(out) is not cclib.parser.data.ccData_optdone_bool:
        if io.check_file(out, 1):
            ccdata = parse_cclib(out)
        else:
            return '{0} not found'.format(out)
    else:
        ccdata = out
    return ccdata.entropy


def getcc_freeenergy(out):
    if type(out) is not cclib.parser.data.ccData_optdone_bool:
        if io.check_file(out, 1):
            ccdata = parse_cclib(out)
        else:
            return '{0} not found'.format(out)
    else:
        ccdata = out
    return ccdata.freeenergy


def getcc_frequencies(out):
    if type(out) is not cclib.parser.data.ccData_optdone_bool:
        if io.check_file(out, 1):
            ccdata = parse_cclib(out)
        else:
            return '{0} not found'.format(out)
    else:
        ccdata = out
    return ccdata.vibfreqs


def getcc_xyz(out):
    if type(out) is not cclib.parser.data.ccData_optdone_bool:
        if io.check_file(out, 1):
            ccdata = parse_cclib(out)
        else:
            return '{0} not found'.format(out)
    else:
        ccdata = out
    return ccdata.writexyz()


def get_periodic_table():
    """
    Return the periodic table as a list.
    Includes elements with atomic number less than 55.
    >>>pt = get_periodic_table()
    >>>print(len(pt))
    >>>55
    """
    pt = ['X' ,
          'H' ,'He',
          'Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne',
          'Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar'
          'K' ,'Ca','Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
          'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe']
    return pt


def get_symbol(atomno):
    """
    Returns the element symbol for a given atomic number.
    Returns 'X' for atomno=0
    >>>print(get_symbol(1))
    >>>H
    """
    pt = get_periodic_table()
    return pt[atomno]


def get_atomno(symbol):
    """
    Return the atomic number for a given element symbol.
    >>>print(get_atomno('H')
    >>>1
    """
    pt = get_periodic_table()
    symbol = symbol.capitalize()
    return pt.index(symbol)


def run(s, parameters, mult=None, trial=0):
    """
    Runs qc, returns a string specifying the status of the calculation.
    """
    maxtrial = 4
    package = parameters['qcpackage']
    overwrite = parameters['overwrite']
    tempdir = parameters['qctemplate']
    task = parameters['qctask']
    recover = parameters['recover']
    slabel  = parameters['slabel']
    tmpdir = io.fix_path(io.join_path(*[parameters['tmpdir'],slabel]))
    qcnproc  = parameters['qcnproc']
    msg = ''
    if trial > maxtrial:
        logging.error('Maximum number of trials reached')
        return 'Maximum number of trials reached'  
    outfile = parameters['qcoutput']
    inpfile = outfile.replace('out','inp')
    runqc   = True
    if io.check_file(outfile, timeout=1):
        if overwrite:
            msg = 'Overwriting previous calculation "{0}"\n'.format(io.get_path(outfile))
            runqc = True
        else:
            out = io.read_file(outfile)
            if task.startswith('tors'):
                if io.check_file('geom.xyz'):
                    msg = 'Skipping calculation, found "{0}"\n'.format(io.get_path('geom.xyz'))
                    runqc = False               
            elif task.startswith('md'):
                if io.check_file('me_files/reac1_hr.me'):
                    if 'Quantum' in io.read_file('me_files/reac1_hr.me'):
                        msg = 'Skipping calculation, found "{0}"\n'.format(io.get_path('me_files/reac1_hr.me'))
                        runqc = False               
            else:
                if check_output(out):
                    logging.info('Successful calculation found "{0}"\n'.format(io.get_path(outfile)))
                    runqc = False
                else:
                    logging.error('Failed calculation found "{0}"\n'.format(io.get_path(outfile)))
                    if recover:
                        logging.info('Renaming failed output and trying to recover')
                        io.mv(outfile, 'failed_{}_{}'.format(outfile,trial))
                        run(s, parameters, mult=mult, trial=trial+1)
                    else:
                        logging.info('Skipping calculation')
                        runqc = False
    else:
        runqc = True
    if task.startswith('tors'):
        package = 'torsscan'
        templatename = task + '_template' + '.txt'
        templatefile = io.join_path(*[tempdir,templatename])
    elif task.startswith('md'):
        package = 'torsscan'
        templatename =  'md_template.txt'
        templatefile = io.join_path(*[tempdir,templatename])
    else:
        if trial > 0:
            templatename = '{0}_{1}_{2}_template.txt'.format(task,package,trial)
        else:
            templatename = '{0}_{1}_template.txt'.format(task,package)
        templatefile =  io.join_path(*[tempdir,templatename])
        if not io.check_file(templatefile):
            if trial > 0:
                templatename = '{0}_{1}_template.txt'.format(package,trial)
            else:
                templatename = '{0}_template.txt'.format(package)
            templatefile =  io.join_path(*[tempdir,templatename])
    if io.check_file(templatefile) and runqc:           
        tmp = io.read_file(templatefile)
        inptext = get_input(s, tmp, parameters)
    else:
        if runqc:
            logging.warning('Set runqc False. Cannot locate template file {}'.format(templatefile))
            runqc = False
        recover = False
    if runqc:
       
        pwd = io.pwd()
        if task.startswith('md'):
            import shutil
            ddir = pwd.split('md')
            if len(ddir) > 2:
                ddir = 'md'.join(ddir[:-1]) + 'torsscan' + ddir[-1]
            else:
                ddir = ddir[0] + 'torsscan' + ddir[1]
            shutil.rmtree(pwd)
            shutil.copytree(ddir, pwd)
            io.cd(pwd)
        inpfile = io.join_path(*[pwd, inpfile])
        io.write_file(inptext, inpfile)
        if package in ['nwchem', 'molpro', 'mopac', 'gaussian', 'torsscan','torsopt','qchem','md' ]:
            if package.startswith('nwc'):
                io.mkdir(tmpdir)
                if parameters['machinefile']:
                    parameters['qcexe'] = 'mpirun -machinefile {0} nwchem'.format(parameters['machinefile'])
                else:
                    parameters['qcexe'] = 'mpirun -n {0} nwchem'.format(qcnproc)
            elif package.startswith('molp'):
                parameters['qcexe'] = '{0} -n {1} -d {2}'.format(parameters['molpro'],qcnproc,parameters['tmpdir'])
            elif package.startswith('qch'):
                #parameters['qcexe'] = '{0} -np {1}'.format(parameters['qchem'],qcnproc)
                parameters['qcexe'] = '{0}'.format(parameters['qchem'])
            elif package.startswith('mopac'):
                parameters['qcexe'] = parameters['mopac']
                mopacdir = io.get_path(parameters['mopac'],directory=True)
                io.set_env_var('MOPAC_LICENSE',mopacdir)
                ldpath= io.get_env('LD_LIBRARY_PATH')
                ldpath = mopacdir + ':' + ldpath
                io.set_env_var('LD_LIBRARY_PATH',ldpath)
            elif task.startswith('tors') or task.startswith('md'):
                parameters['qcexe'] = parameters['torsscan']
            else:
                parameters['qcexe'] = parameters[package]
        if io.check_file(inpfile, timeout=1):
            if package in  ['nwchem', 'torsscan','torsopt', 'md']:
                command = parameters['qcexe'] + ' ' + inpfile
                if package == 'nwchem':
                    io.mkdir('tmp_nwchem')
                logging.info('Running quantum chemistry calculation with {}'.format(command))
                msg += io.execute(command,stdoutfile=outfile,merge=True)
                if package == 'nwchem':
                    io.rmrf('tmp_nwchem')
            elif package == 'molpro':
                inppath = io.get_path(inpfile)
                if len(inppath) > 255:
                    logging.info('Creating {} since path length > 255 (molpro problem)'.format(tmpdir))
                    io.mkdir(tmpdir)
                    io.symlink(tmpdir,'tmp')
                    io.cp(inpfile,tmpdir)
                    io.cd(tmpdir)
                command = parameters['qcexe'] + ' ' + inpfile + ' -o ' + outfile
                logging.info('Running quantum chemistry calculation with {}'.format(command))
                msg += io.execute(command,stdoutfile='stdouterr.txt',merge=True)
                if len(inppath) > 255:
                    logging.info('Copying {} from {}'.format(outfile,tmpdir))
                    io.cp(outfile,pwd)
                    io.cd(pwd)
            else:
                command = parameters['qcexe'] + ' ' + inpfile + ' ' + outfile
                logging.info('Running quantum chemistry calculation with {}'.format(command))
                msg += io.execute(command)
            outfile2 = inpfile + '.out'
            ### MOPAC may create *.inp.out file depending on the version
            if io.check_file(outfile2, timeout=1):
                io.mv(outfile2,outfile)
            ### 
            if io.check_file(outfile, timeout=1):
                msg += ' Output file: "{0}"\n'.format(io.get_path(outfile))
                out = io.read_file(outfile)
                io.rmrf('tmp')
                if not check_output(out):
                    logging.error('Failed calculation "{0}"\n'.format(io.get_path(outfile)))
                    if recover:
                        logging.info('Attempting to recover, trial {}'.format(trial+1))
                        logging.info('Renamed failed output')
                        io.mv(outfile, 'failed_{}_{}'.format(outfile,trial))
                        run(s, parameters, mult=mult, trial=trial+1)
                    else:
                        logging.info('Skipping calculation')
                        runqc = False
            io.rmrf(tmpdir)
        else:
            msg += 'Failed, cannot find input file "{0}"\n'.format(io.get_path(inpfile))
    return msg


def run_composite(parameters):
    """
    Parses qckeyword for composite method. 
    'opt/mp2/cc-pvdz/gaussian,freq/mp2/cc-pvtz/molpro,sp/mp2/cc-pvqz,composite/cbs-dtq/energy=0.1 * E[0] + 0.4 * E[1] + 0.5 * E[2]'
    """
    qckeyword = parameters['qckeyword'] 
    slabel  = parameters['slabel']
    formula = parameters['formula']
    compositeformula = parameters['composite formula'][0]
    method = parameters['qcmethod']
    allresults = parameters['all results']
    logging.info('Composite energy formula = {0}'.format(compositeformula))
    calcindex = parameters['calcindex']
    e = [0.] * calcindex
    energy = None
    enefile = formula + '.ene'  
    inpfile = formula + '_' + method  + '.inp'  
    for i in range(0,calcindex):
        try:
            qlabel = get_qlabel(qckeyword, i)
            e[i] = allresults[slabel][qlabel]['energy']
        except Exception as err:
            e[i] = 0.
            logging.error('Energy not found for {0} {1}. Exception: {2}'.format(slabel,qlabel,err))
    exec(compositeformula)
    if energy:
        io.write_file(str(energy),enefile )
        io.write_file(compositeformula,inpfile )
        logging.info('Composite energy = {} Hartree\n'.format(energy))
        logging.debug('Energy file: "{}"\n'.format(io.get_path(enefile)))
    else:
        energy = 0.
        logging.error('Problem in given composite formula: {}'.format(compositeformula))
    return energy

def use_given_hof(parameters):
    """
    Prints success of hof read from <smiles>_e<hof>
    """
    method = parameters['qcpackage'] 
    slabel  = parameters['slabel']
    hof = None
    if parameters['hof']:
        hof = parameters['hof'] 
    if hof:
        logging.info('Given HoF = {} kcal\n'.format(hof))
    else:
        hof = 0.
        logging.error('Problem in given energy: {}'.format(hof))
    return hof


def run_qcscript(qcscriptpath, inputpath, geopath, multiplicity):
    """
    Submit jobs using Ahren's script.
     Usage: qcscript.pl input.qcscript initialgeometry.geo multiplicity
     Sample geo file:
     C          0.96737       -0.07578        0.02761
     O          2.33437       -0.07578        0.02761
    """
    from subprocess import Popen, PIPE
    process = Popen(['perl',qcscriptpath, inputpath, geopath, str(multiplicity)], stdout=PIPE, stderr=PIPE)
    msg, err = process.communicate()
    if err:
        msg = 'Failed {0}'.format(err)
    return msg


def run_x2z(xyz,exe='x2z', getinfo=False, outputfile=None):
    """
    Runs Yury's x2z, and returns the output text as a string (default) or 
    if getinfo is True returns a dictionary with parsed information.
    """
    result = ''
    if io.check_file(xyz):
        inp = xyz
    else:
        io.write_file(xyz, 'x2z.xyz')
        inp = 'x2z.xyz'
    logging.debug('Running x2z for {}'.format(inp))
    out, err = io.get_stdout_stderr([exe,inp])
    if err:
        logging.info('STDOUT for x2z:\n {}'.format(out))
        logging.error('STDERR for x2z:\n {}'.format(err))
        try:
            if getinfo:
                x2zinfo = get_x2z_info(out)
                result = x2zinfo
            else:
                result = out
            if outputfile:
                io.write_file(out, outputfile)
            return result
        except:
            return ''
    else:
        if getinfo:
            x2zinfo = get_x2z_info(out)
            result = x2zinfo
        else:
            result = out
        if outputfile:
            io.write_file(out, outputfile)
    return result


def get_x2z_nmethyl(out):
    """
    Return the number of methyl rotors based on x2z output
    """
    if io.check_file(out):
        out = io.read_file(out)
    out = out.lower()
    lines = out.splitlines()
    nrotor = 0
    key = 'rotational groups:'
    n = io.get_line_number(key,lines=lines)
    nmethyl = 0
    for line in lines[n+1:]:
        if line:
            groups = line[3:].split()
            if 'c1h3' in groups:
                nmethyl += 1
        else:
            break
    return nmethyl


def get_x2z_nrotor(out):
    """
    Return the number of rotors based on x2z output
    """
    if io.check_file(out):
        out = io.read_file(out)
    out = out.lower()
    lines = out.splitlines()
    nrotor = 0
    key = 'rotational bond dihedral angles:'
    for line in lines:
        if line.startswith(key):
            line = line.replace(key,'').strip()
            if line:
                if ',' in line:
                    nrotor = len(line.split(','))
                else:
                    nrotor = 1
            else:
                nrotor = 0
    return nrotor


def get_x2z_sym(out):
    """
    Return the rotational symmetry number.
    molecule is nonlinear
has enantiomer? yes

rotational symmetry number = 1

Molecular structure: resonantly stabilized (2 resonances)

   A\A    O1   O2   H3   C4   O5   O6 

   O1     X

   O2     1    X

   H3     1    0    X

   C4     0    1    0    X

   O5     0    0    0  1.5    X

   O6     0    0    0  1.5    0    X


Free radical: Radical sites are O4 O5 

Z-matrix atom order:
 0 -->  0
 1 -->  1
 2 -->  2
 3 -->  3
 4 -->  4
 5 -->  5

Z-Matrix:
O
O, 1 , R1 
H, 1 , R2 , 2 , A2 
C, 2 , R3 , 1 , A3 , 3 , D3 
O, 4 , R4 , 2 , A4 , 1 , D4 
O, 4 , R5 , 2 , A5 , 5 , D5 

 R1   =         2.77828
 R2   =         1.82926 A2   =         95.1267
 R3   =         2.58174 A3   =          110.52 D3   =         74.8222
 R4   =          2.3984 A4   =         119.574 D4   =         64.0599
 R5   =         2.39823 A5   =         114.245 D5   =             180

Rotational bond dihedral angles: d3, d4

Beta-scission bonds: r3

Linear bonds: 
    """
    if io.check_file(out):
        out = io.read_file(out)
    out = out.lower()
    lines = out.splitlines()
    sym = 1
    key = 'rotational symmetry number ='
    for line in lines:
        if line.startswith(key):
            try:
                sym = int(line.split()[-1])
            except:
                logging.error('Failed in parsing sym. number in x2z output {}'.format(out))
    return sym

def get_x2z_zmat(out):
    if io.check_file(out):
        out = io.read_file(out)
    zmat = ''
    atoms = []
    measure = []
    if out:
        props, lines = out.lower().split('z-matrix:\n')
        #zmatrix connectivity
        lines = lines.split('\n')
        for i,line in enumerate(lines):
            if line == '':
                break
            atoms.append('  '.join(line.rstrip('\n').replace(' ','').split(',')))     #Main part of ZMAT
        zmat += '\n'.join(atoms)
        #zmatrix parameters
        for j in range(i+1,len(lines)):
            if lines[j]:
                if 'const' in lines[j].lower() or 'rot' in lines[j] or 'beta' in lines[j]:
                    break
            measure.extend(lines[j].replace('=',' ').rstrip('\n').split())   #Gets parameters
        measure = np.array(measure)                     
        if  (len(measure)%2 != 0):
            measure = measure[:-1]
        measure = measure.reshape( len(measure)/2, 2)      #Puts measurements into two columns
        for angle in measure:
            if 'r' in angle[0].lower():
                angle[1] = str(float(angle[1]) * 0.529177) #bohr to angstrom 
            angle[1] = '{:.4f}'.format(float(angle[1]))
        measure = measure.tolist()
        for i, meas in enumerate(measure): 
            measure[i] = ' = '.join(measure[i])    
        zmat += '\nVariables:\n' + '\n'.join(measure) 
    return zmat

def get_x2z_bonds(out):
    """
    Gets bonds from x2z connectivity matrix
    """
    if io.check_file(out):
        out = io.read_file(out)
    out = out.lower().replace('\n\n','\n')
    lines = out.splitlines()
    mat = []
    startkey = 'molecular structure:'
    endkey = 'z-matrix atom order:'
    save = False
    for line in lines:
        if line.startswith(endkey):
            save = False
        if save:
            mat.append(line)
        if line.startswith(startkey):
            save = True
    atoms = get_atom_list(mat[0])
    
    bonds = {}
    for i, line in enumerate(mat[1:]):
        if line == '':
            break
        for j, typ in enumerate(line.split()[1:-1]):
            if float(typ) > 0.0:
                bonded = sorted([atoms[i], atoms[j]])
                bond = '{0}-{2}({1})'.format(bonded[0],typ,bonded[1])
                if bond in bonds:
                    bonds[bond] += 1
                else:
                    bonds[bond]  = 1
    return bonds

def get_atom_list(line):
    """
    splits the number from the atom 
    a\\a   o1   c2   c3   c4   o5   h6   h7   h8   h9
    -> [O, C, C, C, O, H, H, H, H]
    """
    line = line.split()[1:]
    atoms =  []
    for a in line:
        atom = ''
        num  = ''
        for b in a:
           if b.isdigit():
               num += b
           else:
               atom += b
        atoms.append(atom.capitalize())
    return atoms
 
def get_atom_stoich(line):
    line += 'Z'
    atoms = {}
    atom = ''
    num  = ''
    newatom = False
    for b in line:
       if b.isdigit():
           num += b
       else:
           if atom and b.lower() != b:
               newatom = True
           if newatom:
               if num:
                   atoms[atom] = float(num)
               else:
                   atoms[atom] = 1.
               num = ''
               atom = ''
               newatom = False
           atom += b
    return atoms

def get_x2z_info(out):
    """
    Returns a dictionary that contains information based on
    given x2z output string or output filename . 
    """
    if io.check_file(out):
        out = io.read_file(out)
    out = out.lower()
    lines = out.splitlines()
    x2zinfo = {}
    x2zinfo['moltype'] = 'unknown'
    x2zinfo['nbeta'] = 0
    x2zinfo['nlinear'] = 0
    x2zinfo['nradical'] = 0
    x2zinfo['nsym'] = -1
    x2zinfo['nrotor'] = 0
    x2zinfo['nresonance'] = 0
    x2zinfo['bonds'] = get_x2z_bonds(out)
    x2zinfo['zmat'] = get_x2z_zmat(out)
    x2zinfo['nmethyl'] = get_x2z_nmethyl(out)
    x2zinfo['bonds'] = get_x2z_bonds(out)
    for line in lines:
        if 'molecule is' in line:
            x2zinfo['moltype'] = line.replace('molecule is', '').strip()
        elif 'beta-scission' in line:
            line = line.replace('beta-scission bonds:','')
            x2zinfo['nbeta'] = len(line.split(','))
        elif 'linear bonds' in line:
            line = line.replace('linear bonds:','')
            x2zinfo['nlinear'] = len(line.split(','))
        elif 'free radical' in line:
            """
            Free radical: Radical sites are C2 C4
            Free radical: Radical site is C2
            """
            line = line.split('radical')[-1]
            x2zinfo['nradical'] = len(line.split()) - 2
        elif 'rotational symmetry number =' in line:
            x2zinfo['nsym'] = int(line.split()[-1])
        elif 'rotational bond dihedral angles:' in line:
            line = line.replace('rotational bond dihedral angles:','')
            if line.strip():
                x2zinfo['nrotor'] = len(line.split(','))
        elif 'resonantly stabilized' in line:
            line = line.replace('molecular structure:','').replace('(',' ').replace(')',' ')
            for item in line.split():
                if item.isdigit():
                    x2zinfo['nresonance'] = int(item)
    return x2zinfo                


def get_ob_info(s):
    """
    Returns a dictionary containing open-babel generated info for
    a given chemical identifier (smiles, inchi, etc).
    """
    mol = ob.get_mol(s)
    obinfo = {}
    obinfo['nobrotor'] = ob.get_nrotor(mol)
    obinfo['nelec'] = ob.get_nelectron(mol)
    obinfo['natom'] = ob.get_natom(mol)
    obinfo['nheavy'] = ob.get_natom_heavy(mol)
    obinfo['formula'] = ob.get_formula(mol)
    obinfo['smult']  = ob.get_multiplicity(s)
    obinfo['obmult'] = ob.get_multiplicity(mol)
    return obinfo


def get_list_info(slist,sep=' , '):
    info = ''
    header = 'Index' + sep + 'Given' + sep + 'Can.SMILES' + sep + 'XYZ_SMILES'+ sep
    for i,s in enumerate(slist):
        obinfo = get_ob_info(s)
        xyz = ob.get_xyz(s)
        cansmi = ob.get_smiles(s)
        xyzsmi = ob.get_smiles(xyz)
        x2zinfo = run_x2z(xyz,getinfo=True)
        info += str(i+1) + sep + s     + sep + cansmi       + sep + xyzsmi + sep
        for key,value in obinfo.iteritems():
            if i == 0:
                header += key + sep
            info += str(value) + sep
        for key,value in x2zinfo.iteritems():
            if i == 0:
                header += key + sep
            info += str(value) + sep
        if i == 0:
            header += '\n'
        info += '\n'
    final = header + info
    return final


def check_output(s):
    """
    Returns true/false if quantum chemistry calculation completed/failed.
    """
    if "Normal termination of Gaussian" in s:
        completed = True
    elif "== MOPAC DONE ==" in s:
        completed = True
    elif "Kowalski" in s:
        completed = True
    elif "Variable memory released" in s:
        completed = True
    elif "Thank you very much for using Q-Chem" in s:
        completed = True
    elif io.check_file('me_files/reac1_fr.me'):#TorsScan/ES2KTP
        completed = True
    elif io.check_file('geom.xyz'):#Torsopt/ES2KTP
        completed = True
    else:
        completed = False
    return completed


def find_xyzfile(xyzpath,smilesdir):
    """
    Returns the path for xyzfile.
    """
    xyzfile = ''
    if io.check_file(xyzpath):
        xyzfile = xyzpath
    elif io.check_file(io.join_path(*(smilesdir,xyzpath))):
        xyzfile = io.join_path(*(smilesdir,xyzpath))
    elif io.check_dir(xyzpath):
        try:
            xyzfile = io.find_files(xyzpath, '*.xyz')[0]
        except:
            pass
    elif xyzpath and io.check_dir(io.join_path(*(smilesdir,xyzpath))):
        xyzpath = io.join_path(*(smilesdir,xyzpath))
        try:
            xyzfile = io.find_files(xyzpath, '*.xyz')[0]
        except:
            pass
    return xyzfile 

  
def get_output_package(out,filename=False):
    """
    Returns the name of qc package if the calculation is successful.
    Returns None if failed or unknown package.
    """
    if filename:
        out = io.read_file(out,aslines=False)
    if "Gaussian(R)" in out:
        p = 'gaussian'
    if "Normal termination of Gaussian" in out:
        p = 'gaussian'
    elif "== MOPAC DONE ==" in out:
        p = 'mopac'
    elif "Northwest Computational Chemistry Package" in out:
        if "Straatsma" in out:
            p = 'nwchem'
        else:
            p = 'failed_nwchem'
    elif "PROGRAM SYSTEM MOLPRO" in out:
        if "Variable memory released" in out:
            p = 'molpro'
        else:
            p = 'failed_molpro'
    elif "Task: Submitting EStokTP job" in out:
        p = 'torsscan'
    elif "Beta-scission bonds:" in out:
        p = 'x2z'
    elif "thermo properties for species" in out:
        p = 'thermp'
    elif "Thank you very much for using Q-Chem" in out:
        p = 'qchem'
    else:
        p = None
    return p


def get_output_status(out,filename=False):
    """
    Returns the name of qc package if the calculation is successful.
    Returns None if failed or unknown package.
    """
    p = 'NA'
    s = 'NA'
    m = ''
    nchar = 512
    if filename:
        out = io.read_file(out,aslines=False)
    if "Gaussian(R)" in out:
        p = 'gaussian'
        if "Normal termination of Gaussian" in out:
            s = 'OK'
        else:
            s = 'FAILED'
            m = out[-nchar:-1]
    elif "MOPAC" in out:
        p = 'mopac'
        if "== MOPAC DONE ==" in out:
            s = 'OK'
        else:
            s = 'FAILED'
            m = out[-nchar:-1]
    elif "Northwest Computational Chemistry Package" in out:
        p = 'nwchem'
        if "Straatsma" in out:
            s = 'OK'
        else:
            s = 'FAILED'
            m = out[-nchar:-1]
    elif "PROGRAM SYSTEM MOLPRO" in out:
        p = 'molpro'
        if "Variable memory released" in out:
            s = 'OK'
        else:
            s = 'FAILED'
            m = out[-nchar:-1]
    elif "Task: Submitting EStokTP job" in out:
        p = 'torsscan'
    elif "Beta-scission bonds:" in out:
        p = 'x2z'
    elif "thermo properties for species" in out:
        p = 'thermp'
    else:
        p = None
    return p, s, m



def get_package(templatename):
    """
    Return the quantum chemistry package name based on the template file name
    """
    suffix = templatename.split('.')[-1]
    if 'cbs' in templatename or suffix is 'py':
        p = 'extrapolation'
    elif 'nwchem' in templatename or 'nw' in suffix:
        p = 'nwchem'
    elif 'qchem' in templatename or 'qc' in suffix:
        p = 'qchem'
    elif 'gau' in templatename or 'g09' in suffix or 'com' in suffix or 'g03' in suffix:
        p = 'gaussian'
    elif 'mopac' in templatename or 'mop' in suffix:
        p = 'mopac'
    elif 'molpro' in templatename or 'mlp' in suffix:
        p = 'molpro'
    else:
        p = 'unknown'
    return p

    
def get_gaussian_input(x, template, mult=0):
    """
    Returns Gaussian input file based on a given template.
    """
    if type(x) is str:
        mol = ob.get_mol(x)
    else:
        mol = x
    if mult == 0:
        mult = ob.get_multiplicity(mol)
    charge = ob.get_charge(mol)
    geo = ob.get_geo(mol)
    xyz = ob.get_xyz(mol)
    zmat = ob.get_zmat(mol)
    uniquename = ob.get_inchi_key(mol, mult)
    inp = template.replace("QTC(CHARGE)", str(charge))
    inp = inp.replace("QTC(MULTIPLICITY)", str(mult))
    inp = inp.replace("QTC(UNIQUENAME)", uniquename)
    inp = inp.replace("QTC(ZMAT)", zmat)
    inp = inp.replace("QTC(GEO)", geo)
    inp = inp.replace("QTC(XYZ)", xyz)
    if "QTC(" in inp:
        logging.info("Error in template file:\n" + inp)
        return
    return inp


def get_gaussian_natom(lines):
    """
    NAtoms=     30 NQM=       30 NQMF=       0 NMMI=      0 NMMIF=      0
    """
    import iotools as io
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'NAtoms='
    n = io.get_line_number(keyword, lines=lines)
    return int(lines[n].split()[1])


def get_gaussian_basis(lines):
    """
    Standard basis: CC-pVDZ (5D, 7F)
    """
    import iotools as io
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'Standard basis:'
    n = io.get_line_number(keyword, lines=lines)
    return int(lines[n].split()[2])


def get_method(s):
    """
    """
    methods = ['b3lyp','ccsdt(q)','ccsdt','ccsd(t)','ccsd','mp2','mp3','pm3', 'pm6', 'pm7']
    method = ''
    for m in methods:
        if m in s.lower:
            method = m
            break
    method = method.replace('(','p')
    method = method.replace(')','')
    method = method.replace('*','')

    return method


def get_basis(s):
    """
    """
    opts = ['aug-cc-pvdz','aug-cc-pvtz','aug-cc-pvqz','cc-pvdz','cc-pvtz','cc-pvqz', 'sto-3g', '6-31g']
    basis = ''
    for m in opts:
        if m in s.lower():
            basis = m
            break
    basis = basis.replace('(','p')
    basis = basis.replace(')','')
    basis = basis.replace('*','')
    return basis


def get_gaussian_xyz(lines,optimized=True):
    """
                          Input orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          1           0        0.000000    0.000000    0.122819
      2          1           0        0.000000    0.000000    0.877181
 ---------------------------------------------------------------------
    """
    import iotools as io
    if type(lines) == str:
        lines = lines.splitlines()
    natom = get_gaussian_natom(lines)

    keyword = 'Input orientation:'
    n = io.get_line_number(keyword, lines=lines,getlastone=optimized)
    #for i in range(n):

    return


def get_gaussian_zpve(s):
    """
    Parses zero-point vibrational energy from gaussian
    log file.
    Input:
    s: String containing the log file output.
    Returns:
    If successful:
        Float, zpve in kcal/mol
    else:
        A string showing the error.
    Portion of the relevant output:
    Zero-point vibrational energy     194497.1 (Joules/Mol)
                                   46.48591 (Kcal/Mol)
    """
    key = "Zero-point vibrational energy"
    lines = s.splitlines()
    iline = io.get_line_number(key,lines=lines)
    if iline < 0:
        return 'Not found: {0}'.format(key)
    iline += 1
    line = lines[iline]
    return float(line.split()[0])


def get_gaussian_xmatrix(s,nfreq):
    """
    Parses  X matrix from Gaussian log file.
    Input:
    s: String containing the log file output.
    nfreq : number of vibrational frequencies
    Returns:
    If successful:
        Numpy 2D array of size: nfreq x nfreq
        Only lower half triangle is filled.
        Unit of the elements is cm-1.
    else:
        A string showing the error.
    Portion of the relevant output:
 X matrix of Anharmonic Constants (cm-1)
                1             2             3             4             5
      1       -16.952
      2       -20.275       -16.791
      3       -68.264       -19.750       -17.403
      4       -19.762       -67.602       -20.866       -17.232
      5       -39.101       -38.921       -40.903       -40.724        -9.467
      6       -40.164       -39.960       -39.582       -39.374       -37.734
      7        -3.979        -8.048        -3.791        -7.530         4.360
      8        -8.117        -3.911        -7.606        -3.720         4.345
      9        -5.530        -9.506        -4.745        -9.316        -0.471
     10        -9.552        -5.430        -9.359        -4.641        -0.441
     11        -3.195        -3.073        -0.553        -0.426        11.081
     12        -3.864        -3.746        -5.117        -4.998         2.696
     13        -1.869        -3.005        -1.758        -1.604        -1.937
     14        -3.070        -1.768        -1.671        -1.662        -1.899
     15         1.427         1.429         1.853         1.856         0.340
     16        -2.443        -2.611        -2.489        -3.468        -2.054
     17        -2.798        -2.340        -3.661        -2.385        -2.062
     18         0.189         0.483         1.107         1.390         1.489
                6             7             8             9            10
      6        -9.394
      7        -9.658        -4.966
      8        -9.649        -3.176        -4.942
      9       -10.774        -2.897        -2.061        -3.149
     10       -10.737        -2.054        -2.842        -0.368        -3.127
     11         7.094        -2.539        -2.509        -2.613        -2.570
     12         6.172        -0.754        -0.714        -2.250        -2.202
     13        -1.744        -4.336        -4.722        -3.153        -3.462
     14        -1.707        -4.716        -4.248        -3.466        -3.078
     15         0.219        -1.347        -1.343        -1.539        -1.537
     16        -1.895        -1.989        -2.795        -5.621        -6.255
     17        -1.915        -2.837        -1.935        -6.357        -5.550
     18         0.433        -1.192        -1.110        -0.335        -0.089
               11            12            13            14            15
     11        -7.155
     12       -18.994        -3.594
     13        -5.621        -3.171        -1.575
     14        -5.555        -3.100         0.450        -1.551
     15        -7.109        -4.660        -4.614        -4.605        -5.489
     16        -1.602        -0.998        -4.504         1.095        -2.827
     17        -1.560        -0.941         1.012        -4.392        -2.819
     18         1.778        -0.412        -5.035        -4.808         1.097
               16            17            18
     16         2.460
     17         7.634         2.481
     18         8.273         8.287       -13.180

 Resonance Analysis
    """
    xmat = np.zeros((nfreq,nfreq))
    lines = s.splitlines()
    key = 'X matrix of Anharmonic Constants (cm-1)'
    key2 = 'Total Anharmonic X Matrix (in cm^-1)'
    iline = io.get_line_number(key,lines=lines)
    if iline > 0:
        iline += 1
    else:
        iline = io.get_line_number(key2,lines=lines)
        iline += 2
    line = lines[iline]
    if iline < 3:
        return 'Not found: {0}'.format(key)
    try:
        while line.strip():
            cols = line.split()
            icol = int(cols[0])-1
            for irow in range(icol,nfreq):
                iline += 1
                line = lines[iline]
                cols = line.split()
                ncol = len(cols) - 1
                xmat[irow,icol:icol+ncol] = [float(num.replace('D','E')) for num in cols[1:]]
            iline += 1
            line = lines[iline]
        xmat = xmat.tolist()
    except:
            logging.warning('Ignoring anharmonicities -- unexpected length of xmat')
            xmat  = []
    return xmat


def get_gaussian_fundamentals(s,nfreq=None):
    """
    Parses harmonic and anharmonic frequencies from gaussian
    log file.
    Input:
    s: String containing the log file output.
    nfreq : number of vibrational frequencies
    Returns:
    If successful:
        Numpy 2D array of size: nfreq x 2
        1st column for harmonic frequencies in cm-1
        2nd column for anharmonic frequencies in cm-1
    else:
        A string showing the error.
    Portion of the relevant output:
  Fundamental Bands (DE w.r.t. Ground State)
  1(1)           3106.899     2957.812   -0.042978   -0.008787   -0.008920
  2(1)           3106.845     2959.244   -0.042969   -0.008924   -0.008782
  3(1)           3082.636     2934.252   -0.043109   -0.008543   -0.008705
  4(1)           3082.581     2935.702   -0.043101   -0.008709   -0.008539
  5(1)           3028.430     2918.529   -0.048859   -0.008796   -0.008794
  6(1)           3026.064     2926.301   -0.048438   -0.008788   -0.008785
  7(1)           1477.085     1438.911   -0.044573   -0.001097   -0.007855
  8(1)           1477.063     1439.122   -0.044576   -0.007858   -0.001089
  9(1)           1474.346     1432.546   -0.043241    0.000678   -0.007062
 10(1)           1474.318     1432.981   -0.043245   -0.007065    0.000691
 11(1)           1410.843     1377.548   -0.028060   -0.016937   -0.016944
 12(1)           1387.532     1356.818   -0.027083   -0.016001   -0.016001
 13(1)           1205.022     1177.335   -0.029813   -0.010333   -0.011188
 14(1)           1204.977     1177.775   -0.029806   -0.011191   -0.010328
 15(1)           1011.453      988.386   -0.037241   -0.014274   -0.014270
 16(1)            821.858      814.503   -0.025712   -0.008603   -0.010446
 17(1)            821.847      814.500   -0.025693   -0.010449   -0.008599
 18(1)            317.554      296.967   -0.035184   -0.010866   -0.010861
Overtones (DE w.r.t. Ground State)
    """
    if nfreq is None:
        nfreq = get_gaussian_nfreq(s)
    freqs = np.zeros((nfreq,2))
    lines = s.splitlines()
    key = 'Fundamental Bands (DE w.r.t. Ground State)'
    iline = io.get_line_number(key,lines=lines)
    if iline > 0:
        for i in range(nfreq):
            iline += 1
            line = lines[iline]
            cols = line.split()
            freqs[i,:] = [float(cols[-5]),float(cols[-4])]
    return freqs[freqs[:,0].argsort()]


def get_mopac_input(x, method='pm3', keys='precise nosym threads=1 opt', mult=1, dothermo=False):
    """
    Returns mopac input as a string.
    Note: For doctest I had to escape newline characters \n as \\n
    Since it gives EOL error.
    Note2: Doctest is also sensitive to whitespace at the end of lines.
    Hence, I used .strip() to awoid unnecessary whitespace.
    >>> xyz = "2\\n \\n H 0. 0. 0.\\n H 0. 0. 0.9\\n  \\n"
    >>> logging.info get_mopac_input(xyz,method='pm7',dothermo=True)
    pm7 precise nosym threads=1 opt
    <BLANKLINE>
    <BLANKLINE>
    H   0.00000 1  0.00000 1  0.00000 1
    H   0.00000 1  0.00000 1  0.90000 1
    <BLANKLINE>
    pm7 precise nosym threads=1 oldgeo thermo
    """
    if type(x) is str:
        mol = ob.get_mol(x)
    else:
        mol = x
    multDictionary = {
                        1: '',
                        2: 'uhf doublet',
                        3: 'uhf triplet',
                        4: 'uhf quartet',
                        5: 'uhf quintet',
                        6: 'uhf sextet',
                        7: 'uhf septet',
                        8: 'uhf octet',
                        9: 'uhf nonet',
                    }
    keys = method + ' ' + keys + ' ' + multDictionary[mult]
    inp = ob.get_mop(mol, keys=keys.strip())
    if dothermo:
        inp += '\n' + keys.replace('opt', 'oldgeo thermo').strip()
    return inp


def get_gaussian_islinear(s):
    """
    Returns true if the molecule is linear for the given log.
    """
    if "Linear Molecule" in s or get_gaussian_natom(s) == 2:
        return True
    else:
        return False


def get_gaussian_nfreq(s):
    """
    Return the number of vibrational degrees of freedom for
    a given log.
    """
    natom = get_gaussian_natom(s)
    if natom == 1:
        nvdof = 0
    elif get_gaussian_islinear(s):
        nvdof = 3*natom - 5
    else:
        nvdof = 3*natom - 6
    return nvdof


def get_mopac_natom(lines):
    """
    Return the number of atoms from mopac output
    >>> s = io.read_file('test/input.out')
    >>> logging.info get_mopac_natom(s)
    5
    """
    import iotools as io
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'Empirical Formula'
    n = io.get_line_number(keyword, lines=lines)
    natom = int(lines[n].split()[-2])
    return natom


def get_mopac_xyz(lines):
    """
    Returns xyz string from mopac output lines.
    >>> s = io.read_file('test/input.out')
    >>> logging.info get_mopac_xyz(s)
    5
    <BLANKLINE>
    C 0.0000 -0.0000 0.0000
    H 1.0870 -0.0000 0.0000
    H -0.3623 1.0248 0.0000
    H -0.3624 -0.5124 0.8875
    H -0.3624 -0.5124 -0.8875
    <BLANKLINE>
    """
    import iotools as io
    if type(lines) == str:
        lines = lines.splitlines()
    natom = get_mopac_natom(lines)
    keyword = "ORIENTATION OF MOLECULE IN FORCE CALCULATION"
    xyzline = io.get_line_number(keyword, lines=lines) + 4
    xyz = '{0}\n'.format(natom)
    comment = '\n'
    xyz += comment
    for i in range(natom):
        xyz += ' '.join(lines[xyzline + i].split()[1:]) + '\n'
    return xyz


def get_mopac_freq(lines):
    """
    Returns a float list of vibrational frequencies in cm-1.
    >>> s = io.read_file('test/input.out')
    >>> logging.info get_mopac_freq(s)
    [ 1362.21  1362.44  1362.56  1451.04  1451.06  3207.4   3207.46  3207.59
      3310.99]
    """
    import numpy as np
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'FREQ.'
    natom = get_mopac_natom(lines)
    freqs = np.zeros(3 * natom)
    i = 0
    for line in lines:
        if keyword in line:
            freqs[i] = float(line.split()[1])
            i += 1
    return freqs[:i]


def get_mopac_zpe(lines):
    """
    Return zero point energy in kcal/mol from mopac output.
    >>> s = io.read_file('test/input.out')
    """
    from unittools import au2kcal
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'ZERO POINT ENERGY'
    n = io.get_line_number(keyword, lines=lines)
    return float(lines[n].split()[3]) / au2kcal

def get_mopac_energy(lines):
    """
    Return electronic energy in hartree from mopac output.
    >>> s = io.read_file('../test/input.out')
    """
    from unittools import ev2au
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'ELECTRONIC ENERGY       ='
    n = io.get_line_number(keyword, lines=lines, getlastone=True)
    return float(lines[n].split()[3]) * ev2au


def get_mopac_deltaH0(lines):
    """
    Return delta H in kcal/mol from mopac output.
    >>> s = io.read_file('test/input.out')
    >>> print(get_mopac_deltaH(s))
    -13.02534
    """
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'FINAL HEAT OF FORMATION'
    n = io.get_line_number(keyword, lines=lines)
    return float(lines[n].split()[5])


def get_mopac_deltaH298(lines):
    """
    Return delta H in kcal/mol from mopac output.
    >>> s = io.read_file('test/input.out')
    >>> print(get_mopac_deltaH(s))
    -13.02534
    """
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'CALCULATED THERMODYNAMIC PROPERTIES'
    n = io.get_line_number(keyword, lines=lines,getlastone=True)
    return float(lines[n+10].split()[1])


def get_nwchem_xyz(inp,filename=False):
    """
    Returns geometry in xyz format by parsing NWChem output file.
    Sample NWChem output:
 Output coordinates in angstroms (scale by  1.889725989 to convert to a.u.)

  No.       Tag          Charge          X              Y              Z
 ---- ---------------- ---------- -------------- -------------- --------------
    1 C                    6.0000     0.00000000     0.00000000     0.00000922
    2 H                    1.0000     0.00000000     0.00000000     1.09304166
    3 H                    1.0000     0.00000000    -0.94660523    -0.54652544
    4 H                    1.0000     0.00000000     0.94660523    -0.54652544

    """
    if filename:
        lines = io.read_file(inp,aslines=True)
    else:
        if type(inp) is str:
            lines = inp.splitlines()
        else:
            lines = inp
    keyword = 'No.       Tag          Charge          X              Y              Z'
    n = io.get_line_number(keyword, lines, getlastone=True)
    geolines = ''
    natom = 0
    for line in lines[n+2:]:
        items = line.split()
        if len(items) == 6:
            geolines += '{0}    {1}     {2}     {3}\n'.format(items[1],items[3], items[4], items[5])
            natom += 1
        else:
            break
    xyz = '{0}\nParsed by QTC from NWChem output file\n{1}\n\n'.format(natom, geolines)
    return xyz


def get_nwchem_energies(inp, filename=False):
    if filename:
        lines = io.read_file(inp,aslines=True)
    else:
        if type(inp) is str:
            lines = inp.splitlines()
        else:
            lines = inp
    nwdict = {
        'nre'        : 'Effective nuclear repulsion energy (a.u.)',
        'scf'        : 'Total SCF energy',
        'mp2'        : 'Total MP2 energy',
        'mbpt(2)'      : 'MBPT(2) total energy / hartree',
        'mp3'        : 'Total MP3 energy',
        'ccsd'       : 'CCSD total energy / hartree',
        'ccsd(t)'    : 'CCSD(T) total energy / hartree',
        'ccsd(2)_t'  : 'CCSD(2)_T total energy / hartree',
        'ccsd(2)'    : 'CCSD(2) total energy / hartree',
        'ccsdt'      : 'CCSDT total energy / hartree',
        'ccsdt(2)_q' : 'CCSDT(2)_Q total energy / hartree',
        'ccsdtq'     : 'CCSDTQ total energy / hartree'
    }
    energies={}
#    energies = {'unit':'hartree'}
    for key,value in nwdict.iteritems():
        i = io.get_line_number(value,lines=lines,getlastone=True)
        if i >= 0:
            try:
                energies[key] = float(lines[i].split()[-1])
            except:
                logging.info('Cannot parse {0}'.format(value))
    return energies


def get_nwchem_calculation(inp, filename=False):
    if filename:
        inp = io.read_file(inp,aslines=False)
    if 'Optimization converged' in inp:
        calc = 'geometry optimization'
    elif 'P.Frequency' in inp:
        calc = 'frequency analysis'
    else:
        calc = 'single point'
    return calc
    
    
def get_nwchem_method(inp, filename=False):
    if filename:
        inp = io.read_file(inp,aslines=False)
    nwdict = {
        0	:{'nre'        : 'Effective nuclear repulsion energy (a.u.)'},
        1	:{'scf'        : 'Total SCF energy'},
        2	:{'mp2'        : 'Total MP2 energy'},
        3   :{'mbpt(2)'    : 'MBPT(2) total energy / hartree       ='},
        4	:{'mp3'        : 'Total MP3 energy'},
        5	:{'ccsd'       : 'CCSD total energy / hartree'},
        6	:{'ccsd(t)'    : 'CCSD(T) total energy / hartree'},
        7	:{'ccsd(2)_t'  : 'CCSD(2)_T total energy / hartree'},
        8	:{'ccsd(2)'    : 'CCSD(2) total energy / hartree'},
        9	:{'ccsdt'      : 'CCSDT total energy / hartree'},
        10	:{'ccsdt(2)_q' : 'CCSDT(2)_Q total energy / hartree'},
        11	:{'ccsdtq'     : 'CCSDTQ total energy / hartree'}
    }
    method = 'unknown'
    for i in range(11,-1,-1):
        if nwdict[i].values()[0] in inp:
            method = nwdict[i].keys()[0]
            break
    return method.lower()

    
def get_nwchem_basis(inp, filename=False):
    """
------------------------------------------------------------------------------
       Tag                 Description            Shells   Functions and Types
 ---------------- ------------------------------  ------  ---------------------
 C                        aug-cc-pvdz                9       23   4s3p2d
 H                        aug-cc-pvdz                5        9   3s2p

    """
    if filename:
        lines = io.read_file(inp,aslines=True)
    else:
        if type(inp) is str:
            lines = inp.splitlines()
        else:
            lines = inp
    key = 'Tag                 Description            Shells   Functions and Types'
    i = io.get_line_number(key,lines,getlastone=True)
    basis = []
    nbasis = 0
    for line in lines[i+2:]:
        items = line.split()
        if len(items) == 5:
            basis.append(items[1])
            nbasis += int(items[-2])
        else:
            break
    if len(set(basis)) > 1:
        basis = set(basis)
        basis = '_'.join(basis)
    else:
        basis = basis[0]
    return {'basis': basis,'nbasis': nbasis}


def get_nwchem_nbasis(out):
    """
    Returns the number of basis functions by parsing NWChem output.
    Typical output:
    functions       =    24
    """
    keyword = "functions       ="
    lines = out.splitlines()
    for line in lines:
        if keyword in line:
            if len(line.split()) == 3:
                try:
                    nbasis = int(line.split()[2])
                    break
                except:
                    logging.error('Parser error for get_nwchem_nbasis  in line {}'.format(line))
    return nbasis


def get_nwchem_frequencies(inp, filename=False, minfreq=10):
    """
             (Projected Frequencies expressed in cm-1)

                    1           2           3           4           5           6

 P.Frequency        0.00        0.00        0.00        0.00        0.00        0.00

           1     0.00000     0.11409     0.07801     0.21786     0.00000     0.00000
           2    -0.00312     0.00000     0.00000     0.00000     0.00172     0.25797
           3    -0.01627     0.00000     0.00000     0.00000     0.25748    -0.00191
           4     0.00000    -0.45282    -0.30962     0.65355     0.00000     0.00000
           5     0.57079     0.00000     0.00000     0.00000     0.03802     0.26467
           6    -0.01627     0.00000     0.00000     0.00000     0.25748    -0.00191
           7     0.00000     0.79511    -0.30961     0.00000     0.00000     0.00000
           8    -0.29008     0.00000     0.00000     0.00000    -0.01644     0.25462
           9     0.48076     0.00000     0.00000     0.00000     0.28892     0.00389
          10     0.00000     0.00000     0.85326     0.00000     0.00000     0.00000
          11    -0.29008     0.00000     0.00000     0.00000    -0.01644     0.25462
          12    -0.51329     0.00000     0.00000     0.00000     0.22603    -0.00771

                    7           8           9          10          11          12

 P.Frequency      498.18     1406.65     1406.83     3103.34     3292.00     3292.33

           1    -0.12950     0.00000     0.00000     0.00000     0.00000     0.00000
           2     0.00000     0.00000     0.08818     0.00000    -0.09484     0.00000
           3     0.00000    -0.08818     0.00000    -0.00009     0.00000    -0.09484
           4     0.51400     0.00000     0.00000     0.00000     0.00000     0.00000
           5     0.00000     0.00000    -0.77117     0.00000    -0.01518     0.00000
           6     0.00000    -0.07120     0.00000    -0.57437     0.00000     0.76857
           7     0.51398     0.00000     0.00000     0.00000     0.00000     0.00000
           8     0.00000    -0.36472    -0.13940     0.49839     0.57222     0.33868
           9     0.00000     0.56060     0.36475     0.28770     0.33914     0.18033
          10     0.51398     0.00000     0.00000     0.00000     0.00000     0.00000
          11     0.00000     0.36472    -0.13940    -0.49839     0.57222    -0.33868
          12     0.00000     0.56060    -0.36475     0.28770    -0.33914     0.18033 
    """
    if filename:
        lines = io.read_file(inp,aslines=True)
    else:
        if type(inp) is str:
            lines = inp.splitlines()
        else:
            lines = inp
    key = 'P.Frequency'
    nums = io.get_line_numbers(key,lines)
    freqs = []
    if nums is not -1:
        for num in nums:
            line = lines[num]
            for item in line.split()[1:]:
                freq = item.strip()
                if float(freq) > minfreq:
                    freqs.append(freq)
    return freqs


def get_torsscan_info(s):
    """
    Returns electronic prog/method/basis as a string and energy as a float in hartree units.
    Sample torsscan output
Optimized at : g09/b3lyp/sto-3g
Prog  : gaussian
Method: b3lyp
Basis:  sto-3g
Energy: -114.179051157 A.U.
Rotational Constants:120.49000, 23.49807, 22.51838 GHz
    """
    lines  =s.splitlines()
    optlevel = ''
    prog = ''
    method = ''
    basis = ''
    energy = ''
    for line in lines:
        line = line.strip()
        if 'Optimized at' in line:
            optlevel = line.split()[-1]
        elif line.startswith('Prog'):
             s = line.split()[-1]
        elif line.startswith('Method'):
            method = line.split()[-1]
        elif line.startswith('Basis'):
            basis = line.split()[-1]
        elif line.startswith('Energy'):
            energy = float(line.replace('A.U.','').split()[-1])
    return optlevel,  '{}/{}/{}'.format('torsscan',method,basis), energy


def get_mess_xyz(out):
    """
    Return xyz as a string by parsing mess input file.
    Input follows Geometty line
    Sample input file:
    AtomDistanceMin[angstrom] 0.6
    Temperature(step[K],size)        100.   30
    RelativeTemperatureIncrement            0.001
    Species H2O
    RRHO
    Geometry[angstrom] 3 !torsscan/b3lyp/sto-3g/gaussian
     O  0.000000  0.000000  0.000000
     H  0.000000  0.000000  1.027000
     H  1.018909  0.000000  -0.128659
    
    Core RigidRotor
    SymmetryFactor 1
    End
    Frequencies[1/cm] 3 !torsscan/b3lyp/sto-3g/gaussian
    2016.9034 3614.5724 3839.1980
    ZeroEnergy[kcal/mol] 0.0215757824512 ! torsscan/b3lyp/sto-3g/gaussian
    ElectronicLevels[1/cm]  1
    0 1
    End
    """
    lines = out.splitlines()
    n = io.get_line_number('Geometry[angstrom]', lines)
    natom = int(lines[n].split()[-1])
    xyz = str(natom) + '\n\n'
    i = 0
    for line in lines[n+1:]:
        if line.strip() and len(line.split()) == 4:
            xyz += line + '\n'
            i += 1
            if i == natom:
                break
        else:
            logging.warning('get_mess_xyz: Invalid xyz format in mess input')
    return xyz



def get_mess_frequencies(out):
    """
    Return frequencies as a list of floats by parsing mess input file.
    Input follows Frequencies line
    Sample input file:
    AtomDistanceMin[angstrom] 0.6
    Temperature(step[K],size)        100.   30
    RelativeTemperatureIncrement            0.001
    Species H2O
    RRHO
    Geometry[angstrom] 3 !torsscan/b3lyp/sto-3g/gaussian
     O  0.000000  0.000000  0.000000
     H  0.000000  0.000000  1.027000
     H  1.018909  0.000000  -0.128659
    
    Core RigidRotor
    SymmetryFactor 1
    End
    Frequencies[1/cm] 3 !torsscan/b3lyp/sto-3g/gaussian
    2016.9034 3614.5724 3839.1980
    ZeroEnergy[kcal/mol] 0.0215757824512 ! torsscan/b3lyp/sto-3g/gaussian
    ElectronicLevels[1/cm]  1
    0 1
    End
    >>> out = io.read_file('sample_logfiles/CH4O_pf.inp')
    >>> print qc.get_mess_frequencies(out)
    [3690.23, 3425.92, 3353.85, 3218.7, 1684.07, 1662.1, 1604.55, 1582.52, 1178.31, 1152.3, 1081.44, 398.92]
    """
    out = out.lower()
    lines = out.splitlines()
    n = io.get_line_number('frequencies', lines)
    freqs = []
    for line in lines[n+1:]:
        if line.islower(): #Returns false if line has no letters
            break
        elif line.strip():
            items = line.split()
            for item in items:
                try:
                    freqs.append(float(item))
                except:
                    logging.error('Non-numeric string in frequency lines of mess input: {0}'.format(item))
        else:
            pass
    return sorted(freqs)


def print_list(s):
    return s

def remove_dups(seq):
    """
    https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-whilst-preserving-order
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

