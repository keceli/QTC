#!/usr/bin/env python
"""
Quantum chemistry tools.
"""
import obtools as ob
import iotools as io
import numpy as np
import patools as pa
from Cython.Build.Dependencies import package
try:
    import cclib
except:
    pass

__updated__ = "2017-08-31"
_hartree2kcalmol = 627.509 #kcal mol-1

def get_input(x, template, parameters):
    """
    Returns input file text for a qc calculation based on a given template.
    """
    mol = ob.get_mol(x)
    mult = ob.get_multiplicity(mol)
    nopen = mult - 1
    natom = ob.get_natom(mol)
    charge = ob.get_charge(mol)
    formula = ob.get_formula(mol)
    geo = ob.get_geo(mol)
    xyz = ob.get_xyz(mol)
    zmat = ob.get_zmat(mol)
    uniquename = ob.get_inchi_key(mol, mult)
    smilesname = ob.get_smiles_filename(mol)
    smiles = ob.get_smiles(mol)
    package = parameters['qcpackage'] 
    method  = parameters[ 'qcmethod'] 
    basis   = parameters[  'qcbasis']
    if 'results' in parameters.keys():
        results = parameters['results']
        if 'HoF at 0 K' in results.keys():
            heat = results['HoF at 0 K']
        else:
            heat = 0.
    else:
        heat = 0.    
    task    = parameters['qctask']
    nproc   = parameters['qcnproc']
    inp = template
    if task.startswith('tors'):
        if parameters['optdir']:
            xyzfile =  io.join_path(*[parameters['smilesdir'],parameters['optdir'], str(x).strip() + '.xyz'])
            if io.check_file('xyzfile'):
                inp = inp.replace("QTC(OPTDIR)",xyzfile)
            else:
                inp = inp.replace("QTC(OPTDIR)", 'false') 
        else:
            inp = inp.replace("QTC(OPTDIR)", 'false') 
        inp = inp.replace(      "QTC(NMC)", str(3**parameters['nrotor']+2) )
        inp = inp.replace(  "QTC(HFBASIS)", parameters[  'hfbasis'])
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
            elif task == 'energy':
                task = ''
            elif task == 'anharm':
                task = 'freq=(anharm,vibrot)'
        elif package == 'molpro':
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
    inp = inp.replace("QTC(SMILES)", smiles)
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
    if "QTC(" in inp:
        print(66*'#')
        print("Error in template file: \n" + inp)
        print(66*'#')
    return inp

def update_qckeyword(keyword):
    """
    Update parameters['qckeyword']  to minimize the differences
    based on user input.
    """
    keyword = keyword.lower()
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
    keyword = keyword.replace('631g','6-31g')
    keyword = keyword.replace('6-31gs','6-31g*')
    keyword = keyword.replace('torscan/','torsscan/')
    return keyword

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
    runqc
    parseqc
    writefiles
    anharmonic
    optlevel
    freqlevel
    enlevel
    """
    keyword = parameters['qckeyword']
    xyzdir = parameters['xyzpath']
    optdir   = parameters['optdir']
    package = 'nwchem'
    calcs = keyword.split(',')
    currentcalc = calcs[calcindex]
    tokens = currentcalc.split('/')
    task = tokens[0]
    if parameters['natom'] == 1:
        if task in ['opt', 'freq', 'torsscan', 'torsopt', 'anharm']:
            print('Switching to Task: energy, since Task: {0} can not be run for 1 atom'.format(task))
            task = 'energy'
            parameters['task'] = task
    if task.startswith('ext') or task.startswith('cbs') or task.startswith('comp'):
        task = 'composite'
        if len(tokens) > 2:
            method = tokens[1]
            parameters['formula'] = tokens[2:]
        elif len(tokens) == 2:
            method = 'generic'
            parameters['formula'] = tokens[1]
        else:
            print('ERROR! Invalid qckeyword: {0}'.format(tokens))
            return
        qcdirectory = io.fix_path(io.join_path(*[xyzdir,optdir,task,method]))
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
        else:
            print('ERROR! Invalid qckeyword: {0}'.format(tokens))
            return        
        if task.startswith('opt')  or task.startswith('geo') or task.startswith('min'):
            task = 'opt'
            qcdirectory = io.fix_path(io.join_path(*[xyzdir,task,method,basis,package]))
            parameters['optdir'] = qcdirectory
            parameters['optlevel'] = '{}/{}/{}'.format(package,method,basis)
        elif task.startswith('tors'):
            qcdirectory = io.fix_path(io.join_path(*[xyzdir,task,method,basis,package]))
            parameters['optdir'] = qcdirectory
            parameters['optlevel'] = '{}/{}/{}'.format(package,method,basis)
            parameters['freqdir'] = qcdirectory
            parameters['freqlevel'] = '{}/{}/{}'.format(package,method,basis)
        elif task.startswith('freq') or task.startswith('harm') or task.startswith('hrm'):
            task = 'freq'
            qcdirectory = io.fix_path(io.join_path(*[xyzdir,optdir,task,method,basis,package]))
            parameters['freqdir'] = qcdirectory
            parameters['freqlevel'] = '{}/{}/{}'.format(package,method,basis)
        elif task.startswith('anh') or task.startswith('afre'):
            task = 'anharm'
            qcdirectory = io.fix_path(io.join_path(*[xyzdir,optdir,task,method,basis,package]))
            parameters['anharmdir'] = qcdirectory
            parameters['anharmonic'] = True
        elif task.startswith('sp') or task.startswith('single') or task.startswith('ene'):
            task = 'energy'
            qcdirectory = io.fix_path(io.join_path(*[xyzdir,optdir,task,method,basis,package]))
            parameters['energylevel'] = '{}/{}/{}'.format(package,method,basis)
        else:
            print('ERROR! Invalid qckeyword task: {0}'.format(task))
            return      
    parameters['qcdirectory'] = qcdirectory
    parameters['qcpackage'] = package
    parameters['qcmethod'] = method
    parameters['qcbasis'] = basis
    parameters['qctask'] = task
    parameters['runqc'] = True
    parameters['parseqc'] = True
    parameters['writefiles'] = True
    return parameters 


def get_qc_label(natom, qckeyword, calcindex):
    """
    Returns a string that can be used as a label for
    a given quantum chemistry calculation with natom, qckeyword
    and calcindex. The string considers dependencies of the
    corresponding calculation.
    """
    calcs = qckeyword.split(',')
    if natom == 1:
        qckeyword = qckeyword.replace('opt', 'energy')
        qckeyword = qckeyword.replace('freq', 'energy')
        qckeyword = qckeyword.replace('anharm', 'energy')
        qckeyword = qckeyword.replace('torsscan', 'energy')
        qckeyword = qckeyword.replace('torsopt', 'energy')
        calcs = qckeyword.split(',')
        calc  = calcs[calcindex]
        if calc.startswith('compos'):
            label = ','.join(calcs[0:calcindex+1])
        else:
            label = calcs[calcindex]
    else:
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


def parse_results(filename, parameters):
    """
    Returns a dictionary containing results for the given quantum chemistry output.
    package: <str>
    method : <str>
    file: <str>
    hostname: <str>
    xyz: <str>
    xmat: <str>  anharmonicity X matrix from VPT2
    hindered rotor potential: <str> messpf input section for rotors 
    number of basis functions: <int>
    nuclear repulsion energy: <float>
    energy: <float> in atomic units (Hartree)
    zpve: <float> in cm^-1
    anharmonic zpve: <float> in cm^-1
    harmonic frequencies: <[float]> in cm^-1
    anharmonic frequencies: <[float]> in cm^-1
    projected frequencies: <[float]> in cm^-1
    """
    import iotools
    r = {}
    r['file'] = ''
    r['package'] = ''
    r['method'] = ' '
    r['hostname'] = ' '
    r['xyz'] = ' '
    r['xmat'] = ' '
    r['hindered rotor potential'] = ' ' 
    r['number of basis functions'] = 0
    r['nuclear repulsion energy:' ]= 0.
    r['energy' ]= 0.
    r['zpve'] = 0.
    r['anharmonic zpve'] = 0.
    r['harmonic frequencies' ]= []
    r['anharmonic frequencies'] = []
    r['projected frequencies' ]= []
    if io.check_file(filename,timeout=1):
        r['file'] = io.get_path(filename)
        out = io.read_file(filename,asline=False)
        package = get_output_package(out)
        if package:
            r['package'] = package
        else:
            print('"{0}" can not be parsed.'.format(r['file']))
            return r
    else:
        print('"{0}" can not be found at {1}'.format(filename,io.pwd()))
        return r
#     r['xyz'] = get_xyz(out, package)
#     r['xmat'] = get_xyz(out, package)
#     r['hindered rotor potential'] = get_xyz(out, package) 
#     r['number of basis functions'] = get_xyz(out, package)
#     r['nuclear repulsion energy:' ]= get_xyz(out, package)
#     r['energy' ]= get_xyz(out, package)
#     r['zpve'] = get_xyz(out, package)
#     r['anharmonic zpve'] = get_xyz(out, package)
#     r['harmonic frequencies' ]= get_xyz(out, package)
#     r['anharmonic frequencies'] = get_xyz(out, package)
#     r['projected frequencies' ]= get_xyz(out, package)
    return r


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
    return xyz
        
        
def parse_output(s, smilesname, write=False, store=False, optlevel='sp'):
    package = get_output_package(s)
    if type(s) is list:
        lines = s
        s = ''.join(lines)
    elif type(s) is str:
        lines = s.splitlines()
    else:
        print("First parameter in parse_output should be a string or list of strings")
    d = {}
    [method,calculation,xyz,basis] = ['na']*4
    nbasis = 0
    energy = 0
    energies = {}
    freqs = []
    pfreqs = []
    afreqs = []
    xmat= []
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
        basisinfo = get_nwchem_basis(lines)
        basis = basisinfo['basis']
        nbasis = basisinfo['number of basis functions']
        energies = get_nwchem_energies(lines)
        energy = energies[method]
        freqs = get_nwchem_frequencies(lines)
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
        freqs       = pa.molpro_freqs(s)
        parsed = True
    elif package == 'gaussian':
        method, energy = pa.gaussian_energy(s)
        zpve           = pa.gaussian_zpve(s)
        azpve       = pa.gaussian_anzpve(s)
        calculation    = pa.gaussian_calc(s)
        basis          = pa.gaussian_basisset(s)
        zmat           = pa.gaussian_zmat(s)
        xyz            = pa.gaussian_xyz(s)
        geo            = pa.gaussian_geo(s)
        hessian        = pa.gaussian_hessian(s)
        freqs       = pa.gaussian_freqs(s)
        afreqs  = get_gaussian_fundamentals(s)[:,1]
        if sum(afreqs) > 0:
            xmat           = get_gaussian_xmatrix(s, get_gaussian_nfreq(s))
            if type(xmat) == str:
                xmat = []
        parsed = True
    elif package == 'torsscan':
        optlevel, method, energy = get_torsscan_info(s)
        if io.check_file('geom.log'):
                out = io.read_file('geom.log', aslines=False)
                xyz = pa.gaussian_xyz(out)
                method, energy = pa.gaussian_energy(out)
                freqs = pa.gaussian_freqs(out)
                parsed = True
        if io.check_file('hrproj_freq.dat'):
            freqlines = io.read_file('hrproj_freq.dat',aslines=True)
            for line in freqlines:
                if float(line) > 0.:
                    pfreqs.append(float(line))
        if io.check_file('RTproj_freq.dat'):
            freqlines = io.read_file('RTproj_freq.dat',aslines=True)
            for line in freqlines:
                if float(line) > 0.:
                    freqs.append(line)
        fname = 'me_files/reac1_zpe.me'
        if io.check_file(fname):
            zpve = float(io.read_file(fname, aslines=False))                
        fname = 'me_files/reac1_hr.me'
        if io.check_file(fname):
            messhindered = io.read_file(fname, aslines=False)
        fname = 'RPHt_input_data.dat'
        if io.check_file(fname):
            RPHtinput = io.read_file(fname, aslines=False)

    if parsed:
        if write:
            fname = smilesname + '.xyz'
            io.write_file(xyz, fname)
            fname = smilesname + '.ene'
            io.write_file(str(energy), fname)
            if zpve:
                fname = smilesname + '.zpve'
                io.write_file(str(zpve), fname)
            if azpve:
                fname = smilesname + '.anzpve'
                io.write_file(str(azpve), fname)
            if len(freqs) > 0:
                fname = smilesname + '.hrm'
                io.write_file('\n'.join(str(x) for x in freqs), fname )
            if sum(afreqs) > 0:
                fname = smilesname + '.anhrm'
                io.write_file('\n'.join(str(x) for x in afreqs), fname)
        d = {'number of basis functions':nbasis,
               'energy':energy,
               'xyz':xyz,
               'harmonic frequencies' : freqs,
               'anharmonic frequencies' : afreqs,
               'projected frequencies': pfreqs,
               'zpve': zpve,
               'anharmonic zpve': azpve,
               'xmat': xmat,
               'hindered potential': messhindered,
               'RPHt input': RPHtinput,
               'Hessian'   : hessian,
               'HoF at 0 K': hof0,
               'HoF at 298 K': hof298}
#         if calculation == 'geometry optimization':
#             for key,value in energies.iteritems():
#                 if key is not method:
#                     d[optlevel][package][calculation][method][basis]['geometry'].update({
#                         'single point':{key:{basis:{'number of basis functions':nbasis,'energy':value}}}})
#                     if write:
#                         fname = '{0}_{1}.ene'.format(method,smilesname)
#                         io.write_file(str(energy), fname)
        if store:
            if optlevel == 'sp':
                if energy:
                    io.db_store_sp_prop(str(energy), smilesname,  'ene', None, package, method, basis)
                if zpve:
                    io.db_store_sp_prop(str(  zpve), smilesname, 'zpve', None, package, method, basis)
                if len(hrmfreqs) > 0:
                    io.db_store_sp_prop(', '.join(freq for freq in hrmfreqs[::-1]) , smilesname,  'hrm', None, package, method, basis)
            else:
                opt1, opt2, opt3 = optlevel.split('/')
                if energy:
                    io.db_store_sp_prop(str(energy), smilesname,  'ene', None, package, method, basis, opt1, opt2, opt3)
                if zpve:
                    io.db_store_sp_prop(str(  zpve), smilesname, 'zpve', None, package, method, basis, opt1, opt2, opt3)
                if zpve:
                    io.db_store_sp_prop(str(anzpve), smilesname,'anzpve', None, package, method, basis, opt1, opt2, opt3)
                if len(hrmfreqs) > 0:
                    io.db_store_sp_prop(', '.join(freq for freq in hrmfreqs[::-1]) , smilesname,  'hrm', None, package, method, basis, opt1, opt2, opt3)
                if len(xmat) > 0:
                    io.db_store_sp_prop('\n'.join([','.join(['{:4}'.format(x) for x in xma]) for xma in xmat]) , smilesname, 'xmat', None, package, method, basis, opt1, opt2, opt3)
                io.db_store_sp_prop(', '.join(freq for freq in hrmfreqs[::-1]) , smilesname,  'hrm', None, package, method, basis, opt1, opt2, opt3)
            if xyz:
                io.db_store_opt_prop(xyz, smilesname,  'xyz', None, package, method, basis)
            if geo:
                io.db_store_opt_prop(geo, smilesname,  'geo', None, package, method, basis)
            if zmat:
                io.db_store_opt_prop(zmat, smilesname, 'zmat', None, package, method, basis)
    return d


def get_listofstrings(array):
    """
    Return a list of strings from a given array
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


def get_symbol(atomno):
    """
    Returns the element symbol for a given atomic number.
    Returns 'X' for atomno=0
    >>>print get_symbol(1)
    >>>H
    """
    syms = ['X',
            'H', 'He'
            'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
    return syms[atomno]


def run(mol, parameters, mult=None):
    """
    Runs qc, returns a string specifying the status of the calculation.
    """
    package = parameters['qcpackage']
    overwrite = parameters['overwrite']
    template = parameters['qctemplate']
    task = parameters['qctask']
    msg = ''
    if mult is None:
        mult = ob.get_multiplicity(mol)
    else:
        ob.set_mult(mol, mult)
    tmp = io.read_file(template)
    inptext = get_input(mol, tmp, parameters)
    if task.startswith('tors'):
        package = task
    prefix = ob.get_smiles_filename(mol) + '_' + package
    inpfile = prefix + '.inp'
    outfile = prefix + '.out'
    if io.check_file(outfile, timeout=1):
        if overwrite:
            msg = 'Overwriting previous calculation "{0}"\n'.format(io.get_path(outfile))
            run = True
        else:
            out = io.read_file(outfile)
            if check_output(out):
                msg = 'Skipping calculation, found "{0}"\n'.format(io.get_path(outfile))
                run = False
            else: 
                msg = 'Failed output found "{0}", renaming and running a new calculation\n'.format(io.get_path(outfile))
                io.mv(outfile, 'failed_'+outfile)
                run = True
    else:
        run = True
    if run:
        io.write_file(inptext, inpfile)
        if io.check_file(inpfile, timeout=1):
            if package in  ['nwchem', 'torsscan','torsopt']:
                command = parameters['qcexe'] + ' ' + inpfile
                msg += io.execute(command,stdoutfile=outfile,merge=True)
            elif package in  ['molpro']:
                command = parameters['qcexe'] + ' ' + inpfile
                msg += io.execute(command,stdoutfile=outfile,merge=True)
                logfile = prefix + '.log'
                if io.check_file(logfile):
                    io.append_file(io.read_file(logfile),filename=outfile)
                    io.rm(logfile)
            else:
                command = parameters['qcexe'] + ' ' + inpfile + ' ' + outfile
                msg += io.execute(command)
            if io.check_file(outfile, timeout=1):
                msg += ' Output file: "{0}"\n'.format(io.get_path(outfile))
        else:
            msg += 'Failed, cannot find input file "{0}"\n'.format(io.get_path(inpfile))
    return msg


def run_extrapolation(mol,parameters):
    if parameters['qckeyword']:
        msg = run_extrapolation_keyword(mol,parameters)
    elif parameters['qctemplate']:
        msg = run_extrapolation_template(mol, parameters)
    else:
        msg = 'Can not run extrapolation, you need to specify qckeyword with -k or qctemplate with -t. \n'
    return msg


def run_extrapolation_keyword(mol, parameters):
    """
    Parses qckeyword for composite method. 
    'opt/mp2/cc-pvdz/gaussian,freq/mp2/cc-pvtz/molpro,sp/mp2/cc-pvqz,composite/cbs-dtq/energy=0.1 * E[0] + 0.4 * E[1] + 0.5 * E[2]'
    """
    keyword = parameters['qckeyword']
    formula = parameters['formula'][0]
    method = parameters['qcmethod']
    msg = ''
    smilesname = ob.get_smiles_filename(mol)
    calcs = keyword.split(',')
    print('Composite energy formula: {0}\n'.format(formula))
    ncalc = len(calcs) 
    calcindex = parameters['calcindex']
    e = [0.] * ncalc
    energy = None
    enefile = smilesname + '.ene'  
    inpfile = smilesname + '_' + method  + '.inp'  
    for i in range(0,calcindex):
        parameters = parse_qckeyword(parameters, i)
        smilesdir = parameters['smilesdir']
        enedir  = io.join_path(*[smilesdir,parameters['qcdirectory']])
        enepath = io.join_path(*[enedir, enefile])
        if io.check_file(enepath):
            e[i] = float(io.read_file(enepath))
            msg += ('E({0}) from "{1}"  = {2}\n'.format(i,enepath,e[i]))
        else:
            msg += ('E({0}) can not be found at "{1}" \n'.format(i,enepath))
    print(msg)
    msg = ''
    exec(formula)
    if energy:
        io.write_file(str(energy),enefile )
        io.write_file(msg,inpfile )
        print('Composite energy: {}\n'.format(energy))
        print('Energy file: "{}"\n'.format(io.get_path(enefile)))
    parse_qckeyword(parameters, calcindex)
    return msg

def run_extrapolation_template(s, parameters):
    lines = io.read_file(parameters['qctemplate'],aslines=True)
    smilesname = ob.get_smiles_filename(s)
    filename = smilesname + '_cbs.ene'
    qcdir = parameters['qcdirectory']
    directories = []
    msg = ''
    for line in lines:
        if 'directories=' in line:
            exec(line)
            ndir = len(directories)
            energies=[0.]*ndir
            for i, edir in enumerate(directories):
                efile = io.join_path(*[edir,smilesname+'.ene'])
                if io.check_file(efile,verbose=True):        
                    energies[i] = float(io.read_file(efile, aslines=False))
                    print('Reading energy from {0} = {1}'.format(edir,energies[i]))
    for line in lines:
        if 'energy=' in line:
            energy = 0
            exec(line)
            print('Extrapolation based on formula: {0}'.format(line))        
            print('Extrapolated energy = {0}'.format(energy))
        if 'filename=' in line:
            exec(line)
    if len(directories) < 1:
        print('You have to specifies directories as a list in the template file')         
    if energy:
        msg += 'Extrapolation successful'
        if parameters['writefiles']:
            if qcdir:
                filename = io.join_path(*[qcdir,filename])
            io.write_file(str(energy), filename)
            msg += 'Extrapolation enegy file {0}'.format(filename)
    else:
        msg += 'Extrapolation failed'      
    return msg

                    
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
    elif "Task: Submitting EStokTP job..." in s:
        if 'completed' in s:
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
    Returns the name of qc package if the calculations is succesful.
    Returns None if failed or unknown package.
    """
    if filename:
        out = io.read_file(out,aslines=False)
    if "Normal termination of Gaussian" in out:
        p = 'gaussian'
    elif "== MOPAC DONE ==" in out:
        p = 'mopac'
    elif "Straatsma" in out:
        p = 'nwchem'
    elif "Variable memory released" in out:
        p = 'molpro'
    elif "Task: Submitting EStokTP job.." in out:
        p = 'torsscan'
    else:
        p = None
    return p

def get_package(templatename):
    """
    Return the quantum chemistry package name based on the template file name
    """
    suffix = templatename.split('.')[-1]
    if 'cbs' in templatename or suffix is 'py':
        p = 'extrapolation'
    elif 'nwchem' in templatename or 'nw' in suffix:
        p = 'nwchem'
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
        print("Error in template file:\n" + inp)
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
    while line.strip():
        cols = line.split()
        icol = int(cols[0])-1
        for irow in range(icol,nfreq):
            iline += 1
            line = lines[iline]
            cols = line.split()
            ncol = len(cols) - 1
            print cols[1:]
            xmat[irow,icol:icol+ncol] = [float(num.replace('D','E')) for num in cols[1:]]
        iline += 1
        line = lines[iline]
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
    >>> print get_mopac_input(xyz,method='pm7',dothermo=True)
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
    >>> print get_mopac_natom(s)
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
    >>> print get_mopac_xyz(s)
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
    >>> print get_mopac_freq(s)
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
    >>> print get_mopac_zpe(s)
    28.481
    """
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'ZERO POINT ENERGY'
    n = io.get_line_number(keyword, lines=lines)
    return float(lines[n].split()[3])


def get_mopac_deltaH(lines):
    """
    Return delta H in kcal/mol from mopac output.
    >>> s = io.read_file('test/input.out')
    >>> print get_mopac_deltaH(s)
    -13.02534
    """
    if type(lines) == str:
        lines = lines.splitlines()
    keyword = 'FINAL HEAT OF FORMATION'
    n = io.get_line_number(keyword, lines=lines)
    return float(lines[n].split()[5])


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
                print('Cannot parse {0}'.format(value))
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
        3	:{'mp3'        : 'Total MP3 energy'},
        4	:{'ccsd'       : 'CCSD total energy / hartree'},
        5	:{'ccsd(t)'    : 'CCSD(T) total energy / hartree'},
        6	:{'ccsd(2)_t'  : 'CCSD(2)_T total energy / hartree'},
        7	:{'ccsd(2)'    : 'CCSD(2) total energy / hartree'},
        8	:{'ccsdt'      : 'CCSDT total energy / hartree'},
        9	:{'ccsdt(2)_q' : 'CCSDT(2)_Q total energy / hartree'},
        10	:{'ccsdtq'     : 'CCSDTQ total energy / hartree'}
    }
    method = 'unknown'
    for i in range(10,-1,-1):
        if nwdict[i].values()[0] in inp:
            method = nwdict[i].keys()[0]
            break
    return method

    
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
    return {'basis': basis,'number of basis functions': nbasis}


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
                freq = float(item)
                if freq > minfreq:
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
            prog = line.split()[-1]
        elif line.startswith('Method'):
            method = line.split()[-1]
        elif line.startswith('Basis'):
            basis = line.split()[-1]
        elif line.startswith('Energy'):
            energy = float(line.replace('A.U.','').split()[-1])
    return optlevel,  '{}/{}/{}'.format('torsscan',method,basis), energy


def print_list(s):
    return s

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

