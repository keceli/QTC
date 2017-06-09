#!/usr/bin/env python
"""
Quantum chemistry tools.

"""
import subprocess
import obtools as ob
import argparse
import iotools as io
import numpy as np
try:
    import cclib
except:
    pass

__updated__ = "2017-06-09"
_hartree2kcalmol = 627.509 #kcal mol-1


def get_listofstrings(array):
    """
    Return a list of strings from a given array
    """
    n = len(array)
    s = ['']*n
    for i in range(n):
        s [i] = '{0}\n'.format(array[i]) 
    return s


def get_string(array):
    """
    Return a list of strings from a given array
    """
    n = len(array)
    s = ''
    for i in range(n):
        s += '{0}\n'.format(array[i]) 
    return s


def parse_qclog(qclog,qccode='gaussian',anharmonic=False):

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
    lines = s.splitlines()
    if check_logfile(s):
        try:
            if qccode == 'gaussian':
                mol = ob.get_gaussian_mol(qclog)
                zpe = mol.energy
                if cclib:
                    ccdata = parse_cclib(qclog)
                    xyz = ccdata.writexyz()
                    freqs = ccdata.vibfreqs
                    freqs = get_listofstrings(freqs)
                    nfreq = len(freqs)
                    deltaH = ccdata.enthalpy
                    if anharmonic:
                        xmat = ccdata.vibanharms
                        afreqs = get_gaussian_fundamentals(s, nfreq)[:,1]
                        afreqs = get_listofstrings(afreqs)
             
            elif qccode == 'mopac':
                xyz = get_mopac_xyz(lines)
                freqs = get_mopac_freq(lines)
                freqs = get_listofstrings(freqs)
                zpe = get_mopac_zpe(lines)
                deltaH = get_mopac_deltaH(lines)
        except:
            msg = 'Parsing failed for "{0}"\n'.format(io.get_path(qclog))
            return msg,xyz,freqs,zpe,deltaH,afreqs,xmat
    else:
        msg = 'Failed job: "{0}"\n'.format(io.get_path(qclog))
                
    return msg,xyz,freqs,zpe,deltaH,afreqs,xmat


def parse_qclog_as_dict(qclog,qccode='gaussian',anharmonic=False):

    xyz = None
    freqs = None
    zpe = None
    deltaH = None
    xmat = None
    afreqs = None
    basis = None
    method = None
    msg =''
    if io.check_file(qclog, 1):
        s = io.read_file(qclog, aslines=False)
    else:
        msg = 'File not found: "{0}"\n'.format(io.get_path(qclog)) 
        return msg,xyz,freqs,zpe,deltaH,afreqs,xmat   
    lines = s.splitlines()
    if check_logfile(s):
        try:
            if qccode == 'gaussian':
                mol = ob.get_gaussian_mol(qclog)
                zpe = mol.energy
                if cclib:
                    ccdata = parse_cclib(qclog)
                    xyz = ccdata.writexyz()
                    freqs = ccdata.vibfreqs
                    freqs = get_listofstrings(freqs)
                    nfreq = len(freqs)
                    deltaH = ccdata.enthalpy
                    if anharmonic:
                        xmat = ccdata.vibanharms
                        afreqs = get_gaussian_fundamentals(s, nfreq)[:,1]
                        afreqs = get_listofstrings(afreqs)
             
            elif qccode == 'mopac':
                xyz = get_mopac_xyz(lines)
                freqs = get_mopac_freq(lines)
                freqs = get_listofstrings(freqs)
                zpe = get_mopac_zpe(lines)
                deltaH = get_mopac_deltaH(lines)
        except:
            msg = 'Parsing failed for "{0}"\n'.format(io.get_path(qclog))
            return msg,xyz,freqs,zpe,deltaH,afreqs,xmat
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


def execute(inp, exe,outputfile=None):
    """
    Executes a calculation for a given input file inp and executable exe.
    """
    from subprocess import Popen, PIPE
    process = Popen([exe, inp], stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    if outputfile:
        io.write_file(out,outputfile)    
    if err is None or err == '':
        msg = 'Run {0} {1}: Success.\n'.format(exe, inp)
    else:
        errstr = """ERROR in "{0}"\n
        STDOUT:\n{1}\n
        STDERR:\n{2}""".format(inp, out, err)
        errfile = inp + '.err'
        io.write_file(errstr, errfile)
        msg = 'Run {0} {1}: Failed, see "{2}"\n'.format(exe, inp, io.get_path(errfile))
    return msg



def execute_gaussian(inp, exe='g09'):
    """
    Runs gaussian calculation.
    """
    from subprocess import Popen, PIPE
    process = Popen([exe, inp], stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    
    if err is None or err == '':
        msg = 'Run {0} {1}: Success.'.format(exe, inp)
    else:
        errstr = """ERROR in "{0}"\n
        STDOUT:\n{1}\n
        STDERR:\n{2}""".format(inp, out, err)
        errfile = inp + '.err'
        io.write_file(errstr, errfile)
        msg = 'Run {0} {1}: Failed, see "{2}"'.format(exe, inp, io.get_path(errfile))
    return msg


def run_gaussian(s, exe='g09', template='gaussian_template.txt', mult=None, overwrite=False):
    """
    Runs gaussian calculation
    """
    mol = ob.get_mol(s, make3D=True)
    if mult is None:
        mult = ob.get_multiplicity(mol)
    tmp = io.read_file(template)    
    inptext = get_gaussian_input(mol, tmp, mult)
    prefix = ob.get_smiles_filename(s)
    inpfile = prefix + '.g09'  
    outfile = prefix + '.log'
    if io.check_file(outfile, timeout=1):
        if overwrite:
            msg = 'Overwriting previous calculation "{0}"\n'.format(io.get_path(outfile))
            run = True
        else:
            msg = 'Skipping calculation, found "{0}"\n'.format(io.get_path(outfile))
            run = False
    else:
        run = True
    if run:
        if not io.check_file(inpfile, timeout=1):
            io.write_file(inptext, inpfile)
        if io.check_file(inpfile, timeout=1):
            msg = execute(inpfile, exe)
            if io.check_file(outfile, timeout=1):
                msg += ' Output file: "{0}"\n'.format(io.get_path(outfile))
        else:
            msg = 'Failed, cannot find input file "{0}"\n'.format(io.get_path(inpfile))
    return msg


def run_nwchem(s, exe='nwchem', template='nwchem_template.txt', mult=None, overwrite=False):
    """
    Runs nwchem, returns a string specifying the status of the calculation.
    nwchem inp.nw > nw.log
    NWChem writes output and error to stdout.
    """
    mol = ob.get_mol(s, make3D=True)
    if mult is None:
        mult = ob.get_multiplicity(mol)
    tmp = io.read_file(template)    
    inptext = get_qc_input(mol, tmp, mult)
    prefix = ob.get_smiles_filename(s)
    inpfile = prefix + '.nw'  
    outfile = prefix + '_nwchem.log'
    command = [exe, inpfile]
    if io.check_file(outfile, timeout=1):
        if overwrite:
            msg = 'Overwriting previous calculation "{0}"\n'.format(io.get_path(outfile))
            run = True
        else:
            msg = 'Skipping calculation, found "{0}"\n'.format(io.get_path(outfile))
            run = False
    else:
        run = True
    if run:
        if not io.check_file(inpfile, timeout=1) or overwrite:
            io.write_file(inptext, inpfile)
        if io.check_file(inpfile, timeout=1):
            msg = io.execute(command,stdoutfile=outfile,merge=True)
            if io.check_file(outfile, timeout=1):
                msg += ' Output file: "{0}"\n'.format(io.get_path(outfile))
        else:
            msg = 'Failed, cannot find input file "{0}"\n'.format(io.get_path(inpfile))
    return msg


def run_mopac(s, exe='mopac', method='pm7', mopackeys='precise nosym threads=1 opt', mult=0, overwrite=False):
    """
    Runs mopac calculation.
    If overwrite=False, no calculation if output file already exists.
    """
    mol = ob.get_mol(s, make3D=True)
    if mult == 0:
        mult = ob.get_multiplicity(mol)
    inptext = get_mopac_input(mol, method=method, keys=mopackeys, mult=mult, dothermo=True)
    prefix = ob.get_smiles_filename(s)
    inpfile = prefix + '.mop'
    outfile = prefix + '.out'
    if io.check_file(outfile, timeout=1):
        if overwrite:
            msg = 'Overwriting previous calculation "{0}"\n'.format(io.get_path(outfile))
            run = True
        else:
            msg = 'Skipping calculation, found "{0}"\n'.format(io.get_path(outfile))
            run = False
    else:
        run = True
    if run:
        if not io.check_file(inpfile, timeout=1):
            io.write_file(inptext, inpfile)
        if io.check_file(inpfile, timeout=1):
            msg = execute(inpfile, exe)
            if io.check_file(outfile, timeout=1):
                msg += ' Output file: "{0}"\n'.format(io.get_path(outfile))
        else:
            msg = ''
        if not io.check_file(outfile, timeout=1):
            msg = 'Failed, can not find "{0}"\n'.format(io.get_path(outfile))
    return msg


def execute_mopac(inp, exe='mopac'):
    """
    Runs mopac calculation.
    Mopac is a fortran code that does not return an error code
    but writes error to stderr.
    If there is no error stderr = None or ''
    For some keyword errors, still no stderr or stdout is provided,
    so additional checks are required.
    """
    from subprocess import Popen, PIPE
    import iotools as io
#    subprocess.call([mopacexe, inp])
    process = Popen([exe, inp], stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    if err is None or err == '':
        msg = 'Run {0} {1}: Success.'.format(exe, inp)
    else:
        errstr = """ERROR in {0}\n
        STDOUT:\n{1}\n
        STDERR:\n{2}""".format(inp, out, err)
        errfile = inp + '.err'
        io.write_file(errstr, errfile)
        msg = 'Run {0} {1}: Failed, see {2}\n'.format(exe, inp, errfile)
    return msg


def get_input_text(mol=None, s=None, template='qc_template.txt'):
    """
    Returns the text for quantum chemistry input file based on a template.
    """
    if mol is None:
        mol = ob.get_mol(s)
    inchikey = mol.write(format='inchikey').replace('\n', '')
    xyzlines = mol.write(format='xyz').splitlines(True)
    xyzstr = ''.join(xyzlines[2:])
    with open(template, 'r') as f:
        tmp = f.read()
    tmp = tmp.replace('QCPUT(GEO)\n', xyzstr)
    tmp = tmp.replace('QCPUT(UNIQUEKEY)', inchikey)
    tmp = tmp.replace('QCPUT(MULTIPLICITY)', str(mol.spin))
    tmp = tmp.replace('QCPUT(BASIS)', 'cc-pvdz')
    if mol.spin == 1:
        tmp = tmp.replace('QCPUT(MULTIPLICITY_WORD)', 'singlet')
        tmp = tmp.replace('QCPUT(SCFTYPE)', 'rhf')
    elif mol.spin == 2:
        tmp = tmp.replace('QCPUT(MULTIPLICITY_WORD)', 'doublet')
        tmp = tmp.replace('QCPUT(SCFTYPE)', 'uhf')
    elif mol.spin == 3:
        tmp = tmp.replace('QCPUT(MULTIPLICITY_WORD)', 'triplet')
        tmp = tmp.replace('QCPUT(SCFTYPE)', 'uhf')
    return tmp


def check_logfile(s):
    """
    Returns true if gaussian calculation completed succesfully
    """
    if "Normal termination of Gaussian" in s:
        return True
    elif "== MOPAC DONE ==" in s:
        return True
    else:
        return False

def get_qc_input(x, template, mult=None):
    """
    Returns Gaussian input file based on a given template.
    """
    mol = ob.get_mol(x)
    if mult is None:
        mult = ob.get_multiplicity(mol)
    nopen = mult - 1
    charge = ob.get_charge(mol)
    formula = ob.get_formula(mol)
    geo = ob.get_geo(mol)
    xyz = ob.get_xyz(mol)
    zmat = ob.get_zmat(mol)
    uniquekey = ob.get_unique_key(mol, mult)
    smilesname = ob.get_smiles_filename(mol)
    inp = template.replace("QTC(CHARGE)", str(charge))
    inp = inp.replace("QTC(MULTIPLICITY)", str(mult))
    inp = inp.replace("QTC(NOPEN)", str(nopen))
    inp = inp.replace("QTC(UNIQUEKEY)", uniquekey)
    inp = inp.replace("QTC(SMILESNAME)", smilesname)
    inp = inp.replace("QTC(ZMAT)", zmat)
    inp = inp.replace("QTC(GEO)", geo)
    inp = inp.replace("QTC(XYZ)", xyz)
    inp = inp.replace("QTC(FORMULA)", formula)
    if "QTC(" in inp:
        print("Error in template file:\n" + inp)
        return
    return inp


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
    uniquekey = ob.get_unique_key(mol, mult)
    inp = template.replace("QTC(CHARGE)", str(charge))
    inp = inp.replace("QTC(MULTIPLICITY)", str(mult))
    inp = inp.replace("QTC(UNIQUEKEY)", uniquekey)
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
    iline = io.get_line_number(key,lines=lines)
    iline += 1
    line = lines[iline]
    if iline < 0:
        return 'Not found: {0}'.format(key)
    while line.strip():
        cols = line.split()
        icol = int(cols[0])-1 
        for irow in range(icol,nfreq):
            iline += 1
            line = lines[iline] 
            cols = line.split()
            ncol = len(cols) - 1 
            xmat[irow,icol:icol+ncol] = [float(num) for num in cols[1:]]  
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
    if iline < 0:
        return 'Not found: {0}'.format(key)
    for i in range(nfreq):
        iline += 1
        line = lines[iline]
        cols = line.split()
        freqs[i,:] = [float(cols[1]),float(cols[2])]        
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
    if "Linear Molecule" in s:
        return True
    else:
        return False
    

def get_gaussian_nfreq(s):
    """
    Return the number of vibrational degrees of freedom for 
    a given log.
    """
    natom = get_gaussian_natom(s)
    if get_gaussian_islinear(s):
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


def print_list(s):
    return s

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

