Restart at: 0
(0 - beginning, 1 - reac/prod/ts.dat files already created, 
 2 - estoktp/theory files already created, 3 - estoktp level 0 successfully completed
 4 - estoktp level 1 successfully completed, 5 - estoktp 1D successfully completed, 
 6 - estoktp MD successfully completed or not desired, 7 - )

============
Manual species input
============
Reaction type:  
Reactant list (SMILES):  QTC(SMILESNAME)
Product  list (SMILES): 
No. of transition states: 0

============
QTC options 
============
Use QTC xyz: QTC(METHOD)/QTC(BASIS)/QTC(TASK)
             False - make openbabel xyz from smiles
             prog/method/basis - uses xyz in database (recommended)
             logfilename.log - makes torsscan parse out xyz from logfile in cwd
             True - xyz in testchem form in cwd

Use xyz as (start-starting geometry, 0-level0 geometry): 0

==============
BLUES options   
==============
Run on node (0 if on login OR on a node, d to debug): 0
No. of cores high: 24
No. of cores  low:  16


================
EStoKTP options
================
No. MC  sampling points: 5
Scan interval (degrees): 360
No. of steps on the PES: 4

------------------------------------------------ 
   Module      :    Program    :      Theory
------------------------------------------------
Opt            :               : 
Opt_1          : QTC(METHOD)   : QTC(BASIS)/QTC(TASK)
1dTau          : QTC(METHOD)   : QTC(BASIS)/QTC(TASK)
MdTau          :               :
Symm           :               :
HL             :               :
Irc            :               :
------------------------------------------------

=============
THERMO options
=============
Perform all thermochemistry? (default true): true
Anharmonic (0 for level0 theory, 1 for level1 theory, false for off): false
Overwrite anharmonic: false
Basis for heat of formation: auto 
Parse all: True


