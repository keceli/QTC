Restart at: 0
(0 - beginning, 1 - reac/prod/ts.dat files already created, 
 2 - estoktp/theory files already created, 3 - estoktp level 0 successfully completed
 4 - estoktp level 1 successfully completed, 5 - estoktp 1D successfully completed, 
 6 - estoktp MD successfully completed or not desired, 7 - )
Parse all: True

============
Manual species input
============
Reaction type:  
Reactant list (SMILES):  QTC(SMILESNAME)
Product  list (SMILES): 

============
QTC options 
============
Use QTC xyz: False
             False - make openbabel xyz from smiles
             prog/method/basis - uses xyz in database (recommended)
             logfilename.log - makes torsscan parse out xyz from logfile in cwd
             True - xyz in testchem form in cwd

Use xyz as (start-starting geometry, 0-level0 geometry): start

==============
BLUES options   
==============
Run on node (0 if on login OR sshed on node, d to debug): 0
No. of cores high: 12
No. of cores  low:  6


================
EStoKTP options
================
No. MC  sampling points: 5
No. of transition state: 0
Scan interval (degrees): 360
No. of steps on the PES: 4

------------------------------------------------ 
   Module      :    Program    :      Theory
------------------------------------------------
Opt            :      g09      : b3lyp/6-31* 
Opt_1          :      g09      : m062x/6-311+g(d,p)
1dTau          :      g09      : m062x/6-311+g(d,p)
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


