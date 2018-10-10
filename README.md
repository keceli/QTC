# Quantum Thermochemistry Calculator 
[![Build Status](https://travis-ci.org/keceli/QTC.svg?branch=master)](https://travis-ci.org/keceli/QTC)
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/keceli/QTC/master?filepath=https%3A%2F%2Fgithub.com%2Fkeceli%2FQTC%2Fblob%2Fmaster%2Fnotebooks%2Ftutorial.ipynb)

QTC includes modules that integrates open babel with quantum chemistry calculations and generates NASA polynomials in different formats.

It depends on:
  * [Open Babel](http://openbabel.org/) for cheminformatics 
  * MOPAC, NWChem, Gaussian, Molpro for quantum chemistry calculations
  * [MESS](https://github.com/PACChem/MESS) for calculating partititon function
  * [RMG](https://github.com/ReactionMechanismGenerator/RMG-Py) for generating a species list important for combustion chemistry
  * PAC99, thermp for format conversions

## Installation
conda simplifies the installation process
If you don't have it you can install it without root privilages
```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
```
After the installation is completed you can install the dependencies with:
```
conda env create -f environment.yml
```
Run qtc after activating qtc-env environment, i.e.
```
source activate qtc-env
```
Alternatively, you can install them with conda:
```
conda install numpy psutil
conda install -c openbabel openbabel
conda install -c mcs07 cirpy 
```
## Run QTC
```
python src/qtc.py -i O -k 'opt/mp2/dz/nwchem' -Q
```

## Acknowledgment

This work was supported by the U.S. Department of Energy, Office of Basic Energy
Sciences, Division of Chemical Sciences, Geosciences, and Biosciences under DOE
Contract Number DE-AC02-06CH11357 as well as the Exascale Computing Project
(ECP), Project Number: 17-SC-20-SC.  The ECP is a collaborative effort of two
DOE organizations, the Office of Science and the National Nuclear Security
Administration, responsible for the planning and preparation of a capable
exascale ecosystem including software, applications, hardware, advanced system
engineering, and early test bed platforms to support the nation's exascale
computing imperative. 

