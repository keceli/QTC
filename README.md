# Quantum Thermo Chemistry (QTC) 
[![Build Status](https://travis-ci.org/keceli/QTC.svg?branch=master)](https://travis-ci.org/keceli/QTC)

QTC includes modules that integrates open babel with quantum chemistry calculations and generates NASA polynomials in different formats.

It depends on:
  * Open Babel for cheminformatics (read/write chemical identifiers)
  * MOPAC, NWChem, Gaussian, Molpro for quantum chemistry calculations
  * MESS for calculating partititon function
  * RMG for generating a species list important for combustion chemistry
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
conda install numpy psutil
conda install -c openbabel openbabel
conda install -c mcs07 cirpy 

```
You may prefer to create a virtual environment first to avoid conflicts with other modules
