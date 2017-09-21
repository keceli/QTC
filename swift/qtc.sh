#!/bin/bash -l
set -eu

#echo QTC.SH

#echo $PATH

if [[ ${#} != 2 ]]
then
  echo "qtc.sh: Requires 2 arguments!"
  exit 1
fi

MOLECULE=$1
INDEX=$2

soft add @nwchem-6.5
soft add +molpro-2012.1.27
soft add +g09-e.01

# Set up PYTHONPATH for QTC
PYTHONPATH+=:
PYTHONPATH+=/home/keceli/.local/lib/python2.7/site-packages:
PYTHONPATH+=/home/keceli/openbabel-2.4.1/install/lib/pys:
PYTHONPATH+=/home/keceli/openbabel-2.4.1/install/lib/python2.7/site-packages
export PYTHONPATH

#  /home/keceli/backup/QTC/qtc.py

set -x
MOLPRO=$( which molpro )
mkdir -pv /scratch/$USER

# printenv | sort

#env -u PMI_FD -u PMI_RANK -u PMI_SIZE -u HYDI_CONTROL_FD  \
#    $MOLPRO -d /scratch/wozniak "/home/wozniak/N#N_molpro.inp"

#/home/keceli/anaconda2/bin/python -u /home/keceli/qtc/qtc.py \
/home/keceli/anaconda2/bin/python -u /home/keceli/pacc/codes/QTC/qtc.py \
                                   -i $MOLECULE \
                                   -k 'torsopt/hf/sto-3g/gaussian' \
                                    -p 16 \
                                    -t 'templates' \
                                    -f "swift_torsopt-$INDEX"
    # -k energy/mp2/sto-3g/gaussian
    # -k energy/mp2/sto-3g/nwchem -p 2
    # Use +g09-e.01
    # Use +molpro-2012.1.27
    #-k 'torsscan/b2plypd3/cc-pvtz/gaussian,energy/ccsd(t)-f12/cc-pvdz-f12/molpro,energy/mp2-f12/cc-pvtz-f12/molpro,energy/mp2-f12/cc-pvdz-f12/molpro,composite/cbs_f12_dztz/energy=e[1]+e[2]-e[3]' \
