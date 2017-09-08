#!/bin/bash
set -eu

echo QTC.SH

echo $PATH

soft add @nwchem-6.5
soft add +molpro-2012.1.27
soft add +g09-e.01

MOLECULE=$1

#  /home/keceli/backup/QTC/qtc.py

set -x
MOLPRO=$( which molpro )
mkdir -pv /scratch/wozniak

# printenv | sort

env -u PMI_FD -u PMI_RANK -u PMI_SIZE -u HYDI_CONTROL_FD  \
    $MOLPRO -d /scratch/wozniak "/home/wozniak/N#N_molpro.inp"

# /home/keceli/anaconda2/bin/python -u /home/keceli/qtc/qtc.py \
#                                   -i $MOLECULE \
#                                   -k energy/mp2/sto-3g/molpro \
#                                   -p 2
    # -k energy/mp2/sto-3g/gaussian
    # -k energy/mp2/sto-3g/nwchem -p 2
    # Use +g09-e.01
    # Use +molpro-2012.1.27
