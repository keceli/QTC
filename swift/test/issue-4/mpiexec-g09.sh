#!/bin/sh
# set -x
printenv  > env2.txt
mpiexec -n 1 g09 $*
