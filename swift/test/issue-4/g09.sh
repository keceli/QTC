#!/bin/sh
# set -x
printenv  > env2.txt

export LD_LIBRARY_PATH=/software/mvapich2-intel-psm-1.9.5/lib:/soft/intel/13.1.3/lib/intel64:/soft/intel/13.1.3/composer_xe_2013.5.192/mpirt/lib/intel64:/soft/intel/13.1.3/ipp/lib/intel64:/soft/intel/13.1.3/tbb/lib/intel64/gcc4.4:/soft/lcrc/lib:/usr/lib64:/usr/lib

g09 $*
