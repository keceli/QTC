#!/bin/sh
# set -x
printenv > env1.txt
time env -u PMI_FD -u PMI_RANK -u PMI_SIZE -u HYDI_CONTROL_FD ./g09.sh $*
