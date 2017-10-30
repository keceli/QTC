#!/bin/sh

# QTC TASK WRAPPER

ID=$1
echo ID=$ID RANK=$PMI_RANK SIZE=$PMI_SIZE
echo "[$ID]" "hosts file: ${swift_write_hosts}"

PATH=/home/keceli/anaconda2/bin:$HOME/proj/qtc:/home/keceli/bin:$PATH
export PYTHONPATH=/home/keceli/.local/lib/python2.7/site-packages:/home/keceli/bebop/lib/openbabel-2.4.1/install/lib/python2.7/site-packages:/home/keceli/thermo/code/qtc

if [ ${swift_write_hosts:-NOTHING} = NOTHING ]
then
  echo "No hosts file provided!"
  exit 1
fi

qtc.py $* -m $swift_write_hosts
