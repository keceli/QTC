#!/bin/sh

ID=$1
echo ID=$ID RANK=$PMI_RANK SIZE=$PMI_SIZE
echo "[$ID]" "hosts file: ${swift_write_hosts}"
cat ${swift_write_hosts} | sed "s/^/[$ID] /"
