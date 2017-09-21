#!/bin/bash -l
set -eu

# RUN WORKFLOW
# Run this to run the Swift/T workflow

if [[ ${#} != 4 ]]
then
  echo "run-workflow.sh: Requires MODE PROCS CFG INPUT!"
  sleep 5
  exit 0
fi

MODE=$1
PROCS=$2
CFG=$3
INPUT=$4

export THIS=$( cd $( dirname $0 ) ; /bin/pwd )

case $MODE in
  login) 
    MACHINE=""
    ;;
  tcg)
    MACHINE="-t f:hosts.txt"
    ;;
  blues)
    MACHINE="-m pbs"
    ;;
  default)
    echo "unknown MODE=$MODE"
    exit 1
    ;;
esac

# Add Swift/T to PATH
PATH=/home/wozniak/Public/sfw/blues/compute/swift-t/stc/bin:$PATH

soft add +java-1.8

# Run it!
# /home/keceli/qtc/swift/
# -t f:hosts.txt
nice swift-t -l -n $PROCS $MACHINE $THIS/workflow.swift \
        --input=$INPUT --cfg=$CFG
# sorted.txt
#|& tee workflow.out
