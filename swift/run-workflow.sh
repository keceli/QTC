#!/bin/bash -l
set -eu

# RUN WORKFLOW
# Run this to run the Swift/T workflow

if [[ ${#} != 2 ]]
then
  echo "run-workflow.sh: Requires MODE INPUT!"
  sleep 5
  exit 0
fi

export THIS=$( cd $( dirname $0 ) ; /bin/pwd )

MODE=$1
INPUT=$2

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
swift-t -l -n 8 $MACHINE $THIS/workflow.swift --input=$INPUT
# sorted.txt
#|& tee workflow.out
