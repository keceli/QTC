#!/bin/bash -l
set -eu

# RUN WORKFLOW
# Run this to run the Swift/T workflow

if [[ ${#} != 1 ]]
then
  echo "run-workflow.sh: Requires MODE!"
  sleep 5
  exit 0
fi

MODE=$1
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
swift-t -l -n 8 $MACHINE workflow.swift --input=test/syngas.txt
# sorted.txt
#|& tee workflow.out
