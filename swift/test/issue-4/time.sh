#!/bin/bash -l

# Timing information is also available in the log file

INPUT=$1

if [[ ${#} != 1 ]]
then
  echo "Requires INPUT!"
  exit 1
fi

soft add +g09-e.01
time nice g09 $INPUT
