#!/bin/bash -l
set -eu

# Set P to argument 1 if given, else use default of 2
P=${1:-2}

soft add +gcc-5.3.0

PATH=$HOME/Public/sfw/blues/compute/swift-t/pacc/stc/bin:$PATH

# MLS = MPI Launch for Swift
export MLS=/home/wozniak/Public/sfw/blues/mls/src

swift-t -n 10 -I $MLS -r $MLS workflow.swift -P=$P
