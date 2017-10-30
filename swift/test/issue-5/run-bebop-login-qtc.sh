#!/bin/bash -l
set -eu

# Set P to argument 1 if given, else use default of 2
P=${1:-2}

module load intel/17.0.4-74uvhji
module load jdk/8u141-b15-mopj6qr

PATH=$HOME/Public/sfw/bebop/compute/swift-t/stc/bin:$PATH

# MLS = MPI Launch for Swift
export MLS=/home/wozniak/Public/sfw/bebop/mls/src

swift-t -n 10 -I $MLS -r $MLS workflow-qtc.swift -P=$P
