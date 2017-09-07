#!/bin/bash
set -eu

# RUN WORKFLOW BLUES
# Run this to run the Swift/T workflow on Blues

# Add Swift/T to PATH
PATH=/home/wozniak/Public/sfw/blues/compute/swift-t/stc/bin:$PATH

# Set up PYTHONPATH for QTC
PYTHONPATH+=/home/keceli/.local/lib/python2.7/site-packages:
PYTHONPATH+=/home/keceli/openbabel-2.4.1/install/lib/pys:
PYTHONPATH+=/home/keceli/openbabel-2.4.1/install/lib/python2.7/site-packages

# Run it!
swift-t -l -e PYTHONPATH \
        -m pbs \
        workflow.swift --input=$PWD/test/syngas.txt |& tee workflow.out
