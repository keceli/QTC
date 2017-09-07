#!/bin/bash
set -eu

# RUN WORKFLOW
# Run this to run the Swift/T workflow

# Add Swift/T to PATH
PATH=/home/wozniak/Public/sfw/blues/compute/swift-t/stc/bin:$PATH

# Set up PYTHONPATH for QTC
PYTHONPATH+=/home/keceli/.local/lib/python2.7/site-packages:
PYTHONPATH+=/home/keceli/openbabel-2.4.1/install/lib/pys:
PYTHONPATH+=/home/keceli/openbabel-2.4.1/install/lib/python2.7/site-packages

# Run it!
nice swift-t -l workflow.swift --input=test/syngas.txt |& tee workflow.out
