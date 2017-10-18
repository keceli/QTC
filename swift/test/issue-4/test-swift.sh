#!/bin/bash -l
set -eu

PATH=/home/wozniak/Public/sfw/blues/compute/swift-t/stc/bin:$PATH
which swift-t
soft add +java-1.8
soft add +g09-e.01

set -x
time nice swift-t test.swift $*
