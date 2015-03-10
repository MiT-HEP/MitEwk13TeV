#!/bin/bash

    EffType=$1
NoPseudoExp=$2

root -l -q -b lepEffSys.C+\(\"$EffType\", $NoPseudoExp\)

rm *so *d