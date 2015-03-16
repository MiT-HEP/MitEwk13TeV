#!/bin/bash

INFILE=/afs/cern.ch/work/j/jlawhorn/PYTHIA-GEN-SIM/Wme_bacon.root
PDF=CT10

echo $INFILE $PDF > wme_acc.txt

for i in `seq 0 44`;
do
    echo root -l -q accWme_unc.C+\(\"${INFILE}\",\"${PDF}\",${i}\)
    root -l -q accWme_unc.C+\(\"${INFILE}\",\"${PDF}\",${i}\) >> wme_acc.txt 2> screwyou.txt
done