#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=2157.8 # updated removed the bad data in the lumiBin#6
LUMI2=90.5 # this number need to be updated (muon CR is prescaled, electron CR is full lumi)

root -l -q fitWm.C+\(\"Wmunu_default\",${LUMI},${LUMI2},0\) &
sleep 100
root -l -q fitWe.C+\(\"Wenu_default\",${LUMI},${LUMI},0\) #central

#rm *.so *.d
