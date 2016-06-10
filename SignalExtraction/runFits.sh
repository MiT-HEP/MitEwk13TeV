#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=2318

#root -l -q fitWm.C+\(\"Wmunu\",${LUMI},0\)
#root -l -q fitWm.C+\(\"Wmunu_newBaconCorr\",${LUMI},0\)
# root -l -q fitZm.C+\(\"Zmumu_newBacon_lowPt\",${LUMI},0\)
#root -l -q fitWm_mc.C+\(\"Wmunu_mc\",${LUMI},0\)
root -l -q fitZm.C+\(\"Zmumu_newBacon_lowPt\",${LUMI},0\)

#root -l -q fitWe.C+\(\"Wenu\",${LUMI},0\)
#root -l -q fitZee.C+\(\"Zee\",${LUMI},0\)
#root -l -q fitZee2.C+\(\"May23/Zee2\",${LUMI},0\)

# root -l -q  plotZmm.C+\(\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\,0\)
# root -l -q  plotZmm.C+\(\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\,1\)
# root -l -q  plotZmmResScaleUncert.C+\(\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\)
#root -l -q   plotZee.C+\(\"Zee\",${LUMI}\)

#root -l -q postFitWe.C+\(\"Wenu_post\",${LUMI},0\)
#root -l -q postFitWm.C+\(\"Wmunu_post\",${LUMI},0\)
#rm *.so *.d
