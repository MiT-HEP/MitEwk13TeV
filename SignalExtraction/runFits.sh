#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=96
LUMI2=96

# root -l -q fitWm.C+\(\"Wmunu\",${LUMI},0\)
# root -l -q fitWm.C+\(\"WmunuPuppi_WShapeDiff_testBkg\",${LUMI},${LUMI2},0\)
# root -l -q fitWm.C+\(\"AntiWmunu_WshapeDiff2\",${LUMI},${LUMI2},0\)
# root -l -q fitBinsWm.C+\(\"Wmunu_etaBins_wShapes\",${LUMI},${LUMI2},0\)
# root -l -q fitZm.C+\(\"Zmumu_newBacon_lowPt\",${LUMI},0\)
#root -l -q fitWm_mc.C+\(\"Wmunu_mc\",${LUMI},0\)
# root -l -q fitZm.C+\(\"Zmumu_PF_testToys\",${LUMI},0\)

# root -l -q fitWe.C+\(\"Wenu_wShapes\",${LUMI},0\)
#root -l -q fitZee.C+\(\"Zee\",${LUMI},0\)
#root -l -q fitZee2.C+\(\"May23/Zee2\",${LUMI},0\)

# root -l -q  plotZmm.C+\(\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\,0\)
# root -l -q  plotZmm.C+\(\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\,1\)
# root -l -q  plotZmmResScaleUncert.C+\(\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\)
#root -l -q   plotZee.C+\(\"Zee\",${LUMI}\)

# root -l -q postFitWe.C+\(\"Wenu_post_Wshapes\",${LUMI},0\)
root -l -q postFitWm.C+\(\"Wmunu_post_etaBins3\",${LUMI},0\)
# root -l -q postFitWm.C+\(\"WmunuPF_post_WShapeDiff\",${LUMI},0\)
#rm *.so *.d
