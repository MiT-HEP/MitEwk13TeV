#! /bin/bash

HFLUMI=41.8556
PIXLUMI=41.8556
LUMI=${HFLUMI}

#root -l -q fitWm.C+\(\"Wmunu\",${LUMI},0\)
#root -l -q fitWm_mc.C+\(\"Wmunu_mc\",${LUMI},0\)
#root -l -q fitZmm.C+\(\"Zmumu\",${LUMI},0\)

root -l -q fitWe.C+\(\"Wenu\",${LUMI},0\)
#root -l -q fitZee.C+\(\"Zee\",${LUMI},0\)
#root -l -q fitZee2.C+\(\"May23/Zee2\",${LUMI},0\)

#root -l -q  plotZmm.C+\(\"Zmumu\",${LUMI}\)
#root -l -q   plotZee.C+\(\"Zee\",${LUMI}\)

#root -l -q postFitWe.C+\(\"Wenu_post\",${LUMI},0\)
#root -l -q postFitWm.C+\(\"Wmunu_post\",${LUMI},0\)
#rm *.so *.d
