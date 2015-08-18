#! /bin/bash

HFLUMI=40.03
PIXLUMI=40.03
LUMI=${HFLUMI}

root -l -q fitWm.C+\(\"Wmunu\",${LUMI},0\)
#root -l -q fitZmm.C+\(\"Zmumu\",${LUMI},0\)

#root -l -q fitWe.C+\(\"Wenu\",${LUMI},0\)
#root -l -q fitZee.C+\(\"Zee\",${LUMI},0\)
#root -l -q fitZee2.C+\(\"May23/Zee2\",${LUMI},0\)

#root -l -q plotZmm.C+\(\"Zmumu\",${LUMI}\)
#root -l -q  plotZee.C+\(\"Zee\",${LUMI}\)

#rm *.so *.d
