#! /bin/bash

LUMI=2318.3

root -l -q  plotZee.C+\(\"Zee\",${LUMI}\)
#root -l -q  plotZeeSystematics.C+\(\"Zee\",${LUMI}\)
#root -l -q  plotZeeCorrelations.C+\(\"Zee\",${LUMI}\)
#root -l -q  plotZeeStatCorrelations.C+\(\"Zee\",${LUMI}\)
