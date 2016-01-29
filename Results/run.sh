#! /bin/bash

LUMI=2215

root -l -q  plotZmm.C+\(\"Zmumu\",${LUMI}\)
root -l -q  plotZmmSystematics.C+\(\"Zmumu\",${LUMI}\)
root -l -q  plotZmmCorrelations.C+\(\"Zmumu\",${LUMI}\)
root -l -q  plotZmmStatCorrelations.C+\(\"Zmumu\",${LUMI}\)
