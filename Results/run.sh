#! /bin/bash

LUMI=2263

root -l -q  plotZmm.C+\(\"Zmm\",${LUMI}\)
root -l -q  plotZmmSystematics.C+\(\"Zmm\",${LUMI}\)
root -l -q  plotZmmCorrelations.C+\(\"Zmm\",${LUMI}\)
root -l -q  plotZmmStatCorrelations.C+\(\"Zmm\",${LUMI}\)
