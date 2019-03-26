#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================
# format: makeGenMllTemplate*.C+("conf with file location","outputDirectory","outputFileName")
root -l -b -q makeGenMllTemplate_forElectron_v2.C+\(\"zee_13.conf\",\"GenReweightMassHist_v2\",\"zee-amc-pyth-mll.root\"\)
# root -l -b -q makeGenMllTemplate_forElectron_v2.C+\(\"zee_PowPhot.conf\",\"GenReweightMassHist_v2\",\"zee-pow-phot-mll.root\"\)
# root -l -b -q makeGenMllTemplate_forElectron_v2.C+\(\"zee_PowPyth.conf\",\"GenReweightMassHist_v2\",\"zee-pow-pyth-mll.root\"\)

# root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_13.conf\",\"GenReweightMassHist\",\"zmm-amc-pyth-mll.root\"\)
# root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_PowPhot.conf\",\"GenReweightMassHist\",\"zmm-pow-phot-mll.root\"\)
# root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_PowPyth.conf\",\"GenReweightMassHist\",\"zmm-pow-pyth-mll.root\"\)

# rm *.so *.d
