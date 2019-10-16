#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================
# format: makeGenMllTemplate*.C+("conf with file location","outputDirectory","outputFileName")
# root -l -b -q makeGenMllTemplate_forElectron_v2.C+\(\"zee_13.conf\",\"GenReweightMassHist/EleGSF\",\"zee-amc-pyth-mll.root\"\)
# root -l -b -q makeGenMllTemplate_forElectron_v2.C+\(\"zee_PowPhot.conf\",\"GenReweightMassHist/EleGSF\",\"zee-pow-phot-mll.root\"\)
# root -l -b -q makeGenMllTemplate_forElectron_v2.C+\(\"zee_PowPyth.conf\",\"GenReweightMassHist/EleGSF\",\"zee-pow-pyth-mll.root\"\)

 # root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_13.conf\",\"GenReweightMassHist/MuSta\",\"zmm-amc-pyth-mll.root\"\)
 # root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_PowPhot.conf\",\"GenReweightMassHist/MuSta\",\"zmm-pow-phot-mll.root\"\)
 # root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_PowPyth.conf\",\"GenReweightMassHist/MuSta\",\"zmm-pow-pyth-mll.root\"\)
 
  root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_13.conf\",\"GenReweightMassHist/MuSIT\",\"zmm-amc-pyth-mll.root\"\)
 root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_PowPhot.conf\",\"GenReweightMassHist/MuSIT\",\"zmm-pow-phot-mll.root\"\)
 root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_PowPyth.conf\",\"GenReweightMassHist/MuSIT\",\"zmm-pow-pyth-mll.root\"\)

# rm *.so *.d
