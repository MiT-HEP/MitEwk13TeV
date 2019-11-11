#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

OUTDIR="/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/probes/"
 
# format: makeGenMllTemplate*.C+("conf with file location","outputDirectory","outputFileName")
# root -l -b -q makeGenMllTemplate_forElectron_v2.C+\(\"zee_13.conf\",\"${OUTDIR}/GenReweightMassHist/EleGSF\",\"zee-amc-pyth-mll.root\"\)
# root -l -b -q makeGenMllTemplate_forElectron_v2.C+\(\"zee_PowPhot.conf\",\"${OUTDIR}/GenReweightMassHist/EleGSF\",\"zee-pow-phot-mll.root\"\)
root -l -b -q makeGenMllTemplate_forElectron_v2.C+\(\"zee_PowPyth.conf\",\"${OUTDIR}/GenReweightMassHist/EleGSF\",\"zee-pow-pyth-mll.root\"\)

 # root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_13.conf\",\"${OUTDIR}/GenReweightMassHist/MuSta\",\"zmm-amc-pyth-mll.root\"\)
 # root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_PowPhot.conf\",\"${OUTDIR}/GenReweightMassHist/MuSta\",\"zmm-pow-phot-mll.root\"\)
 # root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_PowPyth.conf\",\"${OUTDIR}/GenReweightMassHist/MuSta\",\"zmm-pow-pyth-mll.root\"\)
 

 
  # root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_13.conf\",\"${OUTDIR}/GenReweightMassHist/MuSIT\",\"zmm-amc-pyth-mll.root\"\)
 # root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_PowPhot.conf\",\"${OUTDIR}/GenReweightMassHist/MuSIT\",\"zmm-pow-phot-mll.root\"\)
 # root -l -b -q makeGenMllTemplate_forMuon_v2.C+\(\"zmm_PowPyth.conf\",\"${OUTDIR}/GenReweightMassHist/MuSIT\",\"zmm-pow-pyth-mll.root\"\)

# rm *.so *.d
