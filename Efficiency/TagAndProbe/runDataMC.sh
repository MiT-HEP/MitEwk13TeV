# directory of tag-and-probe output
NTUPLEDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency-results

# integrated luminosity for data
LUMI=40.0

#
# Muon efficiencies - Data vs. MC
#

root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuHLTEff/eff.root\",\"${NTUPLEDIR}/DataZmm_MuHLTEff/eff.root\",\"muHLTEff\",0.01,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"trigger\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuHLTEff_pos/eff.root\",\"${NTUPLEDIR}/DataZmm_MuHLTEff_pos/eff.root\",\"muHLT_pos\",0.0,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"trigger\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuHLTEff_neg/eff.root\",\"${NTUPLEDIR}/DataZmm_MuHLTEff_neg/eff.root\",\"muHLT_neg\",0.0,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"trigger\"\)

root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuHLTEff_finepTbins/eff.root\",\"${NTUPLEDIR}/DataZmm_MuHLTEff_finepTbins/eff.root\",\"muHLTEff_finepTbins\",0.01,1.2,${LUMI},\"mufineptbins.bins\",\"Muon\",\"trigger\"\)

#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuSelEff/eff.root\",\"${NTUPLEDIR}/DataZmm_MuSelEff/eff.root\",\"muSelEff\",0.01,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"selection\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuSelEff_pos/eff.root\",\"${NTUPLEDIR}/DataZmm_MuSelEff_pos/eff.root\",\"muSel_pos\",0.0,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"trigger\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuSelEff_neg/eff.root\",\"${NTUPLEDIR}/DataZmm_MuSelEff_neg/eff.root\",\"muSel_neg\",0.0,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"trigger\"\)

#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuSelEff_finepTbins/eff.root\",\"${NTUPLEDIR}/DataZmm_MuSelEff_finepTbins/eff.root\",\"muSelEff_finepTbins\",0.01,1.2,${LUMI},\"mufineptbins.bins\",\"Muon\",\"selection\"\)

#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuTrkEff/eff.root\",\"${NTUPLEDIR}/DataZmm_MuTrkEff/eff.root\",\"muTrkEff\",0.01,1.2,${LUMI}\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuTrkEff_pos/eff.root\",\"${NTUPLEDIR}/DataZmm_MuTrkEff_pos/eff.root\",\"muTrk_pos\",0.0,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"tracking\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuTrkEff_neg/eff.root\",\"${NTUPLEDIR}/DataZmm_MuTrkEff_neg/eff.root\",\"muTrk_neg\",0.0,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"tracking\"\)

#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuTrkEff_finepTbins/eff.root\",\"${NTUPLEDIR}/DataZmm_MuTrkEff_finepTbins/eff.root\",\"muTrkEff_finepTbins\",0.01,1.2,${LUMI}\)

#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuStaEff_iso/eff.root\",\"${NTUPLEDIR}/DataZmm_MuStaEff_iso/eff.root\",\"muStaEff\",0.01,1.2,${LUMI}\)
#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuStaEff_iso_pos/eff.root\",\"${NTUPLEDIR}/DataZmm_MuStaEff_iso_pos/eff.root\",\"muSta_iso_pos\",0.0,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"stand-alone\"\)
#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuStaEff_iso_neg/eff.root\",\"${NTUPLEDIR}/DataZmm_MuStaEff_iso_neg/eff.root\",\"muSta_iso_neg\",0.0,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"stand-alone\"\)

#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuStaEff_iso_finepTbins/eff.root\",\"${NTUPLEDIR}/DataZmm_MuStaEff_iso_finepTbins/eff.root\",\"muStaEff_finepTbins\",0.01,1.2,${LUMI}\)

#
# Electron efficiencies - Data vs. MC
#

#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff/eff.root\",\"eleHLTEff\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)
#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_pos/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_pos/eff.root\",\"eleHLTEff_pos\",0.80,1.03,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)
#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_neg/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_neg/eff.root\",\"eleHLTEff_neg\",0.80,1.03,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)

#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_finepTbins/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_finepTbins/eff.root\",\"eleHLTEff_finepTbins\",0.01,1.3,${LUMI},\"elfineptbins.bins\",\"Supercluster\",\"trigger\"\)

#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleGsfSelEff/eff.root\",\"${NTUPLEDIR}/DataZee_EleGsfSelEff/eff.root\",\"eleGsfSelEff\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"GSF+ID+Iso\"\)
#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleGsfSelEff_pos/eff.root\",\"${NTUPLEDIR}/DataZee_EleGsfSelEff_pos/eff.root\",\"eleGsfSel_pos\",0.0,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"GSF+ID+Iso\"\)
#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleGsfSelEff_neg/eff.root\",\"${NTUPLEDIR}/DataZee_EleGsfSelEff_neg/eff.root\",\"eleGsfSel_neg\",0.0,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"GSF+ID+Iso\"\)

#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleGsfSelEff_finepTbins/eff.root\",\"${NTUPLEDIR}/DataZee_EleGsfSelEff_finepTbins/eff.root\",\"eleGsfSelEff_finepTbins\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"GSF+ID+Iso\"\)

#
# Debugging trigger efficiency
#
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/DataZee_EleHLTEff/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_afiq/eff.root\",\"eledebugHLT\",0.0,1.25,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)

#root -l -b -q plotAllRunsSingleRun.C+\(\"${NTUPLEDIR}/AllRunsSingleRun\",\"${NTUPLEDIR}/DataZee_EleHLTEff_oldruns_tightID/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_newruns_tightID/eff.root\",\"eleHLT\",0.2,1.1\)
