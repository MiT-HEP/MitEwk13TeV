# directory of tag-and-probe output
NTUPLEDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency-results

# integrated luminosity for data
LUMI=7.3

#
# Muon efficiencies - Data vs. MC
#

#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuHLTEff/eff.root\",\"${NTUPLEDIR}/DataZmm_MuHLTEff/eff.root\",\"muHLT\",0.0,1.2,${LUMI},\"musel.bins\",\"Muon\",\"trigger\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuHLTEff_pos/eff.root\",\"${NTUPLEDIR}/DataZmm_MuHLTEff_pos/eff.root\",\"muHLT_pos\",0.0,1.2,${LUMI},\"muhlt.bins\",\"Muon\",\"trigger\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuHLTEff_neg/eff.root\",\"${NTUPLEDIR}/DataZmm_MuHLTEff_neg/eff.root\",\"muHLT_neg\",0.0,1.2,${LUMI},\"muhlt.bins\",\"Muon\",\"trigger\"\)

#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuHLTEff_finepTbins/eff.root\",\"${NTUPLEDIR}/DataZmm_MuHLTEff_finepTbins/eff.root\",\"muHLT\",0.0,1.2,${LUMI},\"musel.bins\",\"Muon\",\"trigger\"\)

#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuSelEff/eff.root\",\"${NTUPLEDIR}/DataZmm_MuSelEff/eff.root\",\"muSel\",0.0,1.2,${LUMI},\"muhlt.bins\",\"Muon\",\"selection\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuSelEff_pos/eff.root\",\"${NTUPLEDIR}/DataZmm_MuSelEff_pos/eff.root\",\"muSel_pos\",0.0,1.2,${LUMI},\"muhlt.bins\",\"Muon\",\"trigger\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuSelEff_neg/eff.root\",\"${NTUPLEDIR}/DataZmm_MuSelEff_neg/eff.root\",\"muSel_neg\",0.0,1.2,${LUMI},\"muhlt.bins\",\"Muon\",\"trigger\"\)

#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuSelEff_finepTbins/eff.root\",\"${NTUPLEDIR}/DataZmm_MuSelEff_finepTbins/eff.root\",\"muSel\",0.0,1.2,${LUMI},\"muhlt.bins\",\"Muon\",\"selection\"\)

#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuTrkEff/eff.root\",\"${NTUPLEDIR}/DataZmm_MuTrkEff/eff.root\",\"muTrk\",0.0,1.2,${LUMI}\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuTrkEff_pos/eff.root\",\"${NTUPLEDIR}/DataZmm_MuTrkEff_pos/eff.root\",\"muTrk_pos\",0.0,1.2,${LUMI},\"muhlt.bins\",\"Muon\",\"tracking\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuTrkEff_neg/eff.root\",\"${NTUPLEDIR}/DataZmm_MuTrkEff_neg/eff.root\",\"muTrk_neg\",0.0,1.2,${LUMI},\"muhlt.bins\",\"Muon\",\"tracking\"\)

#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuStaEff_iso/eff.root\",\"${NTUPLEDIR}/DataZmm_MuStaEff_iso/eff.root\",\"muSta\",0.0,1.2,${LUMI}\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuStaEff_iso_pos/eff.root\",\"${NTUPLEDIR}/DataZmm_MuStaEff_iso_pos/eff.root\",\"muSta_iso_pos\",0.0,1.2,${LUMI},\"muhlt.bins\",\"Muon\",\"stand-alone\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuStaEff_iso_neg/eff.root\",\"${NTUPLEDIR}/DataZmm_MuStaEff_iso_neg/eff.root\",\"muSta_iso_neg\",0.0,1.2,${LUMI},\"muhlt.bins\",\"Muon\",\"stand-alone\"\)

#
# Electron efficiencies - Data vs. MC
#

#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff/eff.root\",\"eleHLT\",0.0,1.3,${LUMI},\"elgsfsel.bins\",\"Supercluster\",\"trigger\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_pos/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_pos/eff.root\",\"eleHLT_pos\",0.80,1.03,${LUMI},\"elgsfsel.bins\",\"Supercluster\",\"trigger\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_neg/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_neg/eff.root\",\"eleHLT_neg\",0.80,1.03,${LUMI},\"elgsfsel.bins\",\"Supercluster\",\"trigger\"\)

#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_finepTbins/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_finepTbins/eff.root\",\"eleHLT\",0.0,1.3,${LUMI},\"elgsfsel.bins\",\"Supercluster\",\"trigger\"\)

#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleGsfSelEff/eff.root\",\"${NTUPLEDIR}/DataZee_EleGsfSelEff/eff.root\",\"eleGsfSel\",0.0,1.3,${LUMI},\"elgsfsel.bins\",\"Supercluster\",\"GSF+ID+Iso\"\)
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleGsfSelEff_pos/eff.root\",\"${NTUPLEDIR}/DataZee_EleGsfSelEff_pos/eff.root\",\"eleGsfSel_pos\",0.0,1.3,${LUMI},\"elgsfsel.bins\",\"Supercluster\",\"GSF+ID+Iso\"\)
root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleGsfSelEff_neg/eff.root\",\"${NTUPLEDIR}/DataZee_EleGsfSelEff_neg/eff.root\",\"eleGsfSel_neg\",0.0,1.3,${LUMI},\"elgsfsel.bins\",\"Supercluster\",\"GSF+ID+Iso\"\)

#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleGsfSelEff_finepTbins/eff.root\",\"${NTUPLEDIR}/DataZee_EleGsfSelEff_finepTbins/eff.root\",\"eleGsfSel\",0.0,1.3,${LUMI},\"elgsfsel.bins\",\"Supercluster\",\"GSF+ID+Iso\"\)

#
# Debugging trigger efficiency
#
#root -l -b -q plotDataMC.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/DataZee_EleHLTEff/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_afiq/eff.root\",\"eledebugHLT\",0.0,1.25,${LUMI},\"elgsfsel.bins\",\"Supercluster\",\"trigger\"\)

#root -l -b -q plotAllRunsSingleRun.C+\(\"${NTUPLEDIR}/AllRunsSingleRun\",\"${NTUPLEDIR}/DataZee_EleHLTEff_oldruns_tightID/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_newruns_tightID/eff.root\",\"eleHLT\",0.2,1.1\)
