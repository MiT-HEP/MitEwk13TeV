# directory of tag-and-probe output
NTUPLEDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency-results

# integrated luminosity for data
LUMI=42.0

#
# Muon efficiencies - Data vs. MC
#

# trigger efficiency
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuHLTEff/eff.root\",\"${NTUPLEDIR}/DataZmm_MuHLTEff/eff.root\",\"muHLTEff\",0.01,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"trigger\"\)

# tracking efficiency - N_{global} / N_{sta} (use typeBits attribute of baconhep::TMuons)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuTrkEff/eff.root\",\"${NTUPLEDIR}/DataZmm_MuTrkEff/eff.root\",\"muTrkEff\",0.01,1.3,${LUMI},\"mupteta.bins\",\"Muon\",\"tracking\"\)
# tracking efficiency v2 - N_{sta matched to tracks} / N_{sta} (use trkID attribute of baconhep::TMuons)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuTrkEff_v2/eff.root\",\"${NTUPLEDIR}/DataZmm_MuTrkEff_v2/eff.root\",\"muTrkEff_v2\",0.01,1.3,${LUMI},\"mupteta.bins\",\"Muon\",\"tracking\"\)

# stand-alone+id+iso efficiency
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuSelStaEff_iso/eff.root\",\"${NTUPLEDIR}/DataZmm_MuSelStaEff_iso/eff.root\",\"muSelStaEff\",0.01,1.3,${LUMI},\"mupteta.bins\",\"Muon\",\"standalone+ID+Iso\"\)
# compare sf(sel)*sf(sta) with sf(sel+sta)
#root -l -b -q muSelStaEff_compareSFs.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuSelStaEff_iso/eff.root\",\"${NTUPLEDIR}/DataZmm_MuSelStaEff_iso/eff.root\",\"muSelStaEff_ratios\",0.01,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"standalone+ID+Iso\"\)

# id+iso efficiency
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuSelEff/eff.root\",\"${NTUPLEDIR}/DataZmm_MuSelEff/eff.root\",\"muSelEff\",0.01,1.2,${LUMI},\"mupteta.bins\",\"Muon\",\"ID+Iso\"\)

# stand-alone efficiency
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zmm_MuStaEff_iso/eff.root\",\"${NTUPLEDIR}/DataZmm_MuStaEff_iso/eff.root\",\"muStaEff\",0.01,1.3,${LUMI},\"mupteta.bins\",\"Muon\",\"standalone\"\)

#
# Electron efficiencies - Data vs. MC
#

# trigger efficiency
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff/eff.root\",\"eleHLTEff\",0.01,1.2,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)

# with Kevin's binning
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_binsFromKevin/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_binsFromKevin/eff.root\",\"eleHLTEff_binsFromKevin\",0.01,1.2,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)
# for specific runs
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_run251244/eff.root\",\"eleHLTEff_run251244\",0.01,1.2,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_run251251/eff.root\",\"eleHLTEff_run251251\",0.01,1.2,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_run251252/eff.root\",\"eleHLTEff_run251252\",0.01,1.2,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_run251561/eff.root\",\"eleHLTEff_run251561\",0.01,1.2,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_run251562/eff.root\",\"eleHLTEff_run251562\",0.01,1.2,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_run251643/eff.root\",\"eleHLTEff_run251643\",0.01,1.2,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_run251721/eff.root\",\"eleHLTEff_run251721\",0.01,1.2,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleHLTEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleHLTEff_run251883/eff.root\",\"eleHLTEff_run251883\",0.01,1.2,${LUMI},\"elpteta.bins\",\"Supercluster\",\"trigger\"\)

# reco+id+iso efficiency
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleGsfSelEff/eff.root\",\"${NTUPLEDIR}/DataZee_EleGsfSelEff/eff.root\",\"eleGsfSelEff\",0.01,1.2,${LUMI},\"elpteta.bins\",\"Supercluster\",\"GSF+ID+Iso\"\)

# supercluster efficiency
root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSCEff/eff.root\",\"${NTUPLEDIR}/DataZee_EleSCEff/eff.root\",\"eleSCEff\",0.01,1.3,${LUMI},\"elsupercluster.bins\",\"Supercluster\",\"supercluster\",1\)

# id+iso efficiency
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSelEff/eff.root\",\"${NTUPLEDIR}/DataZee_EleSelEff/eff.root\",\"eleSelEff\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"ID+Iso\"\)

# with Kevin's binning
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSelEff_binsFromKevin/eff.root\",\"${NTUPLEDIR}/DataZee_EleSelEff_binsFromKevin/eff.root\",\"eleSelEff_binsFromKevin\",0.01,1.3,${LUMI},\"elsel.bins\",\"Supercluster\",\"ID+Iso\"\)
# with tight id
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSelEff_tightID/eff.root\",\"${NTUPLEDIR}/DataZee_EleSelEff_tightID/eff.root\",\"eleSelEff_tightID\",0.01,1.3,${LUMI},\"elsel.bins\",\"Supercluster\",\"ID+Iso\"\)
# for specific runs
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSelEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleSelEff_run251244/eff.root\",\"eleSelEff_run251244\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"ID+Iso\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSelEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleSelEff_run251251/eff.root\",\"eleSelEff_run251251\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"ID+Iso\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSelEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleSelEff_run251252/eff.root\",\"eleSelEff_run251252\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"ID+Iso\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSelEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleSelEff_run251561/eff.root\",\"eleSelEff_run251561\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"ID+Iso\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSelEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleSelEff_run251562/eff.root\",\"eleSelEff_run251562\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"ID+Iso\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSelEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleSelEff_run251643/eff.root\",\"eleSelEff_run251643\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"ID+Iso\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSelEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleSelEff_run251721/eff.root\",\"eleSelEff_run251721\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"ID+Iso\"\)
#root -l -b -q plotDataMC_singlepTbins.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleSelEff_singleRunBins/eff.root\",\"${NTUPLEDIR}/DataZee_EleSelEff_run251883/eff.root\",\"eleSelEff_run251883\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"ID+Iso\"\)

# gsf efficiency
#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleGsfEff/eff.root\",\"${NTUPLEDIR}/DataZee_EleGsfEff/eff.root\",\"eleGsfEff\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"GSF\"\)

# gsf efficiency with deltaR matching
#root -l -b -q plotDataMC_withRatios.C+\(\"${NTUPLEDIR}/DataMC\",\"${NTUPLEDIR}/Zee_EleGsfEff_dRmatching/eff.root\",\"${NTUPLEDIR}/DataZee_EleGsfEff_dRmatching/eff.root\",\"eleGsfEff_dRmatching\",0.01,1.3,${LUMI},\"elpteta.bins\",\"Supercluster\",\"GSF\"\)

