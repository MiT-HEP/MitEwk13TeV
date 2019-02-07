#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
# # directory of probes ntuples
NTUPLEDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_13TeV_Efficiency_v1/probes
OUTPUTDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_13TeV_Efficiency_v1/results

# NTUPLEDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2_Efficiency
# OUTPUTDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2_Efficiency-results

# integrated luminosity for data
LUMI13=213.1
LUMI5=291.1

#
# Muon efficiencies
#############################################
## 13 TeV ###################################
#############################################
# trigger efficiency
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff_Test/probes.root\",\"${OUTPUTDIR_13}/Zmm/MuHLTEff_Test/Positive\",\"png\",0,0,1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff_Test/probes.root\",\"${OUTPUTDIR_13}/Zmm/MuHLTEff_Test/Negative\",\"png\",0,0,-1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",2,2,1,1,\"${NTUPLEDIR_13}/Zmm/Data/MuHLTEff_Test/probes.root\",\"${OUTPUTDIR_13}/DataZmm/MuHLTEff_Test/Positive\",\"png\",0,0,1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",2,2,1,1,\"${NTUPLEDIR_13}/Zmm/Data/MuHLTEff_Test/probes.root\",\"${OUTPUTDIR_13}/DataZmm/MuHLTEff_Test/Negative\",\"png\",0,0,-1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff_Test/probes.root\"\)

# # id+iso efficiency
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuSelEff_Test/probes.root\",\"${OUTPUTDIR_13}/Zmm/MuSelEff_Test/Positive\",\"png\",0,0,1,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSelEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuSelEff_Test/probes.root\",\"${OUTPUTDIR_13}/Zmm/MuSelEff_Test/Negative\",\"png\",0,0,-1,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSelEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",2,2,1,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSelEff_Test/probes.root\",\"${OUTPUTDIR_13}/DataZmm/MuSelEff_Test/Positive\",\"png\",0,0,1,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSelEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",2,2,1,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSelEff_Test/probes.root\",\"${OUTPUTDIR_13}/DataZmm/MuSelEff_Test/Negative\",\"png\",0,0,-1,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSelEff_Test/probes.root\"\)

# # tracking efficiency
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuTrkEff_Test/probes.root\",\"${OUTPUTDIR_13}/Zmm/MuTrkEff_Test/Positive\",\"png\",0,0,1,\"Muon\",\"tracking\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuTrkEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuTrkEff_Test/probes.root\",\"${OUTPUTDIR_13}/Zmm/MuTrkEff_Test/Negative\",\"png\",0,0,-1,\"Muon\",\"tracking\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuTrkEff_Test/probes.root\"\)
root -l -b -q plotEff.C+\(\"mufineptbins.bins\",2,2,1,1,\"${NTUPLEDIR_13}/Zmm/Data/MuTrkEff_Test/probes.root\",\"${OUTPUTDIR_13}/DataZmm/MuTrkEff_Test/Positive\",\"png\",0,0,1,\"Muon\",\"tracking\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuTrkEff_Test/probes.root\"\)
root -l -b -q plotEff.C+\(\"mufineptbins.bins\",2,2,1,1,\"${NTUPLEDIR_13}/Zmm/Data/MuTrkEff_Test/probes.root\",\"${OUTPUTDIR_13}/DataZmm/MuTrkEff_Test/Negative\",\"png\",0,0,-1,\"Muon\",\"tracking\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuTrkEff_Test/probes.root\"\)

# # standalone efficiency
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_Test/probes.root\",\"${OUTPUTDIR_13}/Zmm/MuStaEff_Test/Positive\",\"png\",0,0,1,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_Test/probes.root\",\"${OUTPUTDIR_13}/Zmm/MuStaEff_Test/Negative\",\"png\",0,0,-1,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",2,2,1,1,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff_Test/probes.root\",\"${OUTPUTDIR_13}/DataZmm/MuStaEff_Test/Positive\",\"png\",0,0,1,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",2,2,1,1,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff_Test/probes.root\",\"${OUTPUTDIR_13}/DataZmm/MuStaEff_Test/Negative\",\"png\",0,0,-1,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_Test/probes.root\"\)

## 5 TeV ######################################
# trigger efficiency
# root -l -b -q plotEff.C+\(\"mupteta.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff_Test/probes.root\",\"${OUTPUTDIR_5}/Zmm/MuHLTEff_Test\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI},\"${NTUPLEDIR}/Zmm/MC/MuHLTEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff_Test/probes.root\",\"${OUTPUTDIR_5}/Zmm/MuHLTEff_Test\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",2,2,4,1,\"${NTUPLEDIR_5}/Zmm/Data/MuHLTEff_Test/probes.root\",\"${OUTPUTDIR_5}/DataZmmFit/MuHLTEff_Test\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/Data/MuHLTEff_Test/probes.root\",\"${OUTPUTDIR_5}/DataZmm/MuHLTEff_Test\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff_Test/probes.root\"\)

# id+iso efficiency
# root -l -b -q plotEff.C+\(\"mupteta.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff_Test/probes.root\",\"${OUTPUTDIR_5}/Zmm/MuSelEff_Test\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI},\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff_Test/probes.root\",\"${OUTPUTDIR_5}/Zmm/MuSelEff_Test\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",2,2,4,1,\"${NTUPLEDIR_5}/Zmm/Data/MuSelEff_Test/probes.root\",\"${OUTPUTDIR_5}/DataZmmFit/MuSelEff_Test\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/Data/MuSelEff_Test/probes.root\",\"${OUTPUTDIR_5}/DataZmm/MuSelEff_Test\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff_Test/probes.root\"\)

# tracking efficiency
# root -l -b -q plotEff.C+\(\"mupteta.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/MC/MuTrkEff_Test/probes.root\",\"${OUTPUTDIR_5}/Zmm/MuTrkEff_Test\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI},\"${NTUPLEDIR_5}/Zmm/MC/MuTrkEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuTrkEff_Test/probes.root\",\"${OUTPUTDIR_5}/Zmm/MuTrkEff_Test\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuTrkEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/Data/MuTrkEff_Test/probes.root\",\"${OUTPUTDIR_5}/DataZmmFit/MuTrkEff_Test\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuTrkEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/Data/MuTrkEff_Test/probes.root\",\"${OUTPUTDIR_5}/DataZmm/MuTrkEff_Test\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuTrkEff_Test/probes.root\"\)

# standalone efficiency
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_Test/probes.root\",\"${OUTPUTDIR_5}/Zmm/MuStaEff_Test\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_Test/probes.root\",\"${OUTPUTDIR_5}/Zmm/MuStaEff_Test\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"muSta_5.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/Data/MuStaEff_Test/probes.root\",\"${OUTPUTDIR_5}/DataZmmFit/MuStaEff_Test\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/Data/MuStaEff_Test/probes.root\",\"${OUTPUTDIR_5}/DataZmm/MuStaEff_Test\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_Test/probes.root\"\)


# standalone+id+iso efficiency
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm/MC/MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm/MuSelStaEff_iso\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"mupteta.bins\",2,1,2,1,\"${NTUPLEDIR}/Zmm/Data/MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/DataZmm/MuSelStaEff_iso\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm/MC/MuSelStaEff_iso/probes.root\"\)
# with a single bin
#root -l -b -q plotEff.C+\(\"musinglebin.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm/MC/MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm/MC/MuSelStaEff_iso_singleBin\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"musinglebin.bins\",2,1,2,1,\"${NTUPLEDIR}/Zmm/Data/MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm/Data/MuSelStaEff_iso_singleBin\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm/MC/MuSelStaEff_iso/probes.root\"\)

# standalone efficiency with deltaR matching
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm/MC/MuStaEff_iso_v2/probes.root\",\"${OUTPUTDIR}/Zmm/MC/MuStaEff_iso_v2\",\"png\",0,0,0,\"Muon\",\"stand-alone-deltaR\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm/Data/MuStaEff_iso_v2/probes.root\",\"${OUTPUTDIR}/Zmm/Data/MuStaEff_iso_v2\",\"png\",0,0,0,\"Muon\",\"stand-alone-deltaR\",0.7,1.02,${LUMI}\)

rm *.so *.d
