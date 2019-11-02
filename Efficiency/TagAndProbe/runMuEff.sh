#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
# # directory of probes ntuples
# NTUPLEDIR_13=../testReweights_v2_2
# OUTPUTDIR_13=../testReweights_v2_2/results


POSTFIX=
v=
NTUPLEDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency${v}/LowPU2017ID_13TeV/probes 
OUTPUTDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency${v}/LowPU2017ID_13TeV/results

# NTUPLEDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2_Efficiency
# OUTPUTDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2_Efficiency-results

# integrated luminosity for data
LUMI13=199.270
LUMI5=291.1


# POSTFIX_CB=_CBxBW_v1
#
# Muon efficiencies
#############################################
## 13 TeV ###################################
#############################################
# trigger efficiency
# root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuHLTEff_aMCxPythia/Positive\",\"png\",0,0,1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuHLTEff_aMCxPythia/Negative\",\"png\",0,0,-1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\"\)

# root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/Data/MuHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuHLTEff_aMCxPythia/Positive\",\"png\",0,0,1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/Data/MuHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuHLTEff_aMCxPythia/Negative\",\"png\",0,0,-1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\"\)


# root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuSITEff_aMCxPythia/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)

# root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuSITEff_aMCxPythia_tagPt/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPt/probes.root\"\)

# root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)

#mc
# root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuStaEff_aMCxPythia/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)

# root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuStaEff_aMCxPythia_tagPt/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPt/probes.root\"\)
# # # amc@nlo x pythia
# root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,6,2,6,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,6,2,6,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)

# # id+iso efficiency
################  mc  #####################
# root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuSITEff_aMCxPythia/Positive\",\"png\",0,0,1,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuSITEff_aMCxPythia/Negative\",\"png\",0,0,-1,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)
############ DATA #######################
## Tag pt change cut
# root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_aMCxPythia_tagPt${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPt/probes.root\"\)
# # # amc@nlo x pythia
# root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)
# # # powheg x pythia
# root -l -b -q plotEff.C+\(\"mu_sit.bins\",5,1,5,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_POWxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)
# # # powheg x photos
# root -l -b -q plotEff.C+\(\"mu_sit.bins\",6,1,6,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_POWxPhotos${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)
# # # minlo 
# root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_minloxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_minlo/probes.root\"\)
# # # # Alternate BKG fit
# root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,7,2,7,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_POWBKG${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)


# # # standalone efficiency
#mc
# root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuStaEff_aMCxPythia/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)
# # # amc@nlo x pythia
root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,6,2,6,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_aMCxPythia_tagPt${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPt/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,6,2,6,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)
# # # # powheg x pythia
# root -l -b -q plotEff.C+\(\"mu_sta.bins\",5,6,5,6,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_POWxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)
# # # # # # powheg x photos
# root -l -b -q plotEff.C+\(\"mu_sta.bins\",6,6,6,6,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_POWxPhotos${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)
# # # # # minlo
# root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,6,2,6,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_minloxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_minlo/probes.root\"\)
# # # bkg power law
# root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,7,2,7,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_POWBKG${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)


# # ## Make a bunch of plots ##
TYPE=_aMCxPythia_tagPt
# # # ## TYPE
 # # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZmmHLT${TYPE}NewTAble/pt/Negative\",\"$OUTPUTDIR_13/Zmm/MC/MuHLTEff${TYPE}/Negative/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuHLTEff${TYPE}/Negative/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"mufineptbins.bins\"\)
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZmmHLT${TYPE}${POSTFIX}/pt/Negative\",\"$OUTPUTDIR_13/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"mu_hlt.bins\"\)
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZmmHLT${TYPE}${POSTFIX}/pt/Positive\",\"$OUTPUTDIR_13/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"mu_hlt.bins\"\)
 
 # root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZmmHLT${TYPE}${POSTFIX}/eta/Negative\",\"$OUTPUTDIR_13/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"mu_hlt.bins\"\)
 # root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZmmHLT${TYPE}${POSTFIX}/eta/Positive\",\"$OUTPUTDIR_13/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"mu_hlt.bins\"\)
  
 # # ## SEL
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZmmSIT${TYPE}${POSTFIX}/pt/Combined\",\"$OUTPUTDIR_13/Zmm/MC/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"mu_sit.bins\"\)
 
 # root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZmmSIT${TYPE}${POSTFIX}/eta/Combined\",\"$OUTPUTDIR_13/Zmm/MC/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"mu_sit.bins\"\)
 
 # # # # # # STA
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZmmSta${TYPE}${POSTFIX}/pt/Combined\",\"$OUTPUTDIR_13/Zmm/MC/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"mu_sta.bins\"\)
 
 # root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZmmSta${TYPE}${POSTFIX}/eta/Combined\",\"$OUTPUTDIR_13/Zmm/MC/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"mu_sta.bins\"\)

## 5 TeV ######################################
############ NEED TO BE REDONE #######################
# trigger efficiency
# root -l -b -q plotEff.C+\(\"mupteta.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuHLTEff\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI},\"${NTUPLEDIR}/Zmm/MC/MuHLTEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuHLTEff\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",2,2,4,1,\"${NTUPLEDIR_5}/Zmm/Data/MuHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/DataFit/MuHLTEff\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/Data/MuHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuHLTEff\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff/probes.root\"\)

# id+iso efficiency
# root -l -b -q plotEff.C+\(\"mupteta.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuSelEff\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI},\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuSelEff\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",2,2,4,1,\"${NTUPLEDIR_5}/Zmm/Data/MuSelEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/DataFit/MuSelEff\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/Data/MuSelEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuSelEff\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSelEff/probes.root\"\)

# tracking efficiency
# root -l -b -q plotEff.C+\(\"mupteta.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuSITEff\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuSITEff\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/DataFit/MuSITEff\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuSITEff\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\"\)

# standalone efficiency
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuStaEff\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuStaEff\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"muSta_5.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/DataFit/MuStaEff\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"mufineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuStaEff\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\"\)

#############################################
# standalone+id+iso efficiency
# root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm/MC/MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm/MC/MuSelStaEff_iso\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI}\)
# root -l -b -q plotEff.C+\(\"mupteta.bins\",2,1,2,1,\"${NTUPLEDIR}/Zmm/Data/MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm/Data/MuSelStaEff_iso\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm/MC/MuSelStaEff_iso/probes.root\"\)
# with a single bin
# root -l -b -q plotEff.C+\(\"musinglebin.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm/MC/MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm/MC/MuSelStaEff_iso_singleBin\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI}\)
# root -l -b -q plotEff.C+\(\"musinglebin.bins\",2,1,2,1,\"${NTUPLEDIR}/Zmm/Data/MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm/Data/MuSelStaEff_iso_singleBin\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm/MC/MuSelStaEff_iso/probes.root\"\)

# standalone efficiency with deltaR matching
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm/MC/MuStaEff_iso_v2/probes.root\",\"${OUTPUTDIR}/Zmm/MC/MuStaEff_iso_v2\",\"png\",0,0,0,\"Muon\",\"stand-alone-deltaR\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm/Data/MuStaEff_iso_v2/probes.root\",\"${OUTPUTDIR}/Zmm/Data/MuStaEff_iso_v2\",\"png\",0,0,0,\"Muon\",\"stand-alone-deltaR\",0.7,1.02,${LUMI}\)

# rm *.so *.d
