#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# # directory of probes ntuples
# NTUPLEDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_13TeV_2017ID_Efficiency_v1/probes


POSTFIX=
POSTFIX2=
# NTUPLEDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_13TeV${POSTFIX2}/probes
NTUPLEDIR_13=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/probes
OUTPUTDIR_13=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results
# NTUPLEDIR_13=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV_biasedMC/probes
# OUTPUTDIR_13=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV_biasedMC/results
# OUTPUTDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_13TeV_2017ID_Efficiency_v1/results
# OUTPUTDIR_13=../LowPU2017ID_13TeV${POSTFIX2}/Plots
# /afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_5TeV

# NTUPLEDIR_5o=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_Efficiency_5TeV
# NTUPLEDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_Efficiency_matchTrig_5TeV
# OUTPUTDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_Efficiency_matchTrig_5TeV-results


NTUPLEDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2_Efficiency
OUTPUTDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2_Efficiency-results

# integrated luminosity for data
LUMI13=199.2
LUMI5=294.4

# POSTFIX=_aMCNLOxPythia
# POSTFIX=_POWHEGxPythia_v2
# POSTFIX=_POWHEGxPhotos
POSTFIX_CB=_CBxBW_v1


#################################################################
# Electron efficiencies
#
# 13 TeV shits
# # # supercluster efficiency
# root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleHLTEff_aMCxPythia${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleHLTEff_aMCxPythia${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/Data/EleHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleHLTEff_aMCxPythia${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/Data/EleHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleHLTEff_aMCxPythia${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\"\)

### Different Tag Pt cut
# root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleHLTEff_aMCxPythia${POSTFIX}_tagPt/Negative\",\"png\",0,0,-1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/HLTEff_tagPt/probes.root\"\)
# root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleHLTEff_aMCxPythia${POSTFIX}_tagPt/Positive\",\"png\",0,0,1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/HLTEff_tagPt/probes.root\"\)
# root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/Data/EleHLTEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleHLTEff_aMCxPythia${POSTFIX}_tagPt/Negative\",\"png\",0,0,-1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/HLTEff_tagPt/probes.root\"\)
# root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/Data/EleHLTEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleHLTEff_aMCxPythia${POSTFIX}_tagPt/Positive\",\"png\",0,0,1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/HLTEff_tagPt/probes.root\"\)

# # # # ####################### gsf+is+iso efficiency
# root -l -b -q plotEff.C+\(\"ele_gsf.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,1,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,7,2,7,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWBKG${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)

### Tag Pt
# root -l -b -q plotEff.C+\(\"ele_gsf.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\"\)
# root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,1,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\"\)


# root -l -b -q plotEff.C+\(\"ele_gsf.bins\",5,1,5,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"ele_gsf.bins\",6,1,6,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWxPhotos${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)

# root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,1,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_minloxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_minlo/probes.root\"\)


# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,7,2,7,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWBKG${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)

###################################################
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,1,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,1,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,7,2,7,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWBKG${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,7,2,7,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWBKG${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)

### Tag Pt
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt/Negative\",\"png\",0,0,-1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt/Positive\",\"png\",0,0,1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,1,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt/Negative\",\"png\",0,0,-1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,1,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt/Positive\",\"png\",0,0,1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\"\)


# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",5,1,5,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWxPythia${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",5,1,5,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWxPythia${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",6,1,6,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWxPhotos${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",6,1,6,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWxPhotos${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)

# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,1,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_minloxPythia${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_minlo/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,1,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_minloxPythia${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_minlo/probes.root\"\)

TYPE=_aMCxPythia
# # # # #### Draw Plots
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/pt/Negative\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"ele_hlt.bins\"\)
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/pt/Positive\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"ele_hlt.bins\"\)
 
 # root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/eta/Negative\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"ele_hlt.bins\"\)
 # root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/eta/Positive\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"ele_hlt.bins\"\)
  
 root -l -b -q plotDataMC.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/etapt/Negative\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"EtaPtBins\",0.00,1.20,$LUMI13,\"ele_hlt.bins\"\)
 root -l -b -q plotDataMC.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/etapt/Positive\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"EtaPtBins\",0.00,1.20,$LUMI13,\"ele_hlt.bins\"\)
 
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZeeGSFSel${TYPE}${POSTFIX}/pt/Combined\",\"$OUTPUTDIR_13/Zee/MC/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"ele_gsf.bins\"\)
 
 # root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZeeGSFSel${TYPE}${POSTFIX}/eta/Combined\",\"$OUTPUTDIR_13/Zee/MC/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"ele_gsf.bins\"\)
 
 root -l -b -q plotDataMC.C+\(\"$OUTPUTDIR_13/Plots/ZeeGSFSel${TYPE}${POSTFIX}/etapt/Combined\",\"$OUTPUTDIR_13/Zee/MC/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaPtBins\",0.00,1.20,$LUMI13,\"ele_gsf.bins\"\)
#####################################################################
# Electron efficiencies
###################### REPLACE THESE / COPY 13TeV DOWN HERE
#
# 5 TeV shits
# supercluster efficiency
# # # # root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,2,1,1,\"${NTUPLEDIR}/Zee_HLTEff/probes.root\",\"${OUTPUTDIR}/Zee/HLTEff\",\"png\",1,0,0,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI},\"${NTUPLEDIR}/Zee_HLTEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zee/MC/HLTEff/probes.root\",\"${OUTPUTDIR_5}/Zee/MC/HLTEff\",\"png\",0,0,0,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/HLTEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins_5.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zee/Data/HLTEff/probes.root\",\"${OUTPUTDIR_5}/Zee/Data/HLTEff\",\"png\",0,0,0,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/HLTEff/probes.root\"\)

# gsf+is+iso efficiency
# # # # root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,2,1,1,\"${NTUPLEDIR}/Zee_GSFSelEff/probes.root\",\"${OUTPUTDIR}/Zee/GSFSelEff\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_GSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zee/MC/GSFSelEff/probes.root\",\"${OUTPUTDIR_5}/Zee/MC/GSFSelEff\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/GSFSelEff/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins_5.bins\",2,2,1,1,\"${NTUPLEDIR_5}/Zee/Data/GSFSelEff/probes.root\",\"${OUTPUTDIR_5}/Zee/Data/GSFSelEff\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/GSFSelEff/probes.root\"\)


 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR/Plots/ZeeHLT\",\"$OUTPUTDIR/Zee/HLTEff/eff.root\",\"$OUTPUTDIR/DataZee/HLTEff/eff.root\",\"EtaBins\",0.00,1.20,$LUMI,\"elfineptbins.bins\"\)
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR/Plots/ZeeGSFSel\",\"$OUTPUTDIR/Zee/GSFSelEff/eff.root\",\"$OUTPUTDIR/Zee/Data/GSFSelEff/eff.root\",\"EtaBins\",0.00,1.20,$LUMI,\"elfineptbins.bins\"\)




# rm *.so *.d
