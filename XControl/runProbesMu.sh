#!/bin/bash

lepton=$1
efftype=$2
xaxisname=$3
yaxisname=$4

workdir="/afs/cern.ch/work/x/xniu/public/Lumi/Ele/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"
aramfiledir="/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Zmumu"

#select probles
#root -l -b << EOF
#gSystem->Load("${workdir}/Efficiency/selectProbes${lepton}Eff_C.so")
#selectProbes${lepton}Eff("${aramfiledir}/ntuples/data_select.root", "${filedir}/${lepton}HLTEff/Data", 0, 0, 0)
#selectProbes${lepton}Eff("${aramfiledir}/ntuples/zmm_select.raw.root", "${filedir}/${lepton}HLTEff/MC", 0, 1, 1)
#selectProbes${lepton}Eff("${aramfiledir}/ntuples/data_select.root", "${filedir}/${lepton}StaEff/Data", 4, 0, 0)
#selectProbes${lepton}Eff("${aramfiledir}/ntuples/zmm_select.raw.root", "${filedir}/${lepton}StaEff/MC", 4, 1, 1)
#selectProbes${lepton}Eff("${aramfiledir}/ntuples/data_select.root", "${filedir}/${lepton}SITEff/Data", 8, 0, 0)
#selectProbes${lepton}Eff("${aramfiledir}/ntuples/zmm_select.raw.root", "${filedir}/${lepton}SITEff/MC", 8, 1, 1)

#do eff nominal
#root -l -b << EOF
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CPlot_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/MitStyleRemix_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser1D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser2D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/plotEff_C.so")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/Data/probes.root","${filedir}/${lepton}HLTEff/MGpositive","png",0,0,1,"${lepton}","HLT",0.7,1.02,2305)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/MC/probes.root","${filedir}/${lepton}HLTEff/CTpositive","png",0,0,1,"${lepton}","HLT",0.7,1.02,2305)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",2,1,2,1,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/MGpositive_FineBin","png",0,0,1,"${lepton}","SIT",0.7,1.02,2305,"${filedir}/${lepton}SITEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",0,0,0,0,"${filedir}/${lepton}SITEff/MC/probes.root","${filedir}/${lepton}SITEff/CTpositive","png",0,0,1,"${lepton}","SIT",0.7,1.02,2305)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,6,2,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/MGpositive","png",0,0,1,"${lepton}","Sta",0.7,1.02,2305,"${filedir}/${lepton}StaEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",0,0,0,0,"${filedir}/${lepton}StaEff/MC/probes.root","${filedir}/${lepton}StaEff/CTpositive","png",0,0,1,"${lepton}","Sta",0.7,1.02,2305)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/Data/probes.root","${filedir}/${lepton}HLTEff/MGnegative","png",0,0,-1,"${lepton}","HLT",0.7,1.02,2305)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/MC/probes.root","${filedir}/${lepton}HLTEff/CTnegative","png",0,0,-1,"${lepton}","HLT",0.7,1.02,2305)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",2,1,2,1,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/MGnegative_FineBin","png",0,0,-1,"${lepton}","SIT",0.7,1.02,2305,"${filedir}/${lepton}SITEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",0,0,0,0,"${filedir}/${lepton}SITEff/MC/probes.root","${filedir}/${lepton}SITEff/CTnegative","png",0,0,-1,"${lepton}","SIT",0.7,1.02,2305)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,6,2,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/MGnegative","png",0,0,-1,"${lepton}","Sta",0.7,1.02,2305,"${filedir}/${lepton}StaEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",0,0,0,0,"${filedir}/${lepton}StaEff/MC/probes.root","${filedir}/${lepton}StaEff/CTnegative","png",0,0,-1,"${lepton}","Sta",0.7,1.02,2305)
#.q
#EOF


#do eff shape systematics:SIT
#root -l -b << EOF
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CPlot_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/MitStyleRemix_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser1D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser2D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/plotChargeDependentEff_C.so")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,1,2,1,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/MGpositive","png",0,0,1,"${lepton}","SIT",0.7,1.02,2305,"${filedir}/${lepton}SITEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,1,2,1,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/MGnegative","png",0,0,-1,"${lepton}","SIT",0.7,1.02,2305,"${filedir}/${lepton}SITEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",1,1,1,1,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/CBpositive","png",0,0,1,"${lepton}","SIT",0.7,1.02,2305,"${filedir}/${lepton}SITEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",1,1,1,1,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/CBnegative","png",0,0,-1,"${lepton}","SIT",0.7,1.02,2305,"${filedir}/${lepton}SITEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,7,2,7,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/PLpositive","png",0,0,1,"${lepton}","SIT",0.7,1.02,2305,"${filedir}/${lepton}SITEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,7,2,7,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/PLnegative","png",0,0,-1,"${lepton}","SIT",0.7,1.02,2305,"${filedir}/${lepton}SITEff/MC/probes.root")
#.q
#EOF


#do eff shape systematics:Sta
#root -l -b << EOF
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CPlot_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/MitStyleRemix_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser1D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser2D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/plotChargeDependentEff_C.so")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,6,2,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/MGpositive","png",0,0,1,"${lepton}","Sta",0.7,1.02,2305,"${filedir}/${lepton}StaEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,6,2,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/MGnegative","png",0,0,-1,"${lepton}","Sta",0.7,1.02,2305,"${filedir}/${lepton}StaEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",1,6,1,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/CBpositive","png",0,0,1,"${lepton}","Sta",0.7,1.02,2305,"${filedir}/${lepton}StaEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",1,6,1,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/CBnegative","png",0,0,-1,"${lepton}","Sta",0.7,1.02,2305,"${filedir}/${lepton}StaEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,7,2,7,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/PLpositive","png",0,0,1,"${lepton}","Sta",0.7,1.02,2305,"${filedir}/${lepton}StaEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,7,2,7,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/PLnegative","png",0,0,-1,"${lepton}","Sta",0.7,1.02,2305,"${filedir}/${lepton}StaEff/MC/probes.root")
#.q
#EOF

#do eff binninb systematics
root -l -b << EOF
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/CPlot_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/MitStyleRemix_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser2D_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/plotChargeDependentEff_C.so")
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/Data/probes.root","${filedir}/${lepton}HLTEff/1MGpositive","png",1,0,1,"${lepton}","HLT",0.7,1.02,2305)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/MC/probes.root","${filedir}/${lepton}HLTEff/1CTpositive","png",1,0,1,"${lepton}","HLT",0.7,1.02,2305)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",2,1,2,1,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/1MGpositive","png",1,0,1,"${lepton}","SIT",0.7,1.02,2305,"${filedir}/${lepton}SITEff/MC/probes.root")
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}SITEff/MC/probes.root","${filedir}/${lepton}SITEff/1CTpositive","png",1,0,1,"${lepton}","SIT",0.7,1.02,2305)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",2,6,2,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/1MGpositive","png",1,0,1,"${lepton}","Sta",0.7,1.02,2305,"${filedir}/${lepton}StaEff/MC/probes.root")
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}StaEff/MC/probes.root","${filedir}/${lepton}StaEff/1CTpositive","png",1,0,1,"${lepton}","Sta",0.7,1.02,2305)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/Data/probes.root","${filedir}/${lepton}HLTEff/1MGnegative","png",1,0,-1,"${lepton}","HLT",0.7,1.02,2305)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/MC/probes.root","${filedir}/${lepton}HLTEff/1CTnegative","png",1,0,-1,"${lepton}","HLT",0.7,1.02,2305)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",2,1,2,1,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/1MGnegative","png",1,0,-1,"${lepton}","SIT",0.7,1.02,2305,"${filedir}/${lepton}SITEff/MC/probes.root")
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}SITEff/MC/probes.root","${filedir}/${lepton}SITEff/1CTnegative","png",1,0,-1,"${lepton}","SIT",0.7,1.02,2305)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",2,6,2,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/1MGnegative","png",1,0,-1,"${lepton}","Sta",0.7,1.02,2305,"${filedir}/${lepton}StaEff/MC/probes.root")
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}StaEff/MC/probes.root","${filedir}/${lepton}StaEff/1CTnegative","png",1,0,-1,"${lepton}","Sta",0.7,1.02,2305)
.q
EOF

