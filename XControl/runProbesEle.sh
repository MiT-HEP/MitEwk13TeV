#!/bin/bash

lepton=$1
efftype=$2
xaxisname=$3
yaxisname=$4

workdir="/afs/cern.ch/work/x/xniu/public/Lumi/Ele/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"
aramfiledir="/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Zee"

#select probles
#root -l -b << EOF
#gSystem->Load("${workdir}/Efficiency/selectProbes${lepton}Eff_C.so")
#selectProbes${lepton}Eff("${aramfiledir}/ntuples/data_select.root", "${filedir}/${lepton}HLTEff/Data", 0, 0, 0)
#selectProbes${lepton}Eff("${aramfiledir}/ntuples/zee_select.root", "${filedir}/${lepton}HLTEff/MC", 0, 1, 1)
#selectProbes${lepton}Eff("${aramfiledir}/ntuples/data_select.root", "${filedir}/${lepton}GsfSelEff/Data", 4, 0, 0)
#selectProbes${lepton}Eff("${aramfiledir}/ntuples/zee_select.root", "${filedir}/${lepton}GsfSelEff/MC", 4, 1, 1)
#.q
#EOF

#do eff. nominal
#for parametrizating trigger eff.
#root -l -b << EOF
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CPlot_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/MitStyleRemix_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser1D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser2D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/printEff_C.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/plotEff_C.so")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}HLTpteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/Data/probes.root","${filedir}/${lepton}HLTEff/MGpositive","png",0,0,1,"${lepton}","HLT",0.7,1.02,2305)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}HLTpteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/MC/probes.root","${filedir}/${lepton}HLTEff/CTpositive","png",0,0,1,"${lepton}","HLT",0.7,1.02,2305)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}HLTpteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/Data/probes.root","${filedir}/${lepton}HLTEff/MGnegative","png",0,0,-1,"${lepton}","HLT",0.7,1.02,2305)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}HLTpteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/MC/probes.root","${filedir}/${lepton}HLTEff/CTnegative","png",0,0,-1,"${lepton}","HLT",0.7,1.02,2305)
#.q
#EOF

#for GsfSel
#root -l -b << EOF
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CPlot_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/MitStyleRemix_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser1D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser2D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/printEff_C.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/plotChargeDependentEff_C.so")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",2,1,2,2,"${filedir}/${lepton}GsfSelEff/Data/probes.root","${filedir}/${lepton}GsfSelEff/MGpositive","png",0,0,1,"${lepton}","GsfSel",0.7,1.02,2305,"${filedir}/${lepton}GsfSelEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",0,0,0,0,"${filedir}/${lepton}GsfSelEff/MC/probes.root","${filedir}/${lepton}GsfSelEff/CTpositive","png",0,0,1,"${lepton}","GsfSel",0.7,1.02,2305)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",2,1,2,2,"${filedir}/${lepton}GsfSelEff/Data/probes.root","${filedir}/${lepton}GsfSelEff/MGnegative","png",0,0,-1,"${lepton}","GsfSel",0.7,1.02,2305,"${filedir}/${lepton}GsfSelEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}finepteta.bins",0,0,0,0,"${filedir}/${lepton}GsfSelEff/MC/probes.root","${filedir}/${lepton}GsfSelEff/CTnegative","png",0,0,-1,"${lepton}","GsfSel",0.7,1.02,2305)
#.q
#EOF

#for shape syst.:
#root -l -b << EOF
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CPlot_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/MitStyleRemix_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser1D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser2D_cc.so")
#gSystem->Load("${workdir}/Efficiency/TagAndProbe/plotChargeDependentEff_C.so")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,1,2,2,"${filedir}/${lepton}GsfSelEff/Data/probes.root","${filedir}/${lepton}GsfSelEff/MGpositive_sys","png",0,0,1,"${lepton}","GsfSel",0.7,1.02,2305,"${filedir}/${lepton}GsfSelEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,1,2,2,"${filedir}/${lepton}GsfSelEff/Data/probes.root","${filedir}/${lepton}GsfSelEff/MGnegative_sys","png",0,0,-1,"${lepton}","GsfSel",0.7,1.02,2305,"${filedir}/${lepton}GsfSelEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",1,1,1,2,"${filedir}/${lepton}GsfSelEff/Data/probes.root","${filedir}/${lepton}GsfSelEff/CBpositive","png",0,0,1,"${lepton}","GsfSel",0.7,1.02,2305,"${filedir}/${lepton}GsfSelEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",1,1,1,2,"${filedir}/${lepton}GsfSelEff/Data/probes.root","${filedir}/${lepton}GsfSelEff/CBnegative","png",0,0,-1,"${lepton}","GsfSel",0.7,1.02,2305,"${filedir}/${lepton}GsfSelEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,7,2,7,"${filedir}/${lepton}GsfSelEff/Data/probes.root","${filedir}/${lepton}GsfSelEff/PLpositive","png",0,0,1,"${lepton}","GsfSel",0.7,1.02,2305,"${filedir}/${lepton}GsfSelEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,7,2,7,"${filedir}/${lepton}GsfSelEff/Data/probes.root","${filedir}/${lepton}GsfSelEff/PLnegative","png",0,0,-1,"${lepton}","GsfSel",0.7,1.02,2305,"${filedir}/${lepton}GsfSelEff/MC/probes.root")
#.q
#EOF


#for binning syst.:
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
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/Data/probes.root","${filedir}/${lepton}HLTEff/1MGnegative","png",1,0,-1,"${lepton}","HLT",0.7,1.02,2305)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/MC/probes.root","${filedir}/${lepton}HLTEff/1CTnegative","png",1,0,-1,"${lepton}","HLT",0.7,1.02,2305)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",2,1,2,2,"${filedir}/${lepton}GsfSelEff/Data/probes.root","${filedir}/${lepton}GsfSelEff/1MGpositive","png",1,0,1,"${lepton}","GsfSel",0.7,1.02,2305,"${filedir}/${lepton}GsfSelEff/MC/probes.root")
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}GsfSelEff/MC/probes.root","${filedir}/${lepton}GsfSelEff/1CTpositive","png",1,0,1,"${lepton}","GsfSel",0.7,1.02,2305)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",2,1,2,2,"${filedir}/${lepton}GsfSelEff/Data/probes.root","${filedir}/${lepton}GsfSelEff/1MGnegative","png",1,0,-1,"${lepton}","GsfSel",0.7,1.02,2305,"${filedir}/${lepton}GsfSelEff/MC/probes.root")
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}GsfSelEff/MC/probes.root","${filedir}/${lepton}GsfSelEff/1CTnegative","png",1,0,-1,"${lepton}","GsfSel",0.7,1.02,2305)
.q
EOF
