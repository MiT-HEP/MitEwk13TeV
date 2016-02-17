#!/bin/bash

lepton=$1
efftype=$2
effnum=$3
#bkgmodel=$3
xaxisname=$4
yaxisname=$5

workdir="/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

#gSystem->Load("${workdir}/Efficiency/selectProbes${lepton}Eff_C.so")
#selectProbes${lepton}Eff("${filedir}/${lepton}/ntuples/data_select.root", "${filedir}/${lepton}HLTEff/Data", 0, 0, 0)
#selectProbes${lepton}Eff("${filedir}/${lepton}/ntuples/zmm_select.raw.root", "${filedir}/${lepton}HLTEff/MC", 0, 1, 1)
#selectProbes${lepton}Eff("${filedir}/${lepton}/ntuples/data_select.root", "${filedir}/${lepton}StaEff/Data", 4, 0, 0)
#selectProbes${lepton}Eff("${filedir}/${lepton}/ntuples/zmm_select.raw.root", "${filedir}/${lepton}StaEff/MC", 4, 1, 1)
#selectProbes${lepton}Eff("${filedir}/${lepton}/ntuples/data_select.root", "${filedir}/${lepton}SITEff/Data", 8, 0, 0)
#selectProbes${lepton}Eff("${filedir}/${lepton}/ntuples/zmm_select.raw.root", "${filedir}/${lepton}SITEff/MC", 8, 1, 1)

root -l -b << EOF
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/CPlot_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/MitStyleRemix_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/CEffUser2D_cc.so")
.L ${workdir}/Efficiency/TagAndProbe/plotEff.C+g
gSystem->Load("${workdir}/Efficiency/TagAndProbe/plotEff_C.so")
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",2,6,2,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/1MG","png",0,0,0,"${lepton}","Sta",0.7,1.02,2263,"${filedir}/${lepton}StaEff/MC/probes.root")

.q
EOF

#gSystem->Load("${workdir}/Efficiency/TagAndProbe/plotDataMC_singlepTbins_C.so")
#plotDataMC_singlepTbins("${filedir}/Results","${filedir}/${lepton}${efftype}Eff/CT/eff.root", "${filedir}/${lepton}${efftype}Eff/MG/eff.root", "${lepton}${efftype}Eff", 0.01,1.2,2263, "${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins", "${xaxisname}", "${yaxisname}")

#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/Data/probes.root","${filedir}/${lepton}HLTEff/MG","png",0,0,0,"${lepton}","HLT",0.7,1.02,2263)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/MC/probes.root","${filedir}/${lepton}HLTEff/CT","png",0,0,0,"${lepton}","HLT",0.7,1.02,2263)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,1,2,1,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/MG","png",0,0,0,"${lepton}","SIT",0.7,1.02,2263,"${filedir}/${lepton}SITEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",0,0,0,0,"${filedir}/${lepton}SITEff/MC/probes.root","${filedir}/${lepton}SITEff/CT","png",0,0,0,"${lepton}","SIT",0.7,1.02,2263)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,6,2,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/MG","png",0,0,0,"${lepton}","Sta",0.7,1.02,2263,"${filedir}/${lepton}StaEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",0,0,0,0,"${filedir}/${lepton}StaEff/MC/probes.root","${filedir}/${lepton}StaEff/CT","png",0,0,0,"${lepton}","Sta",0.7,1.02,2263)

#sys
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",1,1,1,1,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/CB","png",0,0,0,"${lepton}","SIT",0.7,1.02,2263,"${filedir}/${lepton}SITEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,7,2,7,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/PL","png",0,0,0,"${lepton}","SIT",0.7,1.02,2263,"${filedir}/${lepton}SITEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",1,6,1,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/CB","png",0,0,0,"${lepton}","Sta",0.7,1.02,2263,"${filedir}/${lepton}StaEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}pteta.bins",2,7,2,7,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/PL","png",0,0,0,"${lepton}","Sta",0.7,1.02,2263,"${filedir}/${lepton}StaEff/MC/probes.root")

#bin
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/Data/probes.root","${filedir}/${lepton}HLTEff/1MG","png",0,0,0,"${lepton}","HLT",0.7,1.02,2263)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}HLTEff/MC/probes.root","${filedir}/${lepton}HLTEff/1CT","png",0,0,0,"${lepton}","HLT",0.7,1.02,2263)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",2,1,2,1,"${filedir}/${lepton}SITEff/Data/probes.root","${filedir}/${lepton}SITEff/1MG","png",0,0,0,"${lepton}","SIT",0.7,1.02,2263,"${filedir}/${lepton}SITEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}SITEff/MC/probes.root","${filedir}/${lepton}SITEff/1CT","png",0,0,0,"${lepton}","SIT",0.7,1.02,2263)
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",2,6,2,6,"${filedir}/${lepton}StaEff/Data/probes.root","${filedir}/${lepton}StaEff/1MG","png",0,0,0,"${lepton}","Sta",0.7,1.02,2263,"${filedir}/${lepton}StaEff/MC/probes.root")
#plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}1pteta.bins",0,0,0,0,"${filedir}/${lepton}StaEff/MC/probes.root","${filedir}/${lepton}StaEff/1CT","png",0,0,0,"${lepton}","Sta",0.7,1.02,2263)
