#! /bin/bash

# INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/SM/
#/afs/cern.ch/user/s/sabrandt/work/public/SM/newBacon/
# INPUTDIR=/data/t3home000/sabrandt/2018_09_07_Masters_Incl_2
# INPUTDIR=/data/t3home000/sabrandt/2018_09_07_Masters_Incl
INPUTDIR=/data/t3home000/sabrandt/2018_11_11_LowPUFlat_Aram

LUMI=2152
#LUMI=1 #fitting for the fractions

# root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/we_select.root\",3,3,1,1,0,\"puppiU1\",\"puppiU2\",\"puppiMet\",\"WepMCPuppi_test\",${LUMI},0\)  #Zmumu data/ /BKG selection selection
# root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppiMet\",\"ZmmMCPuppi_genBin\",${LUMI},0\)  #Zmumu data/ /BKG selection selection
 root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"u1\",\"u2\",\"met\",\"ZmmMCPF_lowPU13_15GeV\",${LUMI},0\)  #Zmumu data/ /BKG selection selection
# root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"u1\",\"u2\",\"met\",\"ZmmDataPF_bkgTopEWK_lowPU13_15GeV\",${LUMI},0\)  #Zmumu data/ /BKG selection selection

# muons
# root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmDataPuppi_bkgTopEWK_TEST\",${LUMI},0\)  #Zmumu data/ /BKG selection selection
# root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmMCPuppi_puDown\",${LUMI},0\)  #Zmumu data/ /BKG selection selection
# root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmDataPuppi_bkgTopEWK_rap05\",${LUMI},1\)  #Zmumu data/ /BKG selection selection
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmDataPuppi_bkgTopEWK_rap05-1\",${LUMI},2\)  #Zmumu data/ /BKG selection selection
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmDataPuppi_bkgTopEWK_rap1\",${LUMI},3\)  #Zmumu data/ /BKG selection selection
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmMCPuppi_rap05\",${LUMI},1\)  #Zmumu signal MC selection
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmMCPuppi_rap05-1\",${LUMI},2\)  #Zmumu signal MC selection
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmMCPuppi_rap1\",${LUMI},3\)  #Zmumu signal MC selection
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZeeMCPuppi\",${LUMI}\) & #Zmumu signal MC selection

#rm *.so *.d *.pcm
