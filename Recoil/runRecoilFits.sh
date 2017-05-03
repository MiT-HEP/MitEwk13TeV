#! /bin/bash

# INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/SM/
#/afs/cern.ch/user/s/sabrandt/work/public/SM/newBacon/
INPUTDIR=/afs/cern.ch/work/a/arapyan/public/flat_ntuples

LUMI=2318
#LUMI=1 #fitting for the fractions
# muons
root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmDataPuppi_bkg\",${LUMI}\) & #Zmumu data/ /BKG selection selection
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmDataPuppi\",${LUMI}\) #Zmumu data selection
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmMCPuppi\",${LUMI}\) #Zmumu MC selection

# CENTRAL VALUE
#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,-1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmmMCPuppi\",${LUMI},0\) #Wmunu data selection
#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmpMCPuppi\",${LUMI},0\) & #Wmunu data selection

# W/Z recoil differences
###                 int etaBinCategory=0 // 0 is inclusive, 1 is fabs(eta)<=0.5,  2 is fabs(eta)=[0.5,1], 3 is fabs(eta)>=1
#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,-1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmmMCPuppi_rap05\",${LUMI},1\) &  #Wmunu
#sleep 100
#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmpMCPuppi_rap05\",${LUMI},1\) & #Wmunu
#sleep 100
#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,-1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmmMCPuppi_rap05-1\",${LUMI},2\) & #Wmunu
#sleep 100
#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmpMCPuppi_rap05-1\",${LUMI},2\) & #Wmunu
#sleep 100
#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,-1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmmMCPuppi_rap1\",${LUMI},3\) & #Wmunu
#sleep 100
#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmpMCPuppi_rap1\",${LUMI},3\) & #Wmunu

# Z recoil differences
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmMCPuppi_rap05\",${LUMI},1\) & #Zmumu
#sleep 100
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmMCPuppi_rap05-1\",${LUMI},2\) & #Zmumu
#sleep 100
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmMCPuppi_rap1\",${LUMI},3\) & #Zmumu
#sleep 100
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmDataPuppi_bkg_rap05\",${LUMI},1\) & #Zmumu data/ /BKG selection selection
#sleep 100
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmDataPuppi_bkg_rap05-1\",${LUMI},2\) & #Zmumu data/ /BKG selection selection
#sleep 100
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmDataPuppi_bkg_rap1\",${LUMI},3\) & #Zmumu data/ /BKG selection selection

# PileUp systematics
#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,-1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmmMCPuppi_PileupUp\",${LUMI},0\) & #Wmunu data selection
#sleep 100
#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmpMCPuppi_PileupUp\",${LUMI},0\) & #Wmunu data selection
#sleep 100
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmMCPuppi_PileupUp\",${LUMI},0\) & #Zmumu MC selection
#sleep 100
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmDataPuppi_bkg_PileupUp\",${LUMI},0\) & #Zmumu data/ /BKG selection selection
#sleep 100

#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,-1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmmMCPuppi_PileupDown\",${LUMI},0\) & #Wmunu data selection
#sleep 100
#root -l -b -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmpMCPuppi_PileupDown\",${LUMI},0\) & #Wmunu data selection
#sleep 100
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmMCPuppi_PileupDown\",${LUMI},0\) & #Zmumu MC selection
#sleep 100
#root -l -b -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmDataPuppi_bkg_PileupDown\",${LUMI},0\) & #Zmumu data/ /BKG selection selection
#sleep 100


## below old configs from stephanie
#root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,0,\"u1\",\"u2\",\"pf\",\"ZmmDataPF_bkg\",${LUMI}\) #Zmumu data/ /BKG selection selection  
#root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",3,3,1,\"u1\",\"u2\",\"pf\",\"ZmmDataPF\",${LUMI}\) #Zmumu data selection
#root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"u1\",\"u2\",\"pf\",\"ZmmMCPF\",${LUMI}\) #Zmumu MC selection

#   root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",3,3,1,\"puppiU1\",\"puppiU2\",\"puppi\",\"ZmmMCPuppi_newBacon_eta1\"\) #Zmumu data selection
# root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",2,2,1,\"ppU1\",\"ppU2\",\"t1pf\",\"ZmmData\"\) #Zmumu data selection
# root -l -q fitRecoilWm3.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmpMC_newBacon\"\) #Wmunu data selection
# root -l -q fitRecoilWm3.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",3,3,1,-1,0,\"puppiU1\",\"puppiU2\",\"puppi\",\"WmmMC_newBacon\"\) #Wmunu data selection
# root -l -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",2,2,1,-1,0,\"mvaU1\",\"mvaU2\",\"mva\",\"WmmMC\"\) #Wmunu signal MC -

# electrons
#   root -l -q fitRecoilZee.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.raw.root\",2,2,1,\"mvaU1\",\"mvaU2\",\"mva\",\"ZeeData\"\) #Zmumu data selection
# root -l -q fitRecoilZee.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",2,2,1,\"ZeeMC\"\) #Zee signal MC selection
# root -l -q fitRecoilWe.C+\(\"${INPUTDIR}/Wenu/ntuples/\",2,2,1,1,0,\"mvaU1\",\"mvaU2\",\"mva\"\"WepMC\"\) #Wenu +
# root -l -q fitRecoilWe.C+\(\"${INPUTDIR}/Wenu/ntuples/\",2,2,1,-1,0,\"mvaU1\",\"mvaU2\",\"mva\"\"WemMC\"\) #Wenu signal MC -

rm *.so *.d *.pcm
