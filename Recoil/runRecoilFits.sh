#! /bin/bash

INPUTDIR=/data/blue/Bacon/Run2/wz_flat

# muons
# root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",2,2,1,\"u1\",\"u2\",\"pf\",\"ZmmData\"\) #Zmumu data selection
root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",2,2,1,\"mvaU1\",\"u2\",\"mva\",\"ZmmData\"\) #Zmumu data selection
root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",2,2,1,\"ppU1\",\"u2\",\"t1pf\",\"ZmmData\"\) #Zmumu data selection
root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",2,2,1,\"tkU1\",\"u2\",\"tk\",\"ZmmData\"\) #Zmumu data selection
#root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",2,2,1,\"ZmmMC\"\) #Zmumu signal MC selection
#root -l -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",2,2,1,1,0,\"WmpMC\"\) #Wmunu data selection
#root -l -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/\",2,2,1,-1,0,\"WmmMC\"\) #Wmunu signal MC -

# electrons
# root -l -q fitRecoilZee.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",2,2,1,\"ZeeData\"\) #Zee data selection
 #root -l -q fitRecoilZee.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",2,2,1,\"ZeeMC\"\) #Zee signal MC selection
# root -l -q fitRecoilWe.C+\(\"${INPUTDIR}/Wenu/ntuples/\",2,2,1,1,0,\"WepMC\"\) #Wenu +
# root -l -q fitRecoilWe.C+\(\"${INPUTDIR}/Wenu/ntuples/\",2,2,1,-1,0,\"WemMC\"\) #Wenu signal MC -

rm *.so *.d *.pcm
