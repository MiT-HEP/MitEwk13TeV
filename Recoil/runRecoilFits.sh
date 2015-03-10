#! /bin/bash

INPUTDIR=/scratch/klawhorn/EWKAnaStore/8TeV/Selection

root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",2,2,1,\"ZmmData\"\)
#root -l -q fitRecoilZmm.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",2,2,1,\"ZmmMC\"\)
#root -l -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/data_select.root\",2,2,1,1,\"WmpMC\"\)
#root -l -q fitRecoilWm.C+\(\"${INPUTDIR}/Wmunu/ntuples/data_select.root\",2,2,1,-1,\"WmmMC\"\)

root -l -q fitRecoilZee.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",2,2,1,\"ZeeData\"\)
#root -l -q fitRecoilZee.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",2,2,1,\"ZeeMC\"\)
#root -l -q fitRecoilWe.C+\(\"${INPUTDIR}/Wenu/ntuples/data_select.root\",2,2,1,1,\"WepMC\"\)
#root -l -q fitRecoilWe.C+\(\"${INPUTDIR}/Wenu/ntuples/data_select.root\",2,2,1,-1,\"WemMC\"\)

rm *.so *.d
