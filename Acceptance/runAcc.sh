#! /bin/bash

OUTPUTDIR5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2_Acceptance
INPUTDIR5=/eos/cms/store/user/sabrandt/StandardModel/LowPU_5TeV_Try2/
# OUTPUTDIR13=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_Acceptance_Test_13TeV
# INPUTDIR13=/eos/cms/store/user/sabrandt/StandardModel/LowPU_13TeV/
#
# W->munu
#
#root -l -q computeAccGenWm.C+\(\"wmp.conf\",\"Wmunu/plus\",1\)
#root -l -q computeAccGenWm.C+\(\"wmm.conf\",\"Wmunu/minus\",-1\)
#root -l -q computeAccSelWm.C+\(\"wmp.conf\",\"Wmunu/plus\",1\)
#root -l -q computeAccSelWm.C+\(\"wmm.conf\",\"Wmunu/minus\",-1\)

#
# W->enu
#
#root -l -q computeAccGenWe.C+\(\"wep.conf\",\"Wenu/plus\",1\)
#root -l -q computeAccGenWe.C+\(\"wem.conf\",\"Wenu/minus\",-1\)
#root -l -q computeAccSCWe.C+\(\"wep.conf\",\"Wenu/plus\",1\)
#root -l -q computeAccSCWe.C+\(\"wem.conf\",\"Wenu/minus\",-1\)
#root -l -q computeAccSelWe.C+\(\"wep.conf\",\"Wenu/plus\",1\)
#root -l -q computeAccSelWe.C+\(\"wem.conf\",\"Wenu/minus\",-1\)

#
# Z->mumu
#
#root -l -q computeAccGenZmm.C+\(\"zmm.conf\",\"Zmumu\"\)
#root -l -q computeAccSelZmm.C+\(\"zmm.conf\",\"Zmumu\"\)
#root -l -q computeAccSelZmmBinned.C+\(\"zmm.conf\",\"Zmumu\"\)

#
# Z->ee
#
#root -l -q computeAccGenZee.C+\(\"zee.conf\",\"Zee\"\)
#root -l -q computeAccSCZee.C+\(\"zee.conf\",\"Zee\"\)
# root -l -q computeAccSelZee.C+\(\"zee.conf\",\"${OUTPUTDIR}/Zee\"\)
root -l -q computeAccSelZeeBinned_Charge.C+\(\"zee.conf\",\"blah\",\"${INPUTDIR5}\",\"${OUTPUTDIR5}/Zee\",25,2.5,0,1,0,0,0\)
# root -l -q computeAccSelZeeBinned_Charge.C+\(\"zee_13.conf\",\"blah\",\"${INPUTDIR13}\",\"${OUTPUTDIR13}/Zee\",25,2.5,0,1,0,0,1\)

rm *.so *.d
