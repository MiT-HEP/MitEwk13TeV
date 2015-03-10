#! /bin/bash

#
# W->munu
#
#root -l -q computeAccGenWm.C+\(\"wmp.conf\",\"Wmunu/plus\",1\)
#root -l -q computeAccGenWm.C+\(\"wmm.conf\",\"Wmunu/minus\",-1\)
root -l -q computeAccSelWm.C+\(\"wmp.conf\",\"Wmunu/plus\",1\)
root -l -q computeAccSelWm.C+\(\"wmm.conf\",\"Wmunu/minus\",-1\)

#
# W->enu
#
#root -l -q computeAccGenWe.C+\(\"wep.conf\",\"Wenu/plus\",1\)
#root -l -q computeAccGenWe.C+\(\"wem.conf\",\"Wenu/minus\",-1\)
#root -l -q computeAccSCWe.C+\(\"wep.conf\",\"Wenu/plus\",1\)
#root -l -q computeAccSCWe.C+\(\"wem.conf\",\"Wenu/minus\",-1\)
root -l -q computeAccSelWe.C+\(\"wep.conf\",\"Wenu/plus\",1\)
root -l -q computeAccSelWe.C+\(\"wem.conf\",\"Wenu/minus\",-1\)

#
# Z->mumu
#
#root -l -q computeAccGenZmm.C+\(\"zmm.conf\",\"Zmumu\"\)
#root -l -q computeAccSelZmm.C+\(\"zmm.conf\",\"Zmumu\"\)
root -l -q computeAccSelZmmBinned.C+\(\"zmm.conf\",\"Zmumu\"\)

#
# Z->ee
#
#root -l -q computeAccGenZee.C+\(\"zee.conf\",\"Zee\"\)
#root -l -q computeAccSCZee.C+\(\"zee.conf\",\"Zee\"\)
#root -l -q computeAccSelZee.C+\(\"zee.conf\",\"Zee\"\)
root -l -q computeAccSelZeeBinned.C+\(\"zee.conf\",\"Zee\"\)

rm *.so *.d
