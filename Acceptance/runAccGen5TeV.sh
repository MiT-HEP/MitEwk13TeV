#! /bin/bash

OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance/BLAH2

#
# # W->munu
# # #
# root -l -q computeAccGenWm_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WmpGen5TeV_Undressed\",0,1\)
# root -l -q computeAccGenWm_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WmmGen5TeV_Undressed\",0,-1\)
# root -l -q computeAccGenWm_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WmGen5TeV_Undressed\",0,0\)

# root -l -q computeAccGenWm_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WmpGen5TeV_Dressed\",1,1\)
# root -l -q computeAccGenWm_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WmmGen5TeV_Dressed\",1,-1\)
# root -l -q computeAccGenWm_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WmGen5TeV_Dressed\",1,0\)


# # W->enu
# #
# root -l -q computeAccGenWe_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WepGen5TeV_Undressed\",0,1\)
# root -l -q computeAccGenWe_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WemGen5TeV_Undressed\",0,-1\)
# root -l -q computeAccGenWe_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WeGen5TeV_Undressed\",0,0\)

# root -l -q computeAccGenWe_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WepGen5TeV_Dressed\",1,1\)
# root -l -q computeAccGenWe_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WemGen5TeV_Dressed\",1,-1\)
# root -l -q computeAccGenWe_Sys.C+\(\"w_5.conf\",\"${OUTDIR}/WeGen5TeV_Dressed\",1,0\)

#
# Z->mumu
#
# root -l -q computeAccGenZmm_Sys.C+\(\"z_5.conf\",\"${OUTDIR}/ZmmGen5TeV_Undressed\",0\)
# root -l -q computeAccGenZmm_Sys.C+\(\"z_5.conf\",\"${OUTDIR}/ZmmGen5TeV_Dressed\",1\)

#
# Z->ee
#
# root -l -q computeAccGenZee_Sys.C+\(\"z_5.conf\",\"${OUTDIR}/ZeeGen5TeV_Undressed\",0\)
# root -l -q computeAccGenZee_Sys.C+\(\"z_5.conf\",\"${OUTDIR}/ZeeGen5TeV_Dressed\",1\)

# rm *.so *.d
