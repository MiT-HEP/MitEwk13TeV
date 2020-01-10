#! /bin/bash

OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance/GenAcc13_Powheg_Unc

## minlo and powheg
## Undressed
# root -l -q computeAccGenWm_Sys.C+\(\"wmp_13_powheg.conf\",\"${OUTDIR}/WmpGen13TeV_powheg_Undressed\",0,1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wmm_13_powheg.conf\",\"${OUTDIR}/WmmGen13TeV_powheg_Undressed\",0,-1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WmGen13TeV_powheg_Undressed\",0,0\)

## Dressed
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WmpGen13TeV_Dressed\",1,1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WmmGen13TeV_Dressed\",1,-1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WmGen13TeV_Dressed\",1,0\)


# # W->enu
# # #
# root -l -q computeAccGenWe_Sys.C+\(\"wep_13_powheg.conf\",\"${OUTDIR}/WepGen13TeV_powheg_Undressed\",0,1\)
# root -l -q computeAccGenWe_Sys.C+\(\"wem_13_powheg.conf\",\"${OUTDIR}/WemGen13TeV_powheg_Undressed\",0,-1\)
# root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WeGen13TeV_Undressed\",0,0\)

# root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WepGen13TeV_Dressed\",1,1\)
# root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WemGen13TeV_Dressed\",1,-1\)
# root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WeGen13TeV_Dressed\",1,0\)

#
# Z->mm
#

root -l -q computeAccGenZmm_Sys.C+\(\"zmm_13_minlo.conf\",\"${OUTDIR}/ZmmGen13TeV_minlo_Undressed\",0\)
# root -l -q computeAccGenZmm_Sys.C+\(\"zee_13_minlo.conf\",\"${OUTDIR}/ZmmGen13TeV_Minlo_UnDressed_2017\",0\)
# root -l -q computeAccGenZmm_Sys.C+\(\"z_13.conf\",\"${OUTDIR}/ZmmGen13TeV_UnDressed_2017\",0\)
# root -l -q computeAccGenZmm_Sys.C+\(\"zmm_2015.conf\",\"${OUTDIR}/ZmmGen13TeV_UnDressed_2015\",0\)

#
# Z->ee
#

# root -l -q computeAccGenZee_Sys.C+\(\"zee_13_minlo.conf\",\"${OUTDIR}/ZeeGen13TeV_minlo_Undressed\",0\)
# root -l -q computeAccGenZee_Sys.C+\(\"z_13.conf\",\"${OUTDIR}/ZeeGen13TeV_Dressed\",1\)



#########################################
###   regular aMC@nlo
#########################################



#
# W->munu
# #
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WmpGen13TeV_Undressed\",0,1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WmmGen13TeV_Undressed\",0,-1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WmGen13TeV_Undressed\",0,0\)

# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WmpGen13TeV_Dressed\",1,1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WmmGen13TeV_Dressed\",1,-1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WmGen13TeV_Dressed\",1,0\)


# # W->enu
# #
# root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WepGen13TeV_Undressed\",0,1\)
# root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WemGen13TeV_Undressed\",0,-1\)
# root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WeGen13TeV_Undressed\",0,0\)

# root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WepGen13TeV_Dressed\",1,1\)
# root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WemGen13TeV_Dressed\",1,-1\)
# root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/WeGen13TeV_Dressed\",1,0\)

#
# Z->mm
#

# root -l -q computeAccGenZmm_Sys.C+\(\"z_13.conf\",\"${OUTDIR}/ZmmGen13TeV_Undressed\",0\)
# root -l -q computeAccGenZmm_Sys.C+\(\"zee_13_minlo.conf\",\"${OUTDIR}/ZmmGen13TeV_Minlo_UnDressed_2017\",0\)
# root -l -q computeAccGenZmm_Sys.C+\(\"z_13.conf\",\"${OUTDIR}/ZmmGen13TeV_UnDressed_2017\",0\)
# root -l -q computeAccGenZmm_Sys.C+\(\"zmm_2015.conf\",\"${OUTDIR}/ZmmGen13TeV_UnDressed_2015\",0\)

#
# Z->ee
#

# root -l -q computeAccGenZee_Sys.C+\(\"z_13.conf\",\"${OUTDIR}/ZeeGen13TeV_Undressed\",0\)
# root -l -q computeAccGenZee_Sys.C+\(\"z_13.conf\",\"${OUTDIR}/ZeeGen13TeV_Dressed\",1\)


# rm *.so *.d
