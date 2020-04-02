#! /bin/bash

OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance/FINAL_GEN/V


########################################
##             W->munu
########################################

# # ################################################################################
# ##########     aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -q computeAccGenWm.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wmp_13TeV_amcPythia_dressed\",\"ACC_wmp_13_amcPythia_dressed\",1,1\)
root -l -q computeAccGenWm.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wmm_13TeV_amcPythia_dressed\",\"ACC_wmm_13_amcPythia_dressed\",1,-1\)
root -l -q computeAccGenWm.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wm_13TeV_amcPythia_dressed\",\"ACC_wm_13_amcPythia_dressed\",1,0\)

# # ################################################################################
# ##########     aMCnlo+Pythia with pT reweight for RESUMMATION
root -l -q computeAccGenWm.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wmp_13TeV_amcPythia_ptWeight\",\"ACC_wmp_13_amcPythia_ptWeight\",1,1\)
root -l -q computeAccGenWm.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wmm_13TeV_amcPythia_ptWeight\",\"ACC_wmm_13_amcPythia_ptWeight\",1,-1\)
root -l -q computeAccGenWm.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wm_13TeV_amcPythia_ptWeight\",\"ACC_wm_13_amcPythia_ptWeight\",1,0\)

# # ################################################################################
# ##########     Powheg+Pythia for RESUMMATION
root -l -q computeAccGenWm.C+\(\"wm13_powheg.conf\",\"${OUTDIR}/GEN_wmp_13TeV_powPythia_dressed\",\"ACC_wmp_13_powPythia_dressed\",1,1\)
root -l -q computeAccGenWm.C+\(\"wm13_powheg.conf\",\"${OUTDIR}/GEN_wmm_13TeV_powPythia_dressed\",\"ACC_wmm_13_powPythia_dressed\",1,-1\)
root -l -q computeAccGenWm.C+\(\"wm13_powheg.conf\",\"${OUTDIR}/GEN_wm_13TeV_powPythia_dressed\",\"ACC_wm_13_powPythia_dressed\",1,0\)

# # ################################################################################
# ##########     Powheg+Photos for EWK
root -l -q computeAccGenWm.C+\(\"wm13_photos.conf\",\"${OUTDIR}/GEN_wmp_13TeV_powPhotos_undressed\",\"ACC_wmp_13_powPhotos_undressed\",0,1\)
root -l -q computeAccGenWm.C+\(\"wm13_photos.conf\",\"${OUTDIR}/GEN_wmm_13TeV_powPhotos_undressed\",\"ACC_wmm_13_powPhotos_undressed\",0,-1\)
root -l -q computeAccGenWm.C+\(\"wm13_photos.conf\",\"${OUTDIR}/GEN_wm_13TeV_powPhotos_undressed\",\"ACC_wm_13_powPhotos_undressed\",0,0\)


# # ################################################################################
# # ##########     Powheg+Pythia for EWK
root -l -q computeAccGenWm.C+\(\"wm13_pythia.conf\",\"${OUTDIR}/GEN_wmp_13TeV_powPythia_undressed\",\"ACC_wmp_13_powPythia_undressed\",0,1\)
root -l -q computeAccGenWm.C+\(\"wm13_pythia.conf\",\"${OUTDIR}/GEN_wmm_13TeV_powPythia_undressed\",\"ACC_wmm_13_powPythia_undressed\",0,-1\)
root -l -q computeAccGenWm.C+\(\"wm13_pythia.conf\",\"${OUTDIR}/GEN_wm_13TeV_powPythia_undressed\",\"ACC_wm_13_powPythia_undressed\",0,0\)


########################################
##             W->enu
########################################

# ################################################################################
##########     aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -q computeAccGenWe.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wep_13TeV_amcPythia_dressed\",\"ACC_wep_13_amcPythia_dressed\",1,1\)
root -l -q computeAccGenWe.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wem_13TeV_amcPythia_dressed\",\"ACC_wem_13_amcPythia_dressed\",1,-1\)
root -l -q computeAccGenWe.C+\(\"w13.conf\",\"${OUTDIR}/GEN_we_13TeV_amcPythia_dressed\",\"ACC_we_13_amcPythia_dressed\",1,0\)

# ################################################################################
##########     aMCnlo+Pythia for RESUMMATION with pT weight
root -l -q computeAccGenWe.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wep_13TeV_amcPythia_ptWeight\",\"ACC_wep_13_amcPythia_ptWeight\",1,1\)
root -l -q computeAccGenWe.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wem_13TeV_amcPythia_ptWeight\",\"ACC_wem_13_amcPythia_ptWeight\",1,-1\)
root -l -q computeAccGenWe.C+\(\"w13.conf\",\"${OUTDIR}/GEN_we_13TeV_amcPythia_ptWeight\",\"ACC_we_13_amcPythia_ptWeight\",1,0\)

# ################################################################################
##########     Powheg+Pythia for RESUMMATION
root -l -q computeAccGenWe.C+\(\"we13_powheg.conf\",\"${OUTDIR}/GEN_wep_13TeV_powPythia_dressed\",\"ACC_wep_13_powPythia_dressed\",1,1\)
root -l -q computeAccGenWe.C+\(\"we13_powheg.conf\",\"${OUTDIR}/GEN_wem_13TeV_powPythia_dressed\",\"ACC_wem_13_powPythia_dressed\",1,-1\)
root -l -q computeAccGenWe.C+\(\"we13_powheg.conf\",\"${OUTDIR}/GEN_we_13TeV_powPythia_dressed\",\"ACC_we_13_powPythia_dressed\",1,0\)

# ################################################################################
##########     Powheg+Photos for EWK
root -l -q computeAccGenWe.C+\(\"we13_photos.conf\",\"${OUTDIR}/GEN_wep_13TeV_powPhotos_undressed\",\"ACC_wep_13_powPhotos_undressed\",1,1\)
root -l -q computeAccGenWe.C+\(\"we13_photos.conf\",\"${OUTDIR}/GEN_wem_13TeV_powPhotos_undressed\",\"ACC_wem_13_powPhotos_undressed\",1,-1\)
root -l -q computeAccGenWe.C+\(\"we13_photos.conf\",\"${OUTDIR}/GEN_we_13TeV_powPhotos_undressed\",\"ACC_we_13_powPhotos_undressed\",1,0\)

# # ################################################################################
# ##########     Powheg+Pythia for EWK
root -l -q computeAccGenWe.C+\(\"we13_pythia.conf\",\"${OUTDIR}/GEN_wep_13TeV_powPythia_undressed\",\"ACC_wep_13_powPythia_undressed\",1,1\)
root -l -q computeAccGenWe.C+\(\"we13_pythia.conf\",\"${OUTDIR}/GEN_wem_13TeV_powPythia_undressed\",\"ACC_wem_13_powPythia_undressed\",1,-1\)
root -l -q computeAccGenWe.C+\(\"we13_pythia.conf\",\"${OUTDIR}/GEN_we_13TeV_powPythia_undressed\",\"ACC_we_13_powPythia_undressed\",1,0\)


########################################
##             Z->mm
########################################

# ################################################################################
# #########             aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -q computeAccGenZmm.C+\(\"z13.conf\",\"${OUTDIR}/GEN_zmm_13TeV_amcPythia_dressed\",\"ACC_zmm_13_amcPythia_dressed\",1\)
# ################################################################################
# #########             aMCnlo+Pythia for RESUMMATION with pT weights
root -l -q computeAccGenZmm.C+\(\"z13.conf\",\"${OUTDIR}/GEN_zmm_13TeV_amcPythia_ptWeight\",\"ACC_zmm_13_amcPythia_ptWeight\",1\)
# # ################################################################################
# # #########             Powheg+minlo+Pythia for RESUMMATION
root -l -q computeAccGenZmm.C+\(\"zmm13_minlo.conf\",\"${OUTDIR}/GEN_zmm_13TeV_powPythia_dressed\",\"ACC_zmm_13_powPythia_dressed\",1\)
# # ################################################################################
# # ##########            Powheg+Photos for EWK
root -l -q computeAccGenZmm.C+\(\"zmm13_photos.conf\",\"${OUTDIR}/GEN_zmm_13TeV_powPhotos_undressed\",\"ACC_zmm_13_powPhotos_undressed\",0\)
# # ################################################################################
# # #########             Powheg+Pythia for EWK
root -l -q computeAccGenZmm.C+\(\"zmm13_pythia.conf\",\"${OUTDIR}/GEN_zmm_13TeV_powPythia_undressed\",\"ACC_zmm_13_powPythia_undressed\",0\)


########################################
##             Z->ee
########################################
# ################################################################################
# # ##########             aMCnlo+Pythia for RESUMMATION
root -l -q computeAccGenZee.C+\(\"z13.conf\",\"${OUTDIR}/GEN_zee_13TeV_amcPythia_dressed\",\"ACC_zee_13_amcPythia_dressed\",1\)# ################################################################################
# # ##########             aMCnlo+Pythia for RESUMMATION
root -l -q computeAccGenZee.C+\(\"z13.conf\",\"${OUTDIR}/GEN_zee_13TeV_amcPythia_ptWeight\",\"ACC_zee_13_amcPythia_ptWeight\",1\)
# ################################################################################
# # ##########             Powheg+minlo+Pythia for RESUMMATION
root -l -q computeAccGenZee.C+\(\"zee13_minlo.conf\",\"${OUTDIR}/GEN_zee_13TeV_powPythia_dressed\",\"ACC_zee_13_powPythia_dressed\",1\)
################################################################################
# # ##########             Powheg+Photos for EWK
root -l -q computeAccGenZee.C+\(\"zee13_photos.conf\",\"${OUTDIR}/GEN_zee_13TeV_powPhotos_undressed\",\"ACC_zee_13_powPhotos_undressed\",0\)
################################################################################
# ##########             Powheg+Pythia for EWK
root -l -q computeAccGenZee.C+\(\"zee13_pythia.conf\",\"${OUTDIR}/GEN_zee_13TeV_powPythia_undressed\",\"ACC_zee_13_powPythia_undressed\",0\)

# rm *.so *.d
