#! /bin/bash

OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance/FINAL_GEN/V


########################################
##             W->munu
########################################

# # ################################################################################
# ##########     aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -q computeAccGenWm.C+\(\"w5.conf\",\"${OUTDIR}/GEN_wmp_5TeV_amcPythia_dressed\",\"ACC_wmp_5_amcPythia_dressed\",1,1\)
root -l -q computeAccGenWm.C+\(\"w5.conf\",\"${OUTDIR}/GEN_wmm_5TeV_amcPythia_dressed\",\"ACC_wmm_5_amcPythia_dressed\",1,-1\)
root -l -q computeAccGenWm.C+\(\"w5.conf\",\"${OUTDIR}/GEN_wm_5TeV_amcPythia_dressed\",\"ACC_wm_5_amcPythia_dressed\",1,0\)

# # ################################################################################
# ##########     aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -q computeAccGenWm.C+\(\"w5.conf\",\"${OUTDIR}/GEN_wmp_5TeV_amcPythia_ptWeight\",\"ACC_wmp_5_amcPythia_ptWeight\",1,1\)
root -l -q computeAccGenWm.C+\(\"w5.conf\",\"${OUTDIR}/GEN_wmm_5TeV_amcPythia_ptWeight\",\"ACC_wmm_5_amcPythia_ptWeight\",1,-1\)
root -l -q computeAccGenWm.C+\(\"w5.conf\",\"${OUTDIR}/GEN_wm_5TeV_amcPythia_ptWeight\",\"ACC_wm_5_amcPythia_ptWeight\",1,0\)

# # ################################################################################
# ##########     Powheg+Pythia for RESUMMATION
root -l -q computeAccGenWm.C+\(\"wm5_powheg.conf\",\"${OUTDIR}/GEN_wmp_5TeV_powPythia_dressed\",\"ACC_wmp_5_powPythia_dressed\",1,1\)
root -l -q computeAccGenWm.C+\(\"wm5_powheg.conf\",\"${OUTDIR}/GEN_wmm_5TeV_powPythia_dressed\",\"ACC_wmm_5_powPythia_dressed\",1,-1\)
root -l -q computeAccGenWm.C+\(\"wm5_powheg.conf\",\"${OUTDIR}/GEN_wm_5TeV_powPythia_dressed\",\"ACC_wm_5_powPythia_dressed\",1,0\)

# # ################################################################################
# ##########     Powheg+Photos for EWK
root -l -q computeAccGenWm.C+\(\"wm5_photos.conf\",\"${OUTDIR}/GEN_wmp_5TeV_powPhotos_undressed\",\"ACC_wmp_5_powPhotos_undressed\",0,1\)
root -l -q computeAccGenWm.C+\(\"wm5_photos.conf\",\"${OUTDIR}/GEN_wmm_5TeV_powPhotos_undressed\",\"ACC_wmm_5_powPhotos_undressed\",0,-1\)
root -l -q computeAccGenWm.C+\(\"wm5_photos.conf\",\"${OUTDIR}/GEN_wm_5TeV_powPhotos_undressed\",\"ACC_wm_5_powPhotos_undressed\",0,0\)


# # ################################################################################
# ##########     Powheg+Pythia for EWK
root -l -q computeAccGenWm.C+\(\"wm5_pythia.conf\",\"${OUTDIR}/GEN_wmp_5TeV_powPythia_undressed\",\"ACC_wmp_5_powPythia_undressed\",0,1\)
root -l -q computeAccGenWm.C+\(\"wm5_pythia.conf\",\"${OUTDIR}/GEN_wmm_5TeV_powPythia_undressed\",\"ACC_wmm_5_powPythia_undressed\",0,-1\)
root -l -q computeAccGenWm.C+\(\"wm5_pythia.conf\",\"${OUTDIR}/GEN_wm_5TeV_powPythia_undressed\",\"ACC_wm_5_powPythia_undressed\",0,0\)


########################################
##             W->enu
########################################

# ################################################################################
##########     aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -q computeAccGenWe.C+\(\"w5.conf\",\"${OUTDIR}/GEN_wep_5TeV_amcPythia_dressed\",\"ACC_wep_5_amcPythia_dressed\",1,1\)
root -l -q computeAccGenWe.C+\(\"w5.conf\",\"${OUTDIR}/GEN_wem_5TeV_amcPythia_dressed\",\"ACC_wem_5_amcPythia_dressed\",1,-1\)
root -l -q computeAccGenWe.C+\(\"w5.conf\",\"${OUTDIR}/GEN_we_5TeV_amcPythia_dressed\",\"ACC_we_5_amcPythia_dressed\",1,0\)

# ################################################################################
##########     aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -q computeAccGenWe.C+\(\"w5.conf\",\"${OUTDIR}/GEN_wep_5TeV_amcPythia_ptWeight\",\"ACC_wep_5_amcPythia_ptWeight\",1,1\)
root -l -q computeAccGenWe.C+\(\"w5.conf\",\"${OUTDIR}/GEN_wem_5TeV_amcPythia_ptWeight\",\"ACC_wem_5_amcPythia_ptWeight\",1,-1\)
root -l -q computeAccGenWe.C+\(\"w5.conf\",\"${OUTDIR}/GEN_we_5TeV_amcPythia_ptWeight\",\"ACC_we_5_amcPythia_ptWeight\",1,0\)

# ################################################################################
##########     Powheg+Pythia for RESUMMATION
root -l -q computeAccGenWe.C+\(\"we5_powheg.conf\",\"${OUTDIR}/GEN_wep_5TeV_powPythia_dressed\",\"ACC_wep_5_powPythia_dressed\",1,1\)
root -l -q computeAccGenWe.C+\(\"we5_powheg.conf\",\"${OUTDIR}/GEN_wem_5TeV_powPythia_dressed\",\"ACC_wem_5_powPythia_dressed\",1,-1\)
root -l -q computeAccGenWe.C+\(\"we5_powheg.conf\",\"${OUTDIR}/GEN_we_5TeV_powPythia_dressed\",\"ACC_we_5_powPythia_dressed\",1,0\)

# ################################################################################
##########     Powheg+Photos for EWK
root -l -q computeAccGenWe.C+\(\"we5_photos.conf\",\"${OUTDIR}/GEN_wep_5TeV_powPhotos_undressed\",\"ACC_wep_5_powPhotos_undressed\",1,1\)
root -l -q computeAccGenWe.C+\(\"we5_photos.conf\",\"${OUTDIR}/GEN_wem_5TeV_powPhotos_undressed\",\"ACC_wem_5_powPhotos_undressed\",1,-1\)
root -l -q computeAccGenWe.C+\(\"we5_photos.conf\",\"${OUTDIR}/GEN_we_5TeV_powPhotos_undressed\",\"ACC_we_5_powPhotos_undressed\",1,0\)

# ################################################################################
##########     Powheg+Pythia for EWK
root -l -q computeAccGenWe.C+\(\"we5_pythia.conf\",\"${OUTDIR}/GEN_wep_5TeV_powPythia_undressed\",\"ACC_wep_5_powPythia_undressed\",1,1\)
root -l -q computeAccGenWe.C+\(\"we5_pythia.conf\",\"${OUTDIR}/GEN_wem_5TeV_powPythia_undressed\",\"ACC_wem_5_powPythia_undressed\",1,-1\)
root -l -q computeAccGenWe.C+\(\"we5_pythia.conf\",\"${OUTDIR}/GEN_we_5TeV_powPythia_undressed\",\"ACC_we_5_powPythia_undressed\",1,0\)



########################################
##             Z->mm
########################################

# ################################################################################
# #########             aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -q computeAccGenZmm.C+\(\"z5.conf\",\"${OUTDIR}/GEN_zmm_5TeV_amcPythia_dressed\",\"ACC_zmm_5_amcPythia_dressed\",1\)
# ################################################################################
# #########             aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -q computeAccGenZmm.C+\(\"z5.conf\",\"${OUTDIR}/GEN_zmm_5TeV_amcPythia_ptWeight\",\"ACC_zmm_5_amcPythia_ptWeight\",1\)
# # ################################################################################
# # #########             Powheg+minlo+Pythia for RESUMMATION
root -l -q computeAccGenZmm.C+\(\"zmm5_minlo.conf\",\"${OUTDIR}/GEN_zmm_5TeV_powPythia_dressed\",\"ACC_zmm_5_powPythia_dressed\",1\)
# ################################################################################
# # ##########            Powheg+Photos for EWK
root -l -q computeAccGenZmm.C+\(\"zmm5_photos.conf\",\"${OUTDIR}/GEN_zmm_5TeV_powPhotos_undressed\",\"ACC_zmm_5_powPhotos_undressed\",0\)
# # ################################################################################
# # #########             Powheg+Pythia for EWK
root -l -q computeAccGenZmm.C+\(\"zmm5_pythia.conf\",\"${OUTDIR}/GEN_zmm_5TeV_powPythia_undressed\",\"ACC_zmm_5_powPythia_undressed\",0\)


########################################
##             Z->ee
########################################
# ################################################################################
# # ##########             aMCnlo+Pythia for RESUMMATION
root -l -q -b computeAccGenZee.C+\(\"z5.conf\",\"${OUTDIR}/GEN_zee_5TeV_amcPythia_dressed\",\"ACC_zee_5_amcPythia_dressed\",1\)
# ################################################################################
# # ##########             aMCnlo+Pythia for RESUMMATION
root -l -q -b computeAccGenZee.C+\(\"z5.conf\",\"${OUTDIR}/GEN_zee_5TeV_amcPythia_ptWeight\",\"ACC_zee_5_amcPythia_ptWeight\",1\)
# ################################################################################
# # ##########             Powheg+minlo+Pythia for RESUMMATION
root -l -q -b computeAccGenZee.C+\(\"zee5_minlo.conf\",\"${OUTDIR}/GEN_zee_5TeV_powPythia_dressed\",\"ACC_zee_5_powPythia_dressed\",1\)
################################################################################
# # ##########             Powheg+Photos for EWK
root -l -q -b computeAccGenZee.C+\(\"zee5_photos.conf\",\"${OUTDIR}/GEN_zee_5TeV_powPhotos_undressed\",\"ACC_zee_5_powPhotos_undressed\",0\)
################################################################################
# ##########             Powheg+Pythia for EWK
root -l -q -b computeAccGenZee.C+\(\"zee5_pythia.conf\",\"${OUTDIR}/GEN_zee_5TeV_powPythia_undressed\",\"ACC_zee_5_powPythia_undressed\",0\)

# rm *.so *.d
