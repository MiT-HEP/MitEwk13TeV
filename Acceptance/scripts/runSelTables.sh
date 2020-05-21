#!/bin/sh

DIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance/FINAL_SEL/TEST"
TEMP="v0_5"
WORK=${DIR}/${TEMP}
RESULTS="result_"${TEMP}
s=13

# fsr

# mkdir -p ${WORK}
# cp ${DIR}/SEL_z*/ACC*.txt ${WORK}
# # ## Combine any of the W files which were computed separately
# ./combine_Sel.sh 13 ""  ${DIR}

# # ## Make the W+/W-, W+/Z, W-/Z, and W/Z ratios
# mv ${DIR}/ACC_*.txt ${WORK}
# ./selRatio.sh 13  ${DIR} ${TEMP}


# ## start computing different things
# mv ${DIR}/ACC_*.txt ${WORK}
# ./selComp.sh 13  "bkg"   "mu"   ${DIR} ${TEMP}
# ./selComp.sh 13  "fsr"   "mu"   ${DIR} ${TEMP}
# ./selComp.sh 13  "mc"    "mu"   ${DIR} ${TEMP}
# ./selComp.sh 13  "tagpt" "mu"   ${DIR} ${TEMP}
# ./selComp.sh 13  "stat"  "mu"   ${DIR} ${TEMP}

# ./selComp.sh 13  "bkg"   "ele"   ${DIR} ${TEMP}
# ./selComp.sh 13  "fsr"   "ele"   ${DIR} ${TEMP}
# ./selComp.sh 13  "mc"    "ele"   ${DIR} ${TEMP}
# ./selComp.sh 13  "tagpt" "ele"   ${DIR} ${TEMP}
# ./selComp.sh 13  "stat"  "ele"   ${DIR} ${TEMP}

# # ./computePDF.sh
# ./computePDF.sh 13 "mu"   ${DIR} ${TEMP}
# ./computePDF.sh 13 "ele"  ${DIR} ${TEMP}

# # ./computeQCD.sh
# ./computeQCD.sh 13 "mu"   ${DIR} ${TEMP}
# ./computeQCD.sh 13 "ele"  ${DIR} ${TEMP}

# mkdir -p ${DIR}/${RESULTS}
# mv ${DIR}/BKG*.txt ${DIR}/TAG*.txt ${DIR}/MC*.txt ${DIR}/STAT*.txt ${DIR}/FSR*.txt ${DIR}/${RESULTS}
# # ### Make the table
# ./fillSel.sh 13 "mu" "_acc13" ${DIR} ${RESULTS}
# ./fillSel.sh 13 "ele" "_acc13" ${DIR} ${RESULTS}


# mkdir -p ${WORK}
# cp ${DIR}/SEL_*5*/ACC*.txt ${WORK}
# # ## Combine any of the W files which were computed separately
# # ./combine_Sel.sh 5 ""  ${DIR}

# # ## Make the W+/W-, W+/Z, W-/Z, and W/Z ratios
# mv ${DIR}/ACC_*.txt ${WORK}
# ./selRatio.sh 5  ${DIR} ${TEMP}


# ## start computing different things
# mv ${DIR}/ACC_*.txt ${WORK}
# ./selComp.sh 5  "bkg"   "mu"   ${DIR} ${TEMP}
# ./selComp.sh 5  "fsr"   "mu"   ${DIR} ${TEMP}
# ./selComp.sh 5  "mc"    "mu"   ${DIR} ${TEMP}
# ./selComp.sh 5  "tagpt" "mu"   ${DIR} ${TEMP}
# ./selComp.sh 5  "stat"  "mu"   ${DIR} ${TEMP}

# ./selComp.sh 5  "bkg"   "ele"   ${DIR} ${TEMP}
# ./selComp.sh 5  "fsr"   "ele"   ${DIR} ${TEMP}
# ./selComp.sh 5  "mc"    "ele"   ${DIR} ${TEMP}
# ./selComp.sh 5  "tagpt" "ele"   ${DIR} ${TEMP}
# ./selComp.sh 5  "stat"  "ele"   ${DIR} ${TEMP}


# mkdir -p ${DIR}/${RESULTS}
# mv ${DIR}/BKG*.txt ${DIR}/TAG*.txt ${DIR}/MC*.txt ${DIR}/STAT*.txt ${DIR}/FSR*.txt ${DIR}/${RESULTS}
# # ### Make the table
# ./fillSel.sh 5 "mu" "_acc5" ${DIR} ${RESULTS}
./fillSel.sh 5 "ele" "_acc5" ${DIR} ${RESULTS}