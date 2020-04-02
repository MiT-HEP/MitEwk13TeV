#!/bin/sh

DIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance/FINAL_GEN/"
TEMP="v1"
WORK=${DIR}/${TEMP}
RESULTS="result_"${TEMP}
# s=13

# # fsr

# # mkdir -p ${WORK}
# cp ${DIR}/GEN_z*/ACC*.txt ${WORK}
# # # # # ## Combine any of the W files which were computed separately
# # ./combineWjet.sh 13 "dressed"    "amcPythia" ${DIR}
# # ./combineWjet.sh 13 "ptWeight"    "amcPythia" ${DIR}
# # ./combineWjet.sh 13 "dressed"    "powPythia" ${DIR}
# ./combineWjet.sh 13 "undressed"  "powPythia" ${DIR}
# ./combineWjet.sh 13 "undressed"  "powPhotos" ${DIR}


# # # ## Make the W+/W-, W+/Z, W-/Z, and W/Z ratios
# mv ${DIR}/ACC_*.txt ${WORK}
# # ./makeRatios.sh 13 "dressed"    "amcPythia" ${DIR} ${TEMP}
# # ./makeRatios.sh 13 "ptWeight"    "amcPythia" ${DIR} ${TEMP}
# # ./makeRatios.sh 13 "dressed"    "powPythia" ${DIR} ${TEMP}
# ./makeRatios.sh 13 "undressed"  "powPythia" ${DIR} ${TEMP}
# ./makeRatios.sh 13 "undressed"  "powPhotos" ${DIR} ${TEMP}


# ## start computing different things
# mv ${DIR}/ACC_*.txt ${WORK}
# ./directCompare.sh 13 "ewk" "mu"   ${DIR} ${TEMP}
# ./directCompare.sh 13 "res" "mu"   ${DIR} ${TEMP}
# ./directCompare.sh 13 "ewk" "ele"  ${DIR} ${TEMP}
# ./directCompare.sh 13 "res" "ele"  ${DIR} ${TEMP}

# # ./computePDF.sh
# ./computePDF.sh 13 "mu"   ${DIR} ${TEMP}
# ./computePDF.sh 13 "ele"  ${DIR} ${TEMP}

# # ./computeQCD.sh
# ./computeQCD.sh 13 "mu"   ${DIR} ${TEMP}
# ./computeQCD.sh 13 "ele"  ${DIR} ${TEMP}

# mkdir -p ${DIR}/${RESULTS}
# mv ${DIR}/RES*.txt ${DIR}/EWK*.txt ${DIR}/QCD*.txt ${DIR}/PDF*.txt ${DIR}/${RESULTS}
# # # ### Make the table
./fillTable.sh 13 "mu" "_13" ${DIR} ${RESULTS}
./fillTable.sh 13 "ele" "_13" ${DIR} ${RESULTS}

# s=5

# # fsr

# mkdir -p ${WORK}
# cp ${DIR}/GEN_z*/ACC*.txt ${WORK}
# # # ## Combine any of the W files which were computed separately
# ./combineWjet.sh 5 "dressed"    "amcPythia" ${DIR}
# ./combineWjet.sh 5 "ptWeight"    "amcPythia" ${DIR}
# # ./combineWjet.sh 5 "dressed"    "powPythia" ${DIR}
# # ./combineWjet.sh 5 "undressed"  "powPythia" ${DIR}
# # ./combineWjet.sh 5 "undressed"  "powPhotos" ${DIR}


# # # # ## Make the W+/W-, W+/Z, W-/Z, and W/Z ratios
# mv ${DIR}/ACC_*.txt ${WORK}
# ./makeRatios.sh 5 "dressed"    "amcPythia" ${DIR} ${TEMP}
# ./makeRatios.sh 5 "ptWeight"    "amcPythia" ${DIR} ${TEMP}
# # # ./makeRatios.sh 13 "dressed"    "powPythia" ${DIR} ${TEMP}
# # # ./makeRatios.sh 13 "undressed"  "powPythia" ${DIR} ${TEMP}
# # # ./makeRatios.sh 13 "undressed"  "powPhotos" ${DIR} ${TEMP}


# # ## start computing different things
# mv ${DIR}/ACC_*.txt ${WORK}
# # ./directCompare.sh 13 "ewk" "mu"   ${DIR} ${TEMP}
# ./directCompare.sh 5 "res" "mu"   ${DIR} ${TEMP}
# # ./directCompare.sh 13 "ewk" "ele"  ${DIR} ${TEMP}
# ./directCompare.sh 5 "res" "ele"  ${DIR} ${TEMP}

# # # ./computePDF.sh
# ./computePDF.sh 5 "mu"   ${DIR} ${TEMP}
# ./computePDF.sh 5 "ele"  ${DIR} ${TEMP}

# # # ./computeQCD.sh
# ./computeQCD.sh 5 "mu"   ${DIR} ${TEMP}
# ./computeQCD.sh 5 "ele"  ${DIR} ${TEMP}

# mkdir -p ${DIR}/${RESULTS}
# mv ${DIR}/RES*.txt ${DIR}/EWK*.txt ${DIR}/QCD*.txt ${DIR}/PDF*.txt ${DIR}/${RESULTS}
# # # ### Make the table
./fillTable.sh 5 "mu" "_5" ${DIR} ${RESULTS}
./fillTable.sh 5 "ele" "_5" ${DIR} ${RESULTS}

