# workdir="/afs/cern.ch/work/x/xniu/public/Lumi/Ele/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
# FILEDIR="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"
# binnum=$1
# FILEDIR=/afs/cern.ch/work/s/sabrandt/public/LowPU_13TeV_Efficiency_v1/results/
FILEDIR=/afs/cern.ch/work/s/sabrandt/public/SM/LowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/testReweights_v2_2/results
# POSTFIX=_v1
# POSTFIX=_CBxBW_v1
POSTFIX=_POWxPythia_v1
POSTFIX_alt=_POWxPhotos_v1
#
# CMSSW_BASE="/afs/cern.ch/work/x/xniu/public/Lumi/Ele/CMSSW_7_6_3_patch2/src/"
TOP="$PWD"
#
BINVAR=etapt
FOLDER=Zmm
EFFTYPE=MuSITEff #MuHLTEff, MuSelEff, MuStaEff
NBINS=2 #Muons have 63 bins
CHARGE=Negative
NTOYS=30
binnum=1
toynum=0
# FILEDIR1=${FILEDIR}/${FOLDER}/Data/${EFFTYPE}_aMCxPythia${POSTFIX}/${CHARGE}/plots/
OUTPUTDIR=${FILEDIR}/${EFFTYPE}${POSTFIX}${POSTFIX_alt}_origtest/${CHARGE}/
# STAGEDIR=${TOP}/${EFFTYPE}/${CHARGE}/Step2Output/${POSTFIX}v${POSTFIX_alt}/
STAGEDIR=./
SIG_FILE=${FILEDIR}/${FOLDER}/Data/${EFFTYPE}${POSTFIX_alt}/${CHARGE}/plots/
BKG_FILE=${FILEDIR}/${FOLDER}/Data/${EFFTYPE}${POSTFIX}/${CHARGE}/plots/
MAS_FILE=${FILEDIR}/${FOLDER}/Data/${EFFTYPE}${POSTFIX}/${CHARGE}/plots/

# #
mkdir -p ${OUTPUTDIR}/
# while [ ${binnum} -lt ${NBINS} ] 
# do
# # echo ${binnum}
# # echo "HELLO"
# # echo "${FILEDIR}\/Zmm"
# # echo "${FILEDIR}/Zmm/Data/MuHLTEff_v1/Negative/plots/"
# # echo "${FILEDIR}/Zmm/Data/MuHLTEff${POSTFIX}/${CHARGE}/plots/"
# # echo \"${FILEDIR}/Zmm/Data/MuHLTEff${POSTFIX_CB}/${CHARGE}/plots/\"
# # echo \"${BINVAR}_${binnum}\"
# echo " ACTUALLY DOING MACRO NOW" 

root -l -b -q toyGenAndPull.C+\(\"${SIG_FILE}\",\"${BKG_FILE}\",\"${MAS_FILE}\",\"${BINVAR}_${binnum}\",\"${OUTPUTDIR}_asimov\",\"pull_${binnum}\",${binnum},${binnum},${NTOYS}\)
# root -l -b -q makePseudoData.C+\(\"${BKG_FILE}\",\"${BKG_FILE}\",\"${BINVAR}_${binnum}\",\"${STAGEDIR}\",1,${binnum},${NTOYS}\)
# # mkdir ${FILEDIR}/${EFFTYPE}${POSTFIX_CB}
# #### Loop to gather data from each of the toys
# while [ ${toynum} -lt ${NTOYS} ] 
# do
# root -l -q -b doStep3.C+\(\"${STAGEDIR}\",\"${BINVAR}_${binnum}_${toynum}.dat\",\"${BKG_FILE}/${BINVAR}_${binnum}.root\",\"${OUTPUTDIR}\",\"_${BINVAR}_${binnum}\"\)
# # rm ${STAGEDIR}/${BINVAR}_${binnum}_${toynum}.dat
# toynum=$[${toynum}+1]
# done
# end that loop

# binnum=$[${binnum}+1]
# echo $binnum
# done




