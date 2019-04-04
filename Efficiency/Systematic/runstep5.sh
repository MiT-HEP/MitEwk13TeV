#!/bin/bash
lepton=$1
efftype=$2
charge=$3
# lepton
binvar=etapt

plotdir="/afs/cern.ch/work/s/sabrandt/public/SM/LowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/Systematic"

for ((binnum=0; binnum<10; binnum++))
#for binnum in 0 #2 3 5 6 8
do

# if [ ${binnum} -eq 2 ] || [ ${binnum} -eq 9 ] || [ ${binnum} -eq 14 ] || [ ${binnum} -eq 21 ] || [ ${binnum} -eq 26 ] || [ ${binnum} -eq 33 ]; then
	echo "0.0" >>${plotdir}/Results/Sig_Array.txt
	echo "," >>${plotdir}/Results/Sig_Array.txt
	# echo "0.0" >>${plotdir}/Results/Bkg_Array.txt
	# echo "," >>${plotdir}/Results/Bkg_Array.txt
	# continue
# fi

cat ${plotdir}/Results/Sig_${lepton}_${efftype}_${charge}_${binvar}_${binnum}.txt | awk '{print $1 * 100.}' >> ${plotdir}/Results/Sig_Array.txt
# cat ${plotdir}/Results/Bkg_${lepton}_${efftype}_${charge}_${binvar}_${binnum}.txt | awk '{print $1 * 100.}' >> ${plotdir}/Results/Bkg_Array.txt

echo "," >>${plotdir}/Results/Sig_Array.txt
# echo "," >>${plotdir}/Results/Bkg_Array.txt

done

