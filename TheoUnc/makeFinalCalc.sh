#!/bin/bash

ACC_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-acc

for CHAN in wme wpe wmm wpm zee zmm
do
    awk '/Tot.*/ {print $2}' ${ACC_DIR}/${CHAN}_ct10nlo.txt > ${ACC_DIR}/${CHAN}_parse_ct10nlo.txt
    cat ${ACC_DIR}/${CHAN}_ct10nlo_as_01* | awk '/Tot.*/ {print $3}' > ${ACC_DIR}/${CHAN}_parse_ct10nlo_as.txt

    awk '/Tot.*/ {print $2 }' ${ACC_DIR}/${CHAN}_nnpdf23_nlo_as_0119.txt > ${ACC_DIR}/${CHAN}_parse_nnpdf23_nlo.txt
    awk '/Tot.*/ {print $2 }' ${ACC_DIR}/${CHAN}_nnpdf23_nlo_as_0118.txt >> ${ACC_DIR}/${CHAN}_parse_nnpdf23_nlo.txt
    awk '/Tot.*/ {print $2 }' ${ACC_DIR}/${CHAN}_nnpdf23_nlo_as_0117.txt >> ${ACC_DIR}/${CHAN}_parse_nnpdf23_nlo.txt
    awk '/Tot.*/ {print $2 }' ${ACC_DIR}/${CHAN}_nnpdf23_nlo_as_0116.txt >> ${ACC_DIR}/${CHAN}_parse_nnpdf23_nlo.txt
    awk '/Tot.*/ {print $2 }' ${ACC_DIR}/${CHAN}_nnpdf23_nlo_as_0120.txt >> ${ACC_DIR}/${CHAN}_parse_nnpdf23_nlo.txt
    awk '/Tot.*/ {print $2 }' ${ACC_DIR}/${CHAN}_nnpdf23_nlo_as_0121.txt >> ${ACC_DIR}/${CHAN}_parse_nnpdf23_nlo.txt
    awk '/Tot.*/ {print $2 }' ${ACC_DIR}/${CHAN}_nnpdf23_nlo_as_0122.txt >> ${ACC_DIR}/${CHAN}_parse_nnpdf23_nlo.txt

    awk '/Tot.*/ {print $3 }' ${ACC_DIR}/${CHAN}_mstw2008nlo68cl.txt > ${ACC_DIR}/${CHAN}_parse_mstw2008nlo68cl.txt
    awk '/Tot.*/ {print $3 }' ${ACC_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68cl.txt > ${ACC_DIR}/${CHAN}_parse_mstw2008nlo68cl_asmz+68cl.txt
    awk '/Tot.*/ {print $3 }' ${ACC_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68clhalf.txt > ${ACC_DIR}/${CHAN}_parse_mstw2008nlo68cl_asmz+68clhalf.txt
    awk '/Tot.*/ {print $3 }' ${ACC_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68cl.txt > ${ACC_DIR}/${CHAN}_parse_mstw2008nlo68cl_asmz-68cl.txt
    awk '/Tot.*/ {print $3 }' ${ACC_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68clhalf.txt > ${ACC_DIR}/${CHAN}_parse_mstw2008nlo68cl_asmz-68clhalf.txt

done

root -l -q getMSTW2008uncertainties.C+\(\"${ACC_DIR}\"\) > ${ACC_DIR}/mstw_uncert.txt
root -l -q getNNPDF23uncertainties.C+\(\"${ACC_DIR}\"\)  > ${ACC_DIR}/nnpdf_uncert.txt
root -l -q getCT10uncertainties.C+\(\"${ACC_DIR}\"\)     > ${ACC_DIR}/cteq_uncert.txt