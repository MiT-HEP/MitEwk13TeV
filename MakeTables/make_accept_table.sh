#!/bin/bash

PDF_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/NNPDF30_amc

wpmPRE=`grep  Pre-FSR  ${PDF_DIR}/wpm_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`
wpmPOST=`grep Post-FSR ${PDF_DIR}/wpm_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`

echo "&" ${wpmPRE} " &" ${wpmPOST} " &" > temp_acc.txt

wmmPRE=`grep  Pre-FSR  ${PDF_DIR}/wmm_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`
wmmPOST=`grep Post-FSR ${PDF_DIR}/wmm_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`

echo "&" ${wmmPRE} " &" ${wmmPOST} " &" >> temp_acc.txt

zmmPRE=`grep  Pre-FSR  ${PDF_DIR}/zmm_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`
zmmPOST=`grep Post-FSR ${PDF_DIR}/zmm_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`

echo "&" ${zmmPRE} " &" ${zmmPOST} " &" >> temp_acc.txt

wpePRE=`grep  Pre-FSR  ${PDF_DIR}/wpe_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`
wpePOST=`grep Post-FSR ${PDF_DIR}/wpe_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`

echo "&" ${wpePRE} " &" ${wpePOST} " &" >> temp_acc.txt

wmePRE=`grep  Pre-FSR  ${PDF_DIR}/wme_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`
wmePOST=`grep Post-FSR ${PDF_DIR}/wme_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`

echo "&" ${wmePRE} " &" ${wmePOST} " &" >> temp_acc.txt

zeePRE=`grep  Pre-FSR  ${PDF_DIR}/zee_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`
zeePOST=`grep Post-FSR ${PDF_DIR}/zee_nnpdf30_nlo_as_0118.txt | head -n 1 | awk '{ print $3}'`

echo "&" ${zeePRE} " &" ${zeePOST} " &" >> temp_acc.txt

paste headers_for_acceptance.txt temp_acc.txt > acceptance_table.txt

rm temp_acc.txt
