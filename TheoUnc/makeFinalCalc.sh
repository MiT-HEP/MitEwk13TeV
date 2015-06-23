#!/bin/bash

root -l -q getCT14uncertainties.C+ > ct14.txt
#root -l -q getNNPDF30uncertainties.C+ > nnpdf30.txt
#root -l -q getMMHT2014uncertainties.C+ > mmht2014.txt

mkdir for2Dplots

for pdf in mmht2014 ct14 nnpdf30
do
    for chan in zmm zee wpm wmm wpe wme
    do 
	cat ${pdf}.txt |grep ${chan} | grep Nominal | head -n 1 | awk '{ print $4 }' > for2Dplots/${chan}_${pdf}.txt
	cat ${pdf}.txt |grep ${chan} | grep -v Nominal | awk '{ print $3 }' >> for2Dplots/${chan}_${pdf}.txt
    done
done