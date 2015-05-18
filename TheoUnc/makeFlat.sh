#!/bin/bash

#indir=/store/user/jlawhorn/CT10nlo-PYTHIA6-13-v2/
indir=/store/user/jlawhorn/NNPDF30_PWHG2_RW/
outdir=/afs/cern.ch/work/j/jlawhorn/flat-05-14/

for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${indir} | grep bacon`
do
    if [[ -e ${outdir}${file} ]] 
    then 
	continue
    else
	echo root -l -q makeFlat.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}${file}\"\)
	root -l -q makeFlat.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}${file}\"\)
    fi
done