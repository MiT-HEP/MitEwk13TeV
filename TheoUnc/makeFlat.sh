#!/bin/bash

#indir=/store/user/jlawhorn/CT10nlo-PYTHIA6-13-v2/
indir=/store/user/jlawhorn/NNPDF30-GENSIM/
outdir=/afs/cern.ch/work/j/jlawhorn/baseComp-05-06/

#for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${indir} | grep bacon`
for file in wmm-pythia6-powheg1-bacon.root wmm-pythia8-powheg1-bacon.root wmm-tune4C-powheg1-bacon.root wmm-tuneMonash-bacon.root wmm-tuned6t-bacon.root wmmLOEWK-pythia8-powheg2-bacon.root wmmNLOEWK-pythia8-powheg2-bacon.root
do
    if [[ -e ${outdir}${file} ]] 
    then 
	continue
    else
	echo root -l -q makeFlat.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}${file}\"\)
	root -l -q makeFlat.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}${file}\"\)
    fi
done

#for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${indir} | grep bacon | grep -v z`
#do
#    if [[ -e ${outdir}CT10-${file} ]]
#    then
#	continue
#    else 
#	echo root -l -q makeFlatW.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}CT10-${file}\"\)
#	root -l -q makeFlatW.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}CT10-${file}\"\)
#    fi
#done

#indir=/store/user/jlawhorn/aMCNLO/DYJetsToLL_500000/

#for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${indir}`
#do
#    if [[ -e ${outdir}${file} ]] 
#    then 
#	continue
#    else
	#echo root -l -q makeFlatW.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}${file}\"\)
	#root -l -q makeFlatW.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}${file}\"\)
#	echo root -l -q makeFlatZ.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}DYJets-${file}\"\)
#	root -l -q makeFlatZ.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}DYJets-${file}\"\)
#    fi
#done

#indir=/store/user/jlawhorn/aMCNLO/WJetsToLNu_50001/

#for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${indir}`
#do
#    if [[ -e ${outdir}${file} ]]
#    then
#	continue
#    else 
#	echo root -l -q makeFlatW.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}WJets${file}\"\)
#	root -l -q makeFlatW.C+\(\"root://eoscms/${indir}${file}\",\"${outdir}WJets${file}\"\)
#    fi
#done