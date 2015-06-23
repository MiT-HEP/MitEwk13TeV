#!/bin/bash

#indir=/store/user/jlawhorn/CT10nlo-PYTHIA6-13-v2/
indir=/store/user/jlawhorn/NNPDF30_PWHG2_RW/
#indir=/store/user/jlawhorn/NNPDF30_PWHG2_RW/
#indir=/store/user/jlawhorn/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_Zjets_genbacon_7_4_t4/150515_074452/0000/
outdir=/afs/cern.ch/work/j/jlawhorn/

i=0

#for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${indir} | grep root`
for file in ../../BaconProd/crab/Output.root
#/store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50000/7C480AB4-96DC-E411-B0C6-60EB69BAC786.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50000/FC4C561E-98DC-E411-99B8-60EB69BAC930.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/0287622E-A5DC-E411-82B3-0025905506E8.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/265AD41B-9FDC-E411-9D29-0025905505DE.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/28B63FFF-9CDC-E411-B3F2-60EB69BAC82E.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/2C321AF3-99DC-E411-A546-60EB69BAC786.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/46F6F159-A0DC-E411-93DD-60EB69BACA5A.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/68ED51E2-9EDC-E411-B631-008CFA197C10.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/7488F2E0-A3DC-E411-A0EC-00259055C876.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/76AC4AD3-A1DC-E411-A4E6-0025904B6FF6.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/800381F1-9BDC-E411-AA5D-60EB69BACA8C.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/9A69D560-A8DC-E411-8AFA-60EB69BAC78C.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/9C31FAEF-A2DC-E411-B9E6-00259055C876.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/AE193552-A7DC-E411-856D-0025905521D2.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/CC414016-A6DC-E411-A718-002590574512.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/D8319A38-9CDC-E411-BD67-60EB69BAC82E.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/E63511F8-99DC-E411-A7AB-60EB69BACBB4.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/F682B82C-9DDC-E411-A382-60EB69BACBC0.root /store/mc/RunIIWinter15GS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM/MCRUN2_71_V1-v1/50001/F8A9AB35-9BDC-E411-BA7B-60EB69BACA84.root
do
 #   if [[ -e ${outdir}${file} ]] 
 #   then 
#	continue
 #   else
#	echo ${file}
    #echo root -l -q makeFlat.C+\(\"root://cmsxrootd.fnal.gov/${file}\",\"${outdir}${file}\",23\)
    root -l -q makeFlat.C+\(\"${file}\",\"${outdir}file${i}.root\",23\)
    i=$((i+1))
 #   fi
done