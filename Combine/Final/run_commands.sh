#!/bin/bash
# 
# python ../../../HiggsAnalysis/CombinedLimit/scripts/text2workspace.py Wmunu.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:wModelmu -o Wmunu.root
# combine -M MaxLikelihoodFit  --robustFit 1 --plots --rMax 2 --saveNormalizations Wmunu.root
# # combine -M MaxLikelihoodFit --robustFit 1 --rMax 2 --minimizerAlgo=Minuit2 --minimizerStrategy=1  Wmunu.root 
# mv mlfit.root mlfit_munu.root
# python ../../../HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py mlfit_munu.root
# python diffNuisances.py -A -a  -f text mlfit_munu.root > mlfit_munu.txt

# # # 
# python ../../../HiggsAnalysis/CombinedLimit/scripts/text2workspace.py Wenu_p.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:wModel -o Wenup.root
# combine -M MaxLikelihoodFit --robustFit 1 --rMax 3  Wenup.root --plots --saveNormalizations
# mv mlfit.root mlfitp.root
# python ../../../HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py mlfitp.root
# python diffNuisances.py -A -a  -f text mlfitp.root > mlfitp.txt
# # 
# python ../../../HiggsAnalysis/CombinedLimit/scripts/text2workspace.py Wenu_m.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:wModel -o Wenum.root
# combine -M MaxLikelihoodFit --robustFit 1 --rMax 2  Wenum.root --plots --saveNormalizations
# mv mlfit.root mlfitm.root
# python ../../../HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py mlfitm.root
# python diffNuisances.py -A -a  -f text mlfitm.root > mlfitm.txt

python ../../../HiggsAnalysis/CombinedLimit/scripts/text2workspace.py Wmunu_p.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:wModel -o Wmunup.root
combine -M MaxLikelihoodFit --robustFit 1 --rMax 3  Wmunup.root --plots --saveNormalizations
mv mlfit.root mlfitp.root
python ../../../HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py mlfitp.root
python diffNuisances.py -A -a  -f text mlfitp.root > mlfitp.txt
# 
# python ../../../HiggsAnalysis/CombinedLimit/scripts/text2workspace.py Wmunu_m.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:wModel -o Wmunum.root
# combine -M MaxLikelihoodFit --robustFit 1 --rMax 2  Wmunum.root --plots --saveNormalizations
# mv mlfit.root mlfitm.root
# python ../../../HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py mlfitm.root
# python diffNuisances.py -A -a  -f text mlfitm.root > mlfitm.txt

#combine -M MaxLikelihoodFit --preFitValue=1. --X-rtd FITTER_NEW_CROSSING_ALGO --minimizerAlgoForMinos=Minuit2 --minimizerToleranceForMinos=0.1 --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --minimizerAlgo=Minuit2 --minimizerStrategy=0 --minimizerTolerance=0.1 --cminFallbackAlgo \"Minuit2,0:1.\" --robustFit 1 Wenup.root --plots --saveNormalizations

#python ../../HiggsAnalysis/CombinedLimit/scripts/text2workspace.py Wmunu_p.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:wModel -o Wmunu.root
#combine -M MaxLikelihoodFit --robustFit 1 --rMax 1 Wmunu.root --plots --saveNormalizations
#mv mlfit.root mlfit_munu.root
#python ../../HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py mlfit_munu.root
#python diffNuisances.py -A -a  -f text mlfit_munu.root > mlfit_munu.txt

