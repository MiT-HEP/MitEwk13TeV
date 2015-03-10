#!/bin/bash

folder=$1

count=0

for file in `ls ${folder}/analysis/plots/fitresetapt_0.txt`
do

  echo $file $count >> region_map.txt

  massFail=`grep massFail $file | awk '{print $2}'`
  massPass=`grep massPass $file | awk '{print $2}'`
  widthFail=`grep widthFail $file | awk '{print $2}'`
  widthPass=`grep widthPass $file | awk '{print $2}'`
  NbkgFail=`grep NbkgFail $file | awk '{print $3}'`
  NbkgPass=`grep NbkgPass $file | awk '{print $3}'`
  Nsig=`grep Nsig $file | awk '{print $3}'`
  alphaFail=`grep alphaFail $file | awk '{print $3}'`
  alphaPass=`grep alphaPass $file | awk '{print $3}'`
  eff=`grep eff $file | awk '{print $3}'`
  meanFail=`grep meanFail $file | awk '{print $3}'`
  meanPass=`grep meanPass $file | awk '{print $3}'`
  nFail=`grep [[:blank:]]nFail $file | awk '{print $3}'`
  nPass=`grep [[:blank:]]nPass $file | awk '{print $3}'`
  sigmaFail=`grep sigmaFail $file | awk '{print $3}'`
  sigmaPass=`grep sigmaPass $file | awk '{print $3}'`
  tFail=`grep tFail $file | awk '{print $3}'`
  tPass=`grep tPass $file | awk '{print $3}'`
  Passing=`grep Passing $file | awk '{print $2}'`
  Failing=`grep Failing $file | awk '{print $2}'`
  etaLo=`grep xbinLo $file | awk '{print $2}'`
  etaHi=`grep xbinHi $file | awk '{print $2}'`
  ptLo=`grep ybinLo $file | awk '{print $2}'`
  ptHi=`grep ybinHi $file | awk '{print $2}'`
  if grep -q signEtaPT $file
      then
      doAbsEta=0
      else
      doAbsEta=1
  fi

 echo  ./runJob.sh $count $massFail $massPass $widthFail $widthPass $NbkgFail $NbkgPass \
     $Nsig $alphaFail $alphaPass $eff $meanFail $meanPass $nFail $nPass $sigmaFail $sigmaPass $tFail $tPass $Passing $Failing $etaLo $etaHi $ptLo $ptHi $doAbsEta

 ./runJob.sh $count $massFail $massPass $widthFail $widthPass $NbkgFail $NbkgPass $Nsig $alphaFail $alphaPass $eff $meanFail $meanPass $nFail $nPass $sigmaFail $sigmaPass $tFail $tPass $Passing $Failing $etaLo $etaHi $ptLo $ptHi $doAbsEta

  count=$count+1

done