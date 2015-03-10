#!/bin/bash
effType=$1

for x in 9
do
  echo $x

  #ls /scratch/klawhorn/EffSysStore/CB_MuSelEffResults/etapt_$x_*output 
  cat /scratch/klawhorn/EffSysStore/${effType}Results/etapt_${x}_*output | grep eff | awk '{ print $3 }' > ${effType}_BinUncert/etapt_${x}.dat
  #grep "Full" /scratch/klawhorn/EffSysStore/${effType}Results/etapt_${x}_*output |wc

done

#