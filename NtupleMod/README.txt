MitEwk/NtupleMod/README.txt - Stephanie Brandt ... February 10, 2020

# Overview

This module takes W selection flat ntuples and adds several branches containing corrected variables (recoil corrections, efficiency scale factors, etc). 

## Why?
Some of the steps such as recoil corrections are computationally expensive, and were previously done in the SignalExtraction code. I changed the workflow so that these can be pre-computed and the computation also parallelized. This also allows us to have a single piece of code for the final W signal fit setup, 

## Output
Output is another flat ntuple file, with all original branches preserved and several new branches added. Most of the new branches are vectors containing information grouped by type, i.e. metVars contains all corrected MET values corresponding to the recoil corrections and alternative models for systematics 
### MET variables
**metVars** and **metVarsPhi**  are vectors containing MET info after performing recoil corrections and adjusting for lepton energy scale
1. no corrections
1. main correction
1. rapidity-binned correction
1. gaussian kernel fit
1. lepton energy scale uncertainty (high)
1. lepton energy scale uncertainty (low)
1. statistical uncertainty [10 of these for each of the free parameters per pT bin]

### Efficiency
**evtWeight** vector contains the product of 
* the event MC weight
* efficiency for the event
* prefire weight for the event

The elements in the vector correspond to the different systematic uncertainties.
1. central value
1. MC systematic
1. FSR systematic
1. Background model systematic
1. Tag pT cut systematic
1. Efficiency SF statistical uncertainty
1. Prefire Probability statistical uncertainty

## tl;dr
W selection + efficiency SF + recoil corrections -> **NtupleMod** -> Signal Extraction


*MET recoil corrections - does all of the correction and application and stores the corrected MET in new variables
*Lepton scale corrections (as needed) - apply the rochester corretions (electron is applied at selection step)
*calculates the per-event lepton efficiency scale factors and includes it in a new event weight variable



# Code

###Macros 
* eleNtupleMod.C
* muonNtupleMod.C
###Scripts 
* runFitsMu13.sh 
* runFitsEle13.sh
* runFitsMu5.sh
* runFitsEle5.sh 



