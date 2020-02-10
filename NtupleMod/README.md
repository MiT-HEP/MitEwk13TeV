MitEwk13TeV/NtupleMod/README.md - Stephanie Brandt ... February 10, 2020

# Overview

This module takes W selection flat ntuples and adds several branches containing corrected variables (recoil corrections, efficiency scale factors, etc). 

## Why?
Some of the steps such as recoil corrections are computationally expensive, and were previously done in the SignalExtraction code. Having to rerun this and repeatedly perform the calculations was extremely slow, and the corrections don't change often, so I changed the workflow so that these can be parallelized and pre-computed. I have also ensured that the muon and electron ntuples have consistency between variable names so that the final signal extraction fit code is no longer 2 separate macros.

## Input
* Ntuples from selectW*.C or selectAntiW*.C
* Recoil corrections
* Efficiency scale factors & uncertainties
If you don't have a full set of recoil corrections for systematics you can turn off parts of the computation by resetting some the booleans called doKeys, doStat, doEta, etc. in the macros. The default value for those will be MET without recoil corrections applied.

## Output
Output is another flat ntuple file, with all original branches preserved and several new branches added. Most of the new branches are vectors containing information grouped by type, i.e. metVars contains all corrected MET values corresponding to the recoil corrections and alternative models for systematics. The main modifications are listed below: 
### Electron and Muon ntuple uniformity
Muon and electron code differs in the stage at which energy scale / Rochester corrections are applied. Electrons have it applied at the selection step. Muons previously had it applied in the Signal Extraction fit step, but I have moved it to this code. The output ntuples are now consistent with nomenclature regarding the lepton TLorentzVectors: 
* **lep** *TLorentzVector* contains the scale-corrected lepton
* **lep_raw** *TLorentzVector* contains the raw lepton without energy/momentum scale correction

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


# Code
The code for this is fairly contained, just two macros and some scripts to run.
### Macros 
* eleNtupleMod.C
* muonNtupleMod.C
### Scripts 
* runFitsMu13.sh 
* runFitsEle13.sh
* runFitsMu5.sh
* runFitsEle5.sh 



