MitEwk13TeV/Selection/README.md - Stephanie Brandt .... February 10, 2020

# Overview
Skim a flat ntuple from the Bacon files. Current version compatible with CMSSW_9_4_X.

## Selection steps
There are six macros that select events for each of the channels: select[We,Wm,Zee,Zmm].C for the signals and selectAnti[We,Wm].C for the control regions used for the W fits. 

### Processing order
The selection code can be run in any order. You should make sure you have the relevant flat ntuples before proceeding to the next steps:
1. Efficiency (Z only, electron and muon channels)
1. Recoil (W and Z, muon channel only)
1. NtupleMod (W only, requires W ntuple + Efficiency + Recoil as inputs)
1. SignalExtraction (all)

# Code
## Z selection
The Z selection can be done using the macros:
* selectZee.C
* selectZmm.C
These contain a loose selection which allows the ntuples to be used for the tag-and-probe for the lepton efficiency calculation

## W selection
W signal events can be selected using: 
* selectWe.C
* selectWm.C
Selection of a QCD-enriched control region is done using: 
* selectAntiWe.C
* selectAntiWm.C
The current control region selection uses a reversed isolation cut, based on the standard lepton cut-based ID requirements. 

## Configuration Files
** See next section for important info**
The * *.conf * files specify the locations of the Bacon, type of events contained in the Bacon, and the cross-section of the samples (or the JSON for data). These are parsed by ConfParse.hh. JSONs should be located in this directory as well. 

plot*.C can generate plots, but these are not up-to-date with newest versions of code. 

# Using the code
Check one of the scripts for the macro running setup (runSelection.sh or runSelectionZ*.sh)
I have added some extra arguments which allow the ntuple to be split into sections during the selection and allows for the parallelization of the selection and full recoil+efficiency step for the W ntuples. The full setup is included in the runSelection.sh, the **NSEC** and **ITH** variables control the number of sections and which section is being processed. This causes them to be put into different directories with numbers specifying the splitting i.e. ntuple_0_1. 

## Helper code
### Trigger and Lepton ID information
Cut-based ID for the muons and electrons is stored in the LeptonIDCuts.hh file. Cuts are currently the 2017 standard, and the triggers are for the 2017G and 2017H runs. 
MitEwk13TeV/Utils/LeptonIDCuts.hh

### conf files  (zee/zmm/we/wm).conf:

    add new ntuples for each category

    lines starting with a "$" specify sample type and plotting options with following format:

        $ sample_type fill_color line_color legend_label

    lines without a "$" or "#" specify input file locations and info with following format:

        sample_location mc_xsec data_json_file

	for data, mc_xsec = 0	for mc, data_json_file = NONE

    blank lines will break the parser, so be viligant 