MitEwk/Selection/README.txt - Jay Lawhorn 8/6/13

------| INTRODUCTION |------

* This module takes bacon ntuples from the Ntupler folder and outputs selection ntuples (flat trees). 
  Currently runs on the scntuples in /scratch/ksung/EWKAna/bacon. 

* There are 6 different C macros that make selection cuts for the 4 different signal channels
  (Wenu, Wmunu, Zee, Zmumu) and 2 anti-isolation channels (AntiWenu, AntiWmunu) used for QCD background estimation/validation

* The Zee and Zmumu selection macros also identify tag and probe leptons and classify events based on how well the probe lepton is
  reconstructed for use in efficiency calculations

* Selection macros are select().C

* ().conf specify the input ntuples for data, signal and bacground MCs 

* ConfParse.hh parses the .conf files

* plot().C generate plots for the corresponding channel just to get a basic look at the samples


------| RUN |------

* To run on Kevin's Ntuples with no modification, can just do

source runSelection.sh

------| MODIFICATIONS |------

* The luminosity section JSON files used here should contain only lumisecs contained in the data that are ALSO certified lumisecs. 
  (For further info on how to do this see the corresponding twiki page: http://www.cmsaf.mit.edu/twiki/bin/view/CmsHep/WZBosonCrossSection)

* MitEwk/Utils/LeptonIDCuts.hh:

    contains lepton identification requirements if anything of that nature needs updated

* runSelection.sh:

    NTUPDIR: input ntuple directory

    LUMI: integrated luminosity

* (zee/zmm/we/wm).conf:

    add new ntuples for each category

    lines starting with a "$" specify sample type and plotting options with following format:

        $ sample_type fill_color line_color legend_label

    lines without a "$" or "#" specify input file locations and info with following format:

        sample_location mc_xsec data_json_file

	for data, mc_xsec = 0	for mc, data_json_file = NONE

    blank lines will break the parser, so be viligant 

* select().C:

    kinematic cuts are given near top of file (NOTE: these aren't the final kinematics cuts; they're a bit less stringent)

    electron energy scale corrections in same place

    ULong64_t trigger and ULong64_t triggerObj should be updated with your new triggers from Ntupler/interface/EWKAnaDefs.hh