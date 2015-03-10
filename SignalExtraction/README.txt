MitEwk/SignalExtraction/README.txt - Jay Lawhorn 8/8/13

------| INTRODUCTION |------

This module takes selection ntuples and performs a signal extraction to get yield values for the various decay products. 

* For Z channels, the yield is calculated by simply counting the number of events that passed the selection criteria.

* For Wenu, the yield is calculated from a fit to the recoil-corrected signal MC template, an analytic QCD background model, and other background
  templates from MC with a fixed cross section ratio to our signal.

* For Wmunu, the yield is calculated from a fit to the recoil-corrected signal MC template, a QCD background model generated from an anti-isolation
  selection sample, and other background templates from MC with a fixed cross section ratio to our signal.

* The macros used in the previous analysis:

    fitWe.C
    fitWm.C
    plotZee.C
    plotZmm.C

------| RUN |------

* cmsenv

* source /scratch/ksung/ROOT/root/thisroot.sh (latest version of ROOFit again)

* source runFits.sh

------| MODIFICATION |------

* in runFits.sh:

    integrated luminosity is defined here

* in (fit/plot)(We/Wm/Zee/Zmm).C:

    scale and resolution corrections (getScaleCorr and getResCorr)

    pt and eta cut thresholds

    recoil correction files

    input ntuple locations