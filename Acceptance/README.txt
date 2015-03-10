MitEwk/Acceptance/README.txt - Jay Lawhorn 8/8/13

------| INTRODUCTION |------

This module computes acceptance and efficiency-corrected acceptance values for the various signal channels. The acceptance is calculated from 
MC only and takes into account the geometry of the detector and our kinematic cuts (lepton pT). Efficiency-corrected acceptances also take into account
the various efficiencies calculated in EWKAna/Efficiency and the corresponding scale factors between MC and data.

* runAcc.sh is the main shell wrapper for the acceptance computations. 

* The previous analysis used the following scripts for the final acceptance calculations:

    computeAccSelWe.C
    computeAccSelWm.C
    computeAccSelZeeBinned.C
    computeAccSelZmmBinned.C

* The other scripts in this folder compute intermediate acceptances instead of the full selection acceptance.

------| RUN |------

* After setting up cmsenv, run:

    source /scratch/ksung/ROOT/root/thisroot.sh

  which setup a version of ROOT with the most up to date version (as of Summer 2013) of ROOFit which is needed for this module.

* runAcc.sh will execute the various macros for you.

* Each macro takes input in the following pattern:

    root -l -q computeAccSelWm.C+\(\"input.conf\",\"output_dir\"[,charge]\)

    input.conf is a text file that specifies the location of an input bacon ntuple in the following format:

      bacon_file_path histogram_color histogram_line_style @legend_label

    charge can be either -1 or 1 and is only specified for the W channel macros.

------| MODIFICATION |------

* in each of computeAccSel().C:

    kinematic cuts and efficiency scale factor files are specified near the top of each file and will need modification

* in each of ().conf:

    change file paths