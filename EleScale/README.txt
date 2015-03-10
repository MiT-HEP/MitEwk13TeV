MitEwk/EleScale/README.txt - Jay Lawhorn 8/13/13

------| INTRODUCTION |------

This module computes electron energy scale and resolution correction factors for electron data from raw Zee data and Zee MC. Because the Zee dilepton mass peak is well known 
and mostly background free, the energy corrections are calculated to from a comparison of this peak in data and in MC. These corrections are calculated in bins of eta because of differences in the ECAL barrel and endcap regions using simultaneous fits to dileptons in the various combinations of the two regions (EE-EE, EB-EB, EE-EB).

------| RUN |------

cmsenv
source /scratch/ksung/ROOT/root/bin/thisroot.sh
root -l -q EleScale.C+

------| MODIFICATION |------

in EleScale.C:
   * output file name
   * pileup reweighting file name
   * Data and MC file names
   * Dilepton mass window
   * pT cut
   * eta cut
