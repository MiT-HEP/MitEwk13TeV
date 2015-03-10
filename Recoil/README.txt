MitEwk/Recoil/README.txt - Jay Lawhorn 8/14/13

------| INTRODUCTION |------

* This module generates the recoil correction templates used in the SignalExtraction module.

* While there are a number of macros available here for testing/validation purposes, the main macros of interest are
  fitRecoilZee.C and fitRecoilZmm.C

* The input syntax for each macro basically identical, and follows this pattern:

    root -l -q fitRecoilZmm.C+\(\"input_ntuple_dir\",u1_model,u2_model,sig_only,[charge,]\"output_dir\"\)

  u1_model and u2_model take values {1, 2, 3} and correspond to how many gaussians are used in the fit

  sig_only is a boolean flag where 0 = use backgrounds and 1 = don't use backgrounds

  charge is only specified for W channels. takes values {-1, 0, 1} for a specific charge or 0 for both charges

------| RUN |------

* After running cmsenv, run:

    source /scratch/ksung/ROOT/root/bin/thisroot.sh

  to use a version of ROOT with the most up to date (as of summer 13) version of ROOFit, which will be used by this
  module.

* Then run

    source runSelection.sh
void fitRecoilZee(TString infoldername,  // input ntuple folder                                                                                     
                  Int_t   pfu1model,     // u1 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)                             
                  Int_t   pfu2model,     // u2 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)                             
                  Bool_t  sigOnly,       // signal event only?                                                                                      
                  TString outputDir      // output directory
  to execute the needed recoil corrections macros.

------| MODIFICATIONS |------

* You will want to use recoil corrections obtained from Zll data. However, it's good to make sanity checks on your corrections
  by making them with Zll MC as well as Wlnu MC. The corrections obtained from Wlnu MC's should make essentially no change
  to the Wlnu MC distribution. The Zll MC corrections should in principle make next to no changes to the Wlnu MC distribution
  as well. 

* You will need to change the input parameters to each macro in runRecoilFits.sh. The input syntax is given above.

* The data and correction pT binnings for each macro are given near line 145 and will likely need some tinkering.

* The input file strings are given directly below the binnings, followed by the kinematic cuts for our signal selections. These
  will also possibly need changes.

* Starting parameters and parameter limits for the functions being fit to the recoil parameters in boson pT can be set similarly to
  lines 210-228 in fitRecoilZee.C. 