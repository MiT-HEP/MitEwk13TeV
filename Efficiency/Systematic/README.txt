EWKAna/Efficiency/Systematic/README.txt - Jay Lawhorn 8/27/13

Code for calculating lepton efficiency systematic uncertainties. 

0) Use ROOT version from /scratch/ksung/ROOT/root/bin/thisroot.(c)sh. This has newest version of RooFit.

1) Make generator templates. 
   Run makePsExpTemplates.C (examples in getPsExpTemplates.sh) using ALTERNATE signal/background model in question. 
   (So if you'd usually use MC templates for signal and an exponential for background, you'd use a Crystal Ball signal template
   and exponential background OR MC signal template and linear background, etc...) This step saves a ROOT file containing a RooFit 
   workspace with generator fit info as well as binning information and number of events in that bin.

2) Generate pseudo experiments.
   Run makePseudoData.C (examples in getPsExpTemplates.sh). You'll need to generate a good number of pseudo experiments per bin (~1000).
   This step reads in the ROOT file from step one, and generates the desired number of pseudo experiments. It needs to be run once per bin.
   NOTE: This step doesn't do random seeding properly, so you'll either need to generate all the psuedo experiments in one go or take care of
   that to avoid getting the same n events over and over again.

3) Fit pseudo experiements.
   Run submitPseudoExperiments.sh ....insert more detail here...

4) Do math!