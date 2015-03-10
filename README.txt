MitEwk/README.txt - Jay Lawhorn 8/28/13

------| INTRODUCTION |------

This folder contains the various modules needed to go from BAMBU datasets to W/Z production cross section measurements. 

This code base has a twiki at: http://www.cmsaf.mit.edu/twiki/bin/view/CmsHep/WZBosonCrossSection and each module has its
own README.txt included. The general intention is that the TWIKI page should be read first for an explanation of the 
experimental/scientific background the overall analysis and each individual module. The README.txt's explain the technical
minutia of how to actually run the code as well as guidance on where modifications will be needed to use this code on new
data. 

------| MODULES |------

The order of modules is roughly as follows:

Ntupler:	     BAMBU		-> bacon ntuple

Selection:	     bacon ntuple	-> flat ntuple
		     (e scale cor)

EleScale:	     zee ntuple		-> electron scale/res corrections
		     raw sel ntuple

Efficiency:	     flat ntuple	-> lepton efficiency root files

Acceptance:	     flat ntuple	-> acceptance numbers
		     lep. eff. files

Recoil:		     flat ntuple	-> recoil correction root files

Signal Extraction:   flat ntuple	-> signal yield numbers
		     rec. corr. files

Summary Plots:	     signal yields	-> cross sections and plots
		     acc. numbers
		     syst. uncer. %'s

------| TOOLS |------

The Tools folder contains a number of useful macros I needed along the way including ones for merging root files, merging and doing simple logic with json files, etc. I also used tools available at https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiCalc from the official CMS Luminosity calculation tools.

-----| SETUP INFO |------

Commands to execute on start up for each terminal (or add to your ~/.bashrc):

export SCRAM_ARCH=slc5_amd64_gcc462
export VO_CMS_SW_DIR=/server/02a/cmsprod/cmssoft
source /home/cmsprod/cmssoft/cmsset_default.sh
cd (this directory)
cmsenv