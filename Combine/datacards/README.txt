These datacards can be used to find the signal and background yields for all channels.

First, copy the PhysicsModel.py file in this folder into your CombinedLimit/python directory. This will give you the option of constraining the signal and EWK yields for the W channels to the ratio of the corresponding theoretical cross sections, as was done in the 8 TeV analysis.

For example, for the W->enu channel, run these commands from your CombinedLimit directory:

>> python ./scripts/text2workspace.py Wenu.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:wModel -o Wenu.root // Should create a file name Wenu.root that can be run through Combine.

>> combine -M MaxLikelihoodFit Wenu.root --plots --saveNormalizations // The --plots option saves the bkg-only and sig+bkg fits to the file mlfit.root. The --saveNormalizations option allows you to run the next command.

>> python ./test/mlfitNormsToText.py mlfit.root // Should output a table of yields where the first column is the process name, the second column is the yield in the sig+bkg fit, and the third column is the yield in the bkg-only fit.

If you don't want to constrain the signal and EWK yields, take out the "-P HiggsAnalysis.CombinedLimit.PhysicsModel:wModel" in the first command.
