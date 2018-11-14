#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=2305
LUMI2=96

# root -l -q fitWlikeZe.C+\(\"WlikeEle_incl_FP_test_full_ultrafine\",${LUMI},${LUMI},0\)
 #root -l -q fitWlikeZm.C+\(\"WlikeMu_incl_FP_full_ultrafine\",${LUMI},${LUMI},0\)
 # mv ./Wenu_pdfTemplates.root Wenu_newBacon_JUL5recoil_Zee_workspace
 
# root -l -q fitWe_toys_dualinput.C+\(\"toys_Zmm_Vs_Zmm_test_a3free\",${LUMI},${LUMI},0\)

#root -l -q compareEleMuTemplate.C+\(\"2018_03_11_EleToMu\",${LUMI},${LUMI},0\)
#root -l -q compare10Segments.C+\(\"2018_03_13_compareRW0to1_ele\",${LUMI},${LUMI},0\)
#root -l -q combineWorkspaces.C+\(\"testKeys_2\",\"Wenu_lumi0_PUrwMC_dataSets_incl2/\",\"Wenu_lumi0_PUrwMC_dataSets_eta2\",\"eta\",${LUMI},${LUMI},0\)
#root -l -q combineWorkspaces.C+\(\"testKeys_2\",\"Wenu_lumi0_PUrwMC_dataSets_incl/\",\"Wenu_lumi0_PUrwMC_dataSets_eta\",\"eta\",${LUMI},${LUMI},0\)
#root -l -q combineWorkspaces.C+\(\"testKeys_rb\",\"Wenu_allLumi_WkspPdfs_incl/\",\"Wenu_allLumi_WkspPdfs_eta/\",\"eta\",${LUMI},${LUMI},0\)
#root -l -q combineWorkspaces.C+\(\"testKeys_large\",\"Wenu_allLumi_WkspPdfs_incl/\",\"Wenu_allLumi_WkspPdfs_eta/\",\"eta\",${LUMI},${LUMI},0\)

#root -l -q combineBinnedHist.C+\(\"fullStat_test660kKeyWpSig_lumi0_3\",\"Wenu_lumi0_incl_binnedPdfs_moreBin_3/\",\"Wenu_lumi0_eta_binnedPdfs_moreBin_3/\",\"eta\",${LUMI},${LUMI},0\)
#root -l -q combineBinnedHist.C+\(\"fullStat_eta_keys_Key100k\",\"fullStat_eta_Key100k/\",\"Wenu_lumi0_keys_binnedPdfs_moreBin_2/\",\"keys\",${LUMI},${LUMI},0\)

#for it in {0..0}
#do
#OUTFOLDER=shapes_pol4r2_lumi${it}_mu4
#root -l -q combineBinnedHistmu.C+\(\"${OUTFOLDER}\",\"Wmunu_lumi${it}_incl_moreBin_4/\",\"Wmunu_lumi${it}_eta_moreBin_4/\",\"eta\",${LUMI},${LUMI},0\)
#root -l -q combineBinnedHistmu.C+\(\"${OUTFOLDER}\",\"${OUTFOLDER}\",\"Wmunu_lumi${it}_keys_moreBin_4/\",\"keys\",${LUMI},${LUMI},0\)
###root -l -q fitWe_toys_dualinput.C+\(\"toysUncShapes_pol11rb1_downEWK\",${LUMI},${LUMI},0\)
#mv ${OUTFOLDER}/Wmunu_pdfTemplates_binned.root  ${OUTFOLDER}/Wmunu_pdfTemplates_binned_${it}.root
#done

isMuon=true
##root -l -q combineWorkspaces.C+\(\"shapes_pdfs_etaKeys\",\"shapes_pdfs_eta\",\"Wenu_lumi0_PUrwMC_keys_pdfs\",\"keys\",${LUMI},${LUMI},0\)
for it in {0..9}
do
# ################ MUON SHIT #######################
# OUTFOLDER=Wmunu_met_lumi${it}_35GeV
# # OUTFOLDER2=Wmunu_met_lumi${it}_qcdnorm_shapeUnc_all_stat
# root -l -q combineBinnedHistMu.C+\(\"${OUTFOLDER}\",\"Wmunu_lumi${it}_met_pt35GeV_central/\",\"Wmunu_lumi${it}_met_pt35GeV_eta/\",\"eta\",${LUMI},${LUMI},0\)
# # root -l -q combineBinnedHistMu.C+\(\"${OUTFOLDER}\",\"${OUTFOLDER}\",\"Wmunu_lumi${it}_met_keys_fixZxx/\",\"keys\",${LUMI},${LUMI},0\)
# # root -l -q combineBinnedHistMu.C+\(\"${OUTFOLDER}\",\"${OUTFOLDER}\",\"Wmunu_lumi${it}_met_up1_fixZxx/\",\"stat\",${LUMI},${LUMI},0\)
# mv ${OUTFOLDER}/Wmunu_pdfTemplates_binned.root  ${OUTFOLDER}/Wmunu_pdfTemplates_binned_${it}.root
# mv ${OUTFOLDER}/Wmunu_pdfTemplates_binned_${it}.root  ${OUTFOLDER}/Wmunu_pdfTemplates_binned_stat_${it}.root

 # root -l -q add_stat_shapes.C+\(\"${OUTFOLDER}\",\"${OUTFOLDER2}\",\"hWxMet\",\"wxMet\",${it},${isMuon},0\)
 # root -l -q add_stat_shapes.C+\(\"${OUTFOLDER2}\",\"${OUTFOLDER2}\",\"hZxxMet\",\"zxxMet\",${it},${isMuon},0\)
 # root -l -q add_stat_shapes.C+\(\"${OUTFOLDER2}\",\"${OUTFOLDER2}\",\"hEwkMet\",\"ewkMet\",${it},${isMuon},0\)
 # root -l -q add_stat_shapes.C+\(\"${OUTFOLDER2}\",\"${OUTFOLDER2}\",\"hDibMet\",\"dibMet\",${it},${isMuon},0\)
 # root -l -q add_stat_shapes.C+\(\"${OUTFOLDER2}\",\"${OUTFOLDER2}\",\"hTtbMet\",\"ttbMet\",${it},${isMuon},0\)
# mv ${OUTFOLDER2}/Wmunu_pdfTemplates_binned_${it}.root  ${OUTFOLDER2}/Wmunu_pdfTemplates_binned_stat_${it}.root

# cp ${OUTFOLDER}/Wmunu_pdfTemplates_binned_stat_${it}.root /home/sabrandt/SM/InclusiveMaster/Combine/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/datacards/2018_10_12_Muons_35GeV
#cp ${OUTFOLDER2}/Wmunu_pdfTemplates_binned_stat_${it}.root /home/sabrandt/SM/InclusiveMaster/Combine/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/datacards/2018_10_02_ExtractMu_ShapesLowerPol
# mv 
########## ELECTRON SHIT ####################
 # OUTFOLDER=Wenu_met_lumi${it}_35GeV
 OUTFOLDER=Zmm_lumi${it}_Wlike
# # OUTFOLDER2=Wenu_met_lumi${it}_shapeUnc_all_stat
 # # root -l -q combineNoBinsMuons.C+\(\"${OUTFOLDER}\",\"Zmm_lumi${it}_Wlike_incl\",\"Zmm_lumi${it}_Wlike_eta/\",\"eta\",${LUMI},${LUMI},0,0\)
 # root -l -q combineNoBinsMuons.C+\(\"${OUTFOLDER}\",\"${OUTFOLDER}\",\"Zmm_lumi${it}_Wlike_keys/\",\"keys\",${LUMI},${LUMI},0,0\)
 # root -l -q combineNoBinsMuons.C+\(\"${OUTFOLDER}\",\"${OUTFOLDER}\",\"Zmm_lumi${it}_Wlike_diag/\",\"diag\",${LUMI},${LUMI},0,0\)
 # root -l -q combineNoBins.C+\(\"${OUTFOLDER}\",\"Zee_lumi${it}_Wlike_incl\",\"Zee_lumi${it}_Wlike_eta/\",\"eta\",${LUMI},${LUMI},0,0\)
 # root -l -q combineNoBins.C+\(\"${OUTFOLDER}\",\"${OUTFOLDER}\",\"Zee_lumi${it}_Wlike_keys/\",\"keys\",${LUMI},${LUMI},0,0\)
 # root -l -q combineNoBins.C+\(\"${OUTFOLDER}\",\"${OUTFOLDER}\",\"Zee_lumi${it}_Wlike_diag/\",\"diag\",${LUMI},${LUMI},0,0\)
 # # root -l -q combineBinnedHist.C+\(\"${OUTFOLDER}\",\"Wenu_lumi${it}_met_pt35GeV_central\",\"Wenu_lumi${it}_met_pt35GeV_eta/\",\"eta\",${LUMI},${LUMI},0,0\)
 # # root -l -q combineBinnedHist.C+\(\"${OUTFOLDER}\",\"${OUTFOLDER}\",\"Wenu_lumi${it}_met_keys_fixZxx/\",\"keys\",${LUMI},${LUMI},0,0\)
 # # root -l -q combineBinnedHist.C+\(\"${OUTFOLDER}\",\"${OUTFOLDER}\",\"Wenu_lumi${it}_met_up1_fixZxx/\",\"stat\",${LUMI},${LUMI},0,0\)
 # # root -l -q add_stat_shapes.C+\(\"${OUTFOLDER}\",\"${OUTFOLDER2}\",\"hWxMet\",\"wxMet\",${it},${isMuon},0\)
 # # # root -l -q add_stat_shapes.C+\(\"${OUTFOLDER2}\",\"${OUTFOLDER2}\",\"hZxxMet\",\"zxxMet\",${it},${isMuon},0\)
 # # # root -l -q add_stat_shapes.C+\(\"${OUTFOLDER2}\",\"${OUTFOLDER2}\",\"hEwkMet\",\"ewkMet\",${it},${isMuon},0\)
 # # # root -l -q add_stat_shapes.C+\(\"${OUTFOLDER2}\",\"${OUTFOLDER2}\",\"hDibMet\",\"dibMet\",${it},${isMuon},0\)
 # # # root -l -q add_stat_shapes.C+\(\"${OUTFOLDER2}\",\"${OUTFOLDER2}\",\"hTtbMet\",\"ttbMet\",${it},${isMuon},0\)
# mv ${OUTFOLDER}/Wenu_pdfTemplates.root  ${OUTFOLDER}/Wenu_pdfTemplates_${it}.root
mv ${OUTFOLDER}/Wmunu_pdfTemplates.root  ${OUTFOLDER}/Wmunu_pdfTemplates_${it}.root
# mv ${OUTFOLDER}/Wenu_pdfTemplates_binned.root  ${OUTFOLDER}/Wenu_pdfTemplates_binned_${it}.root
# mv ${OUTFOLDER}/Wenu_pdfTemplates_binned_${it}.root  ${OUTFOLDER}/Wenu_pdfTemplates_binned_stat_${it}.root
# # mv ${OUTFOLDER2}/Wenu_pdfTemplates_binned_${it}.root  ${OUTFOLDER2}/Wenu_pdfTemplates_binned_stat_${it}.root
# # cp ${OUTFOLDER2}/Wenu_pdfTemplates_binned_stat_${it}.root /home/sabrandt/SM/InclusiveMaster/Combine/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/datacards/2018_09_24_ExtractEle
cp ${OUTFOLDER}/Wmunu_pdfTemplates_${it}.root /home/sabrandt/SM/InclusiveMaster/Combine/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/datacards/2018_10_29_Z_Wlike
done

#root -l -q combineWorkspaces.C+\(\"shapes_pdfs_etaKeys\",\"shapes_pdfs_eta\",\"Wenu_lumi0_PUrwMC_keys_pdfs\",\"keys\",${LUMI},${LUMI},0\)
#root -l -q combineWorkspaces.C+\(\"data0incl_shapeUnc1\",\"Wenu_lumi0_PUrwMC_fixWksp_incl_2/\",\"Wenu_lumi1_PUrwMC_fixWksp_incl_lumiMC\",\"Pileup\",${LUMI},${LUMI},0\)
#root -l -q combineWorkspaces.C+\(\"toy1fit0_shapeUnc\",\"toy1fit0_shapeUnc\",\"Wenu_lumi0_PUrwMC_fixWksp_keys\",\"keys\",${LUMI},${LUMI},0\)
#root -l -q compare10SegmentsMuons.C+\(\"2018_03_13_compareRW0to1_Muons\",${LUMI},${LUMI},0\)

# root -l -q drawTemplates.C+\(\"Wenu_toys_2_QCDPdf_Bin1.0GeV_pdfs\",${LUMI},${LUMI},0\)
#rm *.so *.d
