#! /bin/bash


# SFX=_v0_SMPV
# SFX=_v1_MuSIT
# SFX=_v2_MuSta_Direct
# SFX=_v2_MuSta_histMean
# SFX=_v2_EleGSF
# SFX=_v0_SITabs

SFX=_v0
# OUTPUTDIR5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2_Acceptance
# INPUTDIR5=/eos/cms/store/user/sabrandt/StandardModel/LowPU_5TeV_Try2/
OUTPUTDIR13=.
# EFFDIR="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/testReweights_v2_2/results"
EFFDIR="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/LowPU2017ID_13TeV/results"
EFFSYSDIR="../Efficiency/Systematic/${SFX}"
#
# W->munu
# #
# root -l -q computeAccGenWm_Sys.C+\(\"wm_13.conf\",\"TEST_Wmp_file0_Gen13TeV${SFX}\",0,1\)
root -l -q computeAccGenWm_Sys.C+\(\"wm_13.conf\",\"TEST_WmpGen13TeV_allFile${SFX}\",0,1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wm_13.conf\",\"TEST_WmmGen13TeV_allFile${SFX}\",0,-1\)

# root -l -q computeAccGenWm_Sys.C+\(\"wm_13.conf\",\"WmGen13TeV${SFX}\",0,0\)
# root -l -q computeAccGenWm_Sys.C+\(\"wm_13.conf\",\"WmpGen13TeV_dressed${SFX}\",1,1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wm_13.conf\",\"WmmGen13TeV_dressed${SFX}\",1,-1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wm_13.conf\",\"WmGen13TeV_dressed${SFX}\",1,0\)
#root -l -q computeAccGenWm.C+\(\"wmm.conf\",\"Wmunu/minus\",-1\)
#root -l -q computeAccSelWm.C+\(\"wmp.conf\",\"Wmunu/plus\",1\)
#root -l -q computeAccSelWm.C+\(\"wmm.conf\",\"Wmunu/minus\",-1\)

#

# W->enu
#
# root -l -q computeAccGenWe_Sys.C+\(\"we_13.conf\",\"TEST_WepGen13TeV_allFile${SFX}\",0,1\)
# root -l -q computeAccGenWe_Sys.C+\(\"we0_13.conf\",\"TEST_WemGen13TeV_allFile${SFX}\",0,-1\)

# root -l -q computeAccGenWe_Sys.C+\(\"we1_13.conf\",\"TEST_Wep_file1_Gen13TeV${SFX}\",0,1\)
# root -l -q computeAccGenWe_Sys.C+\(\"we2_13.conf\",\"TEST_Wep_file2_Gen13TeV${SFX}\",0,1\)
# root -l -q computeAccGenWe_Sys.C+\(\"we_13.conf\",\"WemGen13TeV${SFX}\",0,-1\)
# root -l -q computeAccGenWe_Sys.C+\(\"we_13.conf\",\"WeGen13TeV${SFX}\",0,0\)
# root -l -q computeAccGenWe_Sys.C+\(\"we_13.conf\",\"WepGen13TeV_dressed${SFX}\",1,1\)
# root -l -q computeAccGenWe_Sys.C+\(\"we_13.conf\",\"WemGen13TeV_dressed${SFX}\",1,-1\)
# root -l -q computeAccGenWe_Sys.C+\(\"we_13.conf\",\"WeGen13TeV_dressed${SFX}\",1,0\)
# #root -l -q computeAccGenWe.C+\(\"wep.conf\",\"Wenu/plus\",1\)
#root -l -q computeAccGenWe.C+\(\"wem.conf\",\"Wenu/minus\",-1\)
#root -l -q computeAccSCWe.C+\(\"wep.conf\",\"Wenu/plus\",1\)
#root -l -q computeAccSCWe.C+\(\"wem.conf\",\"Wenu/minus\",-1\)
#root -l -q computeAccSelWe.C+\(\"wep.conf\",\"Wenu/plus\",1\)
#root -l -q computeAccSelWe.C+\(\"wem.conf\",\"Wenu/minus\",-1\)

#
# Z->mumu
#
# root -l -q computeAccGenZmm.C+\(\"zmm.conf\",\"Zmumu\"\)
# root -l -q computeAccGenZmm_Sys.C+\(\"zmm_13.conf\",\"TEST_ZmumuGen13TeV${SFX}\",1\)
# root -l -q computeAccGenZmm_Sys.C+\(\"zmm_5.conf\",\"ZmumuGen5TeV_dressed${SFX}\",1\) 
#root -l -q computeAccSelZmm.C+\(\"zmm.conf\",\"Zmumu\"\)
#root -l -q computeAccSelZmmBinned.C+\(\"zmm.conf\",\"Zmumu\"\)


# root -l -q computeAccSelWm_Sys.C+\(\"wm0_13.conf\",\"${EFFDIR}/Zmm/\",\"TEST_Wmfile0_pos_test${SFX}\",1,0,\"${EFFSYSDIR}/SysUnc_MuSITEff.root\",\"${EFFSYSDIR}/SysUnc_MuStaEff.root\",1\)
# root -l -q computeAccSelWm_Sys.C+\(\"wm1_13.conf\",\"${EFFDIR}/Zmm/\",\"TEST_Wmfile1_pos_test${SFX}\",1,0,\"${EFFSYSDIR}/SysUnc_MuSITEff.root\",\"${EFFSYSDIR}/SysUnc_MuStaEff.root\",1\)
# root -l -q computeAccSelWm_Sys.C+\(\"wm2_13.conf\",\"${EFFDIR}/Zmm/\",\"TEST_Wmfile2_pos_test${SFX}\",1,0,\"${EFFSYSDIR}/SysUnc_MuSITEff.root\",\"${EFFSYSDIR}/SysUnc_MuStaEff.root\",1\)

# root -l -q computeAccSelZmmBinned_Sys.C+\(\"zmm_13.conf\",\"${EFFDIR}/Zmm/\",\"TEST_Zmm${SFX}_30GeV\",0,\"${EFFSYSDIR}/SysUnc_MuSITEff.root\",\"${EFFSYSDIR}/SysUnc_MuStaEff.root\",1\)
# root -l -q computeAccSelZmmBinned_Sys.C+\(\"zmm_13.conf\",\"${EFFDIR}/Zmm/\",\"TEST_Zmm_checkYesMassCutPassing_ignoreSyst_${SFX}\",0,\"${EFFSYSDIR}/SysUnc_MuSITEff.root\",\"${EFFSYSDIR}/SysUnc_MuStaEff.root\",1\)

## We Acceptance
# root -l -q computeAccSelWe_Sys.C+\(\"we0_13.conf\",\"${EFFDIR}/Zee/\",\"TEST_Wefile0_pos_test${SFX}\",1,0,1,0,\"${EFFSYSDIR}/SysUnc_EleGSFSelEff.root\",1\)
# root -l -q computeAccSelWe_Sys.C+\(\"we1_13.conf\",\"${EFFDIR}/Zee/\",\"TEST_Wefile1_pos_test${SFX}\",1,0,1,0,\"${EFFSYSDIR}/SysUnc_EleGSFSelEff.root\",1\)
# root -l -q computeAccSelWe_Sys.C+\(\"we2_13.conf\",\"${EFFDIR}/Zee/\",\"TEST_Wefile2_pos_test${SFX}\",1,0,1,0,\"${EFFSYSDIR}/SysUnc_EleGSFSelEff.root\",1\)

# root -l -q computeAccSelZeeBinned_Sys.C+\(\"zee_13.conf\",\"${EFFDIR}/Zee/\",\"TEST_Zee_checkYesMassCut_inclGap_ignoreSyst_${SFX}\",0,1,0,\"${EFFSYSDIR}/SysUnc_EleGSFSelEff.root\",1\)
# root -l -q computeAccSelZeeBinned_Sys.C+\(\"zee_13.conf\",\"${EFFDIR}/Zee/\",\"TEST_Zee_GSFPOWPhot_ignoreSyst${SFX}\",0,1,0,\"${EFFSYSDIR}/SysUnc_EleGSFSelEff.root\",1\)
#
# Z->ee
#
# root -l -q computeAccGenZee_Sys.C+\(\"zee_13.conf\",\"ZeeGen13TeV_dressed_v1\",1\)
# root -l -q computeAccGenZee_Sys.C+\(\"zee_5.conf\",\"ZeeGen5TeV_dressed_v1\",1\)
#root -l -q computeAccSCZee.C+\(\"zee.conf\",\"Zee\"\)
# root -l -q computeAccSelZee.C+\(\"zee.conf\",\"${OUTPUTDIR}/Zee\"\)
# root -l -q computeAccSelZeeBinned_Charge.C+\(\"zee.conf\",\"blah\",\"${INPUTDIR5}\",\"${OUTPUTDIR5}/Zee\",25,2.5,0,1,0,0,0\)
# root -l -q computeAccSelZeeBinned_Charge.C+\(\"zee_13.conf\",\"blah\",\"${INPUTDIR13}\",\"${OUTPUTDIR13}/Zee\",25,2.5,0,1,0,0,1\)

# rm *.so *.d
