#! /bin/bash

IN13=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_13TeV
# IN13=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_newXsec
# INPUTDIR5=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_5TeV

OUT=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Recoil_Orig

LUMI13=199.2
LUMI5=291.107
#need to update the lumi^
#LUMI=1 #fitting for the fractions

#--------------------------------------------------------------------------------------------------
#-            CENTRAL VALUES 
# -------------------------------------------------------------------------------------------------
# # ## 13 TeV W MC
# root -l -b -q fitRecoilWm.C+\(\"${IN13}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_2G\",${LUMI13},0,0\) 
# root -l -b -q fitRecoilWm.C+\(\"${IN13}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_2G\",${LUMI13},0,0\)  
# # 13 TeV Z MC
root -l -b -q fitRecoilZmm.C+\(\"${IN13}/Zmumu/ntuples/zmm_select.raw.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_2G\",${LUMI13},0,0\)
# # 13 TeV Z Data
# root -l -b -q fitRecoilZmm.C+\(\"${IN13}/Zmumu/ntuples/data_select.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_2G_bkg_fixRoch\",${LUMI13},0,0\)  

###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ETA-BINNED  
# -------------------------------------------------------------------------------------------------
# for i in {3..3}
# do
# # # # # # # # # 13 TeV W MC
# # root -l -b -q fitRecoilWm.C+\(\"${IN13}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0\) 
# root -l -b -q fitRecoilWm.C+\(\"${IN13}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0\) 
# # # # # # # # # # 13 TeV Z MC
# # root -l -b -q fitRecoilZmm.C+\(\"${IN13}/Zmumu/ntuples/zmm_select.raw.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0\)
# # # # # # # # # # 13 TeV Z Data
# # # root -l -b -q fitRecoilZmm.C+\(\"${IN13}/Zmumu/ntuples/data_select.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_Eta${i}\",${LUMI13},${i},0\)  
# done


###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ROOKEYS PDF
# -------------------------------------------------------------------------------------------------
## 13 TeV W MC
# root -l -b -q fitRecoilWm.C+\(\"${IN13}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_Keys\",${LUMI13},0,1\) 
# root -l -b -q fitRecoilWm.C+\(\"${IN13}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_Keys\",${LUMI13},0,1\) 
# # 13 TeV Z MC
# root -l -b -q fitRecoilZmm.C+\(\"${IN13}/Zmumu/ntuples/zmm_select.raw.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_Keys\",${LUMI13},0,1\)
# 13 TeV Z Data
# root -l -b -q fitRecoilZmm.C+\(\"${IN13}/Zmumu/ntuples/data_select.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_Keys\",${LUMI13},0,1\) 

###################################################################################################
###################################################################################################



#rm *.so *.d *.pcm
