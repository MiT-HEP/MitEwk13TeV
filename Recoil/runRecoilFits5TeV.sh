#! /bin/bash

# IN5=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_5TeV
# IN5=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_5TeV_newXsec
IN5=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_5TeV

OUT=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Recoil

# LUMI5=199.2
LUMI5=291.107
#need to update the lumi^
#LUMI=1 #fitting for the fractions

#--------------------------------------------------------------------------------------------------
#-            CENTRAL VALUES 
# -------------------------------------------------------------------------------------------------
# # ## 5 TeV W MC
# root -l -b -q fitRecoilWm.C+\(\"${IN5}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_5TeV_2G\",${LUMI5},0,0\) 
# root -l -b -q fitRecoilWm.C+\(\"${IN5}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_5TeV_2G\",${LUMI5},0,0\)  
# # 5 TeV Z MC
# root -l -b -q fitRecoilZmm.C+\(\"${IN5}/Zmumu/ntuples/zmm_select.raw.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_5TeV_2G\",${LUMI5},0,0\)
# # 5 TeV Z Data
# root -l -b -q fitRecoilZmm.C+\(\"${IN5}/Zmumu/ntuples/data_select.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_5TeV_2G\",${LUMI5},0,0\)  
# root -l -b -q fitRecoilZmm.C+\(\"${IN5}/Zmumu/ntuples/data_select.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_5TeV_2G_bkg_fixRoch\",${LUMI5},0,0\)  

###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ETA-BINNED  
# -------------------------------------------------------------------------------------------------
# for i in {1..3}
# do
# # # # # # # # # # 5 TeV W MC
# # root -l -b -q fitRecoilWm.C+\(\"${IN5}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_5TeV_Eta${i}\",${LUMI5},${i},0\) 
# root -l -b -q fitRecoilWm.C+\(\"${IN5}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_5TeV_Eta${i}\",${LUMI5},${i},0\) 
# # # # # # # # # 5 TeV Z MC
# # root -l -b -q fitRecoilZmm.C+\(\"${IN5}/Zmumu/ntuples/zmm_select.raw.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_5TeV_Eta${i}\",${LUMI5},${i},0\)
# # # # # # # # # 5 TeV Z Data
# # root -l -b -q fitRecoilZmm.C+\(\"${IN5}/Zmumu/ntuples/data_select.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_5TeV_Eta${i}\",${LUMI5},${i},0\)  
# done


###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ROOKEYS PDF
# -------------------------------------------------------------------------------------------------
## 5 TeV W MC
# root -l -b -q fitRecoilWm.C+\(\"${IN5}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_5TeV_Keys\",${LUMI5},0,1\) 
# root -l -b -q fitRecoilWm.C+\(\"${IN5}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_5TeV_Keys\",${LUMI5},0,1\) 
# # 5 TeV Z MC
root -l -b -q fitRecoilZmm.C+\(\"${IN5}/Zmumu/ntuples/zmm_select.raw.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_5TeV_Keys\",${LUMI5},0,1\)
# 5 TeV Z Data
# root -l -b -q fitRecoilZmm.C+\(\"${IN5}/Zmumu/ntuples/data_select.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_5TeV_Keys\",${LUMI5},0,1\) 

###################################################################################################
###################################################################################################



#rm *.so *.d *.pcm
