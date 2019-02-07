#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# # directory of probes ntuples
# NTUPLEDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_Efficiency
# OUTPUTDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_Efficiency-results


NTUPLEDIR_5o=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_Efficiency_5TeV
# NTUPLEDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_Efficiency_matchTrig_5TeV
# OUTPUTDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_Efficiency_matchTrig_5TeV-results


NTUPLEDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2_Efficiency
OUTPUTDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2_Efficiency-results

# integrated luminosity for data
LUMI13=213.1
LUMI5=294.4

#################################################################
# Electron efficiencies
#
# 13 TeV shits
# supercluster efficiency
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,2,1,1,\"${NTUPLEDIR}/Zee_HLTEff_Test/probes.root\",\"${OUTPUTDIR}/Zee/HLTEff_Test\",\"png\",1,0,0,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI},\"${NTUPLEDIR}/Zee_HLTEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_HLTEff_Test/probes.root\",\"${OUTPUTDIR}/Zee/HLTEff_Test\",\"png\",0,0,0,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR}/Zee_HLTEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,2,1,1,\"${NTUPLEDIR}/DataZee_HLTEff_Test/probes.root\",\"${OUTPUTDIR}/DataZee/HLTEff_Test\",\"png\",0,0,0,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR}/Zee_HLTEff_Test/probes.root\"\)

# gsf+is+iso efficiency
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,2,1,1,\"${NTUPLEDIR}/Zee_GSFSelEff_Test/probes.root\",\"${OUTPUTDIR}/Zee/GSFSelEff_Test\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_GSFSelEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_GSFSelEff_Test/probes.root\",\"${OUTPUTDIR}/Zee/GSFSelEff_Test\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR}/Zee_GSFSelEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,2,1,1,\"${NTUPLEDIR}/DataZee_GSFSelEff_Test/probes.root\",\"${OUTPUTDIR}/DataZee/GSFSelEff_Test\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR}/Zee_GSFSelEff_Test/probes.root\"\)

#### Draw Plots
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR/Plots/ZeeHLT\",\"$OUTPUTDIR/Zee/HLTEff_Test/eff.root\",\"$OUTPUTDIR/DataZee/HLTEff_Test/eff.root\",\"EtaBins\",0.00,1.20,$LUMI,\"elfineptbins.bins\"\)
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR/Plots/ZeeGSFSel\",\"$OUTPUTDIR/Zee/GSFSelEff_Test/eff.root\",\"$OUTPUTDIR/DataZee/GSFSelEff_Test/eff.root\",\"EtaBins\",0.00,1.20,$LUMI,\"elfineptbins.bins\"\)

#####################################################################
# Electron efficiencies
#
# 5 TeV shits
# supercluster efficiency
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,2,1,1,\"${NTUPLEDIR}/Zee_HLTEff_Test/probes.root\",\"${OUTPUTDIR}/Zee/HLTEff_Test\",\"png\",1,0,0,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI},\"${NTUPLEDIR}/Zee_HLTEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zee_HLTEff_Test/probes.root\",\"${OUTPUTDIR_5}/Zee/HLTEff_Test\",\"png\",0,0,0,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI5},\"${NTUPLEDIR_5}/Zee_HLTEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins_5.bins\",2,2,1,1,\"${NTUPLEDIR_5}/DataZee_HLTEff_Test/probes.root\",\"${OUTPUTDIR_5}/DataZeeFit/HLTEff_Test\",\"png\",0,0,0,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI5},\"${NTUPLEDIR_5}/Zee_HLTEff_Test/probes.root\"\)

# gsf+is+iso efficiency
# root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,2,1,1,\"${NTUPLEDIR}/Zee_GSFSelEff_Test/probes.root\",\"${OUTPUTDIR}/Zee/GSFSelEff_Test\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_GSFSelEff_Test/probes.root\"\)
# root -l -b -q plotEff.C+\(\"elfineptbins_5.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zee_GSFSelEff_Test/probes.root\",\"${OUTPUTDIR_5}/Zee/GSFSelEff_Test\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee_GSFSelEff_Test/probes.root\"\)
root -l -b -q plotEff.C+\(\"elfineptbins_5.bins\",2,2,1,1,\"${NTUPLEDIR_5}/DataZee_GSFSelEff_Test/probes.root\",\"${OUTPUTDIR_5}/DataZeeFit/GSFSelEff_Test\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee_GSFSelEff_Test/probes.root\"\)


 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR/Plots/ZeeHLT\",\"$OUTPUTDIR/Zee/HLTEff_Test/eff.root\",\"$OUTPUTDIR/DataZee/HLTEff_Test/eff.root\",\"EtaBins\",0.00,1.20,$LUMI,\"elfineptbins.bins\"\)
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR/Plots/ZeeGSFSel\",\"$OUTPUTDIR/Zee/GSFSelEff_Test/eff.root\",\"$OUTPUTDIR/DataZee/GSFSelEff_Test/eff.root\",\"EtaBins\",0.00,1.20,$LUMI,\"elfineptbins.bins\"\)




# id+iso efficiency
# root -l -b -q plotEff.C+\(\"elsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleSelEff\",\"png\",1,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elsel.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleSelEff\",\"png\",1,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\"\)

# id+iso efficiency with tight ID
#root -l -b -q plotEff.C+\(\"elsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleSelEff_cutonLepInfo_oppCharge_noSCmatch/probes.root\",\"${OUTPUTDIR}/Zee_EleSelEff_binsFromKevin\",\"png\",1,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elsel.bins\",2,1,2,2,\"${NTUPLEDIR}/DataZee_EleSelEff_cutonLepInfo_oppCharge_noSCmatch/probes.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_binsFromKevin\",\"png\",1,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleSelEff_cutonLepInfo_oppCharge_noSCmatch/probes.root\"\)

# trigger efficiency for specific runs
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff_singleRunBins\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.5,1.02,${LUMI}\)

#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff_run251244/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251244\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff_run251251/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251251\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff_run251252/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251252\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff_run251561/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251561\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff_run251562/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251562\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff_run251643/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251643\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff_run251721/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251721\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff_run251883/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251883\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.5,1.02,${LUMI}\)

# id+iso efficiency for specific runs
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleSelEff_singleRunBins\",\"png\",0,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI}\)

#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleSelEff_run251244/probes.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251244\",\"png\",0,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\"\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleSelEff_run251251/probes.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251251\",\"png\",0,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\"\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleSelEff_run251252/probes.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251252\",\"png\",0,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\"\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleSelEff_run251561/probes.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251561\",\"png\",0,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\"\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleSelEff_run251562/probes.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251562\",\"png\",0,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\"\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleSelEff_run251643/probes.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251643\",\"png\",0,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\"\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleSelEff_run251721/probes.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251721\",\"png\",0,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\"\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleSelEff_run251883/probes.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251883\",\"png\",0,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\"\)

rm *.so *.d
