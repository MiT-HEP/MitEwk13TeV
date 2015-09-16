#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
NTUPLEDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency
OUTPUTDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency-results

# integrated luminosity for data
LUMI=42.0

#
# Electron efficiencies
#

# trigger efficiency
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)
# synchronization with Kevin
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff_cutonLepInfo_oppCharge_noSCmatch/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff_binsFromKevin\",\"png\",1,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff_cutonLepInfo_oppCharge_noSCmatch/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_binsFromKevin\",\"png\",1,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)

# reco+id+iso efficiency
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)
# with a single bin
#root -l -b -q plotEff.C+\(\"elbarrelendcap.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff_EBEE\",\"png\",1,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elbarrelendcap.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff_EBEE\",\"png\",1,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)

# supercluster efficiency
#root -l -b -q plotEff.C+\(\"elsupercluster.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleSCEff/probes.root\",\"${OUTPUTDIR}/Zee_EleSCEff\",\"png\",1,0,0,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI}\)
root -l -b -q plotEff.C+\(\"elsupercluster.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleSCEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleSCEff\",\"png\",1,0,0,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI},\"${NTUPLEDIR}/Zee_EleSCEff/probes.root\"\)

# gsf efficiency
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfEff\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfEff\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\"\)

# id+iso efficiency
#root -l -b -q plotEff.C+\(\"elsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleSelEff\",\"png\",1,0,0,\"Supercluster\",\"ID+Iso\",0.5,1.02,${LUMI}\)
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
