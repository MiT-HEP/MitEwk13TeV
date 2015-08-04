#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
NTUPLEDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency
OUTPUTDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency-results

# integrated luminosity for data
LUMI=40.0

#
# Electron efficiencies
#

#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)

#root -l -b -q plotEff.C+\(\"elfineptbins.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff_finepTbins\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elfineptbins.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_finepTbins\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)

#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI}\)
root -l -b -q plotEff_tmp.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)

#root -l -b -q plotEff.C+\(\"elfineptbins.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff_finepTbins\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elfineptbins.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff_finepTbins\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)

#
# Electron (+) efficiencies
#
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff_pos\",\"png\",0,0,1,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_pos\",\"png\",0,0,1,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)

#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff_pos\",\"png\",0,0,1,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff_pos\",\"png\",0,0,1,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)

#
# Electron (-) efficiencies
#
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff_neg\",\"png\",0,0,-1,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_neg\",\"png\",0,0,-1,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)

#root -l -b -q plotEff.C+\(\"elpteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff_neg\",\"png\",0,0,-1,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elpteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff_neg\",\"png\",0,0,-1,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)

rm *.so *.d
