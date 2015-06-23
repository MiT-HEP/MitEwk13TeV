#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
NTUPLEDIR=/afs/cern.ch/work/c/cmedlock/wz-efficiency
OUTPUTDIR="."

#
# Electron efficiencies
#
#root -l -b -q plotEff.C+\(\"elhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleL1Eff/probes.root\",\"${OUTPUTDIR}/Zee_EleL1Eff\",\"png\",0,0,0,\"Supercluster\",\"L1 trigger\",0.8,1.03\)

#root -l -b -q plotEff.C+\(\"elhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.8,1.03\)

#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleSelEff\",\"png\",0,0,0,\"Supercluster\",\"ID+Iso\",0.8,1.03\)

#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfEff\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.8,1.03\)

root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02\)

#
# Electron (+) efficiencies
#
#root -l -b -q plotEff.C+\(\"elhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff_pos\",\"png\",0,0,1,\"Supercluster\",\"trigger\"\)

#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff_pos\",\"png\",1,0,1,\"Supercluster\",\"GSF+ID+Iso\"\)

#
# Electron (-) efficiencies
#
#root -l -b -q plotEff.C+\(\"elhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff_neg\",\"png\",0,0,-1,\"Supercluster\",\"trigger\"\)

#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff_neg\",\"png\",1,0,-1,\"Supercluster\",\"GSF+ID+Iso\"\)

rm *.so *.d
