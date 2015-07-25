#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
NTUPLEDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency
OUTPUTDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency-results

# integrated luminosity for data
LUMI=7.3

#
# Electron efficiencies
#

#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)

#root -l -b -q plotEff.C+\(\"elhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff_finepTbins\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_finepTbins\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15,${LUMI}\)

#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)

#
# Debugging trigger efficiency
#
#root -l -b -q plotEff.C+\(\"elhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15\)
#root -l -b -q plotEff_afiq.C+\(\"elhlt.bins\",0,0,0,0,\"/afs/cern.ch/work/a/afiqaize/public/files/forCatherine/tnp_sin_ele23.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_afiq\",\"png\",0,0,0,\"Supercluster\",\"trigger\",0.0,1.15\)

#
# Debugging bad fit from 1.566 < eta < 2.0 in GsfSel efficiency
#
#root -l -b -q plotEff.C+\(\"eltest_1bin.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff_1bin\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.0,1.3,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)
#root -l -b -q plotEff.C+\(\"eltest_2bins.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff_2bins\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.0,1.3,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)

#
# Testing different signal model to compare with Xinmei
#
#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",1,1,1,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff_CB\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)

#root -l -b -q plotEff.C+\(\"elhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff_finepTbins\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elhlt.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff_finepTbins\",\"png\",0,0,0,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)

#
# Electron (+) efficiencies
#
#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff_pos\",\"png\",0,0,1,\"Supercluster\",\"trigger\",0.0,1.15\)
#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_pos\",\"png\",0,0,1,\"Supercluster\",\"trigger\",0.0,1.15\)

#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff_pos\",\"png\",0,0,1,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff_pos\",\"png\",0,0,1,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)

#
# Electron (-) efficiencies
#
#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/Zee_EleHLTEff_neg\",\"png\",0,0,-1,\"Supercluster\",\"trigger\",0.0,1.15\)
#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZee_EleHLTEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_neg\",\"png\",0,0,-1,\"Supercluster\",\"trigger\",0.0,1.15\)

#root -l -b -q plotEff.C+\(\"elgsfsel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff_neg\",\"png\",0,0,-1,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI}\)
root -l -b -q plotEff.C+\(\"elgsfsel.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZee_EleGsfSelEff/probes.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff_neg\",\"png\",0,0,-1,\"Supercluster\",\"GSF+ID+Iso\",0.5,1.02,${LUMI},\"${NTUPLEDIR}/Zee_EleGsfSelEff/probes.root\"\)

rm *.so *.d
