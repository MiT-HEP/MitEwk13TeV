#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
NTUPLEDIR=/data/blue/ksung/TagAndProbeExample

#
# Electron selection efficiencies (for H->ZZ->2l2v) in Summer11 Madgraph sample
#
root -l -q plotEff.C+\(\"example.bins\",0,0,0,0,\"${NTUPLEDIR}/zee_madgraph_probes.root\",\"MC\",\"png\",1,1,0\)

#
# Electron selection efficiencies (for H->ZZ->2l2v) in Run2011A data
#
root -l -q plotEff.C+\(\"example.bins\",2,1,2,2,\"${NTUPLEDIR}/zee_data_probes.root\",\"Data\",\"png\",1,-1,0,\"${NTUPLEDIR}/zee_madgraph_probes.root\"\)


rm *.so *.d
