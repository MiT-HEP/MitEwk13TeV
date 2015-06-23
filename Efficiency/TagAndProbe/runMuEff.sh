#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
NTUPLEDIR=/afs/cern.ch/work/c/cmedlock/wz-efficiency

#
# Muon efficiencies
#
#root -l -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuL1Eff/probes.root\",\"Zmm_MuL1Eff\",\"png\",0,0,0,\"Muon\",\"L1 trigger\",0.7,1.03\)

#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuHLTEff/probes.root\",\"Zmm_MuHLTEff\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03\)

#root -l -b -q plotEff.C+\(\"musel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\",\"Zmm_MuSelEff\",\"png\",0,0,0,\"Muon\",\"selection\",0.7,1.02\)

#root -l -b -q plotEff.C+\(\"mutrk.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\",\"Zmm_MuTrkEff\",\"png\",1,0,0,\"Muon\",\"tracking\",0.7,1.02\)

#root -l -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuStaEff/probes.root\",\"Zmm_MuStaEff\",\"png\",0,0,0,\"Muon\",\"stand-alone (no iso)\",0.7,1.03\)

#root -l -b -q plotEff.C+\(\"musta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuStaEff_iso/probes.root\",\"Zmm_MuStaEff_iso\",\"png\",1,0,0,\"Muon\",\"stand-alone\",0.7,1.02\)

#
# Muon (+) efficiencies
#
#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuHLTEff/probes.root\",\"Zmm_MuHLTEff_pos\",\"png\",0,0,1,\"Muon\",\"trigger\",0.7,1.02\)

#root -l -b -q plotEff.C+\(\"musel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\",\"Zmm_MuSelEff_pos\",\"png\",0,0,1,\"Muon\",\"selection\",0.7,1.02\)

#root -l -b -q plotEff.C+\(\"mutrk.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\",\"Zmm_MuTrkEff_pos\",\"png\",1,0,1,\"Muon\",\"tracking\",0.7,1.02\)

#root -l -q plotEff.C+\(\"musta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuStaEff/probes.root\",\"Zmm_MuStaEff_pos\",\"png\",1,0,1,\"Muon\",\"stand-alone (no iso)\",0.7,1.03\)

#root -l -b -q plotEff.C+\(\"musta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuStaEff_iso/probes.root\",\"Zmm_MuStaEff_iso_pos\",\"png\",1,0,1,\"Muon\",\"stand-alone\",0.7,1.02\)

#
# Muon (-) efficiencies
#
#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuHLTEff/probes.root\",\"Zmm_MuHLTEff_neg\",\"png\",0,0,-1,\"Muon\",\"trigger\",0.7,1.02\)

#root -l -b -q plotEff.C+\(\"musel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\",\"Zmm_MuSelEff_neg\",\"png\",0,0,-1,\"Muon\",\"selection\",0.7,1.02\)

#root -l -b -q plotEff.C+\(\"mutrk.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\",\"Zmm_MuTrkEff_neg\",\"png\",1,0,-1,\"Muon\",\"tracking\",0.7,1.02\)

#root -l -q plotEff.C+\(\"musta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuStaEff/probes.root\",\"Zmm_MuStaEff_neg\",\"png\",1,0,-1,\"Muon\",\"stand-alone (no iso)\",0.7,1.03\)

root -l -b -q plotEff.C+\(\"musta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuStaEff_iso/probes.root\",\"Zmm_MuStaEff_iso_neg\",\"png\",1,0,-1,\"Muon\",\"standalone\",0.7,1.02\)

rm *.so *.d
