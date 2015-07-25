#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
NTUPLEDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency
OUTPUTDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency-results

# integrated luminosity for data
LUMI=7.3

#
# Muon efficiencies
#

#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuHLTEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuHLTEff\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI}\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZmm_MuHLTEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuHLTEff\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI}\)

#root -l -b -q plotEff.C+\(\"musel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuHLTEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuHLTEff_finepTbins\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI}\)
#root -l -b -q plotEff.C+\(\"musel.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZmm_MuHLTEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuHLTEff_finepTbins\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI}\)

#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuSelEff\",\"png\",0,0,0,\"Muon\",\"selection\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuSelEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuSelEff\",\"png\",0,0,0,\"Muon\",\"selection\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\"\)

root -l -b -q plotEff.C+\(\"musel.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuSelEff_finepTbins\",\"png\",0,0,0,\"Muon\",\"selection\",0.7,1.02,${LUMI}\)
root -l -b -q plotEff.C+\(\"musel.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuSelEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuSelEff_finepTbins\",\"png\",0,0,0,\"Muon\",\"selection\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\"\)

#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuTrkEff\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuTrkEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuTrkEff\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\"\)

#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm_MuStaEff_iso\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuStaEff_iso/probes.root\",\"${OUTPUTDIR}/DataZmm_MuStaEff_iso\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuStaEff_iso/probes.root\"\)

#
# Muon (+) efficiencies
#
#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuHLTEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuHLTEff_pos\",\"png\",0,0,1,\"Muon\",\"trigger\",0.7,1.03,${LUMI}\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZmm_MuHLTEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuHLTEff_pos\",\"png\",0,0,1,\"Muon\",\"trigger\",0.7,1.03,${LUMI},\"${NTUPLEDIR}/Zmm_MuHLTEff/probes.root\"\)

#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuSelEff_pos\",\"png\",0,0,1,\"Muon\",\"selection\",0.7,1.02\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuSelEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuSelEff_pos\",\"png\",0,0,1,\"Muon\",\"selection\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\"\)

#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuTrkEff_pos\",\"png\",0,0,1,\"Muon\",\"tracking\",0.7,1.02\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuTrkEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuTrkEff_pos\",\"png\",0,0,1,\"Muon\",\"tracking\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\"\)

#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm_MuStaEff_iso_pos\",\"png\",0,0,1,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuStaEff_iso/probes.root\",\"${OUTPUTDIR}/DataZmm_MuStaEff_iso_pos\",\"png\",0,0,1,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuStaEff_iso/probes.root\"\)

#
# Muon (-) efficiencies
#
#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuHLTEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuHLTEff_neg\",\"png\",0,0,-1,\"Muon\",\"trigger\",0.7,1.03,${LUMI}\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZmm_MuHLTEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuHLTEff_neg\",\"png\",0,0,-1,\"Muon\",\"trigger\",0.7,1.03,${LUMI},\"${NTUPLEDIR}/Zmm_MuHLTEff/probes.root\"\)

#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuSelEff_neg\",\"png\",0,0,-1,\"Muon\",\"selection\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuSelEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuSelEff_neg\",\"png\",0,0,-1,\"Muon\",\"selection\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\"\)

#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuTrkEff_neg\",\"png\",0,0,-1,\"Muon\",\"tracking\",0.7,1.02\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuTrkEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuTrkEff_neg\",\"png\",0,0,-1,\"Muon\",\"tracking\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\"\)

#root -l -b -q plotEff.C+\(\"muhlt.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm_MuStaEff_iso_neg\",\"png\",0,0,-1,\"Muon\",\"standalone\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"muhlt.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuStaEff_iso/probes.root\",\"${OUTPUTDIR}/DataZmm_MuStaEff_iso_neg\",\"png\",0,0,-1,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuStaEff_iso/probes.root\"\)

rm *.so *.d
