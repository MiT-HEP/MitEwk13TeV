#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
NTUPLEDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency
OUTPUTDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency-results

# integrated luminosity for data
LUMI=42.0

#
# Muon efficiencies
#

# trigger efficiency
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuHLTEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuHLTEff\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI}\)
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZmm_MuHLTEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuHLTEff\",\"png\",0,0,0,\"Muon\",\"trigger\",0.7,1.03,${LUMI}\)

# id+iso efficiency
root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuSelEff\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI}\)
root -l -b -q plotEff.C+\(\"mupteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuSelEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuSelEff\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuSelEff/probes.root\"\)

# tracking efficiency
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuTrkEff\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"mupteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuTrkEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuTrkEff\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\"\)
# with a single bin
#root -l -b -q plotEff.C+\(\"musinglebin.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\",\"${OUTPUTDIR}/Zmm_MuTrkEff_singleBin\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"musinglebin.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuTrkEff/probes.root\",\"${OUTPUTDIR}/DataZmm_MuTrkEff_singleBin\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuTrkEff/probes.root\"\)

# tracking efficiency v2
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuTrkEff_v2/probes.root\",\"${OUTPUTDIR}/Zmm_MuTrkEff_v2\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"mupteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuTrkEff_v2/probes.root\",\"${OUTPUTDIR}/DataZmm_MuTrkEff_v2\",\"png\",0,0,0,\"Muon\",\"tracking\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuTrkEff_v2/probes.root\"\)

# standalone efficiency
root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm_MuStaEff_iso\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI}\)
root -l -b -q plotEff.C+\(\"mupteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuStaEff_iso/probes.root\",\"${OUTPUTDIR}/DataZmm_MuStaEff_iso\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuStaEff_iso/probes.root\"\)

# standalone+id+iso efficiency
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm_MuSelStaEff_iso\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"mupteta.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/DataZmm_MuSelStaEff_iso\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuSelStaEff_iso/probes.root\"\)
# with a single bin
#root -l -b -q plotEff.C+\(\"musinglebin.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/Zmm_MuSelStaEff_iso_singleBin\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"musinglebin.bins\",2,1,2,1,\"${NTUPLEDIR}/DataZmm_MuSelStaEff_iso/probes.root\",\"${OUTPUTDIR}/DataZmm_MuSelStaEff_iso_singleBin\",\"png\",0,0,0,\"Muon\",\"stand-alone+ID+Iso\",0.7,1.02,${LUMI},\"${NTUPLEDIR}/Zmm_MuSelStaEff_iso/probes.root\"\)

# standalone efficiency with deltaR matching
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/Zmm_MuStaEff_iso_v2/probes.root\",\"${OUTPUTDIR}/Zmm_MuStaEff_iso_v2\",\"png\",0,0,0,\"Muon\",\"stand-alone-deltaR\",0.7,1.02,${LUMI}\)
#root -l -b -q plotEff.C+\(\"mupteta.bins\",0,0,0,0,\"${NTUPLEDIR}/DataZmm_MuStaEff_iso_v2/probes.root\",\"${OUTPUTDIR}/DataZmm_MuStaEff_iso_v2\",\"png\",0,0,0,\"Muon\",\"stand-alone-deltaR\",0.7,1.02,${LUMI}\)

rm *.so *.d
