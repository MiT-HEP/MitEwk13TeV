#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
NTUPLEDIR=/scratch/klawhorn/EWKAnaR12a/Efficiency

# generate templates
#root -l -q makePsExpTemplates.C+\(\"musel.bins\",1,1,1,1,\"${NTUPLEDIR}/R12a_MuSelEff/probes.root\",\"CB_MuSelEff/analysis\",\"png\",0,0,0\)
#root -l -q makePsExpTemplates.C+\(\"mutrk.bins\",1,1,1,1,\"${NTUPLEDIR}/R12a_MuTrkEff/probes.root\",\"CB_MuTrkEff/analysis\",\"png\",0,0,0\)
#root -l -q makePsExpTemplates.C+\(\"musta.bins\",1,2,1,2,\"${NTUPLEDIR}/R12a_MuStaEff_iso/probes.root\",\"CB_MuStaEff_iso/analysis\",\"png\",0,0,0\)

#root -l -q makePsExpTemplates.C+\(\"musta.bins\",2,1,2,6,\"${NTUPLEDIR}/R12a_MuStaEff_iso/probes.root\",\"BK_MuStaEff_iso\",\"png\",0,0,0,\"/scratch/klawhorn/EWKAnaR12a/Efficiency/Zmm_MuSelEff/probes.root\",\"../muStaBkg/probes.root\"\)

#root -l -q makePsExpTemplates.C+\(\"musel.bins\",2,4,2,4,\"${NTUPLEDIR}/R12a_MuSelEff/probes.root\",\"Lin_MuSelEff/analysis\",\"png\",0,0,0,\"/scratch/klawhorn/EWKAnaR12a/Efficiency/Zmm_MuSelEff/probes.root\"\)
#root -l -q makePsExpTemplates.C+\(\"mutrk.bins\",2,4,2,4,\"${NTUPLEDIR}/R12a_MuTrkEff/probes.root\",\"Lin_MuTrkEff/analysis\",\"png\",0,0,0,\"/scratch/klawhorn/EWKAnaR12a/Efficiency/Zmm_MuTrkEff/probes.root\"\)

#root -l -q makePsExpTemplates.C+\(\"musta.bins\",2,4,2,4,\"${NTUPLEDIR}/R12a_MuStaEff_iso/probes.root\",\"Lin_MuStaEff_iso/analysis\",\"png\",0,0,0,\"/scratch/klawhorn/EWKAnaR12a/Efficiency/Zmm_MuStaEff_iso/probes.root\"\)

#root -l -q makePsExpTemplates.C+\(\"elgsfsel.bins\",2,6,2,6,\"${NTUPLEDIR}/R12a_EleGsfSelEff/probes.root\",\"BK_EleEff\",\"png\",0,0,0,\"/scratch/klawhorn/EWKAnaR12a/Efficiency/Zee_EleGsfSelEff/probes.root\",\"../qcd/probes.root\"\)

#root -l -q makePsExpTemplates.C+\(\"elgsfsel.bins\",2,2,2,2,\"${NTUPLEDIR}/R12a_EleGsfSelEff/probes.root\",\"EleGsfSelEff\",\"png\",0,0,0,\"/scratch/klawhorn/EWKAnaR12a/Efficiency/Zee_EleGsfSelEff/probes.root\"\)

# generate pseudodata from templates

#root -l -q makePseudoData.C+\(\"BK_EleEff/analysis/plots/\",\"etapt_0\",1000,\"/scratch/klawhorn/EffSysStore/BK_EleEff/\"\)
#root -l -q makePseudoData.C+\(\"BK_EleEff/analysis/plots/\",\"etapt_1\",1000,\"/scratch/klawhorn/EffSysStore/BK_EleEff/\"\)
#root -l -q makePseudoData.C+\(\"BK_EleEff/analysis/plots/\",\"etapt_2\",1000,\"/scratch/klawhorn/EffSysStore/BK_EleEff/\"\)
#root -l -q makePseudoData.C+\(\"BK_EleEff/analysis/plots/\",\"etapt_3\",1000,\"/scratch/klawhorn/EffSysStore/BK_EleEff/\"\)
#root -l -q makePseudoData.C+\(\"BK_EleEff/analysis/plots/\",\"etapt_4\",1000,\"/scratch/klawhorn/EffSysStore/BK_EleEff/\"\)
#root -l -q makePseudoData.C+\(\"BK_EleEff/analysis/plots/\",\"etapt_5\",1000,\"/scratch/klawhorn/EffSysStore/BK_EleEff/\"\)
#root -l -q makePseudoData.C+\(\"BK_EleEff/analysis/plots/\",\"etapt_6\",1000,\"/scratch/klawhorn/EffSysStore/BK_EleEff/\"\)
#root -l -q makePseudoData.C+\(\"BK_EleEff/analysis/plots/\",\"etapt_7\",1000,\"/scratch/klawhorn/EffSysStore/BK_EleEff/\"\)
#root -l -q makePseudoData.C+\(\"BK_EleEff/analysis/plots/\",\"etapt_8\",1000,\"/scratch/klawhorn/EffSysStore/BK_EleEff/\"\)
root -l -q makePseudoData.C+\(\"BK_EleEff/analysis/plots/\",\"etapt_9\",1000,\"/scratch/klawhorn/EffSysStore/BK_EleEff/\"\)


#root -l -q makePseudoData.C+\(\"BK_MuStaEff_iso/analysis/plots/\",\"etapt_0\",1000,\"/scratch/klawhorn/EffSysStore/BK_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"BK_MuStaEff_iso/analysis/plots/\",\"etapt_1\",1000,\"/scratch/klawhorn/EffSysStore/BK_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"BK_MuStaEff_iso/analysis/plots/\",\"etapt_2\",1000,\"/scratch/klawhorn/EffSysStore/BK_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"BK_MuStaEff_iso/analysis/plots/\",\"etapt_3\",1000,\"/scratch/klawhorn/EffSysStore/BK_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"BK_MuStaEff_iso/analysis/plots/\",\"etapt_4\",1000,\"/scratch/klawhorn/EffSysStore/BK_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"BK_MuStaEff_iso/analysis/plots/\",\"etapt_5\",1000,\"/scratch/klawhorn/EffSysStore/BK_MuStaEff_iso/\"\)

#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_0\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_1\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_2\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_3\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_4\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_5\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)

#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_6\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_7\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_8\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_9\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_10\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_11\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)

#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_12\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_13\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_14\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_15\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_16\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuSelEff/analysis/plots/\",\"etapt_17\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuSelEff/\"\)


#root -l -q makePseudoData.C+\(\"CB_MuStaEff_iso/analysis/plots/\",\"etapt_0\",1000,\"/scratch/klawhorn/EffSysStore/CB_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"CB_MuStaEff_iso/analysis/plots/\",\"etapt_1\",1000,\"/scratch/klawhorn/EffSysStore/CB_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"CB_MuStaEff_iso/analysis/plots/\",\"etapt_2\",1000,\"/scratch/klawhorn/EffSysStore/CB_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"CB_MuStaEff_iso/analysis/plots/\",\"etapt_3\",1000,\"/scratch/klawhorn/EffSysStore/CB_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"CB_MuStaEff_iso/analysis/plots/\",\"etapt_4\",1000,\"/scratch/klawhorn/EffSysStore/CB_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"CB_MuStaEff_iso/analysis/plots/\",\"etapt_5\",1000,\"/scratch/klawhorn/EffSysStore/CB_MuStaEff_iso/\"\)

#root -l -q makePseudoData.C+\(\"Lin_MuStaEff_iso/analysis/plots/\",\"etapt_0\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuStaEff_iso/analysis/plots/\",\"etapt_1\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuStaEff_iso/analysis/plots/\",\"etapt_2\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuStaEff_iso/analysis/plots/\",\"etapt_3\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuStaEff_iso/analysis/plots/\",\"etapt_4\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuStaEff_iso/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuStaEff_iso/analysis/plots/\",\"etapt_5\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuStaEff_iso/\"\)

#root -l -q makePseudoData.C+\(\"Lin_MuTrkEff/analysis/plots/\",\"etapt_0\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuTrkEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuTrkEff/analysis/plots/\",\"etapt_1\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuTrkEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuTrkEff/analysis/plots/\",\"etapt_2\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuTrkEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuTrkEff/analysis/plots/\",\"etapt_3\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuTrkEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuTrkEff/analysis/plots/\",\"etapt_4\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuTrkEff/\"\)
#root -l -q makePseudoData.C+\(\"Lin_MuTrkEff/analysis/plots/\",\"etapt_5\",1000,\"/scratch/klawhorn/EffSysStore/Lin_MuTrkEff/\"\)

# actually do pseudofits
#root -l doPseudoFits.C+\(\"/scratch/klawhorn/EffSysStore/etapt_0_0000.dat\",\"CB_MuSelEff/analysis/plots/etapt_0.root\",2,1,2,1,\"test\",0,0,\"/scratch/klawhorn/EWKAnaR12a/Efficiency/Zmm_MuSelEff/probes.root\"\)

#rm *so *d