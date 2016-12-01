{    

    gROOT->Macro("../Utils/CPlot.cc++");
    gROOT->Macro("../Utils/MitStyleRemix.cc++");  

    gSystem->AddIncludePath("-I$CMSSW_BASE/src/HiggsAnalysis/GBRLikelihood/interface/");    
    gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisGBRLikelihood.so");
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}


