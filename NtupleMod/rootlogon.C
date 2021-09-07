{    
  if(gSystem->Getenv("CMSSW_VERSION")) {
    gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconAnaDataFormats.so");
    gROOT->Macro("../Utils/PdfDiagonalizer.cc++");  
    gROOT->Macro("../Utils/CEffUser2D.cc+");
  }          
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
