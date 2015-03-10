{    
  if(gSystem->Getenv("CMSSW_VERSION")) {
    TString rfitpath("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include");
    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I";
    path += rfitpath;
    gSystem->SetIncludePath(path.Data());
    
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }
    
    gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libEWKAnaNtupler.so");
    
    gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
    gROOT->Macro("../Utils/CPlot.cc+");
    gROOT->Macro("../Utils/MitStyleRemix.cc+"); 
    gROOT->Macro("CEffUser1D.cc+");
    gROOT->Macro("CEffUser2D.cc+");    
  }
               
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
