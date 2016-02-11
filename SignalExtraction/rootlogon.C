{    

  if(gSystem->Getenv("CMSSW_VERSION")) {
    TString rfitpath("/afs/cern.ch/cms/$SCRAM_ARCH/lcg/roofit/6.02.00-cms/include/");
    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I";
    path += rfitpath;
    gSystem->SetIncludePath(path.Data());

    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }

    gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconAnaDataFormats.so");

    gROOT->Macro("$CMSSW_BASE/src/BaconAna/macros/setRootEnv.C+");

  gROOT->Macro("../Utils/RooVoigtianShape.cc+");
  gROOT->Macro("../Utils/RooCMSShape.cc+");

  gROOT->Macro("../Utils/CPlot.cc++");
  gROOT->Macro("../Utils/MitStyleRemix.cc++");  

  gROOT->Macro("RooVoigtianShape.cc+");
  gROOT->Macro("RooCMSShape.cc+");

  gROOT->Macro("CEffUser1D.cc+");
  gROOT->Macro("CEffUser2D.cc+");
  gROOT->Macro("muresolution_run2.cc+");
  gROOT->Macro("rochcor2015.cc+");
  	
  TString path = gSystem->GetIncludePath();
  path += " -I ../EleScale/"
  gSystem->SetIncludePath(path.Data());
  gROOT->Macro("EnergyScaleCorrection_class.cc+");

  }
               
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
