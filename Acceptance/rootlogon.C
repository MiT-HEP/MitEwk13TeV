{    

  if(gSystem->Getenv("CMSSW_VERSION")) {
    TString rfitpath("/afs/cern.ch/cms/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms9/include/");
    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I";
    path += rfitpath;
    gSystem->SetIncludePath(path.Data());

    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }

    gSystem->Load("$CMSSW_BASE/lib/slc6_amd64_gcc493/libBaconAnaDataFormats.so");

    gROOT->Macro("$CMSSW_BASE/src/BaconAna/macros/setRootEnv.C+");
    gROOT->Macro("../Utils/CPlot.cc+");
    gROOT->Macro("../Utils/MitStyleRemix.cc+");
    gROOT->Macro("CEffUser1D.cc+");
    gROOT->Macro("CEffUser2D.cc+");    

   {
//    TString path = gSystem->GetIncludePath();
//    path += " -I./ ";
//    gSystem->SetIncludePath(path.Data());
    gROOT->Macro("../EleScale/EnergyScaleCorrection_class.cc+");
   }

  }
  
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
