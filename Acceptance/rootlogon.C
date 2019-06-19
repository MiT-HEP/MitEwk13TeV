{    
  if(gSystem->Getenv("CMSSW_VERSION")) {    
    // TString rfitpath("/afs/cern.ch/cms/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms9/include/");
    // TString path = gSystem->GetIncludePath();
    // path += "-I. -I$ROOTSYS/src -I";
    // path += rfitpath;
    // gSystem->SetIncludePath(path.Data());
    // 
    // TString str = gSystem->GetMakeSharedLib();
    // if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
    //   str.ReplaceAll("g++", "g++ -m32");
    //   gSystem->SetMakeSharedLib(str);
    // }
    
    gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconAnaDataFormats.so");
    
    gROOT->Macro("$CMSSW_BASE/src/BaconAna/macros/setRootEnv.C+");
    // gSystem->SetMakeSharedLib(cmd);
    //#endif
  
  
   }

  gROOT->Macro("../Utils/CPlot.cc++");
  gROOT->Macro("../Utils/MitStyleRemix.cc++");
  gROOT->Macro("CCorrUser2D.cc+");
 
  gROOT->Macro("../Utils/CEffUser1D.cc+");
  gROOT->Macro("../Utils/CEffUser2D.cc+");
               
    {
     //TString path = gSystem->GetIncludePath();
     //path += " -I../EleScale/ ";
     //gSystem->SetIncludePath(path.Data());
      gSystem->AddIncludePath("-I../EleScale");
      gInterpreter->AddIncludePath("../EleScale");
      gROOT->SetMacroPath(TString(gROOT->GetMacroPath()) + ":../EleScale");
      // gROOT->Macro("EnergyScaleCorrection_class.cc+"); // commented out to test the new ones
      gROOT->Macro("EnergyScaleCorrection.cc+");
    }
                     

  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
