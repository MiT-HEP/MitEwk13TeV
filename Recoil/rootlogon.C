{    

 if(gSystem->Getenv("CMSSW_VERSION")) {
    //TString rfitpath("/afs/cern.ch/cms/$SCRAM_ARCH/lcg/roofit/6.02.00-cms/include/");
    //TString path = gSystem->GetIncludePath();
    //path += "-I. -I$ROOTSYS/src -I";
    //path += rfitpath;
    //gSystem->SetIncludePath(path.Data());

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
  gROOT->Macro("../Utils/PdfDiagonalizer.cc++");  

  gROOT->Macro("../SignalExtraction/RooVoigtianShape.cc+");
  gROOT->Macro("../SignalExtraction/RooCMSShape.cc+");

  gROOT->Macro("../SignalExtraction/CEffUser1D.cc+");
  gROOT->Macro("../SignalExtraction/CEffUser2D.cc+");

  gROOT->Macro("../SignalExtraction/muresolution_run2r.cc+");
  gROOT->Macro("../SignalExtraction/rochcor2015r.cc+");
 
  //gROOT->Macro("../Utils/RecoilCorrector_asym2.hh++"); 
  {  
    //TString path = gSystem->GetIncludePath();
    //path += " -I../EleScale/ ";
    //gSystem->SetIncludePath(path.Data());
    gSystem->AddIncludePath("-I../EleScale");
    gInterpreter->AddIncludePath("../EleScale");
    gROOT->SetMacroPath(TString(gROOT->GetMacroPath()) + ":../EleScale");
    gROOT->Macro("../SignalExtraction/EnergyScaleCorrection_class.cc+");
  }

  }
               
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
