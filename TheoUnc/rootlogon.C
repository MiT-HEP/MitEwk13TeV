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
    
    gSystem->Load("$CMSSW_BASE/lib/slc6_amd64_gcc481/libBaconAnaDataFormats.so");
    gSystem->Load("$LHAPDFSYS/lib/libLHAPDF.so");
    
    gROOT->Macro("$CMSSW_BASE/src/BaconAna/macros/setRootEnv.C+");
    gROOT->Macro("../Utils/MitStyleRemix.cc++");
    gROOT->Macro("compare.cc++");
    gSystem->AddIncludePath("-I/afs/cern.ch/sw/lcg/external/Boost/1.55.0_python2.7/x86_64-slc6-gcc47-opt/boost-1_55/boost/");
    gSystem->AddIncludePath("-I$LHAPDFSYS/include/");

  } 
}
