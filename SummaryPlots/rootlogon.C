{    
  gROOT->Macro("MitStyleRemix.cc++");  

  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
