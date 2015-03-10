{     
  gROOT->Macro("RooVoigtianShape.cc+");
  gROOT->Macro("RooCMSShape.cc+");
  
  gROOT->Macro("CPlot.cc+");
  gROOT->Macro("MitStyleRemix.cc+");
  gROOT->Macro("CEffUser1D.cc+");
  gROOT->Macro("CEffUser2D.cc+");
              
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
