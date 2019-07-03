{     
//  gROOT->Macro("RooVoigtianShape.cc+");
  gROOT->Macro("RooCMSShape.cc+");
  
  gROOT->Macro("../../CPlot.cc+");
  gROOT->Macro("../../Utils/MitStyleRemix.cc+");
  gROOT->Macro("../../Utils/CEffUser1D.cc+");
  gROOT->Macro("../../Utils/CEffUser2D.cc+");
              
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
