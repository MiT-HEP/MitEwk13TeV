{    

    gROOT->Macro("../Utils/CPlot.cc++");
    gROOT->Macro("../Utils/MitStyleRemix.cc++");  
               
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
