//================================================================================================
//
// Compute cross sections and produce summary plots
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "CSummaryPlot.hh"          // Class to for making cross section summary plots
#endif

void xsec(const TString infilename="input.txt", const TString outputDir=".")
{   
  // theory predictions
  const Double_t theory_xsW   = 12.503;   const Double_t theory_xsWerr   = theory_xsW*0.27/10.44;//0.018;
  const Double_t theory_xsWp  =  7.322;   const Double_t theory_xsWperr  = theory_xsWp*0.17/6.15;//0.015;
  const Double_t theory_xsWm  =  5.181;   const Double_t theory_xsWmerr  = theory_xsWm*0.11/4.29;//0.010;
  const Double_t theory_xsZ   =  1.1323;  const Double_t theory_xsZerr   = theory_xsZ*0.03/0.97;//0.0009;
  
  const Double_t theory_xsWZ  = 11.04;    const Double_t theory_xsWZerr  = theory_xsWZ*0.04/10.74;//0.02;
  const Double_t theory_xsWpm =  1.413;   const Double_t theory_xsWpmerr = theory_xsWpm*0.01/1.43;//0.004;
  
  
  CSummaryPlot::sOutDir = outputDir;
  TString format("png");
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  Double_t nWmp,      nWmm,      nWm,      nZmm,      nWep,      nWem,      nWe,      nZee;
  Double_t nWmpErr,   nWmmErr,   nWmErr,   nZmmErr,   nWepErr,   nWemErr,   nWeErr,   nZeeErr;
  Double_t accWmp,    accWmm,    accWm,    accZmm,    accWep,    accWem,    accWe,    accZee;
  Double_t accWmpErr, accWmmErr, accWmErr, accZmmErr, accWepErr, accWemErr, accWeErr, accZeeErr;
  Double_t lumi, lumiErr;
  
  Double_t errWmp=0, errWmm=0, errWm=0, errZmm=0, errWmpm=0, errWZm=0;
  Double_t errWep=0, errWem=0, errWe=0, errZee=0, errWepm=0, errWZe=0; 
  
  ifstream ifs;
  ifs.open(infilename.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    
    if(state==0) {
      stringstream ss(line);
      ss >> lumi >> lumiErr;
      state++;
      
    } else if(state==1) {
      stringstream ss(line);
      string label;
      ss >> label >> nWmp >> nWmm >> nWm >> nZmm >> nWep >> nWem >> nWe >> nZee;
      state++;
      
    } else if(state==2) {
      stringstream ss(line);
      string label;
      ss >> label >> nWmpErr >> nWmmErr >> nWmErr >> nZmmErr >> nWepErr >> nWemErr >> nWeErr >> nZeeErr;
      state++;
      
    } else if(state==3) {
      stringstream ss(line);
      string label;
      ss >> label >> accWmp >> accWmm >> accWm >> accZmm >> accWep >> accWem >> accWe >> accZee;
      state++;
    
    } else if(state==4) {
      stringstream ss(line);
      string label;
      ss >> label >> accWmpErr >> accWmmErr >> accWmErr >> accZmmErr >> accWepErr >> accWemErr >> accWeErr >> accZeeErr;
      state++;
      errWmp  = accWmpErr*accWmpErr/accWmp/accWmp;
      errWmm  = accWmmErr*accWmmErr/accWmm/accWmm;
      errWm   = accWmErr*accWmErr/accWm/accWm;
      errZmm  = accZmmErr*accZmmErr/accZmm/accZmm;
      errWmpm = errWmp+errWmm;
      errWZm  = errWm+errZmm;
      errWep  = accWepErr*accWepErr/accWep/accWep;
      errWem  = accWemErr*accWemErr/accWem/accWem;
      errWe   = accWeErr*accWeErr/accWe/accWe;
      errZee  = accZeeErr*accZeeErr/accZee/accZee;
      errWepm = errWep+errWem;
      errWZe  = errWe+errZee;
    
    } else {
      stringstream ss(line);
      string label;
      Double_t err[12];
      ss >> label >> err[0] >> err[1] >> err[2] >> err[3] >> err[4] >> err[5] >> err[6] >> err[7] >> err[8] >> err[9] >> err[10] >> err[11];
      errWmp  += err[0] *err[0];
      errWmm  += err[1] *err[1];
      errWm   += err[2] *err[2];
      errZmm  += err[3] *err[3];
      errWmpm += err[4] *err[4];
      errWZm  += err[5] *err[5];
      errWep  += err[6] *err[6];
      errWem  += err[7] *err[7];
      errWe   += err[8] *err[8];
      errZee  += err[9] *err[9];
      errWepm += err[10]*err[10];
      errWZe  += err[11]*err[11];
    }
  }
  ifs.close();
  
  cout << endl;   
  cout << setprecision(3) << fixed; 
  
  Double_t xsWmp = nWmp/accWmp/lumi/1e6;
  cout << "    W+(mu): " << setw(8) << xsWmp  << " +/- " << setw(7) << xsWmp*nWmpErr/nWmp << " (stat) +/- " << setw(7) << sqrt(errWmp)*xsWmp << " (syst) +/- " << setw(7) << xsWmp*lumiErr << " (lumi) nb" << endl;
  
  Double_t xsWmm = nWmm/accWmm/lumi/1e6;
  cout << "    W-(mu): " << setw(8) << xsWmm  << " +/- " << setw(7) << xsWmm*nWmmErr/nWmm << " (stat) +/- " << setw(7) << sqrt(errWmm)*xsWmm << " (syst) +/- " << setw(7) << xsWmm*lumiErr << " (lumi) nb" << endl;
  
  nWm=nWmm+nWmp;
  Double_t xsWm = xsWmm+xsWmp; 
  nWmErr=sqrt(nWmmErr*nWmmErr + nWmpErr*nWmpErr);  
  cout << "     W(mu): " << setw(8) << xsWm   << " +/- " << setw(7) << xsWm*nWmErr/nWm    << " (stat) +/- " << setw(7) << sqrt(errWm)*xsWm   << " (syst) +/- " << setw(7) << xsWm*lumiErr  << " (lumi) nb" <<  endl;
  
  Double_t xsZmm  = nZmm/accZmm/lumi/1e6;
  cout << "     Z(mu): " << setw(8) << xsZmm  << " +/- " << setw(7) << xsZmm*nZmmErr/nZmm << " (stat) +/- " << setw(7) << sqrt(errZmm)*xsZmm << " (syst) +/- " << setw(7) << xsZmm*lumiErr << " (lumi) nb" << endl;
  
  cout << setprecision(3) << fixed;
  
  Double_t xsWmpm = nWmp/nWmm/accWmp*accWmm;
  cout << " W+/W-(mu): " << setw(8) << xsWmpm << " +/- " << setw(7) << sqrt(nWmpErr*nWmpErr/nWmp/nWmp+nWmmErr*nWmmErr/nWmm/nWmm)*xsWmpm << " (stat) +/- " << setw(7) << sqrt(errWmpm)*xsWmpm << " (syst)" << endl;
  
  Double_t xsWZm = xsWm/xsZmm;
  cout << "   W/Z(mu): " << setw(8) << xsWZm  << " +/- " << setw(7) << sqrt(nWmErr*nWmErr/nWm/nWm+nZmmErr*nZmmErr/nZmm/nZmm)*xsWZm      << " (stat) +/- " << setw(7) << sqrt(errWZm)*xsWZm   << " (syst)" << endl;

  cout << endl;
  cout << setprecision(3) << fixed; 
  
  Double_t xsWep = nWep/accWep/lumi/1e6;
  cout << "     W+(e): " << setw(8) << xsWep  << " +/- " << setw(7) << xsWep*nWepErr/nWep   << " (stat) +/- " << setw(7) << sqrt(errWep)*xsWep   << " (syst) +/- " << setw(7) << xsWep*lumiErr << " (lumi) nb" << endl;
  
  Double_t xsWem = nWem/accWem/lumi/1e6;
  cout << "     W-(e): " << setw(8) << xsWem  << " +/- " << setw(7) << xsWem*nWemErr/nWem   << " (stat) +/- " << setw(7) << sqrt(errWem)*xsWem   << " (syst) +/- " << setw(7) << xsWem*lumiErr << " (lumi) nb" << endl;

  nWe=nWem+nWep;
  Double_t xsWe = xsWem+xsWep; 
  nWeErr=sqrt(nWemErr*nWemErr + nWepErr*nWepErr);
  cout << "      W(e): " << setw(8) << xsWe   << " +/- " << setw(7) << xsWe*nWeErr/nWe      << " (stat) +/- " << setw(7) << sqrt(errWe)*xsWe	 << " (syst) +/- " << setw(7) << xsWe*lumiErr  << " (lumi) nb" <<  endl;
  
  Double_t xsZee = nZee/accZee/lumi/1e6;
  cout << "      Z(e): " << setw(8) << xsZee  << " +/- " << setw(7) << xsZee*nZeeErr/nZee   << " (stat) +/- " << setw(7) << sqrt(errZee)*xsZee   << " (syst) +/- " << setw(7) << xsZee*lumiErr << " (lumi) nb" << endl;
  
  cout << setprecision(3) << fixed;
  
  Double_t xsWepm = nWep/nWem/accWep*accWem;
  cout << "  W+/W-(e): " << setw(8) << xsWepm << " +/- " << setw(7) << sqrt(nWepErr*nWepErr/nWep/nWep+nWemErr*nWemErr/nWem/nWem)*xsWepm << " (stat) +/- " << setw(7) << sqrt(errWepm)*xsWepm << " (syst)" << endl;
  
  Double_t xsWZe = xsWe/xsZee;
  cout << "    W/Z(e): " << setw(8) << xsWZe  << " +/- " << setw(7) << sqrt(nWeErr*nWeErr/nWe/nWe+nZeeErr*nZeeErr/nZee/nZee)*xsWZe      << " (stat) +/- " << setw(7) << sqrt(errWZe)*xsWZe   << " (syst)" << endl;

  cout << endl;

  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  
  TCanvas *c = MakeCanvas("c","c",800,600);
  c->SetTickx(1);
  c->SetTicky(0);
  c->SetFrameFillStyle(0);
  c->SetFrameLineWidth(2);
  c->SetFrameBorderMode(0);  
  c->SetLeftMargin(0.07);
  
  gStyle->SetEndErrorSize(8);
  
  enum { kEle=0, kMu=1 };
  
  CSummaryPlot plotWplus("xsWplus",
                         "#sigma(pp#rightarrowW^{+})#timesBR(W^{+}#rightarrowl^{+}#nu) [nb]",
                         "W^{+}#rightarrowe^{+}#nu",
			 "W^{+}#rightarrow#mu^{+}#nu",
			 "W^{+}#rightarrowl^{+}#nu (combined)",
		         lumi*1e3,theory_xsWp,theory_xsWperr,0,9.3);
  plotWplus.SetResults(kEle,xsWep, xsWep*nWepErr/nWep, sqrt(errWep)*xsWep, xsWep*lumiErr);
  plotWplus.SetResults(kMu ,xsWmp, xsWmp*nWmpErr/nWmp, sqrt(errWmp)*xsWmp, xsWmp*lumiErr);
  plotWplus.Draw(c,format);
  
  CSummaryPlot plotWminus("xsWminus",
                          "#sigma(pp#rightarrowW^{-})#timesBR(W^{+}#rightarrowl^{-}#nu) [nb]",
                          "W^{-}#rightarrowe^{-}#nu",
			  "W^{-}#rightarrow#mu^{-}#nu",
			  "W^{-}#rightarrowl^{-}#nu (combined)",
		          lumi*1e3,theory_xsWm,theory_xsWmerr,0,6.5);
  plotWminus.SetResults(kEle,xsWem, xsWem*nWemErr/nWem, sqrt(errWem)*xsWem, xsWem*lumiErr);
  plotWminus.SetResults(kMu ,xsWmm, xsWmm*nWmmErr/nWmm, sqrt(errWmm)*xsWmm, xsWmm*lumiErr);
  plotWminus.Draw(c,format);
  
  CSummaryPlot plotW("xsW",
                     "#sigma(pp#rightarrowW)#timesBR(W#rightarrowl#nu) [nb]",
                     "W#rightarrowe#nu",
		     "W#rightarrow#mu#nu",
		     "W#rightarrowl#nu (combined)",
		     lumi*1e3,theory_xsW,theory_xsWerr,0,16);
  plotW.SetResults(kEle,xsWe, xsWe*nWeErr/nWe, sqrt(errWe)*xsWe, xsWe*lumiErr);
  plotW.SetResults(kMu ,xsWm, xsWm*nWmErr/nWm, sqrt(errWm)*xsWm, xsWm*lumiErr);
  plotW.Draw(c,format);
  
  CSummaryPlot plotWpm("ratioWpm",
                       "R_{+/-} = [ #sigma#timesBR ](W^{+}) / [ #sigma#timesBR ](W^{-})",
                       "W^{+}#rightarrowe^{+}#nu, W^{-}#rightarrowe^{-}#nu",
		       "W^{+}#rightarrow#mu^{+}#nu, W^{-}#rightarrow#mu^{-}#nu",
		       "W^{+}#rightarrowl^{+}#nu, W^{-}#rightarrowl^{-}#nu (combined)",
		       lumi*1e3,theory_xsWpm,theory_xsWpmerr,0,1.8);
  plotWpm.SetResults(kEle,xsWepm,sqrt(nWepErr*nWepErr/nWep/nWep+nWemErr*nWemErr/nWem/nWem)*xsWepm,sqrt(errWepm)*xsWepm);
  plotWpm.SetResults(kMu, xsWmpm,sqrt(nWmpErr*nWmpErr/nWmp/nWmp+nWmmErr*nWmmErr/nWmm/nWmm)*xsWmpm,sqrt(errWmpm)*xsWmpm);
  plotWpm.Draw(c,format);  
  
  CSummaryPlot plotZ("xsZ",
                     "#sigma(pp#rightarrowZ)#timesBR(Z#rightarrowll) [nb]",
                     "Z#rightarrowee",
		     "Z#rightarrow#mu#mu",
		     "Z#rightarrowll (combined)",
		     lumi*1e3,theory_xsZ,theory_xsZerr,0,1.45);
  plotZ.SetResults(kEle,xsZee, xsZee*nZeeErr/nZee, sqrt(errZee)*xsZee, xsZee*lumiErr);
  plotZ.SetResults(kMu ,xsZmm, xsZmm*nZmmErr/nZmm, sqrt(errZmm)*xsZmm, xsZmm*lumiErr);
  plotZ.Draw(c,format);
  
  CSummaryPlot plotWZ("ratioWZ",
                      "R_{W/Z} = [ #sigma#timesBR ](W) / [ #sigma#timesBR ](Z)",
                      "W#rightarrowe#nu, Z#rightarrowee",
		      "W#rightarrow#mu#nu, Z#rightarrow#mu#mu",
		      "W#rightarrowl#nu, Z#rightarrowll (combined)",
		      lumi*1e3,theory_xsWZ,theory_xsWZerr,0,14);
  plotWZ.SetResults(kEle,xsWZe,sqrt(nWeErr*nWeErr/nWe/nWe+nZeeErr*nZeeErr/nZee/nZee)*xsWZe,sqrt(errWZe)*xsWZe);
  plotWZ.SetResults(kMu, xsWZm,sqrt(nWmErr*nWmErr/nWm/nWm+nZmmErr*nZmmErr/nZmm/nZmm)*xsWZm,sqrt(errWZm)*xsWZm);
  plotWZ.Draw(c,format);  
}
