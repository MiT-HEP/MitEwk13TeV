//================================================================================================
// 
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>   
#include <TGraphAsymmErrors.h>   
#include <TColor.h>             // histogram class
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include <TGaxis.h>
#include "TLorentzVector.h"           // 4-vector class

#include "../Utils/MyTools.hh"	          // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing


#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);
TGraphAsymmErrors* TH1TOTGraphAsymmErrors(TH1 *h1);
TGraphAsymmErrors* myMakeBandSymmetric(TGraphAsymmErrors* g0, TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);
void myAddtoBand(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2);
TGraphAsymmErrors* myTGraphErrorsDivide_noErrGraph2(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);

//=== MAIN MACRO ================================================================================================= 

void plotZmmStatCorrelations(const TString  outputDir,   // output directory
             const Double_t lumi         // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotZmmStatCorrelations");
  gStyle->SetTitleOffset(0.75,"Y");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  //
  // input ntuple file names
  //
vector<TFile*> file;
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStar.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRap.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1Pt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2Pt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1Eta.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2Eta.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPt.root", "OPEN"));

  
  // plot output file format
  const TString format("png");

   
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  
   
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  

  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_StatCorrelations.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");

  
  // histograms

  TH2D *ZPT_STATCORR_MATRIX=(TH2D*)(file[0]->Get("hCorr_Bayes"));
  TH2D *PHISTAR_STATCORR_MATRIX=(TH2D*)(file[1]->Get("hCorr_Bayes"));
  TH2D *ZRAP_STATCORR_MATRIX=(TH2D*)(file[2]->Get("hCorr_Bayes"));
  TH2D *LEP1PT_STATCORR_MATRIX=(TH2D*)(file[3]->Get("hCorr_Bayes"));
  TH2D *LEP2PT_STATCORR_MATRIX=(TH2D*)(file[4]->Get("hCorr_Bayes"));
  TH2D *LEP1ETA_STATCORR_MATRIX=(TH2D*)(file[5]->Get("hCorr_Bayes"));
  TH2D *LEP2ETA_STATCORR_MATRIX=(TH2D*)(file[6]->Get("hCorr_Bayes"));
  TH2D *LEPNEGPT_STATCORR_MATRIX=(TH2D*)(file[7]->Get("hCorr_Bayes"));
  TH2D *LEPPOSPT_STATCORR_MATRIX=(TH2D*)(file[8]->Get("hCorr_Bayes"));
 


  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  char xlabel[100];     // string buffer for x-axis label
  char ylabel[100];     // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000);  
  
  TCanvas *c = MakeCanvas("c","c",800,800);
  c->cd();
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.15);
  c->SetLeftMargin(0.15);  
  c->SetRightMargin(0.15);  
  c->SetTickx(1);
  c->SetTicky(1);  
  TGaxis::SetMaxDigits(3);

  gStyle->SetTitleOffset(1.4,"Y");
  gStyle->SetPalette(1);
  gStyle->SetPaintTextFormat("4.2f");

  
  
  //
  // ZPt
  //   

  sprintf(ylabel,"p_{T}^{#mu^{+}#mu^{-}} [GeV]");
  sprintf(xlabel,"p_{T}^{#mu^{+}#mu^{-}} [GeV]");
  ZPT_STATCORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  ZPT_STATCORR_MATRIX->GetZaxis()->SetRangeUser(-1,1);
  CPlot plotZmmPt("zmmPtStatCorrelations","",xlabel,ylabel);
  plotZmmPt.AddHist2D(ZPT_STATCORR_MATRIX,"COLZ");
  plotZmmPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmmPt.AddTextBox("Bayes Correlation Matrix",0.53,0.90,0.86,0.95,0);
  plotZmmPt.SetLogx();
  plotZmmPt.SetLogy();
  plotZmmPt.Draw(c,kTRUE,format);

  //
  // PhiStar
  //   

  sprintf(ylabel,"#phi_{#eta}*");
  sprintf(xlabel,"#phi_{#eta}*");
  PHISTAR_STATCORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  PHISTAR_STATCORR_MATRIX->GetZaxis()->SetRangeUser(-1,1);
  CPlot plotZmmPhiStar("zmmPhiStarStatCorrelations","",xlabel,ylabel);
  plotZmmPhiStar.AddHist2D(PHISTAR_STATCORR_MATRIX,"COLZ");
  plotZmmPhiStar.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmmPhiStar.AddTextBox("Bayes Correlation Matrix",0.53,0.90,0.86,0.95,0);
  plotZmmPhiStar.SetLogx();
  plotZmmPhiStar.SetLogy();
  plotZmmPhiStar.Draw(c,kTRUE,format);

  //
  // ZRap
  //   

  sprintf(ylabel,"|y^{#mu^{+}#mu^{-}}|");
  sprintf(xlabel,"|y^{#mu^{+}#mu^{-}}|");
  ZRAP_STATCORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  ZRAP_STATCORR_MATRIX->GetZaxis()->SetRangeUser(-1,1);
  CPlot plotZmmRap("zmmRapStatCorrelations","",xlabel,ylabel);
  plotZmmRap.AddHist2D(ZRAP_STATCORR_MATRIX,"COLZ");
  plotZmmRap.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmmRap.AddTextBox("Bayes Correlation Matrix",0.53,0.90,0.86,0.95,0);
  plotZmmRap.Draw(c,kTRUE,format);

  //
  // Lep1Pt
  //   

  sprintf(ylabel,"p_{T} (leading muon) [GeV]");
  sprintf(xlabel,"p_{T} (leading muon) [GeV]");
  LEP1PT_STATCORR_MATRIX->GetZaxis()->SetRangeUser(-1,1);
  LEP1PT_STATCORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  CPlot plotZmmLep1Pt("zmmLep1PtStatCorrelations","",xlabel,ylabel);
  plotZmmLep1Pt.AddHist2D(LEP1PT_STATCORR_MATRIX,"COLZ");
  plotZmmLep1Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmmLep1Pt.AddTextBox("Bayes Correlation Matrix",0.53,0.90,0.86,0.95,0);
  plotZmmLep1Pt.SetLogx();
  plotZmmLep1Pt.SetLogy();
  plotZmmLep1Pt.Draw(c,kTRUE,format);

  //
  // Lep2Pt
  //   

  sprintf(ylabel,"p_{T} (2nd leading muon) [GeV]");
  sprintf(xlabel,"p_{T} (2nd leading muon) [GeV]");
  LEP2PT_STATCORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  LEP2PT_STATCORR_MATRIX->GetZaxis()->SetRangeUser(-1,1);
  CPlot plotZmmLep2Pt("zmmLep2PtStatCorrelations","",xlabel,ylabel);
  plotZmmLep2Pt.AddHist2D(LEP2PT_STATCORR_MATRIX,"COLZ");
  plotZmmLep2Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmmLep2Pt.AddTextBox("Bayes Correlation Matrix",0.53,0.90,0.86,0.95,0);
  plotZmmLep2Pt.SetLogx();
  plotZmmLep2Pt.SetLogy();
  plotZmmLep2Pt.Draw(c,kTRUE,format);

  //
  // Lep1Eta
  //   

  sprintf(ylabel,"|#eta| (leading muon)");
  sprintf(xlabel,"|#eta| (leading muon)");
  LEP1ETA_STATCORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  LEP1ETA_STATCORR_MATRIX->GetZaxis()->SetRangeUser(-1,1);
  CPlot plotZmmLep1Eta("zmmLep1EtaStatCorrelations","",xlabel,ylabel);
  plotZmmLep1Eta.AddHist2D(LEP1ETA_STATCORR_MATRIX,"COLZ");
  plotZmmLep1Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmmLep1Eta.AddTextBox("Bayes Correlation Matrix",0.53,0.90,0.86,0.95,0);
  plotZmmLep1Eta.Draw(c,kTRUE,format);

  //
  // Lep2Eta
  //   

  sprintf(ylabel,"|#eta| (2nd leading muon)");
  sprintf(xlabel,"|#eta| (2nd leading muon)");
  LEP2ETA_STATCORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  LEP2ETA_STATCORR_MATRIX->GetZaxis()->SetRangeUser(-1,1);
  CPlot plotZmmLep2Eta("zmmLep2EtaStatCorrelations","",xlabel,ylabel);
  plotZmmLep2Eta.AddHist2D(LEP2ETA_STATCORR_MATRIX,"COLZ");
  plotZmmLep2Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmmLep2Eta.AddTextBox("Bayes Correlation Matrix",0.53,0.90,0.86,0.95,0);
  plotZmmLep2Eta.Draw(c,kTRUE,format);

  //
  // LepNegPt
  //   

  sprintf(ylabel,"p_{T}^{#mu^{-}} [GeV]");
  sprintf(xlabel,"p_{T}^{#mu^{-}} [GeV]");
  LEP1PT_STATCORR_MATRIX->GetZaxis()->SetRangeUser(-1,1);
  LEP1PT_STATCORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  CPlot plotZmmLepNegPt("zmmLepNegPtStatCorrelations","",xlabel,ylabel);
  plotZmmLepNegPt.AddHist2D(LEP1PT_STATCORR_MATRIX,"COLZ");
  plotZmmLepNegPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmmLepNegPt.AddTextBox("Bayes Correlation Matrix",0.53,0.90,0.86,0.95,0);
  plotZmmLepNegPt.SetLogx();
  plotZmmLepNegPt.SetLogy();
  plotZmmLepNegPt.Draw(c,kTRUE,format);

  //
  // LepPosPt
  //   

  sprintf(ylabel,"p_{T}^{#mu^{+}} [GeV]");
  sprintf(xlabel,"p_{T}^{#mu^{+}} [GeV]");
  LEP1PT_STATCORR_MATRIX->GetZaxis()->SetRangeUser(-1,1);
  LEP1PT_STATCORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  CPlot plotZmmLepPosPt("zmmLepPosPtStatCorrelations","",xlabel,ylabel);
  plotZmmLepPosPt.AddHist2D(LEP1PT_STATCORR_MATRIX,"COLZ");
  plotZmmLepPosPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmmLepPosPt.AddTextBox("Bayes Correlation Matrix",0.53,0.90,0.86,0.95,0);
  plotZmmLepPosPt.SetLogx();
  plotZmmLepPosPt.SetLogy();
  plotZmmLepPosPt.Draw(c,kTRUE,format);

  
  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  outFile->cd();
  outFile->Write();
  outFile->Close(); 


  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  cout << endl;
 
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     

  gBenchmark->Show("plotZmmStatCorrelations");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------

TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = (TH1D*)hData->Clone("hDiff");
  hDiff->SetName(name);
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff=0;
    Double_t err=0;
    if(hData->GetBinContent(ibin)!=0)
      {
	diff = hFit->GetBinContent(ibin)/hData->GetBinContent(ibin);
	err = hFit->GetBinError(ibin)/hData->GetBinContent(ibin);
      }
    hDiff->SetBinContent(ibin,diff);
    hDiff->SetBinError(ibin,err);   
  }
  
  hDiff->GetYaxis()->SetTitleOffset(0.55);
  hDiff->GetYaxis()->SetTitleSize(0.13);
  hDiff->GetYaxis()->SetLabelSize(0.10);
  hDiff->GetYaxis()->SetNdivisions(104);
  hDiff->GetYaxis()->CenterTitle();
  hDiff->GetXaxis()->SetTitleOffset(1.2);
  hDiff->GetXaxis()->SetTitleSize(0.13);
  hDiff->GetXaxis()->SetLabelSize(0.12);
  hDiff->GetXaxis()->CenterTitle();
  
  return hDiff;
}

TGraphAsymmErrors* TH1TOTGraphAsymmErrors(TH1 *h1){


 if (!h1) cout << "TH1TOTGraph: histogram not found !" << endl;

 TGraphAsymmErrors* g1= new TGraphAsymmErrors();

 Double_t x, y, exh,exl, eyh,eyl;
 for (Int_t i=0; i<h1->GetNbinsX(); i++) {
   y=h1->GetBinContent(i+1);
   eyl=h1->GetBinError(i+1);
   eyh=h1->GetBinError(i+1);
   x=h1->GetBinCenter(i+1);
   exl=h1->GetBinWidth(i+1)/2;
   exh=h1->GetBinWidth(i+1)/2;

   g1->SetPoint(i,x,y);
   g1->SetPointError(i,exl,exh,eyl,eyh);

 }
 return g1;
}

TGraphAsymmErrors* myMakeBandSymmetric(TGraphAsymmErrors* g0, TGraphAsymmErrors* g1,TGraphAsymmErrors* g2) {
   TGraphAsymmErrors* g3= new TGraphAsymmErrors();

  Double_t  x1=0., y1=0., x2=0., y2=0., y0=0, x3=0.;

  Double_t dum;
  for (Int_t i=0; i<g1->GetN(); i++) {
    g0->GetPoint(i, x1,y0);
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    if (i==g1->GetN()-1) x2=x1;
    else                 g2->GetPoint(i+1,x2,dum);

    if (i==0)            x3=x1;
    else                 g2->GetPoint(i-1,x3,dum);

    Double_t tmp=y2;
    if (y1<y2) {y2=y1; y1=tmp;}
    Double_t y3=y0;
    g3->SetPoint(i,x1,y3);

    Double_t binwl=(x1-x3)/2.;
    Double_t binwh=(x2-x1)/2.;
    if (binwl==0.)  binwl= binwh;
    if (binwh==0.)  binwh= binwl;
    g3->SetPointError(i,g0->GetErrorXlow(i),g0->GetErrorXhigh(i),(fabs(y3-y2)+fabs(y1-y3))/2,(fabs(y3-y2)+fabs(y1-y3))/2);

  }
  return g3;
}

void myAddtoBand(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2) {

  Double_t  x1=0., y1=0.,  y2=0., y0=0;
  
  if (g1->GetN()!=g2->GetN())
    cout << " graphs have not the same # of elements " <<g1->GetN()<<" "<<g2->GetN()<< endl;
  Double_t* EYhigh1 = g1-> GetEYhigh();
  Double_t* EYlow1  = g1-> GetEYlow();
  Double_t* EYhigh2 = g2-> GetEYhigh();
  Double_t* EYlow2  = g2-> GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    if (y1==0) y1=1;
    if (y2==0) y2=1;

    Double_t eyh1=0., eyl1=0.,eyh2=0., eyl2=0.;
  
    eyh1=EYhigh1[i];
    eyh2=EYhigh2[i];
    eyh2=sqrt(eyh1*eyh1+eyh2*eyh2);
    g2->SetPointEYhigh(i,eyh2);
    eyl1=EYlow1[i];
    eyl2=EYlow2[i];
    eyl2=sqrt(eyl1*eyl1+eyl2*eyl2);
    g2->SetPointEYlow (i,eyl2);
  }
  return;

}

TGraphAsymmErrors* myTGraphErrorsDivide_noErrGraph2(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2) {

  const Int_t debug=0; 

  TGraphAsymmErrors* g3= new TGraphAsymmErrors();
  Int_t n1=g1->GetN();
  Int_t n2=g2->GetN();

  if (n1!=n2) {
    printf(" vectors do not have same number of entries !  \n");
   return g3;
  }

  Double_t   x1=0.,   y1=0., x2=0., y2=0.;
  Double_t dx1h=0., dx1l=0.;
  Double_t dy1h=0., dy1l=0.;
  Double_t dy2h=0., dy2l=0.;

  Double_t* X1 = g1->GetX();
  Double_t* Y1 = g1->GetY();
  Double_t* EXhigh1 = g1->GetEXhigh();
  Double_t* EXlow1 =  g1->GetEXlow();
  Double_t* EYhigh1 = g1->GetEYhigh();
  Double_t* EYlow1 =  g1->GetEYlow();

  Double_t* X2 = g2->GetX();
  Double_t* Y2 = g2->GetY();
  Double_t* EXhigh2 = g2->GetEXhigh();
  Double_t* EXlow2 =  g2->GetEXlow();
  Double_t* EYhigh2 = g2->GetEYhigh();
  Double_t* EYlow2 =  g2->GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i,x1,y1);
    g2->GetPoint(i,x2,y2);
    dx1h  = EXhigh1[i];
    dx1l  = EXlow1[i];
    if (y1!=0.) dy1h  = EYhigh1[i]/y1;
    else        dy1h  = 0.;
    if (y2!=0.) dy2h  = EYhigh2[i]/y2;
    else        dy2h  = 0.;
    if (y1!=0.) dy1l  = EYlow1 [i]/y1;
    else        dy1l  = 0.;
    if (y2!=0.) dy2l  = EYlow2 [i]/y2;
    else        dy2l  = 0.;
    
    if (y2!=0.) g3->SetPoint(i, x1,y1/y2);
    else       g3->SetPoint(i, x1,y2);
    Double_t el=0.; Double_t eh=0.;

    if (y2!=0.) el=sqrt(dy1l*dy1l)*(y1/y2);
    if (y2!=0.) eh=sqrt(dy1h*dy1h)*(y1/y2);

    g3->SetPointError(i,dx1h,dx1l,fabs(el),fabs(eh));

  }  
  return g3;

}
