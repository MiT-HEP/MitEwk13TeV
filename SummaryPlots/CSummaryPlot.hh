//================================================================================================
//
// Class to for making cross section summary plots
//
//________________________________________________________________________________________________

#ifndef CSUMMARYPLOT_HH
#define CSUMMARYPLOT_HH

#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TStyle.h>                 // class to handle ROOT plotting styles
#include <TSystem.h>                // interface to OS
#include <TCanvas.h>                // canvas class
#include <TGraphErrors.h>           // graph class
#include <TPaveText.h>              // text box class
#include <TBox.h>                   // box class
#include <TLine.h>                  // line class
#include <TString.h>                // ROOT string class
#include <cassert>                  // C++ library for assert
#include <cmath>                    // C++ math library

#include "MitStyleRemix.hh"         // style settings for drawing

class CSummaryPlot
{
public:
  // Constructor
  CSummaryPlot(TString pname="name",                  // plot name (file format extension automatically appended)
               TString xstring="xlabel",              // x-axis label
               TString labelEle="Z#rightarrowee",     // electron channel label
	       TString labelMu="Z#rightarrow#mu#mu",  // muon channel label
	       TString labelComb="Z#rightarrowll",    // combined channels label
	       Double_t luminosity=0,                 // luminosity label
	       Double_t thyval=0,                     // theory predicted value
	       Double_t thyerr=1,                     // uncertainty on theory prediction
	       Double_t x0=0,                         // x-axis minimum
	       Double_t x1=9999);                     // x-axis maximum
  
  // Destructor
  ~CSummaryPlot();
  
  // Draw function
  void Draw(TCanvas *c,             // pointer to canvas
            TString format="png");  // file format
  
  // Assign results to electron or muon channel
  void SetResults(const UInt_t ichan,         // channel index (0=electron, 1=muon)
                  const Double_t value,       // measured value
		  const Double_t stat,        // statistical uncertainty
		  const Double_t syst,        // systematic uncertainty
		  const Double_t lumierr=0);  // luminosity uncertainty
  
  // Output directory for plots (will be created by the macro if it does not exist)
  static TString sOutDir;
   
protected:
  TString name, xlabel;                                 // plot name, x-axis label
  TString labele, labelm, labell;                       // channel description
  Double_t lumi;                                        // integrated lumi
  Double_t theoryval, theoryerr;                        // theory prediction
  Double_t xmin, xmax;                                  // x-axis range
  Double_t val[2], statErr[2], systErr[2], lumiErr[2];  // measured value and uncertainties per channel
};

//--------------------------------------------------------------------------------------------------
TString CSummaryPlot::sOutDir = ".";

//--------------------------------------------------------------------------------------------------
CSummaryPlot::CSummaryPlot(TString pname, TString xstring, 
                           TString labelEle, TString labelMu, TString labelComb,
			   Double_t luminosity, Double_t thyval, Double_t thyerr, Double_t x0, Double_t x1):
name(pname),
xlabel(xstring),
labele(labelEle),
labelm(labelMu),
labell(labelComb),
lumi(luminosity),
theoryval(thyval),
theoryerr(thyerr),
xmin(x0),
xmax(x1)
{
  for(Int_t i=0; i<2; i++) {
    val[i] = statErr[i] = systErr[i] = lumiErr[i] = 0;
  }
}

//--------------------------------------------------------------------------------------------------
CSummaryPlot::~CSummaryPlot(){}

//--------------------------------------------------------------------------------------------------
void CSummaryPlot::SetResults(const UInt_t ichan, const Double_t value, 
                              const Double_t stat, const Double_t syst, const Double_t lumierr)
{
  assert(ichan<2);
  
  val[ichan]     = value;
  statErr[ichan] = stat;
  systErr[ichan] = syst;
  lumiErr[ichan] = lumierr;
}

//--------------------------------------------------------------------------------------------------
void CSummaryPlot::Draw(TCanvas *c, TString format)
{
  Double_t xval, err, yval;
  
  // electron channel
  xval = val[0];
  yval = 2.8;
  err = sqrt(statErr[0]*statErr[0] + systErr[0]*systErr[0] + lumiErr[0]*lumiErr[0]);  
  TGraphErrors grEle1(1,&xval,&yval,&err,0);
  err = sqrt(statErr[0]*statErr[0] + systErr[0]*systErr[0]); 
  TGraphErrors grEle2(1,&xval,&yval,&err,0);
  err = statErr[0]; 
  TGraphErrors grEle3(1,&xval,&yval,&err,0);
  TGraphErrors grEle4(1,&xval,&yval,0,0);
  
  // muon channel
  xval = val[1];
  yval = 1.8;
  err = sqrt(statErr[1]*statErr[1] + systErr[1]*systErr[1] + lumiErr[1]*lumiErr[1]);  
  TGraphErrors grMu1(1,&xval,&yval,&err,0);
  err = sqrt(statErr[1]*statErr[1] + systErr[1]*systErr[1]); 
  TGraphErrors grMu2(1,&xval,&yval,&err,0);
  err = statErr[1]; 
  TGraphErrors grMu3(1,&xval,&yval,&err,0);
  TGraphErrors grMu4(1,&xval,&yval,0,0);
  
  // combined
  Double_t sum = 1.0/(statErr[0]*statErr[0] + systErr[0]*systErr[0])+1.0/(statErr[1]*statErr[1] + systErr[1]*systErr[1]);
  Double_t wgt[2];
  wgt[0] = 1.0/(statErr[0]*statErr[0] + systErr[0]*systErr[0])/sum;
  wgt[1] = 1.0/(statErr[1]*statErr[1] + systErr[1]*systErr[1])/sum;
  xval = val[0]*wgt[0] + val[1]*wgt[1];
  yval = 0.8;
  err = sqrt(1.0/sum + lumiErr[0]*lumiErr[0]/val[0]/val[0]*xval*xval);
  TGraphErrors grComb1(1,&xval,&yval,&err,0);
  err = 1.0/sqrt(sum); 
  TGraphErrors grComb2(1,&xval,&yval,&err,0);
  err =  sqrt(statErr[0]*statErr[0]*wgt[0]*wgt[0] + statErr[1]*statErr[1]*wgt[1]*wgt[1]);
  TGraphErrors grComb3(1,&xval,&yval,&err,0);
  TGraphErrors grComb4(1,&xval,&yval,0,0);
  
  Double_t valComb = xval;
  Double_t statErrComb = sqrt(statErr[0]*statErr[0]*wgt[0]*wgt[0] + statErr[1]*statErr[1]*wgt[1]*wgt[1]);
  Double_t systErrComb = sqrt(systErr[0]*systErr[0]*wgt[0]*wgt[0] + systErr[1]*systErr[1]*wgt[1]*wgt[1]);
  Double_t lumiErrComb = lumiErr[0]/val[0]*xval;
  
  //--------------------------------------------------------------------------------------------------------------
  // Drawing
  //============================================================================================================== 
  
  // set up axes
  grEle1.SetTitle("");
  grEle1.GetXaxis()->SetTitle(xlabel);
  grEle1.GetXaxis()->SetTitleSize(0.05);
  grEle1.GetYaxis()->SetTitle("");
  grEle1.GetYaxis()->SetRangeUser(0,5);
  grEle1.GetXaxis()->SetLimits(xmin,xmax);
  grEle1.GetXaxis()->SetNdivisions(506);
  grEle1.GetYaxis()->SetNdivisions(0);
  grEle1.Draw("AP");
  
  // theory uncertainty band
  TBox theory_box(theoryval-theoryerr,0,theoryval+theoryerr,3.5);
  theory_box.SetLineColor(796);
  theory_box.SetFillColor(796);
  theory_box.Draw();
  
  // theory prediction line
  TLine theory_line(theoryval,0,theoryval,3.5);
  theory_line.SetLineColor(kBlue);
  theory_line.SetLineStyle(1);
  theory_line.SetLineWidth(2);
  theory_line.Draw();
  
  // electron channel: stat+syst+lumi uncertainties
  grEle1.SetMarkerStyle(kFullCircle);
  grEle1.SetMarkerSize(0);
  grEle1.SetLineWidth(2);
  grEle1.SetMarkerColor(kGreen+2);
  grEle1.SetLineColor(kGreen+2);
  if(lumiErr[0]>0) grEle1.Draw("EPSAME");
  
  // electron channel: stat+syst uncertainties
  grEle2.SetMarkerStyle(kFullCircle);
  grEle2.SetMarkerSize(0);
  grEle2.SetLineWidth(2);
  grEle2.SetMarkerColor(kRed);
  grEle2.SetLineColor(kRed);
  grEle2.Draw("EPSAME");
  
  // electron channel: stat uncertainties
  grEle3.SetMarkerStyle(kFullCircle);
  grEle3.SetMarkerSize(0);
  grEle3.SetLineWidth(2);
  grEle3.SetMarkerColor(kBlack);
  grEle3.SetLineColor(kBlack);
  grEle3.Draw("EPSAME"); 
  
  // electron channel: draw marker
  grEle4.SetMarkerStyle(kFullCircle);
  grEle4.SetMarkerSize(0.9);
  grEle4.SetLineWidth(2);
  grEle4.SetMarkerColor(kBlack);
  grEle4.SetLineColor(kBlack);
  grEle4.Draw("EPSAME"); 
  
  // muon channel: stat+syst+lumi uncertainties
  grMu1.SetMarkerStyle(kFullCircle);
  grMu1.SetMarkerSize(0);
  grMu1.SetLineWidth(2);
  grMu1.SetMarkerColor(kGreen+2);
  grMu1.SetLineColor(kGreen+2);
  if(lumiErr[1]>0) grMu1.Draw("EPSAME");
  
  // muon channel: stat+syst uncertainties
  grMu2.SetMarkerStyle(kFullCircle);
  grMu2.SetMarkerSize(0);
  grMu2.SetLineWidth(2);
  grMu2.SetMarkerColor(kRed);
  grMu2.SetLineColor(kRed);
  grMu2.Draw("EPSAME");
  
  // muon channel: stat uncertainties
  grMu3.SetMarkerStyle(kFullCircle);
  grMu3.SetMarkerSize(0);
  grMu3.SetLineWidth(2);
  grMu3.SetMarkerColor(kBlack);
  grMu3.SetLineColor(kBlack);
  grMu3.Draw("EPSAME"); 
  
  // muon channel: draw marker
  grMu4.SetMarkerStyle(kFullCircle);
  grMu4.SetMarkerSize(0.9);
  grMu4.SetLineWidth(2);
  grMu4.SetMarkerColor(kBlack);
  grMu4.SetLineColor(kBlack);
  grMu4.Draw("EPSAME"); 
  
  // combined: stat+syst+lumi uncertainties
  grComb1.SetMarkerStyle(kFullCircle);
  grComb1.SetMarkerSize(0);
  grComb1.SetLineWidth(2);
  grComb1.SetMarkerColor(kGreen+2);
  grComb1.SetLineColor(kGreen+2);
  if(lumiErr[0]>0 && lumiErr[1]>0) grComb1.Draw("EPSAME");
  
  // combined: stat+syst uncertainties
  grComb2.SetMarkerStyle(kFullCircle);
  grComb2.SetMarkerSize(0);
  grComb2.SetLineWidth(2);
  grComb2.SetMarkerColor(kRed);
  grComb2.SetLineColor(kRed);
  grComb2.Draw("EPSAME");
  
  // combined: stat uncertainties
  grComb3.SetMarkerStyle(kFullCircle);
  grComb3.SetMarkerSize(0);
  grComb3.SetLineWidth(2);
  grComb3.SetMarkerColor(kBlack);
  grComb3.SetLineColor(kBlack);
  grComb3.Draw("EPSAME");
  
  // combined: draw marker
  grComb4.SetMarkerStyle(kFullCircle);
  grComb4.SetMarkerSize(0.9);
  grComb4.SetLineWidth(2);
  grComb4.SetMarkerColor(kBlack);
  grComb4.SetLineColor(kBlack);
  grComb4.Draw("EPSAME");
  
  // CMS label
  TPaveText tb1(0.05,0.93,0.31,0.99,"NDC");
  tb1.SetFillStyle(0);
  tb1.SetBorderSize(0);
  tb1.SetTextAlign(12);
  tb1.AddText("CMS Preliminary");
  tb1.Draw();
  
  char buffer[200]; 
  
  // lumi label
  sprintf(buffer,"%.1f pb^{-1}  at  #sqrt{s} = 8 TeV",lumi);
  TPaveText tb2(0.65,0.93,0.97,0.99,"NDC");
  tb2.SetFillStyle(0);
  tb2.SetBorderSize(0);
  tb2.SetTextAlign(32);
  tb2.AddText(buffer);
  tb2.Draw();
  
  // theory predictions text
  sprintf(buffer,"%.2f #pm %.2f",theoryval,theoryerr);
  if(lumiErr[0]>0) sprintf(buffer,"%s nb",buffer);
  TPaveText tb3(0.61,0.70,0.91,0.75,"NDC");
  tb3.SetFillStyle(0);
  tb3.SetBorderSize(0);
  tb3.SetTextAlign(22);
  tb3.AddText(buffer);
  tb3.Draw();  
  
  // theory source text
  TPaveText tb4(0.45,0.79,0.95,0.89,"NDC");
  tb4.SetFillStyle(0);
  tb4.SetBorderSize(0);
  tb4.SetTextAlign(12);
  tb4.AddText("NNLO, FEWZ+MSTW2008 prediction");
  tb4.AddText("[with MSTW2008 68% CL uncertainty]");
  tb4.Draw();

  // electron channel: label text
  TPaveText texte(0.10,0.58,0.50,0.64,"NDC");
  texte.SetFillStyle(0);
  texte.SetBorderSize(0);
  texte.SetTextAlign(12);
  texte.AddText(labele);
  texte.Draw(); 
  
  // electron channel: result text
  sprintf(buffer,"%.2f #pm %.2f_{stat} #pm %.2f_{syst}",val[0],statErr[0],systErr[0]);
  if(lumiErr[0]>0) sprintf(buffer,"%s #pm %.2f_{lumi} nb",buffer,lumiErr[0]);
  TPaveText resulte(0.12,0.53,0.60,0.58,"NDC");
  resulte.SetFillStyle(0);
  resulte.SetBorderSize(0);
  resulte.SetTextAlign(12);
  resulte.AddText(buffer);
  resulte.Draw(); 
  
  // muon channel: label text
  TPaveText textm(0.10,0.42,0.50,0.48,"NDC");
  textm.SetFillStyle(0);
  textm.SetBorderSize(0);
  textm.SetTextAlign(12);
  textm.AddText(labelm);
  textm.Draw(); 
  
  // muon channel: result text
  sprintf(buffer,"%.2f #pm %.2f_{stat} #pm %.2f_{syst}",val[1],statErr[1],systErr[1]);
  if(lumiErr[1]>0) sprintf(buffer,"%s #pm %.2f_{lumi} nb",buffer,lumiErr[1]);
  TPaveText resultm(0.12,0.37,0.60,0.42,"NDC");
  resultm.SetFillStyle(0);
  resultm.SetBorderSize(0);
  resultm.SetTextAlign(12);
  resultm.AddText(buffer);
  resultm.Draw(); 
  
  // combined: label text
  TPaveText text(0.10,0.26,0.50,0.32,"NDC");
  text.SetFillStyle(0);
  text.SetBorderSize(0);
  text.SetTextAlign(12);
  text.AddText(labell);
  text.Draw(); 
  
  // combined: result text
  sprintf(buffer,"%.2f #pm %.2f_{stat} #pm %.2f_{syst}",valComb,statErrComb,systErrComb);
  if(lumiErrComb>0) sprintf(buffer,"%s #pm %.2f_{lumi} nb",buffer,lumiErrComb);
  TPaveText result(0.12,0.21,0.60,0.26,"NDC");
  result.SetFillStyle(0);
  result.SetBorderSize(0);
  result.SetTextAlign(12);
  result.AddText(buffer);
  result.Draw(); 

  // save plot (create output directory if necessary)
  gSystem->mkdir(sOutDir,kTRUE);
  c->SaveAs(sOutDir+TString("/")+name+TString(".")+format);
}
#endif
