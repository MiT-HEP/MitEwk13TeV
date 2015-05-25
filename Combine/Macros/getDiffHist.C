#include <iostream>

#include <TH1D.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TImage.h>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooCurve.h"

void getDiffHist(const char *inputname,TString channel) {

  // Map channels to axis labels
  vector<TString> vchannels;
  vector<TString> vxlabels;
  vector<TString> vsiglabels;

  vchannels.push_back("Zee");   vchannels.push_back("Zmm");
  vchannels.push_back("Wenu");  vchannels.push_back("Wenu_p");  vchannels.push_back("Wenu_m");
  vchannels.push_back("Wmunu"); vchannels.push_back("Wmunu_p"); vchannels.push_back("Wmunu_m");

  vxlabels.push_back("M(#e^{+}#e^{-}) [GeV/c^{2}]"); vxlabels.push_back("M(#mu^{+}#mu^{-}) [GeV/c^{2}]");
  vxlabels.push_back("#slash{E}_{T} [GeV]"); vxlabels.push_back("#slash{E}_{T} [GeV]"); vxlabels.push_back("#slash{E}_{T} [GeV]");
  vxlabels.push_back("#slash{E}_{T} [GeV]"); vxlabels.push_back("#slash{E}_{T} [GeV]"); vxlabels.push_back("#slash{E}_{T} [GeV]");

  vsiglabels.push_back("Z#rightarrowee"); vsiglabels.push_back("Z#rightarrow#mu#mu");
  vsiglabels.push_back("W#rightarrowe#nu"); vsiglabels.push_back("W^{+}#rightarrowe^{+}#nu"); vsiglabels.push_back("W^{-}#rightarrowe^{-}#nu");
  vsiglabels.push_back("W#rightarrow#mu#nu"); vsiglabels.push_back("W^{+}#rightarrow#mu^{+}#nu"); vsiglabels.push_back("W^{-}#rightarrow#mu^{-}#nu");

  TString xlabel;
  TString siglabel;
  for(int k=0; k<vchannels.size(); k++) {
    if (channel==vchannels[k]) {
      xlabel = vxlabels[k];
      siglabel = vsiglabels[k];
      break;
    }
  }

  // Read in data and fit
  TFile *inputfile = new TFile(inputname);
  RooPlot *fitplot = new RooPlot(); RooPlot *fitplot0 = new RooPlot(); RooPlot *fitplot1 = new RooPlot();
  TGraph  *data    = new TGraph();  TGraph  *data0    = new TGraph();  TGraph  *data1    = new TGraph();
  TGraph  *fit     = new TGraph();  TGraph  *fit0     = new TGraph();  TGraph  *fit1     = new TGraph();
  TGraph  *sig     = new TGraph();  TGraph  *sig0     = new TGraph();  TGraph  *sig1     = new TGraph();
  TGraph  *bkg     = new TGraph();  TGraph  *bkg0     = new TGraph();  TGraph  *bkg1     = new TGraph();
  // Access the plot
  if (channel=="Wmunu" || channel=="Wmunu_p" || channel=="Wmunu_m") {
    fitplot = (RooPlot*)inputfile->Get("sel_fit_s");   data = (TGraph*)fitplot->getHist("h_sel");
    fitplot0 = (RooPlot*)inputfile->Get("anti_fit_s"); data0 = (TGraph*)fitplot0->getHist("h_anti");
  }
  else if (channel=="Zee") {
    fitplot = (RooPlot*)inputfile->Get("NoSel_fit_s"); data = (TGraph*)fitplot->getHist("h_NoSel");
    fitplot0 = (RooPlot*)inputfile->Get("SC_fit_s");   data0 = (TGraph*)fitplot0->getHist("h_SC");
  }
  else if (channel=="Zmm") {
    fitplot = (RooPlot*)inputfile->Get("NoSel_fit_s"); data = (TGraph*)fitplot->getHist("h_NoSel");
    fitplot0 = (RooPlot*)inputfile->Get("Sta_fit_s");  data0 = (TGraph*)fitplot0->getHist("h_Sta");
    fitplot1 = (RooPlot*)inputfile->Get("Trk_fit_s");  data1 = (TGraph*)fitplot1->getHist("h_Trk");
  }
  else {
    fitplot = (RooPlot*)inputfile->Get(channel+"_fit_s");
    data = (TGraph*)fitplot->getHist("h_"+channel);
  }
  // Separate the different curves
  if (channel=="Zee") {
    fit = (TGraph*)fitplot->getCurve("pdf_binNoSel_Norm[m]");
    sig = (TGraph*)fitplot->getCurve("pdf_binNoSel_Norm[m]_Comp[shapeSig*]");
    bkg = (TGraph*)fitplot->getCurve("pdf_binNoSel_Norm[m]_Comp[shapeBkg*]");
    data->GetYaxis()->SetTitle("Events / 1.0 GeV/c^{2}");  data->SetTitle(channel+" - NoSel");
    fit0 = (TGraph*)fitplot0->getCurve("pdf_binSC_Norm[m]"); 
    sig0 = (TGraph*)fitplot0->getCurve("pdf_binSC_Norm[m]_Comp[shapeSig*]");
    bkg0 = (TGraph*)fitplot0->getCurve("pdf_binSC_Norm[m]_Comp[shapeBkg*]");
    data0->GetYaxis()->SetTitle("Events / 1.0 GeV/c^{2}"); data0->SetTitle(channel+" - SC");
  }
  else if (channel=="Zmm") {
    fit = (TGraph*)fitplot->getCurve("pdf_binNoSel_Norm[m]");
    sig = (TGraph*)fitplot->getCurve("pdf_binNoSel_Norm[m]_Comp[shapeSig*]");
    bkg = (TGraph*)fitplot->getCurve("pdf_binNoSel_Norm[m]_Comp[shapeBkg*]");
    data->GetYaxis()->SetTitle("Events / 1.0 GeV/c^{2}");  data->SetTitle(channel+" - NoSel");
    fit0 = (TGraph*)fitplot0->getCurve("pdf_binSta_Norm[m]");
    sig0 = (TGraph*)fitplot0->getCurve("pdf_binSta_Norm[m]_Comp[shapeSig*]");
    bkg0 = (TGraph*)fitplot0->getCurve("pdf_binSta_Norm[m]_Comp[shapeBkg*]");
    data0->GetYaxis()->SetTitle("Events / 1.0 GeV/c^{2}"); data0->SetTitle(channel+" - Sta");
    fit1 = (TGraph*)fitplot1->getCurve("pdf_binTrk_Norm[m]");
    sig1 = (TGraph*)fitplot1->getCurve("pdf_binTrk_Norm[m]_Comp[shapeSig*]");
    bkg1 = (TGraph*)fitplot1->getCurve("pdf_binTrk_Norm[m]_Comp[shapeBkg*]");
    data1->GetYaxis()->SetTitle("Events / 1.0 GeV/c^{2}"); data1->SetTitle(channel+" - Trk");
  }
  else if (channel=="Wmunu" || channel=="Wmunu_p" || channel=="Wmunu_m") {
    fit = (TGraph*) fitplot->getCurve("pdf_binsel_Norm[pfmet]");
    sig = (TGraph*) fitplot->getCurve("pdf_binsel_Norm[pfmet]_Comp[shapeSig*]");
    bkg = (TGraph*) fitplot->getCurve("pdf_binsel_Norm[pfmet]_Comp[shapeBkg*]");
    data->GetYaxis()->SetTitle("Events / 2.0 GeV/c^{2}"); data->SetTitle(channel+" - Sel");
    fit0 = (TGraph*) fitplot0->getCurve("pdf_binanti_Norm[pfmet]");
    sig0 = (TGraph*) fitplot0->getCurve("pdf_binanti_Norm[pfmet]_Comp[shapeSig*]");
    bkg0 = (TGraph*) fitplot0->getCurve("pdf_binanti_Norm[pfmet]_Comp[shapeBkg*]");
    data0->GetYaxis()->SetTitle("Events / 2.0 GeV/c^{2}"); data0->SetTitle(channel+" - Anti");
  }
  else {
    fit = (TGraph*) fitplot->getCurve("pdf_bin"+channel+"_Norm[pfmet]");
    sig = (TGraph*) fitplot->getCurve("pdf_bin"+TString(channel)+"_Norm[pfmet]_Comp[shapeSig*]");
    bkg = (TGraph*) fitplot->getCurve("pdf_bin"+TString(channel)+"_Norm[pfmet]_Comp[shapeBkg*]");
    data->GetYaxis()->SetTitle("Events / 2.0 GeV/c^{2}"); data->SetTitle(channel);
  }

  // Fill difference plot
//  cout << data->GetN() << "," << data->GetXaxis()->GetXmin() << "," << data->GetXaxis()->GetXmax() << endl;
  double xshift = 0;
  if (channel=="Zee" || channel=="Zmm") xshift=60;
  TH1D *hDiff  = new TH1D("hDiff","",data->GetN(),data->GetXaxis()->GetXmin()+xshift,data->GetXaxis()->GetXmax()+xshift);
  TH1D *hDiff0 = new TH1D("hDiff0","",data0->GetN(),data0->GetXaxis()->GetXmin()+xshift,data0->GetXaxis()->GetXmax()+xshift);
  TH1D *hDiff1 = new TH1D("hDiff1","",data1->GetN(),data1->GetXaxis()->GetXmin()+xshift,data1->GetXaxis()->GetXmax()+xshift);
  double x = 0;
  double y = 0;
  for (int jbin=0; jbin<data->GetN(); jbin++) {
    data->GetPoint(jbin,x,y);
    double diff = y-fit->Eval(x);
    double err = TMath::Sqrt(y);
    if(err==0) err==TMath::Sqrt(fit->Eval(x));
    if(err>0) hDiff->SetBinContent(jbin,diff/err);
    else      hDiff->SetBinContent(jbin,0);
    hDiff->SetBinError(jbin,1);
  }

  if (!(channel=="Wenu" || channel=="Wenu_p" || channel=="Wenu_m")) {
    for (int jbin=0; jbin<data0->GetN(); jbin++) {
      data0->GetPoint(jbin,x,y);
      double diff = y-fit0->Eval(x);
      double err = TMath::Sqrt(y);
      if(err==0) err==TMath::Sqrt(fit0->Eval(x));
      if(err>0) hDiff0->SetBinContent(jbin,diff/err);
      else      hDiff0->SetBinContent(jbin,0);
      hDiff0->SetBinError(jbin,1);
    }
  }

  if (channel=="Zmm") {
    for (int jbin=0; jbin<data1->GetN(); jbin++) {
      data0->GetPoint(jbin,x,y);
      double diff = y-fit1->Eval(x);
      double err = TMath::Sqrt(y);
      if(err==0) err==TMath::Sqrt(fit1->Eval(x));
      if(err>0) hDiff1->SetBinContent(jbin,diff/err);
      else      hDiff1->SetBinContent(jbin,0);
      hDiff1->SetBinError(jbin,1);
    }
  }

  hDiff->GetYaxis()->SetTitleOffset(0.42); hDiff0->GetYaxis()->SetTitleOffset(0.42); hDiff1->GetYaxis()->SetTitleOffset(0.42);
  hDiff->GetYaxis()->SetTitleSize(0.13);   hDiff0->GetYaxis()->SetTitleSize(0.13);   hDiff1->GetYaxis()->SetTitleSize(0.13);
  hDiff->GetYaxis()->SetLabelSize(0.10);   hDiff0->GetYaxis()->SetLabelSize(0.10);   hDiff1->GetYaxis()->SetLabelSize(0.10);
  hDiff->GetYaxis()->SetNdivisions(104);   hDiff0->GetYaxis()->SetNdivisions(104);   hDiff1->GetYaxis()->SetNdivisions(104);
  hDiff->GetYaxis()->CenterTitle();        hDiff0->GetYaxis()->CenterTitle();        hDiff1->GetYaxis()->CenterTitle();
  hDiff->GetXaxis()->SetTitleOffset(1.2);  hDiff0->GetXaxis()->SetTitleOffset(1.2);  hDiff1->GetXaxis()->SetTitleOffset(1.2);
  hDiff->GetXaxis()->SetTitleSize(0.13);   hDiff0->GetXaxis()->SetTitleSize(0.13);   hDiff1->GetXaxis()->SetTitleSize(0.13);
  hDiff->GetXaxis()->SetLabelSize(0.12);   hDiff0->GetXaxis()->SetLabelSize(0.12);   hDiff1->GetXaxis()->SetLabelSize(0.12);
  hDiff->GetXaxis()->CenterTitle();        hDiff0->GetXaxis()->CenterTitle();        hDiff1->GetXaxis()->CenterTitle();
  hDiff->SetStats(0);                      hDiff0->SetStats(0);                      hDiff1->SetStats(0);
  hDiff->SetMarkerStyle(7);                hDiff0->SetMarkerStyle(7);                hDiff1->SetMarkerStyle(7);
  hDiff->SetMarkerSize(50);                hDiff0->SetMarkerSize(50);                hDiff1->SetMarkerSize(50);
  hDiff->GetXaxis()->SetTitle(xlabel);     hDiff0->GetXaxis()->SetTitle(xlabel);     hDiff1->GetXaxis()->SetTitle(xlabel);
  hDiff->GetYaxis()->SetTitle("#chi");     hDiff0->GetYaxis()->SetTitle("#chi");     hDiff1->GetYaxis()->SetTitle("#chi");
  hDiff->GetYaxis()->SetRangeUser(-5,5);   hDiff0->GetYaxis()->SetRangeUser(-5,5);   hDiff1->GetYaxis()->SetRangeUser(-5,5);

  // Make plot
  TCanvas *c  = new TCanvas("c", "c", 800,800);
  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c->Divide(1,2,0,0);              c0->Divide(1,2,0,0);              c1->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0); c0->cd(1)->SetPad(0,0.3,1.0,1.0); c1->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);     c0->cd(1)->SetTopMargin(0.1);     c1->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.01); c0->cd(1)->SetBottomMargin(0.01); c1->cd(1)->SetBottomMargin(0.01);
  c->cd(1)->SetLeftMargin(0.15);   c0->cd(1)->SetLeftMargin(0.15);   c1->cd(1)->SetLeftMargin(0.15);  
  c->cd(1)->SetRightMargin(0.07);  c0->cd(1)->SetRightMargin(0.07);  c1->cd(1)->SetRightMargin(0.07);  
  c->cd(1)->SetTickx(1);           c0->cd(1)->SetTickx(1);           c1->cd(1)->SetTickx(1);
  c->cd(1)->SetTicky(1);           c0->cd(1)->SetTicky(1);           c1->cd(1)->SetTicky(1);  
  c->cd(2)->SetPad(0,0,1.0,0.3);   c0->cd(2)->SetPad(0,0,1.0,0.3);   c1->cd(2)->SetPad(0,0,1.0,0.3);
  c->cd(2)->SetTopMargin(0.05);    c0->cd(2)->SetTopMargin(0.05);    c1->cd(2)->SetTopMargin(0.05);
  c->cd(2)->SetBottomMargin(0.45); c0->cd(2)->SetBottomMargin(0.45); c1->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.15);   c0->cd(2)->SetLeftMargin(0.15);   c1->cd(2)->SetLeftMargin(0.15);
  c->cd(2)->SetRightMargin(0.07);  c0->cd(2)->SetRightMargin(0.07);  c1->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetTickx(1);           c0->cd(2)->SetTickx(1);           c1->cd(2)->SetTickx(1);
  c->cd(2)->SetTicky(1);           c0->cd(2)->SetTicky(1);           c1->cd(2)->SetTicky(1);
  c->cd(2)->SetGridy(1);           c0->cd(2)->SetGridy(1);           c1->cd(2)->SetGridy(1);

  TLegend *leg = new TLegend();
  if (channel=="Zee" || channel=="Zmm") leg = new TLegend(0.1746231,0.661362,0.5050251,0.8279053,NULL,"brNDC");
  else                                  leg = new TLegend(0.5665829,0.661362,0.8969849,0.8279053,NULL,"brNDC");
  leg->SetTextFont(62);
  leg->SetTextSize(0.03330866);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(19);
  leg->SetFillStyle(1001);
  leg->AddEntry(data,"data","lp");
  leg->AddEntry(fit,"fit","l");

  leg->AddEntry(sig,siglabel,"l");
  leg->AddEntry(bkg,"EWK+t#bar{t}+QCD","l");

  c->cd(1);
  data->Draw("ap");
  fit->Draw("pl");
  sig->Draw("pl");
  bkg->Draw("pl");
  leg->Draw("same");
  c->cd(2);
  hDiff->Draw("");
  c->Print(channel+"_bin1.png");
//  TImage *img = TImage::Create();
//  img->FromPad(c);
//  img->WriteImage(channel+".png");

  if (!(channel=="Wenu" || channel=="Wenu_p" || channel=="Wenu_m")) {
    c0->cd(1);
    data0->Draw("ap");
    fit0->Draw("pl");
    sig0->Draw("pl");
    bkg0->Draw("pl");
    leg->Draw("same");
    c0->cd(2);
    hDiff0->Draw("");
    c0->Print(channel+"_bin2.png");
  }

  if (channel=="Zmm") {
    c1->cd(1);
    data1->Draw("ap");
    fit1->Draw("pl");
    sig1->Draw("pl");
    bkg1->Draw("pl");
    leg->Draw("same");
    c1->cd(2);
    hDiff1->Draw("");
    c1->Print(channel+"_bin3.png");
  }    
}
