#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <TGraph.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <sstream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TBox.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TLatex.h>
#include <TLorentzVector.h>     // 4-vector class

#include "../Utils/MitStyleRemix.cc"
#include "CorrPlot.hh"

using namespace std;

#endif

enum {CTEQ=1, NNPDF, MSTW};

void comPM(TString p_list, TString m_list, Double_t p_xsec, Double_t m_xsec, std::vector<Double_t> &combined);

void rPlot() {

  // theory
  Double_t we_xs=19.6976; Double_t wm_xs=19.6796; 
  Double_t ze_xs=1.86766; Double_t zm_xs=1.86766;

  //measured cross sections
  Double_t we_xs_meas = 20.072; Double_t wm_xs_meas = 19.867;
  Double_t rw_meas = 1.0/(we_xs_meas/wm_xs_meas);
  Double_t ze_xs_meas = 1.919; Double_t zm_xs_meas = 1.898;
  Double_t rz_meas = 1.0/(ze_xs_meas/zm_xs_meas);

  //lepton reco+id uncertainties
  Double_t we_lep = 1.589151*0.01*we_xs_meas;
  Double_t ze_lep = 2.370844*0.01*ze_xs;

  //background subtraction/modeling
  Double_t we_bkg = 1.372443*0.01*we_xs_meas;
  Double_t ze_bkg = 0.64*0.01*ze_xs;
  
  //theory
  Double_t we_theo = 1.390827*0.01*we_xs_meas;
  Double_t ze_theo = 1.609658*0.01*ze_xs;

   //theory
  Double_t we_ewk = 0.534883*0.01*wm_xs_meas;
  Double_t ze_ewk = 0.813941*0.01*zm_xs_meas;
  
  //lepton reco+id uncertainties
  Double_t wm_lep = 1.676842*0.01*wm_xs_meas;
  Double_t zm_lep = 2.075331*0.01*zm_xs_meas;

  //background subtraction/modeling
  Double_t wm_bkg = 0.631506*0.01*wm_xs_meas;
  Double_t zm_bkg = 0.55*0.01*zm_xs_meas;

  //theory
  Double_t wm_theo = 1.316283*0.01*wm_xs_meas;
  Double_t zm_theo = 1.523023*0.01*zm_xs_meas;

  //theory
  Double_t wm_ewk = 0.534883*0.01*wm_xs_meas;
  Double_t zm_ewk = 0.813941*0.01*zm_xs_meas;


  //sys
  Double_t we_sys = sqrt(we_bkg*we_bkg+we_theo*we_theo); 
  Double_t wm_sys = sqrt(wm_bkg*wm_bkg+wm_theo*wm_theo); 
  //Double_t we_sys = sqrt(we_bkg*we_bkg+we_ewk*we_ewk); 
  //Double_t wm_sys = sqrt(wm_bkg*wm_bkg+wm_ewk*wm_ewk); 
  Double_t rw_unc_sys = rw_meas*sqrt(we_sys*we_sys/we_xs_meas/we_xs_meas + wm_sys*wm_sys/wm_xs_meas/wm_xs_meas);
  Double_t ze_sys = sqrt(ze_bkg*ze_bkg+ze_theo*ze_theo); 
  Double_t zm_sys = sqrt(zm_bkg*zm_bkg+zm_theo*zm_theo); 
  //Double_t ze_sys = sqrt(ze_bkg*ze_bkg+ze_ewk*ze_ewk); 
  //Double_t zm_sys = sqrt(zm_bkg*zm_bkg+zm_ewk*zm_ewk); 
  Double_t rz_unc_sys = rz_meas*sqrt(ze_sys*ze_sys/ze_xs_meas/ze_xs_meas + zm_sys*zm_sys/zm_xs_meas/zm_xs_meas);

  Double_t rw_unc_lep = rw_meas*sqrt(we_lep*we_lep/we_xs_meas/we_xs_meas + wm_lep*wm_lep/wm_xs_meas/wm_xs_meas);
  Double_t rz_unc_lep = rz_meas*sqrt(ze_lep*ze_lep/ze_xs_meas/ze_xs_meas + zm_lep*zm_lep/zm_xs_meas/zm_xs_meas);

  //Double_t rw_unc_theo = rw_meas*sqrt(we_theo*we_theo/we_xs_meas/we_xs_meas + wm_lep*wm_lep/wm_xs_meas/wm_xs_meas);
  //Double_t rz_unc_theo = rz_meas*sqrt(ze_lep*ze_lep/ze_xs_meas/ze_xs_meas + zm_lep*zm_lep/zm_xs_meas/zm_xs_meas);

  
  //lepton reco+id uncertainties
  //Double_t we_lep = 3.2*0.01; Double_t wm_lep = 1.0*0.01;
  //Double_t rw_unc_lep = rw_meas*(we_lep*we_lep+wm_lep*wm_lep);
  //Double_t ze_lep = 6.2*0.01; Double_t zm_lep = 1.9*0.01;
  //Double_t rz_unc_lep = rz_meas*(ze_lep*we_lep+zm_lep*wm_lep);

  //background subtraction/modeling
  //Double_t we_bkg = 1.1*0.01; Double_t wm_bkg = 0.1*0.01;
  //Double_t rw_unc_bkg = rw_meas*(we_bkg*we_bkg+wm_bkg*wm_bkg);
  //Double_t ze_bkg = 0.1*0.01; Double_t zm_bkg = 0.1*0.01;
  //Double_t rz_unc_bkg = rz_meas*(ze_bkg*we_bkg+zm_bkg*wm_bkg);

  //statistical
  Double_t we_stat = 0.124; Double_t wm_stat = 0.082;
  Double_t rw_unc_stat = rw_meas*sqrt(we_stat*we_stat/we_xs_meas/we_xs_meas + wm_stat*wm_stat/wm_xs_meas/wm_xs_meas);
  Double_t ze_stat = 0.016; Double_t zm_stat = 0.012;
  Double_t rz_unc_stat = rz_meas*sqrt(ze_stat*ze_stat/ze_xs_meas/ze_xs_meas + zm_stat*zm_stat/zm_xs_meas/zm_xs_meas);

  //theo
  //Double_t we_theo = 2.49; Double_t wm_theo = 0.09;
  //Double_t rw_unc_theo = rw_meas*(we_theo*we_theo/we_xs_meas/we_xs_meas + wm_theo*wm_theo/wm_xs_meas/wm_xs_meas);
  //Double_t ze_theo = 0.02; Double_t zm_theo = 0.01;
  //Double_t rz_unc_theo = rz_meas*(ze_theo*ze_theo/ze_xs_meas/ze_xs_meas + zm_theo*zm_theo/zm_xs_meas/zm_xs_meas);

  //theory
  //Double_t we_theo = 1.0*0.01; Double_t wm_theo = 1.3*0.01;
  //Double_t we_theo_co = 0.0*0.01; Double_t wm_theo_co = 0.0*0.01;
  //Double_t rw_unc_theo = rw_meas*(we_theo*we_theo+wm_theo*wm_theo);
  //Double_t ze_theo = 1.4*0.01; Double_t zm_theo = 1.1*0.01;
  //Double_t ze_theo_co = 0.0*0.01; Double_t zm_theo_co = 0.0*0.01;
  //Double_t rz_unc_theo = rz_meas*(ze_theo*we_theo+zm_theo*wm_theo);

  //data
  std::vector<Double_t> wUncert;
  std::vector<Double_t> zUncert;
  //wUncert.push_back(0.043); zUncert.push_back(0);
  //wUncert.push_back(0);     zUncert.push_back(0.072);
  wUncert.push_back(rw_unc_stat); zUncert.push_back(0);
  wUncert.push_back(0);           zUncert.push_back(rz_unc_stat);
  wUncert.push_back(rw_unc_sys); zUncert.push_back(0);
  wUncert.push_back(0);           zUncert.push_back(rz_unc_sys);
  wUncert.push_back(rw_unc_lep);           zUncert.push_back(rz_unc_lep);
  

  //wUncert.push_back(rw_unc_lep); zUncert.push_back(0);
  //wUncert.push_back(0);          zUncert.push_back(rz_unc_lep);
  //wUncert.push_back(rw_unc_bkg); zUncert.push_back(0);
  //wUncert.push_back(0);          zUncert.push_back(rz_unc_bkg);
  //wUncert.push_back(rw_unc_theo); zUncert.push_back(rz_unc_theo);
  //wUncert.push_back(0);          zUncert.push_back(rz_unc_theo);

  Double_t x_xs=rw_meas; Double_t y_xs=rz_meas;
  std::vector<Double_t> &x=wUncert; std::vector<Double_t> &y=zUncert;

  TGraph *gr = new TGraph(x.size()-1);
  TGraph *nom = new TGraph(1);
  nom->SetPoint(0,x_xs,y_xs);

  TMatrixDSym covMatrix(2); covMatrix=0;
  for (UInt_t i=0; i<x.size(); i++) {
    gr->SetPoint(i,x[i],y[i]);
    //cout << x[i] << ", " << y[i] << endl;
    Double_t dX=x[i];
    Double_t dY=y[i];
    
    covMatrix(0,0)+=dX*dX;
    covMatrix(1,1)+=dY*dY;
    covMatrix(0,1)+=dX*dY;
    covMatrix(1,0)+=dX*dY;
  }

  //cout << "matrix: " << endl;
  //cout << covMatrix(0,0) << " " << covMatrix(1,0) << endl;
  //cout << covMatrix(0,1) << " " << covMatrix(1,1) << endl;

  TVectorD eigVals(2); eigVals=TMatrixDSymEigen(covMatrix).GetEigenValues();
  TMatrixD eigVecs(2,2); eigVecs=TMatrixDSymEigen(covMatrix).GetEigenVectors();
  
  TLorentzVector vec(0,0,0,0);
  vec.SetPx(eigVecs(0,0));
  vec.SetPy(eigVecs(1,0));

  TEllipse *ell = new TEllipse(x_xs,y_xs,sqrt(eigVals(0)),sqrt(eigVals(1)),0,360,vec.Phi()/TMath::Pi()*180);

  TCanvas *c1 = MakeCanvas("c1", "", 800, 600);

  TLegend *l = new TLegend(0.2,0.2,0.45,0.5);
  l->SetShadowColor(0); l->SetLineColor(0);

  l->SetTextFont(42);
  l->SetTextSize(0.04);

  Double_t xmin=0.90 + 0.001, xmax=1.04 -0.001;
  Double_t ymin=0.92 + 0.001, ymax=1.08 -0.001;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // setGrid draws it on top of everything, but it shouldn't. It's related to some redraw
  Float_t gridX[] = { 0.90,.95,1.00,1.05 };
  Float_t gridY[] = { 0.90,.95,1.00,1.05 };

  TH2D *axis = new TH2D("axis","",1000,xmin,xmax,1000,ymin,ymax);
  axis->GetXaxis()->SetTitle("#sigma^{tot}_{W}xBR(W#rightarrow #mu#nu)/#sigma^{tot}_{W}xBR(W#rightarrow e#nu)");
  axis->GetYaxis()->SetTitle("#sigma^{tot}_{Z}xBR(Z#rightarrow #mu#mu)/#sigma^{tot}_{Z}xBR(Z#rightarrow ee)");
  axis->GetXaxis()->SetLimits(xmin,xmax);
  axis->GetYaxis()->SetRangeUser(ymin,ymax);

  axis->GetXaxis()->SetTitleOffset(1.4);
  axis->GetYaxis()->SetTitleOffset(1.4);

  axis->GetXaxis()->SetDecimals(); // 2 digi decimals
  axis->GetYaxis()->SetDecimals();

  axis->GetXaxis()->SetNdivisions(505);
  axis->GetYaxis()->SetNdivisions(505);
  axis->Draw("AXIS");

  /*for(unsigned i=0;i< sizeof(gridX)/sizeof(Float_t) ;++i)
  {
 	TLine *l = new TLine(  gridX[i], ymin, gridX[i] ,ymax);
	l->SetLineColor(kGray);
	l->SetLineStyle(3);
	l->SetLineWidth(1);
	l->Draw("L SAME");
  }

  for(unsigned i=0;i< sizeof(gridY)/sizeof(Float_t) ;++i)
  {
 	TLine *l = new TLine(  xmin, gridY[i], xmax ,gridY[i]);
	l->SetLineColor(kGray);
	l->SetLineStyle(3);
	l->SetLineWidth(1);
	l->Draw("L SAME");
	}*/

  nom->SetTitle("");
  //nom->GetXaxis()->SetTitle("#sigma^{tot}_{W}xBR(W#rightarrow e#nu)/#sigma^{tot}_{W}xBR(W#rightarrow #mu#nu)");
  //nom->GetYaxis()->SetTitle("#sigma^{tot}_{Z}xBR(Z#rightarrow ee)/#sigma^{tot}_{Z}xBR(Z#rightarrow #mu#mu)");
  //nom->Draw("");
  //nom->GetXaxis()->SetLimits(xmin,xmax);
  //nom->GetYaxis()->SetRangeUser(xmin,xmax);
  //
  // CMS
  char lumitxt[150];
  sprintf(lumitxt,"#scale[0.75]{%i pb^{-1} (13 TeV)}", 43);
  //sprintf(lumitxt,"#bf{%i pb^{-1} (13 TeV)}", 43);
  TPaveText *lumi = new TPaveText(0.65,0.93,0.97,0.99,"NDC");
  lumi->SetFillStyle(0); lumi->SetShadowColor(0); lumi->SetLineColor(0);
  //lumi->SetTextFont(62);
  lumi->AddText(lumitxt);
  TPaveText *prelim = new TPaveText(0.205,0.80,0.465,0.88,"NDC");
  prelim->SetFillStyle(0); prelim->SetShadowColor(0); prelim->SetLineColor(0);
  //prelim->SetTextFont(62);
  prelim->AddText("CMS #scale[0.75]{#bf{#it{Preliminary}}}");
  lumi->Draw();
  prelim->Draw();


  //TLatex *cms=new TLatex();
  //cms->SetNDC();
  //cms->SetTextFont(42);
  //cms->SetTextSize(0.05);
  //cms->SetTextAlign(13);
  //cms->DrawLatex(0.20,.9,"#bf{CMS},#scale[0.75]{ #it{Preliminary}}"); // one-line
  //cms->DrawLatex(0.20,.9,"#bf{CMS} {#scale[0.75]{#it{Preliminary}}}"); // two-lines

  //cms->SetTextSize(0.03);
  //cms->SetTextAlign(31);
  //cms->DrawLatex(.95,.93,"43 pb^{-1} (13 TeV)");

  //pdg
  TBox *b1 = new TBox(xmin,1.0009-0.0028,xmax,1.0009+0.0028);
  b1->SetFillColor(kGreen-7);
  b1->SetLineColor(kGreen+2);
  b1->SetLineWidth(2);
  TLine *l1 = new TLine(xmin,1.0009,xmax,1.0009);
  l1->SetLineWidth(2);
  l1->SetLineColor(kGreen+2);


  //pdg
  TColor *mylightblue = new TColor(2350,  145./255, 212./255.,254./255.);
  TColor *mydarkblue = new TColor(2351, 0, 147./255., 250./255.);

  TBox *b2 = new TBox(0.993-0.019,ymin,0.993+0.019,ymax);
  b2->SetFillColor( 2350 );
  b2->SetLineColor( 2351 );
  b2->SetLineWidth(2);

  TLine *l2 = new TLine(0.993,ymin,0.993,ymax);
  l2->SetLineWidth(2);
  l2->SetLineColor(2351);

  b2->Draw("same f");
  b1->Draw("same f");
  // line draw last
  l2->Draw("same l");
  l1->Draw("same l");

  nom->SetMarkerColor(kBlack);
  nom->SetMarkerStyle(20);
  nom->SetMarkerSize(0.8);

  ell->SetLineWidth(2);
  ell->SetLineColor(kBlack);
  ell->SetFillColor(kBlack);
  ell->SetFillStyle(3004);

  TGraph *sm = new TGraph();
  sm->SetPoint(0,1.0,1.0);
  sm->SetMarkerStyle(21);
  sm->SetMarkerColor(kRed);

  ell->Draw("same l");
  nom->Draw("same p");
  sm->Draw("same p");
  c1->RedrawAxis();

  l->AddEntry(ell, "Data","lf");
  l->AddEntry(b1, "PDG avg R_{Z}", "lf");
  l->AddEntry(b2, "PDG avg R_{W}", "lf");
  l->AddEntry(sm, "SM", "p");

  l->Draw();

  axis->Draw("AXIS X+Y+ SAME"); // Redraw
  axis->Draw("AXIS  SAME"); // Redraw

  c1->SaveAs("rat_el_mu.png");
  c1->SaveAs("rat_el_mu.pdf");

}

void comPM(TString p_list,
           TString m_list,
           Double_t p_xsec,
           Double_t m_xsec,
           std::vector<Double_t> &combined) {

  ifstream ifs;
  ifs.open(p_list);
  assert(ifs.is_open());
  string line;

  vector<Double_t>  vP;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t acc;
    ss >> acc;
    vP.push_back(acc);
  }
  ifs.close();

  ifs.open(m_list);
  assert(ifs.is_open());

  vector<Double_t> vM;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t acc;
    ss >> acc;
    vM.push_back(acc);
  }
  ifs.close();

  for (UInt_t i=0; i<vP.size(); i++) {
    Double_t c = (vP[i]*p_xsec+vM[i]*m_xsec)/(p_xsec+m_xsec);
    combined.push_back(c);
  }

}
