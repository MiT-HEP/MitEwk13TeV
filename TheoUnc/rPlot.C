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
#include <TCanvas.h>
#include <TLorentzVector.h>     // 4-vector class

#include "../Utils/MitStyleRemix.hh"
#include "CorrPlot.hh"

using namespace std;

#endif

enum {CTEQ=1, NNPDF, MSTW};

void comPM(TString p_list, TString m_list, Double_t p_xsec, Double_t m_xsec, std::vector<Double_t> &combined);

void rPlot() {

  // theory
  Double_t we_xs=19.70; Double_t wm_xs=19.70; 
  Double_t ze_xs=1.87; Double_t zm_xs=1.87;

  // measured
  Double_t we_yield=222674; Double_t wm_yield=299939;
  Double_t ze_yield=15288; Double_t zm_yield=23666;

  //measured cross sections
  Double_t we_xs_meas = 20.71; Double_t wm_xs_meas = 20.38;
  Double_t rw_meas = we_xs_meas/wm_xs_meas;
  Double_t ze_xs_meas = 1.96; Double_t zm_xs_meas = 1.94;
  Double_t rz_meas = ze_xs_meas/zm_xs_meas;

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
  Double_t we_stat = 0.15; Double_t wm_stat = 0.09;
  Double_t rw_unc_stat = rw_meas*(we_stat*we_stat/we_xs_meas/we_xs_meas + wm_stat*wm_stat/wm_xs_meas/wm_xs_meas);
  Double_t ze_stat = 0.02; Double_t zm_stat = 0.01;
  Double_t rz_unc_stat = rz_meas*(ze_stat*ze_stat/ze_xs_meas/ze_xs_meas + zm_stat*zm_stat/zm_xs_meas/zm_xs_meas);

  //sys
  Double_t we_sys = 0.78; Double_t wm_sys = 0.35;
  Double_t rw_unc_sys = rw_meas*(we_sys*we_sys/we_xs_meas/we_xs_meas + wm_sys*wm_sys/wm_xs_meas/wm_xs_meas);
  Double_t ze_sys = 0.13; Double_t zm_sys = 0.05;
  Double_t rz_unc_sys = rz_meas*(ze_sys*ze_sys/ze_xs_meas/ze_xs_meas + zm_sys*zm_sys/zm_xs_meas/zm_xs_meas);

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
  wUncert.push_back(0.043); zUncert.push_back(0);
  wUncert.push_back(0);     zUncert.push_back(0.072);
  //wUncert.push_back(rw_unc_stat); zUncert.push_back(0);
  //wUncert.push_back(0);           zUncert.push_back(rz_unc_stat);
  //wUncert.push_back(rw_unc_sys); zUncert.push_back(0);
  //wUncert.push_back(0);           zUncert.push_back(rz_unc_sys);
  //wUncert.push_back(rw_unc_lep); zUncert.push_back(0);
  //wUncert.push_back(0);          zUncert.push_back(rz_unc_lep);
  //wUncert.push_back(rw_unc_bkg); zUncert.push_back(0);
  //wUncert.push_back(0);          zUncert.push_back(rz_unc_bkg);
  //wUncert.push_back(rw_unc_theo); zUncert.push_back(0);
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

  TLegend *l = new TLegend(0.2,0.2,0.5,0.5);
  l->SetShadowColor(0); l->SetLineColor(0);

  Double_t xmin=0.8, xmax=1.1;

  nom->SetTitle("");
  nom->GetXaxis()->SetTitle("#sigma^{tot}_{W}xBR(W#rightarrow e#nu)/#sigma^{tot}_{W}xBR(W#rightarrow #mu#nu)");
  nom->GetYaxis()->SetTitle("#sigma^{tot}_{Z}xBR(Z#rightarrow ee)/#sigma^{tot}_{Z}xBR(Z#rightarrow #mu#mu)");
  nom->Draw("");
  nom->GetXaxis()->SetLimits(xmin,xmax);
  nom->GetYaxis()->SetRangeUser(xmin,xmax);

  //pdg
  TBox *b1 = new TBox(xmin,0.9991-0.0024,xmax,0.9991+0.0024);
  b1->SetFillColor(kGreen-7);
  b1->SetLineColor(kGreen+1);
  TLine *l1 = new TLine(xmin,0.9991,xmax,0.9991);
  l1->SetLineWidth(3);
  l1->SetLineColor(kGreen+1);

  b1->Draw("same f");
  l1->Draw("same l");

  //pdg
  TBox *b2 = new TBox(1.0075-0.0207,xmin,1.0075+0.0207,xmax);
  b2->SetFillColor(kAzure+7);
  b2->SetLineColor(kAzure-1);
  TLine *l2 = new TLine(1.0075,xmin,1.0075,xmax);
  l2->SetLineWidth(3);
  l2->SetLineColor(kAzure-1);

  b2->Draw("same f");
  l2->Draw("same l");

  nom->SetMarkerColor(kBlack);
  nom->SetMarkerStyle(20);
  ell->SetLineWidth(3);
  ell->SetLineColor(797);
  ell->SetFillColor(798);
  ell->SetFillStyle(1001);

  TGraph *sm = new TGraph();
  sm->SetPoint(0,1.0,1.0);
  sm->SetMarkerStyle(24);

  ell->Draw("same l");
  nom->Draw("same p");
  sm->Draw("same p");
  c1->RedrawAxis();

  l->AddEntry(ell, "Data","lf");
  l->AddEntry(b1, "PDG avg R_{Z}", "lf");
  l->AddEntry(b2, "PDG avg R_{W}", "lf");
  l->AddEntry(sm, "SM", "p");

  l->Draw();

  c1->SaveAs("rat_el_mu.png");

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
