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
#include <TCanvas.h>
#include <TLorentzVector.h>     // 4-vector class

#include "../Utils/MitStyleRemix.hh"
#include "CorrPlot.hh"

using namespace std;

#endif

enum {CTEQ=1, NNPDF, MSTW};

void comEM(TString e_list, TString m_list, Double_t e_yield, Double_t m_yield, std::vector<Double_t> &combined);

void addPdf(CorrPlot *plot, Int_t pdf, TString label, Int_t color, Double_t x_xs, Double_t y_xs, std::vector<Double_t> &x, std::vector<Double_t> &y);

void addData(CorrPlot *plot, TString label, Int_t color, Double_t x_xs, Double_t y_xs, std::vector<Double_t> &x, std::vector<Double_t> &y, Int_t doFill);

void wpmCorrPlot() {

  TString folder = "for2Dplots";

<<<<<<< HEAD
  Double_t wme_yield=109451;
  Double_t wmm_yield=122375;
  Double_t wpe_yield=140481;
  Double_t wpm_yield=157017;
=======
  // theory
  Double_t wp_xs_nnpdf=11.33;
  Double_t wm_xs_nnpdf=8.37;

  Double_t wp_xs_cteq=11.50;
  Double_t wm_xs_cteq=8.52;

  Double_t wp_xs_mmht=11.58;
  Double_t wm_xs_mmht=8.59;

  // measured
  Double_t wpe_yield=123267;
  Double_t wpm_yield=169662;
  Double_t wme_yield=99407;
  Double_t wmm_yield=130277;

  //measured cross sections
  Double_t wp_xs_meas = 11.74;
  Double_t wm_xs_meas = 8.69;

  //statistical
  Double_t wp_stat = 0.06;
  Double_t wm_stat = 0.05;

  //lepton reco+id uncertainties
  Double_t wpe_lep = 3.2*0.01; 
  Double_t wme_lep = 3.2*0.01;
  Double_t wpm_lep = 1.0*0.01;
  Double_t wmm_lep = 1.0*0.01;

  //background subtraction/modelling
  Double_t wpe_bkg = 0.9*0.01;
  Double_t wme_bkg = 1.9*0.01;
  Double_t wpm_bkg = 0.1*0.01;
  Double_t wmm_bkg = 0.1*0.01;

  //theory
  Double_t wpe_theo = 1.5*0.01;
  Double_t wme_theo = 1.0*0.01;
  Double_t wpm_theo = 2.0*0.01;
  Double_t wmm_theo = 1.4*0.01;

  //lumi
  Double_t wp_lumi = 1.41;
  Double_t wm_lumi = 1.04;

  CorrPlot plot("cplot","","#sigma^{tot}_{W}xBR(W^{-}#rightarrow l#nu) [nb]","#sigma^{tot}_{W}xBR(W^{+}#rightarrow l#nu) [nb]",7.5,10,10,13.5);

  //data
  std::vector<Double_t> wmUncert;
  std::vector<Double_t> wpUncert;
  wpUncert.push_back(wp_stat); wmUncert.push_back(0);
  wpUncert.push_back(0); wmUncert.push_back(wm_stat);

  wpUncert.push_back(wp_xs_meas*(wpm_yield*wpm_lep)/(wpe_yield+wpm_yield)); wmUncert.push_back(wm_xs_meas*(wmm_yield*wmm_lep)/(wme_yield+wmm_yield));
  wpUncert.push_back(wp_xs_meas*(wpe_yield*wpe_lep)/(wpe_yield+wpm_yield)); wmUncert.push_back(wm_xs_meas*(wme_yield*wme_lep)/(wme_yield+wmm_yield));
>>>>>>> ca01ba804896793bb77be166f140f4f1f8e542d5

  wpUncert.push_back(wp_xs_meas*(wpm_yield*wpm_bkg)/(wpe_yield+wpm_yield)); wmUncert.push_back(0);
  wpUncert.push_back(0);                                                    wmUncert.push_back(wm_xs_meas*(wmm_yield*wmm_bkg)/(wme_yield+wmm_yield));
  wpUncert.push_back(wp_xs_meas*(wpe_yield*wpe_bkg)/(wpe_yield+wpm_yield)); wmUncert.push_back(0);
  wpUncert.push_back(0);                                                    wmUncert.push_back(wm_xs_meas*(wme_yield*wme_bkg)/(wme_yield+wmm_yield));

  wpUncert.push_back(wp_xs_meas*(wpm_yield*wpm_theo)/(wpe_yield+wpm_yield)); wmUncert.push_back(0);
  wpUncert.push_back(0);                                                     wmUncert.push_back(wm_xs_meas*(wmm_yield*wmm_theo)/(wme_yield+wmm_yield));
  wpUncert.push_back(wp_xs_meas*(wpe_yield*wpe_theo)/(wpe_yield+wpm_yield)); wmUncert.push_back(0);
  wpUncert.push_back(0);                                                     wmUncert.push_back(wm_xs_meas*(wme_yield*wme_theo)/(wme_yield+wmm_yield));

  addData(&plot, "Data", kViolet, wm_xs_meas, wp_xs_meas, wmUncert, wpUncert,1);

  wpUncert.push_back(wp_lumi); wmUncert.push_back(wm_lumi);

  addData(&plot, "Data, with Lumi", kBlack, wm_xs_meas, wp_xs_meas, wmUncert, wpUncert,0);

  std::vector<Double_t> ct14_minus;
  comEM(folder+"/wme_ct14.txt",
	folder+"/wmm_ct14.txt",
	//wme_yield,
	//wmm_yield,
	1.0,1.0,
	ct14_minus);

  std::vector<Double_t> ct14_plus;
  comEM(folder+"/wpe_ct14.txt",
	folder+"/wpm_ct14.txt",
	//wpe_yield,
	//wpm_yield,
	1.0,1.0,
	ct14_plus);

  //CT14nlo
  addPdf(&plot, CTEQ, "CT14nlo", kGreen, wm_xs_cteq, wp_xs_cteq, ct14_minus, ct14_plus);

  std::vector<Double_t> nnpdf23_minus;
  comEM(folder+"/wme_nnpdf30.txt",
	folder+"/wmm_nnpdf30.txt",
	//wme_yield,
	//wmm_yield,
	1.0,1.0,
	nnpdf23_minus);

  std::vector<Double_t> nnpdf23_plus;
  comEM(folder+"/wpe_nnpdf30.txt",
	folder+"/wpm_nnpdf30.txt",
	//wpe_yield,
	//wpm_yield,
	1.0,1.0,
	nnpdf23_plus);

  //NNPDF2.3nlo
  addPdf(&plot, NNPDF, "NNPDF3.0nlo", kBlue, wm_xs_nnpdf, wp_xs_nnpdf, nnpdf23_minus, nnpdf23_plus);

  std::vector<Double_t> mstw2008_minus;
  comEM(folder+"/wme_mmht2014.txt",
	folder+"/wmm_mmht2014.txt",
	//wme_yield,
	//wmm_yield,
	1.0,1.0,
	mstw2008_minus);

  std::vector<Double_t> mstw2008_plus;
  comEM(folder+"/wpe_mmht2014.txt",
	folder+"/wpm_mmht2014.txt",
	//wpe_yield,
	//wpm_yield,
	1.0,1.0,
	mstw2008_plus);

  //MSTW2008
  addPdf(&plot, MSTW, "MMMHT2014nlo", kRed, wm_xs_mmht, wp_xs_mmht, mstw2008_minus, mstw2008_plus);

  TCanvas *c1 = MakeCanvas("c1", "", 800, 600);
  plot.Draw(c1, "wpm_inc.png",42);

}

void comEM(TString e_list, 
	   TString m_list,
	   Double_t e_yield,
	   Double_t m_yield,
	   std::vector<Double_t> &combined) {

  ifstream ifs;
  ifs.open(e_list);
  assert(ifs.is_open());
  string line;

  vector<Double_t>  vE;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t acc;
    ss >> acc;
    vE.push_back(acc);
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

  combined.push_back(1.0);

  for (UInt_t i=1; i<vE.size(); i++) {
    Double_t c = (e_yield/(e_yield+m_yield))*(vE[i]/vE[0]) + (m_yield/(e_yield+m_yield))*(vM[i]/vM[0]);
    combined.push_back(c);
  }

}

void addPdf(CorrPlot *plot,
	    Int_t pdf,
	    TString label,
	    Int_t color,
	    Double_t x_xs,
	    Double_t y_xs,
	    std::vector<Double_t> &x,
	    std::vector<Double_t> &y) {

  TGraph *gr = new TGraph(x.size()-1);
  TGraph *nom = new TGraph(1);
  nom->SetPoint(0,x_xs*x[0],y_xs*y[0]);
  
  TMatrixDSym covMatrix(2); covMatrix=0;
  for (UInt_t i=1; i<x.size(); i++) {
    gr->SetPoint(i-1,x_xs*x[i],y_xs*y[i]);
    
    Double_t dX=0, dY=0;
    
    //for mstw08
    if (pdf==MSTW) {
      dX=x_xs*(x[0]-x[i]);
      dY=y_xs*(y[0]-y[i]);
    }
    //for ct14
    else if (pdf==CTEQ) {
      dX=x_xs*(x[0]-x[i])/1.645;
      dY=y_xs*(y[0]-y[i])/1.645;
    }
    //for nnpdf
    else if (pdf==NNPDF) {
      dX=x_xs*(x[0]-x[i])/sqrt(x.size());
      dY=y_xs*(y[0]-y[i])/sqrt(x.size());
    }

    covMatrix(0,0)+=dX*dX;
    covMatrix(1,1)+=dY*dY;
    covMatrix(0,1)+=dX*dY;
    covMatrix(1,0)+=dX*dY;
  }

  cout << "matrix: " << endl;
  cout << covMatrix(0,0) << " " << covMatrix(1,0) << endl;
  cout << covMatrix(0,1) << " " << covMatrix(1,1) << endl;

  TVectorD eigVals(2); eigVals=TMatrixDSymEigen(covMatrix).GetEigenValues();
  TMatrixD eigVecs(2,2); eigVecs=TMatrixDSymEigen(covMatrix).GetEigenVectors();
  
  cout << "eigen-values" << endl;
  cout << eigVals(0) << " " << eigVals(1) << endl;
  cout << "eigen-vectors: " << endl;
  cout << eigVecs(0,0) << " " << eigVecs(1,0) << " " << eigVecs(0,1) << " " << eigVecs(1,1) << endl;
  
  TLorentzVector vec(0,0,0,0);
  vec.SetPx(eigVecs(0,0));
  vec.SetPy(eigVecs(1,0));

  cout << x[0] << " " << y[0] << " " << sqrt(eigVals(0)) << " " << sqrt(eigVals(1)) << " " << vec.Phi() << " " << vec.Phi()/TMath::Pi()*180 << endl;

  TEllipse *ell = new TEllipse(x_xs*x[0],y_xs*y[0],sqrt(eigVals(0)),sqrt(eigVals(1)),0,360,vec.Phi()/TMath::Pi()*180);
  //cout << "00000" << endl;
  plot->AddCorrPlot(nom, ell, label, color);
}

void addData(CorrPlot *plot, 
	     TString label, 
	     Int_t color, 
	     Double_t x_xs, 
	     Double_t y_xs, 
	     std::vector<Double_t> &x, 
	     std::vector<Double_t> &y,
	     Int_t doFill) {
  
  TGraph *gr = new TGraph(x.size()-1);
  TGraph *nom = new TGraph(1);
  nom->SetPoint(0,x_xs,y_xs);

  TMatrixDSym covMatrix(2); covMatrix=0;
  for (UInt_t i=0; i<x.size(); i++) {
    gr->SetPoint(i,x[i],y[i]);
    cout << x[i] << ", " << y[i] << endl;
    Double_t dX=x[i];
    Double_t dY=y[i];
    
    covMatrix(0,0)+=dX*dX;
    covMatrix(1,1)+=dY*dY;
    covMatrix(0,1)+=dX*dY;
    covMatrix(1,0)+=dX*dY;
  }

  cout << "matrix: " << endl;
  cout << covMatrix(0,0) << " " << covMatrix(1,0) << endl;
  cout << covMatrix(0,1) << " " << covMatrix(1,1) << endl;

  TVectorD eigVals(2); eigVals=TMatrixDSymEigen(covMatrix).GetEigenValues();
  TMatrixD eigVecs(2,2); eigVecs=TMatrixDSymEigen(covMatrix).GetEigenVectors();
  
  cout << "eigen-values" << endl;
  cout << eigVals(0) << " " << eigVals(1) << endl;
  cout << "eigen-vectors: " << endl;
  cout << eigVecs(0,0) << " " << eigVecs(1,0) << " " << eigVecs(0,1) << " " << eigVecs(1,1) << endl;
  
  TLorentzVector vec(0,0,0,0);
  vec.SetPx(eigVecs(0,0));
  vec.SetPy(eigVecs(1,0));

  cout << x_xs << " " << y_xs << " " << sqrt(eigVals(0)) << " " << sqrt(eigVals(1)) << " " << vec.Phi() << " " << vec.Phi()/TMath::Pi()*180 << endl;

  TEllipse *ell = new TEllipse(x_xs,y_xs,sqrt(eigVals(0)),sqrt(eigVals(1)),0,360,vec.Phi()/TMath::Pi()*180);
  //cout << "00000" << endl;
  if (doFill) plot->AddCorrPlot(nom, ell, label, color, kFullDotLarge, 1, 1001);
  else 
    plot->AddCorrPlot(nom, ell, label, color);
}

