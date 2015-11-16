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
#include <TColor.h>

#include "../Utils/MitStyleRemix.hh"
#include "CorrPlot.hh"

using namespace std;

#endif

enum {CTEQ=1, NNPDF, MSTW, ABM, HERA};

void toVec(TString e_list, std::vector<Double_t> &e_vec);

void addPdf(CorrPlot *plot, Int_t pdf, TString label, Int_t color, Double_t x_xs, Double_t y_xs, std::vector<Double_t> &x, std::vector<Double_t> &y);

void addData(CorrPlot *plot, TString label, Int_t color, Double_t x_xs, Double_t y_xs, std::vector<Double_t> &x, std::vector<Double_t> &y, Int_t doFill);

void wpmCorrPlotFidM() {

  TString folder = "for2Dplots";
  
  Double_t wp_acc_nnpdf = 0.444214;
  Double_t wm_acc_nnpdf = 0.4586;
  Double_t wp_acc_mmht = 0.441985;
  Double_t wm_acc_mmht = 0.457923;
  Double_t wp_acc_cteq = 0.4425;
  Double_t wm_acc_cteq = 0.459134;
  Double_t wp_acc_abm = 0.442622;
  Double_t wm_acc_abm = 0.467855;
  Double_t wp_acc_hera = 0.43823;
  Double_t wm_acc_hera = 0.461555;

  //Theory
  Double_t wp_xs_nnpdf=11328.8*wp_acc_nnpdf;
  Double_t wm_xs_nnpdf=8369.09*wm_acc_nnpdf;

  cout << wp_xs_nnpdf << ", " << wm_xs_nnpdf << endl;

  Double_t wp_xs_cteq=11501.7*wp_acc_cteq;
  Double_t wm_xs_cteq=8520.06* wm_acc_cteq;

  Double_t wp_xs_mmht=11578.1*wp_acc_mmht;
  Double_t wm_xs_mmht=8588.01*wm_acc_mmht;

  Double_t wp_xs_abm=11725.7*wp_acc_abm;
  Double_t wm_xs_abm=8554.61*wm_acc_abm;

  Double_t wp_xs_hera=11775.9*wp_acc_hera;
  Double_t wm_xs_hera=8703.65*wm_acc_hera;

  //measured cross sections
  Double_t wp_xs_meas = 5043.61;
  Double_t wm_xs_meas = 3904.06;

  //statistical
  Double_t wp_stat = 24.876;
  Double_t wm_stat = 27.0574;

  //lepton reco+id uncertainties
  Double_t wp_lep = 1.784405*0.01; 
  Double_t wm_lep = 1.545671*0.01;

  //background subtraction/modelling
  Double_t wp_bkg = 0.631506*0.01;
  Double_t wm_bkg = 0.631506*0.01;

  //lumi
  Double_t wp_lumi = 242.097;
  Double_t wm_lumi = 187.567;

  CorrPlot plot("cplot","","#sigma^{fid}_{W^{-}}xBR(W^{-}#rightarrow #mu#nu) [pb]","#sigma^{fid}_{W^{+}}xBR(W^{+}#rightarrow #mu#nu) [pb]",3500,4200,4300,5400);

  //data
  std::vector<Double_t> wmUncert;
  std::vector<Double_t> wpUncert;
  wpUncert.push_back(wp_stat); wmUncert.push_back(0);
  wpUncert.push_back(0); wmUncert.push_back(wm_stat);

  wpUncert.push_back(wp_xs_meas*wp_lep); wmUncert.push_back(wm_xs_meas*wm_lep);

  wpUncert.push_back(wp_xs_meas*wp_bkg); wmUncert.push_back(0);
  wpUncert.push_back(0);                  wmUncert.push_back(wm_xs_meas*wm_bkg);
 
  addData(&plot, "Data (stat #oplus sys)", kViolet, wm_xs_meas, wp_xs_meas, wmUncert, wpUncert,1);

  wpUncert.push_back(wp_lumi); wmUncert.push_back(wm_lumi);

  addData(&plot, "Data (stat #oplus sys #oplus lumi)", kBlack, wm_xs_meas, wp_xs_meas, wmUncert, wpUncert,0);

  std::vector<Double_t> ct14_minus;
  toVec(folder+"/wmm_ct14.txt",
	ct14_minus);

  std::vector<Double_t> ct14_plus;
  toVec(folder+"/wpm_ct14.txt",
	ct14_plus);

  //CT14nlo
  addPdf(&plot, CTEQ, "CT14", kGreen+1, wm_xs_cteq, wp_xs_cteq, ct14_minus, ct14_plus);

  std::vector<Double_t> nnpdf23_minus;
  toVec(folder+"/wmm_nnpdf30.txt",
	nnpdf23_minus);

  std::vector<Double_t> nnpdf23_plus;
  toVec(folder+"/wpm_nnpdf30.txt",
	nnpdf23_plus);

  //NNPDF2.3nlo
  addPdf(&plot, NNPDF, "NNPDF3.0", kBlue, wm_xs_nnpdf, wp_xs_nnpdf, nnpdf23_minus, nnpdf23_plus);

  std::vector<Double_t> mstw2008_minus;
  toVec(folder+"/wmm_mmht2014.txt",
	mstw2008_minus);

  std::vector<Double_t> mstw2008_plus;
  toVec(folder+"/wpm_mmht2014.txt",
	mstw2008_plus);

  //MSTW2008
  addPdf(&plot, MSTW, "MMMHT2014", kRed, wm_xs_mmht, wp_xs_mmht, mstw2008_minus, mstw2008_plus);

  std::vector<Double_t> abm_minus;
  toVec(folder+"/wmm_abm.txt",
	abm_minus);
  
  std::vector<Double_t> abm_plus;
  toVec(folder+"/wpm_abm.txt",
	abm_plus);

  //ABM2011
  addPdf(&plot, ABM, "ABM", TColor::GetColor(248,206,104), wm_xs_abm, wp_xs_abm, abm_minus, abm_plus);

  std::vector<Double_t> hera_minus;
  toVec(folder+"/wmm_hera.txt",
	hera_minus);
  
  std::vector<Double_t> hera_plus;
  toVec(folder+"/wpm_hera.txt",
	hera_plus);

  //HERA2015
  addPdf(&plot, HERA, "HERAPDF15", kBlue+2, wm_xs_hera, wp_xs_hera, hera_minus, hera_plus);

  TCanvas *c1 = MakeCanvas("c1", "", 800, 600);
  plot.Draw(c1, "wpm_fid_mu.png",43);
}

void toVec(TString e_list, 
	   std::vector<Double_t> &e_vec) {

  ifstream ifs;
  ifs.open(e_list);
  assert(ifs.is_open());
  string line;

  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t acc;
    ss >> acc;
    e_vec.push_back(acc);
  }
  ifs.close();

  for (UInt_t i=1; i<e_vec.size(); i++) {
    e_vec[i]/=e_vec[0];
  }
  e_vec[0]=1.0;

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
  
  Double_t wp_acc_nnpdf = 0.444214;
  Double_t wm_acc_nnpdf = 0.4586;
  Double_t wp_acc_mmht = 0.441985;
  Double_t wm_acc_mmht = 0.457923;
  Double_t wp_acc_cteq = 0.4425;
  Double_t wm_acc_cteq = 0.459134;
  Double_t wp_acc_abm = 0.442622;
  Double_t wm_acc_abm = 0.467855;
  Double_t wp_acc_hera = 0.43823;
  Double_t wm_acc_hera = 0.461555;

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
    else if (pdf==ABM)
      {
	dX=x_xs*(x[0]-x[i]);
	dY=y_xs*(y[0]-y[i]);
      }
    else if (pdf==HERA)
      {
	dX=x_xs*(x[0]-x[i]);
	dY=y_xs*(y[0]-y[i]);
      }
    covMatrix(0,0)+=dX*dX;
    covMatrix(1,1)+=dY*dY;
    covMatrix(0,1)+=dX*dY;
    covMatrix(1,0)+=dX*dY;
  }

   if (pdf==MSTW) {
     covMatrix(1,1)+=182.83*182.83*wp_acc_mmht*wp_acc_mmht;
     covMatrix(0,0)+=133.45*133.45*wm_acc_mmht*wm_acc_mmht;
     covMatrix(0,1)+=182.83*133.45*wp_acc_mmht*wm_acc_mmht;
     covMatrix(1,0)+=182.83*133.45*wp_acc_mmht*wm_acc_mmht;
   }
   else if (pdf==CTEQ) {
     covMatrix(1,1)+=312.16*312.16*wp_acc_cteq*wp_acc_cteq;
     covMatrix(0,0)+=224.99*224.99*wm_acc_cteq*wm_acc_cteq;
     covMatrix(0,1)+=312.16*224.99*wp_acc_cteq*wm_acc_cteq;
     covMatrix(1,0)+=312.16*224.99*wp_acc_cteq*wm_acc_cteq;
   }
   else if (pdf==NNPDF) {
     covMatrix(1,1)+=258.56*258.56*wp_acc_nnpdf*wp_acc_nnpdf;
     covMatrix(0,0)+=195.89*195.89*wm_acc_nnpdf*wm_acc_nnpdf;
     covMatrix(0,1)+=195.89*258.56*wp_acc_nnpdf*wm_acc_nnpdf;
     covMatrix(1,0)+=195.89*258.56*wp_acc_nnpdf*wm_acc_nnpdf;
   }
   else if (pdf==ABM) {
     covMatrix(1,1)+=103.526*103.526*wp_acc_abm*wp_acc_abm;
     covMatrix(0,0)+=70.2464*70.2464*wm_acc_abm*wm_acc_abm;
     covMatrix(0,1)+=103.526*70.2464*wp_acc_abm*wm_acc_abm;
     covMatrix(1,0)+=103.526*70.2464*wp_acc_abm*wm_acc_abm;
   }
   else if (pdf==HERA) {
     covMatrix(1,1)+=248.2875*248.2875*wp_acc_hera*wp_acc_hera;
     covMatrix(0,0)+=142.66*142.66*wm_acc_hera*wm_acc_hera;
     covMatrix(0,1)+=248.2875*142.66*wp_acc_hera*wm_acc_hera;
     covMatrix(1,0)+=248.2875*142.66*wp_acc_hera*wm_acc_hera;
   }
  //cout << "matrix: " << endl;
  //cout << covMatrix(0,0) << " " << covMatrix(1,0) << endl;
  //cout << covMatrix(0,1) << " " << covMatrix(1,1) << endl;

  TVectorD eigVals(2); eigVals=TMatrixDSymEigen(covMatrix).GetEigenValues();
  TMatrixD eigVecs(2,2); eigVecs=TMatrixDSymEigen(covMatrix).GetEigenVectors();
  
  //cout << "eigen-values" << endl;
  //cout << eigVals(0) << " " << eigVals(1) << endl;
  //cout << "eigen-vectors: " << endl;
  //cout << eigVecs(0,0) << " " << eigVecs(1,0) << " " << eigVecs(0,1) << " " << eigVecs(1,1) << endl;
  
  TLorentzVector vec(0,0,0,0);
  vec.SetPx(eigVecs(0,0));
  vec.SetPy(eigVecs(1,0));

  //cout << x[0] << " " << y[0] << " " << sqrt(eigVals(0)) << " " << sqrt(eigVals(1)) << " " << vec.Phi() << " " << vec.Phi()/TMath::Pi()*180 << endl;

  TEllipse *ell = new TEllipse(x_xs*x[0],y_xs*y[0],sqrt(eigVals(0)),sqrt(eigVals(1)),0,360,vec.Phi()/TMath::Pi()*180);
  ////cout << "00000" << endl;
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
  
  //cout << "eigen-values" << endl;
  //cout << eigVals(0) << " " << eigVals(1) << endl;
  //cout << "eigen-vectors: " << endl;
  //cout << eigVecs(0,0) << " " << eigVecs(1,0) << " " << eigVecs(0,1) << " " << eigVecs(1,1) << endl;
  
  TLorentzVector vec(0,0,0,0);
  vec.SetPx(eigVecs(0,0));
  vec.SetPy(eigVecs(1,0));

  //cout << x_xs << " " << y_xs << " " << sqrt(eigVals(0)) << " " << sqrt(eigVals(1)) << " " << vec.Phi() << " " << vec.Phi()/TMath::Pi()*180 << endl;

  TEllipse *ell = new TEllipse(x_xs,y_xs,sqrt(eigVals(0)),sqrt(eigVals(1)),0,360,vec.Phi()/TMath::Pi()*180);
  if (doFill) plot->AddCorrPlot(nom, ell, label, color, kFullDotLarge, 1, 1001);
  else 
    plot->AddCorrPlot(nom, ell, label, color);
}

