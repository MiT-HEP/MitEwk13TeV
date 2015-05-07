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

#include <../Utils/MitStyleRemix.hh>
#include <CorrPlot.hh>

using namespace std;

#endif

enum {CTEQ=1, NNPDF, MSTW};

void comEM(TString e_list, TString m_list, Double_t e_yield, Double_t m_yield, std::vector<Double_t> &combined);

void addPdf(CorrPlot *plot, Int_t pdf, TString label, Int_t color, Double_t x_xs, Double_t y_xs, std::vector<Double_t> &x, std::vector<Double_t> &y);

void wpmCorrPlot() {

  Double_t wme_yield=30857.0;
  Double_t wmm_yield=33848.0;
  Double_t wpe_yield=44206.0;
  Double_t wpm_yield=47637.0;

  Double_t wp_xs=7.12;
  Double_t wm_xs=5.06;

  CorrPlot plot("cplot","","#sigma^{tot}_{W}xBR(W^{-}#rightarrow l#nu) [nb]","#sigma^{tot}_{W}xBR(W^{+}#rightarrow l#nu) [nb]");

  std::vector<Double_t> ct10_minus;
  comEM("/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wme_parse_ct10nlo.txt",
	"/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wmm_parse_ct10nlo.txt",
	wme_yield,
	wmm_yield,
	ct10_minus);

  std::vector<Double_t> ct10_plus;
  comEM("/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wpe_parse_ct10nlo.txt",
	"/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wpm_parse_ct10nlo.txt",
	wpe_yield,
	wpm_yield,
	ct10_plus);

  //CT10nlo
  addPdf(&plot, CTEQ, "CT10nlo", kGreen, wm_xs, wp_xs, ct10_minus, ct10_plus);

  std::vector<Double_t> nnpdf23_minus;
  comEM("/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wme_corr_nnpdf23_nlo.txt",
	"/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wmm_corr_nnpdf23_nlo.txt",
	wme_yield,
	wmm_yield,
	nnpdf23_minus);

  std::vector<Double_t> nnpdf23_plus;
  comEM("/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wpe_corr_nnpdf23_nlo.txt",
	"/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wpm_corr_nnpdf23_nlo.txt",
	wpe_yield,
	wpm_yield,
	nnpdf23_plus);


  //NNPDF2.3nlo
  addPdf(&plot, NNPDF, "NNPDF2.3nlo", kBlue, wm_xs, wp_xs, nnpdf23_minus, nnpdf23_plus);

  std::vector<Double_t> mstw2008_minus;
  comEM("/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wme_parse_mstw2008nlo68cl.txt",
	"/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wmm_parse_mstw2008nlo68cl.txt",
	wme_yield,
	wmm_yield,
	mstw2008_minus);

  std::vector<Double_t> mstw2008_plus;
  comEM("/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wpe_parse_mstw2008nlo68cl.txt",
	"/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wpm_parse_mstw2008nlo68cl.txt",
	wpe_yield,
	wpm_yield,
	mstw2008_plus);

  //MSTW2008
  addPdf(&plot, MSTW, "MSTW2008nlo", kRed, wm_xs, wp_xs, mstw2008_minus, mstw2008_plus);

  TCanvas *c1 = MakeCanvas("c1", "", 800, 600);
  plot.Draw(c1, "test.png");

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
    Double_t acc, scale;
    ss >> acc >> scale;
    vE.push_back(scale);
  }
  ifs.close();

  ifs.open(m_list);
  assert(ifs.is_open());

  vector<Double_t> vM;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t acc, scale;
    ss >> acc >> scale;
    vM.push_back(scale);
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
    //for ct10
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

  /*  cout << "matrix: " << endl;
  cout << covMatrix(0,0) << " " << covMatrix(1,0) << endl;
  cout << covMatrix(0,1) << " " << covMatrix(1,1) << endl;*/

  TVectorD eigVals(2); eigVals=TMatrixDSymEigen(covMatrix).GetEigenValues();
  TMatrixD eigVecs(2,2); eigVecs=TMatrixDSymEigen(covMatrix).GetEigenVectors();
  
  /*  cout << "eigen-values" << endl;
  cout << eigVals(0) << " " << eigVals(1) << endl;
  cout << "eigen-vectors: " << endl;
  cout << eigVecs(0,0) << " " << eigVecs(1,0) << " " << eigVecs(0,1) << " " << eigVecs(1,1) << endl;*/
  
  TLorentzVector vec(0,0,0,0);
  vec.SetPx(eigVecs(0,0));
  vec.SetPy(eigVecs(1,0));

  //cout << x[0] << " " << y[0] << " " << sqrt(eigVals(0)) << " " << sqrt(eigVals(1)) << " " << vec.Phi() << " " << vec.Phi()/TMath::Pi()*180 << endl;
  TEllipse *ell = new TEllipse(x_xs*x[0],y_xs*y[0],sqrt(eigVals(0)),sqrt(eigVals(1)),0,360,vec.Phi()/TMath::Pi()*180);
  
  plot->AddCorrPlot(nom, ell, label, color);
}
