#if !defined(__CINT__) //|| defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <TRandom3.h>
#include <TGraph.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <sstream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <TChain.h>
#include <TH1.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TEllipse.h>
#include <TCanvas.h>
#include <TLorentzVector.h>     // 4-vector class

using namespace std;

#endif

void correlationPlot(TString wpe_unc="/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wpe_parse_ct10nlo.txt", 
		     TString wme_unc="/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wme_parse_ct10nlo.txt") {

  Double_t wp_xs=7.322;
  Double_t wm_xs=5.181;

  ifstream ifs;
  ifs.open(wpe_unc);
  assert(ifs.is_open());
  string line;

  vector<Double_t>  vWpe;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t acc, scale;
    ss >> acc >> scale;
    vWpe.push_back(scale);
  }
  ifs.close();

  ifs.open(wme_unc);
  assert(ifs.is_open());

  vector<Double_t>  vWme;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t acc, scale;
    ss >> acc >> scale;
    vWme.push_back(scale);
  }
  ifs.close();


  TGraph *gr = new TGraph(vWme.size()-1);
  TGraph *nom = new TGraph(1);
  nom->SetPoint(0,wm_xs*vWme[0],wp_xs*vWpe[0]);
  nom->SetMarkerStyle(8);
  nom->SetMarkerColor(kRed);
  gr->SetMarkerStyle(8);

  TMatrixDSym covMatrix(2); covMatrix=0;
  for (Int_t i=1; i<vWme.size(); i++) {
    gr->SetPoint(i-1,wm_xs*vWme[i],wp_xs*vWpe[i]);
    Double_t dX=wm_xs*(vWme[0]-vWme[i])/1.645;
    Double_t dY=wp_xs*(vWpe[0]-vWpe[i])/1.645;

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

  cout << vWme[0] << " " << vWpe[0] << " " << sqrt(eigVals(0)) << " " << sqrt(eigVals(1)) << " " << vec.Phi() << " " << vec.Phi()/TMath::Pi()*180 << endl;
  TEllipse *ell = new TEllipse(wm_xs*vWme[0],wp_xs*vWpe[0],sqrt(eigVals(0)),sqrt(eigVals(1)),0,360,vec.Phi()/TMath::Pi()*180);

  TCanvas *c1 = new TCanvas("c1", "", 600, 600);

  gr->Draw("ap");
  nom->Draw("same p");

  ell->SetLineColor(kRed);
  ell->SetLineWidth(2);
  ell->SetFillStyle(4000);
  ell->Draw("same s");

}
