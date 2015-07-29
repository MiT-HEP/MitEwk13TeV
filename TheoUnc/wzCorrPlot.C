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

void addPdf(CorrPlot *plot, Int_t pdf, TString x1_unc, TString x2_unc, TString y_unc, TString label, Int_t color, Double_t x_xs, Double_t y_xs);

void wzCorrPlot() {

  CorrPlot plot("cplot","","#sigma^{tot}_{W}xBR(W^{-}#rightarrow l#nu) [nb]","#sigma^{tot}_{W}xBR(W^{+}#rightarrow l#nu) [nb]");
  Double_t wp_xs=7.12;
  Double_t wm_xs=5.06;

  //CT10nlo
  addPdf(&plot, 
	 CTEQ,
	 "/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wme_parse_ct10nlo.txt",
	 "/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wpe_parse_ct10nlo.txt",
	 "/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/zee_parse_ct10nlo.txt",
	 "CT10nlo", kGreen, wm_xs, wp_xs);

  //NNPDF2.3nlo
  addPdf(&plot, 
	 NNPDF,
	 "/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wme_corr_nnpdf23_nlo.txt",
	 "/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wpe_corr_nnpdf23_nlo.txt",
	 "NNPDF2.3nlo", kBlue, wm_xs, wp_xs);

  //NNPDF2.3nlo
  /*  addPdf(&plot, 
	 NNPDF,
	 "nnpdf_phil_wme.txt",
	 "nnpdf_phil_wpe.txt",
	 "NNPDF2.3nlo Phil", kBlue, wm_xs, wp_xs);*/

  //MSTW2008
  addPdf(&plot,
	 MSTW,
	 "/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wme_parse_mstw2008nlo68cl.txt",
	 "/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc/wpe_parse_mstw2008nlo68cl.txt",
	 "MSTW2008nlo", kRed, wm_xs, wp_xs);

  TCanvas *c1 = MakeCanvas("c1", "", 800, 600);
  plot.Draw(c1, "test.png");

}

void addPdf(CorrPlot *plot,
	    Int_t pdf,
	    TString x_unc,
	    TString y_unc,
	    TString label,
	    Int_t color,
	    Double_t x_xs,
	    Double_t y_xs) {

  ifstream ifs;
  ifs.open(y_unc);
  assert(ifs.is_open());
  string line;

  vector<Double_t>  vY;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t acc, scale;
    ss >> acc >> scale;
    vY.push_back(scale);
  }
  ifs.close();

  ifs.open(x_unc);
  assert(ifs.is_open());

  vector<Double_t>  vX;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t acc, scale;
    ss >> acc >> scale;
    vX.push_back(scale);
  }
  ifs.close();

  TGraph *gr = new TGraph(vX.size()-1);
  TGraph *nom = new TGraph(1);
  nom->SetPoint(0,x_xs*vX[0],y_xs*vY[0]);

  TMatrixDSym covMatrix(2); covMatrix=0;
  for (UInt_t i=1; i<vX.size(); i++) {
    gr->SetPoint(i-1,x_xs*vX[i],y_xs*vY[i]);

    Double_t dX=0, dY=0;

    //for mstw08
    if (pdf==MSTW) {
      dX=x_xs*(vX[0]-vX[i]);
      dY=y_xs*(vY[0]-vY[i]);
    }
    //for ct10
    else if (pdf==CTEQ) {
      dX=x_xs*(vX[0]-vX[i])/1.645;
      dY=y_xs*(vY[0]-vY[i])/1.645;
    }
    //for nnpdf
    else if (pdf==NNPDF) {
      dX=x_xs*(vX[0]-vX[i])/sqrt(vX.size());
      dY=y_xs*(vY[0]-vY[i])/sqrt(vX.size());
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

  //cout << vX[0] << " " << vY[0] << " " << sqrt(eigVals(0)) << " " << sqrt(eigVals(1)) << " " << vec.Phi() << " " << vec.Phi()/TMath::Pi()*180 << endl;
  TEllipse *ell = new TEllipse(x_xs*vX[0],y_xs*vY[0],sqrt(eigVals(0)),sqrt(eigVals(1)),0,360,vec.Phi()/TMath::Pi()*180);
  
  plot->AddCorrPlot(nom, ell, label, color);
}
