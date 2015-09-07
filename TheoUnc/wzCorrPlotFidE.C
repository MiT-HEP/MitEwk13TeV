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

void comPM(TString p_list, TString m_list, Double_t p_xsec, Double_t m_xsec, std::vector<Double_t> &combined);

void toVec(TString e_list, std::vector<Double_t> &e_vec);

void addPdf(CorrPlot *plot, Int_t pdf, TString label, Int_t color, Double_t x_xs, Double_t y_xs, std::vector<Double_t> &x, std::vector<Double_t> &y);

void addData(CorrPlot *plot, TString label, Int_t color, Double_t x_xs, Double_t y_xs, std::vector<Double_t> &x, std::vector<Double_t> &y, Int_t doFill);

void wzCorrPlotFidE() {

  TString folder = "for2Dplots";

  CorrPlot plot("cplot","","#sigma^{acc}_{W}xBR(W#rightarrow e#nu) [nb]","#sigma^{acc}_{Z}xBR(Z#rightarrow ee) [nb]",7.5,10.2,0.55,0.77);

  // theory
  Double_t w_xs_nnpdf=19.70*0.433;
  Double_t wp_xs_nnpdf=11.33*0.427;
  Double_t wm_xs_nnpdf=8.37*0.441;
  Double_t z_xs_nnpdf=1.87*0.334;

  Double_t w_xs_mmht=20.17*0.427;
  Double_t wp_xs_mmht=11.58*0.421;
  Double_t wm_xs_mmht=8.59*0.435;
  Double_t z_xs_mmht=1.92*0.336;

  Double_t w_xs_cteq=20.02*0.428;
  Double_t wp_xs_cteq=11.50*0.422;
  Double_t wm_xs_cteq=8.52*0.437;
  Double_t z_xs_cteq=1.91*0.336;

  // measured
  Double_t we_yield=222791;
  Double_t wpe_yield=123342;
  Double_t wme_yield=99449;

  Double_t ze_yield=15387;

  //measured cross sections
  Double_t w_xs_meas = 8.88;
  Double_t z_xs_meas = 0.66;

  //statistical
  Double_t w_stat = 0.06;
  Double_t z_stat = 0.01;

  //lepton reco+id uncertainties
  Double_t we_lep = 3.2*0.01;
  Double_t ze_lep = 6.2*0.01;

  //background subtraction/modeling
  Double_t we_bkg = 1.1*0.01;
  Double_t ze_bkg = 0.1*0.01;

  //theory
  Double_t we_theo = 1.0*0.01;
  Double_t ze_theo = 1.4*0.01;

  //lumi
  Double_t w_lumi = 1.07;
  Double_t z_lumi = 0.08;

  //data
  std::vector<Double_t> wUncert;
  std::vector<Double_t> zUncert;
  wUncert.push_back(w_stat); zUncert.push_back(0);
  wUncert.push_back(0);      zUncert.push_back(z_stat);

  wUncert.push_back(w_xs_meas*we_lep); zUncert.push_back(z_xs_meas*ze_lep);

  wUncert.push_back(w_xs_meas*we_bkg); zUncert.push_back(0);
  wUncert.push_back(0);                zUncert.push_back(z_xs_meas*ze_bkg);

  addData(&plot, "Data", kViolet, w_xs_meas, z_xs_meas, wUncert, zUncert,1);

  wUncert.push_back(w_lumi); zUncert.push_back(z_lumi);

  addData(&plot, "Data, with Lumi", kBlack, w_xs_meas, z_xs_meas, wUncert, zUncert,0);

  // CT14nlo
  std::vector<Double_t> ct14_w;
  comPM(folder+"/wpe_ct14.txt",
	folder+"/wme_ct14.txt",
	wp_xs_cteq, wm_xs_cteq, 
	ct14_w);

  std::vector<Double_t> ct14_z;
  toVec(folder+"/zee_ct14.txt",
	ct14_z);

  cout << w_xs_cteq << ", " << z_xs_cteq << ", " << ct14_w[0] << ", " << ct14_z[0] <<  endl;
  addPdf(&plot, CTEQ, "CT14nlo", kGreen, w_xs_cteq, z_xs_cteq, ct14_w, ct14_z);

  // NNPDF3.0nlo
  std::vector<Double_t> nnpdf30_w;
  comPM(folder+"/wpe_nnpdf30.txt",
	folder+"/wme_nnpdf30.txt",
	wp_xs_nnpdf, wm_xs_nnpdf, 
	nnpdf30_w);

  std::vector<Double_t> nnpdf30_z;
  toVec(folder+"/zee_nnpdf30.txt",
	nnpdf30_z);

  addPdf(&plot, NNPDF, "NNPDF3.0nlo", kBlue, w_xs_nnpdf, z_xs_nnpdf, nnpdf30_w, nnpdf30_z);

  // MMHT2014nlo
  std::vector<Double_t> mmht2014_w;
  comPM(folder+"/wpe_mmht2014.txt",
	folder+"/wme_mmht2014.txt",
	wp_xs_mmht, wm_xs_mmht, 
	mmht2014_w);

  std::vector<Double_t> mmht2014_z;
  toVec(folder+"/zee_mmht2014.txt",
	mmht2014_z);

  addPdf(&plot, MSTW, "MMHT2014nlo", kRed, w_xs_mmht, z_xs_mmht, mmht2014_w, mmht2014_z);

  TCanvas *c1 = MakeCanvas("c1", "", 800, 600);
  plot.Draw(c1, "wz_fid_ele.png",42);

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

  combined.push_back(1.0);

  for (UInt_t i=1; i<vP.size(); i++) {
    Double_t c = (vP[i]*p_xsec+vM[i]*m_xsec)/(vP[0]*p_xsec+vM[0]*m_xsec);
    combined.push_back(c);
  }

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
  //cout << "00000" << endl;                                                                                                        
  if (doFill) plot->AddCorrPlot(nom, ell, label, color, kFullDotLarge, 1, 1001);
  else
    plot->AddCorrPlot(nom, ell, label, color);
}
