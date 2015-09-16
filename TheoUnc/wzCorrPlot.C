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

void comEM(TString e_list, TString m_list, Double_t e_yield, Double_t m_yield, std::vector<Double_t> &combined);

void comEM(std::vector<Double_t> &e_vec, std::vector<Double_t> &m_vec, Double_t e_yield, Double_t m_yield, std::vector<Double_t> &combined);

void addPdf(CorrPlot *plot, Int_t pdf, TString label, Int_t color, Double_t x_xs, Double_t y_xs, std::vector<Double_t> &x, std::vector<Double_t> &y);

void addData(CorrPlot *plot, TString label, Int_t color, Double_t x_xs, Double_t y_xs, std::vector<Double_t> &x, std::vector<Double_t> &y, Int_t doFill);

void wzCorrPlot() {

  TString folder = "for2Dplots";

  CorrPlot plot("cplot","","#sigma^{tot}_{W}xBR(W#rightarrow l#nu) [nb]","#sigma^{tot}_{Z}xBR(Z#rightarrow ll) [nb]",17.5,23.5,1.67,2.22);

  // theory
  Double_t w_xs_nnpdf=19.70;
  Double_t wp_xs_nnpdf=11.33;
  Double_t wm_xs_nnpdf=8.37;
  Double_t z_xs_nnpdf=1.87;

  Double_t w_xs_mmht=20.17;
  Double_t wp_xs_mmht=11.58;
  Double_t wm_xs_mmht=8.59;
  Double_t z_xs_mmht=1.92;

  Double_t w_xs_cteq=20.02;
  Double_t wp_xs_cteq=11.50;
  Double_t wm_xs_cteq=8.52;
  Double_t z_xs_cteq=1.91;

  // measured
  Double_t wpe_yield=123267;
  Double_t wpm_yield=169662;

  Double_t wme_yield=99407;
  Double_t wmm_yield=130277;

  Double_t we_yield=wpe_yield+wme_yield;
  Double_t wm_yield=wpm_yield+wmm_yield;

  Double_t ze_yield=15288;
  Double_t zm_yield=23666;

  //measured cross sections
  Double_t w_xs_meas = 20.43;
  Double_t z_xs_meas = 1.94;

  //statistical
  Double_t w_stat = 0.08;
  Double_t z_stat = 0.01;

  //lepton reco+id uncertainties
  Double_t we_lep = 0.032;
  Double_t wm_lep = 0.010;
  Double_t ze_lep = 0.062;
  Double_t zm_lep = 0.019;

  //background subtraction/modeling
  Double_t we_bkg = 0.011;
  Double_t wm_bkg = 0.001;
  Double_t ze_bkg = 0.001;
  Double_t zm_bkg = 0.001;

  //theory
  Double_t we_theo = 0.01;
  Double_t wm_theo = 0.013;
  Double_t ze_theo = 0.014;
  Double_t zm_theo = 0.011;

  //lumi
  Double_t w_lumi = 2.45;
  Double_t z_lumi = 0.23;

  //data
  std::vector<Double_t> wUncert;
  std::vector<Double_t> zUncert;
  wUncert.push_back(w_stat); zUncert.push_back(0);
  wUncert.push_back(0);      zUncert.push_back(z_stat);

  wUncert.push_back(w_xs_meas*(wm_yield*wm_lep)/(we_yield+wm_yield)); zUncert.push_back(z_xs_meas*(zm_yield*zm_lep)/(ze_yield+zm_yield));
  wUncert.push_back(w_xs_meas*(we_yield*we_lep)/(we_yield+wm_yield)); zUncert.push_back(z_xs_meas*(ze_yield*ze_lep)/(ze_yield+zm_yield));

  wUncert.push_back(w_xs_meas*(wm_yield*wm_bkg)/(we_yield+wm_yield)); zUncert.push_back(0);
  wUncert.push_back(0);                                               zUncert.push_back(z_xs_meas*(zm_yield*zm_bkg)/(ze_yield+zm_yield));
  wUncert.push_back(w_xs_meas*(we_yield*we_bkg)/(we_yield+wm_yield)); zUncert.push_back(0);
  wUncert.push_back(0);                                               zUncert.push_back(z_xs_meas*(ze_yield*ze_bkg)/(ze_yield+zm_yield));

  wUncert.push_back(w_xs_meas*(wm_yield*wm_theo)/(we_yield+wm_yield)); zUncert.push_back(0);
  wUncert.push_back(0);                                                zUncert.push_back(z_xs_meas*(zm_yield*zm_theo)/(ze_yield+zm_yield));
  wUncert.push_back(w_xs_meas*(we_yield*we_theo)/(we_yield+wm_yield)); zUncert.push_back(0);
  wUncert.push_back(0);                                                zUncert.push_back(z_xs_meas*(ze_yield*ze_theo)/(ze_yield+zm_yield));

  addData(&plot, "Data", kViolet, w_xs_meas, z_xs_meas, wUncert, zUncert,1);

  wUncert.push_back(w_lumi); zUncert.push_back(z_lumi);

  addData(&plot, "Data, with Lumi", kBlack, w_xs_meas, z_xs_meas, wUncert, zUncert,0);

  // CT14nlo
  std::vector<Double_t> ct14_we;
  comPM(folder+"/wpe_ct14.txt",
	folder+"/wme_ct14.txt",
	wp_xs_cteq, wm_xs_cteq, ct14_we);

  std::vector<Double_t> ct14_wm;
  comPM(folder+"/wpm_ct14.txt",
	folder+"/wmm_ct14.txt", 
	wp_xs_cteq, wm_xs_cteq, ct14_wm);

  std::vector<Double_t> ct14_w;
  comEM(ct14_we, ct14_wm, we_yield, wm_yield, ct14_w);

  std::vector<Double_t> ct14_z;
  comEM(folder+"/zee_ct14.txt",
	folder+"/zmm_ct14.txt", 
	ze_yield, zm_yield, ct14_z);

  addPdf(&plot, CTEQ, "CT14nlo", kGreen, w_xs_cteq, z_xs_cteq, ct14_w, ct14_z);

  // NNPDF3.0nlo
  std::vector<Double_t> nnpdf30_we;
  comPM(folder+"/wpe_nnpdf30.txt",
	folder+"/wme_nnpdf30.txt",
	wp_xs_nnpdf, wm_xs_nnpdf, nnpdf30_we);

  std::vector<Double_t> nnpdf30_wm;
  comPM(folder+"/wpm_nnpdf30.txt",
	folder+"/wmm_nnpdf30.txt", 
	wp_xs_nnpdf, wm_xs_nnpdf, nnpdf30_wm);

  std::vector<Double_t> nnpdf30_w;
  comEM(nnpdf30_we, nnpdf30_wm, we_yield, wm_yield, nnpdf30_w);

  std::vector<Double_t> nnpdf30_z;
  comEM(folder+"/zee_nnpdf30.txt",
	folder+"/zmm_nnpdf30.txt", 
	ze_yield, zm_yield, nnpdf30_z);

  addPdf(&plot, NNPDF, "NNPDF3.0nlo", kBlue, w_xs_nnpdf, z_xs_nnpdf, nnpdf30_w, nnpdf30_z);

  // MMHT2014nlo
  std::vector<Double_t> mmht2014_we;
  comPM(folder+"/wpe_mmht2014.txt",
	folder+"/wme_mmht2014.txt",
	wp_xs_mmht, wm_xs_mmht, mmht2014_we);

  std::vector<Double_t> mmht2014_wm;
  comPM(folder+"/wpm_mmht2014.txt",
	folder+"/wmm_mmht2014.txt", 
	wp_xs_mmht, wm_xs_mmht, mmht2014_wm);

  std::vector<Double_t> mmht2014_w;
  comEM(mmht2014_we, mmht2014_wm, we_yield, wm_yield, mmht2014_w);

  std::vector<Double_t> mmht2014_z;
  comEM(folder+"/zee_mmht2014.txt",
	folder+"/zmm_mmht2014.txt", 
	ze_yield, zm_yield, mmht2014_z);

  addPdf(&plot, MSTW, "MMHT2014nlo", kRed, w_xs_mmht, z_xs_mmht, mmht2014_w, mmht2014_z);

  TCanvas *c1 = MakeCanvas("c1", "", 800, 600);
  plot.Draw(c1, "wz_inc.png",42);

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

void comEM(std::vector<Double_t> &e_vec, 
	   std::vector<Double_t> &m_vec, 
	   Double_t e_yield, 
	   Double_t m_yield, 
	   std::vector<Double_t> &combined) {

  combined.push_back(1.0);

  for (UInt_t i=1; i<e_vec.size(); i++) {
    Double_t c = (e_yield/(e_yield+m_yield))*(e_vec[i]/e_vec[0]) + (m_yield/(e_yield+m_yield))*(m_vec[i]/m_vec[0]);
    combined.push_back(c);
  }

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

  cout << x[0] << " " << y[0] << " " << sqrt(eigVals(0)) << " " << sqrt(eigVals(1)) << " " << vec.Phi() << " " << vec.Phi()/TMath::\
    Pi()*180 << endl;

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
