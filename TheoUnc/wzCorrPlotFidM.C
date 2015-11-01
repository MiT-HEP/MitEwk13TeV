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
#include <TColor.h>

using namespace std;

#endif

enum {CTEQ=1, NNPDF, MSTW, ABM, HERA};

void comPM(TString p_list, TString m_list, Double_t p_xsec, Double_t m_xsec, std::vector<Double_t> &combined);

void toVec(TString e_list, std::vector<Double_t> &e_vec);

void addPdf(CorrPlot *plot, Int_t pdf, TString label, Int_t color, Double_t x_xs, Double_t y_xs, std::vector<Double_t> &x, std::vector<Double_t> &y);

void addData(CorrPlot *plot, TString label, Int_t color, Double_t x_xs, Double_t y_xs, std::vector<Double_t> &x, std::vector<Double_t> &y, Int_t doFill);

void wzCorrPlotFidM() {

  TString folder = "for2Dplots";

  CorrPlot plot("cplot","","#sigma^{fid}_{W}xBR(W#rightarrow #mu#nu) [pb]","#sigma^{fid}_{Z}xBR(Z#rightarrow #mu#mu) [pb]",8200,9700,620,740);

  Double_t w_acc_nnpdf = 0.45032621;
  Double_t z_acc_nnpdf = 0.362546;
  Double_t w_acc_mmht = 0.44877241;
  Double_t z_acc_mmht = 0.362161;
  Double_t w_acc_cteq = 0.44957843;
  Double_t z_acc_cteq = 0.362458;
  Double_t w_acc_abm = 0.45326575;
  Double_t z_acc_abm = 0.366388;
  Double_t w_acc_hera = 0.44814294;
  Double_t z_acc_hera = 0.359996;

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

  // theory
  Double_t w_xs_nnpdf=19697.6*w_acc_nnpdf;
  Double_t z_xs_nnpdf=1867.66*z_acc_nnpdf;

  Double_t w_xs_cteq=20021.8*w_acc_cteq;
  Double_t z_xs_cteq=1897.28*z_acc_cteq;

  Double_t w_xs_mmht=20166.8*w_acc_mmht;
  Double_t z_xs_mmht=1915.66*z_acc_mmht;

  Double_t w_xs_abm=20280.3*w_acc_abm;
  Double_t z_xs_abm=1920.11*z_acc_abm;

  Double_t w_xs_hera=20479.4*w_acc_hera;
  Double_t z_xs_hera=1929.64*z_acc_hera;

  Double_t wp_xs_nnpdf=11328.8*wp_acc_nnpdf;
  Double_t wm_xs_nnpdf=8369.09*wm_acc_nnpdf;
  Double_t wp_xs_cteq=11501.7*wp_acc_cteq;
  Double_t wm_xs_cteq=8520.06* wm_acc_cteq;
  Double_t wp_xs_mmht=11578.1*wp_acc_mmht;
  Double_t wm_xs_mmht=8588.01*wm_acc_mmht;
  Double_t wp_xs_abm=11725.7*wp_acc_abm;
  Double_t wm_xs_abm=8554.61*wm_acc_abm;
  Double_t wp_xs_hera=11775.9*wp_acc_hera;
  Double_t wm_xs_hera=8703.65*wm_acc_hera;


  //measured cross sections
  Double_t w_xs_meas = 8946.63;
  Double_t z_xs_meas = 688.112;

  //statistical
  Double_t w_stat = 36.9267;
  Double_t z_stat = 4.350555;

  //lepton reco+id uncertainties
  Double_t we_lep = 1.676842*0.01;
  Double_t ze_lep = 2.075331*0.01;

  //background subtraction/modeling
  Double_t we_bkg = 0.631506*0.01;
  Double_t ze_bkg = 0.55*0.01;

  //lumi
  Double_t w_lumi = 429.611;
  Double_t z_lumi = 32.9917;

  //data
  std::vector<Double_t> wUncert;
  std::vector<Double_t> zUncert;
  wUncert.push_back(w_stat); zUncert.push_back(0);
  wUncert.push_back(0);      zUncert.push_back(z_stat);

  wUncert.push_back(w_xs_meas*we_lep); zUncert.push_back(z_xs_meas*ze_lep);

  wUncert.push_back(w_xs_meas*we_bkg); zUncert.push_back(0);
  wUncert.push_back(0);                zUncert.push_back(z_xs_meas*ze_bkg);

  addData(&plot, "Data (stat #oplus sys)", kViolet, w_xs_meas, z_xs_meas, wUncert, zUncert,1);

  wUncert.push_back(w_lumi); zUncert.push_back(z_lumi);

  addData(&plot, "Data (stat #oplus sys #oplus lumi)", kBlack, w_xs_meas, z_xs_meas, wUncert, zUncert,0);

  // CT14nlo
  std::vector<Double_t> ct14_w;
  comPM(folder+"/wpm_ct14.txt",
	folder+"/wmm_ct14.txt",
	wp_xs_cteq, wm_xs_cteq, 
	ct14_w);

  std::vector<Double_t> ct14_z;
  toVec(folder+"/zmm_ct14.txt",
	ct14_z);
 
  cout << w_xs_cteq << ", " << z_xs_cteq << ", " << ct14_w[0] << ", " << ct14_z[0] <<  endl;

  addPdf(&plot, CTEQ, "CT14", kGreen+1, w_xs_cteq, z_xs_cteq, ct14_w, ct14_z);

  // NNPDF3.0nlo
  std::vector<Double_t> nnpdf30_w;
  comPM(folder+"/wpm_nnpdf30.txt",
	folder+"/wmm_nnpdf30.txt",
	wp_xs_nnpdf, wm_xs_nnpdf, 
	nnpdf30_w);

  std::vector<Double_t> nnpdf30_z;
  toVec(folder+"/zmm_nnpdf30.txt",
	nnpdf30_z);

  addPdf(&plot, NNPDF, "NNPDF3.0", kBlue, w_xs_nnpdf, z_xs_nnpdf, nnpdf30_w, nnpdf30_z);

  // MMHT2014nlo
  std::vector<Double_t> mmht2014_w;
  comPM(folder+"/wpm_mmht2014.txt",
	folder+"/wmm_mmht2014.txt",
	wp_xs_mmht, wm_xs_mmht, 
	mmht2014_w);

  std::vector<Double_t> mmht2014_z;
  toVec(folder+"/zmm_mmht2014.txt",
	mmht2014_z);

  addPdf(&plot, MSTW, "MMHT2014", kRed, w_xs_mmht, z_xs_mmht, mmht2014_w, mmht2014_z);

  // ABM
  std::vector<Double_t> abm_w;
  comPM(folder+"/wpm_abm.txt",
	folder+"/wmm_abm.txt",
	wp_xs_abm, wm_xs_abm, 
	abm_w);

  std::vector<Double_t> abm_z;
  toVec(folder+"/zmm_abm.txt",
	abm_z);

  addPdf(&plot, ABM, "ABM", TColor::GetColor(248,206,104), w_xs_abm, z_xs_abm, abm_w, abm_z);
  
   //HERA
  std::vector<Double_t> hera_w;
  comPM(folder+"/wpm_hera.txt",
	folder+"/wmm_hera.txt",
	wp_xs_hera, wm_xs_hera, 
	hera_w);

  std::vector<Double_t> hera_z;
  toVec(folder+"/zmm_hera.txt",
	hera_z);

  addPdf(&plot, HERA, "HERAPDF15", kBlue+2, w_xs_hera, z_xs_hera, hera_w, hera_z);

  TCanvas *c1 = MakeCanvas("c1", "", 800, 600);
  plot.Draw(c1, "wz_fid_mu.png",43);

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

  Double_t w_acc_nnpdf = 0.45032621;
  Double_t z_acc_nnpdf = 0.362546;
  Double_t w_acc_mmht = 0.44877241;
  Double_t z_acc_mmht = 0.362161;
  Double_t w_acc_cteq = 0.44957843;
  Double_t z_acc_cteq = 0.362458;
  Double_t w_acc_abm = 0.45326575;
  Double_t z_acc_abm = 0.366388;
  Double_t w_acc_hera = 0.44814294;
  Double_t z_acc_hera = 0.359996;

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
     covMatrix(1,1)+=29.2049*29.2049*z_acc_mmht*z_acc_mmht;
     covMatrix(0,0)+=309.135*309.135*w_acc_mmht*w_acc_mmht;
     covMatrix(0,1)+=309.135*29.2049*w_acc_mmht*z_acc_mmht;
     covMatrix(1,0)+=309.135*29.2049*w_acc_mmht*z_acc_mmht;
   }
   else if (pdf==CTEQ) {
     covMatrix(1,1)+=47.765*47.765*z_acc_cteq*z_acc_cteq;
     covMatrix(0,0)+=312.1599*312.1599*w_acc_cteq*w_acc_cteq;
     covMatrix(0,1)+=312.1599*47.765*w_acc_cteq*z_acc_cteq;
     covMatrix(1,0)+=312.1599*47.765*w_acc_cteq*z_acc_cteq;
   }
   else if (pdf==NNPDF) {
     covMatrix(1,1)+=39.6931*39.6931*z_acc_nnpdf*z_acc_nnpdf;
     covMatrix(0,0)+=448.38*448.38*w_acc_nnpdf*w_acc_nnpdf;
     covMatrix(0,1)+=448.38*39.6931*w_acc_nnpdf*z_acc_nnpdf;
     covMatrix(1,0)+=448.38*39.6931*w_acc_nnpdf*z_acc_nnpdf;
   }
   else if (pdf==ABM) {
     covMatrix(1,1)+=15.88485*15.88485*z_acc_abm*z_acc_abm;
     covMatrix(0,0)+=172.0875*172.0875*w_acc_abm*w_acc_abm;
     covMatrix(0,1)+=172.0875*15.88485*w_acc_abm*z_acc_abm;
     covMatrix(1,0)+=172.0875*15.88485*w_acc_abm*z_acc_abm;
   }
   else if (pdf==HERA) {
     covMatrix(1,1)+=32.50225*32.50225*z_acc_hera*z_acc_hera;
     covMatrix(0,0)+=384.382*384.382*w_acc_hera*w_acc_hera;
     covMatrix(0,1)+=384.38*32.50225*w_acc_hera*z_acc_hera;
     covMatrix(1,0)+=384.38*32.50225*w_acc_hera*z_acc_hera;
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
