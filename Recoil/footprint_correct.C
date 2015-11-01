//================================================================================================
//
// Perform fits to recoil against Z->ee events
//
//  * Outputs a ROOT file of functions parametrizing the distribution of recoil components
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>                   // standard I/O
#include <fstream>                    // standard I/O
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TF1.h>                      // 1D function
#include <TFitResult.h>               // class to handle fit results
#include <TGraphErrors.h>             // graph class
#include "TLorentzVector.h"           // 4-vector class
#include "TVector2.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TStyle.h"

#include "../Utils/CPlot.hh"          // helper class for plots
#include "../Utils/MitStyleRemix.hh"  // style settings for drawing
//#include "../Utils/RecoilCorrector.hh"    // class to handle recoil corrections for MET
#include "../Utils/RecoilCorrector_htautau_hist.hh"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#endif

using namespace RooFit;

//=== FUNCTION DECLARATIONS ======================================================================================

// generate web page
void makeHTML(const TString outDir,  const Int_t nbins, 
              const Int_t pfu1model, const Int_t pfu2model);

TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {

    hDiff->SetBinContent(ibin,hData->GetBinContent(ibin)/hFit->GetBinContent(ibin));
    hDiff->SetBinError(ibin,0);
  }

  hDiff->GetYaxis()->SetTitleOffset(0.42);
  hDiff->GetYaxis()->SetTitleSize(0.13);
  hDiff->GetYaxis()->SetLabelSize(0.10);
  hDiff->GetYaxis()->SetNdivisions(104);
  hDiff->GetYaxis()->CenterTitle();
  hDiff->GetXaxis()->SetTitleOffset(1.2);
  hDiff->GetXaxis()->SetTitleSize(0.13);
  hDiff->GetXaxis()->SetLabelSize(0.12);
  hDiff->GetXaxis()->CenterTitle();

  return hDiff;
}


//--------------------------------------------------------------------------------------------------
// function to describe Gaussian widths as a function of dilepton pT
/*Double_t sigmaFunc(Double_t *x, Double_t *par) {
  // par[0]: quadratic coefficient
  // par[1]: linear coefficient
  // par[2]: constant term
  
  Double_t a  = par[0];
  Double_t b  = par[1];
  Double_t c  = par[2];
    
  return a*x[0]*x[0] + b*x[0] + c;
  }*/

//--------------------------------------------------------------------------------------------------
// function to describe relative fraction in a double Gaussian based on 
// functions for sigma0, sigma1, and sigma2
/*Double_t frac2Func(Double_t *x, Double_t *par) {
  // par[0..3]:  sigma0
  // par[4..7]:  sigma1
  // par[8..11]: sigma2
  
  TF1 s0("_s0",sigmaFunc,0,7000,4); s0.SetParameters(par[0],par[1],par[2],par[3]);
  TF1 s1("_s1",sigmaFunc,0,7000,4); s1.SetParameters(par[4],par[5],par[6],par[7]);
  TF1 s2("_s2",sigmaFunc,0,7000,4); s2.SetParameters(par[8],par[9],par[10],par[11]);
  
  return (s0.Eval(x[0]) - s1.Eval(x[0]))/(s2.Eval(x[0]) - s1.Eval(x[0]));
  }*/


//--------------------------------------------------------------------------------------------------
/*Double_t dMean(const TF1 *fcn, const Double_t x, const TFitResultPtr fs) {
  Double_t df[2];
  df[0] = 1;
  df[1] = x;
  Double_t err2 = df[0]*df[0]*(fs->GetCovarianceMatrix()[0][0]) 
                  + df[1]*df[1]*(fs->GetCovarianceMatrix()[1][1]) 
		  + 2.0*df[0]*df[1]*(fs->GetCovarianceMatrix()[0][1]);
  assert(err2>=0);
  return sqrt(err2);
  }*/

//--------------------------------------------------------------------------------------------------
/*Double_t dSigma(const TF1 *fcn, const Double_t x, const TFitResultPtr fs) {
  Double_t df[4];
  Double_t a  = fcn->GetParameter(0);
  Double_t b  = fcn->GetParameter(1);
  Double_t c  = fcn->GetParameter(2);
  
  df[0] = x*x;
  df[1] = x;
  df[2] = 1;
  
  Double_t err2=0;
  for(Int_t i=0; i<3; i++) {
    err2 += df[i]*df[i]*(fs->GetCovarianceMatrix()[i][i]);
    for(Int_t j=i+1; j<3; j++) {
      err2 += 2.0*df[i]*df[j]*(fs->GetCovarianceMatrix()[i][j]);
    }
  }
  assert(err2>=0);
  return sqrt(err2);
  }*/


//--------------------------------------------------------------------------------------------------
// perform fit of recoil component
void performFit(const vector<TH1D*> hv, const vector<TH1D*> hbkgv, const Double_t *ptbins, const Int_t nbins,
                const Int_t model, const Bool_t sigOnly,
                TCanvas *c, const char *plabel, const char *xlabel,
                Double_t *meanArr,   Double_t *meanErrArr,
                Double_t *sigma0Arr, Double_t *sigma0ErrArr,
                Double_t *sigma1Arr, Double_t *sigma1ErrArr,
                Double_t *sigma2Arr, Double_t *sigma2ErrArr,
                Double_t *sigma3Arr, Double_t *sigma3ErrArr,
                Double_t *frac2Arr,  Double_t *frac2ErrArr,
                Double_t *frac3Arr,  Double_t *frac3ErrArr);


//=== MAIN MACRO ================================================================================================= 

void footprint_correct(TString infilename="/data/blue/Bacon/Run2/wz_flat_final_read/Wenu/ntuples/we_select.raw.root",  // input ntuple
//void footprint_correct(TString infilename="/data/blue/jlawhorn/Wenu/ntuples/we_select.raw.root",  // input ntuple
		       std::string uparName = "puppiU1",
		       std::string uprpName = "puppiU2",
		       std::string metName = "puppi",
		       TString outputDir="test" // output directory
) {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
 
  
  CPlot::sOutDir = outputDir + TString("/plots");

  vector<TString> fnamev;
  vector<Bool_t> isBkgv;
  fnamev.push_back(infilename); isBkgv.push_back(kFALSE);
  
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.5;
     
 
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  TFile *infile = 0;
  TTree *intree = 0;  
  
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  TLorentzVector *lep=0;
  TLorentzVector *sc=0; 
  TLorentzVector *genLep=0;
  TLorentzVector *genPreLep=0;

  vector<Float_t> upper;
  upper.push_back(2.5);
  upper.push_back(2.1); 
  upper.push_back(1.5);
  upper.push_back(0);
  upper.push_back(-1.5); 
  upper.push_back(-2.1); 
  vector<Float_t> lower; 
  lower.push_back(2.1); 
  lower.push_back(1.5);
  lower.push_back(0);
  lower.push_back(-1.5); 
  lower.push_back(-2.1);
  lower.push_back(-2.5); 

  assert(upper.size()==lower.size());

  Int_t nbins=upper.size();

  char title[50];
  vector<TProfile*> pBias; 
  for (Int_t i=0; i<nbins; i++) {
    sprintf(title,"pBias_%i",i);
    pBias.push_back(new TProfile(title,"",10,-0.2,0.2,"e(J)"));
  }

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);
    intree = (TTree*)infile->Get("Events");

    intree->SetBranchAddress("runNum",   &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    intree->SetBranchAddress("npv",      &npv);       // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);       // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genVy",    &genVy);     // GEN boson rapidity (signal MC)
    intree->SetBranchAddress("genVMass", &genVMass);  // GEN boson mass (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("puppiMet",    &met);       // MET
    intree->SetBranchAddress("puppiMetPhi", &metPhi);    // phi(MET)
    intree->SetBranchAddress("sumEt",    &sumEt);     // Sum ET
    intree->SetBranchAddress("mt",       &mt);        // transverse mass
    intree->SetBranchAddress("puppiU1", &u1);  // parallel component of recoil      
    intree->SetBranchAddress("puppiU2", &u2);  // perpendicular component of recoil
    intree->SetBranchAddress("q",        &q);         // lepton charge
    intree->SetBranchAddress("lep",      &lep);       // lepton 4-vector
    intree->SetBranchAddress("genLep",   &genLep);    // lepton 4-vector
    intree->SetBranchAddress("genLep", &genPreLep);  // lepton 4-vector
    intree->SetBranchAddress("sc",       &sc);        // electron Supercluster 4-vector
  
    //
    // Loop over events
    //
    for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    //for(Int_t ientry=0; ientry<100; ientry++) {
      intree->GetEntry(ientry);
      
      if(sc->Pt()        < PT_CUT)  continue;   
      if(fabs(sc->Eta()) > ETA_CUT) continue;
      if(fabs(lep->Eta()) > ETA_CUT) continue;
      if(lep->Pt()        < PT_CUT)  continue;
      
      TVector2 vMET;    vMET.Set(met*cos(metPhi), met*sin(metPhi));
      TVector2 vBoson;  vBoson.Set(genVPt*cos(genVPhi), genVPt*sin(genVPhi));
      TVector2 vLep;    vLep.Set(lep->Pt()*cos(lep->Phi()),lep->Pt()*sin(lep->Phi()));
      TVector2 vGenLep; vGenLep.Set(genPreLep->Pt()*cos(genPreLep->Phi()),genPreLep->Pt()*sin(genPreLep->Phi()));
      TVector2 vReco = -1.0*(vMET+vLep);

      //cout << ((vBoson.Px())*(vReco.Px())+((vBoson.Py())*(vReco.Py())))/(vBoson.Mod()) << ", " << u1 << endl;
      //cout << ((vBoson.Px())*(vReco.Py())-((vBoson.Py())*(vReco.Px())))/(vBoson.Mod()) << ", " << u2 << endl;
      Double_t u1e=((vLep.Px())*(vReco.Px())+((vLep.Py())*(vReco.Py())))/(vLep.Mod());
      Double_t u2e=((vLep.Px())*(vReco.Py())-((vLep.Py())*(vReco.Px())))/(vLep.Mod());

      //cout << sqrt(u1*u1+u2*u2) << ", " << sqrt(u1e*u1e+u2e*u2e) << ", " << vReco.Mod() << endl;

      for (Int_t i=0; i<nbins; i++) {
	if( sc->Eta()>lower[i] && sc->Eta()<upper[i] ) {
	  pBias[i]->Fill(cos(vLep.DeltaPhi(vBoson)),u1e/lep->Pt());
	}
      }

    }
    
    delete infile;
    infile=0, intree=0;
  }  

  TCanvas *c = new TCanvas("c","",800,600);

  vector<Float_t> correction;
  for (Int_t i=0; i<nbins; i++) {
    TF1 *f = new TF1("f","pol1",-0.1,0.1);
    pBias[i]->Fit(f);
    //cout << lower[i] << " < |eta| < " << upper[i] << ": ";
    correction.push_back(f->Eval(0.0));
    pBias[i]->Draw();
    pBias[i]->GetYaxis()->SetTitle("Bias");
    pBias[i]->GetXaxis()->SetTitle("cos(deltaPhi)");
    f->SetLineColor(kRed);
    f->Draw("same");
    sprintf(title,"pbias_%i.png",i);
    c->SaveAs(title);
    //delete f;
  }

  TGraph *gr = new TGraph();

  for (Int_t i=0; i<nbins; i++) {
    cout << lower[i] << " < |eta| < " << upper[i] << ": " << correction[i] << endl;
    gr->SetPoint(i,0.5*(lower[i]+upper[i]),correction[i]);
  }

  gr->Draw("ap");

}
