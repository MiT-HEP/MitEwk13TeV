//================================================================================================
//
// Perform fits to recoil against Z->ee events
//
//  * Outputs a ROOT file of functions parametrizing the distribution of recoil components
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>                   // standard I/O
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TF1.h>                      // 1D function
#include <TFitResult.h>               // class to handle fit results
#include <TGraphErrors.h>             // graph class
#include "Math/LorentzVector.h"       // 4-vector class

#include "../Utils/CPlot.hh"          // helper class for plots
#include "../Utils/MitStyleRemix.hh"  // style settings for drawing

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

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== FUNCTION DECLARATIONS ======================================================================================

// generate web page
void makeHTML(const TString outDir,  const Int_t nbins, 
              const Int_t pfu1model, const Int_t pfu2model);

//--------------------------------------------------------------------------------------------------
// function to describe Gaussian widths as a function of dilepton pT
Double_t sigmaFunc(Double_t *x, Double_t *par) {
  // par[0]: quadratic coefficient
  // par[1]: linear coefficient
  // par[2]: constant term
  
  Double_t a  = par[0];
  Double_t b  = par[1];
  Double_t c  = par[2];
    
  return a*x[0]*x[0] + b*x[0] + c;
}

//--------------------------------------------------------------------------------------------------
// function to describe relative fraction in a double Gaussian based on 
// functions for sigma0, sigma1, and sigma2
Double_t frac2Func(Double_t *x, Double_t *par) {
  // par[0..3]:  sigma0
  // par[4..7]:  sigma1
  // par[8..11]: sigma2
  
  TF1 s0("_s0",sigmaFunc,0,7000,4); s0.SetParameters(par[0],par[1],par[2],par[3]);
  TF1 s1("_s1",sigmaFunc,0,7000,4); s1.SetParameters(par[4],par[5],par[6],par[7]);
  TF1 s2("_s2",sigmaFunc,0,7000,4); s2.SetParameters(par[8],par[9],par[10],par[11]);
  
  return (s0.Eval(x[0]) - s1.Eval(x[0]))/(s2.Eval(x[0]) - s1.Eval(x[0]));
}


//--------------------------------------------------------------------------------------------------
Double_t dMean(const TF1 *fcn, const Double_t x, const TFitResultPtr fs) {
  Double_t df[2];
  df[0] = 1;
  df[1] = x;
  Double_t err2 = df[0]*df[0]*(fs->GetCovarianceMatrix()[0][0]) 
                  + df[1]*df[1]*(fs->GetCovarianceMatrix()[1][1]) 
		  + 2.0*df[0]*df[1]*(fs->GetCovarianceMatrix()[0][1]);
  assert(err2>=0);
  return sqrt(err2);
}

//--------------------------------------------------------------------------------------------------
Double_t dSigma(const TF1 *fcn, const Double_t x, const TFitResultPtr fs) {
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
}


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

void fitRecoilZmm(TString infilename,  // input ntuple
                  Int_t   pfu1model,   // u1 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
                  Int_t   pfu2model,   // u2 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
	          Bool_t  sigOnly,     // signal event only?
	          TString outputDir    // output directory
) {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
 
  
  CPlot::sOutDir = outputDir + TString("/plots");

  Double_t ptbins[] = {
     0,  5, 10, 15, 20, 25, 30, 35, 40, 50,
    60, 70, 80, 90
  };
  Int_t nbins = sizeof(ptbins)/sizeof(Double_t)-1;

  Double_t corrbins[] = { 0, 10, 30, 50 };
  Int_t ncorrbins = sizeof(corrbins)/sizeof(Double_t)-1;

  TString formulaPFu1mean("pol1");
  TString formulaPFu2mean("pol1");

  vector<TString> fnamev;
  vector<Bool_t> isBkgv;
  fnamev.push_back(infilename); isBkgv.push_back(kFALSE);
  fnamev.push_back("/scratch/ksung/EWKAna/8TeV/Selection/Zmumu/ntuples/top_select.root"); isBkgv.push_back(kTRUE); 
  fnamev.push_back("/scratch/ksung/EWKAna/8TeV/Selection/Zmumu/ntuples/ewk_select.root"); isBkgv.push_back(kTRUE);
  
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.1;
     
 
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  char hname[100];
  vector<TH1D*> hPFu1v,  hPFu1Bkgv;
  vector<TH1D*> hPFu2v,  hPFu2Bkgv;
  for(Int_t ibin=0; ibin<nbins; ibin++) {
    sprintf(hname,"hPFu1_%i",ibin);    hPFu1v.push_back(new TH1D(hname,"",100,-100-ptbins[ibin],100-ptbins[ibin]));    hPFu1v[ibin]->Sumw2();
    sprintf(hname,"hPFu1Bkg_%i",ibin); hPFu1Bkgv.push_back(new TH1D(hname,"",100,-100-ptbins[ibin],100-ptbins[ibin])); hPFu1Bkgv[ibin]->Sumw2();
    
    sprintf(hname,"hPFu2_%i",ibin);    hPFu2v.push_back(new TH1D(hname,"",100,-100,100));    hPFu2v[ibin]->Sumw2();	       
    sprintf(hname,"hPFu2Bkg_%i",ibin); hPFu2Bkgv.push_back(new TH1D(hname,"",100,-100,100)); hPFu2Bkgv[ibin]->Sumw2();
  }

  vector<TH2D*> hPFu1u2v;
  for(Int_t ibin=0; ibin<ncorrbins; ibin++) {
    sprintf(hname,"hPFu1u2_%i",ibin);
    hPFu1u2v.push_back(new TH2D(hname,"",300,-10,10,300,-10,10));
  }
    
  TFitResultPtr fitresPFu1mean;   TF1 *fcnPFu1mean   = new TF1("fcnPFu1mean",formulaPFu1mean,0,7000);
  TFitResultPtr fitresPFu1sigma1; TF1 *fcnPFu1sigma1 = new TF1("fcnPFu1sigma1",sigmaFunc,0,7000,3);
  TFitResultPtr fitresPFu1sigma2; TF1 *fcnPFu1sigma2 = new TF1("fcnPFu1sigma2",sigmaFunc,0,7000,3);  
  TFitResultPtr fitresPFu1sigma0; TF1 *fcnPFu1sigma0 = new TF1("fcnPFu1sigma0",sigmaFunc,0,7000,3);    
  TFitResultPtr fitresPFu1frac2;  TF1 *fcnPFu1frac2  = new TF1("fcnPFu1frac2",frac2Func,0,7000,12);

  TFitResultPtr fitresPFu2mean;   TF1 *fcnPFu2mean   = new TF1("fcnPFu2mean",formulaPFu2mean,0,7000);
  TFitResultPtr fitresPFu2sigma1; TF1 *fcnPFu2sigma1 = new TF1("fcnPFu2sigma1",sigmaFunc,0,7000,3);  
  TFitResultPtr fitresPFu2sigma2; TF1 *fcnPFu2sigma2 = new TF1("fcnPFu2sigma2",sigmaFunc,0,7000,3);  
  TFitResultPtr fitresPFu2sigma0; TF1 *fcnPFu2sigma0 = new TF1("fcnPFu2sigma0",sigmaFunc,0,7000,3);  
  TFitResultPtr fitresPFu2frac2;  TF1 *fcnPFu2frac2  = new TF1("fcnPFu2frac2",frac2Func,0,7000,12);
      
//  fcnPFu1sigma1->SetParameter(0,-5e-4); fcnPFu1sigma1->SetParLimits(0,-1e-3,1e-3);
//  fcnPFu1sigma1->SetParameter(1,0.07);  fcnPFu1sigma1->SetParLimits(1,-1,1);
//  fcnPFu1sigma1->SetParameter(2,5.5);   fcnPFu1sigma1->SetParLimits(2,3,7);
//  fcnPFu1sigma2->SetParameter(0,1e-4);  fcnPFu1sigma2->SetParLimits(0,-2e-2,2e-2);
//  fcnPFu1sigma2->SetParameter(1,0.06);  fcnPFu1sigma2->SetParLimits(1,-1,1);
//  fcnPFu1sigma2->SetParameter(2,9.5);   fcnPFu1sigma2->SetParLimits(2,5,15);
//  fcnPFu1sigma0->SetParameter(0,-1e-4); fcnPFu1sigma0->SetParLimits(0,-1e-2,1e-2);
//  fcnPFu1sigma0->SetParameter(1,0.07);  fcnPFu1sigma0->SetParLimits(1,-1,1);
  fcnPFu1sigma0->SetParameter(2,7);     fcnPFu1sigma0->SetParLimits(2,5,9); 
  
//  fcnPFu2sigma1->SetParameter(0,-1e-4);  fcnPFu2sigma1->SetParLimits(0,-5e-3,5e-3);
//  fcnPFu2sigma1->SetParameter(1,0.05);   fcnPFu2sigma1->SetParLimits(1,-1,1);
//  fcnPFu2sigma1->SetParameter(2,5);      fcnPFu2sigma1->SetParLimits(2,2,7);
//  fcnPFu2sigma2->SetParameter(0,-1e-4);  fcnPFu2sigma2->SetParLimits(0,-5e-3,5e-3);
//  fcnPFu2sigma2->SetParameter(1,0.05);   fcnPFu2sigma2->SetParLimits(1,-1,1);
//  fcnPFu2sigma2->SetParameter(2,9);      fcnPFu2sigma2->SetParLimits(2,5,15);
//  fcnPFu2sigma0->SetParameter(0,-1e-4);  fcnPFu2sigma0->SetParLimits(0,-1e-2,1e-2);
//  fcnPFu2sigma0->SetParameter(1,0.03);   fcnPFu2sigma0->SetParLimits(1,-1,1);
//  fcnPFu2sigma0->SetParameter(2,7);      fcnPFu2sigma0->SetParLimits(2,5,9);
    
  TFile *infile = 0;
  TTree *intree = 0;  
  
  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, u1, u2;
  Int_t   q1, q2;
  LorentzVector *dilep=0, *lep1=0, *lep2=0;
  

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);
    intree = (TTree*)infile->Get("Events");
  
    intree->SetBranchAddress("runNum",   &runNum);     // event run number
    intree->SetBranchAddress("lumiSec",	 &lumiSec);    // event lumi section
    intree->SetBranchAddress("evtNum",	 &evtNum);     // event number
    intree->SetBranchAddress("matchGen", &matchGen);   // event has both leptons matched to MC Z->ll
    intree->SetBranchAddress("category", &category);   // dilepton category
    intree->SetBranchAddress("npv",	 &npv);        // number of primary vertices
    intree->SetBranchAddress("npu",	 &npu);        // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);     // GEN boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);    // GEN boson phi (signal MC)
    intree->SetBranchAddress("genVy",    &genVy);      // GEN boson rapidity (signal MC)
    intree->SetBranchAddress("genVMass", &genVMass);   // GEN boson mass (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);   // event weight per 1/fb (MC)
    intree->SetBranchAddress("met",	 &met);        // MET
    intree->SetBranchAddress("metPhi",	 &metPhi);     // phi(MET)
    intree->SetBranchAddress("sumEt",	 &sumEt);      // Sum ET
    intree->SetBranchAddress("u1",	 &u1);         // parallel component of recoil
    intree->SetBranchAddress("u2",	 &u2);         // perpendicular component of recoil
    intree->SetBranchAddress("q1",	 &q1);         // charge of tag lepton
    intree->SetBranchAddress("q2",	 &q2);         // charge of probe lepton
    intree->SetBranchAddress("dilep",	 &dilep);      // dilepton 4-vector
    intree->SetBranchAddress("lep1",	 &lep1);       // tag lepton 4-vector
    intree->SetBranchAddress("lep2",	 &lep2);       // probe lepton 4-vector 
  
    //
    // Loop over events
    //
    for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
    
      if(category!=1 && category!=2)                                 continue;
      if(dilep->M() < MASS_LOW || dilep->M() > MASS_HIGH)            continue;
      if(lep1->Pt()        < PT_CUT  || lep2->Pt()        < PT_CUT)  continue;
      if(fabs(lep1->Eta()) > ETA_CUT || fabs(lep2->Eta()) > ETA_CUT) continue;
    
      Int_t ipt=-1;
      for(Int_t ibin=0; ibin<nbins; ibin++) {
        if(dilep->Pt() > ptbins[ibin] && dilep->Pt() <= ptbins[ibin+1])
          ipt = ibin;
      }
      if(ipt<0) continue;
    
      if(isBkgv[ifile]) {
        hPFu1Bkgv[ipt]->Fill(u1);
        hPFu2Bkgv[ipt]->Fill(u2);
      
      } else {
        hPFu1v[ipt]->Fill(u1);
        hPFu2v[ipt]->Fill(u2);
      }
    }
    
    delete infile;
    infile=0, intree=0;   
  }  
  
  Double_t xval[nbins], xerr[nbins];
  for(Int_t ibin=0; ibin<nbins; ibin++) {
    xval[ibin] = 0.5*(ptbins[ibin+1]+ptbins[ibin]);
    xerr[ibin] = 0.5*(ptbins[ibin+1]-ptbins[ibin]);
  }
  
  //
  // Arrays and graphs to store fit results
  //
  TGraphErrors *grPFu1mean=0;   Double_t pfu1Mean[nbins],   pfu1MeanErr[nbins];
  TGraphErrors *grPFu1sigma0=0; Double_t pfu1Sigma0[nbins], pfu1Sigma0Err[nbins];
  TGraphErrors *grPFu1sigma1=0; Double_t pfu1Sigma1[nbins], pfu1Sigma1Err[nbins];
  TGraphErrors *grPFu1sigma2=0; Double_t pfu1Sigma2[nbins], pfu1Sigma2Err[nbins];
  TGraphErrors *grPFu1sigma3=0; Double_t pfu1Sigma3[nbins], pfu1Sigma3Err[nbins];
  TGraphErrors *grPFu1frac2=0;  Double_t pfu1Frac2[nbins],  pfu1Frac2Err[nbins];  
  TGraphErrors *grPFu1frac3=0;  Double_t pfu1Frac3[nbins],  pfu1Frac3Err[nbins];

  TGraphErrors *grPFu2mean=0;   Double_t pfu2Mean[nbins],   pfu2MeanErr[nbins];
  TGraphErrors *grPFu2sigma0=0; Double_t pfu2Sigma0[nbins], pfu2Sigma0Err[nbins];
  TGraphErrors *grPFu2sigma1=0; Double_t pfu2Sigma1[nbins], pfu2Sigma1Err[nbins];
  TGraphErrors *grPFu2sigma2=0; Double_t pfu2Sigma2[nbins], pfu2Sigma2Err[nbins];
  TGraphErrors *grPFu2sigma3=0; Double_t pfu2Sigma3[nbins], pfu2Sigma3Err[nbins];
  TGraphErrors *grPFu2frac2=0;  Double_t pfu2Frac2[nbins],  pfu2Frac2Err[nbins];  
  TGraphErrors *grPFu2frac3=0;  Double_t pfu2Frac3[nbins],  pfu2Frac3Err[nbins];
        
  
  TCanvas *c = MakeCanvas("c","c",800,600);

  // Fitting PF-MET u1
  performFit(hPFu1v, hPFu1Bkgv, ptbins, nbins, pfu1model, sigOnly,
             c, "pfu1", "PF u_{1} [GeV]",
	     pfu1Mean,   pfu1MeanErr,
	     pfu1Sigma0, pfu1Sigma0Err,
	     pfu1Sigma1, pfu1Sigma1Err,
	     pfu1Sigma2, pfu1Sigma2Err,
	     pfu1Sigma3, pfu1Sigma3Err,
	     pfu1Frac2,  pfu1Frac2Err,
	     pfu1Frac3,  pfu1Frac3Err);

  // Fitting PF-MET u2
  performFit(hPFu2v, hPFu2Bkgv, ptbins, nbins, pfu2model, sigOnly,
             c, "pfu2", "PF u_{2} [GeV]",
	     pfu2Mean,   pfu2MeanErr,
	     pfu2Sigma0, pfu2Sigma0Err,
	     pfu2Sigma1, pfu2Sigma1Err,
	     pfu2Sigma2, pfu2Sigma2Err,
	     pfu2Sigma3, pfu2Sigma3Err,
	     pfu2Frac2,  pfu2Frac2Err,
	     pfu2Frac3,  pfu2Frac3Err);  

 
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
   
  char fitparam[100];
  char chi2ndf[50]; 
    
//  TGraphErrors *errBand = new TGraphErrors(nbins+2);
  
  //
  // Plotting PF-MET u1 vs. dilepton pT
  //
  grPFu1mean = new TGraphErrors(nbins,xval,pfu1Mean,xerr,pfu1MeanErr);
  grPFu1mean->SetName("grPFu1mean");
  fitresPFu1mean = grPFu1mean->Fit("fcnPFu1mean","QMRN0FBSE");
  sprintf(chi2ndf,"#chi^{2}/ndf = %.2f",(fcnPFu1mean->GetChisquare())/(fcnPFu1mean->GetNDF()));
//  errBand->SetPoint(0,0.002*(xval[nbins-1]),fcnPFu1mean->Eval(0.002*(xval[nbins-1])));
//  errBand->SetPointError(0,0,dMean(fcnPFu1mean,0.002*(xval[nbins-1]),fitresPFu1mean));
//  for(Int_t i=1; i<=nbins; i++) {
//    errBand->SetPoint(i,xval[i-1],fcnPFu1mean->Eval(xval[i-1]));
//    errBand->SetPointError(i,0,dMean(fcnPFu1mean,xval[i-1],fitresPFu1mean));
//  }
//  errBand->SetPoint(nbins+1,1.2*(xval[nbins-1]),fcnPFu1mean->Eval(1.2*(xval[nbins-1])));
//  errBand->SetPointError(nbins+1,0,dMean(fcnPFu1mean,1.2*(xval[nbins-1]),fitresPFu1mean));  
  CPlot plotPFu1mean("pfu1mean","","p_{T}(ll) [GeV/c]","#mu(u_{1}) [GeV]");
  plotPFu1mean.AddGraph(grPFu1mean,"");
//  plotPFu1mean.AddGraph(errBand,"3",kAzure-9,kFullDotSmall);
  plotPFu1mean.AddGraph(grPFu1mean,"",kBlack,kOpenCircle);
  plotPFu1mean.AddFcn(fcnPFu1mean,kRed);
  plotPFu1mean.AddTextBox(chi2ndf,0.65,0.87,0.95,0.82,0,kBlack,-1);
  sprintf(fitparam,"p_{0} = %.3f #pm %.3f",fcnPFu1mean->GetParameter(0),fcnPFu1mean->GetParError(0)); plotPFu1mean.AddTextBox(fitparam,0.65,0.80,0.95,0.75,0,kBlack,-1);
  sprintf(fitparam,"p_{1} = %.3f #pm %.3f",fcnPFu1mean->GetParameter(1),fcnPFu1mean->GetParError(1)); plotPFu1mean.AddTextBox(fitparam,0.65,0.75,0.95,0.70,0,kBlack,-1);
  plotPFu1mean.Draw(c,kTRUE,"png");
  
  grPFu1sigma1 = new TGraphErrors(nbins,xval,pfu1Sigma1,xerr,pfu1Sigma1Err);  
  grPFu1sigma1->SetName("grPFu1sigma1");
  fitresPFu1sigma1 = grPFu1sigma1->Fit("fcnPFu1sigma1","QMRN0SE");
  sprintf(chi2ndf,"#chi^{2}/ndf = %.2f",(fcnPFu1sigma1->GetChisquare())/(fcnPFu1sigma1->GetNDF()));
//  errBand->SetPoint(0,0.002*(xval[nbins-1]),fcnPFu1sigma1->Eval(0.002*(xval[nbins-1])));
//  errBand->SetPointError(0,0,dSigma(fcnPFu1sigma1,0.002*(xval[nbins-1]),fitresPFu1sigma1));
//  for(Int_t i=1; i<=nbins; i++) {
//    errBand->SetPoint(i,xval[i-1],fcnPFu1sigma1->Eval(xval[i-1]));
//    errBand->SetPointError(i,0,dSigma(fcnPFu1sigma1,xval[i-1],fitresPFu1sigma1));
//  }
//  errBand->SetPoint(nbins+1,1.2*(xval[nbins-1]),fcnPFu1sigma1->Eval(1.2*(xval[nbins-1])));
//  errBand->SetPointError(nbins+1,0,dSigma(fcnPFu1sigma1,1.2*(xval[nbins-1]),fitresPFu1sigma1));  
  CPlot plotPFu1sigma1("pfu1sigma1","","p_{T}(ll) [GeV/c]","#sigma_{1}(u_{1}) [GeV]");
  plotPFu1sigma1.AddGraph(grPFu1sigma1,"");
//  plotPFu1sigma1.AddGraph(errBand,"3",kAzure-9,kFullDotSmall);
  plotPFu1sigma1.AddGraph(grPFu1sigma1,"",kBlack,kOpenCircle);
  plotPFu1sigma1.AddFcn(fcnPFu1sigma1,kRed);
  plotPFu1sigma1.AddTextBox(chi2ndf,0.21,0.87,0.41,0.82,0,kBlack,-1);
  sprintf(fitparam,"p_{0} = (%.1f #pm %.1f) #times 10^{-5}",1e5*(fcnPFu1sigma1->GetParameter(0)),1e5*(fcnPFu1sigma1->GetParError(0))); 
  plotPFu1sigma1.AddTextBox(fitparam,0.21,0.80,0.51,0.75,0,kBlack,-1);
  sprintf(fitparam,"p_{1} = %.3f #pm %.3f",fcnPFu1sigma1->GetParameter(1),fcnPFu1sigma1->GetParError(1)); plotPFu1sigma1.AddTextBox(fitparam,0.21,0.75,0.51,0.70,0,kBlack,-1);
  sprintf(fitparam,"p_{2} = %.3f #pm %.3f",fcnPFu1sigma1->GetParameter(2),fcnPFu1sigma1->GetParError(2)); plotPFu1sigma1.AddTextBox(fitparam,0.21,0.70,0.51,0.65,0,kBlack,-1);
//  sprintf(fitparam,"p_{3} = %.3f #pm %.3f",fcnPFu1sigma1->GetParameter(3),fcnPFu1sigma1->GetParError(3)); plotPFu1sigma1.AddTextBox(fitparam,0.21,0.65,0.51,0.60,0,kBlack,-1);
//  if(njetcut==0) plotPFu1sigma1.SetYRange(0,15);
//  else           plotPFu1sigma1.SetYRange(0,30);
  plotPFu1sigma1.Draw(c,kTRUE,"png");
  
  if(pfu1model>=2) {
    grPFu1sigma2 = new TGraphErrors(nbins,xval,pfu1Sigma2,xerr,pfu1Sigma2Err);    
    grPFu1sigma2->SetName("grPFu1sigma2");
    fitresPFu1sigma2 = grPFu1sigma2->Fit("fcnPFu1sigma2","QMRN0SE");
    sprintf(chi2ndf,"#chi^{2}/ndf = %.2f",(fcnPFu1sigma2->GetChisquare())/(fcnPFu1sigma2->GetNDF()));    
//    errBand->SetPoint(0,0.002*(xval[nbins-1]),fcnPFu1sigma2->Eval(0.002*(xval[nbins-1])));
//    errBand->SetPointError(0,0,dSigma(fcnPFu1sigma2,0.002*(xval[nbins-1]),fitresPFu1sigma2));
//    for(Int_t i=1; i<=nbins; i++) {
//      errBand->SetPoint(i,xval[i-1],fcnPFu1sigma2->Eval(xval[i-1]));
//      errBand->SetPointError(i,0,dSigma(fcnPFu1sigma2,xval[i-1],fitresPFu1sigma2));
//    }
//    errBand->SetPoint(nbins+1,1.2*(xval[nbins-1]),fcnPFu1sigma2->Eval(1.2*(xval[nbins-1])));
//    errBand->SetPointError(nbins+1,0,dSigma(fcnPFu1sigma2,1.2*(xval[nbins-1]),fitresPFu1sigma2));  
    CPlot plotPFu1sigma2("pfu1sigma2","","p_{T}(ll) [GeV/c]","#sigma_{2}(u_{1}) [GeV]");
    plotPFu1sigma2.AddGraph(grPFu1sigma2,"");
//    plotPFu1sigma2.AddGraph(errBand,"3",kAzure-9,kFullDotSmall);
    plotPFu1sigma2.AddGraph(grPFu1sigma2,"",kBlack,kOpenCircle);
    plotPFu1sigma2.AddFcn(fcnPFu1sigma2,kRed);
    plotPFu1sigma2.AddTextBox(chi2ndf,0.21,0.87,0.41,0.82,0,kBlack,-1);
    sprintf(fitparam,"p_{0} = (%.1f #pm %.1f) #times 10^{-5}",1e5*(fcnPFu1sigma2->GetParameter(0)),1e5*(fcnPFu1sigma2->GetParError(0))); 
    plotPFu1sigma2.AddTextBox(fitparam,0.21,0.80,0.51,0.75,0,kBlack,-1);
    sprintf(fitparam,"p_{1} = %.3f #pm %.3f",fcnPFu1sigma2->GetParameter(1),fcnPFu1sigma2->GetParError(1)); plotPFu1sigma2.AddTextBox(fitparam,0.21,0.75,0.51,0.70,0,kBlack,-1);
    sprintf(fitparam,"p_{2} = %.3f #pm %.3f",fcnPFu1sigma2->GetParameter(2),fcnPFu1sigma2->GetParError(2)); plotPFu1sigma2.AddTextBox(fitparam,0.21,0.70,0.51,0.65,0,kBlack,-1);
//    sprintf(fitparam,"p_{3} = %.3f #pm %.3f",fcnPFu1sigma2->GetParameter(3),fcnPFu1sigma2->GetParError(3)); plotPFu1sigma2.AddTextBox(fitparam,0.21,0.65,0.51,0.60,0,kBlack,-1);
//    if(njetcut==0) plotPFu1sigma2.SetYRange(0,30);
//    else           plotPFu1sigma2.SetYRange(0,60);
    plotPFu1sigma2.Draw(c,kTRUE,"png");


    grPFu1sigma0 = new TGraphErrors(nbins,xval,pfu1Sigma0,xerr,pfu1Sigma0Err);    
    grPFu1sigma0->SetName("grPFu1sigma0");
    fitresPFu1sigma0 = grPFu1sigma0->Fit("fcnPFu1sigma0","QMRN0SE");
    sprintf(chi2ndf,"#chi^{2}/ndf = %.2f",(fcnPFu1sigma0->GetChisquare())/(fcnPFu1sigma0->GetNDF()));    
//    errBand->SetPoint(0,0.002*(xval[nbins-1]),fcnPFu1sigma0->Eval(0.002*(xval[nbins-1])));
//    errBand->SetPointError(0,0,dSigma(fcnPFu1sigma0,0.002*(xval[nbins-1]),fitresPFu1sigma0));
//    for(Int_t i=1; i<=nbins; i++) {
//      errBand->SetPoint(i,xval[i-1],fcnPFu1sigma0->Eval(xval[i-1]));
//      errBand->SetPointError(i,0,dSigma(fcnPFu1sigma0,xval[i-1],fitresPFu1sigma0));
//    }
//    errBand->SetPoint(nbins+1,1.2*(xval[nbins-1]),fcnPFu1sigma0->Eval(1.2*(xval[nbins-1])));
//    errBand->SetPointError(nbins+1,0,dSigma(fcnPFu1sigma0,1.2*(xval[nbins-1]),fitresPFu1sigma0));    
    CPlot plotPFu1sigma0("pfu1sigma0","","p_{T}(ll) [GeV/c]","#sigma(u_{1}) [GeV]");
    plotPFu1sigma0.AddGraph(grPFu1sigma0,"");
//    plotPFu1sigma0.AddGraph(errBand,"3",kAzure-9,kFullDotSmall);
    plotPFu1sigma0.AddGraph(grPFu1sigma0,"",kBlack,kOpenCircle);
    plotPFu1sigma0.AddFcn(fcnPFu1sigma0,kRed);
    plotPFu1sigma0.AddTextBox(chi2ndf,0.21,0.87,0.41,0.82,0,kBlack,-1);
    sprintf(fitparam,"p_{0} = (%.1f #pm %.1f) #times 10^{-5}",1e5*(fcnPFu1sigma0->GetParameter(0)),1e5*(fcnPFu1sigma0->GetParError(0)));
    plotPFu1sigma0.AddTextBox(fitparam,0.21,0.80,0.51,0.75,0,kBlack,-1);
    sprintf(fitparam,"p_{1} = %.3f #pm %.3f",fcnPFu1sigma0->GetParameter(1),fcnPFu1sigma0->GetParError(1)); plotPFu1sigma0.AddTextBox(fitparam,0.21,0.75,0.51,0.70,0,kBlack,-1);
    sprintf(fitparam,"p_{2} = %.3f #pm %.3f",fcnPFu1sigma0->GetParameter(2),fcnPFu1sigma0->GetParError(2)); plotPFu1sigma0.AddTextBox(fitparam,0.21,0.70,0.51,0.65,0,kBlack,-1);
//    sprintf(fitparam,"p_{3} = %.3f #pm %.3f",fcnPFu1sigma0->GetParameter(3),fcnPFu1sigma0->GetParError(3)); plotPFu1sigma0.AddTextBox(fitparam,0.21,0.65,0.51,0.60,0,kBlack,-1);
//    if(njetcut==0) plotPFu1sigma0.SetYRange(6,15);
//    else           plotPFu1sigma0.SetYRange(6,30);
    plotPFu1sigma0.Draw(c,kTRUE,"png");
    
    grPFu1frac2 = new TGraphErrors(nbins,xval,pfu1Frac2, xerr,pfu1Frac2Err);
    grPFu1frac2->SetName("grPFu1frac2");
    fcnPFu1frac2->SetParameter(0,fcnPFu1sigma0->GetParameter(0));
    fcnPFu1frac2->SetParameter(1,fcnPFu1sigma0->GetParameter(1));
    fcnPFu1frac2->SetParameter(2,fcnPFu1sigma0->GetParameter(2));
    fcnPFu1frac2->SetParameter(3,fcnPFu1sigma0->GetParameter(3));
    fcnPFu1frac2->SetParameter(4,fcnPFu1sigma1->GetParameter(0));
    fcnPFu1frac2->SetParameter(5,fcnPFu1sigma1->GetParameter(1));
    fcnPFu1frac2->SetParameter(6,fcnPFu1sigma1->GetParameter(2));
    fcnPFu1frac2->SetParameter(7,fcnPFu1sigma1->GetParameter(3));
    fcnPFu1frac2->SetParameter(8,fcnPFu1sigma2->GetParameter(0));
    fcnPFu1frac2->SetParameter(9,fcnPFu1sigma2->GetParameter(1));
    fcnPFu1frac2->SetParameter(10,fcnPFu1sigma2->GetParameter(2));
    fcnPFu1frac2->SetParameter(11,fcnPFu1sigma2->GetParameter(3));
    CPlot plotPFu1frac2("pfu1frac2","","p_{T}(ll) [GeV/c]","f_{2}");
    plotPFu1frac2.AddGraph(grPFu1frac2,"",kBlack,kOpenCircle);
    plotPFu1frac2.AddFcn(fcnPFu1frac2,kRed);
    plotPFu1frac2.Draw(c,kTRUE,"png");
  }
  
  if(pfu1model>=3) {  
    grPFu1sigma3 = new TGraphErrors(nbins,xval,pfu1Sigma3,xerr,pfu1Sigma3Err);
    grPFu1sigma3->SetName("grPFu1sigma3");
    CPlot plotPFu1sigma3("pfu1sigma3","","p_{T}(ll) [GeV/c]","#sigma_{3} [GeV]");
    plotPFu1sigma3.AddGraph(grPFu1sigma3,"",kBlack,kOpenCircle);
    plotPFu1sigma3.Draw(c,kTRUE,"png");
  
    grPFu1frac3 = new TGraphErrors(nbins,xval,pfu1Frac3, xerr,pfu1Frac3Err);
    grPFu1frac3->SetName("grPFu1frac3");
    CPlot plotPFu1frac3("pfu1frac3","","p_{T}(ll) [GeV/c]","f_{3}");
    plotPFu1frac3.AddGraph(grPFu1frac3,"",kBlack,kOpenCircle);
    plotPFu1frac3.Draw(c,kTRUE,"png");
  }

  //
  // Plotting PF-MET u2 vs. dilepton pT
  //
  grPFu2mean = new TGraphErrors(nbins,xval,pfu2Mean,xerr,pfu2MeanErr);
  grPFu2mean->SetName("grPFu2mean");
  fitresPFu2mean = grPFu2mean->Fit("fcnPFu2mean","QMRN0FBSE");
  sprintf(chi2ndf,"#chi^{2}/ndf = %.2f",(fcnPFu2mean->GetChisquare())/(fcnPFu2mean->GetNDF()));  
//  errBand->SetPoint(0,0.002*(xval[nbins-1]),fcnPFu2mean->Eval(0.002*(xval[nbins-1])));
//  errBand->SetPointError(0,0,dMean(fcnPFu2mean,0.002*(xval[nbins-1]),fitresPFu2mean));
//  for(Int_t i=1; i<=nbins; i++) {
//    errBand->SetPoint(i,xval[i-1],fcnPFu2mean->Eval(xval[i-1]));
//    errBand->SetPointError(i,0,dMean(fcnPFu2mean,xval[i-1],fitresPFu2mean));
//  }
//  errBand->SetPoint(nbins+1,1.2*(xval[nbins-1]),fcnPFu2mean->Eval(1.2*(xval[nbins-1])));
//  errBand->SetPointError(nbins+1,0,dMean(fcnPFu2mean,1.2*(xval[nbins-1]),fitresPFu2mean));    
  CPlot plotPFu2mean("pfu2mean","","p_{T}(ll) [GeV/c]","#mu(u_{2}) [GeV]");
  plotPFu2mean.AddGraph(grPFu2mean,"");
//  plotPFu2mean.AddGraph(errBand,"3",kAzure-9,kFullDotSmall);
  plotPFu2mean.AddGraph(grPFu2mean,"",kBlack,kOpenCircle);
  plotPFu2mean.AddFcn(fcnPFu2mean,kRed);
  plotPFu2mean.AddTextBox(chi2ndf,0.21,0.87,0.41,0.82,0,kBlack,-1);
  sprintf(fitparam,"p_{0} = %.3f #pm %.3f",fcnPFu2mean->GetParameter(0),fcnPFu2mean->GetParError(0)); plotPFu2mean.AddTextBox(fitparam,0.21,0.80,0.51,0.75,0,kBlack,-1);
  sprintf(fitparam,"p_{1} = %.3f #pm %.3f",fcnPFu2mean->GetParameter(1),fcnPFu2mean->GetParError(1)); plotPFu2mean.AddTextBox(fitparam,0.21,0.75,0.51,0.70,0,kBlack,-1);
  plotPFu2mean.Draw(c,kTRUE,"png");

  
  grPFu2sigma1 = new TGraphErrors(nbins,xval,pfu2Sigma1,xerr,pfu2Sigma1Err);
  grPFu2sigma1->SetName("grPFu2sigma1");
  fitresPFu2sigma1 = grPFu2sigma1->Fit("fcnPFu2sigma1","QMRN0SE");
  sprintf(chi2ndf,"#chi^{2}/ndf = %.2f",(fcnPFu2sigma1->GetChisquare())/(fcnPFu2sigma1->GetNDF()));  
//  errBand->SetPoint(0,0.002*(xval[nbins-1]),fcnPFu2sigma1->Eval(0.002*(xval[nbins-1])));
//  errBand->SetPointError(0,0,dSigma(fcnPFu2sigma1,0.002*(xval[nbins-1]),fitresPFu2sigma1));
//  for(Int_t i=1; i<=nbins; i++) {
//    errBand->SetPoint(i,xval[i-1],fcnPFu2sigma1->Eval(xval[i-1]));
//    errBand->SetPointError(i,0,dSigma(fcnPFu2sigma1,xval[i-1],fitresPFu2sigma1));
//  }
//  errBand->SetPoint(nbins+1,1.2*(xval[nbins-1]),fcnPFu2sigma1->Eval(1.2*(xval[nbins-1])));
//  errBand->SetPointError(nbins+1,0,dSigma(fcnPFu2sigma1,1.2*(xval[nbins-1]),fitresPFu2sigma1));    
  CPlot plotPFu2sigma1("pfu2sigma1","","p_{T}(ll) [GeV/c]","#sigma_{1}(u_{2}) [GeV]");
  plotPFu2sigma1.AddGraph(grPFu2sigma1,"");
//  plotPFu2sigma1.AddGraph(errBand,"3",kAzure-9,kFullDotSmall);
  plotPFu2sigma1.AddGraph(grPFu2sigma1,"",kBlack,kOpenCircle);
  plotPFu2sigma1.AddFcn(fcnPFu2sigma1,kRed);
  plotPFu2sigma1.AddTextBox(chi2ndf,0.21,0.87,0.41,0.82,0,kBlack,-1);
  sprintf(fitparam,"p_{0} = (%.1f #pm %.1f) #times 10^{-5}",1e5*(fcnPFu2sigma1->GetParameter(0)),1e5*(fcnPFu2sigma1->GetParError(0)));
  plotPFu2sigma1.AddTextBox(fitparam,0.21,0.80,0.51,0.75,0,kBlack,-1);
  sprintf(fitparam,"p_{1} = %.3f #pm %.3f",fcnPFu2sigma1->GetParameter(1),fcnPFu2sigma1->GetParError(1)); plotPFu2sigma1.AddTextBox(fitparam,0.21,0.75,0.51,0.70,0,kBlack,-1);
  sprintf(fitparam,"p_{2} = %.3f #pm %.3f",fcnPFu2sigma1->GetParameter(2),fcnPFu2sigma1->GetParError(2)); plotPFu2sigma1.AddTextBox(fitparam,0.21,0.70,0.51,0.65,0,kBlack,-1);
//  sprintf(fitparam,"p_{3} = %.3f #pm %.3f",fcnPFu2sigma1->GetParameter(3),fcnPFu2sigma1->GetParError(3)); plotPFu2sigma1.AddTextBox(fitparam,0.21,0.65,0.51,0.60,0,kBlack,-1);  
//  if(njetcut==0) plotPFu2sigma1.SetYRange(0,15);
//  else           plotPFu2sigma1.SetYRange(0,20);
  plotPFu2sigma1.Draw(c,kTRUE,"png");
  
  if(pfu2model>=2) {
    grPFu2sigma2 = new TGraphErrors(nbins,xval,pfu2Sigma2,xerr,pfu2Sigma2Err);
    grPFu2sigma2->SetName("grPFu2sigma2");
    fitresPFu2sigma2 = grPFu2sigma2->Fit("fcnPFu2sigma2","QMRN0SE");
    sprintf(chi2ndf,"#chi^{2}/ndf = %.2f",(fcnPFu2sigma2->GetChisquare())/(fcnPFu2sigma2->GetNDF()));        
//    errBand->SetPoint(0,0.002*(xval[nbins-1]),fcnPFu2sigma2->Eval(0.002*(xval[nbins-1])));
//    errBand->SetPointError(0,0,dSigma(fcnPFu2sigma2,0.002*(xval[nbins-1]),fitresPFu2sigma2));
//    for(Int_t i=1; i<=nbins; i++) {
//      errBand->SetPoint(i,xval[i-1],fcnPFu2sigma2->Eval(xval[i-1]));
//      errBand->SetPointError(i,0,dSigma(fcnPFu2sigma2,xval[i-1],fitresPFu2sigma2));
//    }
//    errBand->SetPoint(nbins+1,1.2*(xval[nbins-1]),fcnPFu2sigma2->Eval(1.2*(xval[nbins-1])));
//    errBand->SetPointError(nbins+1,0,dSigma(fcnPFu2sigma2,1.2*(xval[nbins-1]),fitresPFu2sigma2));    
    CPlot plotPFu2sigma2("pfu2sigma2","","p_{T}(ll) [GeV/c]","#sigma_{2}(u_{2}) [GeV]");    
    plotPFu2sigma2.AddGraph(grPFu2sigma2,"");
//    plotPFu2sigma2.AddGraph(errBand,"3",kAzure-9,kFullDotSmall);
    plotPFu2sigma2.AddGraph(grPFu2sigma2,"",kBlack,kOpenCircle);
    plotPFu2sigma2.AddFcn(fcnPFu2sigma2,kRed);
    plotPFu2sigma2.AddTextBox(chi2ndf,0.21,0.87,0.41,0.82,0,kBlack,-1);
    sprintf(fitparam,"p_{0} = (%.1f #pm %.1f) #times 10^{-5}",1e5*(fcnPFu2sigma2->GetParameter(0)),1e5*(fcnPFu2sigma2->GetParError(0)));
    plotPFu2sigma2.AddTextBox(fitparam,0.21,0.80,0.51,0.75,0,kBlack,-1);
    sprintf(fitparam,"p_{1} = %.3f #pm %.3f",fcnPFu2sigma2->GetParameter(1),fcnPFu2sigma2->GetParError(1)); plotPFu2sigma2.AddTextBox(fitparam,0.21,0.75,0.51,0.70,0,kBlack,-1);
    sprintf(fitparam,"p_{2} = %.3f #pm %.3f",fcnPFu2sigma2->GetParameter(2),fcnPFu2sigma2->GetParError(2)); plotPFu2sigma2.AddTextBox(fitparam,0.21,0.70,0.51,0.65,0,kBlack,-1);
//    sprintf(fitparam,"p_{3} = %.3f #pm %.3f",fcnPFu2sigma2->GetParameter(3),fcnPFu2sigma2->GetParError(3)); plotPFu2sigma2.AddTextBox(fitparam,0.21,0.65,0.51,0.60,0,kBlack,-1);    
//    if(njetcut==0) plotPFu2sigma2.SetYRange(0,30);
//    else           plotPFu2sigma2.SetYRange(0,40);
    plotPFu2sigma2.Draw(c,kTRUE,"png");


    grPFu2sigma0 = new TGraphErrors(nbins,xval,pfu2Sigma0,xerr,pfu2Sigma0Err);
    grPFu2sigma0->SetName("grPFu2sigma0");
    fitresPFu2sigma0 = grPFu2sigma0->Fit("fcnPFu2sigma0","QMRN0SE");
    sprintf(chi2ndf,"#chi^{2}/ndf = %.2f",(fcnPFu2sigma0->GetChisquare())/(fcnPFu2sigma0->GetNDF()));    
//    errBand->SetPoint(0,0.002*(xval[nbins-1]),fcnPFu2sigma0->Eval(0.002*(xval[nbins-1])));
//    errBand->SetPointError(0,0,dSigma(fcnPFu2sigma0,0.002*(xval[nbins-1]),fitresPFu2sigma0));
//    for(Int_t i=1; i<=nbins; i++) {
//      errBand->SetPoint(i,xval[i-1],fcnPFu2sigma0->Eval(xval[i-1]));
//      errBand->SetPointError(i,0,dSigma(fcnPFu2sigma0,xval[i-1],fitresPFu2sigma0));
//    }
//    errBand->SetPoint(nbins+1,1.2*(xval[nbins-1]),fcnPFu2sigma0->Eval(1.2*(xval[nbins-1])));
//    errBand->SetPointError(nbins+1,0,dSigma(fcnPFu2sigma0,1.2*(xval[nbins-1]),fitresPFu2sigma0));    
    CPlot plotPFu2sigma0("pfu2sigma0","","p_{T}(ll) [GeV/c]","#sigma(u_{2}) [GeV]");
    plotPFu2sigma0.AddGraph(grPFu2sigma0,"");
//    plotPFu2sigma0.AddGraph(errBand,"3",kAzure-9,kFullDotSmall);
    plotPFu2sigma0.AddGraph(grPFu2sigma0,"",kBlack,kOpenCircle);
    plotPFu2sigma0.AddFcn(fcnPFu2sigma0,kRed);
    plotPFu2sigma0.AddTextBox(chi2ndf,0.21,0.87,0.41,0.82,0,kBlack,-1);
    sprintf(fitparam,"p_{0} = (%.1f #pm %.1f) #times 10^{-5}",1e5*(fcnPFu2sigma0->GetParameter(0)),1e5*(fcnPFu2sigma0->GetParError(0)));
    plotPFu2sigma0.AddTextBox(fitparam,0.21,0.80,0.51,0.75,0,kBlack,-1);
    sprintf(fitparam,"p_{1} = %.3f #pm %.3f",fcnPFu2sigma0->GetParameter(1),fcnPFu2sigma0->GetParError(1)); plotPFu2sigma0.AddTextBox(fitparam,0.21,0.75,0.51,0.70,0,kBlack,-1);
    sprintf(fitparam,"p_{2} = %.3f #pm %.3f",fcnPFu2sigma0->GetParameter(2),fcnPFu2sigma0->GetParError(2)); plotPFu2sigma0.AddTextBox(fitparam,0.21,0.70,0.51,0.65,0,kBlack,-1);
//    sprintf(fitparam,"p_{3} = %.3f #pm %.3f",fcnPFu2sigma0->GetParameter(3),fcnPFu2sigma0->GetParError(3)); plotPFu2sigma0.AddTextBox(fitparam,0.21,0.65,0.51,0.60,0,kBlack,-1);    
//    if(njetcut==0) plotPFu2sigma0.SetYRange(6,15);
//    else           plotPFu2sigma0.SetYRange(6,20);
    plotPFu2sigma0.Draw(c,kTRUE,"png");

    
    grPFu2frac2 = new TGraphErrors(nbins,xval,pfu2Frac2, xerr,pfu2Frac2Err);
    grPFu2frac2->SetName("grPFu2frac2");
    fcnPFu2frac2->SetParameter(0,fcnPFu2sigma0->GetParameter(0));
    fcnPFu2frac2->SetParameter(1,fcnPFu2sigma0->GetParameter(1));
    fcnPFu2frac2->SetParameter(2,fcnPFu2sigma0->GetParameter(2));
    fcnPFu2frac2->SetParameter(3,fcnPFu2sigma0->GetParameter(3));
    fcnPFu2frac2->SetParameter(4,fcnPFu2sigma1->GetParameter(0));
    fcnPFu2frac2->SetParameter(5,fcnPFu2sigma1->GetParameter(1));
    fcnPFu2frac2->SetParameter(6,fcnPFu2sigma1->GetParameter(2));
    fcnPFu2frac2->SetParameter(7,fcnPFu2sigma1->GetParameter(3));
    fcnPFu2frac2->SetParameter(8,fcnPFu2sigma2->GetParameter(0));
    fcnPFu2frac2->SetParameter(9,fcnPFu2sigma2->GetParameter(1));
    fcnPFu2frac2->SetParameter(10,fcnPFu2sigma2->GetParameter(2));
    fcnPFu2frac2->SetParameter(11,fcnPFu2sigma2->GetParameter(3));
    CPlot plotPFu2frac2("pfu2frac2","","p_{T}(ll) [GeV/c]","f_{2}");
    plotPFu2frac2.AddFcn(fcnPFu2frac2,kRed);
    plotPFu2frac2.AddGraph(grPFu2frac2,"",kBlack,kOpenCircle);
    plotPFu2frac2.Draw(c,kTRUE,"png");
  }
  
  if(pfu2model>=3) {  
    grPFu2sigma3 = new TGraphErrors(nbins,xval,pfu2Sigma3,xerr,pfu2Sigma3Err);
    grPFu2sigma3->SetName("grPFu2sigma3");
    CPlot plotPFu2sigma3("pfu2sigma3","","p_{T}(ll) [GeV/c]","#sigma_{3} [GeV]");
    plotPFu2sigma3.AddGraph(grPFu2sigma3,"",kBlack,kOpenCircle);
    plotPFu2sigma3.Draw(c,kTRUE,"png");
  
    grPFu2frac3 = new TGraphErrors(nbins,xval,pfu2Frac3, xerr,pfu2Frac3Err);
    grPFu2frac3->SetName("grPFu2frac3");
    CPlot plotPFu2frac3("pfu2frac3","","p_{T}(ll) [GeV/c]","f_{3}");
    plotPFu2frac3.AddGraph(grPFu2frac3,"",kBlack,kOpenCircle);
    plotPFu2frac3.Draw(c,kTRUE,"png");
  }

  //
  // Correlations (N.B. kind of whack and not used)
  //
  cout << endl;
  cout << " Assessing correlations..." << endl;  
  cout << endl;
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    if(isBkgv[ifile]) continue;
    
    infile = new TFile(fnamev[ifile]);
    intree = (TTree*)infile->Get("Events");
  
    intree->SetBranchAddress("runNum",   &runNum);     // event run number
    intree->SetBranchAddress("lumiSec",	 &lumiSec);    // event lumi section
    intree->SetBranchAddress("evtNum",	 &evtNum);     // event number
    intree->SetBranchAddress("matchGen", &matchGen);   // event has both leptons matched to MC Z->ll
    intree->SetBranchAddress("category", &category);   // dilepton category
    intree->SetBranchAddress("npv",	 &npv);        // number of primary vertices
    intree->SetBranchAddress("npu",	 &npu);        // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);     // GEN boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);    // GEN boson phi (signal MC)
    intree->SetBranchAddress("genVy",    &genVy);      // GEN boson rapidity (signal MC)
    intree->SetBranchAddress("genVMass", &genVMass);   // GEN boson mass (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);   // event weight per 1/fb (MC)
    intree->SetBranchAddress("met",	 &met);        // MET
    intree->SetBranchAddress("metPhi",	 &metPhi);     // phi(MET)
    intree->SetBranchAddress("sumEt",	 &sumEt);      // Sum ET
    intree->SetBranchAddress("u1",	 &u1);         // parallel component of recoil
    intree->SetBranchAddress("u2",	 &u2);         // perpendicular component of recoil
    intree->SetBranchAddress("q1",	 &q1);         // charge of tag lepton
    intree->SetBranchAddress("q2",	 &q2);         // charge of probe lepton
    intree->SetBranchAddress("dilep",	 &dilep);      // dilepton 4-vector
    intree->SetBranchAddress("lep1",	 &lep1);       // tag lepton 4-vector
    intree->SetBranchAddress("lep2",	 &lep2);       // probe lepton 4-vector 
    
    for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
    
      if(category!=1 && category!=2)                                 continue;
      if(dilep->M() < MASS_LOW || dilep->M() > MASS_HIGH)            continue;
      if(lep1->Pt()        < PT_CUT  || lep2->Pt()        < PT_CUT)  continue;
      if(fabs(lep1->Eta()) > ETA_CUT || fabs(lep2->Eta()) > ETA_CUT) continue;

      Int_t ipt=-1;
      for(Int_t ibin=0; ibin<nbins; ibin++) {
        if(dilep->Pt() > ptbins[ibin] && dilep->Pt() <= ptbins[ibin+1]) ipt = ibin;
      }
      if(ipt<0) continue;
    
      Double_t zpfu1 = (u1 - pfu1Mean[ipt])/(pfu1Sigma0[ipt]);
      Double_t zpfu2 = (u2 - pfu2Mean[ipt])/(pfu2Sigma0[ipt]);
    
      for(Int_t ibin=0; ibin<ncorrbins; ibin++) {
        if(dilep->Pt() > corrbins[ibin] && dilep->Pt() <= corrbins[ibin+1]) {
          hPFu1u2v[ibin]->Fill(fabs(zpfu1),fabs(zpfu2));	
        }
      }
      if(dilep->Pt() > corrbins[ncorrbins]) {
        hPFu1u2v[ncorrbins-1]->Fill(fabs(zpfu1),fabs(zpfu2));    
      }
    }
  }
  delete infile;
  infile=0, intree=0;

  TH1D *hCorrPFu1u2 = new TH1D("hCorrPFu1u2","",ncorrbins,corrbins);
  for(Int_t ibin=0; ibin<ncorrbins; ibin++) { hCorrPFu1u2->SetBinContent(ibin+1,hPFu1u2v[ibin]->GetCorrelationFactor()); }
  CPlot plotCorrPFu1u2("corrpfu1u2","","p_{T}(ll) [GeV/c]","corr(PF u_{1}, PF u_{2})");
  plotCorrPFu1u2.AddHist1D(hCorrPFu1u2,"");
  plotCorrPFu1u2.Draw(c,kTRUE,"png");
  

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  char outfname[100];
  sprintf(outfname,"%s/fits.root",outputDir.Data());
  TFile *outfile = new TFile(outfname,"RECREATE");
  
  if(grPFu1mean)    grPFu1mean->Write();
  if(grPFu1sigma0)  grPFu1sigma0->Write();
  if(grPFu1sigma1)  grPFu1sigma1->Write();
  if(grPFu1sigma2)  grPFu1sigma2->Write();
  if(grPFu1sigma3)  grPFu1sigma3->Write();
  if(grPFu1frac2)   grPFu1frac2->Write();
  if(grPFu1frac3)   grPFu1frac3->Write();
  if(fcnPFu1mean)   fcnPFu1mean->Write();
  if(fcnPFu1sigma0) fcnPFu1sigma0->Write();
  if(fcnPFu1sigma1) fcnPFu1sigma1->Write();
  if(fcnPFu1sigma2) fcnPFu1sigma2->Write();
  fitresPFu1mean->SetName("fitresPFu1mean");     fitresPFu1mean->Write();
  fitresPFu1sigma0->SetName("fitresPFu1sigma0"); fitresPFu1sigma0->Write();
  fitresPFu1sigma1->SetName("fitresPFu1sigma1"); fitresPFu1sigma1->Write();
  fitresPFu1sigma2->SetName("fitresPFu1sigma2"); fitresPFu1sigma2->Write();

  if(grPFu2mean)    grPFu2mean->Write();
  if(grPFu2sigma0)  grPFu2sigma0->Write();
  if(grPFu2sigma1)  grPFu2sigma1->Write();
  if(grPFu2sigma2)  grPFu2sigma2->Write();
  if(grPFu2sigma3)  grPFu2sigma3->Write();
  if(grPFu2frac2)   grPFu2frac2->Write();
  if(grPFu2frac3)   grPFu2frac3->Write();
  if(fcnPFu2mean)   fcnPFu2mean->Write();
  if(fcnPFu2sigma0) fcnPFu2sigma0->Write();
  if(fcnPFu2sigma1) fcnPFu2sigma1->Write();
  if(fcnPFu2sigma2) fcnPFu2sigma2->Write();
  fitresPFu2mean->SetName("fitresPFu2mean");     fitresPFu2mean->Write();
  fitresPFu2sigma0->SetName("fitresPFu2sigma0"); fitresPFu2sigma0->Write();
  fitresPFu2sigma1->SetName("fitresPFu2sigma1"); fitresPFu2sigma1->Write();
  fitresPFu2sigma2->SetName("fitresPFu2sigma2"); fitresPFu2sigma2->Write();
  
  hCorrPFu1u2->Write();
    
  outfile->Close();
  delete outfile;
  
  makeHTML(outputDir,nbins,pfu1model,pfu2model);
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir,  const Int_t nbins,
              const Int_t pfu1model, const Int_t pfu2model)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/plots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Recoil Fits</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  if(pfu1model==1) {
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1mean.png\"><img src=\"plots/pfu1mean.png\" alt=\"plots/pfu1mean.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1sigma1.png\"><img src=\"plots/pfu1sigma1.png\" alt=\"plots/pfu1sigma1.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
  
  } else if(pfu1model==2) {
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1mean.png\"><img src=\"plots/pfu1mean.png\" alt=\"plots/pfu1mean.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1sigma0.png\"><img src=\"plots/pfu1sigma0.png\" alt=\"plots/pfu1sigma0.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1sigma1.png\"><img src=\"plots/pfu1sigma1.png\" alt=\"plots/pfu1sigma1.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1sigma2.png\"><img src=\"plots/pfu1sigma2.png\" alt=\"plots/pfu1sigma2.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1frac2.png\"><img src=\"plots/pfu1frac2.png\" alt=\"plots/pfu1frac2.png\" width=\"100%\"></a></td>" << endl;    
    htmlfile << "</tr>" << endl;
    
  } else if(pfu1model==3) {
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1mean.png\"><img src=\"plots/pfu1mean.png\" alt=\"plots/pfu1mean.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1sigma0.png\"><img src=\"plots/pfu1sigma0.png\" alt=\"plots/pfu1sigma0.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1sigma1.png\"><img src=\"plots/pfu1sigma1.png\" alt=\"plots/pfu1sigma1.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1sigma2.png\"><img src=\"plots/pfu1sigma2.png\" alt=\"plots/pfu1sigma2.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1frac2.png\"><img src=\"plots/pfu1frac2.png\" alt=\"plots/pfu1frac2.png\" width=\"100%\"></a></td>" << endl;      
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1sigma3.png\"><img src=\"plots/pfu1sigma3.png\" alt=\"plots/pfu1sigma3.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1frac3.png\"><img src=\"plots/pfu1frac3.png\" alt=\"plots/pfu1frac3.png\" width=\"100%\"></a></td>" << endl;      
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "PF u1 fits:";
  htmlfile << " <a target=\"_blank\" href=\"pfu1fits.html\">[linear scale]</a>";
  htmlfile << " <a target=\"_blank\" href=\"pfu1fitslog.html\">[log scale]</a>" << endl;
  htmlfile << "<hr />" << endl;
  
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  if(pfu2model==1) {
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2mean.png\"><img src=\"plots/pfu2mean.png\" alt=\"plots/pfu2mean.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2sigma1.png\"><img src=\"plots/pfu2sigma1.png\" alt=\"plots/pfu2sigma1.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
  
  } else if(pfu2model==2) {
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2mean.png\"><img src=\"plots/pfu2mean.png\" alt=\"plots/pfu2mean.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2sigma0.png\"><img src=\"plots/pfu2sigma0.png\" alt=\"plots/pfu2sigma0.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2sigma1.png\"><img src=\"plots/pfu2sigma1.png\" alt=\"plots/pfu2sigma1.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2sigma2.png\"><img src=\"plots/pfu2sigma2.png\" alt=\"plots/pfu2sigma2.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2frac2.png\"><img src=\"plots/pfu2frac2.png\" alt=\"plots/pfu2frac2.png\" width=\"100%\"></a></td>" << endl;    
    htmlfile << "</tr>" << endl;
    
  } else if(pfu2model==3) {
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2mean.png\"><img src=\"plots/pfu2mean.png\" alt=\"plots/pfu2mean.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2sigma0.png\"><img src=\"plots/pfu2sigma0.png\" alt=\"plots/pfu2sigma0.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2sigma1.png\"><img src=\"plots/pfu2sigma1.png\" alt=\"plots/pfu2sigma1.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2sigma2.png\"><img src=\"plots/pfu2sigma2.png\" alt=\"plots/pfu2sigma2.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2frac2.png\"><img src=\"plots/pfu2frac2.png\" alt=\"plots/pfu2frac2.png\" width=\"100%\"></a></td>" << endl;      
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2sigma3.png\"><img src=\"plots/pfu2sigma3.png\" alt=\"plots/pfu2sigma3.png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2frac3.png\"><img src=\"plots/pfu2frac3.png\" alt=\"plots/pfu2frac3.png\" width=\"100%\"></a></td>" << endl;      
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "PF u2 fits:";
  htmlfile << " <a target=\"_blank\" href=\"pfu2fits.html\">[linear scale]</a>";
  htmlfile << " <a target=\"_blank\" href=\"pfu2fitslog.html\">[log scale]</a>" << endl;
  htmlfile << "<hr />" << endl;

        
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/corrpfu1u2.png\"><img src=\"plots/corrpfu1u2.png\" alt=\"plots/corrpfu1u2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"20%\"></td>" << endl;
  htmlfile << "<td width=\"20%\"></td>" << endl;
  htmlfile << "<td width=\"20%\"></td>" << endl;
  htmlfile << "<td width=\"20%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
  
  Int_t ibin=0;
  
  //
  // PF u1 fits page
  //
  sprintf(htmlfname,"%s/pfu1fits.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;   
  for(ibin=0; ibin<nbins; ibin++) {
    if(ibin%5 == 0) htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1fit_" << ibin << ".png\"><img src=\"plots/pfu1fit_" << ibin << ".png\" alt=\"plots/pfu1fit_" << ibin << ".png\" width=\"100%\"></a></td>" << endl;
    if(ibin%5 == 4) htmlfile << "</tr>" << endl;
  }
  if(ibin%5 != 4) {
    while(ibin%5 != 4) {
      htmlfile << "<td width=\"20%\"></td>" << endl;
      ibin++;
    }
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();  

  sprintf(htmlfname,"%s/pfu1fitslog.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;   
  for(ibin=0; ibin<nbins; ibin++) {
    if(ibin%5 == 0) htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu1fitlog_" << ibin << ".png\"><img src=\"plots/pfu1fitlog_" << ibin << ".png\"alt=\"plots/pfu1fitlog_" << ibin << ".png\" width=\"100%\"></a></td>" << endl;
    if(ibin%5 == 4) htmlfile << "</tr>" << endl;
  }
  if(ibin%5 != 4) {
    while(ibin%5 != 4) {
      htmlfile << "<td width=\"20%\"></td>" << endl;
      ibin++;
    }
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close(); 
  
  //
  // PF u2 fits page
  //
  sprintf(htmlfname,"%s/pfu2fits.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
  for(ibin=0; ibin<nbins; ibin++) {
    if(ibin%5 == 0) htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/pfu2fit_" << ibin << ".png\"><img src=\"plots/pfu2fit_" << ibin << ".png\"alt=\"plots/pfu2fit_" << ibin << ".png\" width=\"100%\"></a></td>" << endl;
    if(ibin%5 == 4) htmlfile << "</tr>" << endl;
  }
  if(ibin%5 != 4) {
    while(ibin%5 != 4) {
      htmlfile << "<td width=\"20%\"></td>" << endl;
      ibin++;
    }
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close(); 

  sprintf(htmlfname,"%s/pfu2fitslog.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
  for(ibin=0; ibin<nbins; ibin++) {
    if(ibin%5 == 0) htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plotslog/pfu2fit_" << ibin << ".png\"><img src=\"plots/pfu2fitlog_" << ibin << ".png\"alt=\"plots/pfu2fitlog_" << ibin << ".png\" width=\"100%\"></a></td>" << endl;
    if(ibin%5 == 4) htmlfile << "</tr>" << endl;
  }
  if(ibin%5 != 4) {
    while(ibin%5 != 4) {
      htmlfile << "<td width=\"20%\"></td>" << endl;
      ibin++;
    }
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close(); 
}

//--------------------------------------------------------------------------------------------------
void performFit(const vector<TH1D*> hv, const vector<TH1D*> hbkgv, const Double_t *ptbins, const Int_t nbins,
                const Int_t model, const Bool_t sigOnly,
                TCanvas *c, const char *plabel, const char *xlabel,
		Double_t *meanArr,   Double_t *meanErrArr,
		Double_t *sigma0Arr, Double_t *sigma0ErrArr,
		Double_t *sigma1Arr, Double_t *sigma1ErrArr,
		Double_t *sigma2Arr, Double_t *sigma2ErrArr,
		Double_t *sigma3Arr, Double_t *sigma3ErrArr,
		Double_t *frac2Arr,  Double_t *frac2ErrArr,
		Double_t *frac3Arr,  Double_t *frac3ErrArr
) {
  char pname[50];
  char ylabel[50];
  char binlabel[50];
  char nsigtext[50];
  char nbkgtext[50];
  
  char meantext[50];
  char sig0text[50];
  char sig1text[50];
  char sig2text[50];
  char sig3text[50];
  
  for(Int_t ibin=0; ibin<nbins; ibin++) {
    RooRealVar u("u","u",hv[ibin]->GetXaxis()->GetXmin(),hv[ibin]->GetXaxis()->GetXmax());
    u.setBins(100);
    RooDataHist dataHist("dataHist","dataHist",RooArgSet(u),hv[ibin]);
    
    //
    // Set up background histogram templates
    //
    RooDataHist bkgHist("bkgHist","bkgHist",RooArgSet(u),hbkgv[ibin]);
    RooHistPdf bkg("bkg","bkg",u,bkgHist,0);
    
    //
    // Set up fit parameters
    //
    RooRealVar mean("mean","mean",
                    hv[ibin]->GetMean(),
                    hv[ibin]->GetXaxis()->GetXmin(),
                    hv[ibin]->GetXaxis()->GetXmax());
    RooRealVar sigma1("sigma1","sigma1",0.8*(hv[ibin]->GetRMS()),0,20);
    RooRealVar sigma2("sigma2","sigma2",1.2*(hv[ibin]->GetRMS()),0,50);
    RooRealVar sigma3("sigma3","sigma3",10.*(hv[ibin]->GetRMS()),0,200); 
    RooRealVar frac2("frac2","frac2",0.5,0.1,0.9);
    RooRealVar frac3("frac3","frac3",0.05,0,0.15);
    RooGaussian gauss1("gauss1","gauss1",u,mean,sigma1);
    RooGaussian gauss2("gauss2","gauss2",u,mean,sigma2);
    RooGaussian gauss3("gauss3","gauss3",u,mean,sigma3);

    if(ibin>0) {
      mean.setVal(meanArr[ibin-1]);
      sigma1.setMin(0.6*(sigma0Arr[ibin-1]));
      sigma1.setMax(1.0*(sigma0Arr[ibin-1]));
      sigma1.setVal(0.8*(sigma0Arr[ibin-1]));
      sigma2.setMin(0.7*(sigma0Arr[ibin-1]));
      sigma2.setMax(1.5*(sigma0Arr[ibin-1]));
      sigma2.setVal(1.2*(sigma0Arr[ibin-1]));
    }
    
    //
    // Define formula for overall width (sigma0)
    //
    char formula[100];
    RooArgList params;
    if(model==1) {
      sprintf(formula,"sigma1");
    
    } else if(model==2) {
      sprintf(formula,"(1.-frac2)*sigma1 + frac2*sigma2");
      params.add(frac2);
      params.add(sigma1);
      params.add(sigma2);

    } else if(model==3) {
      sprintf(formula,"(1.-frac2-frac3)*sigma1 + frac2*sigma2 + frac3*sigma3");
      params.add(frac2);
      params.add(frac3);
      params.add(sigma1);
      params.add(sigma2);
      params.add(sigma3);
    }       
    RooFormulaVar sigma0("sigma0",formula,params);
    
    //
    // Construct fit model
    //
    RooArgList shapes;
    if(model>=3) shapes.add(gauss3);
    if(model>=2) shapes.add(gauss2);
    shapes.add(gauss1);
  
    RooArgList fracs;
    if(model>=3) fracs.add(frac3);
    if(model>=2) fracs.add(frac2);
    
    RooAddPdf sig("sig","sig",shapes,fracs);
    
    RooArgList parts;
    parts.add(sig);
    if(!sigOnly) parts.add(bkg);
    
    RooArgList yields;
    RooRealVar nsig("nsig","nsig",0.98*(hv[ibin]->Integral()),0.,hv[ibin]->Integral());
    yields.add(nsig);
    RooRealVar nbkg("nbkg","nbkg",0.01*(hv[ibin]->Integral()),0.,0.50*(hv[ibin]->Integral()));
    if(!sigOnly) yields.add(nbkg);
    else         nbkg.setVal(0);
    
    RooAddPdf modelpdf("modelpdf","modelpdf",parts,yields);
        
    //
    // Perform fit
    //
    RooFitResult *fitResult=0;
    fitResult = modelpdf.fitTo(dataHist,
                               //RooFit::Minos(),
			       RooFit::Strategy(2),
	                       RooFit::Save());
    
    if(sigma1.getVal() > sigma2.getVal()) {
      Double_t wide = sigma1.getVal();
      Double_t thin = sigma2.getVal();
      sigma1.setVal(thin);
      sigma2.setVal(wide);
      frac2.setVal(1.0-frac2.getVal());
      fitResult = modelpdf.fitTo(dataHist,
                                 //RooFit::Minos(),
			         RooFit::Strategy(2),
	                         RooFit::Save());
    }
    
    meanArr[ibin]      = mean.getVal();
    meanErrArr[ibin]   = mean.getError();
    sigma1Arr[ibin]    = sigma1.getVal();
    sigma1ErrArr[ibin] = sigma1.getError();
    if(model>=2) {
      sigma0Arr[ibin]    = sigma0.getVal();
      sigma0ErrArr[ibin] = sigma0.getPropagatedError(*fitResult);
      sigma2Arr[ibin]    = sigma2.getVal();
      sigma2ErrArr[ibin] = sigma2.getError();
      frac2Arr[ibin]     = frac2.getVal();
      frac2ErrArr[ibin]  = frac2.getError();
    }
    if(model>=3) {    
      sigma3Arr[ibin]    = sigma3.getVal();
      sigma3ErrArr[ibin] = sigma3.getError();
      frac3Arr[ibin]     = frac3.getVal();
      frac3ErrArr[ibin]  = frac3.getError();
    }
    
    //
    // Plot fit results
    //
    RooPlot *frame = u.frame(Bins(100));
    dataHist.plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelpdf.plotOn(frame);
    if(!sigOnly) modelpdf.plotOn(frame,Components("bkg"),LineStyle(kDotted),LineColor(kMagenta+2));
    if(model>=2) sig.plotOn(frame,Components("gauss1"),LineStyle(kDashed),LineColor(kRed));
    if(model>=2) sig.plotOn(frame,Components("gauss2"),LineStyle(kDashed),LineColor(kCyan+2));
    if(model>=3) sig.plotOn(frame,Components("gauss3"),LineStyle(kDashed),LineColor(kGreen+2));
    
    sprintf(pname,"%sfit_%i",plabel,ibin);
    sprintf(ylabel,"Events / %.1f GeV",hv[ibin]->GetBinWidth(1));
    sprintf(binlabel,"%i < p_{T} < %i",(Int_t)ptbins[ibin],(Int_t)ptbins[ibin+1]);    
    if(sigOnly) {
      sprintf(nsigtext,"N_{evts} = %i",(Int_t)hv[ibin]->Integral());
    } else {
      sprintf(nsigtext,"N_{sig} = %.1f #pm %.1f",nsig.getVal(),nsig.getError());
      sprintf(nbkgtext,"N_{bkg} = %.1f #pm %.1f",nbkg.getVal(),nbkg.getError());
    }
    sprintf(meantext,"#mu = %.1f #pm %.1f",meanArr[ibin],meanErrArr[ibin]);
    sprintf(sig1text,"#sigma = %.1f #pm %.1f",sigma1Arr[ibin],sigma1ErrArr[ibin]);
    if(model>=2) {
      sprintf(sig0text,"#sigma = %.1f #pm %.1f",sigma0Arr[ibin],sigma0ErrArr[ibin]);
      sprintf(sig1text,"#sigma_{1} = %.1f #pm %.1f",sigma1Arr[ibin],sigma1ErrArr[ibin]);          
      sprintf(sig2text,"#sigma_{2} = %.1f #pm %.1f",sigma2Arr[ibin],sigma2ErrArr[ibin]);
    }
    if(model>3)
      sprintf(sig3text,"#sigma_{3} = %.1f #pm %.1f",sigma3Arr[ibin],sigma3ErrArr[ibin]);
    
    CPlot plot(pname,frame,"",xlabel,ylabel);
    plot.AddTextBox(binlabel,0.21,0.80,0.51,0.85,0,kBlack,-1);
    if(sigOnly) plot.AddTextBox(nsigtext,0.21,0.78,0.51,0.73,0,kBlack,-1);
    else        plot.AddTextBox(0.21,0.78,0.51,0.68,0,kBlack,-1,2,nsigtext,nbkgtext);
    if(model==1)      plot.AddTextBox(0.70,0.90,0.95,0.80,0,kBlack,-1,2,meantext,sig1text);
    else if(model==2) plot.AddTextBox(0.70,0.90,0.95,0.70,0,kBlack,-1,4,meantext,sig0text,sig1text,sig2text);
    else if(model==3) plot.AddTextBox(0.70,0.90,0.95,0.65,0,kBlack,-1,5,meantext,sig0text,sig1text,sig2text,sig3text);
    plot.Draw(c,kTRUE,"png");
    
    sprintf(pname,"%sfitlog_%i",plabel,ibin);
    plot.SetName(pname);
    plot.SetLogy();
    plot.Draw(c,kTRUE,"png");        
  }
}
