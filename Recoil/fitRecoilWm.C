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
#include <sstream>
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TF1.h>                      // 1D function
#include <TFitResult.h>               // class to handle fit results
#include <TGraphErrors.h>             // graph class
#include "TLorentzVector.h"           // 4-vector class

#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"          // helper class for plots
#include "../Utils/MitStyleRemix.hh"  // style settings for drawing

#include <../RochesterCorr/RoccoR.cc>
#include "TRandom.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooFormulaVar.h"
#include "RooRealIntegral.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "RooConstVar.h"
#endif

using namespace RooFit;
using namespace std;

bool doElectron=false;
// bool do5TeV=true;

//=== FUNCTION DECLARATIONS ======================================================================================

// generate web page
// void makeHTML(const TString outDir,  const Int_t nbins, 
              // const Int_t pfu1model, const Int_t pfu2model);
//--------------------------------------------------------------------------------------------------
// perform fit of recoil component
void performFit(const vector<TH1D*> hv, const vector<TH1D*> hbkgv, const Double_t *ptbins, const Int_t nbins,
                const Int_t model, const Bool_t sigOnly,
		const vector<RooDataSet> lDataSet, const vector<RooRealVar> lVar,
                TCanvas *c, const char *plabel, const char *xlabel,
                Double_t *mean1Arr,   Double_t *mean1ErrArr,
                Double_t *mean2Arr,   Double_t *mean2ErrArr,
                Double_t *mean3Arr,   Double_t *mean3ErrArr,
                Double_t *sigma0Arr, Double_t *sigma0ErrArr,
                Double_t *sigma1Arr, Double_t *sigma1ErrArr,
                Double_t *sigma2Arr, Double_t *sigma2ErrArr,
                Double_t *sigma3Arr, Double_t *sigma3ErrArr,
                Double_t *frac2Arr,  Double_t *frac2ErrArr,
                Double_t *frac3Arr,  Double_t *frac3ErrArr,
                RooWorkspace *workspace,
		int etaBinCategory, bool do_keys);

//=== MAIN MACRO ================================================================================================= 

void fitRecoilWm(TString infoldername,  // input ntuple
		 Int_t   pfu1model,     // u1 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
		 Int_t   pfu2model,     // u2 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
		 Bool_t  sigOnly,       // signal event only?
		 Int_t   charge,        // charge requirement
		 Bool_t  useData,       // use Data? (0 = signal MC, 1 = data)
     std::string metVar = "met", // variable name to pull
     std::string metPhiVar = "metPhi", // variable name for MET phi 
		 std::string metName = "pf", // for recordkeeping?
		 TString outputDir ="./",     // output directory
		 Double_t lumi=1, // lumi value for the data fits
		 int etaBinCategory=0, // 0 is inclusive, 1 is fabs(eta)<=0.5,  2 is fabs(eta)=[0.5,1], 3 is fabs(eta)>=1
     bool do_keys=0, // 0 for regular function, 1 for RooKeysPDF
     bool do5TeV=0 // 0 for 13 TeV, 1 for 5 TeV
) {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
 
  
  CPlot::sOutDir = outputDir + TString("/plots");


  // preserving the fine binning at low pT but the higher-pT bins (>75 GeV have been adjusted to be slightly wider)
   // Double_t ptbins[] = {0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,90,100,125,150,1000}; 
   Double_t ptbins[] = {0,25,50,75,100,125,150,1000}; 
   Int_t nbins = sizeof(ptbins)/sizeof(Double_t)-1;

  vector<TString> fnamev;
  vector<Bool_t> isBkgv;

  if (useData == 0){
	  if(doElectron){
      fnamev.push_back(infoldername+"/Wenu/ntuples/we0_select.root"); isBkgv.push_back(kFALSE);
      fnamev.push_back(infoldername+"/Wenu/ntuples/we1_select.root"); isBkgv.push_back(kFALSE);
      fnamev.push_back(infoldername+"/Wenu/ntuples/we2_select.root"); isBkgv.push_back(kFALSE);
	  } else {
      if(!do5TeV){
        fnamev.push_back(infoldername+"/Wmunu/ntuples/wm0_select.raw.root"); isBkgv.push_back(kFALSE);
        fnamev.push_back(infoldername+"/Wmunu/ntuples/wm1_select.raw.root"); isBkgv.push_back(kFALSE);
        fnamev.push_back(infoldername+"/Wmunu/ntuples/wm2_select.raw.root"); isBkgv.push_back(kFALSE);
      } else {
        fnamev.push_back(infoldername+"/Wmunu/ntuples/wm_select.raw.root"); isBkgv.push_back(kFALSE);
      }
    }
  } else if (useData == 1){
      fnamev.push_back(TString(infoldername) + TString("/Wmunu/ntuples/data_select.root")); isBkgv.push_back(kFALSE);
  } else {
    cout << "useData value doesn't make sense" << endl;
  }

  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  const Double_t mu_MASS = 0.1057;
  //Setting up rochester corrections
  RoccoR  rc("../RochesterCorr/RoccoR2017.txt");
 
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  char hname[100];
  vector<TH1D*> hPFu1v,  hPFu1Bkgv;
  vector<TH1D*> hPFu2v,  hPFu2Bkgv;

  vector<RooDataSet> lDataSetU1;
  vector<RooDataSet> lDataSetU2;

  vector<RooRealVar> vu1Var;
  vector<RooRealVar> vu2Var;

  RooWorkspace pdfsU1("pdfsU1");
  RooWorkspace pdfsU2("pdfsU2");

  for(Int_t ibin=0; ibin<nbins; ibin++) {

    int range=100;
    if(ptbins[ibin]>80) range=125;
    if(ptbins[ibin]>150) range=150;
    sprintf(hname,"hPFu1_%i",ibin);    hPFu1v.push_back(new TH1D(hname,"",100,-range-ptbins[ibin],range-ptbins[ibin]));    hPFu1v[ibin]->Sumw2();
    sprintf(hname,"hPFu1Bkg_%i",ibin); hPFu1Bkgv.push_back(new TH1D(hname,"",50,-range-ptbins[ibin],range-ptbins[ibin])); hPFu1Bkgv[ibin]->Sumw2();

    sprintf(hname,"hPFu2_%i",ibin);    hPFu2v.push_back(new TH1D(hname,"",100,-range,range));    hPFu2v[ibin]->Sumw2();
    sprintf(hname,"hPFu2Bkg_%i",ibin); hPFu2Bkgv.push_back(new TH1D(hname,"",50,-range,range)); hPFu2Bkgv[ibin]->Sumw2();

    std::stringstream name;
    name << "u_" << ibin;

    RooRealVar u1Var(name.str().c_str(),name.str().c_str(), 0, -range-ptbins[ibin], range-ptbins[ibin]);
    RooRealVar u2Var(name.str().c_str(),name.str().c_str(), 0, -range, range);

    vu1Var.push_back(u1Var);
    vu2Var.push_back(u2Var);

    sprintf(hname,"hDataSetU1_%i",ibin);  RooDataSet dataSetU1(hname,hname,RooArgSet(u1Var)); lDataSetU1.push_back(dataSetU1);
    sprintf(hname,"hDataSetU2_%i",ibin);  RooDataSet dataSetU2(hname,hname,RooArgSet(u2Var)); lDataSetU2.push_back(dataSetU2);
    
  }
   
  TFile *infile = 0;
  TTree *intree = 0;  
  
  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  // UInt_t  npv, npu;
  Float_t genVPt, genVPhi, genVy;
  Float_t genMuonPt;
  Float_t scale1fb, puWeight;//, scale1fbUp, scale1fbDown;
  Float_t met, metPhi;//, mt, u1, u2;
  Int_t   q;
  UInt_t nTkLayers; // for roch corr
  TLorentzVector *lep=0, *lep_raw=0, *genV = 0, *genLep =0;

//   Float_t puWeight;
//   Float_t scale1fb;

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);
    intree = (TTree*)infile->Get("Events");

    intree->SetBranchAddress("genMuonPt",   &genMuonPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)   
    intree->SetBranchAddress("genVy",    &genVy);     // GEN W boson rapidity (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
    intree->SetBranchAddress(metVar.c_str(),        &met);        // Uncorrected PF MET
    intree->SetBranchAddress(metPhiVar.c_str(),     &metPhi);     // phi(MET)
    intree->SetBranchAddress("q",        &q);         // lepton charge
    intree->SetBranchAddress("lep",      &lep);       // lepton 4-vector 
    intree->SetBranchAddress("nTkLayers",   &nTkLayers);       // lepton 4-vector
	  if(doElectron) intree->SetBranchAddress("lep_raw",         &lep_raw);       // probe lepton 4-vector
    
    TH1D* hGenWeights; double totalNorm = 1.0;
    hGenWeights = (TH1D*)infile->Get("hGenWeights");
    totalNorm = hGenWeights->Integral();
    
    //
    // Loop over events
    //
    int iterator=20;
    // for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    for(Int_t ientry=0; ientry<intree->GetEntries(); ientry+=iterator) {
      intree->GetEntry(ientry);

      // apply rochester correction
      TLorentzVector mu;  mu.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      double SF1=1;
      if(infoldername.Contains("data_")) {
        SF1 = rc.kScaleDT(q, mu.Pt(), mu.Eta(), mu.Phi());
      } else if (!doElectron) {
        if(genMuonPt > 0) SF1 = rc.kSpreadMC(q, mu.Pt(), mu.Eta(), mu.Phi(), genMuonPt);
        else SF1 = rc.kSmearMC(q, mu.Pt(),  mu.Eta(), mu.Phi(), nTkLayers, gRandom->Uniform(1));
      } else { SF1 = 1; } // set to 1 if electrons
      
      mu*=SF1;

      if(charge== 1 && q<0) continue;
      if(charge==-1 && q>0) continue;
      if(mu.Pt()          < PT_CUT)  continue;  
      if(fabs(lep->Eta()) > ETA_CUT) continue;
      
      // 0 is inclusive, 1 is fabs(eta)<=0.5,  2 is fabs(eta)=[0.5,1], 3 is fabs(eta)>=1
      if(etaBinCategory==1 && fabs(genVy)>0.5) continue;
      if(etaBinCategory==2 && (fabs(genVy)<=0.5 || fabs(genVy)>=1 )) continue;
      if(etaBinCategory==3 && fabs(genVy)<1) continue;
    
      Int_t ipt=-1;
      for(Int_t ibin=0; ibin<nbins; ibin++) {
        if(genVPt > ptbins[ibin] && genVPt <= ptbins[ibin+1])
          ipt = ibin;
      }
      if(ipt<0) continue;

    double pU1=0.;
    double pU2=0.;
    TVector2 vLepRaw1;
	  if(doElectron){
		  vLepRaw1.Set((lep_raw->Pt())*cos(lep_raw->Phi()),(lep_raw->Pt())*sin(lep_raw->Phi()));
	  } else {
		  vLepRaw1.Set((lep->Pt())*cos(lep->Phi()),(lep->Pt())*sin(lep->Phi()));
	  }		  
    TVector2 vLepCor1((mu.Pt())*cos(lep->Phi()),(mu.Pt())*sin(lep->Phi()));

    TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
    Double_t corrMetWithLepton = (vMetCorr + vLepRaw1 - vLepCor1).Mod();
    Double_t corrMetWithLeptonPhi = (vMetCorr + vLepRaw1 - vLepCor1).Phi();
    double pUX  = corrMetWithLepton*cos(corrMetWithLeptonPhi) + mu.Pt()*cos(lep->Phi());
    double pUY  = corrMetWithLepton*sin(corrMetWithLeptonPhi) + mu.Pt()*sin(lep->Phi());
    double pU   = sqrt(pUX*pUX+pUY*pUY);
    // projected on the corrected W 
    double pCos = - (pUX*cos(genVPhi) + pUY*sin(genVPhi))/pU;
    double pSin =   (pUX*sin(genVPhi) - pUY*cos(genVPhi))/pU;
    pU1   = pU*pCos; // U1 in data
    pU2   = pU*pSin; // U2 in data

    vu1Var[ipt].setVal(pU1);
    vu2Var[ipt].setVal(pU2);
    lDataSetU1[ipt].add(RooArgSet(vu1Var[ipt])); // need to add the weights
    lDataSetU2[ipt].add(RooArgSet(vu2Var[ipt]));

      if(isBkgv[ifile]) {
        hPFu1Bkgv[ipt]->Fill(pU1,scale1fb*lumi/totalNorm);
        hPFu2Bkgv[ipt]->Fill(pU2,scale1fb*lumi/totalNorm);
      } else {
        hPFu1v[ipt]->Fill(pU1,scale1fb*lumi/totalNorm);
        hPFu2v[ipt]->Fill(pU2,scale1fb*lumi/totalNorm);
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
  

  Double_t pfu1Mean[nbins],   pfu1MeanErr[nbins];
  Double_t pfu1Mean2[nbins],  pfu1Mean2Err[nbins];
  Double_t pfu1Mean3[nbins],  pfu1Mean3Err[nbins];
  Double_t pfu1Sigma0[nbins], pfu1Sigma0Err[nbins];
  Double_t pfu1Sigma1[nbins], pfu1Sigma1Err[nbins];
  Double_t pfu1Sigma2[nbins], pfu1Sigma2Err[nbins];
  Double_t pfu1Sigma3[nbins], pfu1Sigma3Err[nbins];
  Double_t pfu1Frac2[nbins],  pfu1Frac2Err[nbins];  
  Double_t pfu1Frac3[nbins],  pfu1Frac3Err[nbins];

  Double_t pfu2Mean[nbins],   pfu2MeanErr[nbins];
  Double_t pfu2Mean2[nbins],  pfu2Mean2Err[nbins];
  Double_t pfu2Mean3[nbins],  pfu2Mean3Err[nbins];
  Double_t pfu2Sigma0[nbins], pfu2Sigma0Err[nbins];
  Double_t pfu2Sigma1[nbins], pfu2Sigma1Err[nbins];
  Double_t pfu2Sigma2[nbins], pfu2Sigma2Err[nbins];
  Double_t pfu2Sigma3[nbins], pfu2Sigma3Err[nbins];
  Double_t pfu2Frac2[nbins],  pfu2Frac2Err[nbins];  
  Double_t pfu2Frac3[nbins],  pfu2Frac3Err[nbins];
  
  TCanvas *c = MakeCanvas("c","c",800,800);

  // Fitting PF-MET u1
  performFit(hPFu1v, hPFu1Bkgv, ptbins, nbins, pfu1model, sigOnly,
	     lDataSetU1, vu1Var,
             c, "pfu1", "u_{#parallel  } [GeV]",
	     pfu1Mean,   pfu1MeanErr,
	     pfu1Mean2,  pfu1Mean2Err,
	     pfu1Mean3,  pfu1Mean3Err,
	     pfu1Sigma0, pfu1Sigma0Err,
	     pfu1Sigma1, pfu1Sigma1Err,
	     pfu1Sigma2, pfu1Sigma2Err,
	     pfu1Sigma3, pfu1Sigma3Err,
	     pfu1Frac2,  pfu1Frac2Err,
	     pfu1Frac3,  pfu1Frac3Err,
	     &pdfsU1,
	     etaBinCategory, do_keys);

          
  std::cout << "writing" << std::endl;

  char outpdfname[50];

  sprintf(outpdfname,"%s/%s.root",outputDir.Data(),"pdfsU1");
  pdfsU1.writeToFile(outpdfname);
  
  // Fitting PF-MET u2
  performFit(hPFu2v, hPFu2Bkgv, ptbins, nbins, pfu2model, sigOnly,
	     lDataSetU2, vu2Var,
             c, "pfu2", "u_{#perp  } [GeV]",
	     pfu2Mean,   pfu2MeanErr,
	     pfu2Mean2,  pfu2Mean2Err,
	     pfu2Mean3,  pfu2Mean3Err,
	     pfu2Sigma0, pfu2Sigma0Err,
	     pfu2Sigma1, pfu2Sigma1Err,
	     pfu2Sigma2, pfu2Sigma2Err,
	     pfu2Sigma3, pfu2Sigma3Err,
	     pfu2Frac2,  pfu2Frac2Err,
	     pfu2Frac3,  pfu2Frac3Err,
	     &pdfsU2,
	     etaBinCategory, do_keys);
         
  sprintf(outpdfname,"%s/%s.root",outputDir.Data(),"pdfsU2");
  pdfsU2.writeToFile(outpdfname);

  delete infile;
  infile=0, intree=0;

  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;
}

//--------------------------------------------------------------------------------------------------
void performFit(const vector<TH1D*> hv, const vector<TH1D*> hbkgv, const Double_t *ptbins, const Int_t nbins,
                const Int_t model, const Bool_t sigOnly,
		vector<RooDataSet> lDataSet, vector<RooRealVar> lVar,
                TCanvas *c, const char *plabel, const char *xlabel,
		Double_t *mean1Arr,   Double_t *mean1ErrArr,
		Double_t *mean2Arr,   Double_t *mean2ErrArr,
		Double_t *mean3Arr,   Double_t *mean3ErrArr,
		Double_t *sigma0Arr, Double_t *sigma0ErrArr,
		Double_t *sigma1Arr, Double_t *sigma1ErrArr,
		Double_t *sigma2Arr, Double_t *sigma2ErrArr,
		Double_t *sigma3Arr, Double_t *sigma3ErrArr,
		Double_t *frac2Arr,  Double_t *frac2ErrArr,
		Double_t *frac3Arr,  Double_t *frac3ErrArr,
		RooWorkspace *wksp,
		int etaBinCategory, bool do_keys
) {
  char pname[50];
  char lumi[50];
  char ylabel[50];
  char binlabel[50];
  char binYlabel[50];
  char nsigtext[50];
  char nbkgtext[50];
  
  char mean1text[50];
  char mean2text[50];
  char mean3text[50];
  char sig0text[50];
  char sig1text[50];
  char sig2text[50];
  char sig3text[50];
  
  for(Int_t ibin=0; ibin<nbins; ibin++) {
    
    std::stringstream name;
    // unfortunately have to give each variable individual names for each bin
    name << "u_" << ibin;
    RooRealVar u(name.str().c_str(),name.str().c_str(),hv[ibin]->GetXaxis()->GetXmin(),hv[ibin]->GetXaxis()->GetXmax());name.str(""); 
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
    name.str(""); name << "mean1_" << ibin;
    RooRealVar mean1(name.str().c_str(),name.str().c_str(),
                    hv[ibin]->GetMean(),
                    hv[ibin]->GetXaxis()->GetXmin()+50,
                    hv[ibin]->GetXaxis()->GetXmax()-50);
    name.str(""); name << "mean2_" << ibin;
    RooRealVar mean2(name.str().c_str(),name.str().c_str(),
                    hv[ibin]->GetMean(),
                    hv[ibin]->GetXaxis()->GetXmin()+50,
                    hv[ibin]->GetXaxis()->GetXmax()-50);
    name.str(""); name << "mean3_" << ibin;
    RooRealVar mean3(name.str().c_str(),name.str().c_str(),
                    hv[ibin]->GetMean()*0.85,
                    hv[ibin]->GetXaxis()->GetXmin()+50,
                    hv[ibin]->GetXaxis()->GetXmax()-50);
    name.str(""); name << "sigma1_" << ibin;
    RooRealVar sigma1(name.str().c_str(),name.str().c_str(),0.3*(hv[ibin]->GetRMS()),0.,2.3*(hv[ibin]->GetRMS()));
    name.str(""); name << "sigma2_" << ibin;
    RooRealVar sigma2(name.str().c_str(),name.str().c_str(),1.0*(hv[ibin]->GetRMS()),0.,4.5*(hv[ibin]->GetRMS()));
    name.str(""); name << "sigma3_" << ibin;
    RooRealVar sigma3(name.str().c_str(),name.str().c_str(),2.0*(hv[ibin]->GetRMS()),0.,9*(hv[ibin]->GetRMS()));
    name.str(""); name << "frac2_" << ibin;
    RooRealVar frac2(name.str().c_str(),name.str().c_str(),0.5,0.15,0.85);
    name.str(""); name << "frac3_" << ibin;
    RooRealVar frac3(name.str().c_str(),name.str().c_str(),0.05,0,0.15);
    
    if(string(plabel)==string("pfu2")) {

      mean1.setVal(0); mean1.setConstant(kTRUE);
      mean2.setVal(0); mean2.setConstant(kTRUE);
      mean3.setVal(0); mean3.setConstant(kTRUE);

    }

    name.str(""); name << "gauss1_" << ibin;
    RooGaussian gauss1(name.str().c_str(),name.str().c_str(),u,mean1,sigma1);
    name.str(""); name << "gauss2_" << ibin;
    RooGaussian gauss2(name.str().c_str(),name.str().c_str(),u,mean2,sigma2);name.str(""); 
    name.str(""); name << "gauss3_" << ibin;
    RooGaussian gauss3(name.str().c_str(),name.str().c_str(),u,mean3,sigma3);name.str("");
    
    RooGaussian constGauss1("constGauss1","constGauss1",mean1,RooConst(hv[ibin]->GetMean()),RooConst(0.15*hv[ibin]->GetRMS()));
    RooGaussian constGauss2("constGauss2","constGauss2",mean2,RooConst(hv[ibin]->GetMean()),RooConst(0.15*hv[ibin]->GetRMS()));
    RooGaussian constGauss3("constGauss3","constGauss3",mean3,RooConst(hv[ibin]->GetMean()),RooConst(0.15*hv[ibin]->GetRMS()));

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
    
    // sig pdfsU1
    name.str(""); name << "sig_" << ibin;
    RooAddPdf sig(name.str().c_str(),name.str().c_str(),shapes,fracs); name.str(""); 
    
    RooArgList parts;
    parts.add(sig);
    if(!sigOnly) parts.add(bkg);
    
    
    RooArgList yields;
    name.str(""); name << "nsig_" << ibin;
    RooRealVar nsig(name.str().c_str(),name.str().c_str(),0.98*(hv[ibin]->Integral()),0.,1.1*hv[ibin]->Integral()); // just to be sure that doesn't it the boundary
    yields.add(nsig);
    name.str(""); name << "nbkg_" << ibin;
    RooRealVar nbkg(name.str().c_str(),name.str().c_str(),0.01*(hv[ibin]->Integral()),0.,0.50*(hv[ibin]->Integral()));
    if(!sigOnly) yields.add(nbkg);
    else         nbkg.setVal(0);
    
//     std::stringstream name;
    name.str("");  name << "modelpdf_" << ibin << std::endl;
    RooAddPdf modelpdf(name.str().c_str(),name.str().c_str(),parts,yields);name.str(""); 
    
    std::cout << "name = " << name.str().c_str() << std::endl;
    
    std::cout << "on bin # " << ibin << std::endl;
    std::cout << "signal evts = " << nsig.getVal() << std::endl;
    std::cout << "total events = " << hv[ibin]->Integral() << std::endl;
    std::cout << "sigma1 = " << sigma1.getVal() << std::endl;
    std::cout << "sigma2 = " << sigma2.getVal() << std::endl;

    if(ibin>0) std::cout << "sigma max = " << 1.5*hv[ibin-1]->GetRMS() << " " << 1.8*hv[ibin-1]->GetRMS() << std::endl;
    
    //
    // Perform fit
    //
    RooFitResult *fitResult=0;
    fitResult = modelpdf.fitTo(dataHist,
			       NumCPU(4),
			       Minimizer("Minuit2","minimize"),
			       ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                               //RooFit::Minos(),
			       RooFit::Strategy(2),
	                       RooFit::Save());



    int nTries = 0;

    do {
      // if(ibin==22||ibin==27)break;
      fitResult = modelpdf.fitTo(dataHist,
         NumCPU(4),
         Minimizer("Minuit2","scan"),
         ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
         RooFit::Minos(),
         RooFit::Strategy(2),
         RooFit::Save());
                 
        // nTries++;

      fitResult = modelpdf.fitTo(dataHist,
				 NumCPU(4),
				 Minimizer("Minuit2","migrad"),
				 ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 RooFit::Hesse(),
				 RooFit::Strategy(2),
				 RooFit::Save());
         
      fitResult = modelpdf.fitTo(dataHist,
				 NumCPU(4),
				 Minimizer("Minuit2","improve"),
				 ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 RooFit::Minos(),
				 RooFit::Strategy(2),
				 RooFit::Save());
				 
      fitResult = modelpdf.fitTo(dataHist,
         NumCPU(4),
         Minimizer("Minuit2","minimize"),
         ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
         RooFit::Minos(),
         RooFit::Strategy(2),
         RooFit::Save());
         
        nTries++;
    } while((fitResult->status()>0 || fitResult->covQual()<3)&&nTries < 10);

    c->SetFillColor(kWhite);
    if(fitResult->status()>0) c->SetFillColor(kYellow);

    wksp->import(u);
    wksp->import(modelpdf);
    RooRealVar* myXd = (RooRealVar*) wksp->var("u_0");
    
    mean1Arr[ibin]      = mean1.getVal();
    mean1ErrArr[ibin]   = mean1.getError();
    sigma1Arr[ibin]    = sigma1.getVal();
    sigma1ErrArr[ibin] = sigma1.getError();
    if(model>=2) {
      mean2Arr[ibin]      = mean2.getVal();
      mean2ErrArr[ibin]   = mean2.getError();
      sigma0Arr[ibin]    = sigma0.getVal();
      sigma0ErrArr[ibin] = sigma0.getPropagatedError(*fitResult);
      sigma2Arr[ibin]    = sigma2.getVal();
      sigma2ErrArr[ibin] = sigma2.getError();
      frac2Arr[ibin]     = frac2.getVal();
      frac2ErrArr[ibin]  = frac2.getError();
    }
    if(model>=3) {    
      mean3Arr[ibin]      = mean3.getVal();
      mean3ErrArr[ibin]   = mean3.getError();
      sigma3Arr[ibin]    = sigma3.getVal();
      sigma3ErrArr[ibin] = sigma3.getError();
      frac3Arr[ibin]     = frac3.getVal();
      frac3ErrArr[ibin]  = frac3.getError();
    }
    
    std::cout << "Plot Fit results " << std::endl;
    //
    // Plot fit results
    //
    RooPlot *frame = u.frame(Bins(100));
    dataHist.plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelpdf.plotOn(frame);
    if(!sigOnly) modelpdf.plotOn(frame,Components("bkg"),LineStyle(kDotted),LineColor(kMagenta+2));
    name.str(""); name << "gauss1_" << ibin ;
    if(model>=2) sig.plotOn(frame,Components(name.str().c_str()),LineStyle(kDashed),LineColor(kRed));
    name.str(""); name << "gauss2_" << ibin ;
    if(model>=2) sig.plotOn(frame,Components(name.str().c_str()),LineStyle(kDashed),LineColor(kMagenta));
    name.str(""); name << "gauss3_" << ibin ;
    if(model>=3) sig.plotOn(frame,Components(name.str().c_str()),LineStyle(kDashed),LineColor(kGreen+2));

    // draw the curve
    sig.plotOn(frame,FillColor(7),VisualizeError(*fitResult,1),RooFit::Components(sig)); // 1 sigma band
    sig.plotOn(frame,RooFit::LineColor(kBlue));

    // redraw the data
    dataHist.plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    
    if(do_keys) {
      // rookeys

      lDataSet[ibin].Print();
      name.str(""); name << "key_" << ibin ;
      RooKeysPdf pdf_keys(name.str().c_str(),name.str().c_str(),lVar[ibin], lDataSet[ibin], RooKeysPdf::NoMirror, 2);

      RooPlot* xframe  = lVar[ibin].frame(Title(Form("%s Wp_{T}=%d",plabel,ibin))) ;
      lDataSet[ibin].plotOn(xframe) ;
      TCanvas* c = new TCanvas("validatePDF","validatePDF",800,400) ;
      c->cd();
      pdf_keys.plotOn(xframe,LineColor(kBlue)) ;
      xframe->Draw() ;

      c->SaveAs(Form("%s_%d_datasetW.png",plabel,ibin));

      name.str("");

      pdf_keys.plotOn(frame,LineColor(kRed)) ;

      wksp->import(pdf_keys);
      wksp->Print();

    }

    sprintf(lumi,"CMS                               2.3 fb^{-1} (13 TeV)");
    sprintf(pname,"%sfit_%i",plabel,ibin);
    sprintf(ylabel,"Events / %.1f GeV",hv[ibin]->GetBinWidth(1));
    //    sprintf(binlabel,"%i < p_{T}(W) < %i",(Int_t)ptbins[ibin],(Int_t)ptbins[ibin+1]);
    sprintf(binlabel,"p_{T}(W) = %.1f - %.1f GeV ",ptbins[ibin],ptbins[ibin+1]);

    if(etaBinCategory==1) sprintf(binYlabel,"|y| < 0.5");
    if(etaBinCategory==2) sprintf(binYlabel,"0.5 < |y| < 1");
    if(etaBinCategory==3) sprintf(binYlabel,"|y| > 1");
    if(sigOnly) {
      sprintf(nsigtext,"N_{evts} = %i",(Int_t)hv[ibin]->Integral());
    } else {
      sprintf(nsigtext,"N_{sig} = %.1f #pm %.1f",nsig.getVal(),nsig.getError());
      sprintf(nbkgtext,"N_{bkg} = %.1f #pm %.1f",nbkg.getVal(),nbkg.getError());
    }
    sprintf(mean1text,"#mu_{1} = %.1f #pm %.1f",mean1Arr[ibin],mean1ErrArr[ibin]);
    sprintf(sig1text,"#sigma = %.1f #pm %.1f",sigma1Arr[ibin],sigma1ErrArr[ibin]);
    if(model>=2) {
      sprintf(mean2text,"#mu_{2} = %.1f #pm %.1f",mean2Arr[ibin],mean2ErrArr[ibin]);
      // sprintf(sig0text,"#sigma = %.1f #pm %.1f",sigma0Arr[ibin],sigma0ErrArr[ibin]);
      sprintf(sig1text,"#sigma_{1} = %.1f #pm %.1f",sigma1Arr[ibin],sigma1ErrArr[ibin]);          
      sprintf(sig2text,"#sigma_{2} = %.1f #pm %.1f",sigma2Arr[ibin],sigma2ErrArr[ibin]);
    }
    if(model>=3){
      sprintf(mean3text,"#mu_{3} = %.1f #pm %.1f",mean3Arr[ibin],mean3ErrArr[ibin]);
      sprintf(sig3text,"#sigma_{3} = %.1f #pm %.1f",sigma3Arr[ibin],sigma3ErrArr[ibin]);
    }
    
    CPlot plot(pname,frame,"",xlabel,ylabel);
    plot.AddTextBox(lumi,0.05,0.92,0.95,0.97,0,kBlack,-1);
    plot.AddTextBox(binlabel,0.21,0.80,0.51,0.85,0,kBlack,-1);
    if(etaBinCategory!=0) plot.AddTextBox(binYlabel,0.21,0.85,0.51,0.9,0,kBlack,-1);
    if(sigOnly) plot.AddTextBox(nsigtext,0.21,0.78,0.51,0.73,0,kBlack,-1);
    else        plot.AddTextBox(0.21,0.78,0.51,0.68,0,kBlack,-1,2,nsigtext,nbkgtext);
    if(model==1)      plot.AddTextBox(0.70,0.90,0.95,0.80,0,kBlack,-1,2,mean1text,sig1text);
    else if(model==2) plot.AddTextBox(0.70,0.90,0.95,0.70,0,kBlack,-1,5,mean1text,mean2text,sig0text,sig1text,sig2text);
    //    else if(model==3) plot.AddTextBox(0.70,0.90,0.95,0.65,0,kBlack,-1,7,mean1text,mean2text,mean3text,sig0text,sig1text,sig2text,sig3text);
    else if(model==3) plot.AddTextBox(0.70,0.90,0.95,0.65,0,kBlack,-1,6,mean1text,mean2text,mean3text,sig1text,sig2text,sig3text);
    plot.Draw(c,kTRUE,"png");
    plot.Draw(c,kTRUE,"pdf");
    
    sprintf(pname,"%sfitlog_%i",plabel,ibin);
    plot.SetYRange(0.1,10*hv[ibin]->GetMaximum());
    plot.SetName(pname);
    plot.SetLogy();
    plot.Draw(c,kTRUE,"png");        
    plot.Draw(c,kTRUE,"pdf");        

    // reset color canvas
    c->SetFillColor(kWhite);

  }
}