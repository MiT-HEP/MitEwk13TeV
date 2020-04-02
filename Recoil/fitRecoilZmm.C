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
#include <TLatex.h>
#include <TFitResult.h>               // class to handle fit results
#include <TGraphErrors.h>             // graph class
#include "TLorentzVector.h"           // 4-vector class

#include "../Utils/CPlot.hh"          // helper class for plots
#include "../Utils/MitStyleRemix.hh"  // style settings for drawing

#include "../RochesterCorr/RoccoR.cc"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooFormulaVar.h"
#include "RooRealIntegral.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "RooConstVar.h"
#include "RooMinuit.h"

#endif

using namespace RooFit;
using namespace std;

bool doLikelihoodScan=false;
bool doLog=false; // true for data; false for MC
bool doElectron=false;

//=== FUNCTION DECLARATIONS ======================================================================================
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
		Double_t *chi2Arr,   Double_t *chi2ErrArr,
                RooWorkspace *workspace,
		const char *outputname,
		int etaBinCategory, bool do_keys);

//=== MAIN MACRO ================================================================================================= 

void fitRecoilZmm(TString indir="/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/Production_13TeV",  // input ntuple
                  TString fname="data_select.root",
                  Int_t   pfu1model=2,   // u1 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
                  Int_t   pfu2model=2,   // u2 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
                  Bool_t  sigOnly=1,     // signal event only?
                  std::string metVar = "met",
                  std::string metPhiVar = "metPhi",
                  std::string metName = "pf",
                  TString outputDir="./", // output directory
                  Double_t lumi=1,
		  int etaBinCategory=0, // 0 is inclusive, 1 is fabs(eta)<=0.5,  2 is fabs(eta)=[0.5,1], 3 is fabs(eta)>=1
      bool do_keys=0,
      bool do_5TeV=0
) {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
 
  bool doFootprint = false;
  CPlot::sOutDir = outputDir + TString("/plots");
  gSystem->mkdir(outputDir,kTRUE);
  
  // preserving the fine binning at low pT but the higher-pT bins (>75 GeV have been adjusted to be slightly wider)
  Double_t ptbins[] = {0,5.0,25,100,1000}; 
  // Double_t ptbins[] = {0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,90,100,125,150,1000}; 
  Int_t nbins = sizeof(ptbins)/sizeof(Double_t)-1;
  TString infilename = indir+fname;


  // Double_t corrbins[] = { 0, 10, 30, 50 };
  // Int_t ncorrbins = sizeof(corrbins)/sizeof(Double_t)-1;

  TString formulaPFu1mean("pol2");
  TString formulaPFu2mean("pol2");
  TString formulaPFu1meanScale("pol2");
  TString formulaPFu2meanScale("pol2");
  
  vector<TString> fnamev;
  vector<Bool_t> isBkgv;
  fnamev.push_back(infilename); isBkgv.push_back(kFALSE);
  
  // // comment this part out until we have the correct MC Samples for the background
  // // 13 TeV low pu
  // // not doin any of the backgrounds
   if(!do_5TeV) {
     if(doElectron) {
       // fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU_Fixed/Zee/ntuples/top_select.root"); isBkgv.push_back(kTRUE);
       // fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU_Fixed/Zee/ntuples/ewk_select.root"); isBkgv.push_back(kTRUE);
     } else {
       fnamev.push_back(indir+"top1_select.raw.root"); isBkgv.push_back(kTRUE);
       fnamev.push_back(indir+"top2_select.raw.root"); isBkgv.push_back(kTRUE);
       fnamev.push_back(indir+"top3_select.raw.root"); isBkgv.push_back(kTRUE);
       fnamev.push_back(indir+"zxx_select.raw.root"); isBkgv.push_back(kTRUE);
       fnamev.push_back(indir+"wx0_select.raw.root"); isBkgv.push_back(kTRUE);
       fnamev.push_back(indir+"wx1_select.raw.root"); isBkgv.push_back(kTRUE);
       fnamev.push_back(indir+"wx2_select.raw.root"); isBkgv.push_back(kTRUE);
       fnamev.push_back(indir+"ww_select.raw.root"); isBkgv.push_back(kTRUE);
       fnamev.push_back(indir+"wz_select.raw.root"); isBkgv.push_back(kTRUE);
       fnamev.push_back(indir+"zz_select.raw.root"); isBkgv.push_back(kTRUE);

     }
   }
   if(do_5TeV) {
     fnamev.push_back(indir+"top_select.raw.root"); isBkgv.push_back(kTRUE);
     fnamev.push_back(indir+"wx_select.raw.root"); isBkgv.push_back(kTRUE);
     fnamev.push_back(indir+"zxx_select.raw.root"); isBkgv.push_back(kTRUE);
     fnamev.push_back(indir+"ww_select.raw.root"); isBkgv.push_back(kTRUE);
     fnamev.push_back(indir+"zz_select.raw.root"); isBkgv.push_back(kTRUE);
     fnamev.push_back(indir+"wz_select.raw.root"); isBkgv.push_back(kTRUE);
   }
  
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.4;
  Double_t mu_MASS = 0.1057;

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

    // Puppi-PF
    int range=100;
    if(ptbins[ibin]>80) range=125;
    if(ptbins[ibin]>150) range=150;

    sprintf(hname,"hPFu1_%i",ibin);    hPFu1v.push_back(new TH1D(hname,"",100,-range-ptbins[ibin],range-ptbins[ibin]));    hPFu1v[ibin]->Sumw2();
    sprintf(hname,"hPFu1Bkg_%i",ibin); hPFu1Bkgv.push_back(new TH1D(hname,"",100,-range-ptbins[ibin],range-ptbins[ibin]));  hPFu1Bkgv[ibin]->Sumw2();
    
    sprintf(hname,"hPFu2_%i",ibin);    hPFu2v.push_back(new TH1D(hname,"",100,-range,range));    hPFu2v[ibin]->Sumw2();
    sprintf(hname,"hPFu2Bkg_%i",ibin); hPFu2Bkgv.push_back(new TH1D(hname,"",100,-range,range));  hPFu2Bkgv[ibin]->Sumw2();

    std::stringstream name;
    name << "u_" << ibin;

    RooRealVar u1Var(name.str().c_str(),name.str().c_str(), 0, -range-ptbins[ibin], range-ptbins[ibin]);
    RooRealVar u2Var(name.str().c_str(),name.str().c_str(), 0, -range, range);
    u1Var.setBins(100);
    u2Var.setBins(100);

    vu1Var.push_back(u1Var);
    vu2Var.push_back(u2Var);

    sprintf(hname,"hDataSetU1_%i",ibin);  RooDataSet dataSetU1(hname,hname,RooArgSet(u1Var)); lDataSetU1.push_back(dataSetU1);
    sprintf(hname,"hDataSetU2_%i",ibin);  RooDataSet dataSetU2(hname,hname,RooArgSet(u2Var)); lDataSetU2.push_back(dataSetU2);
  }

  TFile *infile = 0;
  TTree *intree = 0;  

  UInt_t  category;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t scale1fb;
  Float_t met, metPhi;
  Int_t   q1, q2;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0, *lep1_raw=0, *lep2_raw=0, *genlep1=0, *genlep2=0;

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);
    intree = (TTree*)infile->Get("Events");
  
    intree->SetBranchAddress("category", &category);   // dilepton category
    intree->SetBranchAddress("genVPt",   &genVPt);     // GEN boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);    // GEN boson phi (signal MC)
    intree->SetBranchAddress("genVy",    &genVy);      // GEN boson rapidity (signal MC)
    intree->SetBranchAddress("genVMass", &genVMass);   // GEN boson mass (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);   // event weight per 1/fb (MC)
    intree->SetBranchAddress(metVar.c_str(),        &met);        // Uncorrected PF MET
    intree->SetBranchAddress(metPhiVar.c_str(),     &metPhi);     // phi(MET)
    intree->SetBranchAddress("q1",	 &q1);         // charge of tag lepton
    intree->SetBranchAddress("q2",	 &q2);         // charge of probe lepton
    intree->SetBranchAddress("dilep",	 &dilep);      // dilepton 4-vector
    intree->SetBranchAddress("lep1",	 &lep1);       // tag lepton 4-vector
    intree->SetBranchAddress("lep2",	 &lep2);       // probe lepton 4-vector 
    intree->SetBranchAddress("genlep1",	 &genlep1);       // tag lepton 4-vector
    intree->SetBranchAddress("genlep2",	 &genlep2);       // probe lepton 4-vector 
    if(doElectron) intree->SetBranchAddress("lep1_raw",         &lep1_raw);       // tag lepton 4-vector
    if(doElectron) intree->SetBranchAddress("lep2_raw",         &lep2_raw);       // probe lepton 4-vector

    TH1D* hGenWeights;
    double totalNorm = 1.0;
    cout << "Hello " << endl;
    if(!infilename.Contains("data_")) {
      hGenWeights = (TH1D*)infile->Get("hGenWeights");
      totalNorm = hGenWeights->Integral();
    }
    //
    // Loop over events
    //
    int iterator=10;
    // for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    for(Int_t ientry=0; ientry<intree->GetEntries(); ientry+=iterator) {
    // for(Int_t ientry=0; ientry<1000; ientry++) {
      intree->GetEntry(ientry);
    
      TLorentzVector mu1, mu2;
      if(doElectron) mu_MASS=0.000511;
      mu1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),mu_MASS);
      mu2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),mu_MASS);
	  
      double SF1=1, SF2=1;
      if(infilename.Contains("data_")) {
        SF1 = rc.kScaleDT(q1, mu1.Pt(), mu1.Eta(), mu1.Phi());//, s=0, m=0);
        SF2 = rc.kScaleDT(q2, mu2.Pt(), mu2.Eta(), mu2.Phi());//s=0, m=0);
      } else if (!doElectron) {
        SF1 = rc.kSpreadMC(q1, mu1.Pt(), mu1.Eta(), mu1.Phi(), genlep1->Pt());//, s=0, m=0);
        SF2 = rc.kSpreadMC(q2, mu2.Pt(), mu2.Eta(), mu2.Phi(), genlep2->Pt());//, s=0, m=0);
        if(genlep1->Pt()==0 && genlep2->Pt()==0) {SF1=1; SF2=1;}
      } else {SF1=1; SF2=1;}
      
      
      mu1*=SF1;
      mu2*=SF2;
	  
      TLorentzVector dl;
      dl=mu1+mu2;
      double mll=dl.M();
      double etall=dl.Rapidity();
      double ptll=dl.Pt();
	  
      if(!isBkgv[ifile]) {
        if(category!=1 && category!=2 && category != 3)                continue;
        if(mll < MASS_LOW || mll > MASS_HIGH)            continue;
      }
      if(mu1.Pt()          < PT_CUT  || mu2.Pt()          < PT_CUT)  continue;
      if(fabs(mu1.Eta()) > ETA_CUT || fabs(mu2.Eta()) > ETA_CUT) continue;

      if(!isBkgv[ifile]) {
        if(!infilename.Contains("data_")) {
          if(etaBinCategory==1 && fabs(genVy)>0.5) continue;
          if(etaBinCategory==2 && (fabs(genVy)<=0.5 || fabs(genVy)>=1 )) continue;
          if(etaBinCategory==3 && fabs(genVy)<1) continue;
        } else {
          if(etaBinCategory==1 && fabs(etall)>0.5) continue;
          if(etaBinCategory==2 && (fabs(etall)<=0.5 || fabs(etall)>=1 )) continue;
          if(etaBinCategory==3 && fabs(etall)<1) continue;
        }
      }
      Int_t ipt=-1;
      for(Int_t ibin=0; ibin<nbins; ibin++) {
        if(ptll > ptbins[ibin] && ptll <= ptbins[ibin+1])
          ipt = ibin;
      }
      if(ipt<0) continue;

      /////////
      /// RECO filling the Zee or the Zmm
      //
      double pU1=0.;
      double pU2=0.;

      TVector2 vLepRaw1, vLepRaw2;
      if(doElectron) {
        vLepRaw1.Set((lep1_raw->Pt())*cos(lep1_raw->Phi()),(lep1_raw->Pt())*sin(lep1_raw->Phi()));
        vLepRaw2.Set((lep2_raw->Pt())*cos(lep2_raw->Phi()),(lep2_raw->Pt())*sin(lep2_raw->Phi()));
      } else {
        vLepRaw1.Set((lep1->Pt())*cos(lep1->Phi()),(lep1->Pt())*sin(lep1->Phi()));
        vLepRaw2.Set((lep2->Pt())*cos(lep2->Phi()),(lep2->Pt())*sin(lep2->Phi()));
      }
      TVector2 vLepCor1((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
      TVector2 vLepCor2((mu2.Pt())*cos(mu2.Phi()),(mu2.Pt())*sin(mu2.Phi()));

      TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
      Double_t corrMetWithLepton = (vMetCorr + vLepRaw1 + vLepRaw2 - vLepCor1 - vLepCor2).Mod();
      Double_t corrMetWithLeptonPhi = (vMetCorr + vLepRaw1 + vLepRaw2 - vLepCor1 - vLepCor2).Phi();
      double pUX  = corrMetWithLepton*cos(corrMetWithLeptonPhi) + dl.Pt()*cos(dl.Phi());
      double pUY  = corrMetWithLepton*sin(corrMetWithLeptonPhi) + dl.Pt()*sin(dl.Phi());
      double pU   = sqrt(pUX*pUX+pUY*pUY);
      double pCos = - (pUX*cos(dl.Phi()) + pUY*sin(dl.Phi()))/pU;
      double pSin =   (pUX*sin(dl.Phi()) - pUY*cos(dl.Phi()))/pU;
      pU1   = pU*pCos; // U1 in data
      pU2   = pU*pSin; // U2 in data
	  
      if(isBkgv[ifile]) {
        hPFu1Bkgv[ipt]->Fill(pU1,scale1fb*lumi/totalNorm);
        hPFu2Bkgv[ipt]->Fill(pU2,scale1fb*lumi/totalNorm);
      } else {
        if(infilename.Contains("data_")) {
          scale1fb = 1.0;
          lumi = 1.0;
        }

        hPFu1v[ipt]->Fill(pU1,scale1fb*lumi/totalNorm);
        hPFu2v[ipt]->Fill(pU2,scale1fb*lumi/totalNorm);

        // this is the dataset for the RooKey
        // clean the under/overflow
        int range=100;
        if(ptbins[ipt]>80) range=125;
        if(ptbins[ipt]>150) range=150;

        if(pU1<(-range-ptbins[ipt])) continue;
        if(pU1>(range-ptbins[ipt])) continue;
        if(pU2<(-range)) continue;
        if(pU2>(range)) continue;

        vu1Var[ipt].setVal(pU1);
        vu2Var[ipt].setVal(pU2);

        lDataSetU1[ipt].add(RooArgSet(vu1Var[ipt])); // need to add the weights
        lDataSetU2[ipt].add(RooArgSet(vu2Var[ipt]));
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
  Double_t pfu1Mean[nbins],   pfu1MeanErr[nbins];
  Double_t pfu1Mean2[nbins],  pfu1Mean2Err[nbins];
  Double_t pfu1Mean3[nbins],  pfu1Mean3Err[nbins];
  Double_t pfu1MeanScale[nbins],   pfu1MeanErrScale[nbins];
  Double_t pfu1Mean2Scale[nbins],  pfu1Mean2ErrScale[nbins];
  Double_t pfu1Mean3Scale[nbins],  pfu1Mean3ErrScale[nbins];
  Double_t pfu1Sigma0[nbins], pfu1Sigma0Err[nbins];
  Double_t pfu1Sigma1[nbins], pfu1Sigma1Err[nbins];
  Double_t pfu1Sigma2[nbins], pfu1Sigma2Err[nbins];
  Double_t pfu1Sigma3[nbins], pfu1Sigma3Err[nbins];
  Double_t pfu1Frac2[nbins],  pfu1Frac2Err[nbins];  
  Double_t pfu1Frac3[nbins],  pfu1Frac3Err[nbins];
  Double_t pfu1chi2[nbins],   pfu1chi2Err[nbins];;

  Double_t pfu2Mean[nbins],   pfu2MeanErr[nbins];
  Double_t pfu2Mean2[nbins],  pfu2Mean2Err[nbins];
  Double_t pfu2Mean3[nbins],  pfu2Mean3Err[nbins];
  Double_t pfu2Sigma0[nbins], pfu2Sigma0Err[nbins];
  Double_t pfu2Sigma1[nbins], pfu2Sigma1Err[nbins];
  Double_t pfu2Sigma2[nbins], pfu2Sigma2Err[nbins];
  Double_t pfu2Sigma3[nbins], pfu2Sigma3Err[nbins];
  Double_t pfu2Frac2[nbins],  pfu2Frac2Err[nbins];  
  Double_t pfu2Frac3[nbins],  pfu2Frac3Err[nbins];
  Double_t pfu2chi2[nbins],   pfu2chi2Err[nbins];
  
  
  TCanvas *c = MakeCanvas("c","c",800,800);

  char outpdfname[50];
  sprintf(outpdfname,"%s/%s.root",outputDir.Data(),"pdfsU1");

  // Fitting PF-MET u1
  performFit(hPFu1v, hPFu1Bkgv, ptbins, nbins, pfu1model, sigOnly,
	     lDataSetU1, vu1Var,
             c, "pfu1", "u_{#parallel} [GeV]",
	     pfu1Mean,   pfu1MeanErr,
	     pfu1Mean2,  pfu1Mean2Err,
	     pfu1Mean3,  pfu1Mean3Err,
	     pfu1Sigma0, pfu1Sigma0Err,
	     pfu1Sigma1, pfu1Sigma1Err,
	     pfu1Sigma2, pfu1Sigma2Err,
	     pfu1Sigma3, pfu1Sigma3Err,
	     pfu1Frac2,  pfu1Frac2Err,
	     pfu1Frac3,  pfu1Frac3Err,
	     pfu1chi2,   pfu1chi2Err,
	     &pdfsU1,
	     outpdfname,
	     etaBinCategory, do_keys);
          
  pdfsU1.writeToFile(outpdfname,kFALSE);

  sprintf(outpdfname,"%s/%s.root",outputDir.Data(),"pdfsU2");

  // Fitting PF-MET u2         
  performFit(hPFu2v, hPFu2Bkgv, ptbins, nbins, pfu2model, sigOnly,
	     lDataSetU2, vu2Var,
             c, "pfu2", "u_{#perp  } [GeV/c]",
	     pfu2Mean,   pfu2MeanErr,
	     pfu2Mean2,  pfu2Mean2Err,
	     pfu2Mean3,  pfu2Mean3Err,
	     pfu2Sigma0, pfu2Sigma0Err,
	     pfu2Sigma1, pfu2Sigma1Err,
	     pfu2Sigma2, pfu2Sigma2Err,
	     pfu2Sigma3, pfu2Sigma3Err,
	     pfu2Frac2,  pfu2Frac2Err,
	     pfu2Frac3,  pfu2Frac3Err,
	     pfu2chi2,   pfu2chi2Err,
	     &pdfsU2,
	     outpdfname,
	     etaBinCategory, do_keys);

  pdfsU2.writeToFile(outpdfname,kFALSE);

  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;

  return;
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
		Double_t *chi2Arr,   Double_t *chi2ErrArr,
		RooWorkspace *wksp,
		const char *outputname,
		int etaBinCategory, bool do_keys
) {
  char lumi[50];
  char pname[50];
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


  double frac2_ini=0;
  double frac3_ini=0;

  double sigma1_ini=0;
  double sigma2_ini=0;
  double sigma3_ini=0;

  for(Int_t ibin=0; ibin<nbins; ibin++) {
    
    TH1D* hvLOG= (TH1D*) hv[ibin]->Clone();
    for(Int_t i=0; i<hv[ibin]->GetNbinsX()+1; i++) {
      if(hv[ibin]->GetBinContent(i)==0) {
	hvLOG->SetBinContent(i,0);
	hvLOG->SetBinError(i,0);
      } else {
	hvLOG->SetBinContent(i,TMath::Log10(hv[ibin]->GetBinContent(i)));
	hvLOG->SetBinError(i,0);
      }
    }

    TH1D* hbkgvLOG= (TH1D*) hbkgv[ibin]->Clone();
    for(Int_t i=0; i<hbkgv[ibin]->GetNbinsX()+1; i++) {
      if(hbkgv[ibin]->GetBinContent(i)==0) {
	hbkgvLOG->SetBinContent(i,0);
	hbkgvLOG->SetBinError(i,0);
      } else {
	hbkgvLOG->SetBinContent(i,TMath::Log10(hbkgv[ibin]->GetBinContent(i)) );
	hbkgvLOG->SetBinError(i,0 );
      }
    }


    std::stringstream name;
    // unfortunately have to give each variable individual names for each bin
    name << "u_" << ibin;
    RooRealVar u(name.str().c_str(),name.str().c_str(),hv[ibin]->GetXaxis()->GetXmin(),hv[ibin]->GetXaxis()->GetXmax());name.str(""); 
    u.setBins(100);

    //    HERE THE LOG histo
    RooDataHist dataHist("dataHist","dataHist",RooArgSet(u),hv[ibin]);
    RooDataHist dataHistLog("dataHistLog","dataHistLog",RooArgSet(u),hvLOG);

    //
    // Set up background histogram templates
    //
    RooDataHist bkgHist("bkgHist","bkgHist",RooArgSet(u),hbkgv[ibin]);
    RooDataHist bkgHistLog("bkgHistLog","bkgHistLog",RooArgSet(u),hbkgvLOG);
    //    RooHistPdf bkg("bkg","bkg",u,bkgHist,0);
    name.str("");  name << "bkg_" << ibin << std::endl;
    RooHistPdf bkg(name.str().c_str(),name.str().c_str(),u,bkgHist,0);name.str("");
    name.str("");  name << "bkgLog_" << ibin << std::endl;
    RooHistPdf bkgLog(name.str().c_str(),name.str().c_str(),u,bkgHistLog,0);name.str("");

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
    RooRealVar sigma1(name.str().c_str(),name.str().c_str(),0.3*(hv[ibin]->GetRMS()),0.1*hv[ibin]->GetRMS(),2.5*(hv[ibin]->GetRMS()));
    name.str(""); name << "sigma2_" << ibin;
    RooRealVar sigma2(name.str().c_str(),name.str().c_str(),1.0*(hv[ibin]->GetRMS()),0.1*hv[ibin]->GetRMS(),4.5*(hv[ibin]->GetRMS()));
    name.str(""); name << "sigma3_" << ibin;
    RooRealVar sigma3(name.str().c_str(),name.str().c_str(),2.0*(hv[ibin]->GetRMS()),0.1*hv[ibin]->GetRMS(),9*hv[ibin]->GetRMS());
    // fraction
    name.str(""); name << "frac2_" << ibin;
    RooRealVar frac2(name.str().c_str(),name.str().c_str(),0.5,0.15,0.85);
    name.str(""); name << "frac3_" << ibin;
    RooRealVar frac3(name.str().c_str(),name.str().c_str(),0.05,0.,0.15);

    if(model==2) {
      frac2.setVal(0.5);
      frac2.setRange(0.,1.);
    }

    if(string(plabel)==string("pfu2")) {
      mean1.setVal(0); mean1.setRange(-5.,5.); mean1.setConstant(kTRUE);
      mean2.setVal(0); mean2.setRange(-5.,5.); mean2.setConstant(kTRUE);
      mean3.setVal(0); mean3.setRange(-5.,5.); mean3.setConstant(kTRUE);
    }

    name.str(""); name << "gauss1_" << ibin;
    RooGaussian gauss1(name.str().c_str(),name.str().c_str(),u,mean1,sigma1);
    name.str(""); name << "gauss2_" << ibin;
    RooGaussian gauss2(name.str().c_str(),name.str().c_str(),u,mean2,sigma2);name.str(""); 
    name.str(""); name << "gauss3_" << ibin;
    RooGaussian gauss3(name.str().c_str(),name.str().c_str(),u,mean3,sigma3);name.str("");
    
    RooGaussian constGauss1("constGauss1","constGauss1",mean1,RooConst(hv[ibin]->GetMean()),RooConst(0.15*hv[ibin]->GetRMS()));
    RooGaussian constGauss2("constGauss2","constGauss2",mean2,RooConst(hv[ibin]->GetMean()),RooConst(0.15*hv[ibin]->GetRMS()));
    RooGaussian constGauss3("constGauss3","constGauss3",mean3,RooConst(0.85*hv[ibin]->GetMean()),RooConst(0.15*hv[ibin]->GetRMS()));
    

    
    //
    // Define formula for overall width (sigma0)
    //
    char formula[100];
    RooArgList params;
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
    //    if(!sigOnly) parts.add(bkg);
    if(!sigOnly && !doLog) parts.add(bkg);
    if(!sigOnly && doLog) parts.add(bkgLog);

    RooArgList yields;

    name.str(""); name << "nsig_" << ibin;
    RooRealVar nsig(name.str().c_str(),name.str().c_str(),0.98*(hv[ibin]->Integral()),0.,1.1*hv[ibin]->Integral()); // just to be sure that doesn't it the boundary
    name.str(""); name << "nbkg_" << ibin;
    RooRealVar nbkg(name.str().c_str(),name.str().c_str(),(hbkgv[ibin]->Integral()),0.,0.25*(hv[ibin]->Integral()));
    RooRealVar *lAbkgFrac =new RooRealVar("AbkgFrac","AbkgFrac",(hv[ibin]->Integral()/(hv[ibin]->Integral()+hbkgv[ibin]->Integral())),0.9,1.0);

    if(sigOnly) {
      yields.add(nsig);  nbkg.setVal(0);
    } else {
      RooFormulaVar * sigbkgFrac= new RooFormulaVar("bkgfrac","@0",RooArgSet(*lAbkgFrac));
      yields.add(*sigbkgFrac);
    }


    name.str("");  name << "modelpdf_" << ibin << std::endl;
    RooAddPdf modelpdf(name.str().c_str(),name.str().c_str(),parts,yields);name.str("");

    RooFitResult *fitResultLOG=0;
    if(doLog) {
      fitResultLOG = modelpdf.fitTo(dataHistLog,
				    NumCPU(4),
				    Minimizer("Minuit2","minimize"),
				    ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				    //				    ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				    RooFit::Minos(),
				    RooFit::Strategy(2),
				    RooFit::Save());
    }


    RooFitResult *fitResult=0;
    fitResult = modelpdf.fitTo(dataHist,
			       NumCPU(4),
			       Minimizer("Minuit2","minimize"),
			       ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
			       //			       ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
			       RooFit::Minos(),
			       RooFit::Strategy(2),
	               RooFit::Save());
    int nTries = 0;

    do {
      fitResult = modelpdf.fitTo(dataHist,
				 NumCPU(4),
				 //				 Minimizer("Minuit2","minimize"),
				 Minimizer("Minuit2","scan"),
				 ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 RooFit::Minos(),
				 RooFit::Strategy(2),
				 RooFit::Save());
                 
        // nTries++;

        fitResult = modelpdf.fitTo(dataHist,
				 NumCPU(4),
				 //				 Minimizer("Minuit2","minimize"),
				 Minimizer("Minuit2","migrad"),
				 ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 RooFit::Hesse(),
				 RooFit::Strategy(2),
				 RooFit::Save());
         fitResult = modelpdf.fitTo(dataHist,
				 NumCPU(4),
				 //				 Minimizer("Minuit2","minimize"),
				 Minimizer("Minuit2","improve"),
				 ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 RooFit::Minos(),
				 RooFit::Strategy(2),
				 RooFit::Save());
				 
         fitResult = modelpdf.fitTo(dataHist,
			       NumCPU(4),
			       Minimizer("Minuit2","minimize"),
			       ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
			       //			       ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
			       RooFit::Minos(),
			       RooFit::Strategy(2),
	               RooFit::Save());
        nTries++;
    }while((fitResult->status()>0 || fitResult->covQual()<3)&&nTries < 10);



  /*
  //// Add the likelihood
  */

    if(doLikelihoodScan && (ibin==5 || ibin==15 || ibin==25) ) {

      TCanvas* c_like = new TCanvas("c_like","c_like",800,800);
      RooAbsReal* nll = modelpdf.createNLL(dataHist,NumCPU(4)) ;
      RooMinuit(*nll).migrad() ;
      // mean1
      RooPlot* frame_mean1 = mean1.frame(Bins(100),Range(hv[ibin]->GetXaxis()->GetXmin()+50,hv[ibin]->GetXaxis()->GetXmax()-50),"mean1");
      if(string(plabel)==string("pfu2")) frame_mean1->SetAxisRange(-5., 5., "X");
      //      RooPlot* frame_mean1 = mean1.frame(Bins(1000),Range(0.5*mean1.getVal(),1.5*mean1.getVal()),"mean1");
      RooAbsReal* nll_mean1 = nll->createProfile(mean1) ;
      nll_mean1->plotOn(frame_mean1,LineColor(kRed)) ; frame_mean1->SetMinimum(0); frame_mean1->SetMaximum(100);

      TString plotname;
      if(string(plabel)==string("pfu1")) plotname=Form("LogLikelihood_mean1_Zpt%.1f_u1.png",ptbins[ibin]);
      if(string(plabel)==string("pfu2")) plotname=Form("LogLikelihood_mean1_Zpt%.1f_u2.png",ptbins[ibin]);
      c_like->Clear();frame_mean1->Draw();c_like->Update();c_like->SaveAs(plotname);

      // mean2
      RooPlot* frame_mean2 = mean2.frame(Bins(100),Range(hv[ibin]->GetXaxis()->GetXmin()+30,hv[ibin]->GetXaxis()->GetXmax()-30),"mean2");
      if(string(plabel)==string("pfu2")) frame_mean2->SetAxisRange(-5., 5., "X");
      //      RooPlot* frame_mean2 = mean2.frame(Bins(1000),Range(0.5*mean2.getVal(),1.5*mean2.getVal()),"mean2");
      RooAbsReal* nll_mean2 = nll->createProfile(mean2) ;
      nll_mean2->plotOn(frame_mean2,LineColor(kRed)) ;   frame_mean2->SetMinimum(0); frame_mean2->SetMaximum(100);

      if(string(plabel)==string("pfu1")) plotname=Form("LogLikelihood_mean2_Zpt%.1f_u1.png",ptbins[ibin]);
      if(string(plabel)==string("pfu2")) plotname=Form("LogLikelihood_mean2_Zpt%.1f_u2.png",ptbins[ibin]);
      c_like->Clear();frame_mean2->Draw();c_like->Update();c_like->SaveAs(plotname);

      // mean3
      RooPlot* frame_mean3 = mean3.frame(Bins(100),Range(hv[ibin]->GetXaxis()->GetXmin()+50,hv[ibin]->GetXaxis()->GetXmax()-50),"mean3");
      if(string(plabel)==string("pfu2")) frame_mean3->SetAxisRange(-5., 5., "X");
      //      RooPlot* frame_mean3 = mean2.frame(Bins(1000),Range(0.5*mean2.getVal(),1.5*mean2.getVal()),"mean3");
      RooAbsReal* nll_mean3 = nll->createProfile(mean3) ;
      nll_mean3->plotOn(frame_mean3,LineColor(kRed)) ;   frame_mean3->SetMinimum(0); frame_mean3->SetMaximum(100);

      if(string(plabel)==string("pfu1")) plotname=Form("LogLikelihood_mean3_Zpt%.1f_u1.png",ptbins[ibin]);
      if(string(plabel)==string("pfu2")) plotname=Form("LogLikelihood_mean3_Zpt%.1f_u2.png",ptbins[ibin]);
      c_like->Clear();frame_mean3->Draw();c_like->Update();c_like->SaveAs(plotname);

      ////
      ////
      ////
      ////

      // sigma 1
      RooPlot* frame_sigma1 = sigma1.frame(Bins(100),Range(0.01*(hv[ibin]->GetRMS()),2.5*(hv[ibin]->GetRMS())),"sigma1");
      //      RooPlot* frame_sigma1 = sigma1.frame(Bins(1000),Range(0,5*sigma1.getVal()),"sigma1");
      RooAbsReal* nll_sigma1 = nll->createProfile(sigma1) ;
      nll_sigma1->plotOn(frame_sigma1,LineColor(kRed)) ; frame_sigma1->SetMinimum(0); frame_sigma1->SetMaximum(100);

      if(string(plabel)==string("pfu1")) plotname=Form("LogLikelihood_sigma1_Zpt%.1f_u1.png",ptbins[ibin]);
      if(string(plabel)==string("pfu2")) plotname=Form("LogLikelihood_sigma1_Zpt%.1f_u2.png",ptbins[ibin]);
      c_like->Clear();frame_sigma1->Draw();c_like->Update();c_like->SaveAs(plotname);

      // sigma 2
      RooPlot* frame_sigma2 = sigma2.frame(Bins(100),Range(0.01*(hv[ibin]->GetRMS()),4.5*(hv[ibin]->GetRMS())),"sigma2");
      //      RooPlot* frame_sigma2 = sigma2.frame(Bins(1000),Range(0,10*sigma2.getVal()),"sigma2");
      RooAbsReal* nll_sigma2 = nll->createProfile(sigma2) ;
      nll_sigma2->plotOn(frame_sigma2,LineColor(kRed)) ;  frame_sigma2->SetMinimum(0); frame_sigma2->SetMaximum(100);

      if(string(plabel)==string("pfu1")) plotname=Form("LogLikelihood_sigma2_Zpt%.1f_u1.png",ptbins[ibin]);
      if(string(plabel)==string("pfu2")) plotname=Form("LogLikelihood_sigma2_Zpt%.1f_u2.png",ptbins[ibin]);
      c_like->Clear();frame_sigma2->Draw();c_like->Update();c_like->SaveAs(plotname);

      // sigma 3
      RooPlot* frame_sigma3 = sigma3.frame(Bins(100),Range(0.5*(hv[ibin]->GetRMS()),9.*(hv[ibin]->GetRMS())),"sigma3");
      //      RooPlot* frame_sigma3 = sigma3.frame(Bins(1000),Range(0,10*sigma3.getVal()),"sigma3");
      RooAbsReal* nll_sigma3 = nll->createProfile(sigma3) ;
      nll_sigma3->plotOn(frame_sigma3,LineColor(kRed)) ; frame_sigma3->SetMinimum(0);  frame_sigma3->SetMaximum(100);

      if(string(plabel)==string("pfu1")) plotname=Form("LogLikelihood_sigma3_Zpt%.1f_u1.png",ptbins[ibin]);
      if(string(plabel)==string("pfu2")) plotname=Form("LogLikelihood_sigma3_Zpt%.1f_u2.png",ptbins[ibin]);
      c_like->Clear();frame_sigma3->Draw();c_like->Update();c_like->SaveAs(plotname);

      ////
      ////
      ////
      ////

      // fraction 2
      RooPlot* frame_frac2 = frac2.frame(Bins(100),Range(0.,1.),"frac2");
      RooAbsReal* nll_frac2 = nll->createProfile(frac2) ;
      nll_frac2->plotOn(frame_frac2,LineColor(kRed)) ;  frame_frac2->SetMinimum(0); frame_frac2->SetMaximum(100);

      if(string(plabel)==string("pfu1")) plotname=Form("LogLikelihood_frac2_Zpt%.1f_u1.png",ptbins[ibin]);
      if(string(plabel)==string("pfu2")) plotname=Form("LogLikelihood_frac2_Zpt%.1f_u2.png",ptbins[ibin]);
      c_like->Clear();frame_frac2->Draw();c_like->Update();c_like->SaveAs(plotname);

      // fraction 3
      RooPlot* frame_frac3 = frac3.frame(Bins(100),Range(0.,1.),"frac3");
      RooAbsReal* nll_frac3 = nll->createProfile(frac3) ;
      nll_frac3->plotOn(frame_frac3,LineColor(kRed)) ;  frame_frac3->SetMinimum(0); frame_frac3->SetMaximum(100);

      if(string(plabel)==string("pfu1")) plotname=Form("LogLikelihood_frac3_Zpt%.1f_u1.png",ptbins[ibin]);
      if(string(plabel)==string("pfu2")) plotname=Form("LogLikelihood_frac3_Zpt%.1f_u2.png",ptbins[ibin]);
      c_like->Clear();frame_frac3->Draw();c_like->Update();c_like->SaveAs(plotname);

    }

  /*
  //// end the likelihood
  */

    char rname[100];
    if(string(plabel)==string("pfu1")) sprintf(rname,"fitResultU1_%i",ibin);
    if(string(plabel)==string("pfu2")) sprintf(rname,"fitResultU2_%i",ibin);

    fitResult->SetName(rname);

    TFile *lFile = TFile::Open(outputname,"UPDATE");
    fitResult->Write();
    lFile->Write();
    lFile->Close();

    c->SetFillColor(kWhite);
    if(fitResult->status()>0) c->SetFillColor(kYellow);

    wksp->import(u);
    wksp->import(modelpdf);
    wksp->import(sig);
    if(!doLog) wksp->import(bkg);
    if(doLog) wksp->import(bkgLog);

    frac2_ini=frac2.getVal();
    frac3_ini=frac3.getVal();
    sigma1_ini=sigma1.getVal();
    sigma2_ini=sigma2.getVal();
    sigma3_ini=sigma3.getVal();


    mean1Arr[ibin]      = mean1.getVal();
    mean1ErrArr[ibin]   = mean1.getError();
    sigma1Arr[ibin]     = sigma1.getVal();
    sigma1ErrArr[ibin]  = sigma1.getError();
    if(model>=2) {
      mean2Arr[ibin]      = mean2.getVal();
      mean2ErrArr[ibin]   = mean2.getError();
      sigma0Arr[ibin]     = sigma0.getVal();
      sigma0ErrArr[ibin]  = sigma0.getPropagatedError(*fitResult);
      sigma2Arr[ibin]     = sigma2.getVal();
      sigma2ErrArr[ibin]  = sigma2.getError();
      frac2Arr[ibin]      = frac2.getVal();
      frac2ErrArr[ibin]   = frac2.getError();
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

    if(!doLog) name.str(""); name << "bkg_" << ibin ;
    if(doLog) name.str(""); name << "bkgLog_" << ibin ;
    if(!sigOnly) modelpdf.plotOn(frame,Components(bkg),FillColor(kRed), DrawOption("F"));
    name.str(""); name << "gauss1_" << ibin ;
    if(model>=2) sig.plotOn(frame,Components(name.str().c_str()),LineStyle(kDashed),LineColor(kRed));
    name.str(""); name << "gauss2_" << ibin ;
    if(model>=2) sig.plotOn(frame,Components(name.str().c_str()),LineStyle(kDashed),LineColor(kMagenta));
    name.str(""); name << "gauss3_" << ibin ;
    if(model>=3) sig.plotOn(frame,Components(name.str().c_str()),LineStyle(kDashed),LineColor(kGreen+2));
    if(!doLog) name.str(""); name << "bkg_" << ibin ;
    if(doLog) name.str(""); name << "bkgLog_" << ibin ;

    // draw the curve
    if(!sigOnly) modelpdf.plotOn(frame,FillColor(kGray),VisualizeError(*fitResult,1),RooFit::Components(modelpdf)); // 1 sigma band
    if(!sigOnly) modelpdf.plotOn(frame,RooFit::LineColor(kGray+2));
    sig.plotOn(frame,FillColor(7),VisualizeError(*fitResult,1),RooFit::Components(sig)); // 1 sigma band
    sig.plotOn(frame,RooFit::LineColor(kBlue));
    
    // redraw the data
    dataHist.plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));

    if(do_keys) {

      //      lDataSet[ibin].Print();
      name.str(""); name << "key_" << ibin;
      //      RooKeysPdf * pdf_keys = new RooKeysPdf(name.str().c_str(),name.str().c_str(), u, dataHist, RooKeysPdf::NoMirror, 2);
      RooKeysPdf pdf_keys(name.str().c_str(),name.str().c_str(),lVar[ibin], lDataSet[ibin], RooKeysPdf::NoMirror, 2);

      RooPlot* xframe  = lVar[ibin].frame(Title(Form("%s Zp_{T}=%.1f - %.1f GeV/c ",plabel,ptbins[ibin],ptbins[ibin+1])));

      lDataSet[ibin].plotOn(xframe) ;
      TCanvas* c = new TCanvas("validatePDF","validatePDF",800,800);
      c->cd();
      pdf_keys.plotOn(xframe,LineColor(kBlue));
      xframe->Draw() ;

      c->SaveAs(Form("%s_%d_dataset.png",plabel,ibin));

      name.str("");

      pdf_keys.plotOn(frame,LineColor(kRed)) ;

      //      wksp->import(lDataSet[ibin]);
      wksp->import(pdf_keys);
      //      wksp->import(lVar[ibin],RooFit::RecycleConflictNodes(),RooFit::Silence());
      //      wksp->import(pdf_keys,RooFit::RecycleConflictNodes(),RooFit::Silence());

      wksp->Print();

    }

    int sizeParam=0;
    if(string(plabel)==string("pfu1")) sizeParam=8; // 3 means + 3 sigma + 2 frac
    if(string(plabel)==string("pfu2")) sizeParam=5; // 0 means + 3 sigma + 2 frac

    TString nameRooHist=Form("h_%s",dataHist.GetName());
    TString nameRooCurve=Form("sig_%d_Norm[u_%d]",ibin,ibin);

    RooHist* hist = frame->getHist(nameRooHist.Data());
    RooCurve* fitCurve = frame->getCurve(nameRooCurve.Data());

    RooHist* hist_pull = hist->makePullHist(*fitCurve);
    hist_pull->SetTitle("");
    hist_pull->GetYaxis()->SetTitle("pull");
    hist_pull->GetYaxis()->SetRangeUser(-5.,5.);
    hist_pull->SetMarkerColor(kAzure);
    hist_pull->SetLineColor(kAzure);
    hist_pull->SetFillColor(kAzure);
    //    hist_pull->GetYaxis()->SetTitleFont(42);
    //    hist_pull->GetXaxis()->SetTitleFont(42);
    hist_pull->GetYaxis()->SetTitleSize  (0.055);
    hist_pull->GetYaxis()->SetTitleOffset(1.600);
    hist_pull->GetYaxis()->SetLabelOffset(0.014);
    hist_pull->GetYaxis()->SetLabelSize  (0.050);
    hist_pull->GetYaxis()->SetLabelFont  (42);
    hist_pull->GetXaxis()->SetTitleSize  (0.055);
    hist_pull->GetXaxis()->SetTitleOffset(1.300);
    hist_pull->GetXaxis()->SetLabelOffset(0.014);
    hist_pull->GetXaxis()->SetLabelSize  (0.050);
    hist_pull->GetXaxis()->SetLabelFont  (42);

    chi2Arr[ibin]  = frame->chiSquare(nameRooCurve.Data(),nameRooHist.Data(),sizeParam);
    chi2ErrArr[ibin]  = 0 ;
    if(chi2Arr[ibin] > 10) { chi2Arr[ibin]=0; chi2ErrArr[ibin]=200; } // just a larger number so that is easy to notice on the plot
    //    cout << " chi2Arr[ibin]=" << chi2Arr[ibin] << " chi2ErrArr[ibin]=" << chi2ErrArr[ibin] << endl;

    // if(do_5TeV) sprintf(lumi,"CMS                               27.4 pb^{-1} (5 TeV)");
    // if(!do_5TeV) sprintf(lumi,"CMS                               2.2 fb^{-1} (13 TeV)");
    sprintf(pname,"%sfit_%i",plabel,ibin);
    sprintf(ylabel,"Events / %.1f GeV",hv[ibin]->GetBinWidth(1));
    sprintf(binlabel,"p_{T}(Z) = %.1f - %.1f GeV/c",ptbins[ibin],ptbins[ibin+1]);

    if(etaBinCategory==1) sprintf(binYlabel,"|y| < 0.5");
    if(etaBinCategory==2) sprintf(binYlabel,"0.5 < |y| < 1");
    if(etaBinCategory==3) sprintf(binYlabel,"|y| > 1");

    if(sigOnly) {
      sprintf(nsigtext,"N_{evts} = %i",(Int_t)hv[ibin]->Integral());
    } else {
      sprintf(nsigtext,"N_{sig}/(N_{bkg}+N_{sig}) = %.3f #pm %.3f",lAbkgFrac->getVal(),lAbkgFrac->getError());
      //      sprintf(nsigtext,"N_{sig} = %.1f #pm %.1f",nsig.getVal(),nsig.getError());
      //      sprintf(nbkgtext,"N_{bkg} = %.1f #pm %.1f",nbkg.getVal(),nbkg.getError());
    }
    sprintf(mean1text,"#mu_{1} = %.1f #pm %.1f",mean1Arr[ibin],mean1ErrArr[ibin]);
    sprintf(sig1text,"#sigma = %.1f #pm %.1f",sigma1Arr[ibin],sigma1ErrArr[ibin]);
    if(model>=2) {
      sprintf(mean2text,"#mu_{2} = %.1f #pm %.1f",mean2Arr[ibin],mean2ErrArr[ibin]);
      //      sprintf(mean2text,"#mu_{2} = #mu_{1} ");
      sprintf(sig0text,"#sigma = %.1f #pm %.1f",sigma0Arr[ibin],sigma0ErrArr[ibin]);
      sprintf(sig1text,"#sigma_{1} = %.1f #pm %.1f",sigma1Arr[ibin],sigma1ErrArr[ibin]);          
      sprintf(sig2text,"#sigma_{2} = %.1f #pm %.1f",sigma2Arr[ibin],sigma2ErrArr[ibin]);
    }
    if(model>=3){
      sprintf(mean3text,"#mu_{3} = %.1f #pm %.1f",mean3Arr[ibin],mean3ErrArr[ibin]);
      sprintf(sig3text,"#sigma_{3} = %.1f #pm %.1f",sigma3Arr[ibin],sigma3ErrArr[ibin]);
    }
    
    ///////////
    /////////// Draw Linear
    ///////////
    ///////////

    TCanvas *cLin = MakeCanvas("cLin","cLin",800,800);
    cLin->Divide(1,2,0,0);
    cLin->cd(1)->SetPad(0,0.3,1.0,1.0);
    cLin->cd(1)->SetTopMargin(0.1);
    cLin->cd(1)->SetBottomMargin(0.01);
    cLin->cd(1)->SetLeftMargin(0.15);
    cLin->cd(1)->SetRightMargin(0.07);
    cLin->cd(1)->SetTickx(1);
    cLin->cd(1)->SetTicky(1);
    cLin->cd(2)->SetPad(0,0,1.0,0.3);
    cLin->cd(2)->SetTopMargin(0.05);
    cLin->cd(2)->SetBottomMargin(0.45);
    cLin->cd(2)->SetLeftMargin(0.15);
    cLin->cd(2)->SetRightMargin(0.07);
    cLin->cd(2)->SetTickx(1);
    cLin->cd(2)->SetTicky(1);


    CPlot plot(pname,frame,"",xlabel,ylabel);
    //    pad1->cd();
    plot.AddTextBox(lumi,0.1,0.92,0.95,0.97,0,kBlack,-1);
    plot.AddTextBox(binlabel,0.21,0.80,0.51,0.85,0,kBlack,-1);
    if(etaBinCategory!=0) plot.AddTextBox(binYlabel,0.21,0.85,0.51,0.9,0,kBlack,-1);
    if(sigOnly) plot.AddTextBox(nsigtext,0.21,0.78,0.51,0.73,0,kBlack,-1);
    //    else        plot.AddTextBox(0.21,0.78,0.51,0.68,0,kBlack,-1,2,nsigtext,nbkgtext);
    else        plot.AddTextBox(nsigtext,0.21,0.78,0.51,0.73,0,kBlack,-1); // this print the fraction now
    if(model==1)      plot.AddTextBox(0.70,0.90,0.95,0.80,0,kBlack,-1,2,mean1text,sig1text);
    else if(model==2) plot.AddTextBox(0.70,0.90,0.95,0.70,0,kBlack,-1,5,mean1text,mean2text,sig0text,sig1text,sig2text);
    //    else if(model==3) plot.AddTextBox(0.70,0.90,0.95,0.65,0,kBlack,-1,7,mean1text,mean2text,mean3text,sig0text,sig1text,sig2text,sig3text);
    else if(model==3) plot.AddTextBox(0.70,0.90,0.95,0.65,0,kBlack,-1,6,mean1text,mean2text,mean3text,sig1text,sig2text,sig3text);
    plot.Draw(cLin,kFALSE,"png",1);
    plot.Draw(cLin,kFALSE,"pdf",1);

    cLin->cd(2);
    hist_pull->Draw("A3 L ");
    TLine *lineZero = new TLine(hv[ibin]->GetXaxis()->GetXmin(), 0,  hv[ibin]->GetXaxis()->GetXmax(), 0);
    lineZero->SetLineColor(kBlack);
    lineZero->Draw("same");
    TLine *lineZero1SigmaM = new TLine(hv[ibin]->GetXaxis()->GetXmin(), 0,  hv[ibin]->GetXaxis()->GetXmax(), 0);
    lineZero1SigmaM->SetLineColor(11);
    lineZero1SigmaM->Draw("same");
    TLine *lineZero1SigmaP = new TLine(hv[ibin]->GetXaxis()->GetXmin(), 0,  hv[ibin]->GetXaxis()->GetXmax(), 0);
    lineZero1SigmaP->SetLineColor(11);
    lineZero1SigmaP->Draw("same");

    plot.Draw(cLin,kTRUE,"png",1);
    plot.Draw(cLin,kTRUE,"pdf",1);

    TCanvas *c1 = MakeCanvas("c1","c1",800,800);
    c1->Divide(1,2,0,0);
    c1->cd(1)->SetPad(0,0.3,1.0,1.0);
    c1->cd(1)->SetTopMargin(0.1);
    c1->cd(1)->SetBottomMargin(0.01);
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetRightMargin(0.07);
    c1->cd(1)->SetTickx(1);
    c1->cd(1)->SetTicky(1);
    c1->cd(2)->SetPad(0,0,1.0,0.3);
    c1->cd(2)->SetTopMargin(0.05);
    c1->cd(2)->SetBottomMargin(0.45);
    c1->cd(2)->SetLeftMargin(0.15);
    c1->cd(2)->SetRightMargin(0.07);
    c1->cd(2)->SetTickx(1);
    c1->cd(2)->SetTicky(1);

    sprintf(pname,"%sfitlog_%i",plabel,ibin);
    plot.SetYRange(0.1,10*hv[ibin]->GetMaximum());
    plot.SetName(pname);
    plot.SetLogy();
    plot.Draw(c1,kFALSE,"png",1);
    plot.Draw(c1,kFALSE,"pdf",1);

    c1->cd(2);
    hist_pull->SetTitle("");
    hist_pull->GetYaxis()->SetTitle("pull");
    hist_pull->GetYaxis()->SetRangeUser(-5.,5.);
    hist_pull->SetMarkerColor(kAzure);
    hist_pull->SetLineColor(kAzure);
    hist_pull->SetFillColor(kAzure);
    hist_pull->Draw("A3 L ");
    lineZero->SetLineColor(kBlack);
    lineZero->Draw("same");
    lineZero1SigmaM->SetLineColor(11);
    lineZero1SigmaM->Draw("same");
    lineZero1SigmaP->SetLineColor(11);
    lineZero1SigmaP->Draw("same");
    plot.Draw(c1,kTRUE,"png",1);
    plot.Draw(c1,kTRUE,"pdf",1);

    // reset color canvas
    c->SetFillColor(kWhite);

  }

  return;

}