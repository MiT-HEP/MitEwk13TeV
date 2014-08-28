//================================================================================================
//
// Perform fit to extract W->enu signal
//
//  * outputs plots and fit results summary
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>                         // histogram class
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include <TGaxis.h>
#include "Math/LorentzVector.h"           // 4-vector class

#include "../Utils/MyTools.hh"	          // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../Utils/RecoilCorrector2.hh"    // class to handle recoil corrections for MET

// RooFit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

// print correlations of fitted parameters
void printCorrelations(ostream& os, RooFitResult *res);

// print chi2 test and KS test results
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
			   const Double_t ksprob, const Double_t ksprobpe);

// make webpage
void makeHTML(const TString outDir);

// make datacard for Combine
void makeDatacard(const TString outputDir, const TString channel, const TString ewk_option, double ewkSF, double nQCD_exp, double lumi);

Double_t getScaleCorr(const Double_t eta)
{
  if     (fabs(eta) < 0.4)    { return 1.0/1.00243; }
  else if(fabs(eta) < 0.8)    { return 1.0/1.00549; }
  else if(fabs(eta) < 1.2)    { return 1.0/1.00634; }
  else if(fabs(eta) < 1.4442) { return 1.0/1.00669; }
  else if(fabs(eta) < 1.566)  { return 1.0/0.998547; }
  else                        { return 1.0/0.983544; }
}

Double_t getResCorr(const Double_t eta)
{
  if     (fabs(eta) < 0.4)    { return 0.464061;   }
  else if(fabs(eta) < 0.8)    { return 0.00985329; }
  else if(fabs(eta) < 1.2)    { return 0.822958;   }
  else if(fabs(eta) < 1.4442) { return 0.71369;    }
  else if(fabs(eta) < 1.566)  { return 1.35987;    }
  else                        { return 1.27686;    }
}


//=== MAIN MACRO ================================================================================================= 

void fitWe(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)
	   const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
) {
  gBenchmark->Start("fitWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  // MET histogram binning and range
  const Int_t    NBINS  = 50;
  const Double_t METMAX = 100;

  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.5;

  // file format for output plots
  const TString format("png"); 
    
  // recoil correction
  RecoilCorrector recoilCorr("../Recoil/ZeeDataMay23_2/fits.root");//, (!) uncomment to perform corrections to recoil from W-MC/Z-MC
                             //"../Recoil/WepMC_2/fits.root",
			     //"../Recoil/WemMC_2/fits.root",
			     //"../Recoil/ZeeMC_2/fits.root");
  
  // Phil's bias correction
  TFile philCorrFile("/scratch/cmedlock/EWKAnaStore/8TeV/Utils/Scale.root");
  TH1D *hPhilCorr = (TH1D*)philCorrFile.Get("Scale");
   
  // NNLO boson pT k-factors
  TFile nnloCorrFile("/scratch/cmedlock/EWKAnaStore/8TeV/Utils/Ratio.root");
  TH1D *hNNLOCorr = (TH1D*)nnloCorrFile.Get("RpT_B");
  
  //
  // input ntuple file names
  //
  enum { eData, eWenu, eEWK };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;

  // Roughly 4 pileup on average
  fnamev.push_back("/scratch/klawhorn/EWKAnaStore/8TeV/Selection/Wenu/ntuples/data_select.root"); typev.push_back(eData);
  fnamev.push_back("/scratch/klawhorn/EWKAnaStore/8TeV/Selection/Wenu/ntuples/we_select.root");   typev.push_back(eWenu);
  fnamev.push_back("/scratch/klawhorn/EWKAnaStore/8TeV/Selection/Wenu/ntuples/ewk_select.root");  typev.push_back(eEWK);
  fnamev.push_back("/scratch/klawhorn/EWKAnaStore/8TeV/Selection/Wenu/ntuples/top_select.root");  typev.push_back(eEWK);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;
  
  //
  // Declare MET histograms
  //
  TH1D *hDataMet  = new TH1D("hDataMet", "",NBINS,0,METMAX); hDataMet->Sumw2();
  TH1D *hDataMetm = new TH1D("hDataMetm","",NBINS,0,METMAX); hDataMetm->Sumw2();  
  TH1D *hDataMetp = new TH1D("hDataMetp","",NBINS,0,METMAX); hDataMetp->Sumw2();
  TH1D *hWenuMet  = new TH1D("hWenuMet", "",NBINS,0,METMAX); hWenuMet->Sumw2();
  TH1D *hWenuMetp = new TH1D("hWenuMetp","",NBINS,0,METMAX); hWenuMetp->Sumw2();
  TH1D *hWenuMetm = new TH1D("hWenuMetm","",NBINS,0,METMAX); hWenuMetm->Sumw2();
  TH1D *hEWKMet   = new TH1D("hEWKMet",  "",NBINS,0,METMAX); hEWKMet->Sumw2();
  TH1D *hEWKMetp  = new TH1D("hEWKMetp", "",NBINS,0,METMAX); hEWKMetp->Sumw2();
  TH1D *hEWKMetm  = new TH1D("hEWKMetm", "",NBINS,0,METMAX); hEWKMetm->Sumw2();
  TH1D *hQCDMet   = new TH1D("hQCDMet",  "",NBINS,0,METMAX); hQCDMet->Sumw2();
  TH1D *hQCDMetp  = new TH1D("hQCDMetp", "",NBINS,0,METMAX); hQCDMetp->Sumw2();
  TH1D *hQCDMetm  = new TH1D("hQCDMetm", "",NBINS,0,METMAX); hQCDMetm->Sumw2();

  TH1D *hWenuMet_5GeV  = new TH1D("hWenuMet_5GeV", "",NBINS,0,METMAX); hWenuMet_5GeV->Sumw2();
  TH1D *hWenuMetp_5GeV = new TH1D("hWenuMetp_5GeV","",NBINS,0,METMAX); hWenuMetp_5GeV->Sumw2();
  TH1D *hWenuMetm_5GeV = new TH1D("hWenuMetm_5GeV","",NBINS,0,METMAX); hWenuMetm_5GeV->Sumw2();

  TH1D *hWenuMet_10GeV  = new TH1D("hWenuMet_10GeV", "",NBINS,0,METMAX); hWenuMet_10GeV->Sumw2();
  TH1D *hWenuMetp_10GeV = new TH1D("hWenuMetp_10GeV","",NBINS,0,METMAX); hWenuMetp_10GeV->Sumw2();
  TH1D *hWenuMetm_10GeV = new TH1D("hWenuMetm_10GeV","",NBINS,0,METMAX); hWenuMetm_10GeV->Sumw2();

  TH1D *hWenuMet_15GeV  = new TH1D("hWenuMet_15GeV", "",NBINS,0,METMAX); hWenuMet_15GeV->Sumw2();
  TH1D *hWenuMetp_15GeV = new TH1D("hWenuMetp_15GeV","",NBINS,0,METMAX); hWenuMetp_15GeV->Sumw2();
  TH1D *hWenuMetm_15GeV = new TH1D("hWenuMetm_15GeV","",NBINS,0,METMAX); hWenuMetm_15GeV->Sumw2();

  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi;
  Float_t genLepPt, genLepPhi;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  LorentzVector *lep=0;
  LorentzVector *sc=0;
    
  TFile *infile=0;
  TTree *intree=0;

  //
  // Loop over files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);	  assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);

    intree->SetBranchAddress("runNum",   &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    intree->SetBranchAddress("npv",      &npv);       // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);       // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genLepPt",   &genLepPt);    // GEN lepton pT (signal MC)
    intree->SetBranchAddress("genLepPhi",  &genLepPhi);   // GEN lepton phi (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("met",      &met);       // MET
    intree->SetBranchAddress("metPhi",   &metPhi);    // phi(MET)
    intree->SetBranchAddress("sumEt",    &sumEt);     // Sum ET
    intree->SetBranchAddress("mt",       &mt);        // transverse mass
    intree->SetBranchAddress("u1",       &u1);        // parallel component of recoil
    intree->SetBranchAddress("u2",       &u2);        // perpendicular component of recoil
    intree->SetBranchAddress("q",        &q);         // lepton charge
    intree->SetBranchAddress("lep",      &lep);       // lepton 4-vector
    intree->SetBranchAddress("sc",       &sc);        // electron Supercluster 4-vector

    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      
      if(sc->Pt()        < PT_CUT)  continue;	
      if(fabs(sc->Eta()) > ETA_CUT) continue;

      if(typev[ifile]==eData) {
        hDataMet->Fill(met);
	if(q>0) { hDataMetp->Fill(met); } 
	else    { hDataMetm->Fill(met); }
      
      } else {
        Double_t weight = 1;
        weight *= scale1fb*lumi;
	
	if(typev[ifile]==eWenu) {
          Double_t corrMet=met, corrMetPhi=metPhi;
          Double_t corrMet_5GeV=met, corrMet_10GeV=met, corrMet_15GeV=met;
        
	  Double_t philcorr=1;
          for(Int_t ibin=1; ibin<=hPhilCorr->GetNbinsX(); ibin++) {
            if(fabs(sc->Eta()) >= hPhilCorr->GetBinLowEdge(ibin) &&
               fabs(sc->Eta()) < (hPhilCorr->GetBinLowEdge(ibin)+hPhilCorr->GetBinWidth(ibin)))
              philcorr = hPhilCorr->GetBinContent(ibin);
          }

	  // apply recoil corrections to W MC
	  Double_t lepPt = philcorr*(lep->Pt());
	  //Double_t lepPt = philcorr*(gRandom->Gaus((lep->Pt())*getScaleCorr(sc->Eta()),getResCorr(sc->Eta())));  // (!) uncomment to apply scale/res corrections to MC
	  recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),nsigma,q);
          // smear MET
          Double_t METres = 0.0;
          Double_t smearMETx_5GeV=0,  smearMETy_5GeV=0;
          Double_t smearMETx_10GeV=0, smearMETy_10GeV=0;
          Double_t smearMETx_15GeV=0, smearMETy_15GeV=0;
          TVector2 smearMET_5GeV, smearMET_10GeV, smearMET_15GeV;

          smearMETx_5GeV = gRandom->Gaus(corrMet*cos(corrMetPhi),5.0);
          smearMETy_5GeV = gRandom->Gaus(corrMet*sin(corrMetPhi),5.0);
          smearMET_5GeV.Set(smearMETx_5GeV,smearMETy_5GeV);
          corrMet_5GeV  = smearMET_5GeV.Mod();

          smearMETx_10GeV = gRandom->Gaus(corrMet*cos(corrMetPhi),10.0);
          smearMETy_10GeV = gRandom->Gaus(corrMet*sin(corrMetPhi),10.0);
          smearMET_10GeV.Set(smearMETx_10GeV,smearMETy_10GeV);
          corrMet_10GeV = smearMET_10GeV.Mod();

          smearMETx_15GeV = gRandom->Gaus(corrMet*cos(corrMetPhi),15.0);
          smearMETy_15GeV = gRandom->Gaus(corrMet*sin(corrMetPhi),15.0);
          smearMET_15GeV.Set(smearMETx_15GeV,smearMETy_15GeV);
          corrMet_15GeV = smearMET_15GeV.Mod();

          Double_t nnlocorr=1;
          for(Int_t ibin=1; ibin<=hNNLOCorr->GetNbinsX(); ibin++) {
            if(genVPt >= hNNLOCorr->GetBinLowEdge(ibin) &&
               genVPt < (hNNLOCorr->GetBinLowEdge(ibin)+hNNLOCorr->GetBinWidth(ibin)))
              nnlocorr = hNNLOCorr->GetBinContent(ibin);
          }
	  //weight *= nnlocorr;  // (!) uncomment to apply NNLO corrections

	  hWenuMet->Fill(corrMet,weight);
          hWenuMet_5GeV->Fill(corrMet_5GeV,weight);
          hWenuMet_10GeV->Fill(corrMet_10GeV,weight);
          hWenuMet_15GeV->Fill(corrMet_15GeV,weight);
	  if(q>0) { hWenuMetp->Fill(corrMet,weight);} 
	  else    { hWenuMetm->Fill(corrMet,weight); }
	  if(q>0) { hWenuMetp_5GeV->Fill(corrMet_5GeV,weight);} 
	  else    { hWenuMetm_5GeV->Fill(corrMet_5GeV,weight); }
	  if(q>0) { hWenuMetp_10GeV->Fill(corrMet_10GeV,weight);} 
	  else    { hWenuMetm_10GeV->Fill(corrMet_10GeV,weight); }
	  if(q>0) { hWenuMetp_15GeV->Fill(corrMet_15GeV,weight);} 
	  else    { hWenuMetm_15GeV->Fill(corrMet_15GeV,weight); }
        }
        if(typev[ifile]==eEWK) {
          hEWKMet->Fill(met,weight);
          if(q>0) { hEWKMetp->Fill(met,weight); }
          else    { hEWKMetm->Fill(met,weight); }
        }
      }
    }
  }  
  delete infile;
  infile=0, intree=0;   
  
  cout << " hWenuMet->Integral(): " << hWenuMet->Integral() << " and hWenuMet->GetRMS(): " << hWenuMet->GetRMS() << endl;
  cout << "hWenuMetp->Integral(): " << hWenuMetp->Integral() << endl;
  cout << "hWenuMetm->Integral(): " << hWenuMetm->Integral() << endl;

  cout << " hEWKMet->Integral(): " << hEWKMet->Integral() << endl;
  cout << "hEWKMetp->Integral(): " << hEWKMetp->Integral() << endl;
  cout << "hEWKMetm->Integral(): " << hEWKMetm->Integral() << endl;

  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  // Remove this constraint by commenting out the "cewk(p,m).setConstant(kTRUE)" lines
  //

  RooRealVar nSig("nSig","nSig",0.7*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar nQCD("nQCD","nQCD",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar cewk("cewk","cewk",0.1,0,5) ;
  cewk.setVal(hEWKMet->Integral()/hWenuMet->Integral());
  cewk.setConstant(kTRUE);
  RooFormulaVar nEWK("nEWK","nEWK","cewk*nSig",RooArgList(nSig,cewk));

  RooRealVar nSigp("nSigp","nSigp",0.7*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nQCDp("nQCDp","nQCDp",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar cewkp("cewkp","cewkp",0.1,0,5) ;
  cewkp.setVal(hEWKMetp->Integral()/hWenuMetp->Integral());
  cewkp.setConstant(kTRUE);
  RooFormulaVar nEWKp("nEWKp","nEWKp","cewkp*nSigp",RooArgList(nSigp,cewkp));

  RooRealVar nSigm("nSigm","nSigm",0.7*(hDataMetm->Integral()),0,hDataMetm->Integral());
  RooRealVar nQCDm("nQCDm","nQCDm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  RooRealVar cewkm("cewkm","cewkm",0.1,0,5) ;
  cewkm.setVal(hEWKMetm->Integral()/hWenuMetm->Integral());
  cewkm.setConstant(kTRUE);
  RooFormulaVar nEWKm("nEWKm","nEWKm","cewkm*nSigm",RooArgList(nSigm,cewkm));

  //
  // Construct PDFs for fitting
  //
  RooRealVar pfmet("pfmet","pfmet",0,METMAX);
  pfmet.setBins(NBINS);
   
  // Signal PDFs
  RooDataHist wenuMet ("wenuMET", "wenuMET", RooArgSet(pfmet),hWenuMet);  RooHistPdf pdfWe ("we", "we", pfmet,wenuMet, 1);
  RooDataHist wenuMetp("wenuMETp","wenuMETp",RooArgSet(pfmet),hWenuMetp); RooHistPdf pdfWep("wep","wep",pfmet,wenuMetp,1);
  RooDataHist wenuMetm("wenuMETm","wenuMETm",RooArgSet(pfmet),hWenuMetm); RooHistPdf pdfWem("wem","wem",pfmet,wenuMetm,1); 

  RooDataHist wenuMet_5GeV ("wenuMET_5GeV", "wenuMET_5GeV", RooArgSet(pfmet),hWenuMet_5GeV);
  RooDataHist wenuMetp_5GeV("wenuMETp_5GeV","wenuMETp_5GeV",RooArgSet(pfmet),hWenuMetp_5GeV);
  RooDataHist wenuMetm_5GeV("wenuMETm_5GeV","wenuMETm_5GeV",RooArgSet(pfmet),hWenuMetm_5GeV);
  RooHistPdf pdfWe_5GeV ("we_5GeV", "we_5GeV", pfmet,wenuMet_5GeV, 1);
  RooHistPdf pdfWep_5GeV("wep_5GeV","wep_5GeV",pfmet,wenuMetp_5GeV,1);
  RooHistPdf pdfWem_5GeV("wem_5GeV","wem_5GeV",pfmet,wenuMetm_5GeV,1); 

  RooDataHist wenuMet_10GeV ("wenuMET_10GeV", "wenuMET_10GeV", RooArgSet(pfmet),hWenuMet_10GeV);
  RooDataHist wenuMetp_10GeV("wenuMETp_10GeV","wenuMETp_10GeV",RooArgSet(pfmet),hWenuMetp_10GeV);
  RooDataHist wenuMetm_10GeV("wenuMETm_10GeV","wenuMETm_10GeV",RooArgSet(pfmet),hWenuMetm_10GeV);
  RooHistPdf pdfWe_10GeV ("we_10GeV", "we_10GeV", pfmet,wenuMet_10GeV, 1);
  RooHistPdf pdfWep_10GeV("wep_10GeV","wep_10GeV",pfmet,wenuMetp_10GeV,1);
  RooHistPdf pdfWem_10GeV("wem_10GeV","wem_10GeV",pfmet,wenuMetm_10GeV,1); 

  RooDataHist wenuMet_15GeV ("wenuMET_15GeV", "wenuMET_15GeV", RooArgSet(pfmet),hWenuMet_15GeV);
  RooDataHist wenuMetp_15GeV("wenuMETp_15GeV","wenuMETp_15GeV",RooArgSet(pfmet),hWenuMetp_15GeV);
  RooDataHist wenuMetm_15GeV("wenuMETm_15GeV","wenuMETm_15GeV",RooArgSet(pfmet),hWenuMetm_15GeV);
  RooHistPdf pdfWe_15GeV ("we_15GeV", "we_15GeV", pfmet,wenuMet_15GeV, 1);
  RooHistPdf pdfWep_15GeV("wep_15GeV","wep_15GeV",pfmet,wenuMetp_15GeV,1);
  RooHistPdf pdfWem_15GeV("wem_15GeV","wem_15GeV",pfmet,wenuMetm_15GeV,1); 
  
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", pfmet,ewkMet, 1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",pfmet,ewkMetp,1); 
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",pfmet,ewkMetm,1); 
  
  // QCD Pdfs
  CPepeModel1 qcd("qcd",pfmet);
  CPepeModel1 qcdp("qcdp",pfmet);
  CPepeModel1 qcdm("qcdm",pfmet);
/*
  for (int jbin=0;jbin<NBINS;jbin++) {
    hQCDMet->SetBinContent(jbin,1);
    hQCDMetp->SetBinContent(jbin,1);
    hQCDMetm->SetBinContent(jbin,1);
  }

  RooDataHist qcdMet ("qcdMET", "qcdMET", RooArgSet(pfmet),hQCDMet);  RooHistPdf qcd("pepe1Pdf_qcd","pepe1Pdf_qcd",pfmet,qcdMet,1);
  RooDataHist qcdMetp("qcdMETp","qcdMETp",RooArgSet(pfmet),hQCDMetp); RooHistPdf qcdp("pepe1Pdf_qcdp","pepe1Pdf_qcdp",pfmet,qcdMetp,1);
  RooDataHist qcdMetm("qcdMETm","qcdMETm",RooArgSet(pfmet),hQCDMetm); RooHistPdf qcdm("pepe1Pdf_qcdm","pepe1Pdf_qcdm",pfmet,qcdMetm,1);
*/  
  // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWe,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWep,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWem,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));

  // For Combine

  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet),hDataMet);
  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(pfmet),hDataMetp);
  RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(pfmet),hDataMetm);

  RooWorkspace combine_workspace("combine_workspace");

  combine_workspace.import(dataMet);
  combine_workspace.import(dataMetp);
  combine_workspace.import(dataMetm);

  combine_workspace.import(pdfWe);  combine_workspace.import(pdfWe_5GeV); combine_workspace.import(pdfWe_10GeV); combine_workspace.import(pdfWe_15GeV);
  combine_workspace.import(pdfWep);combine_workspace.import(pdfWep_5GeV);combine_workspace.import(pdfWep_10GeV);combine_workspace.import(pdfWep_15GeV);
  combine_workspace.import(pdfWem);combine_workspace.import(pdfWem_5GeV);combine_workspace.import(pdfWem_10GeV);combine_workspace.import(pdfWem_15GeV);
  combine_workspace.import(pdfEWK);
  combine_workspace.import(pdfEWKp);
  combine_workspace.import(pdfEWKm);
  combine_workspace.import(*(qcd.model));
  combine_workspace.import(*(qcdp.model));
  combine_workspace.import(*(qcdm.model));

  combine_workspace.writeToFile("Wenu_pdfTemplates.root");

  RooAbsPdf *we       = combine_workspace.pdf("we");
  RooAbsPdf *we_5GeV  = combine_workspace.pdf("we_5GeV");
  RooAbsPdf *we_10GeV = combine_workspace.pdf("we_10GeV");
  RooAbsPdf *we_15GeV = combine_workspace.pdf("we_15GeV");
  RooPlot *plot = combine_workspace.var("pfmet")->frame();  plot->SetTitle("Wenu - MC Shape Templates");
  we->plotOn(plot,RooFit::LineColor(1));
  we_5GeV->plotOn(plot,RooFit::LineColor(2));
  we_10GeV->plotOn(plot,RooFit::LineColor(3));
  we_15GeV->plotOn(plot,RooFit::LineColor(4));
  TCanvas *wecanvas = new TCanvas();  wecanvas->cd();  plot->Draw();
  wecanvas->Print("/home/cmedlock/public_html/WeShapePlot.png");

  RooAbsPdf *wep = combine_workspace.pdf("wep");
  RooAbsPdf *wep_5GeV  = combine_workspace.pdf("wep_5GeV");
  RooAbsPdf *wep_10GeV = combine_workspace.pdf("wep_10GeV");
  RooAbsPdf *wep_15GeV = combine_workspace.pdf("wep_15GeV");
  RooPlot *plotp = combine_workspace.var("pfmet")->frame();  plotp->SetTitle("Wenu_p - MC Shape Templates");
  wep->plotOn(plotp,RooFit::LineColor(1));
  wep_5GeV->plotOn(plotp,RooFit::LineColor(2));
  wep_10GeV->plotOn(plotp,RooFit::LineColor(3));
  wep_15GeV->plotOn(plotp,RooFit::LineColor(4));
  TCanvas *wepcanvas = new TCanvas();  wepcanvas->cd();  plotp->Draw();
  wepcanvas->Print("/home/cmedlock/public_html/WepShapePlot.png");

  RooAbsPdf *wem = combine_workspace.pdf("wem");
  RooAbsPdf *wem_5GeV  = combine_workspace.pdf("wem_5GeV");
  RooAbsPdf *wem_10GeV = combine_workspace.pdf("wem_10GeV");
  RooAbsPdf *wem_15GeV = combine_workspace.pdf("wem_15GeV");
  RooPlot *plotm = combine_workspace.var("pfmet")->frame();  plotm->SetTitle("Wenu_m - MC Shape Templates");
  wem->plotOn(plotm,RooFit::LineColor(1));
  wem_5GeV->plotOn(plotm,RooFit::LineColor(2));
  wem_10GeV->plotOn(plotm,RooFit::LineColor(3));
  wem_15GeV->plotOn(plotm,RooFit::LineColor(4));
  TCanvas *wemcanvas = new TCanvas();  wemcanvas->cd();  plotm->Draw();
  wemcanvas->Print("/home/cmedlock/public_html/WemShapePlot.png");

  /*
  // Make datacards

  makeDatacard(outputDir,"Wenu",  "tied",cewk.getVal(), nQCD.getVal(), lumi);
  makeDatacard(outputDir,"Wenu_p","tied",cewkp.getVal(),nQCDp.getVal(),lumi);
  makeDatacard(outputDir,"Wenu_m","tied",cewkm.getVal(),nQCDm.getVal(),lumi);
  */
/*
  //
  // Perform fits
  //
//  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet),hDataMet);
  RooFitResult *fitRes = pdfMet.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE));
  
//  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(pfmet),hDataMetp);
  RooFitResult *fitResp = pdfMetp.fitTo(dataMetp,Extended(),Minos(kTRUE),Save(kTRUE));
  
//  RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(pfmet),hDataMetm);
  RooFitResult *fitResm = pdfMetm.fitTo(dataMetm,Extended(),Minos(kTRUE),Save(kTRUE));
    
  //
  // Use histogram version of fitted PDFs to make ratio plots
  // (Will also use PDF histograms later for Chi^2 and KS tests)
  //
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", pfmet));
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  hMetDiff->SetMarkerStyle(kFullCircle);
  hMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", pfmet));
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet));
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle); 
  hMetmDiff->SetMarkerSize(0.9);
   
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  
  TCanvas *c = MakeCanvas("c","c",800,800);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.01);
  c->cd(1)->SetLeftMargin(0.15);  
  c->cd(1)->SetRightMargin(0.07);  
  c->cd(1)->SetTickx(1);
  c->cd(1)->SetTicky(1);  
  c->cd(2)->SetPad(0,0,1.0,0.3);
  c->cd(2)->SetTopMargin(0.05);
  c->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.15);
  c->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetTickx(1);
  c->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.100,"Y");
  TGaxis::SetMaxDigits(3);
  
  char ylabel[100];  // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  if(lumi<0.1) sprintf(lumitext,"L = %.1f pb^{-1},  #sqrt{s} = 8 TeV",lumi*1000.);
  else         sprintf(lumitext,"L = %.2f fb^{-1},  #sqrt{s} = 8 TeV",lumi);
  
  // plot colors
  Int_t linecolorW   = kOrange-3;
  Int_t fillcolorW   = kOrange-2;
  Int_t linecolorEWK = kOrange+10;
  Int_t fillcolorEWK = kOrange+7;
  Int_t linecolorQCD = kViolet+2;
  Int_t fillcolorQCD = kViolet-5;
  Int_t ratioColor   = kGray+2;
  
  //
  // Dummy histograms for TLegend
  // (I can't figure out how to properly pass RooFit objects...)
  //
  TH1D *hDummyData = new TH1D("hDummyData","",0,0,10);
  hDummyData->SetMarkerStyle(kFullCircle);
  hDummyData->SetMarkerSize(0.9);
  
  TH1D *hDummyW = new TH1D("hDummyW","",0,0,10);
  hDummyW->SetLineColor(linecolorW);
  hDummyW->SetFillColor(fillcolorW);
  hDummyW->SetFillStyle(1001);
  
  TH1D *hDummyEWK = new TH1D("hDummyEWK","",0,0,10);
  hDummyEWK->SetLineColor(linecolorEWK);
  hDummyEWK->SetFillColor(fillcolorEWK);
  hDummyEWK->SetFillStyle(1001);
  
  TH1D *hDummyQCD = new TH1D("hDummyQCD","",0,0,10);
  hDummyQCD->SetLineColor(linecolorQCD);
  hDummyQCD->SetFillColor(fillcolorQCD);
  hDummyQCD->SetFillStyle(1001);
   
  //
  // W MET plot
  //
  RooPlot *weframe = pfmet.frame(Bins(NBINS));
  weframe->GetYaxis()->SetNdivisions(505);
  dataMet.plotOn(weframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMet.plotOn(weframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMet.plotOn(weframe,LineColor(linecolorW));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfEWK,*(qcd.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfEWK,*(qcd.model))),LineColor(linecolorEWK));
  pdfMet.plotOn(weframe,Components(RooArgSet(*(qcd.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMet.plotOn(weframe,Components(RooArgSet(*(qcd.model))),LineColor(linecolorQCD));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfWe)),LineColor(linecolorW),LineStyle(2));
  dataMet.plotOn(weframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMet->GetBinWidth(1));
  CPlot plotMet("fitmet",weframe,"","",ylabel);
  plotMet.SetLegend(0.68,0.57,0.93,0.77);
  plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMet.GetLegend()->AddEntry(hDummyW,"W#rightarrowe#nu","F");
  plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMet.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotMet.AddTextBox("CMS",0.82,0.91,0.92,0.98,0);
//  plotMet.SetYRange(0.1,1.1*(hDataMet->GetMaximum()));
plotMet.SetYRange(0.1,1e4);
  plotMet.Draw(c,kTRUE,format,1);

  CPlot plotMetDiff("fitmet","","#slash{E}_{T} [GeV]","#chi");
  plotMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
  plotMetDiff.SetYRange(-8,8);
  plotMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotMetDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotMetDiff.Draw(c,kTRUE,format,2);
  
  plotMet.SetName("fitmetlog");
  plotMet.SetLogy();
  plotMet.SetYRange(1e-3*(hDataMet->GetMaximum()),10*(hDataMet->GetMaximum()));
  plotMet.Draw(c,kTRUE,format,1);
    
  //
  // W+ MET plot
  //
  RooPlot *wepframe = pfmet.frame(Bins(NBINS));    
  wepframe->GetYaxis()->SetNdivisions(505);
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetp.plotOn(wepframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetp.plotOn(wepframe,LineColor(linecolorW));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),LineColor(linecolorEWK));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(*(qcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(*(qcdp.model))),LineColor(linecolorQCD));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfWep)),LineColor(linecolorW),LineStyle(2));
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
  CPlot plotMetp("fitmetp",wepframe,"","",ylabel);
  plotMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrowe^{+}#nu","F");
  plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetp.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotMetp.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
//  plotMetp.SetYRange(0.1,1.1*(hDataMetp->GetMaximum()));
plotMetp.SetYRange(0.1,5000);
  plotMetp.Draw(c,kFALSE,format,1);

  CPlot plotMetpDiff("fitmetp","","#slash{E}_{T} [GeV]","#chi");
  plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
  plotMetpDiff.SetYRange(-8,8);
  plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetpDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotMetpDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotMetpDiff.Draw(c,kTRUE,format,2);
  
  plotMetp.SetName("fitmetplog");
  plotMetp.SetLogy();
  plotMetp.SetYRange(1e-3*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
  plotMetp.Draw(c,kTRUE,format,1);
  
  //
  // W- MET plot
  //
  RooPlot *wemframe = pfmet.frame(Bins(NBINS)); 
  wemframe->GetYaxis()->SetNdivisions(505);
  dataMetm.plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetm.plotOn(wemframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetm.plotOn(wemframe,LineColor(linecolorW));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),LineColor(linecolorEWK));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(*(qcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(*(qcdm.model))),LineColor(linecolorQCD));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfWem)),LineColor(linecolorW),LineStyle(2));
  dataMetm.plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
  CPlot plotMetm("fitmetm",wemframe,"","",ylabel);
  plotMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrowe^{-}#bar{#nu}","F");
  plotMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetm.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotMetm.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
//  plotMetm.SetYRange(0.1,1.1*(hDataMetm->GetMaximum()));
plotMetm.SetYRange(0.1,5000);
  plotMetm.Draw(c,kFALSE,format,1);

  CPlot plotMetmDiff("fitmetm","","#slash{E}_{T} [GeV]","#chi");
  plotMetmDiff.AddHist1D(hMetmDiff,"EX0",ratioColor);
  plotMetmDiff.SetYRange(-8,8);
  plotMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetmDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotMetmDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotMetmDiff.Draw(c,kTRUE,format,2);
  
  plotMetm.SetName("fitmetmlog");
  plotMetm.SetLogy();
  plotMetm.SetYRange(1e-3*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
  plotMetm.Draw(c,kTRUE,format,1);

    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  
  //
  // Write fit results
  //
  ofstream txtfile;
  char txtfname[100];    
  
  ios_base::fmtflags flags;
  
  Double_t chi2prob, chi2ndf;
  Double_t ksprob, ksprobpe;
  
  chi2prob = hDataMet->Chi2Test(hPdfMet,"PUW");
  chi2ndf  = hDataMet->Chi2Test(hPdfMet,"CHI2/NDFUW");
  ksprob   = hDataMet->KolmogorovTest(hPdfMet);
  ksprobpe = hDataMet->KolmogorovTest(hPdfMet,"DX");
  sprintf(txtfname,"%s/fitresWe.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << hDataMet->Integral() << endl;
  txtfile << "  Signal: " << nSig.getVal() << " +/- " << nSig.getPropagatedError(*fitRes) << endl;
  txtfile << "     QCD: " << nQCD.getVal() << " +/- " << nQCD.getPropagatedError(*fitRes) << endl;
  txtfile << "   Other: " << nEWK.getVal() << " +/- " << nEWK.getPropagatedError(*fitRes) << endl;
  txtfile << endl; 
  txtfile.flags(flags);
  
  fitRes->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitRes);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();
  
  chi2prob = hDataMetp->Chi2Test(hPdfMetp,"PUW");
  chi2ndf  = hDataMetp->Chi2Test(hPdfMetp,"CHI2/NDFUW");
  ksprob   = hDataMetp->KolmogorovTest(hPdfMetp);
  ksprobpe = hDataMetp->KolmogorovTest(hPdfMetp,"DX");  
  sprintf(txtfname,"%s/fitresWep.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << hDataMetp->Integral() << endl;
  txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResp) << endl;
  txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResp) << endl;
  txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResp) << endl;
  txtfile << endl;  
  txtfile.flags(flags);
  
  fitResp->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResp);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();

  chi2prob = hDataMetm->Chi2Test(hPdfMetm,"PUW");
  chi2ndf  = hDataMetm->Chi2Test(hPdfMetm,"CHI2/NDFUW");
  ksprob   = hDataMetm->KolmogorovTest(hPdfMetm);
  ksprobpe = hDataMetm->KolmogorovTest(hPdfMetm,"DX");  
  sprintf(txtfname,"%s/fitresWem.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << hDataMetm->Integral() << endl;
  txtfile << "  Signal: " << nSigm.getVal() << " +/- " << nSigm.getPropagatedError(*fitResm) << endl;
  txtfile << "     QCD: " << nQCDm.getVal() << " +/- " << nQCDm.getPropagatedError(*fitResm) << endl;
  txtfile << "   Other: " << nEWKm.getVal() << " +/- " << nEWKm.getPropagatedError(*fitResm) << endl;
  txtfile << endl;
  txtfile.flags(flags);
  
  fitResm->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResm);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();

  makeHTML(outputDir);
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  */
  gBenchmark->Show("fitWe");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    
    Double_t err = sqrt(hData->GetBinContent(ibin));
    if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    
    if(err>0) hDiff->SetBinContent(ibin,diff/err);
    else      hDiff->SetBinContent(ibin,0);
    hDiff->SetBinError(ibin,1);   
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
void printCorrelations(ostream& os, RooFitResult *res)
{
  ios_base::fmtflags flags = os.flags();
  const RooArgList parlist = res->floatParsFinal();
  
  os << "  Correlation Matrix" << endl;
  os << " --------------------" << endl;
  for(Int_t i=0; i<parlist.getSize(); i++) {
    for(Int_t j=0; j<parlist.getSize(); j++) 
      os << "  " << setw(7) << setprecision(4) << fixed << res->correlationMatrix()(i,j);    
    os << endl;
  }
  os.flags(flags);
}

//--------------------------------------------------------------------------------------------------
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
			   const Double_t ksprob, const Double_t ksprobpe)
{
  ios_base::fmtflags flags = os.flags();
  
  os << "  Chi2 Test" << endl;
  os << " -----------" << endl;
  os << "       prob = " << chi2prob << endl;
  os << "   chi2/ndf = " << chi2ndf << endl;
  os << endl;
  os << "  KS Test" << endl;
  os << " ---------" << endl;
  os << "   prob = " << ksprob << endl;
  os << "   prob = " << ksprobpe << " with 1000 pseudo-experiments" << endl;
  os << endl;
 
  os.flags(flags);
}
			   
//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/WenuFitPlots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Wenu</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmet.png\"><img src=\"fitmet.png\" alt=\"fitmet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetp.png\"><img src=\"fitmetp.png\" alt=\"fitmetp.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetm.png\"><img src=\"fitmetm.png\" alt=\"fitmetm.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetlog.png\"><img src=\"fitmetlog.png\" alt=\"fitmetlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetplog.png\"><img src=\"fitmetplog.png\" alt=\"fitmetplog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetmlog.png\"><img src=\"fitmetmlog.png\" alt=\"fitmetmlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();  
}

void makeDatacard(const TString outputDir, const TString channel, const TString ewk_option, double ewkSF, double nQCD_exp, double lumi)
{
  //
  // Write datacard
  //
  ofstream txtfile;
  char txtfname[100];    
  
  ios_base::fmtflags flags;
  
  // Expected signal yields calculated from theoretical predictions of cross sections
  Double_t nWenu_exp  = lumi*1e6*12.21*0.479*0.695;
  Double_t nWenup_exp = lumi*1e6*7.11*0.484*0.687;
  Double_t nWenum_exp = lumi*1e6*5.09*0.471*0.708;
  Double_t nWmunu_exp  = lumi*1e6*12.21*0.440*0.839;
  Double_t nWmunup_exp = lumi*1e6*7.11*0.441*0.843;
  Double_t nWmunum_exp = lumi*1e6*5.09*0.439*0.831;
  // Expected signal yields for anti-isolation samples are the initial values of nAntiSig, nAntiSigp, and nAntiSigm
  Double_t nAntiWmunu_exp  = 1231.100;
  Double_t nAntiWmunup_exp = 627.250;
  Double_t nAntiWmunum_exp = 603.800;
  Double_t nSig_exp = 0.0;

  // Names of PDF's
  TString PdfFileName, dataPdfname, sigPdfname, ewkPdfname, qcdPdfname;

  // Systematic uncertainties
  Double_t lepton_eff, pT_scale_res, met_scale_res, bkg_model, lumi_8TeV;
  Double_t ewk_constraint;
  if (ewk_option=="tied") ewk_constraint = 1.15;
  else if (ewk_option=="floating") ewk_constraint = 0.0;
  else cout << "Invalid value for ewk_option" << endl;
  if (channel=="Wenu") {
    nSig_exp = nWenu_exp;
    PdfFileName = "Wenu"; dataPdfname = "dataMet"; sigPdfname = "we"; ewkPdfname = "ewk"; qcdPdfname = "pepe1Pdf_qcd";
    lepton_eff = 1.025; pT_scale_res = 1.005; met_scale_res = 1.008; bkg_model = 1.003; lumi_8TeV = 1.026; }
  else if (channel=="Wenu_p") {
    nSig_exp = nWenup_exp;
    PdfFileName = "Wenu"; dataPdfname = "dataMetp"; sigPdfname = "wep"; ewkPdfname = "ewkp"; qcdPdfname = "pepe1Pdf_qcdp";
    lepton_eff = 1.028; pT_scale_res = 1.004; met_scale_res = 1.008; bkg_model = 1.002; lumi_8TeV = 1.026; }
  else if (channel=="Wenu_m") {
    nSig_exp = nWenum_exp;
    PdfFileName = "Wenu"; dataPdfname = "dataMetm"; sigPdfname = "wem"; ewkPdfname = "ewkm"; qcdPdfname = "pepe1Pdf_qcdm";
    lepton_eff = 1.025; pT_scale_res = 1.007; met_scale_res = 1.007; bkg_model = 1.003; lumi_8TeV = 1.026; }
  else if (channel=="Wmunu_sel" || channel=="Wmunu_anti") {
    if (channel=="Wmunu_sel") {
      nSig_exp = nWmunu_exp;
      PdfFileName = "Wmunu"; dataPdfname = "dataMet"; sigPdfname = "wm"; ewkPdfname = "ewk"; qcdPdfname = "pepe1Pdf_qcd"; }
    else {
      nSig_exp = nAntiWmunu_exp;
      PdfFileName = "Wmunu"; dataPdfname = "antiMet"; sigPdfname = "awm"; ewkPdfname = "aewk"; qcdPdfname = "pepe1Pdf_aqcd"; }
    lepton_eff = 1.010; pT_scale_res = 1.003; met_scale_res = 1.005; bkg_model = 1.001; lumi_8TeV = 1.026; }
  else if (channel=="Wmunup_sel" || channel=="Wmunup_anti") {
    if (channel=="Wmunup_sel") {
      nSig_exp = nWmunup_exp;
      PdfFileName = "Wmunu"; dataPdfname = "dataMetp"; sigPdfname = "wmp"; ewkPdfname = "ewkp"; qcdPdfname = "pepe1Pdf_qcdp"; }
    else {
      nSig_exp = nAntiWmunup_exp;
      PdfFileName = "Wmunu"; dataPdfname = "antiMetp"; sigPdfname = "awmp"; ewkPdfname = "aewkp"; qcdPdfname = "pepe1Pdf_aqcdp"; }
    lepton_eff = 1.010; pT_scale_res = 1.003; met_scale_res = 1.005; bkg_model = 1.002; lumi_8TeV = 1.026; }
  else if (channel=="Wmunum_sel" || channel=="Wmunum_anti") {
    if (channel=="Wmunum_sel") {
      nSig_exp = nWmunum_exp;
      PdfFileName = "Wmunu"; dataPdfname = "dataMetm"; sigPdfname = "wmm"; ewkPdfname = "ewkm"; qcdPdfname = "pepe1Pdf_qcdm"; }
    else {
      nSig_exp = nAntiWmunum_exp;
      PdfFileName = "Wmunu"; dataPdfname = "antiMetm"; sigPdfname = "awmm"; ewkPdfname = "aewkm"; qcdPdfname = "pepe1Pdf_aqcdm"; }
    lepton_eff = 1.009; pT_scale_res = 1.003; met_scale_res = 1.005; bkg_model = 1.001; lumi_8TeV = 1.026; }
  else {
    cout << "Not a valid channel name" << endl; }

  // EWK yield is either tied to the Wenu yield
  // or left to float freely in the fit
  Double_t nEWK_exp = nSig_exp*ewkSF;

  sprintf(txtfname,"/home/cmedlock/cms/cmssw/CMSSW_6_1_2/src/HiggsAnalysis/CombinedLimit/%s/"+channel+".txt",outputDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(5);
  txtfile << "imax 1 number of bins" << endl;
  txtfile << "jmax 2 number of processes" << endl;
  txtfile << "kmax * number of nuisance parameters" << endl;
  txtfile << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
  txtfile << "shapes data_obs" << setw(12) << channel << setw(58) << PdfFileName+"_pdfTemplates.root  combine_workspace:" << dataPdfname << endl;
  txtfile << "shapes Sig" << setw(17) << channel << setw(58) << PdfFileName+"_pdfTemplates.root  combine_workspace:" << sigPdfname << endl;
  txtfile << "shapes EWK" << setw(17) << channel << setw(58) << PdfFileName+"_pdfTemplates.root  combine_workspace:" << ewkPdfname << endl;
  txtfile << "shapes QCD" << setw(17) << channel << setw(58) << PdfFileName+"_pdfTemplates.root  combine_workspace:" << qcdPdfname << endl;
  txtfile << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
  txtfile << "bin" << setw(24) << channel << endl;
  txtfile << "observation            -1.0" << endl;
  txtfile << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
  txtfile << "bin" << setw(34) << channel << setw(12) << channel << setw(12) << channel << endl;
  txtfile << "process                             0           1           2" << endl;
  txtfile << "process                           Sig         EWK         QCD" << endl;
  txtfile << "rate" << setw(33) << nSig_exp << setw(12) << nEWK_exp << setw(12) << nQCD_exp << endl;
  txtfile << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
  txtfile << setprecision(4);
  txtfile << "lepton_eff              lnN" << setw(10) << lepton_eff << setw(12) << lepton_eff << setw(12) << "-" << endl;
  txtfile << "pT_scale_res            lnN" << setw(10) << pT_scale_res << setw(12) << pT_scale_res << setw(12) << "-" << endl;
  txtfile << "met_scale_res           lnN" << setw(10) << met_scale_res << setw(12) << met_scale_res << setw(12) << "-" << endl;
  txtfile << "bkg_model               lnN" << setw(10) << bkg_model << setw(12) << bkg_model << setw(12) << "-" << endl;
  txtfile << "lumi_8TeV               lnN" << setw(10) << lumi_8TeV << setw(12) << lumi_8TeV << setw(12) << "-" << endl;
  txtfile << "ewk_constraint          lnN" << setw(10) << "-" << setw(12) << ewk_constraint << setw(12) << "-" << endl;
  txtfile << "qcd_constraint          lnN" << setw(10) << "-" << setw(12) << "-" << setw(12) << "2.0" << endl;
  txtfile << endl;
  txtfile.flags(flags);
  txtfile.close();  
}
