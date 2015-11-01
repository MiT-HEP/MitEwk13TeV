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
#include <TMath.h> // ROOT math library
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include <TGaxis.h>
#include "TLorentzVector.h"           // 4-vector class

#include "../Utils/MyTools.hh"	          // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../Utils/RecoilCorrector_v2.hh"    // class to handle recoil corrections for MET

#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections
// #include "ZBackgrounds.hh"
#include "RooCategory.h"

// RooFit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

// print correlations of fitted parameters
void printCorrelations(ostream& os, RooFitResult *res);

// print chi2 test and KS test results
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
			   const Double_t ksprob, const Double_t ksprobpe);

void apply_postfit_shape(TH1D *histogram, TH1D *up1, TH1D *down1, TH1D *up2, TH1D *down2, double shift, double shift_err);

// make webpage
void makeHTML(const TString outDir);


//=== MAIN MACRO ================================================================================================= 

void postFitWe(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)
	   const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
) {
  gBenchmark->Start("postFitWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  // MET histogram binning and range
  const Int_t    NBINS  = 75;
  const Double_t METMAX = 150;

  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.5;

  TString pufname = "../Tools/pileup_weights_2015B.root";

  // file format for output plots
  const TString format("png"); 

  // recoil correction
  // RecoilCorrector recoilCorr("../Recoil/ZeeData/fits_puppi.root");//, (!) uncomment to perform corrections to recoil from W-MC/Z-MC
                             //"../Recoil/WepMC/fits.root",
			     //"../Recoil/WemMC/fits.root",

  TFile *pufile = new TFile(pufname); assert(pufile);
  TH1D  *puWeights = (TH1D*)pufile->Get("npv_rw");

  
  // recoil correctors

  RecoilCorrector *recoilCorr = new  RecoilCorrector("../Recoil/WmpMC/fits.root","fcnPF");// PUPPI!
    recoilCorr->addMCFile("../Recoil/ZmmMC_default/fits_puppi_ff.root"); // PUPPI
  recoilCorr->addDataFile("../Recoil/ZmumuData_puppi_lin_10_12_small/fits_puppi.root");// PUPPI
  recoilCorr->addMCTrueFile("../Recoil/ZmmData/fits_puppi_new.root");// PUPPI
//     RecoilCorrector *recoilCorr = new  RecoilCorrector("../Recoil/WmunuPlus_MC_mvaFixed_2015_09_22/fits.root","fcnPF");// MVA!
//   recoilCorr->addMCFile("../Recoil/Zmumu_MC_mvaFixed_2015_09_22/fits_mva.root"); // MVA!
//   recoilCorr->addDataFile("../Recoil/Zmumu_Data_mvaFixed_2015_09_22/fits_mva.root"); // MVA!
//   recoilCorr->addMCTrueFile("../Recoil/Zmumu_Data_mvaFixed_2015_09_22/fits_mva.root"); // MVA!

  
  RecoilCorrector *recoilCorrm = new  RecoilCorrector("../Recoil/WmmMC/fits.root","fcnPF"); // PUPPI!
  recoilCorrm->addMCFile("../Recoil/ZmmMC_default/fits_puppi_ff.root");//puppi
  recoilCorrm->addDataFile("../Recoil/ZmumuData_puppi_lin_10_12_small/fits_puppi.root"); //puppi
  recoilCorrm->addMCTrueFile("../Recoil/ZmmData/fits_puppi_new.root"); //puppi
//    RecoilCorrector *recoilCorrm = new  RecoilCorrector("../Recoil/WmunuMinus_MC_mvaFixed_2015_09_22/fits.root","fcnPF"); // MVA!
//   recoilCorrm->addMCFile("../Recoil/Zmumu_MC_mvaFixed_2015_09_22/fits_mva.root"); // MVA
//   recoilCorrm->addDataFile("../Recoil/Zmumu_Data_mvaFixed_2015_09_22/fits_mva.root"); // MVA
//   recoilCorrm->addMCTrueFile("../Recoil/Zmumu_Data_mvaFixed_2015_09_22/fits_mva.root"); // MVA

  enum { eData, eWenu, eEWK };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  
  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_testMVA2/Wenu/ntuples/data_select.root"); typev.push_back(eData);
  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_testMVA2/Wenu/ntuples/we_select.raw.root");   typev.push_back(eWenu);
  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_testMVA2/Wenu/ntuples/ewk_select.raw.root");  typev.push_back(eEWK);
  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_testMVA2/Wenu/ntuples/top_select.raw.root");  typev.push_back(eEWK);
  
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
  TH1D *hWenuMet_RecoilUp  = new TH1D("hWenuMet_RecoilUp", "",NBINS,0,METMAX); hWenuMet_RecoilUp->Sumw2();
  TH1D *hWenuMetp_RecoilUp = new TH1D("hWenuMetp_RecoilUp","",NBINS,0,METMAX); hWenuMetp_RecoilUp->Sumw2();
  TH1D *hWenuMetm_RecoilUp = new TH1D("hWenuMetm_RecoilUp","",NBINS,0,METMAX); hWenuMetm_RecoilUp->Sumw2();
  TH1D *hWenuMet_RecoilDown  = new TH1D("hWenuMet_RecoilDown", "",NBINS,0,METMAX); hWenuMet_RecoilDown->Sumw2();
  TH1D *hWenuMetp_RecoilDown = new TH1D("hWenuMetp_RecoilDown","",NBINS,0,METMAX); hWenuMetp_RecoilDown->Sumw2();
  TH1D *hWenuMetm_RecoilDown = new TH1D("hWenuMetm_RecoilDown","",NBINS,0,METMAX); hWenuMetm_RecoilDown->Sumw2();
  TH1D *hWenuMet_ScaleUp  = new TH1D("hWenuMet_ScaleUp", "",NBINS,0,METMAX); hWenuMet_ScaleUp->Sumw2();
  TH1D *hWenuMetp_ScaleUp = new TH1D("hWenuMetp_ScaleUp","",NBINS,0,METMAX); hWenuMetp_ScaleUp->Sumw2();
  TH1D *hWenuMetm_ScaleUp = new TH1D("hWenuMetm_ScaleUp","",NBINS,0,METMAX); hWenuMetm_ScaleUp->Sumw2();
  TH1D *hWenuMet_ScaleDown  = new TH1D("hWenuMet_ScaleDown", "",NBINS,0,METMAX); hWenuMet_ScaleDown->Sumw2();
  TH1D *hWenuMetp_ScaleDown = new TH1D("hWenuMetp_ScaleDown","",NBINS,0,METMAX); hWenuMetp_ScaleDown->Sumw2();
  TH1D *hWenuMetm_ScaleDown = new TH1D("hWenuMetm_ScaleDown","",NBINS,0,METMAX); hWenuMetm_ScaleDown->Sumw2();


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
  TLorentzVector *lep=0;
  TLorentzVector *sc=0;
    
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
    intree->SetBranchAddress("puppiMet",      &met);       // MET
    intree->SetBranchAddress("puppiMetPhi",   &metPhi);    // phi(MET)
//     intree->SetBranchAddress("mvaMet",      &met);       // MET
//     intree->SetBranchAddress("mvaMetPhi",   &metPhi);    // phi(MET)
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

      double pU1         = 0;  //--
      double pU2         = 0;  //--
    
      //if(sc->Pt()        < PT_CUT)  continue;	
      if(fabs(lep->Eta()) > ETA_CUT) continue;
  
      mt     = sqrt( 2.0 * (lep->Pt()) * (met) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metPhi))) );

      if(typev[ifile]==eData) {
	if(lep->Pt()        < PT_CUT)  continue;
        hDataMet->Fill(met);
	if(q>0) { hDataMetp->Fill(met); } 
	else    { hDataMetm->Fill(met); }
      
      } else {
        Double_t weight = 1;
        weight *= scale1fb*lumi;
	weight *=puWeights->GetBinContent(npv+1);
	if(typev[ifile]==eWenu) {
          Double_t corrMet=met, corrMetPhi=metPhi;
	  
	  // apply recoil corrections to W MC
	  //Double_t lepPt = lep->Pt();
	  Double_t lepPt = (gRandom->Gaus((lep->Pt())*getEleScaleCorr(lep->Eta(),0),getEleResCorr(lep->Eta(),0)));  // (!) uncomment to apply scale/res corrections to MC
	  Double_t lepPtup = (gRandom->Gaus((lep->Pt())*getEleScaleCorr(lep->Eta(),1),getEleResCorr(lep->Eta(),1)));  // (!) uncomment to apply scale/res corrections to MC
	  Double_t lepPtdown = (gRandom->Gaus((lep->Pt())*getEleScaleCorr(lep->Eta(),-1),getEleResCorr(lep->Eta(),-1)));  // (!) uncomment to apply scale/res corrections to MC
	    if(lepPt        > PT_CUT) 
	    {  
	      //recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),nsigma,q);
//           recoilCorr->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
	      hWenuMet->Fill(corrMet,weight);
	      if(q>0) 
		{
		  recoilCorr->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
		  hWenuMetp->Fill(corrMet,weight); 
		} 
	      else    
		{ 
		  recoilCorrm->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
		  hWenuMetm->Fill(corrMet,weight); 
		}
	      corrMet=met, corrMetPhi=metPhi;
	      //recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),1,q);
//           recoilCorr->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,2,2);
	      hWenuMet_RecoilUp->Fill(corrMet,weight);
	      if(q>0) 
		{
		  recoilCorr->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,2,2);
		  hWenuMetp_RecoilUp->Fill(corrMet,weight); 
		} 
	      else    
		{ 
		  recoilCorrm->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,2,2);
		  hWenuMetm_RecoilUp->Fill(corrMet,weight); 
		}
	      corrMet=met, corrMetPhi=metPhi;
	      //recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),-1,q);
//           recoilCorr->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,-2,-2);
	      hWenuMet_RecoilDown->Fill(corrMet,weight);
	      if(q>0) 
		{ 
		  recoilCorr->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,-2,-2);
		  hWenuMetp_RecoilDown->Fill(corrMet,weight); 
		} 
	      else    
		{ 
		  recoilCorrm->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,-2,-2);
		  hWenuMetm_RecoilDown->Fill(corrMet,weight); 
		}
	    }
	    if(lepPtup        > PT_CUT)
	      {
		corrMet=met, corrMetPhi=metPhi;
		//recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPtup,lep->Phi(),0,q);
//         recoilCorr->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
		hWenuMet_ScaleUp->Fill(corrMet,weight);
		if(q>0) 
		  {
		    recoilCorr->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPtup,lep->Phi(),pU1,pU2,0);
		    hWenuMetp_ScaleUp->Fill(corrMet,weight); 
		  } 
		else   
		  { 
		    recoilCorrm->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPtup,lep->Phi(),pU1,pU2,0);
		    hWenuMetm_ScaleUp->Fill(corrMet,weight); 
		  }
	      }
	    if(lepPtdown        > PT_CUT)
	    {
	      corrMet=met, corrMetPhi=metPhi;
	      //recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPtdown,lep->Phi(),0,q);
//           recoilCorr->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
	      hWenuMet_ScaleDown->Fill(corrMet,weight);
	      if(q>0) 
		{ 
		  recoilCorr->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPtdown,lep->Phi(),pU1,pU2,0);
		  hWenuMetp_ScaleDown->Fill(corrMet,weight); 
		} 
	      else   
		{ 
		  recoilCorrm->CorrectType2(corrMet,corrMetPhi,genVPt,genVPhi,lepPtdown,lep->Phi(),pU1,pU2,0);
		  hWenuMetm_ScaleDown->Fill(corrMet,weight); 
		}
	    }
        }
	if(typev[ifile]==eEWK) {
	  if(lep->Pt()        < PT_CUT)  continue;
          hEWKMet->Fill(met,weight);
	  if(q>0) { hEWKMetp->Fill(met,weight); }
	  else    { hEWKMetm->Fill(met,weight); }
	}
      }
    }
  }  
  delete infile;
  infile=0, intree=0;   
  
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  //
  RooRealVar nSig("nSig","nSig",0.7*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar nQCD("nQCD","nQCD",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar cewk("cewk","cewk",0.1,0,5);
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
//   apply_postfit_shape(hWenuMetp, hWenuMetp_RecoilUp, hWenuMetp_RecoilDown,hWenuMetp_ScaleUp, hWenuMetp_ScaleDown, 0.0471 , -0.8684); // mva
//   apply_postfit_shape(hWenuMetm, hWenuMetm_RecoilUp, hWenuMetm_RecoilDown,hWenuMetm_ScaleUp, hWenuMetm_ScaleDown,  -0.1430,-0.4277); //mva 
    apply_postfit_shape(hWenuMetp, hWenuMetp_RecoilUp, hWenuMetp_RecoilDown,hWenuMetp_ScaleUp, hWenuMetp_ScaleDown, -0.5326, -0.9921); // puppi
  apply_postfit_shape(hWenuMetm, hWenuMetm_RecoilUp, hWenuMetm_RecoilDown,hWenuMetm_ScaleUp, hWenuMetm_ScaleDown,  -0.5506, -0.3567); //puppi 
  RooDataHist wenuMet ("wenuMET", "wenuMET", RooArgSet(pfmet),hWenuMet);  RooHistPdf pdfWe ("we", "we", pfmet,wenuMet, 1);
  RooDataHist wenuMetp("wenuMETp","wenuMETp",RooArgSet(pfmet),hWenuMetp); RooHistPdf pdfWep("wep","wep",pfmet,wenuMetp,1);
  RooDataHist wenuMetm("wenuMETm","wenuMETm",RooArgSet(pfmet),hWenuMetm); RooHistPdf pdfWem("wem","wem",pfmet,wenuMetm,1); 
  
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", pfmet,ewkMet, 1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",pfmet,ewkMetp,1); 
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",pfmet,ewkMetm,1); 
  
  // QCD Pdfs
  CPepeModel1 qcd("qcd",pfmet);
  CPepeModel1 qcdp("qcdp",pfmet);
  CPepeModel1 qcdm("qcdm",pfmet);  

  // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWe,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWep,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWem,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm)); 
  
  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet),hDataMet);
  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(pfmet),hDataMetp);
  RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(pfmet),hDataMetm);


  //Now lets get those lovely fit parameters
  std::cout  << "Selectedp: " << hDataMetp->Integral() << endl;
  std::cout  << "Selectedm: " << hDataMetm->Integral() << endl;
  // MVA MET #'s
//   qcdp.a1->setVal(0.2177);
//   qcdp.sigma->setVal(15.9620);
//   nQCDp.setVal(85515.3812);
//   nSigp.setVal(0.9815*119270.0174); //118876.426343
//   
//   qcdm.a1->setVal(0.2216);
//   qcdm.sigma->setVal(15.2119);
//   nQCDm.setVal(78369.3676);
//   nSigm.setVal(0.9614*98640.6); // 100718.15014
  
  // nfinal
  qcdp.a1->setVal(0.2268);
  qcdp.sigma->setVal(14.1877);
  nQCDp.setVal(79665.5222);
  nSigp.setVal(0.9893*123643.3131); 
  
  qcdm.a1->setVal(0.2190);
  qcdm.sigma->setVal(14.2202);
  nQCDm.setVal(74456.3999);
  nSigm.setVal(1.0036*97913.78449); 

  
   
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
  c->cd(1)->SetBottomMargin(0.02);
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
  if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  (8 TeV)",lumi*1000.);
  else         sprintf(lumitext,"%.0f pb^{-1}  (13 TeV)",lumi);
  
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
  CPlot plotMet("fitmet",weframe,"","mT [GeV]",ylabel);
  plotMet.SetLegend(0.68,0.57,0.93,0.77);
  plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMet.GetLegend()->AddEntry(hDummyW,"W#rightarrowe#nu","F");
  plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  //plotMet.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  //plotMet.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  plotMet.AddTextBox("CMS Preliminary",0.55,0.80,0.90,0.86,0);
  plotMet.AddTextBox(lumitext,0.63,0.92,0.95,0.99,0);
//  plotMet.SetYRange(0.1,1.1*(hDataMet->GetMaximum()));
//plotMet.SetYRange(0.1,1e4);
  plotMet.Draw(c,kTRUE,format,1);

  CPlot plotMetDiff("fitmet","","#slash{E}_{T} [GeV]","#chi");
  plotMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
  plotMetDiff.SetYRange(-0.1,0.1);
  plotMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetDiff.AddLine(0, 0.05,METMAX, 0.05,kBlack,3);
  plotMetDiff.AddLine(0,-0.05,METMAX,-0.05,kBlack,3);
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
  wepframe->GetXaxis()->SetLabelOffset(2.0);
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
  plotMetp.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.62,0.80,0.88,0.88,0);
  plotMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
//  plotMetp.SetYRange(0.1,1.1*(hDataMetp->GetMaximum()));
//plotMetp.SetYRange(0.1,5000);
  plotMetp.Draw(c,kFALSE,format,1);

  CPlot plotMetpDiff("fitmetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetpDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
  plotMetpDiff.SetYRange(-0.20,0.20);
  plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotMetpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotMetpDiff.Draw(c,kTRUE,format,2);
  plotMetpDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetp.SetName("fitmetplog");
  plotMetp.SetLogy();
  plotMetp.SetYRange(1e-3*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
  plotMetp.Draw(c,kTRUE,format,1);
  plotMetp.Draw(c,kTRUE,"pdf",1);
  
  //
  // W- MET plot
  //
  RooPlot *wemframe = pfmet.frame(Bins(NBINS)); 
  wemframe->GetYaxis()->SetNdivisions(505);
  wemframe->GetXaxis()->SetLabelOffset(2.0);
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
  plotMetm.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.62,0.80,0.88,0.88,0);
  plotMetm.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
//  plotMetm.SetYRange(0.1,1.1*(hDataMetm->GetMaximum()));
//plotMetm.SetYRange(0.1,5000);
  plotMetm.Draw(c,kFALSE,format,1);

  CPlot plotMetmDiff("fitmetm","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetmDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetmDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetmDiff.AddHist1D(hMetmDiff,"EX0",ratioColor);
  plotMetmDiff.SetYRange(-0.2,0.2);
  plotMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetmDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotMetmDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotMetmDiff.Draw(c,kTRUE,format,2);
  plotMetmDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetm.SetName("fitmetmlog");
  plotMetm.SetLogy();
  plotMetm.SetYRange(1e-3*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
  plotMetm.Draw(c,kTRUE,format,1);
  plotMetm.Draw(c,kTRUE,"pdf",1);

  double chi2prob = hDataMetp->Chi2Test(hPdfMetp,"PUW");
  double chi2ndf  = hDataMetp->Chi2Test(hPdfMetp,"CHI2/NDFUW");
  std::cout << chi2prob << " " << chi2ndf << std::endl;
  chi2prob = hDataMetm->Chi2Test(hPdfMetm,"PUW");
  chi2ndf  = hDataMetm->Chi2Test(hPdfMetm,"CHI2/NDFUW");
  std::cout << chi2prob << " " << chi2ndf << std::endl;
  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  
  makeHTML(outputDir);
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  gBenchmark->Show("postFitWe");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    Double_t diff = diff0/hData->GetBinContent(ibin);
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((1.0/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    //Double_t err = sqrt(hData->GetBinContent(ibin));
    //if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    //if(err>0) hDiff->SetBinContent(ibin,diff/err);
    //else      hDiff->SetBinContent(ibin,0);
    hDiff->SetBinContent(ibin,diff);
    hDiff->SetBinError(ibin,err);   
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

void apply_postfit_shape(TH1D *histogram, TH1D *up1, TH1D *down1,TH1D *up2, TH1D *down2, double shift1, double shift2)
{
  for(int bin=1; bin<=histogram->GetNbinsX(); bin++)
    {
      double upper1 = up1->GetBinContent(bin);
      double lower1 = down1->GetBinContent(bin);
      double upper2 = up2->GetBinContent(bin);
      double lower2 = down2->GetBinContent(bin);
      double central = histogram->GetBinContent(bin);
      double value1 = 1, value2=1;
      if (down1->GetBinContent(bin) < 0.01)
	lower1 = 0;
      if (down2->GetBinContent(bin) < 0.01)
	lower2 = 0;
      if(shift1 > 0 && central > 0 && upper1 != central)
	value1 = shift1 * (upper1 / central - 1) + 1;
      else if(shift1 < 0 && central > 0 && lower1 != central)
	value1 = abs(shift1) * (lower1 / central - 1) + 1;
      if(shift2 > 0 && central > 0 && upper2 != central)
	value2 = shift2 * (upper2 / central - 1) + 1;
      else if(shift2 < 0 && central > 0 && lower2 != central)
	value2 = abs(shift2) * (lower2 / central - 1) + 1;
      if(value1 != 1 && value2 != 1)
	{
	  //central *= value1;
	  central = central*value1+central*value2-central;
	  histogram->SetBinContent(bin, central);
	}
      /*  double uncertainty = 1;
      if (lower && central)
	{
	  uncertainty = TMath::Max(uncertainty,central/lower);
	  uncertainty = TMath::Max(uncertainty,lower/central);
	}
      if (upper && central)
	{
	  uncertainty = TMath::Max(uncertainty,upper / central);
	  uncertainty = TMath::Max(uncertainty,central / upper);
	}
      double temp = uncertainty-1;
      uncertainty = shift_err * TMath::Min(2.0,temp);
      double error = sqrt(histogram->GetBinError(bin)*histogram->GetBinError(bin)+ (histogram->GetBinContent(bin) * uncertainty)*(histogram->GetBinContent(bin) * uncertainty));
      histogram->SetBinError(bin, error);*/
    }
}
