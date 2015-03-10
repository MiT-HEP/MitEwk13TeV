//================================================================================================
//
// Make plots of various distributions after Wenu selection 
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include "Math/LorentzVector.h"           // 4-vector class

#include "ConfParse.hh"                   // input conf file parser
#include "../Utils/CSample.hh"            // helper class to handle samples
#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hMC, const TString name);

// make webpage
void makeHTML(const TString outDir);


//=== MAIN MACRO ================================================================================================= 

void plotWe(const TString  conf,            // input file
            const TString  inputDir,        // input directory
	    const TString  outputDir,       // output directory
	    const Double_t lumi             // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  const TString format("png");

  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.5;

  const Double_t ETA_BARREL = 1.4442;
  const Double_t ETA_ENDCAP = 1.566;
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  //
  // parse .conf file
  //
  confParse(conf, snamev, samplev);
  const Bool_t hasData = (samplev[0]->fnamev.size()>0);

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir + TString("/plots");

  
  //
  // Create histograms
  //
  enum { eInc, ePos, eNeg };
  vector<TH1D*> hPtv[3], hPtBv[3], hPtEv[3];
  vector<TH1D*> hPt2v[3], hPtB2v[3], hPtE2v[3];
  vector<TH1D*> hEtav[3];
  vector<TH1D*> hPhiv[3], hPhiBv[3], hPhiEv[3];
  vector<TH1D*> hMetv[3], hMetBv[3], hMetEv[3];
  vector<TH1D*> hMet2v[3], hMetB2v[3], hMetE2v[3];
  vector<TH1D*> hMetPhiv[3], hMetPhiBv[3], hMetPhiEv[3];
  vector<TH1D*> hMtv[3], hMtBv[3], hMtEv[3];
  vector<TH1D*> hMt2v[3], hMtB2v[3], hMtE2v[3];
  vector<TH1D*> hNPVv[3];
   
  TH1D *hPtMC[3], *hPtBMC[3], *hPtEMC[3];
  TH1D *hPt2MC[3], *hPtB2MC[3], *hPtE2MC[3];
  TH1D *hEtaMC[3];
  TH1D *hPhiMC[3], *hPhiBMC[3], *hPhiEMC[3];
  TH1D *hMetMC[3], *hMetBMC[3], *hMetEMC[3];
  TH1D *hMet2MC[3], *hMetB2MC[3], *hMetE2MC[3];
  TH1D *hMetPhiMC[3], *hMetPhiBMC[3], *hMetPhiEMC[3];
  TH1D *hMtMC[3], *hMtBMC[3], *hMtEMC[3];
  TH1D *hMt2MC[3], *hMtB2MC[3], *hMtE2MC[3];
  TH1D *hNPVMC[3];
  
  char hname[100];
  for(UInt_t ich=0; ich<3; ich++) {
    for(UInt_t isam=0; isam<samplev.size(); isam++) {
      sprintf(hname,"hPt_ch%i_%i",ich,isam);       hPtv[ich].push_back(new TH1D(hname,"",30,20,80));        hPtv[ich][isam]->Sumw2();
      sprintf(hname,"hPtB_ch%i_%i",ich,isam);      hPtBv[ich].push_back(new TH1D(hname,"",30,20,80));       hPtBv[ich][isam]->Sumw2();
      sprintf(hname,"hPtE_ch%i_%i",ich,isam);      hPtEv[ich].push_back(new TH1D(hname,"",30,20,80));       hPtEv[ich][isam]->Sumw2();
      sprintf(hname,"hPt2_ch%i_%i",ich,isam);      hPt2v[ich].push_back(new TH1D(hname,"",40,20,180));       hPt2v[ich][isam]->Sumw2();
      sprintf(hname,"hPtB2_ch%i_%i",ich,isam);     hPtB2v[ich].push_back(new TH1D(hname,"",40,20,180));      hPtB2v[ich][isam]->Sumw2();
      sprintf(hname,"hPtE2_ch%i_%i",ich,isam);     hPtE2v[ich].push_back(new TH1D(hname,"",40,20,180));      hPtE2v[ich][isam]->Sumw2();
      sprintf(hname,"hEta_ch%i_%i",ich,isam);      hEtav[ich].push_back(new TH1D(hname,"",40,-3,3));         hEtav[ich][isam]->Sumw2();
      sprintf(hname,"hPhi_ch%i_%i",ich,isam);      hPhiv[ich].push_back(new TH1D(hname,"",40,-3.2,3.2));     hPhiv[ich][isam]->Sumw2();
      sprintf(hname,"hPhiB_ch%i_%i",ich,isam);     hPhiBv[ich].push_back(new TH1D(hname,"",40,-3.2,3.2));    hPhiBv[ich][isam]->Sumw2();
      sprintf(hname,"hPhiE_ch%i_%i",ich,isam);     hPhiEv[ich].push_back(new TH1D(hname,"",40,-3.2,3.2));    hPhiEv[ich][isam]->Sumw2();
      sprintf(hname,"hMet_ch%i_%i",ich,isam);      hMetv[ich].push_back(new TH1D(hname,"",40,0,80));        hMetv[ich][isam]->Sumw2();
      sprintf(hname,"hMetB_ch%i_%i",ich,isam);     hMetBv[ich].push_back(new TH1D(hname,"",40,0,80));       hMetBv[ich][isam]->Sumw2();
      sprintf(hname,"hMetE_ch%i_%i",ich,isam);     hMetEv[ich].push_back(new TH1D(hname,"",40,0,80));       hMetEv[ich][isam]->Sumw2();
      sprintf(hname,"hMet2_ch%i_%i",ich,isam);     hMet2v[ich].push_back(new TH1D(hname,"",40,0,200));       hMet2v[ich][isam]->Sumw2();
      sprintf(hname,"hMetB2_ch%i_%i",ich,isam);    hMetB2v[ich].push_back(new TH1D(hname,"",40,0,200));      hMetB2v[ich][isam]->Sumw2();
      sprintf(hname,"hMetE2_ch%i_%i",ich,isam);    hMetE2v[ich].push_back(new TH1D(hname,"",40,0,200));      hMetE2v[ich][isam]->Sumw2();
      sprintf(hname,"hMetPhi_ch%i_%i",ich,isam);   hMetPhiv[ich].push_back(new TH1D(hname,"",40,-3.2,3.2));  hMetPhiv[ich][isam]->Sumw2();
      sprintf(hname,"hMetPhiB_ch%i_%i",ich,isam);  hMetPhiBv[ich].push_back(new TH1D(hname,"",40,-3.2,3.2)); hMetPhiBv[ich][isam]->Sumw2();
      sprintf(hname,"hMetPhiE_ch%i_%i",ich,isam);  hMetPhiEv[ich].push_back(new TH1D(hname,"",40,-3.2,3.2)); hMetPhiEv[ich][isam]->Sumw2();
      sprintf(hname,"hMt_ch%i_%i",ich,isam);       hMtv[ich].push_back(new TH1D(hname,"",40,0,120));         hMtv[ich][isam]->Sumw2();
      sprintf(hname,"hMtB_ch%i_%i",ich,isam);      hMtBv[ich].push_back(new TH1D(hname,"",40,0,120));        hMtBv[ich][isam]->Sumw2();
      sprintf(hname,"hMtE_ch%i_%i",ich,isam);      hMtEv[ich].push_back(new TH1D(hname,"",40,0,120));        hMtEv[ich][isam]->Sumw2();
      sprintf(hname,"hMt2_ch%i_%i",ich,isam);      hMt2v[ich].push_back(new TH1D(hname,"",40,0,200));        hMt2v[ich][isam]->Sumw2();
      sprintf(hname,"hMtB2_ch%i_%i",ich,isam);     hMtB2v[ich].push_back(new TH1D(hname,"",40,0,200));       hMtB2v[ich][isam]->Sumw2();
      sprintf(hname,"hMtE2_ch%i_%i",ich,isam);     hMtE2v[ich].push_back(new TH1D(hname,"",40,0,200));       hMtE2v[ich][isam]->Sumw2();
      sprintf(hname,"hNPV_ch%i_%i",ich,isam);      hNPVv[ich].push_back(new TH1D(hname,"",30,-0.5,29.5));    hNPVv[ich][isam]->Sumw2();
    }
    sprintf(hname,"hPtMC_ch%i",ich);      hPtMC[ich] = new TH1D(hname,"",30,20,80);         hPtMC[ich]->Sumw2();
    sprintf(hname,"hPtBMC_ch%i",ich);	  hPtBMC[ich] = new TH1D(hname,"",30,20,80);	    hPtBMC[ich]->Sumw2();
    sprintf(hname,"hPtEMC_ch%i",ich);	  hPtEMC[ich] = new TH1D(hname,"",30,20,80);	    hPtEMC[ich]->Sumw2();
    sprintf(hname,"hPt2MC_ch%i",ich);	  hPt2MC[ich] = new TH1D(hname,"",40,20,180);	    hPt2MC[ich]->Sumw2();
    sprintf(hname,"hPtB2MC_ch%i",ich);	  hPtB2MC[ich] = new TH1D(hname,"",40,20,180);	    hPtB2MC[ich]->Sumw2();
    sprintf(hname,"hPtE2MC_ch%i",ich);	  hPtE2MC[ich] = new TH1D(hname,"",40,20,180);	    hPtE2MC[ich]->Sumw2();
    sprintf(hname,"hEtaMC_ch%i",ich);	  hEtaMC[ich] = new TH1D(hname,"",40,-3,3);	    hEtaMC[ich]->Sumw2();
    sprintf(hname,"hPhiMC_ch%i",ich);	  hPhiMC[ich] = new TH1D(hname,"",40,-3.2,3.2);	    hPhiMC[ich]->Sumw2();
    sprintf(hname,"hPhiBMC_ch%i",ich);	  hPhiBMC[ich] = new TH1D(hname,"",40,-3.2,3.2);    hPhiBMC[ich]->Sumw2();
    sprintf(hname,"hPhiEMC_ch%i",ich);	  hPhiEMC[ich] = new TH1D(hname,"",40,-3.2,3.2);    hPhiEMC[ich]->Sumw2();
    sprintf(hname,"hMetMC_ch%i",ich);	  hMetMC[ich] = new TH1D(hname,"",40,0,80);	    hMetMC[ich]->Sumw2();
    sprintf(hname,"hMetBMC_ch%i",ich);	  hMetBMC[ich] = new TH1D(hname,"",40,0,80);	    hMetBMC[ich]->Sumw2();
    sprintf(hname,"hMetEMC_ch%i",ich);	  hMetEMC[ich] = new TH1D(hname,"",40,0,80);	    hMetEMC[ich]->Sumw2();
    sprintf(hname,"hMet2MC_ch%i",ich);	  hMet2MC[ich] = new TH1D(hname,"",40,0,200);	    hMet2MC[ich]->Sumw2();
    sprintf(hname,"hMetB2MC_ch%i",ich);	  hMetB2MC[ich] = new TH1D(hname,"",40,0,200);	    hMetB2MC[ich]->Sumw2();
    sprintf(hname,"hMetE2MC_ch%i",ich);	  hMetE2MC[ich] = new TH1D(hname,"",40,0,200);	    hMetE2MC[ich]->Sumw2();
    sprintf(hname,"hMetPhiMC_ch%i",ich);  hMetPhiMC[ich] = new TH1D(hname,"",40,-3.2,3.2);  hMetPhiMC[ich]->Sumw2();
    sprintf(hname,"hMetPhiBMC_ch%i",ich); hMetPhiBMC[ich] = new TH1D(hname,"",40,-3.2,3.2); hMetPhiBMC[ich]->Sumw2();
    sprintf(hname,"hMetPhiEMC_ch%i",ich); hMetPhiEMC[ich] = new TH1D(hname,"",40,-3.2,3.2); hMetPhiEMC[ich]->Sumw2();
    sprintf(hname,"hMtMC_ch%i",ich);	  hMtMC[ich] = new TH1D(hname,"",40,0,120);	    hMtMC[ich]->Sumw2();
    sprintf(hname,"hMtBMC_ch%i",ich);	  hMtBMC[ich] = new TH1D(hname,"",40,0,120);	    hMtBMC[ich]->Sumw2();
    sprintf(hname,"hMtEMC_ch%i",ich);	  hMtEMC[ich] = new TH1D(hname,"",40,0,120);	    hMtEMC[ich]->Sumw2();
    sprintf(hname,"hMt2MC_ch%i",ich);	  hMt2MC[ich] = new TH1D(hname,"",40,0,200);	    hMt2MC[ich]->Sumw2();
    sprintf(hname,"hMtB2MC_ch%i",ich);	  hMtB2MC[ich] = new TH1D(hname,"",40,0,200);	    hMtB2MC[ich]->Sumw2();
    sprintf(hname,"hMtE2MC_ch%i",ich);	  hMtE2MC[ich] = new TH1D(hname,"",40,0,200);	    hMtE2MC[ich]->Sumw2();
    sprintf(hname,"hNPVMC_ch%i",ich);	  hNPVMC[ich] = new TH1D(hname,"",30,-0.5,29.5);    hNPVMC[ich]->Sumw2();
  }
  
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t genLepPt, genLepPhi;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  LorentzVector *lep=0;
  ///// electron specific /////
  Float_t trkIso, emIso, hadIso;
  Float_t pfChIso, pfGamIso, pfNeuIso, pfCombIso;
  Float_t sigieie, hovere, eoverp, fbrem, ecalE;
  Float_t dphi, deta;
  Float_t d0, dz;
  UInt_t  isConv, nexphits, typeBits;
  LorentzVector *sc=0;
  
  TFile *infile=0;
  TTree *intree=0;
        
  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if(isam==0 && !hasData) continue;    
    
    // Read input file and get the TTrees
    TString infilename = inputDir + TString("/") + snamev[isam] + TString("_select.root");
    cout << "Processing " << infilename << "..." << endl;
    infile = new TFile(infilename);	    assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);

    intree->SetBranchAddress("runNum",   &runNum);     // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);    // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);     // event number
    intree->SetBranchAddress("npv",      &npv);        // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);        // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);     // GEN boson pT
    intree->SetBranchAddress("genVPhi",  &genVPhi);    // GEN boson phi
    intree->SetBranchAddress("genVy",    &genVy);      // GEN boson rapidity
    intree->SetBranchAddress("genVMass", &genVMass);   // GEN boson mass
    intree->SetBranchAddress("genLepPt", &genLepPt);   // GEN lepton pT
    intree->SetBranchAddress("genLepPhi",&genLepPhi);  // GEN lepton phi
    intree->SetBranchAddress("scale1fb", &scale1fb);   // event weight per 1/fb (MC)
    intree->SetBranchAddress("met",      &met);        // MET
    intree->SetBranchAddress("metPhi",   &metPhi);     // phi(MET)
    intree->SetBranchAddress("sumEt",    &sumEt);      // Sum ET
    intree->SetBranchAddress("mt",       &mt);         // transverse mass
    intree->SetBranchAddress("u1",       &u1);         // parallel component of recoil
    intree->SetBranchAddress("u2",       &u2);         // perpendicular component of recoil
    intree->SetBranchAddress("q",        &q);          // lepton charge
    intree->SetBranchAddress("lep",      &lep);        // lepton 4-vector
    ///// electron specific /////
    intree->SetBranchAddress("trkIso",    &trkIso);     // track isolation of tag lepton
    intree->SetBranchAddress("emIso",     &emIso);      // ECAL isolation of tag lepton
    intree->SetBranchAddress("hadIso",    &hadIso);     // HCAL isolation of tag lepton
    intree->SetBranchAddress("pfChIso",   &pfChIso);    // PF charged hadron isolation of lepton
    intree->SetBranchAddress("pfGamIso",  &pfGamIso);	// PF photon isolation of lepton
    intree->SetBranchAddress("pfNeuIso",  &pfNeuIso);	// PF neutral hadron isolation of lepton
    intree->SetBranchAddress("pfCombIso", &pfCombIso);  // PF combined isolation of electron
    intree->SetBranchAddress("sigieie",   &sigieie);    // sigma-ieta-ieta of electron
    intree->SetBranchAddress("hovere",    &hovere);	// H/E of electron
    intree->SetBranchAddress("eoverp",    &eoverp);	// E/p of electron
    intree->SetBranchAddress("fbrem",     &fbrem);      // brem fraction of electron
    intree->SetBranchAddress("dphi",      &dphi);	// GSF track - ECAL dphi of electron
    intree->SetBranchAddress("deta",      &deta);	// GSF track - ECAL deta of electron
    intree->SetBranchAddress("ecalE",     &ecalE);      // ECAL energy of electron
    intree->SetBranchAddress("d0",        &d0); 	// transverse impact parameter of electron
    intree->SetBranchAddress("dz",        &dz); 	// longitudinal impact parameter of electron
    intree->SetBranchAddress("isConv",    &isConv);     // conversion filter flag of electron
    intree->SetBranchAddress("nexphits",  &nexphits);	// number of missing expected inner hits of electron
    intree->SetBranchAddress("typeBits",  &typeBits);	// electron type of electron
    intree->SetBranchAddress("sc",        &sc);         // electron Supercluster 4-vector

    
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      
      if(sc->Pt()        < PT_CUT)  continue;	
      if(fabs(sc->Eta()) > ETA_CUT) continue;

      Double_t weight = 1;
      if(isam!=0) {
        weight *= scale1fb*lumi;
      }
      
      for(UInt_t ich=0; ich<3; ich++) {
        if(ich==eInc || (ich==ePos && q>0) || (ich==eNeg && q<0)) {
	  hPtv[ich][isam]    ->Fill(lep->Pt(), weight);
	  hPt2v[ich][isam]   ->Fill(lep->Pt(), weight);
	  hEtav[ich][isam]   ->Fill(lep->Eta(),weight);
	  hPhiv[ich][isam]   ->Fill(lep->Phi(),weight);
	  hMetv[ich][isam]   ->Fill(met,       weight);
	  hMet2v[ich][isam]  ->Fill(met,       weight);
	  hMetPhiv[ich][isam]->Fill(metPhi,    weight);
	  hMtv[ich][isam]    ->Fill(mt,        weight);
	  hMt2v[ich][isam]   ->Fill(mt,        weight);
	  hNPVv[ich][isam]   ->Fill(npv,       weight);
	  if(isam!=0) {
	    hPtMC[ich]    ->Fill(lep->Pt(), weight);
	    hPt2MC[ich]   ->Fill(lep->Pt(), weight);
	    hEtaMC[ich]   ->Fill(lep->Eta(),weight);
	    hPhiMC[ich]   ->Fill(lep->Phi(),weight);
	    hMetMC[ich]   ->Fill(met,       weight);
	    hMet2MC[ich]  ->Fill(met,       weight);
	    hMetPhiMC[ich]->Fill(metPhi,    weight);
	    hMtMC[ich]    ->Fill(mt,        weight);
	    hMt2MC[ich]   ->Fill(mt,        weight);
	    hNPVMC[ich]   ->Fill(npv,       weight);
	  }
	  
	  if(fabs(lep->Eta())<ETA_BARREL) {
	    hPtBv[ich][isam]    ->Fill(lep->Pt(), weight);
	    hPtB2v[ich][isam]   ->Fill(lep->Pt(), weight);
	    hPhiBv[ich][isam]   ->Fill(lep->Phi(),weight);
	    hMetBv[ich][isam]   ->Fill(met,       weight);
	    hMetB2v[ich][isam]  ->Fill(met,       weight);
	    hMetPhiBv[ich][isam]->Fill(metPhi,    weight);
	    hMtBv[ich][isam]    ->Fill(mt,        weight);
	    hMtB2v[ich][isam]   ->Fill(mt,        weight);
	    if(isam!=0) {
	      hPtBMC[ich]    ->Fill(lep->Pt(), weight);
	      hPtB2MC[ich]   ->Fill(lep->Pt(), weight);
	      hPhiBMC[ich]   ->Fill(lep->Phi(),weight);
	      hMetBMC[ich]   ->Fill(met,       weight);
	      hMetB2MC[ich]  ->Fill(met,       weight);
	      hMetPhiBMC[ich]->Fill(metPhi,    weight);
	      hMtBMC[ich]    ->Fill(mt,        weight);
	      hMtB2MC[ich]   ->Fill(mt,        weight);
	    }
	    
	  } else if(fabs(lep->Eta())>ETA_ENDCAP) {
	    hPtEv[ich][isam]    ->Fill(lep->Pt(), weight);
	    hPtE2v[ich][isam]   ->Fill(lep->Pt(), weight);
	    hPhiEv[ich][isam]   ->Fill(lep->Phi(),weight);
	    hMetEv[ich][isam]   ->Fill(met,       weight);
	    hMetE2v[ich][isam]  ->Fill(met,       weight);
	    hMetPhiEv[ich][isam]->Fill(metPhi,    weight);
	    hMtEv[ich][isam]    ->Fill(mt,        weight);
	    hMtE2v[ich][isam]   ->Fill(mt,        weight);
	    if(isam!=0) {
	      hPtEMC[ich]    ->Fill(lep->Pt(), weight);
	      hPtE2MC[ich]   ->Fill(lep->Pt(), weight);
	      hPhiEMC[ich]   ->Fill(lep->Phi(),weight);
	      hMetEMC[ich]   ->Fill(met,       weight);
	      hMetE2MC[ich]  ->Fill(met,       weight);
	      hMetPhiEMC[ich]->Fill(metPhi,    weight);
	      hMtEMC[ich]    ->Fill(mt,        weight);
	      hMtE2MC[ich]   ->Fill(mt,        weight);
	    }
	  }
	}
      }
    }
    delete infile;
    infile=0, intree=0;    
  }  

  //
  // Make ratio graphs
  //
  vector<TH1D*> hPtDiffv, hPtBDiffv, hPtEDiffv;
  vector<TH1D*> hPt2Diffv, hPtB2Diffv, hPtE2Diffv;
  vector<TH1D*> hEtaDiffv;
  vector<TH1D*> hPhiDiffv, hPhiBDiffv, hPhiEDiffv;
  vector<TH1D*> hMetDiffv, hMetBDiffv, hMetEDiffv;
  vector<TH1D*> hMet2Diffv, hMetB2Diffv, hMetE2Diffv;
  vector<TH1D*> hMetPhiDiffv, hMetPhiBDiffv, hMetPhiEDiffv;
  vector<TH1D*> hMtDiffv, hMtBDiffv, hMtEDiffv;
  vector<TH1D*> hMt2Diffv, hMtB2Diffv, hMtE2Diffv;
  vector<TH1D*> hNPVDiffv;
  for(Int_t ich=0; ich<3; ich++) {
    sprintf(hname,"hPtDiff_ch%i",ich);      hPtDiffv.push_back(makeDiffHist(hPtv[ich][0],hPtMC[ich],hname));
    sprintf(hname,"hPtBDiff_ch%i",ich);     hPtBDiffv.push_back(makeDiffHist(hPtBv[ich][0],hPtBMC[ich],hname));
    sprintf(hname,"hPtEDiff_ch%i",ich);     hPtEDiffv.push_back(makeDiffHist(hPtEv[ich][0],hPtEMC[ich],hname));
    sprintf(hname,"hPt2Diff_ch%i",ich);     hPt2Diffv.push_back(makeDiffHist(hPt2v[ich][0],hPt2MC[ich],hname));
    sprintf(hname,"hPtB2Diff_ch%i",ich);    hPtB2Diffv.push_back(makeDiffHist(hPtB2v[ich][0],hPtB2MC[ich],hname));
    sprintf(hname,"hPtE2Diff_ch%i",ich);    hPtE2Diffv.push_back(makeDiffHist(hPtE2v[ich][0],hPtE2MC[ich],hname));
    sprintf(hname,"hEtaDiff_ch%i",ich);     hEtaDiffv.push_back(makeDiffHist(hEtav[ich][0],hEtaMC[ich],hname));
    sprintf(hname,"hPhiDiff_ch%i",ich);     hPhiDiffv.push_back(makeDiffHist(hPhiv[ich][0],hPhiMC[ich],hname));
    sprintf(hname,"hPhiBDiff_ch%i",ich);    hPhiBDiffv.push_back(makeDiffHist(hPhiBv[ich][0],hPhiBMC[ich],hname));
    sprintf(hname,"hPhiEDiff_ch%i",ich);    hPhiEDiffv.push_back(makeDiffHist(hPhiEv[ich][0],hPhiEMC[ich],hname));
    sprintf(hname,"hMetDiff_ch%i",ich);     hMetDiffv.push_back(makeDiffHist(hMetv[ich][0],hMetMC[ich],hname));
    sprintf(hname,"hMetBDiff_ch%i",ich);    hMetBDiffv.push_back(makeDiffHist(hMetBv[ich][0],hMetBMC[ich],hname));
    sprintf(hname,"hMetEDiff_ch%i",ich);    hMetEDiffv.push_back(makeDiffHist(hMetEv[ich][0],hMetEMC[ich],hname));
    sprintf(hname,"hMet2Diff_ch%i",ich);    hMet2Diffv.push_back(makeDiffHist(hMet2v[ich][0],hMet2MC[ich],hname));
    sprintf(hname,"hMetB2Diff_ch%i",ich);   hMetB2Diffv.push_back(makeDiffHist(hMetB2v[ich][0],hMetB2MC[ich],hname));
    sprintf(hname,"hMetE2Diff_ch%i",ich);   hMetE2Diffv.push_back(makeDiffHist(hMetE2v[ich][0],hMetE2MC[ich],hname));
    sprintf(hname,"hMetPhiDiff_ch%i",ich);  hMetPhiDiffv.push_back(makeDiffHist(hMetPhiv[ich][0],hMetPhiMC[ich],hname));
    sprintf(hname,"hMetPhiBDiff_ch%i",ich); hMetPhiBDiffv.push_back(makeDiffHist(hMetPhiBv[ich][0],hMetPhiBMC[ich],hname));
    sprintf(hname,"hMetPhiEDiff_ch%i",ich); hMetPhiEDiffv.push_back(makeDiffHist(hMetPhiEv[ich][0],hMetPhiEMC[ich],hname));
    sprintf(hname,"hMtDiff_ch%i",ich);      hMtDiffv.push_back(makeDiffHist(hMtv[ich][0],hMtMC[ich],hname));
    sprintf(hname,"hMtBDiff_ch%i",ich);     hMtBDiffv.push_back(makeDiffHist(hMtBv[ich][0],hMtBMC[ich],hname));
    sprintf(hname,"hMtEDiff_ch%i",ich);     hMtEDiffv.push_back(makeDiffHist(hMtEv[ich][0],hMtEMC[ich],hname));
    sprintf(hname,"hMt2Diff_ch%i",ich);     hMt2Diffv.push_back(makeDiffHist(hMt2v[ich][0],hMt2MC[ich],hname));
    sprintf(hname,"hMtB2Diff_ch%i",ich);    hMtB2Diffv.push_back(makeDiffHist(hMtB2v[ich][0],hMtB2MC[ich],hname));
    sprintf(hname,"hMtE2Diff_ch%i",ich);    hMtE2Diffv.push_back(makeDiffHist(hMtE2v[ich][0],hMtE2MC[ich],hname));
    sprintf(hname,"hNPVDiff_ch%i",ich);     hNPVDiffv.push_back(makeDiffHist(hNPVv[ich][0],hNPVMC[ich],hname));
  }

  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  // string buffers
  char pname[100];       // plot name
  char xlabel[100];      // x-axis label
  char ylabel[100];      // y-axis label
  char lumitext[100];    // lumi label
  char barreltext[100];  // barrel label
  char endcaptext[100];  // endcap label
  char eventtext[100];   // event label

  if(lumi<0.1)    sprintf(lumitext,"#int#font[12]{L}dt = %.1f pb^{-1}",1000.*lumi);
  else if(lumi<1) sprintf(lumitext,"#int#font[12]{L}dt = %.0f pb^{-1}",1000.*lumi);
  else            sprintf(lumitext,"#int#font[12]{L}dt = %.2f fb^{-1}",lumi);
  
  sprintf(barreltext,"|#eta| < %.1f",ETA_BARREL);
  sprintf(endcaptext,"|#eta| > %.1f",ETA_ENDCAP);
  
  TCanvas *c = MakeCanvas("c","c",800,800);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.01);
  c->cd(1)->SetLeftMargin(0.18);  
  c->cd(1)->SetRightMargin(0.07);  
  c->cd(1)->SetTickx(1);
  c->cd(1)->SetTicky(1);  
  c->cd(2)->SetPad(0,0,1.0,0.3);
  c->cd(2)->SetTopMargin(0.05);
  c->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.18);
  c->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetTickx(1);
  c->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.400,"Y");
  
  Int_t ratioColor = kBlue;

  for(UInt_t ich=0; ich<3; ich++) {    
    
    //
    // Lepton pT
    //
    if(ich==eInc) sprintf(xlabel,"p_{T}(e^{#pm}) [GeV/c]");
    if(ich==ePos) sprintf(xlabel,"p_{T}(e^{+}) [GeV/c]");
    if(ich==eNeg) sprintf(xlabel,"p_{T}(e^{-}) [GeV/c]");
    
    sprintf(ylabel,"Events / %.1f GeV/c",hPtv[ich][0]->GetBinWidth(1));
    sprintf(pname,"lpt_ch%i",ich);
    CPlot plotPt(pname,"",xlabel,ylabel);
    plotPt.AddHist1D(hPtv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPt.AddToStack(hPtv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPt.SetLegend(0.75,0.6,0.95,0.9); 
    plotPt.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);   
    plotPt.Draw(c,kTRUE,format,1);
    
    CPlot plotPtGraph(pname,"",xlabel,"#chi");
    plotPtGraph.AddHist1D(hPtDiffv[ich],"E",ratioColor);
    plotPtGraph.SetYRange(-8,8);
    plotPtGraph.AddLine(hPtDiffv[ich]->GetXaxis()->GetXmin(), 0,hPtDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPtGraph.AddLine(hPtDiffv[ich]->GetXaxis()->GetXmin(), 5,hPtDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPtGraph.AddLine(hPtDiffv[ich]->GetXaxis()->GetXmin(),-5,hPtDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotPtGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c",hPtBv[ich][0]->GetBinWidth(1));
    sprintf(pname,"lpt_barrel_ch%i",ich);
    CPlot plotPtB(pname,"",xlabel,ylabel);
    plotPtB.AddHist1D(hPtBv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPtB.AddToStack(hPtBv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPtB.SetLegend(0.75,0.6,0.95,0.9); 
    plotPtB.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0); 
    plotPtB.AddTextBox(barreltext,0.50,0.75,0.70,0.69,0);  
    plotPtB.Draw(c,kFALSE,format,1);
    
    CPlot plotPtBGraph(pname,"",xlabel,"#chi");
    plotPtBGraph.AddHist1D(hPtBDiffv[ich],"E",ratioColor);
    plotPtBGraph.SetYRange(-8,8);
    plotPtBGraph.AddLine(hPtBDiffv[ich]->GetXaxis()->GetXmin(), 0,hPtBDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPtBGraph.AddLine(hPtBDiffv[ich]->GetXaxis()->GetXmin(), 5,hPtBDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPtBGraph.AddLine(hPtBDiffv[ich]->GetXaxis()->GetXmin(),-5,hPtBDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotPtBGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c",hPtEv[ich][0]->GetBinWidth(1));
    sprintf(pname,"lpt_endcap_ch%i",ich);
    CPlot plotPtE(pname,"",xlabel,ylabel);
    plotPtE.AddHist1D(hPtEv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPtE.AddToStack(hPtEv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPtE.SetLegend(0.75,0.6,0.95,0.9); 
    plotPtE.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);
    plotPtE.AddTextBox(endcaptext,0.50,0.75,0.70,0.69,0);
    plotPtE.Draw(c,kFALSE,format,1);
    
    CPlot plotPtEGraph(pname,"",xlabel,"#chi");
    plotPtEGraph.AddHist1D(hPtEDiffv[ich],"E",ratioColor);
    plotPtEGraph.SetYRange(-8,8);
    plotPtEGraph.AddLine(hPtEDiffv[ich]->GetXaxis()->GetXmin(), 0,hPtEDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPtEGraph.AddLine(hPtEDiffv[ich]->GetXaxis()->GetXmin(), 5,hPtEDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPtEGraph.AddLine(hPtEDiffv[ich]->GetXaxis()->GetXmin(),-5,hPtEDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotPtEGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c",hPt2v[ich][0]->GetBinWidth(1));
    sprintf(pname,"lptlog_ch%i",ich);
    CPlot plotPt2(pname,"",xlabel,ylabel);
    plotPt2.AddHist1D(hPt2v[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPt2.AddToStack(hPt2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPt2.SetLegend(0.75,0.6,0.95,0.9); 
    plotPt2.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);
    plotPt2.SetLogy();   
    plotPt2.Draw(c,kFALSE,format,1);
    
    CPlot plotPt2Graph(pname,"",xlabel,"#chi");
    plotPt2Graph.AddHist1D(hPt2Diffv[ich],"E",ratioColor);
    plotPt2Graph.SetYRange(-8,8);
    plotPt2Graph.AddLine(hPt2Diffv[ich]->GetXaxis()->GetXmin(), 0,hPt2Diffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPt2Graph.AddLine(hPt2Diffv[ich]->GetXaxis()->GetXmin(), 5,hPt2Diffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPt2Graph.AddLine(hPt2Diffv[ich]->GetXaxis()->GetXmin(),-5,hPt2Diffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotPt2Graph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c",hPtB2v[ich][0]->GetBinWidth(1));
    sprintf(pname,"lptlog_barrel_ch%i",ich);
    CPlot plotPtB2(pname,"",xlabel,ylabel);
    plotPtB2.AddHist1D(hPtB2v[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPtB2.AddToStack(hPtB2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPtB2.SetLegend(0.75,0.6,0.95,0.9); 
    plotPtB2.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0); 
    plotPtB2.AddTextBox(barreltext,0.50,0.75,0.70,0.69,0);
    plotPtB2.SetLogy();
    plotPtB2.Draw(c,kFALSE,format,1);
    
    CPlot plotPtB2Graph(pname,"",xlabel,"#chi");
    plotPtB2Graph.AddHist1D(hPtB2Diffv[ich],"E",ratioColor);
    plotPtB2Graph.SetYRange(-8,8);
    plotPtB2Graph.AddLine(hPtB2Diffv[ich]->GetXaxis()->GetXmin(), 0,hPtB2Diffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPtB2Graph.AddLine(hPtB2Diffv[ich]->GetXaxis()->GetXmin(), 5,hPtB2Diffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPtB2Graph.AddLine(hPtB2Diffv[ich]->GetXaxis()->GetXmin(),-5,hPtB2Diffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotPtB2Graph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c",hPtE2v[ich][0]->GetBinWidth(1));
    sprintf(pname,"lptlog_endcap_ch%i",ich);
    CPlot plotPtE2(pname,"",xlabel,ylabel);
    plotPtE2.AddHist1D(hPtE2v[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPtE2.AddToStack(hPtE2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPtE2.SetLegend(0.75,0.6,0.95,0.9); 
    plotPtE2.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);
    plotPtE2.AddTextBox(endcaptext,0.50,0.75,0.70,0.69,0);
    plotPtE2.SetLogy();
    plotPtE2.Draw(c,kFALSE,format,1);
    
    CPlot plotPtE2Graph(pname,"",xlabel,"#chi");
    plotPtE2Graph.AddHist1D(hPtE2Diffv[ich],"E",ratioColor);
    plotPtE2Graph.SetYRange(-8,8);
    plotPtE2Graph.AddLine(hPtE2Diffv[ich]->GetXaxis()->GetXmin(), 0,hPtE2Diffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPtE2Graph.AddLine(hPtE2Diffv[ich]->GetXaxis()->GetXmin(), 5,hPtE2Diffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPtE2Graph.AddLine(hPtE2Diffv[ich]->GetXaxis()->GetXmin(),-5,hPtE2Diffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotPtE2Graph.Draw(c,kTRUE,format,2);
    
    //
    // lepton eta
    //
    if(ich==eInc) sprintf(xlabel,"#eta(e^{#pm})");
    if(ich==ePos)  sprintf(xlabel,"#eta(e^{+})");
    if(ich==eNeg)  sprintf(xlabel,"#eta(e^{-})");
    
    sprintf(ylabel,"Events / %.2f",hEtav[ich][0]->GetBinWidth(1));
    sprintf(pname,"leta_ch%i",ich);
    CPlot plotEta(pname,"",xlabel,ylabel);
    plotEta.AddHist1D(hEtav[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotEta.AddToStack(hEtav[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotEta.SetLegend(0.75,0.6,0.95,0.9); 
    plotEta.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);   
    plotEta.SetYRange(0,1.8*(hEtav[ich][0]->GetMaximum()));
    plotEta.Draw(c,kFALSE,format,1);
    
    CPlot plotEtaGraph(pname,"",xlabel,"#chi");
    plotEtaGraph.AddHist1D(hEtaDiffv[ich],"E",ratioColor);
    plotEtaGraph.SetYRange(-8,8);
    plotEtaGraph.AddLine(hEtaDiffv[ich]->GetXaxis()->GetXmin(), 0,hEtaDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotEtaGraph.AddLine(hEtaDiffv[ich]->GetXaxis()->GetXmin(), 5,hEtaDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotEtaGraph.AddLine(hEtaDiffv[ich]->GetXaxis()->GetXmin(),-5,hEtaDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotEtaGraph.Draw(c,kTRUE,format,2);
    
    //
    // Lepton phi
    //
    if(ich==eInc) sprintf(xlabel,"#phi(e^{#pm})");
    if(ich==ePos)  sprintf(xlabel,"#phi(e^{+})");
    if(ich==eNeg)  sprintf(xlabel,"#phi(e^{-})");
    
    sprintf(ylabel,"Events / %.2f",hPhiv[ich][0]->GetBinWidth(1));
    sprintf(pname,"lphi_ch%i",ich);
    CPlot plotPhi(pname,"",xlabel,ylabel);
    plotPhi.AddHist1D(hPhiv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPhi.AddToStack(hPhiv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPhi.SetLegend(0.75,0.6,0.95,0.9); 
    plotPhi.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);   
    plotPhi.SetYRange(0,1.8*(hPhiv[ich][0]->GetMaximum()));
    plotPhi.Draw(c,kFALSE,format,1);
    
    CPlot plotPhiGraph(pname,"",xlabel,"#chi");
    plotPhiGraph.AddHist1D(hPhiDiffv[ich],"E",ratioColor);
    plotPhiGraph.SetYRange(-8,8);
    plotPhiGraph.AddLine(hPhiDiffv[ich]->GetXaxis()->GetXmin(), 0,hPhiDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPhiGraph.AddLine(hPhiDiffv[ich]->GetXaxis()->GetXmin(), 5,hPhiDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPhiGraph.AddLine(hPhiDiffv[ich]->GetXaxis()->GetXmin(),-5,hPhiDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotPhiGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.2f",hPhiBv[ich][0]->GetBinWidth(1));
    sprintf(pname,"lphi_barrel_ch%i",ich);
    CPlot plotPhiB(pname,"",xlabel,ylabel);
    plotPhiB.AddHist1D(hPhiBv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPhiB.AddToStack(hPhiBv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPhiB.SetLegend(0.75,0.6,0.95,0.9); 
    plotPhiB.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0); 
    plotPhiB.AddTextBox(barreltext,0.50,0.75,0.70,0.69,0);  
    plotPhiB.SetYRange(0,1.8*(hPhiBv[ich][0]->GetMaximum()));
    plotPhiB.Draw(c,kFALSE,format,1);
    
    CPlot plotPhiBGraph(pname,"",xlabel,"#chi");
    plotPhiBGraph.AddHist1D(hPhiBDiffv[ich],"E",ratioColor);
    plotPhiBGraph.SetYRange(-8,8);
    plotPhiBGraph.AddLine(hPhiBDiffv[ich]->GetXaxis()->GetXmin(), 0,hPhiBDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPhiBGraph.AddLine(hPhiBDiffv[ich]->GetXaxis()->GetXmin(), 5,hPhiBDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPhiBGraph.AddLine(hPhiBDiffv[ich]->GetXaxis()->GetXmin(),-5,hPhiBDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotPhiBGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.2f",hPhiEv[ich][0]->GetBinWidth(1));
    sprintf(pname,"lphi_endcap_ch%i",ich);
    CPlot plotPhiE(pname,"",xlabel,ylabel);
    plotPhiE.AddHist1D(hPhiEv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPhiE.AddToStack(hPhiEv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPhiE.SetLegend(0.75,0.6,0.95,0.9); 
    plotPhiE.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);
    plotPhiE.AddTextBox(endcaptext,0.50,0.75,0.70,0.69,0);
    plotPhiE.SetYRange(0,1.8*(hPhiEv[ich][0]->GetMaximum()));
    plotPhiE.Draw(c,kFALSE,format,1);
    
    CPlot plotPhiEGraph(pname,"",xlabel,"#chi");
    plotPhiEGraph.AddHist1D(hPhiEDiffv[ich],"E",ratioColor);
    plotPhiEGraph.SetYRange(-8,8);
    plotPhiEGraph.AddLine(hPhiEDiffv[ich]->GetXaxis()->GetXmin(), 0,hPhiEDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPhiEGraph.AddLine(hPhiEDiffv[ich]->GetXaxis()->GetXmin(), 5,hPhiEDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPhiEGraph.AddLine(hPhiEDiffv[ich]->GetXaxis()->GetXmin(),-5,hPhiEDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotPhiEGraph.Draw(c,kTRUE,format,2);
    
    //
    // MET
    //
    sprintf(xlabel,"#slash{E}_{T} [GeV]");
    if(ich==eInc) sprintf(eventtext,"e^{#pm} events");
    if(ich==ePos) sprintf(eventtext,"e^{+} events");
    if(ich==eNeg) sprintf(eventtext,"e^{-} events");
    
    sprintf(ylabel,"Events / %.1f GeV/c",hMetv[ich][0]->GetBinWidth(1));
    sprintf(pname,"met_ch%i",ich);
    CPlot plotMet(pname,"",xlabel,ylabel);
    plotMet.AddHist1D(hMetv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMet.AddToStack(hMetv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMet.SetLegend(0.75,0.6,0.95,0.9); 
    plotMet.AddTextBox(lumitext,0.21,0.85,0.41,0.79,0);   
    plotMet.AddTextBox(eventtext,0.21,0.75,0.41,0.69,0);
    plotMet.Draw(c,kFALSE,format,1);
    
    CPlot plotMetGraph(pname,"",xlabel,"#chi");
    plotMetGraph.AddHist1D(hMetDiffv[ich],"E",ratioColor);
    plotMetGraph.SetYRange(-8,8);
    plotMetGraph.AddLine(hMetDiffv[ich]->GetXaxis()->GetXmin(), 0,hMetDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMetGraph.AddLine(hMetDiffv[ich]->GetXaxis()->GetXmin(), 5,hMetDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMetGraph.AddLine(hMetDiffv[ich]->GetXaxis()->GetXmin(),-5,hMetDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMetGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c",hMetBv[ich][0]->GetBinWidth(1));
    sprintf(pname,"met_barrel_ch%i",ich);
    CPlot plotMetB(pname,"",xlabel,ylabel);
    plotMetB.AddHist1D(hMetBv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMetB.AddToStack(hMetBv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMetB.SetLegend(0.75,0.6,0.95,0.9); 
    plotMetB.AddTextBox(lumitext,0.21,0.85,0.41,0.79,0); 
    plotMetB.AddTextBox(barreltext,0.21,0.75,0.41,0.69,0);  
    plotMetB.AddTextBox(eventtext,0.21,0.65,0.41,0.59,0);
    plotMetB.Draw(c,kFALSE,format,1);
    
    CPlot plotMetBGraph(pname,"",xlabel,"#chi");
    plotMetBGraph.AddHist1D(hMetBDiffv[ich],"E",ratioColor);
    plotMetBGraph.SetYRange(-8,8);
    plotMetBGraph.AddLine(hMetBDiffv[ich]->GetXaxis()->GetXmin(), 0,hMetBDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMetBGraph.AddLine(hMetBDiffv[ich]->GetXaxis()->GetXmin(), 5,hMetBDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMetBGraph.AddLine(hMetBDiffv[ich]->GetXaxis()->GetXmin(),-5,hMetBDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMetBGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c",hMetEv[ich][0]->GetBinWidth(1));
    sprintf(pname,"met_endcap_ch%i",ich);
    CPlot plotMetE(pname,"",xlabel,ylabel);
    plotMetE.AddHist1D(hMetEv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMetE.AddToStack(hMetEv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMetE.SetLegend(0.75,0.6,0.95,0.9); 
    plotMetE.AddTextBox(lumitext,0.21,0.85,0.41,0.79,0);
    plotMetE.AddTextBox(endcaptext,0.21,0.75,0.41,0.69,0);
    plotMetE.AddTextBox(eventtext,0.21,0.65,0.41,0.59,0);
    plotMetE.Draw(c,kFALSE,format,1);
    
    CPlot plotMetEGraph(pname,"",xlabel,"#chi");
    plotMetEGraph.AddHist1D(hMetEDiffv[ich],"E",ratioColor);
    plotMetEGraph.SetYRange(-8,8);
    plotMetEGraph.AddLine(hMetEDiffv[ich]->GetXaxis()->GetXmin(), 0,hMetEDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMetEGraph.AddLine(hMetEDiffv[ich]->GetXaxis()->GetXmin(), 5,hMetEDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMetEGraph.AddLine(hMetEDiffv[ich]->GetXaxis()->GetXmin(),-5,hMetEDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMetEGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c",hMet2v[ich][0]->GetBinWidth(1));
    sprintf(pname,"metlog_ch%i",ich);
    CPlot plotMet2(pname,"",xlabel,ylabel);
    plotMet2.AddHist1D(hMet2v[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMet2.AddToStack(hMet2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMet2.SetLegend(0.75,0.6,0.95,0.9); 
    plotMet2.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);
    plotMet2.AddTextBox(eventtext,0.50,0.75,0.70,0.69,0);
    plotMet2.SetLogy();   
    plotMet2.Draw(c,kFALSE,format,1);
    
    CPlot plotMet2Graph(pname,"",xlabel,"#chi");
    plotMet2Graph.AddHist1D(hMet2Diffv[ich],"E",ratioColor);
    plotMet2Graph.SetYRange(-8,8);
    plotMet2Graph.AddLine(hMet2Diffv[ich]->GetXaxis()->GetXmin(), 0,hMet2Diffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMet2Graph.AddLine(hMet2Diffv[ich]->GetXaxis()->GetXmin(), 5,hMet2Diffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMet2Graph.AddLine(hMet2Diffv[ich]->GetXaxis()->GetXmin(),-5,hMet2Diffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMet2Graph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c",hMetB2v[ich][0]->GetBinWidth(1));
    sprintf(pname,"metlog_barrel_ch%i",ich);
    CPlot plotMetB2(pname,"",xlabel,ylabel);
    plotMetB2.AddHist1D(hMetB2v[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMetB2.AddToStack(hMetB2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMetB2.SetLegend(0.75,0.6,0.95,0.9); 
    plotMetB2.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0); 
    plotMetB2.AddTextBox(barreltext,0.50,0.75,0.70,0.69,0);
    plotMetB2.AddTextBox(eventtext,0.50,0.65,0.70,0.59,0);
    plotMetB2.SetLogy();
    plotMetB2.Draw(c,kFALSE,format,1);
    
    CPlot plotMetB2Graph(pname,"",xlabel,"#chi");
    plotMetB2Graph.AddHist1D(hMetB2Diffv[ich],"E",ratioColor);
    plotMetB2Graph.SetYRange(-8,8);
    plotMetB2Graph.AddLine(hMetB2Diffv[ich]->GetXaxis()->GetXmin(), 0,hMetB2Diffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMetB2Graph.AddLine(hMetB2Diffv[ich]->GetXaxis()->GetXmin(), 5,hMetB2Diffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMetB2Graph.AddLine(hMetB2Diffv[ich]->GetXaxis()->GetXmin(),-5,hMetB2Diffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMetB2Graph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c",hMetE2v[ich][0]->GetBinWidth(1));
    sprintf(pname,"metlog_endcap_ch%i",ich);
    CPlot plotMetE2(pname,"",xlabel,ylabel);
    plotMetE2.AddHist1D(hMetE2v[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMetE2.AddToStack(hMetE2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMetE2.SetLegend(0.75,0.6,0.95,0.9); 
    plotMetE2.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);
    plotMetE2.AddTextBox(endcaptext,0.50,0.75,0.70,0.69,0);
    plotMetE2.AddTextBox(eventtext,0.50,0.65,0.70,0.59,0);
    plotMetE2.SetLogy();
    plotMetE2.Draw(c,kFALSE,format,1);
    
    CPlot plotMetE2Graph(pname,"",xlabel,"#chi");
    plotMetE2Graph.AddHist1D(hMetE2Diffv[ich],"E",ratioColor);
    plotMetE2Graph.SetYRange(-8,8);
    plotMetE2Graph.AddLine(hMetE2Diffv[ich]->GetXaxis()->GetXmin(), 0,hMetE2Diffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMetE2Graph.AddLine(hMetE2Diffv[ich]->GetXaxis()->GetXmin(), 5,hMetE2Diffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMetE2Graph.AddLine(hMetE2Diffv[ich]->GetXaxis()->GetXmin(),-5,hMetE2Diffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMetE2Graph.Draw(c,kTRUE,format,2);
  
    //
    // phi(MET)
    //
    sprintf(xlabel,"#phi(#slash{E}_{T})");
    if(ich==eInc) sprintf(eventtext,"e^{#pm} events");
    if(ich==ePos)  sprintf(eventtext,"e^{+} events");
    if(ich==eNeg)  sprintf(eventtext,"e^{-} events");
    
    sprintf(ylabel,"Events / %.2f",hMetPhiv[ich][0]->GetBinWidth(1));
    sprintf(pname,"metphi_ch%i",ich);
    CPlot plotMetPhi(pname,"",xlabel,ylabel);
    plotMetPhi.AddHist1D(hMetPhiv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMetPhi.AddToStack(hMetPhiv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMetPhi.SetLegend(0.75,0.6,0.95,0.9); 
    plotMetPhi.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);   
    plotMetPhi.AddTextBox(eventtext,0.50,0.75,0.70,0.69,0);
    plotMetPhi.SetYRange(0,1.8*(hMetPhiv[ich][0]->GetMaximum()));
    plotMetPhi.Draw(c,kFALSE,format,1);
    
    CPlot plotMetPhiGraph(pname,"",xlabel,"#chi");
    plotMetPhiGraph.AddHist1D(hMetPhiDiffv[ich],"E",ratioColor);
    plotMetPhiGraph.SetYRange(-8,8);
    plotMetPhiGraph.AddLine(hMetPhiDiffv[ich]->GetXaxis()->GetXmin(), 0,hMetPhiDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMetPhiGraph.AddLine(hMetPhiDiffv[ich]->GetXaxis()->GetXmin(), 5,hMetPhiDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMetPhiGraph.AddLine(hMetPhiDiffv[ich]->GetXaxis()->GetXmin(),-5,hMetPhiDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMetPhiGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.2f",hMetPhiBv[ich][0]->GetBinWidth(1));
    sprintf(pname,"metphi_barrel_ch%i",ich);
    CPlot plotMetPhiB(pname,"",xlabel,ylabel);
    plotMetPhiB.AddHist1D(hMetPhiBv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMetPhiB.AddToStack(hMetPhiBv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMetPhiB.SetLegend(0.75,0.6,0.95,0.9); 
    plotMetPhiB.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0); 
    plotMetPhiB.AddTextBox(barreltext,0.50,0.75,0.70,0.69,0);
    plotMetPhiB.AddTextBox(eventtext,0.50,0.65,0.70,0.59,0);
    plotMetPhiB.SetYRange(0,1.8*(hMetPhiBv[ich][0]->GetMaximum()));
    plotMetPhiB.Draw(c,kFALSE,format,1);
    
    CPlot plotMetPhiBGraph(pname,"",xlabel,"#chi");
    plotMetPhiBGraph.AddHist1D(hMetPhiBDiffv[ich],"E",ratioColor);
    plotMetPhiBGraph.SetYRange(-8,8);
    plotMetPhiBGraph.AddLine(hMetPhiBDiffv[ich]->GetXaxis()->GetXmin(), 0,hMetPhiBDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMetPhiBGraph.AddLine(hMetPhiBDiffv[ich]->GetXaxis()->GetXmin(), 5,hMetPhiBDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMetPhiBGraph.AddLine(hMetPhiBDiffv[ich]->GetXaxis()->GetXmin(),-5,hMetPhiBDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMetPhiBGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.2f",hMetPhiEv[ich][0]->GetBinWidth(1));
    sprintf(pname,"metphi_endcap_ch%i",ich);
    CPlot plotMetPhiE(pname,"",xlabel,ylabel);
    plotMetPhiE.AddHist1D(hMetPhiEv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMetPhiE.AddToStack(hMetPhiEv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMetPhiE.SetLegend(0.75,0.6,0.95,0.9); 
    plotMetPhiE.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);
    plotMetPhiE.AddTextBox(endcaptext,0.50,0.75,0.70,0.69,0);
    plotMetPhiE.AddTextBox(eventtext,0.50,0.65,0.70,0.59,0);
    plotMetPhiE.SetYRange(0,1.8*(hMetPhiEv[ich][0]->GetMaximum()));
    plotMetPhiE.Draw(c,kFALSE,format,1);
    
    CPlot plotMetPhiEGraph(pname,"",xlabel,"#chi");
    plotMetPhiEGraph.AddHist1D(hMetPhiEDiffv[ich],"E",ratioColor);
    plotMetPhiEGraph.SetYRange(-8,8);
    plotMetPhiEGraph.AddLine(hMetPhiEDiffv[ich]->GetXaxis()->GetXmin(), 0,hMetPhiEDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMetPhiEGraph.AddLine(hMetPhiEDiffv[ich]->GetXaxis()->GetXmin(), 5,hMetPhiEDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMetPhiEGraph.AddLine(hMetPhiEDiffv[ich]->GetXaxis()->GetXmin(),-5,hMetPhiEDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMetPhiEGraph.Draw(c,kTRUE,format,2);

    //
    // MT
    //
    sprintf(xlabel,"m_{T} [GeV/c^{2}]");
    if(ich==eInc) sprintf(eventtext,"e^{#pm} events");
    if(ich==ePos)  sprintf(eventtext,"e^{+} events");
    if(ich==eNeg)  sprintf(eventtext,"e^{-} events");
    
    sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMtv[ich][0]->GetBinWidth(1));
    sprintf(pname,"mt_ch%i",ich);
    CPlot plotMt(pname,"",xlabel,ylabel);
    plotMt.AddHist1D(hMtv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMt.AddToStack(hMtv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMt.SetLegend(0.75,0.6,0.95,0.9); 
    plotMt.AddTextBox(lumitext,0.21,0.85,0.41,0.79,0);   
    plotMt.AddTextBox(eventtext,0.21,0.75,0.41,0.69,0);
    plotMt.Draw(c,kFALSE,format,1);
    
    CPlot plotMtGraph(pname,"",xlabel,"#chi");
    plotMtGraph.AddHist1D(hMtDiffv[ich],"E",ratioColor);
    plotMtGraph.SetYRange(-8,8);
    plotMtGraph.AddLine(hMtDiffv[ich]->GetXaxis()->GetXmin(), 0,hMtDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMtGraph.AddLine(hMtDiffv[ich]->GetXaxis()->GetXmin(), 5,hMtDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMtGraph.AddLine(hMtDiffv[ich]->GetXaxis()->GetXmin(),-5,hMtDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMtGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMtBv[ich][0]->GetBinWidth(1));
    sprintf(pname,"mt_barrel_ch%i",ich);
    CPlot plotMtB(pname,"",xlabel,ylabel);
    plotMtB.AddHist1D(hMtBv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMtB.AddToStack(hMtBv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMtB.SetLegend(0.75,0.6,0.95,0.9); 
    plotMtB.AddTextBox(lumitext,0.21,0.85,0.41,0.79,0); 
    plotMtB.AddTextBox(barreltext,0.21,0.75,0.41,0.69,0);  
    plotMtB.AddTextBox(eventtext,0.21,0.65,0.41,0.59,0);
    plotMtB.Draw(c,kFALSE,format,1);
    
    CPlot plotMtBGraph(pname,"",xlabel,"#chi");
    plotMtBGraph.AddHist1D(hMtBDiffv[ich],"E",ratioColor);
    plotMtBGraph.SetYRange(-8,8);
    plotMtBGraph.AddLine(hMtBDiffv[ich]->GetXaxis()->GetXmin(), 0,hMtBDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMtBGraph.AddLine(hMtBDiffv[ich]->GetXaxis()->GetXmin(), 5,hMtBDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMtBGraph.AddLine(hMtBDiffv[ich]->GetXaxis()->GetXmin(),-5,hMtBDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMtBGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMtEv[ich][0]->GetBinWidth(1));
    sprintf(pname,"mt_endcap_ch%i",ich);
    CPlot plotMtE(pname,"",xlabel,ylabel);
    plotMtE.AddHist1D(hMtEv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMtE.AddToStack(hMtEv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMtE.SetLegend(0.75,0.6,0.95,0.9); 
    plotMtE.AddTextBox(lumitext,0.21,0.85,0.41,0.79,0);
    plotMtE.AddTextBox(endcaptext,0.21,0.75,0.41,0.69,0);
    plotMtE.AddTextBox(eventtext,0.21,0.65,0.41,0.59,0);
    plotMtE.Draw(c,kFALSE,format,1);
    
    CPlot plotMtEGraph(pname,"",xlabel,"#chi");
    plotMtEGraph.AddHist1D(hMtEDiffv[ich],"E",ratioColor);
    plotMtEGraph.SetYRange(-8,8);
    plotMtEGraph.AddLine(hMtEDiffv[ich]->GetXaxis()->GetXmin(), 0,hMtEDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMtEGraph.AddLine(hMtEDiffv[ich]->GetXaxis()->GetXmin(), 5,hMtEDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMtEGraph.AddLine(hMtEDiffv[ich]->GetXaxis()->GetXmin(),-5,hMtEDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMtEGraph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMt2v[ich][0]->GetBinWidth(1));
    sprintf(pname,"mtlog_ch%i",ich);
    CPlot plotMt2(pname,"",xlabel,ylabel);
    plotMt2.AddHist1D(hMt2v[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMt2.AddToStack(hMt2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMt2.SetLegend(0.75,0.6,0.95,0.9); 
    plotMt2.AddTextBox(lumitext,0.55,0.85,0.75,0.79,0);
    plotMt2.AddTextBox(eventtext,0.55,0.75,0.75,0.69,0);
    plotMt2.SetLogy();   
    plotMt2.Draw(c,kFALSE,format,1);
    
    CPlot plotMt2Graph(pname,"",xlabel,"#chi");
    plotMt2Graph.AddHist1D(hMt2Diffv[ich],"E",ratioColor);
    plotMt2Graph.SetYRange(-8,8);
    plotMt2Graph.AddLine(hMt2Diffv[ich]->GetXaxis()->GetXmin(), 0,hMt2Diffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMt2Graph.AddLine(hMt2Diffv[ich]->GetXaxis()->GetXmin(), 5,hMt2Diffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMt2Graph.AddLine(hMt2Diffv[ich]->GetXaxis()->GetXmin(),-5,hMt2Diffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMt2Graph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMtB2v[ich][0]->GetBinWidth(1));
    sprintf(pname,"mtlog_barrel_ch%i",ich);
    CPlot plotMtB2(pname,"",xlabel,ylabel);
    plotMtB2.AddHist1D(hMtB2v[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMtB2.AddToStack(hMtB2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMtB2.SetLegend(0.75,0.6,0.95,0.9); 
    plotMtB2.AddTextBox(lumitext,0.55,0.85,0.75,0.79,0); 
    plotMtB2.AddTextBox(barreltext,0.55,0.75,0.75,0.69,0);
    plotMtB2.AddTextBox(eventtext,0.21,0.85,0.33,0.79,0);
    plotMtB2.SetLogy();
    plotMtB2.Draw(c,kFALSE,format,1);
    
    CPlot plotMtB2Graph(pname,"",xlabel,"#chi");
    plotMtB2Graph.AddHist1D(hMtB2Diffv[ich],"E",ratioColor);
    plotMtB2Graph.SetYRange(-8,8);
    plotMtB2Graph.AddLine(hMtB2Diffv[ich]->GetXaxis()->GetXmin(), 0,hMtB2Diffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMtB2Graph.AddLine(hMtB2Diffv[ich]->GetXaxis()->GetXmin(), 5,hMtB2Diffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMtB2Graph.AddLine(hMtB2Diffv[ich]->GetXaxis()->GetXmin(),-5,hMtB2Diffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMtB2Graph.Draw(c,kTRUE,format,2);
    
    
    sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMtE2v[ich][0]->GetBinWidth(1));
    sprintf(pname,"mtlog_endcap_ch%i",ich);
    CPlot plotMtE2(pname,"",xlabel,ylabel);
    plotMtE2.AddHist1D(hMtE2v[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMtE2.AddToStack(hMtE2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMtE2.SetLegend(0.75,0.6,0.95,0.9); 
    plotMtE2.AddTextBox(lumitext,0.55,0.85,0.75,0.79,0);
    plotMtE2.AddTextBox(endcaptext,0.55,0.75,0.75,0.69,0);
    plotMtE2.AddTextBox(eventtext,0.21,0.85,0.33,0.79,0);
    plotMtE2.SetLogy();
    plotMtE2.Draw(c,kFALSE,format,1);
    
    CPlot plotMtE2Graph(pname,"",xlabel,"#chi");
    plotMtE2Graph.AddHist1D(hMtE2Diffv[ich],"E",ratioColor);
    plotMtE2Graph.SetYRange(-8,8);
    plotMtE2Graph.AddLine(hMtE2Diffv[ich]->GetXaxis()->GetXmin(), 0,hMtE2Diffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMtE2Graph.AddLine(hMtE2Diffv[ich]->GetXaxis()->GetXmin(), 5,hMtE2Diffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMtE2Graph.AddLine(hMtE2Diffv[ich]->GetXaxis()->GetXmin(),-5,hMtE2Diffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotMtE2Graph.Draw(c,kTRUE,format,2);

    //
    // PV multiplicity
    //
    sprintf(xlabel,"N_{PV}");
    if(ich==eInc) sprintf(eventtext,"e^{#pm} events");
    if(ich==ePos) sprintf(eventtext,"e^{+} events");
    if(ich==eNeg) sprintf(eventtext,"e^{-} events");
    
    sprintf(ylabel,"Events");
    sprintf(pname,"npv_ch%i",ich);
    CPlot plotNPV(pname,"",xlabel,ylabel);
    plotNPV.AddHist1D(hNPVv[ich][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotNPV.AddToStack(hNPVv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotNPV.SetLegend(0.75,0.6,0.95,0.9); 
    plotNPV.AddTextBox(lumitext,0.50,0.85,0.70,0.79,0);  
    plotNPV.AddTextBox(eventtext,0.50,0.75,0.70,0.69,0);
    plotNPV.Draw(c,kFALSE,format,1);
    
    CPlot plotNPVGraph(pname,"",xlabel,"#chi");
    plotNPVGraph.AddHist1D(hNPVDiffv[ich],"E",ratioColor);
    plotNPVGraph.SetYRange(-8,8);
    plotNPVGraph.AddLine(hNPVDiffv[ich]->GetXaxis()->GetXmin(), 0,hNPVDiffv[ich]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotNPVGraph.AddLine(hNPVDiffv[ich]->GetXaxis()->GetXmin(), 5,hNPVDiffv[ich]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotNPVGraph.AddLine(hNPVDiffv[ich]->GetXaxis()->GetXmin(),-5,hNPVDiffv[ich]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotNPVGraph.Draw(c,kTRUE,format,2);
  }


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
  makeHTML(outputDir);
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("plotWe"); 
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
TH1D *makeDiffHist(TH1D* hData, TH1D* hMC, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff = (hData->GetBinContent(ibin)-hMC->GetBinContent(ibin));
    
    Double_t err = sqrt(hData->GetBinContent(ibin));
    if(err==0) err= sqrt(hMC->GetBinContent(ibin));
    
    if(err>0) hDiff->SetBinContent(ibin,diff/err);
    else      hDiff->SetBinContent(ibin,0);
    hDiff->SetBinError(ibin,1);   
  }
  
  hDiff->GetYaxis()->SetTitleOffset(0.48);
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
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/plots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Wenu</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
      
  for(UInt_t ich=0; ich<3; ich++) {
    if(ich==0) htmlfile << "<h3 style=\"text-align:left; color:DD6600;\"><a target=\"_blank\" href=\"w_all.html\">W<sup>&plusmn;</sup></a></h3>" << endl;
    if(ich==1) htmlfile << "<h3 style=\"text-align:left; color:DD6600;\"><a target=\"_blank\" href=\"w_pos.html\">W<sup>+</sup></a></h3>" << endl;
    if(ich==2) htmlfile << "<h3 style=\"text-align:left; color:DD6600;\"><a target=\"_blank\" href=\"w_neg.html\">W<sup>-</sup></a></h3>" << endl;
    
    htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lpt_ch" << ich << ".png\"><img src=\"plots/lpt_ch" << ich << ".png\" alt=\"plots/lpt_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/met_ch" << ich << ".png\"><img src=\"plots/met_ch" << ich << ".png\" alt=\"plots/met_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mt_ch" << ich << ".png\"><img src=\"plots/mt_ch" << ich << ".png\" alt=\"plots/mt_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "</table>" << endl;
    htmlfile << "<hr />" << endl;
  }
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
  
  for(UInt_t ich=0; ich<3; ich++) {
    if(ich==0) sprintf(htmlfname,"%s/w_all.html",outDir.Data());
    if(ich==1) sprintf(htmlfname,"%s/w_pos.html",outDir.Data());
    if(ich==2) sprintf(htmlfname,"%s/w_neg.html",outDir.Data());
    
    htmlfile.open(htmlfname);
    htmlfile << "<!DOCTYPE html" << endl;
    htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
    htmlfile << "<html>" << endl;
    htmlfile << "<head><title>Wenu</title></head>" << endl;
    htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  
    htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lpt_ch" << ich << ".png\"><img src=\"plots/lpt_ch" << ich << ".png\" alt=\"plots/lpt_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/met_ch" << ich << ".png\"><img src=\"plots/met_ch" << ich << ".png\" alt=\"plots/met_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mt_ch" << ich << ".png\"><img src=\"plots/mt_ch" << ich << ".png\" alt=\"plots/mt_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lptlog_ch" << ich << ".png\"><img src=\"plots/lptlog_ch" << ich << ".png\" alt=\"plots/lptlog_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metlog_ch" << ich << ".png\"><img src=\"plots/metlog_ch" << ich << ".png\" alt=\"plots/metlog_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mtlog_ch" << ich << ".png\"><img src=\"plots/mtlog_ch" << ich << ".png\" alt=\"plots/mtlog_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/leta_ch" << ich << ".png\"><img src=\"plots/leta_ch" << ich << ".png\" alt=\"plots/leta_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lphi_ch" << ich << ".png\"><img src=\"plots/lphi_ch" << ich << ".png\" alt=\"plots/lphi_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metphi_ch" << ich << ".png\"><img src=\"plots/metphi_ch" << ich << ".png\" alt=\"plots/metphi_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/npv_ch" << ich << ".png\"><img src=\"plots/npv_ch" << ich << ".png\" alt=\"plots/npv_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "</table>" << endl;
    htmlfile << "<hr />" << endl;
    
    htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">Barrel</h3>" << endl;
    htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lpt_barrel_ch" << ich << ".png\"><img src=\"plots/lpt_barrel_ch" << ich << ".png\" alt=\"plots/lpt_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/met_barrel_ch" << ich << ".png\"><img src=\"plots/met_barrel_ch" << ich << ".png\" alt=\"plots/met_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mt_barrel_ch" << ich << ".png\"><img src=\"plots/mt_barrel_ch" << ich << ".png\" alt=\"plots/mt_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lptlog_barrel_ch" << ich << ".png\"><img src=\"plots/lptlog_barrel_ch" << ich << ".png\" alt=\"plots/lptlog_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metlog_barrel_ch" << ich << ".png\"><img src=\"plots/metlog_barrel_ch" << ich << ".png\" alt=\"plots/metlog_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mtlog_barrel_ch" << ich << ".png\"><img src=\"plots/mtlog_barrel_ch" << ich << ".png\" alt=\"plots/mtlog_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lphi_barrel_ch" << ich << ".png\"><img src=\"plots/lphi_barrel_ch" << ich << ".png\" alt=\"plots/lphi_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metphi_barrel_ch" << ich << ".png\"><img src=\"plots/metphi_barrel_ch" << ich << ".png\" alt=\"plots/metphi_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "</table>" << endl;
    htmlfile << "<hr />" << endl;
    
    htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">Endcap</h3>" << endl;
    htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lpt_endcap_ch" << ich << ".png\"><img src=\"plots/lpt_endcap_ch" << ich << ".png\" alt=\"plots/lpt_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/met_endcap_ch" << ich << ".png\"><img src=\"plots/met_endcap_ch" << ich << ".png\" alt=\"plots/met_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mt_endcap_ch" << ich << ".png\"><img src=\"plots/mt_endcap_ch" << ich << ".png\" alt=\"plots/mt_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lptlog_endcap_ch" << ich << ".png\"><img src=\"plots/lptlog_endcap_ch" << ich << ".png\" alt=\"plots/lptlog_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metlog_endcap_ch" << ich << ".png\"><img src=\"plots/metlog_endcap_ch" << ich << ".png\" alt=\"plots/metlog_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mtlog_endcap_ch" << ich << ".png\"><img src=\"plots/mtlog_endcap_ch" << ich << ".png\" alt=\"plots/mtlog_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lphi_endcap_ch" << ich << ".png\"><img src=\"plots/lphi_endcap_ch" << ich << ".png\" alt=\"plots/lphi_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metphi_endcap_ch" << ich << ".png\"><img src=\"plots/metphi_endcap_ch" << ich << ".png\" alt=\"plots/metphi_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "</table>" << endl;
    htmlfile << "<hr />" << endl;
    
    htmlfile << "</body>" << endl;
    htmlfile << "</html>" << endl;
    htmlfile.close();  
  }  
}
