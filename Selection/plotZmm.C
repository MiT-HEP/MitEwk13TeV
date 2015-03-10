//================================================================================================
//
// Make plots of various distributions after Zmumu selection
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

void plotZmm(const TString  conf,            // input file
             const TString  inputDir,        // input directory
	     const TString  outputDir,       // output directory
	     const Double_t lumi             // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotZmm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  const TString format("png");
  
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.1;
  
  const Double_t ETA_BARREL = 1.2;
  const Double_t ETA_ENDCAP = 1.2;
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  enum { eMuMu2HLT=1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk };  // event category enum

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
  vector<TH1D*> hMassv[6], hPt1v[6], hPt2v[6], hyv[6], hPhiv[6], hDPhiv[6], hDRv[6];
  vector<TH1D*> hTagPt1v[6], hTagPt2v[6], hTagEtav[6], hTagPhiv[6];
  vector<TH1D*> hProbePt1v[6], hProbePt2v[6], hProbeEtav[6], hProbePhiv[6];
  vector<TH1D*> hNPVv[6];

  TH1D *hMassMC[6], *hPt1MC[6], *hPt2MC[6], *hyMC[6], *hPhiMC[6], *hDPhiMC[6], *hDRMC[6];
  TH1D *hTagPt1MC[6], *hTagPt2MC[6], *hTagEtaMC[6], *hTagPhiMC[6];
  TH1D *hProbePt1MC[6], *hProbePt2MC[6], *hProbeEtaMC[6], *hProbePhiMC[6];
  TH1D *hNPVMC[6];
  
  char hname[100];
  for(UInt_t icat=0; icat<6; icat++) {
    for(UInt_t isam=0; isam<samplev.size(); isam++) {
      sprintf(hname,"hMass_cat%i_%i",icat,isam);     hMassv[icat].push_back(new TH1D(hname,"",30,60,120));       hMassv[icat][isam]->Sumw2();
      sprintf(hname,"hPt1_cat%i_%i",icat,isam);      hPt1v[icat].push_back(new TH1D(hname,"",30,0,60));          hPt1v[icat][isam]->Sumw2();
      sprintf(hname,"hPt2_cat%i_%i",icat,isam);      hPt2v[icat].push_back(new TH1D(hname,"",40,0,160));         hPt2v[icat][isam]->Sumw2();
      sprintf(hname,"hy_cat%i_%i",icat,isam);        hyv[icat].push_back(new TH1D(hname,"",30,-3,3));            hyv[icat][isam]->Sumw2();
      sprintf(hname,"hPhi_cat%i_%i",icat,isam);      hPhiv[icat].push_back(new TH1D(hname,"",40,-3.2,3.2));      hPhiv[icat][isam]->Sumw2();
      sprintf(hname,"hDPhi_cat%i_%i",icat,isam);     hDPhiv[icat].push_back(new TH1D(hname,"",40,0,3.2));        hDPhiv[icat][isam]->Sumw2();
      sprintf(hname,"hDR_cat%i_%i",icat,isam);       hDRv[icat].push_back(new TH1D(hname,"",30,0,9));            hDRv[icat][isam]->Sumw2();      
      sprintf(hname,"hTagPt1_cat%i_%i",icat,isam);   hTagPt1v[icat].push_back(new TH1D(hname,"",30,20,80));      hTagPt1v[icat][isam]->Sumw2();
      sprintf(hname,"hTagPt2_cat%i_%i",icat,isam);   hTagPt2v[icat].push_back(new TH1D(hname,"",40,20,180));     hTagPt2v[icat][isam]->Sumw2();
      sprintf(hname,"hTagEta_cat%i_%i",icat,isam);   hTagEtav[icat].push_back(new TH1D(hname,"",30,-3,3));       hTagEtav[icat][isam]->Sumw2();
      sprintf(hname,"hTagPhi_cat%i_%i",icat,isam);   hTagPhiv[icat].push_back(new TH1D(hname,"",40,-3.2,3.2));   hTagPhiv[icat][isam]->Sumw2();
      sprintf(hname,"hProbePt1_cat%i_%i",icat,isam); hProbePt1v[icat].push_back(new TH1D(hname,"",30,20,80));    hProbePt1v[icat][isam]->Sumw2();
      sprintf(hname,"hProbePt2_cat%i_%i",icat,isam); hProbePt2v[icat].push_back(new TH1D(hname,"",40,20,180));   hProbePt2v[icat][isam]->Sumw2();
      sprintf(hname,"hProbeEta_cat%i_%i",icat,isam); hProbeEtav[icat].push_back(new TH1D(hname,"",30,-3,3));     hProbeEtav[icat][isam]->Sumw2();
      sprintf(hname,"hProbePhi_cat%i_%i",icat,isam); hProbePhiv[icat].push_back(new TH1D(hname,"",40,-3.2,3.2)); hProbePhiv[icat][isam]->Sumw2();      
      sprintf(hname,"NPV_cat%i_%i",icat,isam);       hNPVv[icat].push_back(new TH1D(hname,"",30,-0.5,29.5));     hNPVv[icat][isam]->Sumw2();
    }
    sprintf(hname,"hMassMC_cat%i",icat);     hMassMC[icat]     = new TH1D(hname,"",30,60,120);    hMassMC[icat]->Sumw2();
    sprintf(hname,"hPt1MC_cat%i",icat);	     hPt1MC[icat]      = new TH1D(hname,"",30,0,60);      hPt1MC[icat]->Sumw2();
    sprintf(hname,"hPt2MC_cat%i",icat);	     hPt2MC[icat]      = new TH1D(hname,"",40,0,160);     hPt2MC[icat]->Sumw2();
    sprintf(hname,"hyMC_cat%i",icat);	     hyMC[icat]        = new TH1D(hname,"",30,-3,3);	  hyMC[icat]->Sumw2();
    sprintf(hname,"hPhiMC_cat%i",icat);	     hPhiMC[icat]      = new TH1D(hname,"",40,-3.2,3.2);  hPhiMC[icat]->Sumw2();
    sprintf(hname,"hDPhiMC_cat%i",icat);     hDPhiMC[icat]     = new TH1D(hname,"",40,0,3.2);     hDPhiMC[icat]->Sumw2();
    sprintf(hname,"hDRMC_cat%i",icat);	     hDRMC[icat]       = new TH1D(hname,"",30,0,9);	  hDRMC[icat]->Sumw2();      
    sprintf(hname,"hTagPt1MC_cat%i",icat);   hTagPt1MC[icat]   = new TH1D(hname,"",30,20,80);     hTagPt1MC[icat]->Sumw2();
    sprintf(hname,"hTagPt2MC_cat%i",icat);   hTagPt2MC[icat]   = new TH1D(hname,"",40,20,180);    hTagPt2MC[icat]->Sumw2();
    sprintf(hname,"hTagEtaMC_cat%i",icat);   hTagEtaMC[icat]   = new TH1D(hname,"",30,-3,3);	  hTagEtaMC[icat]->Sumw2();
    sprintf(hname,"hTagPhiMC_cat%i",icat);   hTagPhiMC[icat]   = new TH1D(hname,"",40,-3.2,3.2);  hTagPhiMC[icat]->Sumw2();
    sprintf(hname,"hProbePt1MC_cat%i",icat); hProbePt1MC[icat] = new TH1D(hname,"",30,20,80);     hProbePt1MC[icat]->Sumw2();
    sprintf(hname,"hProbePt2MC_cat%i",icat); hProbePt2MC[icat] = new TH1D(hname,"",40,20,180);    hProbePt2MC[icat]->Sumw2();
    sprintf(hname,"hProbeEtaMC_cat%i",icat); hProbeEtaMC[icat] = new TH1D(hname,"",30,-3,3);	  hProbeEtaMC[icat]->Sumw2();
    sprintf(hname,"hProbePhiMC_cat%i",icat); hProbePhiMC[icat] = new TH1D(hname,"",40,-3.2,3.2);  hProbePhiMC[icat]->Sumw2();	     
    sprintf(hname,"NPVMC_cat%i",icat);	     hNPVMC[icat]      = new TH1D(hname,"",30,-0.5,29.5); hNPVMC[icat]->Sumw2();
  }
  
  vector<TH1D*> hMassBBv, hMassBEv, hMassEEv;
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    sprintf(hname,"hMassBB_%i",isam); hMassBBv.push_back(new TH1D(hname,"",30,60,120)); hMassBBv[isam]->Sumw2();
    sprintf(hname,"hMassBE_%i",isam); hMassBEv.push_back(new TH1D(hname,"",30,60,120)); hMassBEv[isam]->Sumw2();
    sprintf(hname,"hMassEE_%i",isam); hMassEEv.push_back(new TH1D(hname,"",30,60,120)); hMassEEv[isam]->Sumw2();
  }
  TH1D* hMassBBMC = new TH1D("hMassBBMC","",30,60,120); hMassBBMC->Sumw2();
  TH1D* hMassBEMC = new TH1D("hMassBEMC","",30,60,120); hMassBEMC->Sumw2();
  TH1D* hMassEEMC = new TH1D("hMassEEMC","",30,60,120); hMassEEMC->Sumw2();
  
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
  ///// muon specific /////
  Float_t trkIso1, emIso1, hadIso1, trkIso2, emIso2, hadIso2;
  Float_t pfChIso1, pfGamIso1, pfNeuIso1, pfCombIso1, pfChIso2, pfGamIso2, pfNeuIso2, pfCombIso2;
  Float_t d01, dz1, d02, dz2;
  Float_t muNchi21,  muNchi22;
  UInt_t nPixHits1, nTkLayers1, nPixHits2, nTkLayers2;
  UInt_t nValidHits1, nMatch1, nValidHits2, nMatch2;
  UInt_t typeBits1, typeBits2;
  LorentzVector *sta1=0, *sta2=0;

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

    intree->SetBranchAddress("runNum",      &runNum);        // event run number
    intree->SetBranchAddress("lumiSec",     &lumiSec);       // event lumi section
    intree->SetBranchAddress("evtNum",      &evtNum);        // event number
    intree->SetBranchAddress("matchGen",    &matchGen);      // event has both leptons matched to MC Z->ll
    intree->SetBranchAddress("category",    &category);      // dilepton category
    intree->SetBranchAddress("npv",         &npv);	     // number of primary vertices
    intree->SetBranchAddress("npu",         &npu);	     // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",      &genVPt);        // GEN boson pT
    intree->SetBranchAddress("genVPhi",     &genVPhi);       // GEN boson phi
    intree->SetBranchAddress("genVy",       &genVy);         // GEN boson rapidity
    intree->SetBranchAddress("genVMass",    &genVMass);      // GEN boson mass
    intree->SetBranchAddress("scale1fb",    &scale1fb);      // event weight per 1/fb (MC)
    intree->SetBranchAddress("met",         &met);           // MET
    intree->SetBranchAddress("metPhi",      &metPhi);        // phi(MET)
    intree->SetBranchAddress("sumEt",       &sumEt);         // Sum ET
    intree->SetBranchAddress("u1",          &u1);	     // parallel component of recoil
    intree->SetBranchAddress("u2",          &u2);	     // perpendicular component of recoil
    intree->SetBranchAddress("q1",          &q1);	     // charge of tag lepton
    intree->SetBranchAddress("q2",          &q2);            // charge of probe lepton
    intree->SetBranchAddress("dilep",       &dilep);         // dilepton 4-vector
    intree->SetBranchAddress("lep1",        &lep1);          // tag lepton 4-vector
    intree->SetBranchAddress("lep2",        &lep2);          // probe lepton 4-vector    
    ///// muon specific /////
    intree->SetBranchAddress("trkIso1",     &trkIso1);       // track isolation of tag lepton
    intree->SetBranchAddress("trkIso2",     &trkIso2);       // track isolation of probe lepton
    intree->SetBranchAddress("emIso1",      &emIso1);	     // ECAL isolation of tag lepton
    intree->SetBranchAddress("emIso2",      &emIso2);	     // ECAL isolation of probe lepton
    intree->SetBranchAddress("hadIso1",     &hadIso1);       // HCAL isolation of tag lepton
    intree->SetBranchAddress("hadIso2",     &hadIso2);       // HCAL isolation of probe lepton
    intree->SetBranchAddress("pfChIso1",    &pfChIso1);      // PF charged hadron isolation of tag lepton
    intree->SetBranchAddress("pfChIso2",    &pfChIso2);      // PF charged hadron isolation of probe lepton
    intree->SetBranchAddress("pfGamIso1",   &pfGamIso1);     // PF photon isolation of tag lepton
    intree->SetBranchAddress("pfGamIso2",   &pfGamIso2);     // PF photon isolation of probe lepton
    intree->SetBranchAddress("pfNeuIso1",   &pfNeuIso1);     // PF neutral hadron isolation of tag lepton
    intree->SetBranchAddress("pfNeuIso2",   &pfNeuIso2);     // PF neutral hadron isolation of probe lepton
    intree->SetBranchAddress("pfCombIso1",  &pfCombIso1);    // PF combined isolation of tag lepton
    intree->SetBranchAddress("pfCombIso2",  &pfCombIso2);    // PF combined isolation of probe lepton    
    intree->SetBranchAddress("d01",         &d01);	     // transverse impact parameter of tag lepton
    intree->SetBranchAddress("d02",         &d02);	     // transverse impact parameter of probe lepton	 
    intree->SetBranchAddress("dz1",         &dz1);	     // longitudinal impact parameter of tag lepton
    intree->SetBranchAddress("dz2",         &dz2);	     // longitudinal impact parameter of probe lepton	 
    intree->SetBranchAddress("muNchi21",    &muNchi21);      // muon fit normalized chi^2 of tag lepton
    intree->SetBranchAddress("muNchi22",    &muNchi22);      // muon fit normalized chi^2 of probe lepton	     
    intree->SetBranchAddress("nPixHits1",   &nPixHits1);     // number of pixel hits of tag muon
    intree->SetBranchAddress("nPixHits2",   &nPixHits2);     // number of pixel hits of probe muon
    intree->SetBranchAddress("nTkLayers1",  &nTkLayers1);    // number of tracker layers of tag muon
    intree->SetBranchAddress("nTkLayers2",  &nTkLayers2);    // number of tracker layers of probe muon
    intree->SetBranchAddress("nMatch1",     &nMatch1);       // number of matched segments of tag muon
    intree->SetBranchAddress("nMatch2",     &nMatch2);       // number of matched segments of probe muon    
    intree->SetBranchAddress("nValidHits1", &nValidHits1);   // number of valid muon hits of tag muon
    intree->SetBranchAddress("nValidHits2", &nValidHits2);   // number of valid muon hits of probe muon
    intree->SetBranchAddress("typeBits1",   &typeBits1);     // muon type of tag muon
    intree->SetBranchAddress("typeBits2",   &typeBits2);     // muon type of probe muon
    intree->SetBranchAddress("sta1",        &sta1);	     // tag STA muon 4-vector
    intree->SetBranchAddress("sta2",        &sta2);	     // probe STA muon 4-vector   
    
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);

      if(dilep->M()        < MASS_LOW)  continue;
      if(dilep->M()        > MASS_HIGH) continue;
      if(lep1->Pt()        < PT_CUT)    continue;
      if(lep2->Pt()        < PT_CUT)    continue;
      if(fabs(lep1->Eta()) > ETA_CUT)   continue;      
      if(fabs(lep2->Eta()) > ETA_CUT)   continue;
      
      Double_t weight = 1;
      if(isam!=0) {
        weight *= scale1fb*lumi;
      }
      
      hMassv[category][isam]    ->Fill(dilep->M(),weight);
      hPt1v[category][isam]     ->Fill(dilep->Pt(),weight);
      hPt2v[category][isam]     ->Fill(dilep->Pt(),weight);
      hyv[category][isam]       ->Fill(dilep->Rapidity(),weight);
      hPhiv[category][isam]     ->Fill(dilep->Phi(),weight);
      hDPhiv[category][isam]    ->Fill(toolbox::deltaPhi(lep1->Phi(), lep2->Phi()),weight);
      hDRv[category][isam]      ->Fill(toolbox::deltaR(lep1->Eta(),lep1->Phi(),lep2->Eta(),lep2->Phi()),weight);	
      hTagPt1v[category][isam]  ->Fill(lep1->Pt(),weight);
      hTagPt2v[category][isam]  ->Fill(lep1->Pt(),weight);
      hTagEtav[category][isam]  ->Fill(lep1->Eta(),weight);
      hTagPhiv[category][isam]  ->Fill(lep1->Phi(),weight);
      hProbePt1v[category][isam]->Fill(lep2->Pt(),weight);
      hProbePt2v[category][isam]->Fill(lep2->Pt(),weight);
      hProbeEtav[category][isam]->Fill(lep2->Eta(),weight);
      hProbePhiv[category][isam]->Fill(lep2->Phi(),weight);
      hNPVv[category][isam]     ->Fill(npv,weight);
      if(isam!=0) {
        hMassMC[category]    ->Fill(dilep->M(),weight);
        hPt1MC[category]     ->Fill(dilep->Pt(),weight);
        hPt2MC[category]     ->Fill(dilep->Pt(),weight);
        hyMC[category]       ->Fill(dilep->Rapidity(),weight);
        hPhiMC[category]     ->Fill(dilep->Phi(),weight);
        hDPhiMC[category]    ->Fill(toolbox::deltaPhi(lep1->Phi(), lep2->Phi()),weight);
        hDRMC[category]      ->Fill(toolbox::deltaR(lep1->Eta(),lep1->Phi(),lep2->Eta(),lep2->Phi()),weight);	
        hTagPt1MC[category]  ->Fill(lep1->Pt(),weight);
        hTagPt2MC[category]  ->Fill(lep1->Pt(),weight);
        hTagEtaMC[category]  ->Fill(lep1->Eta(),weight);
        hTagPhiMC[category]  ->Fill(lep1->Phi(),weight);
        hProbePt1MC[category]->Fill(lep2->Pt(),weight);
        hProbePt2MC[category]->Fill(lep2->Pt(),weight);
        hProbeEtaMC[category]->Fill(lep2->Eta(),weight);
        hProbePhiMC[category]->Fill(lep2->Phi(),weight);
        hNPVMC[category]     ->Fill(npv,weight);
      }

      if(category==eMuMu2HLT || category==eMuMu1HLT) {
        hMassv[0][isam]    ->Fill(dilep->M(),weight);
        hPt1v[0][isam]     ->Fill(dilep->Pt(),weight);
        hPt2v[0][isam]     ->Fill(dilep->Pt(),weight);
        hyv[0][isam]       ->Fill(dilep->Rapidity(),weight);
        hPhiv[0][isam]     ->Fill(dilep->Phi(),weight);
        hDPhiv[0][isam]    ->Fill(toolbox::deltaPhi(lep1->Phi(), lep2->Phi()),weight);
        hDRv[0][isam]      ->Fill(toolbox::deltaR(lep1->Eta(),lep1->Phi(),lep2->Eta(),lep2->Phi()),weight);	
        hTagPt1v[0][isam]  ->Fill((lep1->Pt()>lep2->Pt()) ? lep1->Pt()  : lep2->Pt(),weight);
        hTagPt2v[0][isam]  ->Fill((lep1->Pt()>lep2->Pt()) ? lep1->Pt()  : lep2->Pt(),weight);
        hTagEtav[0][isam]  ->Fill((lep1->Pt()>lep2->Pt()) ? lep1->Eta() : lep2->Eta(),weight);
        hTagPhiv[0][isam]  ->Fill((lep1->Pt()>lep2->Pt()) ? lep1->Phi() : lep2->Phi(),weight);
        hProbePt1v[0][isam]->Fill((lep1->Pt()>lep2->Pt()) ? lep2->Pt()  : lep1->Pt(),weight);
        hProbePt2v[0][isam]->Fill((lep1->Pt()>lep2->Pt()) ? lep2->Pt()  : lep1->Pt(),weight);
        hProbeEtav[0][isam]->Fill((lep1->Pt()>lep2->Pt()) ? lep2->Eta() : lep1->Eta(),weight);
        hProbePhiv[0][isam]->Fill((lep1->Pt()>lep2->Pt()) ? lep2->Phi() : lep1->Phi(),weight);
        hNPVv[0][isam]     ->Fill(npv,weight);
	
	if(fabs(lep1->Eta()) < ETA_BARREL && fabs(lep2->Eta()) < ETA_BARREL) {
	  hMassBBv[isam]->Fill(dilep->M(),weight);
	} else if(fabs(lep1->Eta()) > ETA_ENDCAP && fabs(lep2->Eta()) > ETA_ENDCAP) {
	  hMassEEv[isam]->Fill(dilep->M(),weight);
	} else {
	  hMassBEv[isam]->Fill(dilep->M(),weight);
	}
	
	if(isam!=0) {
	  hMassMC[0]    ->Fill(dilep->M(),weight);
          hPt1MC[0]     ->Fill(dilep->Pt(),weight);
          hPt2MC[0]     ->Fill(dilep->Pt(),weight);
          hyMC[0]       ->Fill(dilep->Rapidity(),weight);
          hPhiMC[0]     ->Fill(dilep->Phi(),weight);
          hDPhiMC[0]    ->Fill(toolbox::deltaPhi(lep1->Phi(), lep2->Phi()),weight);
          hDRMC[0]      ->Fill(toolbox::deltaR(lep1->Eta(),lep1->Phi(),lep2->Eta(),lep2->Phi()),weight);	
          hTagPt1MC[0]  ->Fill((lep1->Pt()>lep2->Pt()) ? lep1->Pt()  : lep2->Pt(),weight);
          hTagPt2MC[0]  ->Fill((lep1->Pt()>lep2->Pt()) ? lep1->Pt()  : lep2->Pt(),weight);
          hTagEtaMC[0]  ->Fill((lep1->Pt()>lep2->Pt()) ? lep1->Eta() : lep2->Eta(),weight);
          hTagPhiMC[0]  ->Fill((lep1->Pt()>lep2->Pt()) ? lep1->Phi() : lep2->Phi(),weight);
          hProbePt1MC[0]->Fill((lep1->Pt()>lep2->Pt()) ? lep2->Pt()  : lep1->Pt(),weight);
          hProbePt2MC[0]->Fill((lep1->Pt()>lep2->Pt()) ? lep2->Pt()  : lep1->Pt(),weight);
          hProbeEtaMC[0]->Fill((lep1->Pt()>lep2->Pt()) ? lep2->Eta() : lep1->Eta(),weight);
          hProbePhiMC[0]->Fill((lep1->Pt()>lep2->Pt()) ? lep2->Phi() : lep1->Phi(),weight);
          hNPVMC[0]     ->Fill(npv,weight);
	  
	  if(fabs(lep1->Eta()) < ETA_BARREL && fabs(lep2->Eta()) < ETA_BARREL) {
	    hMassBBMC->Fill(dilep->M(),weight);
	  } else if(fabs(lep1->Eta()) > ETA_ENDCAP && fabs(lep2->Eta()) > ETA_ENDCAP) {
	    hMassEEMC->Fill(dilep->M(),weight);
	  } else {
	    hMassBEMC->Fill(dilep->M(),weight);
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
  vector<TH1D*> hMassDiffv, hPt1Diffv, hPt2Diffv, hyDiffv, hPhiDiffv, hDPhiDiffv, hDRDiffv;
  vector<TH1D*> hTagPt1Diffv, hTagPt2Diffv, hTagEtaDiffv, hTagPhiDiffv;
  vector<TH1D*> hProbePt1Diffv, hProbePt2Diffv, hProbeEtaDiffv, hProbePhiDiffv;
  vector<TH1D*> hNPVDiffv;
  TH1D* hMassBBDiff = makeDiffHist(hMassBBv[0],hMassBBMC,"hMassBBDiff");
  TH1D* hMassBEDiff = makeDiffHist(hMassBEv[0],hMassBEMC,"hMassBEDiff");
  TH1D* hMassEEDiff = makeDiffHist(hMassEEv[0],hMassEEMC,"hMassEEDiff");
  for(Int_t icat=0; icat<6; icat++) {
    sprintf(hname,"hMassDiff_cat%i",icat);      hMassDiffv.push_back(makeDiffHist(hMassv[icat][0],hMassMC[icat],hname));
    sprintf(hname,"hPt1Diff_cat%i",icat);       hPt1Diffv.push_back(makeDiffHist(hPt1v[icat][0],hPt1MC[icat],hname));
    sprintf(hname,"hPt2Diff_cat%i",icat);       hPt2Diffv.push_back(makeDiffHist(hPt2v[icat][0],hPt2MC[icat],hname));
    sprintf(hname,"hyDiff_cat%i",icat);         hyDiffv.push_back(makeDiffHist(hyv[icat][0],hyMC[icat],hname));
    sprintf(hname,"hPhiDiff_cat%i",icat);       hPhiDiffv.push_back(makeDiffHist(hPhiv[icat][0],hPhiMC[icat],hname));
    sprintf(hname,"hDPhiDiff_cat%i",icat);      hDPhiDiffv.push_back(makeDiffHist(hDPhiv[icat][0],hDPhiMC[icat],hname));
    sprintf(hname,"hDRDiff_cat%i",icat);        hDRDiffv.push_back(makeDiffHist(hDRv[icat][0],hDRMC[icat],hname));
    sprintf(hname,"hTagPt1Diff_cat%i",icat);    hTagPt1Diffv.push_back(makeDiffHist(hTagPt1v[icat][0],hTagPt1MC[icat],hname));
    sprintf(hname,"hTagPt2Diff_cat%i",icat);    hTagPt2Diffv.push_back(makeDiffHist(hTagPt2v[icat][0],hTagPt2MC[icat],hname));
    sprintf(hname,"hTagEtaDiff_cat%i",icat);    hTagEtaDiffv.push_back(makeDiffHist(hTagEtav[icat][0],hTagEtaMC[icat],hname));
    sprintf(hname,"hTagPhiDiff_cat%i",icat);    hTagPhiDiffv.push_back(makeDiffHist(hTagPhiv[icat][0],hTagPhiMC[icat],hname));
    sprintf(hname,"hProbePt1Diff_cat%i",icat);  hProbePt1Diffv.push_back(makeDiffHist(hProbePt1v[icat][0],hProbePt1MC[icat],hname));
    sprintf(hname,"hProbePt2Diff_cat%i",icat);  hProbePt2Diffv.push_back(makeDiffHist(hProbePt2v[icat][0],hProbePt2MC[icat],hname));
    sprintf(hname,"hProbeEtaDiff_cat%i",icat);  hProbeEtaDiffv.push_back(makeDiffHist(hProbeEtav[icat][0],hProbeEtaMC[icat],hname));
    sprintf(hname,"hProbePhiDiff_cat%i",icat);  hProbePhiDiffv.push_back(makeDiffHist(hProbePhiv[icat][0],hProbePhiMC[icat],hname));
    sprintf(hname,"hNPVDiff_cat%i",icat);       hNPVDiffv.push_back(makeDiffHist(hNPVv[icat][0],hNPVMC[icat],hname));
  }

  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  // string buffers
  char pname[100];       // plot name
  char xlabel[100];      // x-axis label
  char ylabel[100];      // y-axis label
  char lumitext[100];    // lumi label

  if(lumi<0.01)     sprintf(lumitext,"%.2f pb^{-1}  at  #sqrt{s} = 8 TeV",1000.*lumi);
  else if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 8 TeV",1000.*lumi);
  else              sprintf(lumitext,"%.0f fb^{-1}  at  #sqrt{s} = 8 TeV",lumi);
    
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

  for(UInt_t icat=0; icat<6; icat++) {
    
    //
    // Z mass
    //
    sprintf(xlabel,"m(#mu^{+}#mu^{-}) [GeV/c^{2}]");
    
    sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassv[icat][0]->GetBinWidth(1));
    sprintf(pname,"zmass_cat%i",icat);
    CPlot plotMass(pname,"",xlabel,ylabel);
    plotMass.AddHist1D(hMassv[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotMass.AddToStack(hMassv[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotMass.SetLegend(0.70,0.55,0.90,0.85); 
    plotMass.AddTextBox(lumitext,0.21,0.79,0.51,0.85,0);
    plotMass.Draw(c,kFALSE,format,1);
    
    CPlot plotMassGraph(pname,"",xlabel,"#chi");
    plotMassGraph.AddHist1D(hMassDiffv[icat],"E",ratioColor);
    plotMassGraph.SetYRange(-8,8);
    plotMassGraph.AddLine(hMassDiffv[icat]->GetXaxis()->GetXmin(), 0,hMassDiffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotMassGraph.AddLine(hMassDiffv[icat]->GetXaxis()->GetXmin(), 5,hMassDiffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotMassGraph.AddLine(hMassDiffv[icat]->GetXaxis()->GetXmin(),-5,hMassDiffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotMassGraph.Draw(c,kTRUE,format,2);

    sprintf(pname,"zmasslog_cat%i",icat);
    plotMass.SetName(pname);
    plotMass.SetLogy();   
    plotMass.Draw(c,kTRUE,format,1);
    
    //
    // Z pT
    //
    sprintf(xlabel,"p_{T}(#mu^{+}#mu^{-}) [GeV/c]");
    
    sprintf(ylabel,"Events / %.1f GeV/c",hPt1v[icat][0]->GetBinWidth(1));
    sprintf(pname,"zpt_cat%i",icat);
    CPlot plotPt1(pname,"",xlabel,ylabel);
    plotPt1.AddHist1D(hPt1v[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPt1.AddToStack(hPt1v[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPt1.SetLegend(0.70,0.55,0.90,0.85); 
    plotPt1.AddTextBox(lumitext,0.35,0.85,0.65,0.79,0);
    plotPt1.Draw(c,kFALSE,format,1);
    
    CPlot plotPt1Graph(pname,"",xlabel,"#chi");
    plotPt1Graph.AddHist1D(hPt1Diffv[icat],"E",ratioColor);
    plotPt1Graph.SetYRange(-8,8);
    plotPt1Graph.AddLine(hPt1Diffv[icat]->GetXaxis()->GetXmin(), 0,hPt1Diffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPt1Graph.AddLine(hPt1Diffv[icat]->GetXaxis()->GetXmin(), 5,hPt1Diffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPt1Graph.AddLine(hPt1Diffv[icat]->GetXaxis()->GetXmin(),-5,hPt1Diffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotPt1Graph.Draw(c,kTRUE,format,2);
    
    sprintf(ylabel,"Events / %.1f GeV/c",hPt2v[icat][0]->GetBinWidth(1));
    sprintf(pname,"zptlog_cat%i",icat);
    CPlot plotPt2(pname,"",xlabel,ylabel);
    plotPt2.AddHist1D(hPt2v[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPt2.AddToStack(hPt2v[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPt2.SetLegend(0.70,0.55,0.90,0.85); 
    plotPt2.AddTextBox(lumitext,0.35,0.85,0.65,0.79,0);
    plotPt2.SetLogy();   
    plotPt2.Draw(c,kFALSE,format,1);
    
    CPlot plotPt2Graph(pname,"",xlabel,"#chi");
    plotPt2Graph.AddHist1D(hPt2Diffv[icat],"E",ratioColor);
    plotPt2Graph.SetYRange(-8,8);
    plotPt2Graph.AddLine(hPt2Diffv[icat]->GetXaxis()->GetXmin(), 0,hPt2Diffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPt2Graph.AddLine(hPt2Diffv[icat]->GetXaxis()->GetXmin(), 5,hPt2Diffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPt2Graph.AddLine(hPt2Diffv[icat]->GetXaxis()->GetXmin(),-5,hPt2Diffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotPt2Graph.Draw(c,kTRUE,format,2);
    
    //
    // Z rapidity
    //
    sprintf(xlabel,"y(#mu^{+}#mu^{-}");
    
    sprintf(ylabel,"Events / %.2f",hyv[icat][0]->GetBinWidth(1));
    sprintf(pname,"zy_cat%i",icat);
    CPlot ploty(pname,"",xlabel,ylabel);
    ploty.AddHist1D(hyv[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      ploty.AddToStack(hyv[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    ploty.SetLegend(0.70,0.55,0.90,0.85); 
    ploty.AddTextBox(lumitext,0.21,0.79,0.51,0.85,0);
    ploty.SetYRange(0,1.8*(hyv[icat][0]->GetMaximum()));
    ploty.Draw(c,kFALSE,format,1);
    
    CPlot plotyGraph(pname,"",xlabel,"#chi");
    plotyGraph.AddHist1D(hyDiffv[icat],"E",ratioColor);
    plotyGraph.SetYRange(-8,8);
    plotyGraph.AddLine(hyDiffv[icat]->GetXaxis()->GetXmin(), 0,hyDiffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotyGraph.AddLine(hyDiffv[icat]->GetXaxis()->GetXmin(), 5,hyDiffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotyGraph.AddLine(hyDiffv[icat]->GetXaxis()->GetXmin(),-5,hyDiffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotyGraph.Draw(c,kTRUE,format,2);
    
    //
    // Z phi
    //
    sprintf(xlabel,"#phi(#mu^{+}#mu^{-})");
    
    sprintf(ylabel,"Events / %.2f",hPhiv[icat][0]->GetBinWidth(1));
    sprintf(pname,"zphi_cat%i",icat);
    CPlot plotPhi(pname,"",xlabel,ylabel);
    plotPhi.AddHist1D(hPhiv[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotPhi.AddToStack(hPhiv[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotPhi.SetLegend(0.70,0.55,0.90,0.85); 
    plotPhi.AddTextBox(lumitext,0.21,0.79,0.51,0.85,0);
    plotPhi.SetYRange(0,1.8*(hPhiv[icat][0]->GetMaximum()));
    plotPhi.Draw(c,kFALSE,format,1);
    
    CPlot plotPhiGraph(pname,"",xlabel,"#chi");
    plotPhiGraph.AddHist1D(hPhiDiffv[icat],"E",ratioColor);
    plotPhiGraph.SetYRange(-8,8);
    plotPhiGraph.AddLine(hPhiDiffv[icat]->GetXaxis()->GetXmin(), 0,hPhiDiffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotPhiGraph.AddLine(hPhiDiffv[icat]->GetXaxis()->GetXmin(), 5,hPhiDiffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotPhiGraph.AddLine(hPhiDiffv[icat]->GetXaxis()->GetXmin(),-5,hPhiDiffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotPhiGraph.Draw(c,kTRUE,format,2);    
    
    //
    // dilepton dphi
    //
    sprintf(xlabel,"#Delta^{}#phi^{}(e^{+},e^{-})");
    
    sprintf(ylabel,"Events / %.2f",hDPhiv[icat][0]->GetBinWidth(1));
    sprintf(pname,"zdphi_cat%i",icat);
    CPlot plotDPhi(pname,"",xlabel,ylabel);
    plotDPhi.AddHist1D(hDPhiv[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotDPhi.AddToStack(hDPhiv[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotDPhi.SetLegend(0.23,0.55,0.43,0.85); 
    plotDPhi.AddTextBox(lumitext,0.45,0.85,0.75,0.79,0);
    plotDPhi.Draw(c,kFALSE,format,1);
    
    CPlot plotDPhiGraph(pname,"",xlabel,"#chi");
    plotDPhiGraph.AddHist1D(hDPhiDiffv[icat],"E",ratioColor);
    plotDPhiGraph.SetYRange(-8,8);
    plotDPhiGraph.AddLine(hDPhiDiffv[icat]->GetXaxis()->GetXmin(), 0,hDPhiDiffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotDPhiGraph.AddLine(hDPhiDiffv[icat]->GetXaxis()->GetXmin(), 5,hDPhiDiffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotDPhiGraph.AddLine(hDPhiDiffv[icat]->GetXaxis()->GetXmin(),-5,hDPhiDiffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotDPhiGraph.Draw(c,kTRUE,format,2);    
    
    //
    // dilepton dR
    //
    sprintf(xlabel,"#Delta^{}(e^{+},e^{-})");
    
    sprintf(ylabel,"Events / %.2f",hDRv[icat][0]->GetBinWidth(1));
    sprintf(pname,"zdr_cat%i",icat);
    CPlot plotDR(pname,"",xlabel,ylabel);
    plotDR.AddHist1D(hDRv[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotDR.AddToStack(hDRv[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotDR.SetLegend(0.70,0.48,0.90,0.78); 
    plotDR.AddTextBox(lumitext,0.60,0.85,0.90,0.79,0);
    plotDR.Draw(c,kFALSE,format,1);
    
    CPlot plotDRGraph(pname,"",xlabel,"#chi");
    plotDRGraph.AddHist1D(hDRDiffv[icat],"E",ratioColor);
    plotDRGraph.SetYRange(-8,8);
    plotDRGraph.AddLine(hDRDiffv[icat]->GetXaxis()->GetXmin(), 0,hDRDiffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotDRGraph.AddLine(hDRDiffv[icat]->GetXaxis()->GetXmin(), 5,hDRDiffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotDRGraph.AddLine(hDRDiffv[icat]->GetXaxis()->GetXmin(),-5,hDRDiffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotDRGraph.Draw(c,kTRUE,format,2);    
    
    //
    // Tag lepton pT
    //
    if(icat==0) sprintf(xlabel,"leading muon p_{T} [GeV/c]");
    else        sprintf(xlabel,"tag muon p_{T} [GeV/c]");
    
    sprintf(ylabel,"Events / %.1f GeV/c",hTagPt1v[icat][0]->GetBinWidth(1));
    sprintf(pname,"tagpt_cat%i",icat);
    CPlot plotTagPt1(pname,"",xlabel,ylabel);
    plotTagPt1.AddHist1D(hTagPt1v[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotTagPt1.AddToStack(hTagPt1v[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotTagPt1.SetLegend(0.70,0.48,0.90,0.78); 
    plotTagPt1.AddTextBox(lumitext,0.60,0.85,0.90,0.79,0);
    plotTagPt1.Draw(c,kFALSE,format,1);
    
    CPlot plotTagPt1Graph(pname,"",xlabel,"#chi");
    plotTagPt1Graph.AddHist1D(hTagPt1Diffv[icat],"E",ratioColor);
    plotTagPt1Graph.SetYRange(-8,8);
    plotTagPt1Graph.AddLine(hTagPt1Diffv[icat]->GetXaxis()->GetXmin(), 0,hTagPt1Diffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotTagPt1Graph.AddLine(hTagPt1Diffv[icat]->GetXaxis()->GetXmin(), 5,hTagPt1Diffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotTagPt1Graph.AddLine(hTagPt1Diffv[icat]->GetXaxis()->GetXmin(),-5,hTagPt1Diffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotTagPt1Graph.Draw(c,kTRUE,format,2);

    
    sprintf(ylabel,"Events / %.1f GeV/c",hTagPt2v[icat][0]->GetBinWidth(1));
    sprintf(pname,"tagptlog_cat%i",icat);
    CPlot plotTagPt2(pname,"",xlabel,ylabel);
    plotTagPt2.AddHist1D(hTagPt2v[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotTagPt2.AddToStack(hTagPt2v[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotTagPt2.SetLegend(0.70,0.48,0.90,0.78); 
    plotTagPt2.AddTextBox(lumitext,0.60,0.85,0.90,0.79,0);
    plotTagPt2.SetLogy();   
    plotTagPt2.Draw(c,kFALSE,format,1);
    
    CPlot plotTagPt2Graph(pname,"",xlabel,"#chi");
    plotTagPt2Graph.AddHist1D(hTagPt2Diffv[icat],"E",ratioColor);
    plotTagPt2Graph.SetYRange(-8,8);
    plotTagPt2Graph.AddLine(hTagPt2Diffv[icat]->GetXaxis()->GetXmin(), 0,hTagPt2Diffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotTagPt2Graph.AddLine(hTagPt2Diffv[icat]->GetXaxis()->GetXmin(), 5,hTagPt2Diffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotTagPt2Graph.AddLine(hTagPt2Diffv[icat]->GetXaxis()->GetXmin(),-5,hTagPt2Diffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotTagPt2Graph.Draw(c,kTRUE,format,2);
    
    //
    // Tag lepton eta
    //
    if(icat==0) sprintf(xlabel,"leading muon #eta");
    else        sprintf(xlabel,"tag muon #eta");
    
    sprintf(ylabel,"Events / %.2f",hTagEtav[icat][0]->GetBinWidth(1));
    sprintf(pname,"tageta_cat%i",icat);
    CPlot plotTagEta(pname,"",xlabel,ylabel);
    plotTagEta.AddHist1D(hTagEtav[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotTagEta.AddToStack(hTagEtav[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotTagEta.SetLegend(0.70,0.55,0.90,0.85); 
    plotTagEta.AddTextBox(lumitext,0.21,0.85,0.51,0.79,0);
    plotTagEta.SetYRange(0,1.8*(hTagEtav[icat][0]->GetMaximum()));
    plotTagEta.Draw(c,kFALSE,format,1);
    
    CPlot plotTagEtaGraph(pname,"",xlabel,"#chi");
    plotTagEtaGraph.AddHist1D(hTagEtaDiffv[icat],"E",ratioColor);
    plotTagEtaGraph.SetYRange(-8,8);
    plotTagEtaGraph.AddLine(hTagEtaDiffv[icat]->GetXaxis()->GetXmin(), 0,hTagEtaDiffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotTagEtaGraph.AddLine(hTagEtaDiffv[icat]->GetXaxis()->GetXmin(), 5,hTagEtaDiffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotTagEtaGraph.AddLine(hTagEtaDiffv[icat]->GetXaxis()->GetXmin(),-5,hTagEtaDiffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotTagEtaGraph.Draw(c,kTRUE,format,2);
    
    //
    // Tag lepton phi
    //
    if(icat==0) sprintf(xlabel,"leading muon #phi");
    else        sprintf(xlabel,"tag muon #phi");
    
    sprintf(ylabel,"Events / %.2f",hTagPhiv[icat][0]->GetBinWidth(1));
    sprintf(pname,"tagphi_cat%i",icat);
    CPlot plotTagPhi(pname,"",xlabel,ylabel);
    plotTagPhi.AddHist1D(hTagPhiv[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotTagPhi.AddToStack(hTagPhiv[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotTagPhi.SetLegend(0.70,0.55,0.90,0.85); 
    plotTagPhi.AddTextBox(lumitext,0.21,0.85,0.51,0.79,0);
    plotTagPhi.SetYRange(0,1.8*(hTagPhiv[icat][0]->GetMaximum()));
    plotTagPhi.Draw(c,kFALSE,format,1);
    
    CPlot plotTagPhiGraph(pname,"",xlabel,"#chi");
    plotTagPhiGraph.AddHist1D(hTagPhiDiffv[icat],"E",ratioColor);
    plotTagPhiGraph.SetYRange(-8,8);
    plotTagPhiGraph.AddLine(hTagPhiDiffv[icat]->GetXaxis()->GetXmin(), 0,hTagPhiDiffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotTagPhiGraph.AddLine(hTagPhiDiffv[icat]->GetXaxis()->GetXmin(), 5,hTagPhiDiffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotTagPhiGraph.AddLine(hTagPhiDiffv[icat]->GetXaxis()->GetXmin(),-5,hTagPhiDiffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotTagPhiGraph.Draw(c,kTRUE,format,2);
    
    //
    // Probe lepton pT
    //
    if(icat==0) sprintf(xlabel,"trailing muon p_{T} [GeV/c]");
    else        sprintf(xlabel,"probe muon p_{T} [GeV/c]");
    
    sprintf(ylabel,"Events / %.1f GeV/c",hProbePt1v[icat][0]->GetBinWidth(1));
    sprintf(pname,"probept_cat%i",icat);
    CPlot plotProbePt1(pname,"",xlabel,ylabel);
    plotProbePt1.AddHist1D(hProbePt1v[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotProbePt1.AddToStack(hProbePt1v[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotProbePt1.SetLegend(0.70,0.48,0.90,0.78); 
    plotProbePt1.AddTextBox(lumitext,0.60,0.85,0.90,0.79,0);
    plotProbePt1.Draw(c,kFALSE,format,1);
    
    CPlot plotProbePt1Graph(pname,"",xlabel,"#chi");
    plotProbePt1Graph.AddHist1D(hProbePt1Diffv[icat],"E",ratioColor);
    plotProbePt1Graph.SetYRange(-8,8);
    plotProbePt1Graph.AddLine(hProbePt1Diffv[icat]->GetXaxis()->GetXmin(), 0,hProbePt1Diffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotProbePt1Graph.AddLine(hProbePt1Diffv[icat]->GetXaxis()->GetXmin(), 5,hProbePt1Diffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotProbePt1Graph.AddLine(hProbePt1Diffv[icat]->GetXaxis()->GetXmin(),-5,hProbePt1Diffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotProbePt1Graph.Draw(c,kTRUE,format,2);

    
    sprintf(ylabel,"Events / %.1f GeV/c",hProbePt2v[icat][0]->GetBinWidth(1));
    sprintf(pname,"probeptlog_cat%i",icat);
    CPlot plotProbePt2(pname,"",xlabel,ylabel);
    plotProbePt2.AddHist1D(hProbePt2v[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotProbePt2.AddToStack(hProbePt2v[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotProbePt2.SetLegend(0.70,0.48,0.90,0.78); 
    plotProbePt2.AddTextBox(lumitext,0.60,0.85,0.90,0.79,0);
    plotProbePt2.SetLogy();   
    plotProbePt2.Draw(c,kFALSE,format,1);
    
    CPlot plotProbePt2Graph(pname,"",xlabel,"#chi");
    plotProbePt2Graph.AddHist1D(hProbePt2Diffv[icat],"E",ratioColor);
    plotProbePt2Graph.SetYRange(-8,8);
    plotProbePt2Graph.AddLine(hProbePt2Diffv[icat]->GetXaxis()->GetXmin(), 0,hProbePt2Diffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotProbePt2Graph.AddLine(hProbePt2Diffv[icat]->GetXaxis()->GetXmin(), 5,hProbePt2Diffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotProbePt2Graph.AddLine(hProbePt2Diffv[icat]->GetXaxis()->GetXmin(),-5,hProbePt2Diffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotProbePt2Graph.Draw(c,kTRUE,format,2);
    
    //
    // Probe lepton eta
    //
    if(icat==0) sprintf(xlabel,"trailing muon #eta");
    else        sprintf(xlabel,"probe muon #eta");
    
    sprintf(ylabel,"Events / %.2f",hProbeEtav[icat][0]->GetBinWidth(1));
    sprintf(pname,"probeeta_cat%i",icat);
    CPlot plotProbeEta(pname,"",xlabel,ylabel);
    plotProbeEta.AddHist1D(hProbeEtav[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotProbeEta.AddToStack(hProbeEtav[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotProbeEta.SetLegend(0.70,0.55,0.90,0.85); 
    plotProbeEta.AddTextBox(lumitext,0.21,0.85,0.51,0.79,0);
    plotProbeEta.SetYRange(0,1.8*(hProbeEtav[icat][0]->GetMaximum()));
    plotProbeEta.Draw(c,kFALSE,format,1);
    
    CPlot plotProbeEtaGraph(pname,"",xlabel,"#chi");
    plotProbeEtaGraph.AddHist1D(hProbeEtaDiffv[icat],"E",ratioColor);
    plotProbeEtaGraph.SetYRange(-8,8);
    plotProbeEtaGraph.AddLine(hProbeEtaDiffv[icat]->GetXaxis()->GetXmin(), 0,hProbeEtaDiffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotProbeEtaGraph.AddLine(hProbeEtaDiffv[icat]->GetXaxis()->GetXmin(), 5,hProbeEtaDiffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotProbeEtaGraph.AddLine(hProbeEtaDiffv[icat]->GetXaxis()->GetXmin(),-5,hProbeEtaDiffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotProbeEtaGraph.Draw(c,kTRUE,format,2);
    
    //
    // Probe lepton phi
    //
    if(icat==0) sprintf(xlabel,"trailing muon #phi");
    else        sprintf(xlabel,"probe muon #phi");
    
    sprintf(ylabel,"Events / %.2f",hProbePhiv[icat][0]->GetBinWidth(1));
    sprintf(pname,"probephi_cat%i",icat);
    CPlot plotProbePhi(pname,"",xlabel,ylabel);
    plotProbePhi.AddHist1D(hProbePhiv[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotProbePhi.AddToStack(hProbePhiv[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotProbePhi.SetLegend(0.70,0.55,0.90,0.85); 
    plotProbePhi.AddTextBox(lumitext,0.21,0.85,0.51,0.79,0);
    plotProbePhi.SetYRange(0,1.8*(hProbePhiv[icat][0]->GetMaximum()));
    plotProbePhi.Draw(c,kFALSE,format,1);
    
    CPlot plotProbePhiGraph(pname,"",xlabel,"#chi");
    plotProbePhiGraph.AddHist1D(hProbePhiDiffv[icat],"E",ratioColor);
    plotProbePhiGraph.SetYRange(-8,8);
    plotProbePhiGraph.AddLine(hProbePhiDiffv[icat]->GetXaxis()->GetXmin(), 0,hProbePhiDiffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotProbePhiGraph.AddLine(hProbePhiDiffv[icat]->GetXaxis()->GetXmin(), 5,hProbePhiDiffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotProbePhiGraph.AddLine(hProbePhiDiffv[icat]->GetXaxis()->GetXmin(),-5,hProbePhiDiffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);    
    plotProbePhiGraph.Draw(c,kTRUE,format,2);  

    //
    // PV multiplicity
    //
    sprintf(xlabel,"N_{PV}");
    sprintf(ylabel,"Events");
    sprintf(pname,"npv_cat%i",icat);
    CPlot plotNPV(pname,"",xlabel,ylabel);
    plotNPV.AddHist1D(hNPVv[icat][0],samplev[0]->label,"E");
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      plotNPV.AddToStack(hNPVv[icat][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
    }
    plotNPV.SetLegend(0.70,0.55,0.90,0.85); 
    plotNPV.AddTextBox(lumitext,0.35,0.85,0.65,0.79,0);
    plotNPV.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
    plotNPV.Draw(c,kFALSE,format,1);
    
    CPlot plotNPVGraph(pname,"",xlabel,"#chi");
    plotNPVGraph.AddHist1D(hNPVDiffv[icat],"E",ratioColor);
    plotNPVGraph.SetYRange(-8,8);
    plotNPVGraph.AddLine(hNPVDiffv[icat]->GetXaxis()->GetXmin(), 0,hNPVDiffv[icat]->GetXaxis()->GetXmax(), 0,kBlack,1);
    plotNPVGraph.AddLine(hNPVDiffv[icat]->GetXaxis()->GetXmin(), 5,hNPVDiffv[icat]->GetXaxis()->GetXmax(), 5,kBlack,3);
    plotNPVGraph.AddLine(hNPVDiffv[icat]->GetXaxis()->GetXmin(),-5,hNPVDiffv[icat]->GetXaxis()->GetXmax(),-5,kBlack,3);
    plotNPVGraph.Draw(c,kTRUE,format,2);  
  }


  sprintf(pname,"zmassBB");
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassBBv[0]->GetBinWidth(1));
  CPlot plotMassBB(pname,"","m(#mu^{+}#mu^{-}) [GeV/c^{2}]",ylabel);
  plotMassBB.AddHist1D(hMassBBv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    plotMassBB.AddToStack(hMassBBv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  plotMassBB.SetLegend(0.70,0.55,0.90,0.85); 
  plotMassBB.AddTextBox(lumitext,0.21,0.79,0.51,0.85,0);
  plotMassBB.AddTextBox("barrel-barrel",0.21,0.72,0.51,0.78,0);
  plotMassBB.Draw(c,kFALSE,format,1);
  
  CPlot plotMassBBGraph(pname,"","m(#mu^{+}#mu^{-}) [GeV/c^{2}]","#chi");
  plotMassBBGraph.AddHist1D(hMassBBDiff,"E",ratioColor);
  plotMassBBGraph.SetYRange(-8,8);
  plotMassBBGraph.AddLine(hMassBBDiff->GetXaxis()->GetXmin(), 0,hMassBBDiff->GetXaxis()->GetXmax(), 0,kBlack,1);
  plotMassBBGraph.AddLine(hMassBBDiff->GetXaxis()->GetXmin(), 5,hMassBBDiff->GetXaxis()->GetXmax(), 5,kBlack,3);
  plotMassBBGraph.AddLine(hMassBBDiff->GetXaxis()->GetXmin(),-5,hMassBBDiff->GetXaxis()->GetXmax(),-5,kBlack,3);    
  plotMassBBGraph.Draw(c,kTRUE,format,2);

  
  sprintf(pname,"zmassBE");
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassBEv[0]->GetBinWidth(1));
  CPlot plotMassBE(pname,"","m(#mu^{+}#mu^{-}) [GeV/c^{2}]",ylabel);
  plotMassBE.AddHist1D(hMassBEv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    plotMassBE.AddToStack(hMassBEv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  plotMassBE.SetLegend(0.70,0.55,0.90,0.85); 
  plotMassBE.AddTextBox(lumitext,0.21,0.79,0.51,0.85,0);
  plotMassBE.AddTextBox("barrel-endcap",0.21,0.72,0.51,0.78,0);
  plotMassBE.Draw(c,kFALSE,format,1);
  
  CPlot plotMassBEGraph(pname,"","m(#mu^{+}#mu^{-}) [GeV/c^{2}]","#chi");
  plotMassBEGraph.AddHist1D(hMassBEDiff,"E",ratioColor);
  plotMassBEGraph.SetYRange(-8,8);
  plotMassBEGraph.AddLine(hMassBEDiff->GetXaxis()->GetXmin(), 0,hMassBEDiff->GetXaxis()->GetXmax(), 0,kBlack,1);
  plotMassBEGraph.AddLine(hMassBEDiff->GetXaxis()->GetXmin(), 5,hMassBEDiff->GetXaxis()->GetXmax(), 5,kBlack,3);
  plotMassBEGraph.AddLine(hMassBEDiff->GetXaxis()->GetXmin(),-5,hMassBEDiff->GetXaxis()->GetXmax(),-5,kBlack,3);    
  plotMassBEGraph.Draw(c,kTRUE,format,2);  

  
  sprintf(pname,"zmassEE");  
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassEEv[0]->GetBinWidth(1));
  CPlot plotMassEE(pname,"","m(#mu^{+}#mu^{-}) [GeV/c^{2}]",ylabel);
  plotMassEE.AddHist1D(hMassEEv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    plotMassEE.AddToStack(hMassEEv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  plotMassEE.SetLegend(0.70,0.55,0.90,0.85); 
  plotMassEE.AddTextBox(lumitext,0.21,0.79,0.51,0.85,0);
  plotMassEE.AddTextBox("endcap-endcap",0.21,0.72,0.51,0.78,0);
  plotMassEE.Draw(c,kFALSE,format,1);
  
  CPlot plotMassEEGraph(pname,"","m(#mu^{+}#mu^{-}) [GeV/c^{2}]","#chi");
  plotMassEEGraph.AddHist1D(hMassEEDiff,"E",ratioColor);
  plotMassEEGraph.SetYRange(-8,8);
  plotMassEEGraph.AddLine(hMassEEDiff->GetXaxis()->GetXmin(), 0,hMassEEDiff->GetXaxis()->GetXmax(), 0,kBlack,1);
  plotMassEEGraph.AddLine(hMassEEDiff->GetXaxis()->GetXmin(), 5,hMassEEDiff->GetXaxis()->GetXmax(), 5,kBlack,3);
  plotMassEEGraph.AddLine(hMassEEDiff->GetXaxis()->GetXmin(),-5,hMassEEDiff->GetXaxis()->GetXmax(),-5,kBlack,3);    
  plotMassEEGraph.Draw(c,kTRUE,format,2);  
  
    
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
      
  gBenchmark->Show("plotZmm"); 
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
  htmlfile << "<head><title>Zmm</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zmass_cat0.png\"><img src=\"plots/zmass_cat0.png\" alt=\"plots/zmass_cat0.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zmassBB.png\"><img src=\"plots/zmassBB.png\" alt=\"plots/zmassBB.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zmassBE.png\"><img src=\"plots/zmassBE.png\" alt=\"plots/zmassBE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zmassEE.png\"><img src=\"plots/zmassEE.png\" alt=\"plots/zmassEE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"cat0.html\">Selected Zmm</a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
    
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/zmass_cat1.png\"><img src=\"plots/zmass_cat1.png\" alt=\"plots/zmass_cat1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/zmass_cat2.png\"><img src=\"plots/zmass_cat2.png\" alt=\"plots/zmass_cat2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/zmass_cat3.png\"><img src=\"plots/zmass_cat3.png\" alt=\"plots/zmass_cat3.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/zmass_cat4.png\"><img src=\"plots/zmass_cat4.png\" alt=\"plots/zmass_cat4.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"plots/zmass_cat5.png\"><img src=\"plots/zmass_cat5.png\" alt=\"plots/zmass_cat5.png\" width=\"100%\"></a></td>" << endl;  
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"cat1.html\">MuMu2HLT</a></td>" << endl;
  htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"cat2.html\">MuMu1HLT</a></td>" << endl;
  htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"cat3.html\">MuMuNoSel</a></td>" << endl;
  htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"cat4.html\">MuSta</a></td>" << endl;
  htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"cat5.html\">MuTrk</a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
  
  for(UInt_t icat=0; icat<6; icat++) {
    sprintf(htmlfname,"%s/cat%i.html",outDir.Data(),icat);
    htmlfile.open(htmlfname);
    htmlfile << "<!DOCTYPE html" << endl;
    htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
    htmlfile << "<html>" << endl;
    htmlfile << "<head><title>Zmm</title></head>" << endl;
    htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  
    htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zmass_cat"    << icat << ".png\"><img src=\"plots/zmass_cat"    << icat << ".png\" alt=\"plots/zmass_cat"    << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zmasslog_cat" << icat << ".png\"><img src=\"plots/zmasslog_cat" << icat << ".png\" alt=\"plots/zmasslog_cat" << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zpt_cat"      << icat << ".png\"><img src=\"plots/zpt_cat"      << icat << ".png\" alt=\"plots/zpt_cat"      << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zptlog_cat"   << icat << ".png\"><img src=\"plots/zptlog_cat"   << icat << ".png\" alt=\"plots/zptlog_cat"   << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zy_cat"    << icat << ".png\"><img src=\"plots/zy_cat"    << icat << ".png\" alt=\"plots/zy_cat"    << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zphi_cat"  << icat << ".png\"><img src=\"plots/zphi_cat"  << icat << ".png\" alt=\"plots/zphi_cat"  << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zdphi_cat" << icat << ".png\"><img src=\"plots/zdphi_cat" << icat << ".png\" alt=\"plots/zdphi_cat" << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/zdr_cat"   << icat << ".png\"><img src=\"plots/zdr_cat"   << icat << ".png\" alt=\"plots/zdr_cat"   << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/npv_cat" << icat << ".png\"><img src=\"plots/npv_cat" << icat << ".png\" alt=\"plots/npv_cat" << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "</table>" << endl;
    htmlfile << "<hr />" << endl;
    
    htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/tagpt_cat"    << icat << ".png\"><img src=\"plots/tagpt_cat"    << icat << ".png\" alt=\"plots/tagpt_cat"    << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/tagptlog_cat" << icat << ".png\"><img src=\"plots/tagptlog_cat" << icat << ".png\" alt=\"plots/tagptlog_cat" << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/tageta_cat"   << icat << ".png\"><img src=\"plots/tageta_cat"   << icat << ".png\" alt=\"plots/tageta_cat"   << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/tagphi_cat"   << icat << ".png\"><img src=\"plots/tagphi_cat"   << icat << ".png\" alt=\"plots/tagphi_cat"   << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/probept_cat"    << icat << ".png\"><img src=\"plots/probept_cat"    << icat << ".png\" alt=\"plots/probept_cat"    << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/probeptlog_cat" << icat << ".png\"><img src=\"plots/probeptlog_cat" << icat << ".png\" alt=\"plots/probeptlog_cat" << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/probeeta_cat"   << icat << ".png\"><img src=\"plots/probeeta_cat"   << icat << ".png\" alt=\"plots/probeeta_cat"   << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/probephi_cat"   << icat << ".png\"><img src=\"plots/probephi_cat"   << icat << ".png\" alt=\"plots/probephi_cat"   << icat << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "</table>" << endl;
    htmlfile << "<hr />" << endl;
  
    htmlfile << "</body>" << endl;
    htmlfile << "</html>" << endl;
    htmlfile.close();
  }
}
