//================================================================================================
//
// Select W->munu candidates
//
//  * outputs ROOT files of events passing selection
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include "TLorentzVector.h"         // 4-vector class
#include "TH1D.h"
#include "TRandom.h"

#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples
#include "../Utils/LeptonCorr.hh"   // muon scale and resolution corrections

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions
#endif

//=== MAIN MACRO ================================================================================================= 

void selectWmGen(const TString conf="wm.conf", // input file
              const TString outputDir=".",  // output directory
	      const Bool_t  doScaleCorr=0   // apply energy scale corrections?
) {
  gBenchmark->Start("selectWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT    = 20;
  const Double_t ETA_CUT   = 2.4;
  const Double_t MUON_MASS = 0.105658369;

  const Double_t VETO_PT   = 10;
  const Double_t VETO_ETA  = 2.4;

  const Int_t BOSON_ID  = 24;
  const Int_t LEPTON_ID = 13;

  // load trigger menu
  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/pileup_rw_baconDY.root", "read");

  TH1D *h_rw = (TH1D*) f_rw->Get("h_rw_golden");
  TH1D *h_rw_up = (TH1D*) f_rw->Get("h_rw_up_golden");
  TH1D *h_rw_down = (TH1D*) f_rw->Get("h_rw_down_golden");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  //
  // parse .conf file
  //
  confParse(conf, snamev, samplev);

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  const TString ntupDir = outputDir + TString("/ntuples");
  gSystem->mkdir(ntupDir,kTRUE);
  
  //
  // Declare output ntuple variables
  //
  
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  npv, npu;
  UInt_t  triggerDec;
  UInt_t  goodPV;
  UInt_t  matchTrigger;
  UInt_t  ngenlep;
  TLorentzVector *genlep=0;
  TLorentzVector *genV=0;
  Int_t   genq;
  UInt_t nlep;
  TLorentzVector *lep=0;
  Int_t   q;
  Float_t scale1fbGen,scale1fb,scale1fbUp,scale1fbDown;

  std::vector<float> *lheweight = new std::vector<float>();
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info  = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen  = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr    = new TClonesArray("baconhep::TMuon");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;
  
  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {

    // Assume signal sample is given name "wm" -- flag to store GEN W kinematics
    Bool_t isSignal = (snamev[isam].CompareTo("wm",TString::kIgnoreCase)==0);
    // flag to reject W->mnu events when selecting wrong flavor background events
    Bool_t isWrongFlavor = (snamev[isam].CompareTo("wx",TString::kIgnoreCase)==0);
    
    CSample* samp = samplev[isam];

    //
    // Set up output ntuple
    //
    TString outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.raw.root");
    cout << outfilename << endl;
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");
    outTree->Branch("runNum",     &runNum,     "runNum/i");      // event run number
    outTree->Branch("lumiSec",    &lumiSec,    "lumiSec/i");     // event lumi section
    outTree->Branch("evtNum",     &evtNum,     "evtNum/i");      // event number
    outTree->Branch("matchGen",   &matchGen,   "matchGen/i");    // event has both leptons matched to MC Z->ll
    outTree->Branch("npv",        &npv,        "npv/i");         // number of primary vertices
    outTree->Branch("npu",        &npu,        "npu/i");         // number of in-time PU events (MC)
    
    outTree->Branch("triggerDec",   &triggerDec,   "triggerDec/i");    // event pass the trigger
    outTree->Branch("goodPV",   &goodPV,   "goodPV/i");    // event has a good PV
    outTree->Branch("matchTrigger",   &matchTrigger,   "matchTrigger/i");    // event has at least one lepton matched to the trigger
    outTree->Branch("ngenlep",     &ngenlep,     "ngenlep/i");      // number of gen leptons
    outTree->Branch("genlep",   "TLorentzVector",  &genlep);     // gen lepton1 4-vector
    outTree->Branch("genq",          &genq,         "genq/I");          // charge of lepton1
    outTree->Branch("genV",   "TLorentzVector",  &genV);     // gen lepton1 4-vector
//     outTree->Branch("gen",          &genV,         "genq/I");          // charge of lepton1
    outTree->Branch("nlep",     &nlep,     "nlep/i");      // number of leptons
    outTree->Branch("lep",       "TLorentzVector",  &lep);     // lepton1 4-vector
    outTree->Branch("q",          &q,         "q/I");          // charge of lepton1
    outTree->Branch("scale1fbGen",   &scale1fbGen,   "scale1fbGen/F");    // event weight per 1/fb (MC)
    outTree->Branch("scale1fb",   &scale1fb,   "scale1fb/F");    // event weight per 1/fb (MC)
    outTree->Branch("scale1fbUp",    &scale1fbUp,   "scale1fbUp/F");    // event weight per 1/fb (MC)
    outTree->Branch("scale1fbDown",    &scale1fbDown,   "scale1fbDown/F");    // event weight per 1/fb (MC)
    outTree->Branch("lheweight",  &lheweight);


    //
    // loop through files
    //
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      
      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... "; cout.flush();

      infile = TFile::Open(samp->fnamev[ifile]); 
      assert(infile);

      eventTree = (TTree*)infile->Get("Events");
      assert(eventTree);  
      eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");
      Bool_t hasGen = eventTree->GetBranchStatus("GenEvtInfo");
      TBranch *genBr=0, *genPartBr=0;
      if(hasGen) {
        eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
        eventTree->SetBranchAddress("GenParticle",&genPartArr); genPartBr = eventTree->GetBranch("GenParticle");
      }

      // Compute MC event weight per 1/fb
      const Double_t xsec = samp->xsecv[ifile];
      Double_t totalWeightGen=0;
      Double_t totalWeight=0;
      Double_t totalWeightUp=0;
      Double_t totalWeightDown=0;
      Double_t puWeight=0;
      Double_t puWeightUp=0;
      Double_t puWeightDown=0;

      if (hasGen) {
        for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
          infoBr->GetEntry(ientry);
          genBr->GetEntry(ientry);
          puWeight = h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));
          puWeightUp = h_rw_up->GetBinContent(h_rw_up->FindBin(info->nPUmean));
          puWeightDown = h_rw_down->GetBinContent(h_rw_down->FindBin(info->nPUmean));
          totalWeightGen+=gen->weight;
          totalWeight+=gen->weight*puWeight;
          totalWeightUp+=gen->weight*puWeightUp;
          totalWeightDown+=gen->weight*puWeightDown;
        }
      }

      //
      // loop over events
      //
      Double_t nsel=0, nselvar=0;

     // for(UInt_t ientry=0; ientry<1000; ientry++) {
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
        infoBr->GetEntry(ientry);

        if(ientry%1000000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;

        Double_t weightGen=1;
        Double_t weight=1;
        Double_t weightUp=1;
        Double_t weightDown=1;
        if(xsec>0 && totalWeightGen>0) weightGen = xsec/totalWeightGen;
        if(xsec>0 && totalWeight>0) weight = xsec/totalWeight;
        if(xsec>0 && totalWeightUp>0) weightUp = xsec/totalWeightUp;
        if(xsec>0 && totalWeightDown>0) weightDown = xsec/totalWeightDown;
        if(hasGen) {
          genPartArr->Clear();
          genBr->GetEntry(ientry);
          genPartBr->GetEntry(ientry);
          puWeight = h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));
          puWeightUp = h_rw_up->GetBinContent(h_rw_up->FindBin(info->nPUmean));
          puWeightDown = h_rw_down->GetBinContent(h_rw_down->FindBin(info->nPUmean));
          weightGen*=gen->weight;
          weight*=gen->weight*puWeight;
          weightUp*=gen->weight*puWeightUp;
          weightDown*=gen->weight*puWeightDown;
        }


        // veto w -> xv decays for signal and w -> mv for bacground samples (needed for inclusive WToLNu sample)
        if (isWrongFlavor && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;
        else if (isSignal && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue; 

        // trigger requirement               
        // if (!isMuonTrigger(triggerMenu, info->triggerBits)) continue;
        Bool_t passMuTrigger = kFALSE;
	    if (isMuonTrigger(triggerMenu, info->triggerBits)) passMuTrigger=kTRUE;
	
        // good vertex requirement
        // Bool_t hasGoodPV = kFALSE;
      
        // good vertex requirement
        Bool_t hasGoodPV = kFALSE;
        if((info->hasGoodPV)) hasGoodPV=kTRUE;
           
        //
        // SELECTION PROCEDURE:
        //  (1) Look for 1 good muon matched to trigger
        //  (2) Reject event if another muon is present passing looser cuts
        //
        muonArr->Clear();
        muonBr->GetEntry(ientry);

        Int_t nLooseLep=0;
        const baconhep::TMuon *goodMuon=0;
        Bool_t passSel=kFALSE;
        Bool_t hasTriggerMatch=kFALSE;
        for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[i]);

          if(fabs(mu->eta) > ETA_CUT)         continue;  // lepton |eta| cut
	      if(mu->pt     < PT_CUT)          continue;  // lepton pT cut   
          if(!passMuonID(mu))                 continue;  // lepton selection
          if(!isMuonTriggerObj(triggerMenu, mu->hltMatchBits, kFALSE)) continue;
          if(isMuonTriggerObj(triggerMenu, mu->hltMatchBits, kFALSE)) hasTriggerMatch=kTRUE;

          passSel=kTRUE;
          goodMuon = mu;
        }

        if(passSel) {
          /******** We have a W candidate! HURRAY! ********/
          nsel+=weight;
          nselvar+=weight*weight;
          TLorentzVector vlep; 
          vlep.SetPtEtaPhiM(goodMuon->pt, goodMuon->eta, goodMuon->phi, MUON_MASS); 
          

          Int_t glepq1=-99;
          Int_t glepq2=-99;
          TLorentzVector *gvec=new TLorentzVector(0,0,0,0);
          TLorentzVector *glep1=new TLorentzVector(0,0,0,0);
          TLorentzVector *glep2=new TLorentzVector(0,0,0,0);
          genV      = new TLorentzVector(0,0,0,0);
          genlep    = new TLorentzVector(0,0,0,0);
          Bool_t hasGenMatch = kFALSE;
          Bool_t hasTriggerMatch = kFALSE;
            
          if(isSignal && hasGen) {
            toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2,&glepq1,&glepq2,1);
            
            hasGenMatch = ( ((glep1) && toolbox::deltaR(vlep.Eta(), vlep.Phi(), glep1->Eta(), glep1->Phi())<0.5) );
            if (gvec && glep1) {
              genV->SetPtEtaPhiM(gvec->Pt(),gvec->Eta(),gvec->Phi(),gvec->M());
              genlep->SetPtEtaPhiM(glep1->Pt(),glep1->Eta(),glep1->Phi(),glep1->M());
              //genV->Print();
              //genlep->Print(); 
              // genVPt    = gvec->Pt();
              // genVPhi   = gvec->Phi();
              // genVy     = gvec->Rapidity();
              // genVMass  = gvec->M();
              // genLepPt  = glep1->Pt();
              // genLepPhi = glep1->Phi();
            }
            // id_1      = gen->id_1;
            // id_2      = gen->id_2;
            // x_1       = gen->x_1;
            // x_2       = gen->x_2;
            // xPDF_1    = gen->xPDF_1;
            // xPDF_2    = gen->xPDF_2;
            // scalePDF  = gen->scalePDF;
            // weightPDF = gen->weight;
            delete gvec;
            delete glep1;
            delete glep2;
            gvec=0; glep1=0; glep2=0;
          }
          //genV->Print();
          //genlep->Print();
          //
          // Fill tree
          //
          runNum    = info->runNum;
          lumiSec   = info->lumiSec;
          evtNum    = info->evtNum;
          if (hasGenMatch) matchGen=1;
          else matchGen=0;
          if (passMuTrigger) triggerDec=1;
          else triggerDec=0;
          if (hasGoodPV) goodPV=1;
          else goodPV=0;
          if (hasTriggerMatch) matchTrigger=1;
          else matchTrigger=0;
          
          vertexArr->Clear();
          vertexBr->GetEntry(ientry);

          npv       = vertexArr->GetEntries();
          npu	    = info->nPUmean;
          //genV      = new TLorentzVector(0,0,0,0);
          //genlep    = new TLorentzVector(0,0,0,0);
          // lep = &vlep;
          // genVPt    = -999;
          // genVPhi   = -999;
          // genVy     = -999;
          // genVMass  = -999;
          // genLepPt  = -999;
          // genLepPhi = -999;
          // id_1      = -999;
          // id_2      = -999;
          // x_1       = -999;
          // // x_2       = -999;
          // xPDF_1    = -999;
          // xPDF_2    = -999;
          // scalePDF  = -999;
          // weightPDF = -999;
          
          scale1fb     = weight;
          scale1fbGen  = weightGen;
          scale1fbUp   = weightUp;
          scale1fbDown = weightDown;
          
          lheweight->clear();
          for (int j = 0; j<=110; j++){
            lheweight->push_back(gen->lheweight[j]);
          }
      
          // sumEt    = 0;
          q        = goodMuon->q;
          lep      = &vlep;
          

          outTree->Fill();
          delete genV;
          delete genlep;
          genV=0, genlep=0, lep=0;
        }
      }
      delete infile;
      infile=0, eventTree=0;    

      cout << nsel  << " +/- " << sqrt(nselvar);
      if(isam!=0) cout << " per 1/pb";
      cout << endl;
    }
    outFile->Write();
    outFile->Close();
  }
  delete h_rw;
  delete h_rw_up;
  delete h_rw_down;
  delete f_rw;
  delete info;
  delete gen;
  delete genPartArr;
  delete muonArr;
  delete vertexArr;
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " W -> mu nu" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectWm"); 
}
