//================================================================================================
//
// Select W->enu candidates
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
#include "TCanvas.h"
#include "TRandom.h"
#include "TGraph.h"


// define structures to read in ntuple

#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

// lumi section selection with JSON files
#include "BaconAna/Utils/interface/RunLumiRangeMap.hh"

#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples
#include "../Utils/LeptonCorr.hh"   // electron scale and resolution corrections
#include "../EleScale/EnergyScaleCorrection.h" //EGMSmear

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions
#include "../Utils/PrefiringEfficiency.cc"      // prefiring efficiency functions
#endif


//=== MAIN MACRO ================================================================================================= 

void selectWe(const TString  conf        ="we.conf", // input file
const TString  outputDir   =".",  // output directory
const Bool_t   doScaleCorr =0,   // apply energy scale corrections?
const Int_t    sigma       =0,
const Bool_t   doPU        =0,
const Bool_t   is13TeV     =1,
const Int_t    NSEC        =1,
const Int_t    ITH         =0 ) {
  gBenchmark->Start("selectWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  const Double_t PT_CUT   = 25;
  const Double_t ETA_CUT  = 2.4;

  const Double_t VETO_PT   = 10;
  const Double_t VETO_ETA  = 2.4;

  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  const Double_t ELE_MASS = 0.000511;
  const Int_t BOSON_ID  = 24;
  const Int_t LEPTON_ID = 11;

  // load trigger menu
  const baconhep::TTrigger triggerMenu("/afs/cern.ch/work/s/sabrandt/public/SM/LowPU/CMSSW_9_4_12/src/BaconAna/DataFormats/data/HLT_50nsGRun");

  const TString prefireFileName = "../Utils/All2017Gand2017HPrefiringMaps.root";
  PrefiringEfficiency pfire( prefireFileName.Data() , (is13TeV ? "2017H" : "2017G"));

  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/puWeights_76x.root", "read");
  TH1D *h_rw      = (TH1D*) f_rw->Get("puWeights");
  TH1D *h_rw_up   = (TH1D*) f_rw->Get("puWeightsUp");
  TH1D *h_rw_down = (TH1D*) f_rw->Get("puWeightsDown");

  const TString corrFiles = "../EleScale/Run2017_LowPU_v2";
  EnergyScaleCorrection eleCorr( corrFiles.Data(), EnergyScaleCorrection::ECALELF); 

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  // parse .conf file
  confParse(conf, snamev, samplev);
  const Bool_t hasData = (samplev[0]->fnamev.size()>0);

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  const TString ntupDir = outputDir + TString("/ntuples_") + Form("%d",ITH) + TString("_") + Form("%d",NSEC);
  gSystem->mkdir(ntupDir,kTRUE);

  // Declare output ntuple variables
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  TLorentzVector *genV=0, *genLep=0, *genNu = 0;
  Float_t scale1fb, scale1fbUp, scale1fbDown;
  Float_t prefireWeight=1, prefireUp=1,    prefireDown=1;
  Float_t prefirePhoton=1, prefirePhotUp=1, prefirePhotDown=1;
  Float_t prefireJet=1,    prefireJetUp=1,  prefireJetDown=1;
  Float_t met, metPhi;//, mt, u1, u2;
  Float_t puppiMet, puppiMetPhi;//, puppiMt, puppiU1, puppiU2;
  Int_t   q;
  TLorentzVector *lep=0, *lep_raw=0;
  Float_t lepError=0;
  Float_t pfCombIso;
  TLorentzVector *sc=0;

  // Data structures to store info from TTrees
  baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr     = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *electronArr    = new TClonesArray("baconhep::TElectron");
  TClonesArray *muonArr        = new TClonesArray("baconhep::TMuon");
  TClonesArray *scArr          = new TClonesArray("baconhep::TPhoton");
  TClonesArray *vertexArr      = new TClonesArray("baconhep::TVertex");
  TClonesArray *jetArr         = new TClonesArray("baconhep::TJet");

  TFile *infile=0;
  TTree *eventTree=0;

  // loop over samples
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    // Assume data sample is first sample in .conf file
    // If sample is empty (i.e. contains no ntuple files), skip to next sample
    Bool_t isData=kFALSE;
    if(isam==0 && !hasData) continue;
    else if (isam==0) isData=kTRUE;

    Bool_t isSignal        = (snamev[isam].Contains("we"));
    Bool_t isWrongFlavor   = (snamev[isam].Contains("wx"));

    Bool_t isRecoil = (snamev[isam].Contains("zxx")||isSignal||isWrongFlavor);
    Bool_t noGen    = (snamev[isam].Contains("zz")||snamev[isam].Contains("wz")||snamev[isam].Contains("ww"));
    CSample* samp = samplev[isam];

    // Set up output ntuple
    TString outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    if(isam!=0 && !doScaleCorr) outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.raw.root");
    cout << outfilename << endl;
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");
    outTree->Branch("runNum",     &runNum,     "runNum/i");      // event run number
    outTree->Branch("lumiSec",    &lumiSec,    "lumiSec/i");     // event lumi section
    outTree->Branch("evtNum",     &evtNum,     "evtNum/i");      // event number
    outTree->Branch("npv",        &npv,        "npv/i");         // number of primary vertices
    outTree->Branch("npu",        &npu,        "npu/i");         // number of in-time PU events (MC)
    outTree->Branch("genV",       "TLorentzVector", &genV);      // GEN boson 4-vector (signal MC)
    outTree->Branch("genLep",     "TLorentzVector", &genLep);    // GEN lepton 4-vector (signal MC)
    outTree->Branch("genNu",      "TLorentzVector", &genNu);    // GEN lepton 4-vector (signal MC)
    outTree->Branch("prefireWeight", &prefireWeight, "prefireWeight/F");
    outTree->Branch("prefireUp",     &prefireUp,     "prefireUp/F");
    outTree->Branch("prefireDown",   &prefireDown,   "prefireDown/F");
    outTree->Branch("prefirePhoton", &prefirePhoton, "prefirePhoton/F");
    outTree->Branch("prefirePhotUp",     &prefirePhotUp,     "prefirePhotUp/F");
    outTree->Branch("prefirePhotDown",   &prefirePhotDown,   "prefirePhotDown/F");
    outTree->Branch("prefireJet",    &prefireJet,    "prefireJet/F");
    outTree->Branch("prefireJetUp",  &prefireJetUp,  "prefireJetUp/F");
    outTree->Branch("prefireJetDown",&prefireJetDown,"prefireJetDown/F");
    outTree->Branch("scale1fb",   &scale1fb,   "scale1fb/F");    // event weight per 1/fb (MC)
    outTree->Branch("scale1fbUp",   &scale1fbUp,   "scale1fbUp/F");    // event weight per 1/fb (MC)
    outTree->Branch("scale1fbDown",   &scale1fbDown,   "scale1fbDown/F");    // event weight per 1/fb (MC)
    outTree->Branch("met",        &met,        "met/F");         // MET
    outTree->Branch("metPhi",     &metPhi,     "metPhi/F");      // phi(MET)
    outTree->Branch("puppiMet",    &puppiMet,   "puppiMet/F");      // Puppi MET
    outTree->Branch("puppiMetPhi", &puppiMetPhi,"puppiMetPhi/F");   // phi(Puppi MET)
    outTree->Branch("q",          &q,          "q/I");           // lepton charge
    outTree->Branch("lep",       "TLorentzVector", &lep);        // lepton 4-vector
    outTree->Branch("lep_raw",       "TLorentzVector", &lep_raw);        // lepton 4-vector
    outTree->Branch("lepError",   &lepError,   "lepError/F");      // track isolation of tag lepton
    outTree->Branch("pfCombIso",  &pfCombIso,  "pfCombIso/F");   // PF combined isolation of electron
    outTree->Branch("sc",        "TLorentzVector", &sc);         // supercluster 4-vector

    TH1D* hGenWeights = new TH1D("hGenWeights","hGenWeights",10,-10.,10.);
    //
    // loop through files
    //
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... "; cout.flush();      
      infile = TFile::Open(samp->fnamev[ifile]); 
      assert(infile);

      Bool_t hasJSON = kFALSE;
      baconhep::RunLumiRangeMap rlrm;
      if(!samp->jsonv[ifile].Contains("NONE")) {
        hasJSON = kTRUE;
        rlrm.addJSONFile(samp->jsonv[ifile].Data()); 
      }

      eventTree = (TTree*)infile->Get("Events");
      assert(eventTree);
      Bool_t hasJet = eventTree->GetBranchStatus("AK4");
      eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Muon", &muonArr);         TBranch *muonBr     = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PV",   &vertexArr);       TBranch *vertexBr = eventTree->GetBranch("PV");
      eventTree->SetBranchAddress("Photon",   &scArr);       TBranch *scBr       = eventTree->GetBranch("Photon");
      if(hasJet)eventTree->SetBranchAddress("AK4",      &jetArr     ); TBranch *jetBr      = eventTree->GetBranch("AK4");

      Bool_t hasGen = (eventTree->GetBranchStatus("GenEvtInfo")&&!noGen);
      TBranch *genBr=0, *genPartBr=0;
      if(hasGen) {
        eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
        eventTree->SetBranchAddress("GenParticle",&genPartArr); genPartBr = eventTree->GetBranch("GenParticle");
      }
      // Compute MC event weight per 1/fb
      const Double_t xsec = samp->xsecv[ifile];
      Double_t puWeight=0, puWeightUp=0, puWeightDown=0;

      //////////////////////////////////////////////////////////
      // Real selection loop
      //
      // loop over events
      //
      double frac = 1.0/NSEC;
      UInt_t IBEGIN = frac* ITH   *eventTree->GetEntries();
      UInt_t IEND   = frac*(ITH+1)*eventTree->GetEntries();

      Double_t nsel=0, nselvar=0;
      for(UInt_t ientry=IBEGIN; ientry < IEND; ientry++) {
        infoBr->GetEntry(ientry);

        int printIndex = (int)(eventTree->GetEntries()*0.01);
        if(ientry%printIndex==0) cout << "Processing event " << ientry << ". " << (int)(100*(ientry/(double)eventTree->GetEntries())) << " percent done with this file." << endl;

        Double_t weight=xsec, weightUp=xsec, weightDown=xsec;      
        if(hasGen) {
          genPartArr->Clear();
          genBr->GetEntry(ientry);
          genPartBr->GetEntry(ientry);
          puWeight     = doPU ? h_rw     ->GetBinContent(h_rw     ->FindBin(info->nPUmean)) : 1.;
          puWeightUp   = doPU ? h_rw_up  ->GetBinContent(h_rw_up  ->FindBin(info->nPUmean)) : 1.;
          puWeightDown = doPU ? h_rw_down->GetBinContent(h_rw_down->FindBin(info->nPUmean)) : 1.;
          hGenWeights->Fill(0.0,gen->weight);
          weight    *=gen->weight*puWeight;
          weightUp  *=gen->weight*puWeightUp;
          weightDown*=gen->weight*puWeightDown;
        } else {
          hGenWeights->Fill(0.0,1.0);
        }
        // veto w -> xv decays for signal and w -> mv for bacground samples (needed for inclusive WToLNu sample)
        if (isWrongFlavor && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;
        else if (isSignal && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
        // check for certified lumi (if applicable)
        baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
        if(hasJSON && !rlrm.hasRunLumi(rl)) continue;  
        // trigger requirement               
        if (!isEleTrigger(triggerMenu, info->triggerBits, isData,is13TeV)) continue;
        // good vertex requirement
        if(!(info->hasGoodPV)) continue;
        //
        // SELECTION PROCEDURE:
        //  (1) Look for 1 good electron matched to trigger
        //  (2) Reject event if another electron is present passing looser cuts
        //
        electronArr->Clear();           electronBr->GetEntry(ientry);
        scArr      ->Clear();           scBr      ->GetEntry(ientry);
        muonArr    ->Clear();           muonBr    ->GetEntry(ientry);
        jetArr     ->Clear(); if(hasJet)jetBr     ->GetEntry(ientry);

        Int_t nLooseLep=0;
        const baconhep::TElectron *goodEle=0;
        TLorentzVector vEle(0,0,0,0);
        TLorentzVector vGoodEle(0,0,0,0);
        Bool_t passSel=kFALSE;
        double eleRamdom = gRandom->Gaus(0,1);

        for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[i]);
          if(fabs(mu->eta) > VETO_ETA) continue; // loose lepton |eta| cut
          if(mu->pt        < VETO_PT ) continue; // loose lepton pT cut
          if(passMuonLooseID(mu)) nLooseLep++;   // loose lepton selection
          if(nLooseLep>0) {  // extra lepton veto
            passSel=kFALSE;
            break;
          }
        }

        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
          const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[i]);
          vEle.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, ELE_MASS);
          // check ECAL gap
          // if(fabs(vEle.Eta())>=ECAL_GAP_LOW && fabs(vEle.Eta())<=ECAL_GAP_HIGH) continue;

          if(doScaleCorr && (ele->r9 < 1.)){
            float eleAbsEta   = fabs(vEle.Eta());
            double eTregress = ele->ecalEnergy/cosh(fabs(ele->eta));
            if(snamev[isam].Contains("data")) {//Data
              int runNumber = is13TeV ? info->runNum : 306936 ;
              float eleScale = eleCorr.scaleCorr(runNumber, eTregress, eleAbsEta, ele->r9);
              float eleError = eleCorr.scaleCorrUncert(runNumber, eTregress, eleAbsEta, ele->r9);
              (vEle) *= eleScale * (1 + sigma*eleError);
              lepError = eleError;
            } else {//MC
              float eleSmear = eleCorr.smearingSigma(info->runNum, eTregress, eleAbsEta, ele->r9, 12, sigma, 0.);
              (vEle) *= 1. + eleSmear * eleRamdom;
              float eleSmearEP = eleCorr.smearingSigma(info->runNum, eTregress, eleAbsEta,  ele->r9, 12, 1., 0.);
              float eleSmearEM = eleCorr.smearingSigma(info->runNum, eTregress, eleAbsEta,  ele->r9, 12, -1., 0.);
              lepError = eleRamdom * std::hypot(eleSmearEP - eleSmear, eleSmearEM - eleSmear);
            }
          }

          // apply scale and resolution corrections to MC
          if(fabs(vEle.Eta())    > VETO_ETA) continue;
          if(vEle.Pt()           < VETO_PT)  continue; 
          if(passEleLooseID(ele,vEle, info->rhoIso)) nLooseLep++;
          if(nLooseLep>1) {  // extra lepton veto
            passSel=kFALSE;
            break;
          }
          if(vEle.Pt()           < PT_CUT)     continue;  // lepton pT cut
          if(fabs(vEle.Eta())    > ETA_CUT)    continue;  // lepton |eta| cut
          if(!passEleMediumID(ele, vEle, info->rhoIso))     continue;  // lepton selection
          if(!isEleTriggerObj(triggerMenu, ele->hltMatchBits, kFALSE, isData, is13TeV)) continue;
          passSel=kTRUE;
          goodEle = ele;  
          vGoodEle = vEle;
        }

        if(passSel) {
          //******* We have a W candidate! HURRAY! ********
          nsel   +=isData ? 1 : weight;
          nselvar+=isData ? 1 : weight*weight;

          if(!isData){
            pfire.setObjects(scArr,jetArr);
            pfire.computePhotonsOnly(prefirePhoton, prefirePhotUp, prefirePhotDown);
            pfire.computeJetsOnly   (prefireJet   , prefireJetUp , prefireJetDown );
            pfire.computeFullPrefire(prefireWeight, prefireUp    , prefireDown    );
          }

          TLorentzVector vLep(0,0,0,0); TLorentzVector vSC(0,0,0,0); TLorentzVector vLep_raw(0,0,0,0);
          vLep = vGoodEle;
          vLep_raw.SetPtEtaPhiM(goodEle->pt,goodEle->eta,goodEle->phi,ELE_MASS);

          //
          // Fill tree
          //
          runNum    = info->runNum;
          lumiSec   = info->lumiSec;
          evtNum    = info->evtNum;

          vertexArr->Clear();
          vertexBr->GetEntry(ientry);

          npv      = vertexArr->GetEntries();
          npu	    = info->nPUmean;
          genV      = new TLorentzVector(0,0,0,0);
          genLep    = new TLorentzVector(0,0,0,0);
          genNu    = new TLorentzVector(0,0,0,0);

          if(isRecoil && hasGen) {
            Int_t glepq1=-99;
            Int_t glepq2=-99;
            TLorentzVector *gvec=new TLorentzVector(0,0,0,0);
            TLorentzVector *glep1=new TLorentzVector(0,0,0,0);
            TLorentzVector *glep2=new TLorentzVector(0,0,0,0);
            toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2,&glepq1,&glepq2,1);
            if(snamev[isam].Contains("zxx")){ // DY only
              toolbox::fillGen(genPartArr, 23, gvec, glep1, glep2,&glepq1,&glepq2,1);
            }

            TLorentzVector tvec=*glep1+*glep2;
            genV=new TLorentzVector(0,0,0,0);
            genV->SetPtEtaPhiM(tvec.Pt(), tvec.Eta(), tvec.Phi(), tvec.M());
            if (gvec && glep1) {
              genLep    = new TLorentzVector(0,0,0,0);
              if(toolbox::flavor(genPartArr, BOSON_ID)*glepq1<0){
                genLep->SetPtEtaPhiM(glep1->Pt(),glep1->Eta(),glep1->Phi(),glep1->M());
                genNu->SetPtEtaPhiM(glep2->Pt(),glep2->Eta(),glep2->Phi(),glep2->M());
              }
              if(toolbox::flavor(genPartArr, BOSON_ID)*glepq2<0){
                genLep->SetPtEtaPhiM(glep2->Pt(),glep2->Eta(),glep2->Phi(),glep2->M());
                genNu->SetPtEtaPhiM(glep1->Pt(),glep1->Eta(),glep1->Phi(),glep1->M());
              }
            }

            delete gvec;
            delete glep1;
            delete glep2;
            gvec=0; glep1=0; glep2=0;
          }
          scale1fb     = weight;
          scale1fbUp   = weightUp;
          scale1fbDown = weightDown;

          met	     = info->pfMETC;
          metPhi    = info->pfMETCphi;
          puppiMet    = info->puppET;
          puppiMetPhi = info->puppETphi;
          q        = goodEle->q;
          lep       = &vLep;
          lep_raw   = &vLep_raw;

          ///// electron specific /////
          sc       = &vSC;
          pfCombIso = goodEle->chHadIso + TMath::Max(goodEle->neuHadIso + goodEle->gammaIso - 
          (info->rhoIso)*getEffAreaEl(goodEle->eta), 0.);

          outTree->Fill();
          delete genV; 
          delete genLep;
          delete genNu;
          genV=0, genLep=0, lep=0, genNu=0, sc=0;
          prefirePhoton=1; prefirePhotUp=1; prefirePhotDown=1;
          prefireJet   =1; prefireJetUp =1; prefireJetDown =1;
          prefireWeight=1; prefireUp    =1; prefireDown    =1;
        }
      }
      delete infile;
      infile=0, eventTree=0;    

      cout << nsel  << " +/- " << sqrt(nselvar);
      if(isam!=0) cout << " per 1/pb";
      cout << endl;
    }
    outFile->cd();
    hGenWeights->Write();
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
  delete electronArr;
  delete vertexArr;

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " W -> e nu" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  if(doScaleCorr)
  cout << "  *** Scale corrections applied ***" << endl;
  cout << endl;

  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  

  gBenchmark->Show("selectWe"); 
}
