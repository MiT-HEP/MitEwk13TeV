//================================================================================================
//
// Select Z->ee candidates
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

// Utilities belonging to this package
#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples
#include "../EleScale/EnergyScaleCorrection.h" //EGMSmear
#include "../Utils/LeptonIDCuts.hh"             // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"                  // various helper functions
#include "../Utils/PrefiringEfficiency.cc"      // prefiring efficiency functions

#include <TRandom3.h>
#endif

//=== MAIN MACRO ================================================================================================= 

void selectZee(const TString conf        ="zee.conf", // input file
               const TString outputDir   =".",   // output directory
               const Bool_t  doScaleCorr =0,    // apply energy scale corrections?
               const Int_t   sigma       =0,
               const Bool_t  doPU        =0,
               const Bool_t  is13TeV     =1,
               const Int_t   NSEC        =1,
               const Int_t   ITH         =0  ) {
  gBenchmark->Start("selectZee");
  cout << "Currently doing selecting samples with sqrt(s) = "  << ( is13TeV ? "13 TeV" : "5 TeV" ) << std::endl;
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  // probably don't need random
  TRandom3 *rand = new TRandom3();
  rand->SetSeed(1313131313);

  const Double_t MASS_LOW  = 40;
  const Double_t MASS_HIGH = 200;
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.4;
  const Double_t ELE_MASS  = 0.000511;
  
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 11;
  
  // load trigger menu
  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

  // Set up electron energy scale/smear
  const TString corrFiles = "../EleScale/Run2017_LowPU_v2";
  EnergyScaleCorrection ec( corrFiles.Data(), EnergyScaleCorrection::ECALELF);

  // Set up Prefiring Efficiencies
  const TString prefireFileName = "../Utils/All2017Gand2017HPrefiringMaps.root";
  PrefiringEfficiency pfire( prefireFileName.Data() , (is13TeV ? "2017H" : "2017G"));
  
  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/puWeights_76x.root", "read");

  // for systematics we need 3
  TH1D *h_rw = (TH1D*) f_rw->Get("puWeights");
  TH1D *h_rw_up = (TH1D*) f_rw->Get("puWeightsUp");
  TH1D *h_rw_down = (TH1D*) f_rw->Get("puWeightsDown");

  if (h_rw==NULL) cout<<"WARNIG h_rw == NULL"<<endl;
  if (h_rw_up==NULL) cout<<"WARNIG h_rw == NULL"<<endl;
  if (h_rw_down==NULL) cout<<"WARNIG h_rw == NULL"<<endl;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  enum { eEleEle2HLT=1, eEleEle1HLT1L1, eEleEle1HLT, eEleEleNoSel, eEleSC };  // event category enum
  
  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  //
  // parse .conf file
  //
  confParse(conf, snamev, samplev);
  const Bool_t hasData = (samplev[0]->fnamev.size()>0);

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  // const TString ntupDir = outputDir + TString("/ntuples");
  const TString ntupDir = outputDir + TString("/ntuples_") + Form("%d",ITH) + TString("_") + Form("%d",NSEC);
  gSystem->mkdir(ntupDir,kTRUE);
  
  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  TLorentzVector *genV=0;
  Float_t scale1fb,scale1fbUp,scale1fbDown;
  Float_t prefireWeight=1, prefireUp=1,    prefireDown=1;
  Float_t prefirePhoton=1, prefirePhotUp=1, prefirePhotDown=1;
  Float_t prefireJet=1,    prefireJetUp=1,  prefireJetDown=1;
  Float_t met, metPhi, u1, u2;
  Float_t puppiMet, puppiMetPhi, puppiU1, puppiU2;
  Int_t   q1, q2;
  Int_t   glepq1 = -99;
  Int_t   glepq2 = -99;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0, *lep1_raw=0, *lep2_raw=0;
  TLorentzVector *genlep1=0, *genlep2=0;
  
  ///// electron specific /////
  Float_t trkIso1, trkIso2, pfCombIso1, pfCombIso2;
  Float_t r91,r92;
  TLorentzVector *sc1=0, *sc2=0;
  Float_t lep1error, lep2error;
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr     = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *electronArr    = new TClonesArray("baconhep::TElectron");
  TClonesArray *scArr          = new TClonesArray("baconhep::TPhoton");
  TClonesArray *vertexArr      = new TClonesArray("baconhep::TVertex");
  TClonesArray *jetArr         = new TClonesArray("baconhep::TJet");

  TFile *infile=0;
  TTree *eventTree=0;
  
  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    
    // Assume data sample is first sample in .conf file
    // If sample is empty (i.e. contains no ntuple files), skip to next sample
    Bool_t isData=kFALSE;
    if(isam==0 && !hasData) continue;
    else if (isam==0) isData=kTRUE;
    
    Bool_t isSignal        = (snamev[isam].Contains("zee"));
    Bool_t isWboson        = (snamev[isam].Contains("wx"));
    Bool_t isWrongFlavor   = (snamev[isam].Contains("zxx"));
    Bool_t isRecoil = (isWboson||isSignal||isWrongFlavor);
    Bool_t noGen    = (snamev[isam].Contains("zz")||snamev[isam].Contains("wz")||snamev[isam].Contains("ww"));
    
    CSample* samp = samplev[isam];
  
    //
    // Set up output ntuple
    //
    TString outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    if(isam!=0 && !doScaleCorr) outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.raw.root");
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");
    outTree->Branch("runNum",     &runNum,     "runNum/i");      // event run number
    outTree->Branch("lumiSec",    &lumiSec,    "lumiSec/i");     // event lumi section
    outTree->Branch("evtNum",     &evtNum,     "evtNum/i");      // event number
    outTree->Branch("matchGen",   &matchGen,   "matchGen/i");    // event has both leptons matched to MC Z->ll
    outTree->Branch("category",   &category,   "category/i");    // dilepton category
    outTree->Branch("npv",        &npv,        "npv/i");         // number of primary vertices
    outTree->Branch("npu",        &npu,        "npu/i");         // number of in-time PU events (MC)
    outTree->Branch("genV",      "TLorentzVector",  &genV);      // GEN boson 4-vector
    outTree->Branch("prefireWeight", &prefireWeight, "prefireWeight/F");
    outTree->Branch("prefireUp",     &prefireUp,     "prefireUp/F");
    outTree->Branch("prefireDown",   &prefireDown,   "prefireDown/F");
    outTree->Branch("prefirePhoton", &prefirePhoton, "prefirePhoton/F");
    outTree->Branch("prefirePhotUp",     &prefirePhotUp,     "prefirePhotUp/F");
    outTree->Branch("prefirePhotDown",   &prefirePhotDown,   "prefirePhotDown/F");
    outTree->Branch("prefireJet",    &prefireJet,    "prefireJet/F");
    outTree->Branch("prefireJetUp",  &prefireJetUp,  "prefireJetUp/F");
    outTree->Branch("prefireJetDown",&prefireJetDown,"prefireJetDown/F");
    outTree->Branch("scale1fb",      &scale1fb,   "scale1fb/F");    // event weight per 1/fb (MC)
    outTree->Branch("scale1fbUp",    &scale1fbUp,   "scale1fbUp/F");    // event weight per 1/fb (MC)
    outTree->Branch("scale1fbDown",  &scale1fbDown,   "scale1fbDown/F");    // event weight per 1/fb (MC)
    outTree->Branch("met",        &met,        "met/F");         // MET
    outTree->Branch("metPhi",     &metPhi,     "metPhi/F");      // phi(MET)
    outTree->Branch("u1",         &u1,         "u1/F");          // parallel component of recoil
    outTree->Branch("u2",         &u2,         "u2/F");          // perpendicular component of recoil
    outTree->Branch("puppiMet",    &puppiMet,   "puppiMet/F");      // Puppi MET
    outTree->Branch("puppiMetPhi", &puppiMetPhi,"puppiMetPhi/F");   // phi(Puppi MET)
    outTree->Branch("puppiU1",     &puppiU1,    "puppiU1/F");       // parallel component of recoil (Puppi MET)
    outTree->Branch("puppiU2",     &puppiU2,    "puppiU2/F");       // perpendicular component of recoil (Puppi MET)
    outTree->Branch("q1",         &q1,         "q1/I");          // charge of tag lepton
    outTree->Branch("q2",         &q2,         "q2/I");          // charge of probe lepton
    outTree->Branch("glepq1",         &glepq1,         "glepq1/I");          // charge of tag lepton
    outTree->Branch("glepq2",         &glepq2,         "glepq2/I");          // charge of probe lepton
    outTree->Branch("dilep",         "TLorentzVector",  &dilep);    // di-lepton 4-vector
    outTree->Branch("lep1",          "TLorentzVector",  &lep1);     // tag lepton 4-vector
    outTree->Branch("lep2",          "TLorentzVector",  &lep2);     // probe lepton 4-vector
    outTree->Branch("genlep1",       "TLorentzVector",  &genlep1);     // tag lepton 4-vector
    outTree->Branch("genlep2",       "TLorentzVector",  &genlep2);     // probe lepton 4-vector
    outTree->Branch("lep1_raw",      "TLorentzVector",  &lep1_raw);     // tag lepton 4-vector
    outTree->Branch("lep2_raw",      "TLorentzVector",  &lep2_raw);     // probe lepton 4-vector
    ///// electron specific /////
    outTree->Branch("trkIso1",    &trkIso1,    "trkIso1/F");     // track isolation of tag lepton
    outTree->Branch("trkIso2",    &trkIso2,    "trkIso2/F");     // track isolation of probe lepton
    outTree->Branch("pfCombIso1", &pfCombIso1, "pfCombIso1/F");  // PF combine isolation of tag lepton
    outTree->Branch("pfCombIso2", &pfCombIso2, "pfCombIso2/F");  // PF combined isolation of probe lepton    
    outTree->Branch("sc1",        "TLorentzVector",  &sc1);       // tag supercluster 4-vector
    outTree->Branch("sc2",        "TLorentzVector",  &sc2);       // probe supercluster 4-vector
    outTree->Branch("r91",        &r91,        "r91/F");	 // transverse impact parameter of tag
    outTree->Branch("r92",        &r92,        "r92/F");	 // transverse impact parameter of probe	  
    outTree->Branch("lep1error",  &lep1error,  "lep1error/F");   // scale and smear correction uncertainty for tag lepton
    outTree->Branch("lep2error",  &lep2error,  "lep2error/F");   // scale and smear correction uncertainty for probe leptom

    TH1D* hGenWeights = new TH1D("hGenWeights","hGenWeights",10,-10.,10.);
    //
    // loop through files
    //
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {

      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... " << endl; cout.flush();
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
      eventTree->SetBranchAddress("Info",     &info       ); TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Photon",   &scArr      ); TBranch *scBr       = eventTree->GetBranch("Photon");
      eventTree->SetBranchAddress("PV",       &vertexArr  ); TBranch *vertexBr   = eventTree->GetBranch("PV");
      if(hasJet) eventTree->SetBranchAddress("AK4",&jetArr); TBranch *jetBr      = eventTree->GetBranch("AK4");
      Bool_t hasGen = (eventTree->GetBranchStatus("GenEvtInfo")&&!noGen);
      TBranch *genBr=0, *genPartBr=0;
      if(hasGen) {
        eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
        eventTree->SetBranchAddress("GenParticle",&genPartArr); genPartBr = eventTree->GetBranch("GenParticle");
      }

      // Compute MC event weight per 1/fb
      const Double_t xsec = samp->xsecv[ifile];
      Double_t puWeight=0, puWeightUp=0, puWeightDown=0;
      
      TLorentzVector *gvec=new TLorentzVector(0,0,0,0);
      TLorentzVector *glep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *glep2=new TLorentzVector(0,0,0,0);
      TLorentzVector *glep3 =new TLorentzVector(0,0,0,0);
      TLorentzVector *glep4 =new TLorentzVector(0,0,0,0);
      TLorentzVector *gph=new TLorentzVector(0,0,0,0);

      
      //
      // loop over events
      //
      // cout << "n sections " << NSEC << endl;
      double frac = 1.0/NSEC;
      // cout << "n sections " << NSEC << "  frac " << frac << endl;
      UInt_t IBEGIN = frac*ITH*eventTree->GetEntries();
      UInt_t IEND = frac*(ITH+1)*eventTree->GetEntries();
      // cout << "start, end " << IBEGIN << " " << IEND << endl;
      // UInt_t NTEST = (UInt_t)(eventTree->GetEntries()*0.001);
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
          puWeight = doPU ? h_rw->GetBinContent(h_rw->FindBin(info->nPUmean)) : 1.;
          puWeightUp = doPU ? h_rw_up->GetBinContent(h_rw_up->FindBin(info->nPUmean)) : 1.;
          puWeightDown = doPU ? h_rw_down->GetBinContent(h_rw_down->FindBin(info->nPUmean)) : 1.;
          hGenWeights->Fill(0.0,gen->weight);
          weight*=gen->weight*puWeight;
          weightUp*=gen->weight*puWeightUp;
          weightDown*=gen->weight*puWeightDown;
        } else {
          hGenWeights->Fill(0.0,1.0);
        }
	
        // veto z -> xx decays for signal and z -> ee for bacground samples (needed for inclusive DYToLL sample)
        if (isWrongFlavor && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;
        else if (isSignal && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;

        // check for certified lumi (if applicable)
        baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        if(hasJSON && !rlrm.hasRunLumi(rl)) continue;  

        // trigger requirement
        if (!isEleTrigger(triggerMenu, info->triggerBits, isData, is13TeV)) continue;

        // good vertex requirement
        if(!(info->hasGoodPV)) continue;
        
        electronArr->Clear();
        electronBr->GetEntry(ientry);
        scArr->Clear();
        scBr->GetEntry(ientry);
        jetArr->Clear();
        if(hasJet) jetBr->GetEntry(ientry);

        TLorentzVector vTag(0,0,0,0);
        TLorentzVector vTag_raw(0,0,0,0);
        TLorentzVector vTagfinal(0,0,0,0);
        TLorentzVector vTagSCfinal(0,0,0,0);
        TLorentzVector vTagSC(0,0,0,0);
        Double_t tagPt=0;
        Double_t Pt1=0;
        Double_t Pt2=0;
        Int_t itag=-1;
        Int_t tagscID=-1;
        double tagRandom = rand->Gaus(0,1);
        double probeRandom = rand->Gaus(0,1);
        double eleProbeRandom = rand->Gaus(0,1);
        double eleProbeSCRandom = rand->Gaus(0,1);
        // random = tagRandom;
        for(Int_t i1=0; i1<electronArr->GetEntriesFast(); i1++) {
          const baconhep::TElectron *tag = (baconhep::TElectron*)((*electronArr)[i1]);
          double tagEcalE = tag->ecalEnergy;
          double eTregress = tagEcalE/cosh(fabs(tag->eta));
          vTag.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, ELE_MASS);
          vTagSC.SetPtEtaPhiM(tag->scEt, tag->scEta, tag->scPhi, ELE_MASS);
   
          // if(fabs(vTag.Eta())>=ECAL_GAP_LOW && fabs(vTag.Eta())<=ECAL_GAP_HIGH) continue;
            
          float tagError = 0.;
          float tagSCError = 0.;
          if(doScaleCorr && (tag->r9 < 1.)){
            // set up variable and apply scale and smear correction to tag
            // only apply correction to electron pt, not SC pt, since we never use tag SC info.
            float tagSmear = 0.;
            float tagScale = 1.;
            float tagSCSmear = 0.;
            float tagSCScale = 1.;

            float tagAbsEta   = fabs(vTag.Eta());
            float tagEt       = vTag.E() / cosh(tagAbsEta);
            bool  tagisBarrel = tagAbsEta < 1.4442;
            float tagSCAbsEta   = fabs(vTagSC.Eta());
            float tagSCEt       = vTagSC.E() / cosh(tagSCAbsEta);

            if(snamev[isam].Contains("data")){//Data
              int runNumber = is13TeV ? info->runNum : 306936 ;
              tagScale   = ec.scaleCorr(runNumber, eTregress, tagAbsEta  , tag->r9);
              tagSCScale = ec.scaleCorr(runNumber, tagSCEt  , tagSCAbsEta, tag->r9);
              
              tagError   = ec.scaleCorrUncert(runNumber, eTregress, tagAbsEta  , tag->r9,12,1);
              tagSCError = ec.scaleCorrUncert(runNumber, tagSCEt  , tagSCAbsEta, tag->r9,12,1);
              
              (vTag)*=tagScale*(1+sigma*tagError);
              (vTagSC)*=tagSCScale*(1+sigma*tagSCError);

            } else {//MC
              tagSmear = ec.smearingSigma(info->runNum, eTregress, tagAbsEta, tag->r9, 12, sigma, 0.);
              tagSCSmear = ec.smearingSigma(info->runNum, tagSCEt, tagSCAbsEta, tag->r9, 12, sigma, 0.);           
              (vTag) *= (1. + tagSmear*tagRandom);
              (vTagSC) *= 1. + tagSCSmear * tagRandom;
              double tagSmearEP = ec.smearingSigma(info->runNum, eTregress, tagAbsEta, tag->r9, 12,  1., 0.);
              double tagSmearEM = ec.smearingSigma(info->runNum, eTregress, tagAbsEta, tag->r9, 12, -1., 0.);
              tagError = tagRandom * std::hypot(tagSmearEP - tagSmear, tagSmearEM - tagSmear); 
            } 
          }
      
          // Check Kinematic Cuts & Electron ID
          if(vTag.Pt()	         < PT_CUT)                  continue;  // lepton pT cut
          if(fabs(vTag.Eta())    > ETA_CUT)                 continue;  // lepton |eta| cut
          if(!passEleMediumID(tag, vTag, info->rhoIso))     continue;  // lepton selection
          

          double El_Pt=0;
          El_Pt = vTag.Pt();

          if( El_Pt > Pt1 ) {
	    Pt2=Pt1;
	    Pt1=El_Pt;
	  } else if ( El_Pt > Pt2 && El_Pt < Pt1 ) {
	    Pt2=El_Pt;
	  }

          if(!isEleTriggerObj(triggerMenu, tag->hltMatchBits, kFALSE, isData, is13TeV)) continue;
          if(El_Pt<tagPt) continue;
      
          tagPt=El_Pt;
          itag=i1;
          tagscID=tag->scID;
          vTagfinal = vTag;
          vTagSCfinal = vTagSC;
          vTag_raw.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, ELE_MASS);

          trkIso1    = tag->trkIso;
          pfCombIso1 = tag->chHadIso + TMath::Max(tag->neuHadIso + tag->gammaIso - (info->rhoIso)*getEffAreaEl(vTag.Eta()), 0.);
          q1         = tag->q;
          r91        = tag->r9;
          lep1error  = tagError;

        }
  
        if(tagPt<Pt2) continue;

        TLorentzVector vProbe(0,0,0,0); TLorentzVector vProbeSC(0,0,0,0);
        TLorentzVector vProbefinal(0,0,0,0), vProbe_raw(0,0,0,0);
        float probeErrorfinal;
        float probeSCErrorfinal; 
        Double_t probePt=0;
        Int_t iprobe=-1;
        Int_t passID=false;
        UInt_t icat=0;
	
        const baconhep::TElectron *eleProbe=0;

        for(Int_t j=0; j<scArr->GetEntriesFast(); j++) {
          const baconhep::TPhoton *scProbe = (baconhep::TPhoton*)((*scArr)[j]);
          vProbe.SetPtEtaPhiM(scProbe->pt, scProbe->eta, scProbe->phi, ELE_MASS);
          if(scProbe->scID == tagscID) continue;

          // check ECAL gap
          // if(fabs(vProbe.Eta())>=ECAL_GAP_LOW && fabs(vProbe.Eta())<=ECAL_GAP_HIGH) continue;

          float probeError = 0.;
          if(doScaleCorr && (scProbe->r9 < 1.)){
            // set up variable and apply scale and smear correction to probe
            float probeSmear = 0.;
            float probeScale = 1.;

            float probeAbsEta   = fabs(vProbe.Eta());
            float probeEt       = vProbe.E() / cosh(probeAbsEta);
            bool  probeisBarrel = probeAbsEta < 1.4442;
            
            if(snamev[isam].Contains("data")){//Data
              int runNumber = is13TeV ? info->runNum : 306936 ;
              probeScale   = ec.scaleCorr(runNumber, probeEt, probeAbsEta  , scProbe->r9);
              probeError   = ec.scaleCorrUncert(runNumber, probeEt, probeAbsEta  , scProbe->r9,12,1);
              (vProbe) *= probeScale * (1 + sigma*probeError);
            } else {//MC
              probeSmear = ec.smearingSigma(info->runNum, probeEt, probeAbsEta, scProbe->r9, 12, sigma, 0.);

              (vProbe) *= (1. + probeSmear*probeRandom);
              double probeSmearEP = ec.smearingSigma(info->runNum, probeEt, probeAbsEta, scProbe->r9, 12,  1., 0.);
              double probeSmearEM = ec.smearingSigma(info->runNum, probeEt, probeAbsEta, scProbe->r9, 12, -1., 0.);
              probeError = probeRandom * std::hypot(probeSmearEP - probeSmear, probeSmearEM - probeSmear);
            }
          }
          double probeEcalEnergy_tmp = 0;
          if(fabs(vProbe.Eta())  > ETA_CUT) continue;
          for(Int_t i2=0; i2<electronArr->GetEntriesFast(); i2++) {
            if(itag==i2) continue;
            const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[i2]);
            if(scProbe->scID==ele->scID) {
              eleProbe = ele; 
              iprobe   = i2;
              break; 
            }
          }
          double El_Pt=0;
          TLorentzVector vEleProbe(0,0,0,0), vEleProbeSC(0,0,0,0);
          if(eleProbe){
            double probeEcalE = eleProbe->ecalEnergy;
            double eTregress = eleProbe->ecalEnergy/cosh(fabs(eleProbe->eta));
            vEleProbe.SetPtEtaPhiM(eleProbe->pt, eleProbe->eta, eleProbe->phi, ELE_MASS);
            vEleProbeSC.SetPtEtaPhiM(eleProbe->scEt, eleProbe->scEta, eleProbe->scPhi, ELE_MASS);
            // if(fabs(vEleProbe.Eta())>=ECAL_GAP_LOW && fabs(vEleProbe.Eta())<=ECAL_GAP_HIGH) continue;

            float eleProbeError = 0., eleProbeSCError = 0.;
            
            if(doScaleCorr && (eleProbe->r9 < 1.)){
              float eleProbeSmear   = 0., eleProbeScale   = 1.;
              float eleProbeSCSmear = 0., eleProbeSCScale = 1.;

              float eleProbeAbsEta   = fabs(vEleProbe.Eta());
              float eleProbeEt       = vEleProbe.E() / cosh(eleProbeAbsEta);

              float eleProbeSCAbsEta   = fabs(vEleProbeSC.Eta());
              float eleProbeSCEt       = vEleProbeSC.E() / cosh(eleProbeSCAbsEta);
  
              if(snamev[isam].Contains("data")){//Data
                int runNumber = is13TeV ? info->runNum : 306936 ;
                eleProbeScale   = ec.scaleCorr(runNumber, eTregress   , eleProbeAbsEta  , eleProbe->r9);
                eleProbeSCScale = ec.scaleCorr(runNumber, eleProbeSCEt, eleProbeSCAbsEta, eleProbe->r9);
                
                eleProbeError   = ec.scaleCorrUncert(runNumber, eTregress   , eleProbeAbsEta  , eleProbe->r9, 12, 1);
                eleProbeSCError = ec.scaleCorrUncert(runNumber, eleProbeSCEt, eleProbeSCAbsEta, eleProbe->r9, 12, 1);
                
                (vEleProbe)   *= eleProbeScale   * (1 + sigma*eleProbeError);
                (vEleProbeSC) *= eleProbeSCScale * (1 + sigma*eleProbeSCError);
  
              } else {//MC

		float eleProbeR9Prime = eleProbe->r9; // no r9 after 2016

		eleProbeSmear   = ec.smearingSigma(info->runNum, eTregress   , eleProbeAbsEta  , eleProbeR9Prime, 12, sigma, 0.);
		eleProbeSCSmear = ec.smearingSigma(info->runNum, eleProbeSCEt, eleProbeSCAbsEta, eleProbeR9Prime, 12, sigma, 0.);
              
		float eleProbeSmearEP = ec.smearingSigma(info->runNum, eTregress, eleProbeAbsEta, eleProbeR9Prime, 12, 1., 0.);
		float eleProbeSmearEM = ec.smearingSigma(info->runNum, eTregress, eleProbeAbsEta, eleProbeR9Prime, 12, -1., 0.);

		float eleProbeSCSmearEP = ec.smearingSigma(info->runNum, eleProbeSCEt, eleProbeSCAbsEta, eleProbeR9Prime, 12, 1., 0.);
		float eleProbeSCSmearEM = ec.smearingSigma(info->runNum, eleProbeSCEt, eleProbeSCAbsEta, eleProbeR9Prime, 12, -1., 0.);

		(vEleProbe) *= (1.+ eleProbeSmear*eleProbeRandom);
		(vEleProbeSC) *= 1. + eleProbeSCSmear * eleProbeSCRandom;

		eleProbeError = eleProbeRandom * std::hypot(eleProbeSmearEP - eleProbeSmear, eleProbeSmearEM - eleProbeSmear);
		eleProbeSCError = eleProbeSCRandom * std::hypot(eleProbeSCSmearEP - eleProbeSCSmear, eleProbeSCSmearEM - eleProbeSCSmear);
	      }
	    }
	    probeEcalEnergy_tmp = probeEcalE;
	    El_Pt = vEleProbe.Pt();
	    probeErrorfinal = eleProbeError;
	    probeSCErrorfinal = eleProbeSCError;
	  }else{
	    El_Pt = vProbe.Pt();
	    // if(fabs(vProbe.Eta())>=ECAL_GAP_LOW && fabs(vProbe.Eta())<=ECAL_GAP_HIGH) continue;
	    probeErrorfinal = probeError;
	    probeSCErrorfinal = probeError;
	  }
	  if(El_Pt < PT_CUT) continue;
	  if(passID&&eleProbe&&passEleMediumID(eleProbe,vEleProbe,info->rhoIso)&&El_Pt<probePt) continue;
	  if(passID&&eleProbe&&!passEleMediumID(eleProbe,vEleProbe,info->rhoIso)) continue;
	  if(passID&&!eleProbe) continue;
	  if(!passID&&eleProbe&&!passEleMediumID(eleProbe,vEleProbe,info->rhoIso)&&El_Pt<probePt) continue;
	  if(!passID&&!eleProbe&&El_Pt<probePt) continue;
	  if(!passID&&eleProbe&&passEleMediumID(eleProbe,vEleProbe,info->rhoIso)) passID=true;

	  probePt=El_Pt;
	  vProbefinal = (eleProbe) ?  vEleProbe : vProbe ;
	  if(eleProbe) vProbe_raw.SetPtEtaPhiM(eleProbe->pt, eleProbe->eta, eleProbe->phi, ELE_MASS);
	  vProbeSC = (eleProbe) ? vEleProbeSC : vProbe ;
	  // lep2EcalE = probeEcalEnergy_tmp;

	  trkIso2    = (eleProbe) ? eleProbe->trkIso        : -1;
	  pfCombIso2 = (eleProbe) ? 
	    eleProbe->chHadIso + TMath::Max(eleProbe->neuHadIso + eleProbe->gammaIso - 
					    (info->rhoIso)*getEffAreaEl(vEleProbe.Eta()), 0.) :  -1;
	  q2         = (eleProbe) ? eleProbe->q : -q1;
	  lep2error  = probeErrorfinal;
	  // sc2error   = probeSCErrorfinal;


	  // determine event category
	  if(eleProbe) {
	    if(passEleMediumID(eleProbe,vEleProbe,info->rhoIso)) {
	      if(isEleTriggerObj(triggerMenu, eleProbe->hltMatchBits, kFALSE, isData, is13TeV)) {
		icat=eEleEle2HLT;
	      } else if(isEleTriggerObj(triggerMenu, eleProbe->hltMatchBits, kFALSE, isData, is13TeV)) {
		icat=eEleEle1HLT1L1; // does this ever get used
	      } else { icat=eEleEle1HLT; }
	    } else { icat=eEleEleNoSel; }
	  } else { icat=eEleSC; }
	}

	// if(q1 == q2)         continue;  // opposite charge requirement
	// mass window
	TLorentzVector vDilep = vTagfinal + vProbefinal;
	// TLorentzVector vDilepSC = vTagSCfinal + vProbeSC;
	if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
	if(icat==0) continue;
      
      
	// do the prefiring weights
	if(!isData){
	  pfire.setObjects(scArr,jetArr);
	  pfire.computePhotonsOnly(prefirePhoton, prefirePhotUp, prefirePhotDown);
	  pfire.computeJetsOnly   (prefireJet   , prefireJetUp , prefireJetDown );
	  pfire.computeFullPrefire(prefireWeight, prefireUp    , prefireDown    );
	}

	//******** We have a Z candidate! HURRAY! ********
	nsel+=isData ? 1 : weight;
	nselvar+=isData ? 1 : weight*weight;

	// Perform matching of dileptons to GEN leptons from Z decay
	Bool_t hasGenMatch = kFALSE;
	if(isRecoil && hasGen) {
	  toolbox::fillGenBorn(genPartArr, BOSON_ID, gvec, glep3, glep4, glep1, glep2);
        
	  Bool_t match1 = ( ((glep1) && toolbox::deltaR(vTagfinal.Eta(), vTagfinal.Phi(), glep1->Eta(), glep1->Phi())<0.3) || 
			    ((glep2) && toolbox::deltaR(vTagfinal.Eta(), vTagfinal.Phi(), glep2->Eta(), glep2->Phi())<0.3) );
        
	  Bool_t match2 = ( ((glep1) && toolbox::deltaR(vProbefinal.Eta(), vProbefinal.Phi(), glep1->Eta(), glep1->Phi())<0.3) || 
			    ((glep2) && toolbox::deltaR(vProbefinal.Eta(), vProbefinal.Phi(), glep2->Eta(), glep2->Phi())<0.3) );
	  TLorentzVector tvec=*glep1+*glep2;
	  genV=new TLorentzVector(0,0,0,0);
	  genV->SetPtEtaPhiM(tvec.Pt(), tvec.Eta(), tvec.Phi(), tvec.M());
	  genlep1=new TLorentzVector(0,0,0,0);
	  genlep2=new TLorentzVector(0,0,0,0);
	  genlep1->SetPtEtaPhiM(glep1->Pt(),glep1->Eta(),glep1->Phi(),glep1->M());
	  genlep2->SetPtEtaPhiM(glep2->Pt(),glep2->Eta(),glep2->Phi(),glep2->M());
        
	  if(match1 && match2) hasGenMatch = kTRUE;
	}
      
	//
	// Fill tree
	//
	runNum   = info->runNum;
	lumiSec  = info->lumiSec;
	evtNum   = info->evtNum;

	if (hasGenMatch) matchGen=1;
	else matchGen=0;


	category = icat;

	vertexArr->Clear();
	vertexBr->GetEntry(ientry);

	npv      = vertexArr->GetEntries();
	npu      = info->nPUmean;
	scale1fb = weight;
	scale1fbUp   = weightUp;
	scale1fbDown = weightDown;
	met       = info->pfMETC;
	metPhi    = info->pfMETCphi;
      
	puppiMet    = info->puppET;
	puppiMetPhi = info->puppETphi;
	// puppiSumEt = 0;
	lep1       = &vTagfinal;
	lep2       = &vProbefinal;
	lep1_raw   = &vTag_raw;
	lep2_raw   = &vProbe_raw;

	dilep      = &vDilep;
	sc1        = &vTagSCfinal;
	sc2        = &vProbeSC;
	TVector2 vZPt((vDilep.Pt())*cos(vDilep.Phi()),(vDilep.Pt())*sin(vDilep.Phi()));
	TVector2 vMet((info->pfMETC)*cos(info->pfMETCphi), (info->pfMETC)*sin(info->pfMETCphi));
	TVector2 vU = -1.0*(vMet+vZPt);
	u1 = ((vDilep.Px())*(vU.Px()) + (vDilep.Py())*(vU.Py()))/(vDilep.Pt());  // u1 = (pT . u)/|pT|
	u2 = ((vDilep.Px())*(vU.Py()) - (vDilep.Py())*(vU.Px()))/(vDilep.Pt());  // u2 = (pT x u)/|peleProbe	
      
	TVector2 vPuppiMet((info->puppET)*cos(info->puppETphi), (info->puppET)*sin(info->puppETphi));
	TVector2 vPuppiU = -1.0*(vPuppiMet+vZPt);
	puppiU1 = ((vDilep.Px())*(vPuppiU.Px()) + (vDilep.Py())*(vPuppiU.Py()))/(vDilep.Pt());  // u1 = (pT . u)/|pT|
	puppiU2 = ((vDilep.Px())*(vPuppiU.Py()) - (vDilep.Py())*(vPuppiU.Px()))/(vDilep.Pt());  // u2 = (pT x u)/|pT|
	outTree->Fill();
      
	////   -------         RESET EVERYTHING     -------------
	delete genV;
	delete genlep1;
	delete genlep2;
	genV=0, dilep=0, lep1=0, lep2=0, sc1=0, sc2=0, lep1_raw=0, lep2_raw=0, genlep1=0, genlep2=0;
	prefirePhoton=1; prefirePhotUp=1; prefirePhotDown=1;
	prefireJet   =1; prefireJetUp =1; prefireJetDown =1;
	prefireWeight=1; prefireUp    =1; prefireDown    =1;
      }
      delete infile;
      infile=0, eventTree=0;    
      
      cout << nsel  << " +/- " << sqrt(nselvar);
      if(!isData) cout << " per 1/fb";
      cout << endl;
    }
    outFile->cd();
    hGenWeights->Write();
    outFile->Write();
    outFile->Close(); 
  }
  delete h_rw;
  delete f_rw;
  delete info;
  delete gen;
  delete genPartArr;
  delete electronArr;
  delete scArr;
  delete vertexArr;
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " Z -> e e" << endl;
  cout << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectZee"); 
}
