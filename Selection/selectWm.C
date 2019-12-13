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
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

// lumi section selection with JSON files
#include "BaconAna/Utils/interface/RunLumiRangeMap.hh"
#include "CCorrUser2D.hh"

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions
#endif

//=== MAIN MACRO ================================================================================================= 

void selectWm(const TString conf="wm.conf", // input file
              const TString outputDir=".",  // output directory
	          const Bool_t  doScaleCorr=0,   // apply energy scale corrections?
              const Bool_t  doPU=0,
              const Bool_t is13TeV=1,
                const Int_t NSEC = 1,
                const Int_t ITH = 0
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
  
  const Int_t NPDF = 100;
  const Int_t NQCD = 6;

  // load trigger menu
  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

  const TString prefireFileName = "../Utils/All2017Gand2017HPrefiringMaps.root";
  TFile *prefireFile = new TFile(prefireFileName);
  CCorrUser2D prefirePhotonCorr, prefireJetCorr;
  if(!is13TeV){
    prefirePhotonCorr.loadCorr((TH2D*)prefireFile->Get("L1prefiring_photonpt_2017G")); // 5 TeV photon prefire
    prefireJetCorr.loadCorr((TH2D*)prefireFile->Get("L1prefiring_jetpt_2017G")); // 5 TeV jet prefire
  } else if(is13TeV){
    prefirePhotonCorr.loadCorr((TH2D*)prefireFile->Get("L1prefiring_photonpt_2017H")); // 13 TeV photon prefire
    prefireJetCorr.loadCorr((TH2D*)prefireFile->Get("L1prefiring_jetpt_2017H")); // 13 TeV jet prefire
  }
  
  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/puWeights_76x.root", "read");

  TH1D *h_rw = (TH1D*) f_rw->Get("puWeights");
  TH1D *h_rw_up = (TH1D*) f_rw->Get("puWeightsUp");
  TH1D *h_rw_down = (TH1D*) f_rw->Get("puWeightsDown");

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
  // const TString ntupDir = outputDir + TString("/ntuples");
  const TString ntupDir = outputDir + TString("/ntuples_") + Form("%d",ITH) + TString("_") + Form("%d",NSEC);
  gSystem->mkdir(ntupDir,kTRUE);
  
  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  UInt_t  id_1, id_2;
  Double_t x_1, x_2, xPDF_1, xPDF_2;
  Double_t scalePDF, weightPDF;
  TLorentzVector *genV=0, *genLep=0;//, *genMuonMatch=0;
  Float_t genMuonPt;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown;
  Float_t prefireWeight=1, prefireUp=1,    prefireDown=1;
  Float_t prefirePhoton=1, prefirePhotUp=1, prefirePhotDown=1;
  Float_t prefireJet=1,    prefireJetUp=1,  prefireJetDown=1;
  // Float_t prefireWeight, prefireUp, prefireDown;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Float_t metDJee, metPhiDJee, sumEtDJee, u1DJee, u2DJee;
  Float_t tkMet, tkMetPhi, tkSumEt, tkMt, tkU1, tkU2;
  Float_t mvaMet, mvaMetPhi, mvaSumEt, mvaMt, mvaU1, mvaU2;
  Float_t puppiMet, puppiMetPhi, puppiSumEt, puppiMt, puppiU1, puppiU2;
  Int_t   q;
  TLorentzVector *lep=0;
  Int_t lepID;
  vector<Double_t> lheweight(NPDF+NQCD,0);
  ///// muon specific /////
  Float_t trkIso, emIso, hadIso;
  Float_t pfChIso, pfGamIso, pfNeuIso, pfCombIso;
  Float_t d0, dz;
  Float_t muNchi2;
  UInt_t nPixHits, nTkLayers, nValidHits, nMatch, typeBits;
  // Bool_t passHLT;

  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info  = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen  = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr    = new TClonesArray("baconhep::TMuon");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  TClonesArray *scArr          = new TClonesArray("baconhep::TPhoton");
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

    // std::cout << snamev[isam].CompareTo("wm",TString::kIgnoreCase) << std::endl;
    // std::cout << snamev[isam].CompareTo("wmn",TString::kIgnoreCase) << std::endl;
    // std::cout << snamev[isam].CompareTo("wx",TString::kIgnoreCase) << std::endl;
    // std::cout << snamev[isam].CompareTo("wx0",TString::kIgnoreCase) << std::endl;
    // std::cout << snamev[isam].CompareTo("wxn",TString::kIgnoreCase) << std::endl;
    // Assume signal sample is given name "wm" -- flag to store GEN W kinematics
    Bool_t isSignal=false;
    Bool_t isWrongFlavor=false;
    if(is13TeV){
      isSignal = (snamev[isam].CompareTo("wm0",TString::kIgnoreCase)==0||snamev[isam].CompareTo("wm1",TString::kIgnoreCase)==0||snamev[isam].CompareTo("wm2",TString::kIgnoreCase)==0);
      isWrongFlavor = (snamev[isam].CompareTo("wx0",TString::kIgnoreCase)==0||snamev[isam].CompareTo("wx1",TString::kIgnoreCase)==0||snamev[isam].CompareTo("wx2",TString::kIgnoreCase)==0);
    } else {
      isSignal = (snamev[isam].CompareTo("wm",TString::kIgnoreCase)==0);
      isWrongFlavor = (snamev[isam].CompareTo("wx",TString::kIgnoreCase)==0);
    }
    // flag to reject W->mnu events when selecting wrong flavor background events
    // Bool_t isWrongFlavor = (snamev[isam].CompareTo("wx0",TString::kIgnoreCase)==0||snamev[isam].CompareTo("wx1",TString::kIgnoreCase)==0||snamev[isam].CompareTo("wx2",TString::kIgnoreCase)==0);
    //flag to save the info for recoil corrections
    Bool_t isRecoil = (isSignal||(snamev[isam].CompareTo("zxx",TString::kIgnoreCase)==0)||isWrongFlavor);
    Bool_t noGen = (snamev[isam].CompareTo("zz",TString::kIgnoreCase)==0||snamev[isam].CompareTo("wz",TString::kIgnoreCase)==0||snamev[isam].CompareTo("ww",TString::kIgnoreCase)==0);
    
    CSample* samp = samplev[isam];

    //
    // Set up output ntuple
    //
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
    outTree->Branch("id_1",       &id_1,       "id_1/i");        // PDF info -- parton ID for parton 1
    outTree->Branch("id_2",       &id_2,       "id_2/i");        // PDF info -- parton ID for parton 2
    outTree->Branch("x_1",        &x_1,        "x_1/d");         // PDF info -- x for parton 1
    outTree->Branch("x_2",        &x_2,        "x_2/d");         // PDF info -- x for parton 2
    outTree->Branch("xPDF_1",     &xPDF_1,     "xPDF_1/d");      // PDF info -- x*F for parton 1
    outTree->Branch("xPDF_2",     &xPDF_2,     "xPDF_2/d");      // PDF info -- x*F for parton 2
    outTree->Branch("scalePDF",   &scalePDF,   "scalePDF/d");    // PDF info -- energy scale of parton interaction
    outTree->Branch("weightPDF",  &weightPDF,  "weightPDF/d");   // PDF info -- PDF weight
    outTree->Branch("genV",       "TLorentzVector", &genV);      // GEN boson 4-vector (signal MC)
    outTree->Branch("genLep",     "TLorentzVector", &genLep);    // GEN lepton 4-vector (signal MC)
    // outTree->Branch("genMuonMatch","TLorentzVector", &genMuonMatch);    // GEN lepton 4-vector (signal MC)
    outTree->Branch("genMuonPt",  &genMuonPt,     "genMuonPt/F");      // GEN boson pT (signal MC)
    outTree->Branch("genVPt",     &genVPt,     "genVPt/F");      // GEN boson pT (signal MC)
    outTree->Branch("genVPhi",    &genVPhi,    "genVPhi/F");     // GEN boson phi (signal MC)
    outTree->Branch("genVy",      &genVy,      "genVy/F");       // GEN boson rapidity (signal MC)
    outTree->Branch("genVMass",   &genVMass,   "genVMass/F");    // GEN boson mass (signal MC)
    outTree->Branch("genLepPt",   &genLepPt,   "genLepPt/F");    // GEN lepton pT (signal MC)
    outTree->Branch("genLepPhi",  &genLepPhi,  "genLepPhi/F");   // GEN lepton phi (signal MC)
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
    outTree->Branch("sumEt",      &sumEt,      "sumEt/F");       // Sum ET
    outTree->Branch("mt",         &mt,         "mt/F");          // transverse mass
    outTree->Branch("u1",         &u1,         "u1/F");          // parallel component of recoil
    outTree->Branch("u2",         &u2,         "u2/F");          // perpendicular component of recoil
    outTree->Branch("metDJee",        &metDJee,        "metDJee/F");         // MET
    outTree->Branch("metPhiDJee",     &metPhiDJee,     "metPhiDJee/F");      // phi(MET)
    outTree->Branch("sumEtDJee",      &sumEtDJee,      "sumEtDJee/F");       // Sum ET
    outTree->Branch("u1DJee",         &u1DJee,         "u1DJee/F");          // parallel component of recoil
    outTree->Branch("u2DJee",         &u2DJee,         "u2DJee/F");          // perpendicular component of recoil
    outTree->Branch("tkMet",      &tkMet,      "tkMet/F");       // MET (track MET)
    outTree->Branch("tkMetPhi",   &tkMetPhi,   "tkMetPhi/F");    // phi(MET) (track MET)
    outTree->Branch("tkSumEt",    &tkSumEt,    "tkSumEt/F");     // Sum ET (track MET)
    outTree->Branch("tkMt",       &tkMt,       "tkMt/F");        // transverse mass (track MET)
    outTree->Branch("tkU1",       &tkU1,       "tkU1/F");        // parallel component of recoil (track MET)
    outTree->Branch("tkU2",       &tkU2,       "tkU2/F");        // perpendicular component of recoil (track MET)
    outTree->Branch("mvaMet",     &mvaMet,     "mvaMet/F");      // MVA MET
    outTree->Branch("mvaMetPhi",  &mvaMetPhi,  "mvaMetPhi/F");   // phi(MVA MET)
    outTree->Branch("mvaSumEt",   &mvaSumEt,   "mvaSumEt/F");    // Sum ET (mva MET)
    outTree->Branch("mvaMt",      &mvaMt,      "mvaMt/F");      // transverse mass (MVA MET)
    outTree->Branch("mvaU1",      &mvaU1,      "mvaU1/F");       // parallel component of recoil (mva MET)
    outTree->Branch("mvaU2",      &mvaU2,      "mvaU2/F");       // perpendicular component of recoil (mva MET)
    outTree->Branch("puppiMet",    &puppiMet,   "puppiMet/F");      // Puppi MET
    outTree->Branch("puppiMetPhi", &puppiMetPhi,"puppiMetPhi/F");   // phi(Puppi MET)
    outTree->Branch("puppiSumEt",  &puppiSumEt, "puppiSumEt/F");    // Sum ET (Puppi MET)
    outTree->Branch("puppiU1",     &puppiU1,    "puppiU1/F");       // parallel component of recoil (Puppi MET)
    outTree->Branch("puppiU2",     &puppiU2,    "puppiU2/F");       // perpendicular component of recoil (Puppi MET)
    outTree->Branch("q",          &q,          "q/I");           // lepton charge
    outTree->Branch("lep",        "TLorentzVector", &lep);       // lepton 4-vector
    outTree->Branch("lheweight",  "vector<double>", &lheweight);       // lepton 4-vector
    outTree->Branch("lepID",      &lepID,      "lepID/I");       // lepton PDG ID
    ///// muon specific /////
    outTree->Branch("trkIso",     &trkIso,     "trkIso/F");       // track isolation of lepton
    outTree->Branch("emIso",      &emIso,      "emIso/F");        // ECAL isolation of lepton
    outTree->Branch("hadIso",     &hadIso,     "hadIso/F");       // HCAL isolation of lepton
    outTree->Branch("pfChIso",    &pfChIso,    "pfChIso/F");      // PF charged hadron isolation of lepton
    outTree->Branch("pfGamIso",   &pfGamIso,   "pfGamIso/F");     // PF photon isolation of lepton
    outTree->Branch("pfNeuIso",   &pfNeuIso,   "pfNeuIso/F");     // PF neutral hadron isolation of lepton
    outTree->Branch("pfCombIso",  &pfCombIso,  "pfCombIso/F");    // PF combined isolation of lepton
    outTree->Branch("d0",         &d0,         "d0/F");           // transverse impact parameter of lepton
    outTree->Branch("dz",         &dz,         "dz/F");           // longitudinal impact parameter of lepton
    outTree->Branch("muNchi2",    &muNchi2,    "muNchi2/F");      // muon fit normalized chi^2 of lepton
    outTree->Branch("nPixHits",   &nPixHits,   "nPixHits/i");	  // number of pixel hits of muon
    outTree->Branch("nTkLayers",  &nTkLayers,  "nTkLayers/i");	  // number of tracker layers of muon
    outTree->Branch("nMatch",     &nMatch,     "nMatch/i");	  // number of matched segments of muon	 
    outTree->Branch("nValidHits", &nValidHits, "nValidHits/i");   // number of valid muon hits of muon 
    outTree->Branch("typeBits",   &typeBits,   "typeBits/i");     // number of valid muon hits of muon 
    // outTree->Branch("passHLT", &passHLT, "passHLT/b");
    
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
      if(samp->jsonv[ifile].CompareTo("NONE")!=0) { 
        hasJSON = kTRUE;
        rlrm.addJSONFile(samp->jsonv[ifile].Data()); 
      }

      eventTree = (TTree*)infile->Get("Events");
      assert(eventTree);  
      
      
      Bool_t hasJet = eventTree->GetBranchStatus("AK4");
      eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");
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
      Double_t totalWeight=0;
      Double_t totalWeightUp=0;
      Double_t totalWeightDown=0;
      Double_t puWeight=1;
      Double_t puWeightUp=1;
      Double_t puWeightDown=1;

      // if (hasGen) {
      //   for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //     if(ientry%1000000==0) cout << "Pre-processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
      //   // for(UInt_t ientry=0; ientry<(uint)(0.001*eventTree->GetEntries()); ientry++) {
          
      //     // for(int i=0; i < NPDF+NQCD; i++) lheweight.push_back(0);
      //     infoBr->GetEntry(ientry);
      //     genBr->GetEntry(ientry);
      //     // puWeight = doPU ? h_rw->GetBinContent(h_rw->FindBin(info->nPUmean)) : 1.;
      //     // puWeightUp = doPU ? h_rw_up->GetBinContent(h_rw_up->FindBin(info->nPUmean)) : 1.;
      //     // puWeightDown = doPU ? h_rw_down->GetBinContent(h_rw_down->FindBin(info->nPUmean)) : 1.;
      //     totalWeight+=gen->weight*puWeight; // mine has pu and gen separated
      //     totalWeightUp+=gen->weight*puWeightUp;
      //     totalWeightDown+=gen->weight*puWeightDown;
      //   }
      // }
      // else if (not isData){
      //   for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //     if(ientry%1000000==0) cout << "Pre-processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
      //   // for(UInt_t ientry=0; ientry<(uint)(0.001*eventTree->GetEntries()); ientry++) {
      //     // puWeight = doPU ? h_rw->GetBinContent(h_rw->FindBin(info->nPUmean)) : 1.;
      //     // puWeightUp = doPU ? h_rw_up->GetBinContent(h_rw_up->FindBin(info->nPUmean)) : 1.;
      //     // puWeightDown = doPU ? h_rw_down->GetBinContent(h_rw_down->FindBin(info->nPUmean)) : 1.;
      //     totalWeight+= 1.0*puWeight;
      //     totalWeightUp+= 1.0*puWeightUp;
      //     totalWeightDown+= 1.0*puWeightDown;
      //   }

      // }

cout << "n sections " << NSEC << endl;
      double frac = 1.0/NSEC;
      cout << "n sections " << NSEC << "  frac " << frac << endl;
      UInt_t IBEGIN = frac*ITH*eventTree->GetEntries();
      UInt_t IEND = frac*(ITH+1)*eventTree->GetEntries();
      cout << "start, end " << IBEGIN << " " << IEND << endl;
    //
    // loop over events
    //
    Double_t nsel=0, nselvar=0;
    for(UInt_t ientry=IBEGIN; ientry < IEND; ientry++) {
    // for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    // for(UInt_t ientry=0; ientry<(uint)(0.001*eventTree->GetEntries()); ientry++) {
        infoBr->GetEntry(ientry);

        if(ientry%1000000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
        Double_t weight=xsec;
        Double_t weightUp=xsec;
        Double_t weightDown=xsec;
        // passHLT = true;
        // if(xsec>0 && totalWeight>0) weight = xsec/totalWeight;
        // if(xsec>0 && totalWeightUp>0) weightUp = xsec/totalWeightUp;
        // if(xsec>0 && totalWeightDown>0) weightDown = xsec/totalWeightDown;

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
        }else {
          hGenWeights->Fill(0.0,1.0);
        }
        scArr->Clear();
        scBr->GetEntry(ientry);
        
       /* Double_t weight=1;
        if(xsec>0 && totalWeight>0) weight = xsec/totalWeight;
	if(hasGen) {
	  genPartArr->Clear();
	  genBr->GetEntry(ientry);
          genPartBr->GetEntry(ientry);
	  weight*=gen->weight;
	}*/
        
	    // veto w -> xv decays for signal and w -> mv for bacground samples (needed for inclusive WToLNu sample)
        // std::cout << "isWrongFlavor " << isWrongFlavor << " isSignal " << isSignal << " hasGen " << hasGen << " flav " <<  fabs(toolbox::flavor(genPartArr, BOSON_ID)) << std::endl;
        if (isWrongFlavor && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;
        else if (isSignal && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
        // std::cout << "veto1" << std::endl;
	    // check for certified lumi (if applicable)
        baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        if(hasJSON && !rlrm.hasRunLumi(rl)) continue;  
         // std::cout << "veto2" << std::endl;
        // trigger requirement               
        if (!isMuonTrigger(triggerMenu, info->triggerBits,isData,is13TeV)) continue;
      // std::cout << "veto3" << std::endl;
        // good vertex requirement
        if(!(info->hasGoodPV)) continue;
        // std::cout << "veto4" << std::endl;
           
        //
	// SELECTION PROCEDURE:
	//  (1) Look for 1 good muon matched to trigger
	//  (2) Reject event if another muon is present passing looser cuts
	//
      muonArr->Clear();
      muonBr->GetEntry(ientry);
    
        jetArr->Clear();
        if(hasJet)jetBr->GetEntry(ientry);

	Int_t nLooseLep=0;
	const baconhep::TMuon *goodMuon=0;
	Bool_t passSel=kFALSE;

        for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[i]);

          // apply scale and resolution corrections to MC
          Double_t mupt_corr = mu->pt;
          if(doScaleCorr && snamev[isam].CompareTo("data",TString::kIgnoreCase)!=0)
            mupt_corr = gRandom->Gaus(mu->pt*getMuScaleCorr(mu->eta,0),getMuResCorr(mu->eta,0));

          if(fabs(mu->eta) > VETO_ETA) continue; // loose lepton |eta| cut
          if(mupt_corr     < VETO_PT)  continue; // loose lepton pT cut
          // std::cout << "veto5" << std::endl;
          if(passMuonLooseID(mu)) nLooseLep++;   // loose lepton selection
          if(nLooseLep>1) {  // extra lepton veto
            passSel=kFALSE;
            break;
          }
          // std::cout << "veto6" << std::endl;
          
          if(fabs(mu->eta) > ETA_CUT)         continue;  // lepton |eta| cut
          if(mupt_corr     < PT_CUT)          continue;  // lepton pT cut   
          // std::cout << "veto7" << std::endl;
          if(!passMuonID(mu))                 continue;  // lepton selection
          if(!isMuonTriggerObj(triggerMenu, mu->hltMatchBits,isData,is13TeV)) continue;
          // std::cout << "veto8" << std::endl;

	  passSel=kTRUE;
	  goodMuon = mu;
	}



	if(passSel) {
	  /******** We have a W candidate! HURRAY! ********/
	  nsel+=weight;
    nselvar+=weight*weight;
      
    // Loop through the photons to determine the Prefiring scale factor
    // prefireWeight=1; prefireUp=1; prefireDown=1;
    // for(Int_t ip=0; ip<scArr->GetEntriesFast(); ip++) {
      // const baconhep::TPhoton *photon = (baconhep::TPhoton*)((*scArr)[ip]);
      // prefireWeight *= (1.-prefirePhotonCorr.getCorr(photon->eta, photon->pt));
      // prefireUp     *= TMath::Max((1.-(1.2*prefirePhotonCorr.getCorr(photon->eta, photon->pt))),0.0);
      // prefireDown   *= TMath::Max((1.-(0.8*prefirePhotonCorr.getCorr(photon->eta, photon->pt))),0.0);
      // // std::cout << "photon eta " << photon->eta << "  photon pT " << photon->pt << "  prefire weight " << prefireWeight << std::endl;
    // } 

// Loop through Jets
      // set up the met variable, default is PF met
      metDJee      = info->pfMETC;
      metPhiDJee   = info->pfMETCphi;
      sumEtDJee = 0;
      // if(category==1||category==2||category==3) cout << "Selected! " << endl;
      if(hasJet){
        TVector2 vMetEE((info->pfMETC)*cos(info->pfMETCphi),(info->pfMETC)*sin(info->pfMETCphi));
        for(Int_t ip=0; ip<jetArr->GetEntriesFast(); ip++) {
          const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[ip]);
          if(fabs(jet->eta) < 2.65 || fabs(jet->eta) > 3.139) continue;
          if(jet->pt > 50) continue;
          TVector2 vJet((jet->pt)*cos(jet->phi),(jet->pt)*sin(jet->phi));
          vMetEE += vJet;
          // if(category==1||category==2||category==3)cout << "Removing a jet from MET " << jet->pt << " " << jet->eta << " old met " << info->pfMETC << " new met " << vMetEE.Mod() << endl;
        } 
        metDJee      = vMetEE.Mod();
        metPhiDJee   = vMetEE.Phi();
      }
      // sumEtDJee = 0;
        
        if(!isData){
          // Loop through the photons to determine the Prefiring scale factor
          prefirePhoton=1; prefirePhotUp=1; prefirePhotDown=1;
          for(Int_t ip=0; ip<scArr->GetEntriesFast(); ip++) {
            const baconhep::TPhoton *photon = (baconhep::TPhoton*)((*scArr)[ip]);
            if(fabs(photon->eta) < 2 || fabs(photon->eta) > 5) continue;
            prefirePhoton *= 1. - TMath::Max( (double)prefirePhotonCorr.getCorr(photon->eta, photon->pt) , 0.0 );
          } 
          prefirePhotUp = max(prefirePhoton+(1-prefirePhoton)*0.20,1.0);
          prefirePhotDown = max(prefirePhoton-(1-prefirePhoton)*0.20,1.0);
          
        
          prefireJet=1; prefireJetUp=1; prefireJetDown=1;
          if(hasJet){
            for(Int_t ip=0; ip<jetArr->GetEntriesFast(); ip++) {
              const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[ip]);          
              if(fabs(jet->eta) < 2 || fabs(jet->eta) > 5) continue;
              prefireJet*= 1. - TMath::Max((double)prefireJetCorr.getCorr(jet->eta, jet->pt),0.);
            }
          }
          prefireJetUp = max(prefireJet+(1-prefireJet)*0.20,1.0);
          prefireJetDown = max(prefireJet-(1-prefireJet)*0.20,1.0);
          // loop through photons and jets
          // overlap is anything within deltaR < 0.4.
          // take max prefire prob for any overlap cases
          //toolbox::deltaR(jet->eta, jet->phi, photon->eta, photon->phi))<0.4
          // total prefire probability = product of all (1-prob) for photons,jets, & remove the overlap
          prefireWeight=prefireJet*prefirePhoton;
          prefireUp=prefireJetUp*prefirePhotUp;
          prefireDown=prefireJetDown*prefirePhotDown;
          if(hasJet){
            for(Int_t ip=0; ip<scArr->GetEntriesFast(); ip++) {
              const baconhep::TPhoton *photon = (baconhep::TPhoton*)((*scArr)[ip]);
              if(fabs(photon->eta) < 2 || fabs(photon->eta) > 5) continue;
              // now loop through jets:
              double rmP = 1;

              for(Int_t ip=0; ip<jetArr->GetEntriesFast(); ip++) {
                const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[ip]);
                if(fabs(jet->eta) < 2 || fabs(jet->eta) > 5) continue;
                // check if the jet and photon overlap: 
                if(toolbox::deltaR(jet->eta, jet->phi, photon->eta, photon->phi)>0.4) continue;
                // photon & jet overlap, now get min to divide out 
                  rmP = min(TMath::Max( (double)prefirePhotonCorr.getCorr(photon->eta, photon->pt) , 0.0 ), TMath::Max((double)prefireJetCorr.getCorr(jet->eta, jet->pt),0.));
              }
              // divide out the lesser of the two probabilities
              if(rmP<1.0)prefireWeight = prefireWeight / (1 - rmP);
            }
          }
          
          prefireUp = min(prefireWeight+(1-prefireWeight)*0.20,1.0);
          prefireDown = min(prefireWeight-(1-prefireWeight)*0.20,1.0);
        }

    // apply scale and resolution corrections to MC
    Double_t goodMuonpt_corr = goodMuon->pt;
    if(doScaleCorr && snamev[isam].CompareTo("data",TString::kIgnoreCase)!=0)
    goodMuonpt_corr = gRandom->Gaus(goodMuon->pt*getMuScaleCorr(goodMuon->eta,0),getMuResCorr(goodMuon->eta,0));

	  TLorentzVector vLep; 
	  vLep.SetPtEtaPhiM(goodMuonpt_corr, goodMuon->eta, goodMuon->phi, MUON_MASS); 
	  

	  //
	  // Fill tree
	  //
	  runNum    = info->runNum;
	  lumiSec   = info->lumiSec;
	  evtNum    = info->evtNum;
	  
	  vertexArr->Clear();
	  vertexBr->GetEntry(ientry);

	  npv       = vertexArr->GetEntries();
	  npu	    = info->nPUmean;
	  genV      = new TLorentzVector(0,0,0,0);
	  genLep    = new TLorentzVector(0,0,0,0);
	  genVPt    = -999;
	  genVPhi   = -999;
	  genVy     = -999;
	  genVMass  = -999;
	  genLepPt  = -999;
	  genLepPhi = -999;
	  u1        = -999;
	  u2        = -999;
	  tkU1      = -999;
          tkU2      = -999;
	  mvaU1     = -999;
          mvaU2     = -999;
          id_1      = -999;
          id_2      = -999;
          x_1       = -999;
          x_2       = -999;
          xPDF_1    = -999;
          xPDF_2    = -999;
          scalePDF  = -999;
          weightPDF = -999;

    genMuonPt = 0;
    if(hasGen){
      
      // genMuonMatch
      // use the function to get the 
      // std::cout << "---------" << std::endl;
      genMuonPt = toolbox::getGenLep(genPartArr, vLep);
      // std::cout << "gen muon pt " << genMuonPt << std::endl;
      // std::cout << "boson charge " << toolbox::flavor(genPartArr, BOSON_ID) << std::endl;
      
      if(isRecoil&&!isSignal&&!isWrongFlavor){
        // std::cout <<"Filling the Zxx lheweight" << std::endl;
        lheweight[0]=gen->lheweight[0];
        lheweight[1]=gen->lheweight[1];
        lheweight[2]=gen->lheweight[2];
        lheweight[3]=gen->lheweight[3];
        lheweight[4]=gen->lheweight[5];
        lheweight[5]=gen->lheweight[7];
        for(int npdf=0; npdf<NPDF; npdf++) lheweight[npdf]=gen->lheweight[8+npdf];
      }else{
        // std::cout << "filling the lheweight" << std::endl;
        lheweight[0]=gen->lheweight[1];
        lheweight[1]=gen->lheweight[2];
        lheweight[2]=gen->lheweight[3];
        lheweight[3]=gen->lheweight[4];
        lheweight[4]=gen->lheweight[6];
        lheweight[5]=gen->lheweight[8];
        for(int npdf=0; npdf<NPDF; npdf++) lheweight[npdf+NQCD]=gen->lheweight[9+npdf];
        // std::cout << lheweight[0] << "  "  << gen->lheweight[1] << std::endl;
        // std::cout << lheweight[1] << "  "  << gen->lheweight[2] << std::endl;
        // std::cout << lheweight[6] << "  "  << gen->lheweight[9] << std::endl;
      }
    }

	  if(isRecoil && hasGen) {
      Int_t glepq1=-99;
      Int_t glepq2=-99;
	    TLorentzVector *gvec=new TLorentzVector(0,0,0,0);
      TLorentzVector *glep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *glep2=new TLorentzVector(0,0,0,0);
	    toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2,&glepq1,&glepq2,1);
      if((snamev[isam].CompareTo("zxx",TString::kIgnoreCase)==0)){ // DY only
          toolbox::fillGen(genPartArr, 23, gvec, glep1, glep2,&glepq1,&glepq2,1);
      }
   
      TLorentzVector tvec=*glep1+*glep2;
      genV=new TLorentzVector(0,0,0,0);
      genV->SetPtEtaPhiM(tvec.Pt(), tvec.Eta(), tvec.Phi(), tvec.M());
      genVPt   = tvec.Pt();
      genVPhi  = tvec.Phi();
      genVy    = tvec.Rapidity();
      genVMass = tvec.M();
      
      // std::cout << "construct genVPt " << genVPt << std::endl;
      // std::cout << "from " << glep1->Pt() << " and " << glep2->Pt() << std::endl;

          if (gvec && glep1) {
      //genV      = new TLorentzVector(0,0,0,0);
      //genV->SetPtEtaPhiM(gvec->Pt(),gvec->Eta(),gvec->Phi(),gvec->M());
      // std::cout << "boson id "  << BOSON_ID << "  genlepq1 " << glepq1 << " genlepq2 " << glepq2 << std::endl;
      genLep    = new TLorentzVector(0,0,0,0);
      if(toolbox::flavor(genPartArr, BOSON_ID)*glepq1<0){
        // std::cout << "lep1 ! " << BOSON_ID*glepq1 << std::endl;
        genLep->SetPtEtaPhiM(glep1->Pt(),glep1->Eta(),glep1->Phi(),glep1->M());
      }
      if(toolbox::flavor(genPartArr, BOSON_ID)*glepq2<0){
        // std::cout << "lep2 ! " << BOSON_ID*glepq2 << std::endl;
        genLep->SetPtEtaPhiM(glep2->Pt(),glep2->Eta(),glep2->Phi(),glep2->M());
      }
        // genVPt    = gvec->Pt();
        // genVPhi   = gvec->Phi();
        // genVy     = gvec->Rapidity();
        // genVMass  = gvec->M();
        genLepPt  = genLep->Pt();
        // std::cout << "other gen lep pt " << genLepPt << std::endl;
        genLepPhi = genLep->Phi();
        // if(fabs(genLepPt-genMuonPt)/genMuonPt > 0.1) {
          // // std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
          // // std::cout << "lepton code " << toolbox::flavor(genPartArr, BOSON_ID) << std::endl;
          // // std::cout << "genlep1 pT " << genLepPt << "  genmuonpt " << genMuonPt << std::endl;
          // // std::cout << "genlep1 phi " << genLepPhi << "  lepPhi " << vLep.Phi() << std::endl;
          // // std::cout << toolbox::getGenLep(genPartArr, vLep) << std::endl;
          // // std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
        // }
  
        TVector2 vWPt((genVPt)*cos(genVPhi),(genVPt)*sin(genVPhi));
        TVector2 vLepPt(vLep.Px(),vLep.Py());

        TVector2 vMet((info->pfMETC)*cos(info->pfMETCphi), (info->pfMETC)*sin(info->pfMETCphi));
        TVector2 vU = -1.0*(vMet+vLepPt);
        u1 = ((vWPt.Px())*(vU.Px()) + (vWPt.Py())*(vU.Py()))/(genVPt);  // u1 = (pT . u)/|pT|
        u2 = ((vWPt.Px())*(vU.Py()) - (vWPt.Py())*(vU.Px()))/(genVPt);  // u2 = (pT x u)/|pT|
        
        TVector2 vMetDJ((metDJee)*cos(metPhiDJee), (metDJee)*sin(metPhiDJee));
        TVector2 vUDJ = -1.0*(vMetDJ+vLepPt);
        u1DJee = ((vWPt.Px())*(vUDJ.Px()) + (vWPt.Py())*(vUDJ.Py()))/(genVPt);  // u1 = (pT . u)/|pT|
        u2DJee = ((vWPt.Px())*(vUDJ.Py()) - (vWPt.Py())*(vUDJ.Px()))/(genVPt);  // u2 = (pT x u)/|peleProbe	

        TVector2 vTkMet((info->trkMET)*cos(info->trkMETphi), (info->trkMET)*sin(info->trkMETphi));
        TVector2 vTkU = -1.0*(vTkMet+vLepPt);
        tkU1 = ((vWPt.Px())*(vTkU.Px()) + (vWPt.Py())*(vTkU.Py()))/(genVPt);  // u1 = (pT . u)/|pT|
        tkU2 = ((vWPt.Px())*(vTkU.Py()) - (vWPt.Py())*(vTkU.Px()))/(genVPt);  // u2 = (pT x u)/|pT|

        TVector2 vMvaMet((info->mvaMET)*cos(info->mvaMETphi), (info->mvaMET)*sin(info->mvaMETphi));
        TVector2 vMvaU = -1.0*(vMvaMet+vLepPt);
        mvaU1 = ((vWPt.Px())*(vMvaU.Px()) + (vWPt.Py())*(vMvaU.Py()))/(genVPt);  // u1 = (pT . u)/|pT|
        mvaU2 = ((vWPt.Px())*(vMvaU.Py()) - (vWPt.Py())*(vMvaU.Px()))/(genVPt);  // u2 = (pT x u)/|pT|
        
        TVector2 vPuppiMet((info->puppET)*cos(info->puppETphi), (info->puppET)*sin(info->puppETphi));
        TVector2 vPuppiU = -1.0*(vPuppiMet+vLepPt);
        puppiU1 = ((vWPt.Px())*(vPuppiU.Px()) + (vWPt.Py())*(vPuppiU.Py()))/(genVPt);  // u1 = (pT . u)/|pT|
        puppiU2 = ((vWPt.Px())*(vPuppiU.Py()) - (vWPt.Py())*(vPuppiU.Px()))/(genVPt);  // u2 = (pT x u)/|pT|
  
      }
      id_1      = gen->id_1;
      id_2      = gen->id_2;
      x_1       = gen->x_1;
      x_2       = gen->x_2;
      xPDF_1    = gen->xPDF_1;
      xPDF_2    = gen->xPDF_2;
      scalePDF  = gen->scalePDF;
      weightPDF = gen->weight;
      

      

      
      delete gvec;
      delete glep1;
      delete glep2;
      gvec=0; glep1=0; glep2=0;
	  }
	  scale1fb = weight;
    scale1fbUp = weightUp;
    scale1fbDown = weightDown;
    met	   = info->pfMETC;
	  metPhi   = info->pfMETCphi;
	  sumEt    = 0;
	  mt       = sqrt( 2.0 * (vLep.Pt()) * (info->pfMETC) * (1.0-cos(toolbox::deltaPhi(vLep.Phi(),info->pfMETCphi))) );
	  tkMet    = info->trkMET;
	  tkMetPhi = info->trkMETphi;
	  tkSumEt  = 0;
	  tkMt     = sqrt( 2.0 * (vLep.Pt()) * (info->trkMET) * (1.0-cos(toolbox::deltaPhi(vLep.Phi(),info->trkMETphi))) );
	  mvaMet   = info->mvaMET;
	  mvaMetPhi = info->mvaMETphi;
	  mvaSumEt  = 0;
	  mvaMt     = sqrt( 2.0 * (vLep.Pt()) * (info->mvaMET) * (1.0-cos(toolbox::deltaPhi(vLep.Phi(),info->mvaMETphi))) );
          puppiMet = info->puppET;
          puppiMetPhi = info->puppETphi;
	  puppiSumEt  = 0;
	  puppiMt     = sqrt( 2.0 * (vLep.Pt()) * (info->puppET) * (1.0-cos(toolbox::deltaPhi(vLep.Phi(),info->puppETphi))) );
	  q        = goodMuon->q;
	  lep      = &vLep;
	  
	  ///// muon specific /////
	  trkIso     = goodMuon->trkIso;
	  emIso      = goodMuon->ecalIso;
	  hadIso     = goodMuon->hcalIso;
	  pfChIso    = goodMuon->chHadIso;
	  pfGamIso   = goodMuon->gammaIso;
	  pfNeuIso   = goodMuon->neuHadIso;
	  pfCombIso  = goodMuon->chHadIso + TMath::Max(goodMuon->neuHadIso + goodMuon->gammaIso -
						   0.5*(goodMuon->puIso),Double_t(0));
	  d0         = goodMuon->d0;
	  dz         = goodMuon->dz;
	  muNchi2    = goodMuon->muNchi2;
	  nPixHits   = goodMuon->nPixHits;
	  nTkLayers  = goodMuon->nTkLayers;
	  nMatch     = goodMuon->nMatchStn;
	  nValidHits = goodMuon->nValidHits;
	  typeBits   = goodMuon->typeBits;
	  outTree->Fill();
	  delete genV;
	  delete genLep;
	  genV=0, genLep=0, lep=0;
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
