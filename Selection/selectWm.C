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

// lumi section selection with JSON files
#include "BaconAna/Utils/interface/RunLumiRangeMap.hh"

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions
#endif

//=== MAIN MACRO ================================================================================================= 

void selectWm(const TString conf="wm.conf", // input file
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
  const TString ntupDir = outputDir + TString("/ntuples");
  gSystem->mkdir(ntupDir,kTRUE);
  
  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  UInt_t  id_1, id_2;
  Double_t x_1, x_2, xPDF_1, xPDF_2;
  Double_t scalePDF, weightPDF;
  TLorentzVector *genV=0, *genLep=0;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown, puWeight,puWeightUp,puWeightDown;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Float_t tkMet, tkMetPhi, tkSumEt, tkMt, tkU1, tkU2;
  Float_t mvaMet, mvaMetPhi, mvaSumEt, mvaMt, mvaU1, mvaU2;
  Float_t puppiMet, puppiMetPhi, puppiSumEt, puppiMt, puppiU1, puppiU2;
  Int_t   q;
  TLorentzVector *lep=0;
  Int_t lepID;
  ///// muon specific /////
  Float_t trkIso, emIso, hadIso;
  Float_t pfChIso, pfGamIso, pfNeuIso, pfCombIso;
  Float_t d0, dz;
  Float_t muNchi2;
  UInt_t nPixHits, nTkLayers, nValidHits, nMatch, typeBits;
  
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
    
    // Assume data sample is first sample in .conf file
    // If sample is empty (i.e. contains no ntuple files), skip to next sample
    Bool_t isData=kFALSE;
    if(isam==0 && !hasData) continue;
    else if (isam==0) isData=kTRUE;

    // Assume signal sample is given name "wm" -- flag to store GEN W kinematics
    Bool_t isSignal = (snamev[isam].CompareTo("wm",TString::kIgnoreCase)==0);
    //flag to save the info for recoil corrections
    Bool_t isRecoil = ((snamev[isam].CompareTo("wm",TString::kIgnoreCase)==0)||(snamev[isam].CompareTo("zxx",TString::kIgnoreCase)==0)||(snamev[isam].CompareTo("wx",TString::kIgnoreCase)==0));
    // flag to reject W->mnu events when selecting wrong flavor background events
    Bool_t isWrongFlavor = (snamev[isam].CompareTo("wx",TString::kIgnoreCase)==0);
    
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
    outTree->Branch("genVPt",     &genVPt,     "genVPt/F");      // GEN boson pT (signal MC)
    outTree->Branch("genVPhi",    &genVPhi,    "genVPhi/F");     // GEN boson phi (signal MC)
    outTree->Branch("genVy",      &genVy,      "genVy/F");       // GEN boson rapidity (signal MC)
    outTree->Branch("genVMass",   &genVMass,   "genVMass/F");    // GEN boson mass (signal MC)
    outTree->Branch("genLepPt",   &genLepPt,   "genLepPt/F");    // GEN lepton pT (signal MC)
    outTree->Branch("genLepPhi",  &genLepPhi,  "genLepPhi/F");   // GEN lepton phi (signal MC)
    outTree->Branch("scale1fb",   &scale1fb,   "scale1fb/F");    // event weight per 1/fb (MC)
    outTree->Branch("scale1fbUp",   &scale1fbUp,   "scale1fbUp/F");    // event weight per 1/fb (MC)
    outTree->Branch("scale1fbDown",   &scale1fbDown,   "scale1fbDown/F");    // event weight per 1/fb (MC)
    outTree->Branch("met",        &met,        "met/F");         // MET
    outTree->Branch("metPhi",     &metPhi,     "metPhi/F");      // phi(MET)
    outTree->Branch("sumEt",      &sumEt,      "sumEt/F");       // Sum ET
    outTree->Branch("mt",         &mt,         "mt/F");          // transverse mass
    outTree->Branch("u1",         &u1,         "u1/F");          // parallel component of recoil
    outTree->Branch("u2",         &u2,         "u2/F");          // perpendicular component of recoil
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
      Double_t totalWeight=0;
      Double_t totalWeightUp=0;
      Double_t totalWeightDown=0;

//       if (hasGen) {
// 	TH1D *hall = new TH1D("hall", "", 1,0,1);
// 	eventTree->Draw("0.5>>hall", "GenEvtInfo->weight");
// 	totalWeight=hall->Integral();
// 	delete hall;
// 	hall=0;
//       }

      if (hasGen) {
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      infoBr->GetEntry(ientry);
      genBr->GetEntry(ientry);
      puWeight = h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));
      puWeightUp = h_rw_up->GetBinContent(h_rw_up->FindBin(info->nPUmean));
      puWeightDown = h_rw_down->GetBinContent(h_rw_down->FindBin(info->nPUmean));
      totalWeight+=gen->weight*puWeight; // mine has pu and gen separated
      totalWeightUp+=gen->weight*puWeightUp;
      totalWeightDown+=gen->weight*puWeightDown;
    }
      }
      else if (not isData){
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      puWeight = h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));
      puWeightUp = h_rw_up->GetBinContent(h_rw_up->FindBin(info->nPUmean));
      puWeightDown = h_rw_down->GetBinContent(h_rw_down->FindBin(info->nPUmean));
      totalWeight+= 1.0*puWeight;
      totalWeightUp+= 1.0*puWeightUp;
      totalWeightDown+= 1.0*puWeightDown;
    }

      }

      //
      // loop over events
      //
      Double_t nsel=0, nselvar=0;
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
        infoBr->GetEntry(ientry);

        if(ientry%1000000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
	Double_t weight=1;
	Double_t weightUp=1;
	Double_t weightDown=1;
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
	  weight*=gen->weight*puWeight;
	  weightUp*=gen->weight*puWeightUp;
	  weightDown*=gen->weight*puWeightDown;
	}
       /* Double_t weight=1;
        if(xsec>0 && totalWeight>0) weight = xsec/totalWeight;
	if(hasGen) {
	  genPartArr->Clear();
	  genBr->GetEntry(ientry);
          genPartBr->GetEntry(ientry);
	  weight*=gen->weight;
	}*/

	// veto w -> xv decays for signal and w -> mv for bacground samples (needed for inclusive WToLNu sample)
        if (isWrongFlavor && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;
        else if (isSignal && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
        
	// check for certified lumi (if applicable)
        baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        if(hasJSON && !rlrm.hasRunLumi(rl)) continue;  

        // trigger requirement               
        if (!isMuonTrigger(triggerMenu, info->triggerBits)) continue;
      
        // good vertex requirement
        if(!(info->hasGoodPV)) continue;
           
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

        for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[i]);

          // apply scale and resolution corrections to MC
          Double_t mupt_corr = mu->pt;
          if(doScaleCorr && snamev[isam].CompareTo("data",TString::kIgnoreCase)!=0)
            mupt_corr = gRandom->Gaus(mu->pt*getMuScaleCorr(mu->eta,0),getMuResCorr(mu->eta,0));

          if(fabs(mu->eta) > VETO_ETA) continue; // loose lepton |eta| cut
          if(mupt_corr     < VETO_PT)  continue; // loose lepton pT cut
          if(passMuonLooseID(mu)) nLooseLep++;   // loose lepton selection
          if(nLooseLep>1) {  // extra lepton veto
            passSel=kFALSE;
            break;
          }
          
          if(fabs(mu->eta) > ETA_CUT)         continue;  // lepton |eta| cut
	  if(mupt_corr     < PT_CUT)          continue;  // lepton pT cut   
          if(!passMuonID(mu))                 continue;  // lepton selection
          if(!isMuonTriggerObj(triggerMenu, mu->hltMatchBits, kFALSE)) continue;

	  passSel=kTRUE;
	  goodMuon = mu;
	}

	if(passSel) {
	  /******** We have a W candidate! HURRAY! ********/
	  nsel+=weight;
          nselvar+=weight*weight;
	  
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

	  if(isRecoil && hasGen) {
            Int_t glepq1=-99;
            Int_t glepq2=-99;
	    TLorentzVector *gvec=new TLorentzVector(0,0,0,0);
            TLorentzVector *glep1=new TLorentzVector(0,0,0,0);
            TLorentzVector *glep2=new TLorentzVector(0,0,0,0);
	    toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2,&glepq1,&glepq2,1);
	   
            TLorentzVector tvec=*glep1+*glep2;
            genV=new TLorentzVector(0,0,0,0);
            genV->SetPtEtaPhiM(tvec.Pt(), tvec.Eta(), tvec.Phi(), tvec.M());
            genVPt   = tvec.Pt();
            genVPhi  = tvec.Phi();
            genVy    = tvec.Rapidity();
            genVMass = tvec.M();

            if (gvec && glep1) {
	      //genV      = new TLorentzVector(0,0,0,0);
	      //genV->SetPtEtaPhiM(gvec->Pt(),gvec->Eta(),gvec->Phi(),gvec->M());
	      genLep    = new TLorentzVector(0,0,0,0);
	      if(BOSON_ID*glepq1>0)
                genLep->SetPtEtaPhiM(glep1->Pt(),glep1->Eta(),glep1->Phi(),glep1->M());
              if(BOSON_ID*glepq2>0)
                genLep->SetPtEtaPhiM(glep2->Pt(),glep2->Eta(),glep2->Phi(),glep2->M());
              //genVPt    = gvec->Pt();
              //genVPhi   = gvec->Phi();
              //genVy     = gvec->Rapidity();
              //genVMass  = gvec->M();
              genLepPt  = genLep->Pt();
              genLepPhi = genLep->Phi();
	      
              TVector2 vWPt((genVPt)*cos(genVPhi),(genVPt)*sin(genVPhi));
              TVector2 vLepPt(vLep.Px(),vLep.Py());

              TVector2 vMet((info->pfMETC)*cos(info->pfMETCphi), (info->pfMETC)*sin(info->pfMETCphi));
              TVector2 vU = -1.0*(vMet+vLepPt);
              u1 = ((vWPt.Px())*(vU.Px()) + (vWPt.Py())*(vU.Py()))/(genVPt);  // u1 = (pT . u)/|pT|
              u2 = ((vWPt.Px())*(vU.Py()) - (vWPt.Py())*(vU.Px()))/(genVPt);  // u2 = (pT x u)/|pT|

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
