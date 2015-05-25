//================================================================================================
//
// Select Z->mumu candidates
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
#include "TLorentzVector.h"     // 4-vector class

#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

// lumi section selection with JSON files
//#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions
#endif


//=== MAIN MACRO ================================================================================================= 

void selectZmm(const TString conf="zmm.conf", // input file
               const TString outputDir="."    // output directory
) {
  gBenchmark->Start("selectZmm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW  = 40;
  const Double_t MASS_HIGH = 200;
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.4;
  const Double_t MUON_MASS = 0.105658369;

  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 13;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  enum { eMuMu2HLT=1, eMuMu1HLT, eMuMu1HLT1L1, eMuMuNoSel, eMuSta, eMuTrk };  // event category enum
  
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
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  UInt_t  id_1, id_2;
  Double_t x_1, x_2, xPDF_1, xPDF_2;
  Double_t scalePDF, weightPDF;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, u1, u2;
  Float_t tkMet, tkMetPhi, tkSumEt, tkU1, tkU2;
  Int_t   q1, q2;
  Float_t dilep_pt, dilep_eta, dilep_phi, dilep_m;
  Float_t lep1_pt, lep1_eta, lep1_phi, lep1_m;
  Float_t lep2_pt, lep2_eta, lep2_phi, lep2_m;
  ///// muon specific /////
  Float_t trkIso1, emIso1, hadIso1, trkIso2, emIso2, hadIso2;
  Float_t pfChIso1, pfGamIso1, pfNeuIso1, pfCombIso1, pfChIso2, pfGamIso2, pfNeuIso2, pfCombIso2;
  Float_t d01, dz1, d02, dz2;
  Float_t muNchi21,  muNchi22;
  UInt_t nPixHits1, nTkLayers1, nPixHits2, nTkLayers2;
  UInt_t nValidHits1, nMatch1, nValidHits2, nMatch2;
  UInt_t typeBits1, typeBits2;
  Float_t sta1_pt, sta1_eta, sta1_phi, sta1_m;
  Float_t sta2_pt, sta2_eta, sta2_phi, sta2_m;
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr    = new TClonesArray("baconhep::TMuon");
  TClonesArray *pvArr      = new TClonesArray("baconhep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;
  
  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    
    // Assume data sample is first sample in .conf file
    // If sample is empty (i.e. contains no ntuple files), skip to next sample
    if(isam==0 && !hasData) continue;

    // Assume signal sample is given name "zmm"
    // If it's the signal sample, toggle flag to store GEN W kinematics
    Bool_t isSignal = (snamev[isam].CompareTo("zmm",TString::kIgnoreCase)==0);  
    
    CSample* samp = samplev[isam];
  
    //
    // Set up output ntuple
    //
    TString outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");

    outTree->Branch("runNum",      &runNum,     "runNum/i");      // event run number
    outTree->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
    outTree->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
    outTree->Branch("matchGen",    &matchGen,   "matchGen/i");    // event has both leptons matched to MC Z->ll
    outTree->Branch("category",    &category,   "category/i");    // dilepton category
    outTree->Branch("id_1",        &id_1,       "id_1/i");        // PDF info -- parton ID for parton 1
    outTree->Branch("id_2",        &id_2,       "id_2/i");        // PDF info -- parton ID for parton 2
    outTree->Branch("x_1",         &x_1,        "x_1/d");         // PDF info -- x for parton 1
    outTree->Branch("x_2",         &x_2,        "x_2/d");         // PDF info -- x for parton 2
    outTree->Branch("xPDF_1",      &xPDF_1,     "xPDF_1/d");      // PDF info -- x*F for parton 1
    outTree->Branch("xPDF_2",      &xPDF_2,     "xPDF_2/d");      // PDF info -- x*F for parton 2
    outTree->Branch("scalePDF",    &scalePDF,   "scalePDF/d");    // PDF info -- energy scale of parton interaction
    outTree->Branch("weightPDF",   &weightPDF,  "weightPDF/d");   // PDF info -- PDF weight
    outTree->Branch("npv",         &npv,        "npv/i");         // number of primary vertices
    outTree->Branch("npu",         &npu,        "npu/i");         // number of in-time PU events (MC)
    outTree->Branch("genVPt",      &genVPt,     "genVPt/F");      // GEN boson pT (signal MC)
    outTree->Branch("genVPhi",     &genVPhi,    "genVPhi/F");     // GEN boson phi (signal MC)
    outTree->Branch("genVy",       &genVy,      "genVy/F");       // GEN boson rapidity (signal MC)
    outTree->Branch("genVMass",    &genVMass,   "genVMass/F");    // GEN boson mass (signal MC)
    outTree->Branch("scale1fb",    &scale1fb,   "scale1fb/F");    // event weight per 1/fb (MC)
    outTree->Branch("met",         &met,        "met/F");         // MET
    outTree->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)
    outTree->Branch("sumEt",       &sumEt,      "sumEt/F");       // Sum ET
    outTree->Branch("u1",          &u1,         "u1/F");          // parallel component of recoil
    outTree->Branch("u2",          &u2,         "u2/F");          // perpendicular component of recoil
    outTree->Branch("tkMet",       &tkMet,      "tkMet/F");       // MET (track MET)
    outTree->Branch("tkMetPhi",    &tkMetPhi,   "tkMetPhi/F");    // phi(MET) (track MET)
    outTree->Branch("tkSumEt",     &tkSumEt,    "tkSumEt/F");     // Sum ET (track MET)
    outTree->Branch("tkU1",        &tkU1,       "tkU1/F");        // parallel component of recoil (track MET)
    outTree->Branch("tkU2",        &tkU2,       "tkU2/F");        // perpendicular component of recoil (track MET)
    outTree->Branch("q1",          &q1,         "q1/I");          // charge of tag lepton
    outTree->Branch("q2",          &q2,         "q2/I");          // charge of probe lepton
    outTree->Branch("dilep_pt",    &dilep_pt,   "dilep_pt/F");    // pt of di-lepton
    outTree->Branch("dilep_eta",   &dilep_eta,  "dilep_eta/F");   // eta of di-lepton
    outTree->Branch("dilep_phi",   &dilep_phi,  "dilep_phi/F");   // phi of di-lepton
    outTree->Branch("dilep_m",     &dilep_m,    "dilep_m/F");     // m of di-lepton
    outTree->Branch("lep1_pt",     &lep1_pt,    "lep1_pt/F");     // pt of tag lepton
    outTree->Branch("lep1_eta",    &lep1_eta,   "lep1_eta/F");    // eta of tag lepton
    outTree->Branch("lep1_phi",    &lep1_phi,   "lep1_phi/F");    // phi of tag lepton
    outTree->Branch("lep1_m",      &lep1_m,     "lep1_m/F");      // m of tag lepton
    outTree->Branch("lep2_pt",     &lep2_pt,    "lep2_pt/F");     // pt of probe lepton
    outTree->Branch("lep2_eta",    &lep2_eta,   "lep2_eta/F");    // eta of probe lepton
    outTree->Branch("lep2_phi",    &lep2_phi,   "lep2_phi/F");    // phi of probe lepton
    outTree->Branch("lep2_m",      &lep2_m,     "lep2_m/F");      // m of probe lepton      
    ///// muon specific /////
    outTree->Branch("trkIso1",     &trkIso1,     "trkIso1/F");       // track isolation of tag lepton
    outTree->Branch("trkIso2",     &trkIso2,     "trkIso2/F");       // track isolation of probe lepton
    outTree->Branch("emIso1",      &emIso1,      "emIso1/F");        // ECAL isolation of tag lepton
    outTree->Branch("emIso2",      &emIso2,      "emIso2/F");        // ECAL isolation of probe lepton
    outTree->Branch("hadIso1",     &hadIso1,     "hadIso1/F");       // HCAL isolation of tag lepton
    outTree->Branch("hadIso2",     &hadIso2,     "hadIso2/F");       // HCAL isolation of probe lepton
    outTree->Branch("pfChIso1",    &pfChIso1,    "pfChIso1/F");      // PF charged hadron isolation of tag lepton
    outTree->Branch("pfChIso2",    &pfChIso2,    "pfChIso2/F");      // PF charged hadron isolation of probe lepton
    outTree->Branch("pfGamIso1",   &pfGamIso1,   "pfGamIso1/F");     // PF photon isolation of tag lepton
    outTree->Branch("pfGamIso2",   &pfGamIso2,   "pfGamIso2/F");     // PF photon isolation of probe lepton
    outTree->Branch("pfNeuIso1",   &pfNeuIso1,   "pfNeuIso1/F");     // PF neutral hadron isolation of tag lepton
    outTree->Branch("pfNeuIso2",   &pfNeuIso2,   "pfNeuIso2/F");     // PF neutral hadron isolation of probe lepton
    outTree->Branch("pfCombIso1",  &pfCombIso1,  "pfCombIso1/F");    // PF combined isolation of tag lepton
    outTree->Branch("pfCombIso2",  &pfCombIso2,  "pfCombIso2/F");    // PF combined isolation of probe lepton    
    outTree->Branch("d01",         &d01,         "d01/F");           // transverse impact parameter of tag lepton
    outTree->Branch("d02",         &d02,         "d02/F");           // transverse impact parameter of probe lepton	 
    outTree->Branch("dz1",         &dz1,         "dz1/F");           // longitudinal impact parameter of tag lepton
    outTree->Branch("dz2",         &dz2,         "dz2/F");           // longitudinal impact parameter of probe lepton	 
    outTree->Branch("muNchi21",    &muNchi21,    "muNchi21/F");      // muon fit normalized chi^2 of tag lepton
    outTree->Branch("muNchi22",    &muNchi22,    "muNchi22/F");      // muon fit normalized chi^2 of probe lepton
    outTree->Branch("nPixHits1",   &nPixHits1,	 "nPixHits1/i");     // number of pixel hits of tag muon
    outTree->Branch("nPixHits2",   &nPixHits2,	 "nPixHits2/i");     // number of pixel hits of probe muon
    outTree->Branch("nTkLayers1",  &nTkLayers1,  "nTkLayers1/i");    // number of tracker layers of tag muon
    outTree->Branch("nTkLayers2",  &nTkLayers2,  "nTkLayers2/i");    // number of tracker layers of probe muon
    outTree->Branch("nMatch1",     &nMatch1,	 "nMatch1/i");       // number of matched segments of tag muon
    outTree->Branch("nMatch2",     &nMatch2,	 "nMatch2/i");       // number of matched segments of probe muon    
    outTree->Branch("nValidHits1", &nValidHits1, "nValidHits1/i");   // number of valid muon hits of tag muon
    outTree->Branch("nValidHits2", &nValidHits2, "nValidHits2/i");   // number of valid muon hits of probe muon
    outTree->Branch("typeBits1",   &typeBits1,   "typeBits1/i");     // muon type of tag muon
    outTree->Branch("typeBits2",   &typeBits2,   "typeBits2/i");     // muon type of probe muon
    outTree->Branch("sta1_pt",     &sta1_pt,     "sta1_pt/F");       // pt of tag standalone muon
    outTree->Branch("sta1_eta",    &sta1_eta,    "sta1_eta/F");      // eta of tag standalone muon
    outTree->Branch("sta1_phi",    &sta1_phi,    "sta1_phi/F");      // phi of tag standalone muon
    outTree->Branch("sta1_m",      &sta1_m,      "sta1_m/F");        // m of tag standalone muon
    outTree->Branch("sta2_pt",     &sta2_pt,     "sta2_pt/F");       // pt of probe standalone muon
    outTree->Branch("sta2_eta",    &sta2_eta,    "sta2_eta/F");      // eta of probe standalone muon
    outTree->Branch("sta2_phi",    &sta2_phi,    "sta2_phi/F");      // phi of probe standalone muon
    outTree->Branch("sta2_m",      &sta2_m,      "sta2_m/F");        // m of probe standalone muon
    
    //
    // loop through files
    //
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {  

      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... "; cout.flush();
      infile = TFile::Open(samp->fnamev[ifile]); 
      assert(infile);

      const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");
      UInt_t trigger    = triggerMenu.getTriggerBit("HLT_IsoMu20_v*");
      UInt_t trigObjL1  = 6;//triggerMenu.getTriggerObjectBit("HLT_IsoMu20_v*", "hltL1sL1SingleMu16");
      UInt_t trigObjHLT = 7;//triggerMenu.getTriggerObjectBit("HLT_IsoMu20_v*", 
      //"hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09");

      //Bool_t hasJSON = kFALSE;
      //baconhep::RunLumiRangeMap rlrm;
      //if(samp->jsonv[ifile].CompareTo("NONE")!=0) { 
      //  hasJSON = kTRUE;
      //rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
      //}
  
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
      eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
      Bool_t hasGen = eventTree->GetBranchStatus("GenEvtInfo");
      TBranch *genBr=0, *genPartBr=0;
      if(hasGen) {
        eventTree->SetBranchAddress("GenEvtInfo", &gen);
	genBr = eventTree->GetBranch("GenEvtInfo");
	eventTree->SetBranchAddress("GenParticle",&genPartArr);
        genPartBr = eventTree->GetBranch("GenParticle");
      }

      Bool_t hasVer = eventTree->GetBranchStatus("Vertex");
      TBranch *pvBr=0;
      if (hasVer) {
        eventTree->SetBranchAddress("Vertex",       &pvArr);
        pvBr = eventTree->GetBranch("Vertex");
      }
    
      // Compute MC event weight per 1/fb
      Double_t weight = 1;
      const Double_t xsec = samp->xsecv[ifile];
      if(xsec>0) weight = 1000.*xsec/(Double_t)eventTree->GetEntries();     

      //
      // loop over events
      //
      Double_t nsel=0, nselvar=0;

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //for(UInt_t ientry=0; ientry<100; ientry++) {
        infoBr->GetEntry(ientry);
	
	if(genBr) {
	  genBr->GetEntry(ientry);
	  genPartArr->Clear();
          genPartBr->GetEntry(ientry);
	}
     
        // check for certified lumi (if applicable)
        //baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        //if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  

        // trigger requirement               
        if(!(info->triggerBits[trigger])) continue;

        // good vertex requirement
        if(!(info->hasGoodPV)) continue;
        if (hasVer) {
	  pvArr->Clear();
	  pvBr->GetEntry(ientry);
	}

        //
	// SELECTION PROCEDURE:
	//  (1) Find a good muon matched to trigger -> this will be the "tag"
	//  (2) Pair the tag with various probe types which form a tag+probe mass inside 
	//      the Z window and divide candidates into exclusive categories as follows:
	//      (a) if probe is a good muon matched to trigger                  -> MuMu2HLT category
	//      (b) if probe is a good muon not matched to trigger              -> MuMu1HLT category
	//      (c) if probe is a muon failing selection cuts                   -> MuMuNoSel category
	//      (d) if probe is a standalone muon but not global                -> MuSta category
	//      (e) if probe is a tracker muon or non-muon track but not global -> MuTrk category
	//
	muonArr->Clear();
        muonBr->GetEntry(ientry);
        for(Int_t i1=0; i1<muonArr->GetEntriesFast(); i1++) {
          const baconhep::TMuon *tag = (baconhep::TMuon*)((*muonArr)[i1]);
	
	  if(tag->pt        < PT_CUT)        continue;  // lepton pT cut
	  if(fabs(tag->eta) > ETA_CUT)       continue;  // lepton |eta| cut
	  if(!passMuonID(tag))               continue;  // lepton selection
	  if(!(tag->hltMatchBits[trigObjHLT])) continue;  // check trigger matching
	  TLorentzVector vTag;    vTag.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, MUON_MASS);
	  TLorentzVector vTagSta; vTagSta.SetPtEtaPhiM(tag->staPt, tag->staEta, tag->staPhi, MUON_MASS);
	
	  for(Int_t i2=0; i2<muonArr->GetEntriesFast(); i2++) {
	    if(i1==i2) continue;
	    const baconhep::TMuon *probe = (baconhep::TMuon*)((*muonArr)[i2]);
	  
	    if(tag->q == probe->q)         continue;  // opposite charge requirement
	    if(probe->pt        < PT_CUT)  continue;  // lepton pT cut
	    if(fabs(probe->eta) > ETA_CUT) continue;  // lepton |eta| cut
	    
	    TLorentzVector vProbe; vProbe.SetPtEtaPhiM(probe->pt, probe->eta, probe->phi, MUON_MASS);
	    TLorentzVector vProbeSta(0,0,0,0);
	    if(probe->typeBits & baconhep::EMuType::kStandalone)
	      vProbeSta.SetPtEtaPhiM(probe->staPt, probe->staEta, probe->staPhi, MUON_MASS);

	    // mass window
	    TLorentzVector vDilep = vTag + vProbe;
	    if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
	    
	    // determine event category
	    UInt_t icat=0;
	    if(passMuonID(probe)) {
	      if(probe->hltMatchBits[trigObjHLT]) {
	        if(i1>i2) continue;  // make sure we don't double count MuMu2HLT category
		icat=eMuMu2HLT;
		
	      } else if (probe->hltMatchBits[trigObjL1]) {
		icat=eMuMu1HLT1L1;
	      } else {
	        icat=eMuMu1HLT;
	      }
	    } 
	    else if(probe->typeBits & baconhep::EMuType::kGlobal)    { icat=eMuMuNoSel; } 
	    else if(probe->typeBits==baconhep::EMuType::kStandalone) { icat=eMuSta; } 
	    else { icat=eMuTrk; }
	    if(icat==0) continue;
	    
	    /******** We have a Z candidate! HURRAY! ********/
	    
	    nsel+=weight;
            nselvar+=weight*weight;
	    
	    // Perform matching of dileptons to GEN leptons from Z decay
	    Bool_t hasGenMatch = kFALSE;
	    if(isSignal) {
	      TLorentzVector *vec=0, *fvec=0, *lep1=0, *lep2=0;
	      toolbox::fillGen(genPartArr, BOSON_ID, LEPTON_ID, vec, fvec, lep1, lep2);
              Bool_t match1 = ( ((lep1) && toolbox::deltaR(tag->eta, tag->phi, lep1->Eta(), lep1->Phi())<0.5) ||
                                ((lep2) && toolbox::deltaR(tag->eta, tag->phi, lep2->Eta(), lep2->Phi())<0.5) );
	      
              Bool_t match2 = ( ((lep1) && toolbox::deltaR(vProbe.Eta(), vProbe.Phi(), lep1->Eta(), lep1->Phi())<0.5) ||
                                ((lep2) && toolbox::deltaR(vProbe.Eta(), vProbe.Phi(), lep2->Eta(), lep2->Phi())<0.5) );
              if(match1 && match2) {
                hasGenMatch = kTRUE;
                genVPt   = fvec->Pt();
                genVPhi  = fvec->Phi();
                genVy    = fvec->Rapidity();
                genVMass = fvec->M();
              }
              else {
                genVPt   = -999;
                genVPhi  = -999;
                genVy    = -999;
                genVMass = -999;
              }
	    }

	    if (gen) {
              id_1      = gen->id_1;
              id_2      = gen->id_2;
              x_1       = gen->x_1;
              x_2       = gen->x_2;
              xPDF_1    = gen->xPDF_1;
              xPDF_2    = gen->xPDF_2;
              scalePDF  = gen->scalePDF;
              weightPDF = gen->weight;
            }
            else {
              id_1      = -999;
              id_2      = -999;
              x_1       = -999;
              x_2       = -999;
              xPDF_1    = -999;
              xPDF_2    = -999;
              scalePDF  = -999;
              weightPDF = -999;
            }

	    //
	    // Fill tree
	    //
	    runNum   = info->runNum;
	    lumiSec  = info->lumiSec;
	    evtNum   = info->evtNum;
	    matchGen = hasGenMatch ? 1 : 0;
	    category = icat;
	    npv      = hasVer ? pvArr->GetEntriesFast() : 0;
	    npu      = info->nPU;
	    scale1fb = weight;
	    met      = info->pfMET;
	    metPhi   = info->pfMETphi;
	    sumEt    = 0;
	    tkMet    = info->trkMET;
            tkMetPhi = info->trkMETphi;
            tkSumEt  = 0;
	    lep1_pt  = vTag.Pt();
            lep1_eta = vTag.Eta();
            lep1_phi = vTag.Phi();
            lep1_m   = vTag.M();
	    q1       = tag->q;
	    lep2_pt  = vProbe.Pt();
            lep2_eta = vProbe.Eta();
            lep2_phi = vProbe.Phi();
            lep2_m   = vProbe.M();
	    q2       = probe->q;
	    dilep_pt = vDilep.Pt();
            dilep_eta= vDilep.Eta();
            dilep_phi= vDilep.Phi();
            dilep_m  = vDilep.M();
	    
	    TVector2 vZPt((vDilep.Pt())*cos(vDilep.Phi()),(vDilep.Pt())*sin(vDilep.Phi()));        
            TVector2 vMet((info->pfMET)*cos(info->pfMETphi), (info->pfMET)*sin(info->pfMETphi));        
            TVector2 vU = -1.0*(vMet+vZPt);
            u1 = ((vDilep.Px())*(vU.Px()) + (vDilep.Py())*(vU.Py()))/(vDilep.Pt());  // u1 = (pT . u)/|pT|
            u2 = ((vDilep.Px())*(vU.Py()) - (vDilep.Py())*(vU.Px()))/(vDilep.Pt());  // u2 = (pT x u)/|pT|

	    TVector2 vTkMet((info->trkMET)*cos(info->trkMETphi), (info->trkMET)*sin(info->trkMETphi));
            TVector2 vTkU = -1.0*(vTkMet+vZPt);
            tkU1 = ((vDilep.Px())*(vTkU.Px()) + (vDilep.Py())*(vTkU.Py()))/(vDilep.Pt());  // u1 = (pT . u)/|pT|
            tkU2 = ((vDilep.Px())*(vTkU.Py()) - (vDilep.Py())*(vTkU.Px()))/(vDilep.Pt());  // u2 = (pT x u)/|pT|
	  
	    ///// muon specific /////
	    sta1_pt     = vTagSta.Pt();
            sta1_eta    = vTagSta.Eta();
            sta1_phi    = vTagSta.Phi();
            sta1_m      = vTagSta.M();
	    trkIso1     = tag->trkIso;
	    emIso1      = tag->ecalIso;	    
	    hadIso1     = tag->hcalIso;
	    pfChIso1    = tag->chHadIso;
	    pfGamIso1   = tag->gammaIso;
	    pfNeuIso1   = tag->neuHadIso;
	    pfCombIso1  = tag->chHadIso + TMath::Max(tag->neuHadIso + tag->gammaIso - 
						     0.5*(tag->puIso),Double_t(0));
	    d01         = tag->d0;
	    dz1         = tag->dz;
	    muNchi21    = tag->muNchi2;
	    nPixHits1   = tag->nPixHits;
	    nTkLayers1  = tag->nTkLayers;
	    nMatch1     = tag->nMatchStn;
	    nValidHits1 = tag->nValidHits;
	    typeBits1   = tag->typeBits;

	    sta2_pt     = vProbeSta.Pt();
            sta2_eta    = vProbeSta.Eta();
            sta2_phi    = vProbeSta.Phi();
            sta2_m      = vProbeSta.M();	    
	    trkIso2     = probe->trkIso;
	    emIso2      = probe->ecalIso;
	    hadIso2     = probe->hcalIso;
	    pfChIso2    = probe->chHadIso;
	    pfGamIso2   = probe->gammaIso;
	    pfNeuIso2   = probe->neuHadIso;
	    pfCombIso2  = probe->chHadIso + TMath::Max(probe->neuHadIso + probe->gammaIso - 
						       0.5*(probe->puIso),Double_t(0));
	    d02         = probe->d0;
	    dz2         = probe->dz;
	    muNchi22    = probe->muNchi2;
	    nPixHits2   = probe->nPixHits;
	    nTkLayers2  = probe->nTkLayers;
	    nMatch2     = probe->nMatchStn;
	    nValidHits2 = probe->nValidHits;
	    typeBits2   = probe->typeBits;
	    
	    outTree->Fill();
	  }
        }
      }
      delete infile;
      infile=0, eventTree=0;
      
      cout << nsel  << " +/- " << sqrt(nselvar);
      if(isam!=0) cout << " per 1/fb";
      cout << endl;   
    }
    outFile->Write();
    outFile->Close(); 
  }
  delete info;
  delete gen;
  delete muonArr;
  delete pvArr;
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " Z -> mu mu" << endl;
  cout << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectZmm"); 
}
