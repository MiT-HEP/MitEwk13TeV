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
#include "Math/LorentzVector.h"     // 4-vector class

#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples
#include "../Utils/MyTools.hh"      // various helper functions

// define structures to read in ntuple
#include "../Ntupler/interface/EWKAnaDefs.hh"
#include "../Ntupler/interface/TEventInfo.hh"
#include "../Ntupler/interface/TGenInfo.hh"
#include "../Ntupler/interface/TMuon.hh"
#include "../Ntupler/interface/TVertex.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// helper functions for lepton ID selection
#include "../Utils/LeptonIDCuts.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void selectWm(const TString conf,      // input file
              const TString outputDir  // output directory
) {
  gBenchmark->Start("selectWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT    = 20;
  const Double_t ETA_CUT   = 2.4;
  const Double_t MUON_MASS = 0.105658369;


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
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t genLepPt, genLepPhi;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  LorentzVector *lep=0;
  ///// muon specific /////
  Float_t trkIso, emIso, hadIso;
  Float_t pfChIso, pfGamIso, pfNeuIso, pfCombIso;
  Float_t d0, dz;
  Float_t muNchi2;
  UInt_t nPixHits, nTkLayers, nValidHits, nMatch, typeBits;
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info = new mithep::TEventInfo();
  mithep::TGenInfo   *gen  = new mithep::TGenInfo();
  TClonesArray *muonArr    = new TClonesArray("mithep::TMuon");
  TClonesArray *pvArr      = new TClonesArray("mithep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;
  
  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    
    // Assume data sample is first sample in .conf file
    // If sample is empty (i.e. contains no ntuple files), skip to next sample
    if(isam==0 && !hasData) continue;
  
    CSample* samp = samplev[isam];
  
    //
    // Set up output ntuple
    //
    TString outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");

    outTree->Branch("runNum",   &runNum,   "runNum/i");     // event run number
    outTree->Branch("lumiSec",  &lumiSec,  "lumiSec/i");    // event lumi section
    outTree->Branch("evtNum",   &evtNum,   "evtNum/i");     // event number
    outTree->Branch("npv",      &npv,      "npv/i");        // number of primary vertices
    outTree->Branch("npu",      &npu,      "npu/i");        // number of in-time PU events (MC)
    outTree->Branch("genVPt",   &genVPt,   "genVPt/F");     // GEN boson pT (signal MC)
    outTree->Branch("genVPhi",  &genVPhi,  "genVPhi/F");    // GEN boson phi (signal MC)
    outTree->Branch("genVy",    &genVy,    "genVy/F");      // GEN boson rapidity (signal MC)
    outTree->Branch("genVMass", &genVMass, "genVMass/F");   // GEN boson mass (signal MC)
    outTree->Branch("genLepPt", &genLepPt, "genLepPt/F");   // GEN lepton pT (signal MC)
    outTree->Branch("genLepPhi",&genLepPhi,"genLepPhi/F");  // GEN lepton phi (signal MC)
    outTree->Branch("scale1fb", &scale1fb, "scale1fb/F");   // event weight per 1/fb (MC)
    outTree->Branch("met",      &met,      "met/F");        // MET
    outTree->Branch("metPhi",   &metPhi,   "metPhi/F");     // phi(MET)
    outTree->Branch("sumEt",    &sumEt,    "sumEt/F");      // Sum ET
    outTree->Branch("mt",       &mt,       "mt/F");         // transverse mass
    outTree->Branch("u1",       &u1,       "u1/F");         // parallel component of recoil
    outTree->Branch("u2",       &u2,       "u2/F");         // perpendicular component of recoil
    outTree->Branch("q",        &q,        "q/I");          // lepton charge
    outTree->Branch("lep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep);   // lepton 4-vector
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
      infile = new TFile(samp->fnamev[ifile]); 
      assert(infile);
      
      Bool_t hasJSON = kFALSE;
      mithep::RunLumiRangeMap rlrm;
      if(samp->jsonv[ifile].CompareTo("NONE")!=0) { 
        hasJSON = kTRUE;
        rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
      }
  
      eventTree = (TTree*)infile->Get("Events");
      assert(eventTree);  
      eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PV",   &pvArr);   TBranch *pvBr   = eventTree->GetBranch("PV");
      Bool_t hasGen = eventTree->GetBranchStatus("Gen");
      TBranch *genBr=0;
      if(hasGen) {
        eventTree->SetBranchAddress("Gen", &gen);
	genBr = eventTree->GetBranch("Gen");
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
        infoBr->GetEntry(ientry);
	
	if(genBr) genBr->GetEntry(ientry);
        
	// check for certified lumi (if applicable)
        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  

        // trigger requirement               
        ULong64_t trigger = kHLT_Mu15_eta2p1;
	ULong64_t trigObj = kHLT_Mu15_eta2p1_MuObj;   
        if(!(info->triggerBits & trigger)) continue;      
      
        // good vertex requirement
        if(!(info->hasGoodPV)) continue;
        pvArr->Clear();
        pvBr->GetEntry(ientry);
           
        //
	// SELECTION PROCEDURE:
	//  (1) Look for 1 good muon matched to trigger
	//  (2) Reject event if another muon is present passing looser cuts
	//
	muonArr->Clear();
        muonBr->GetEntry(ientry);
	Int_t nLooseLep=0;
	const mithep::TMuon *goodMuon=0;
	Bool_t passSel=kFALSE;
        for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);

          if(fabs(mu->eta) > 2.4) continue;      // loose lepton |eta| cut
          if(mu->pt        < 10)  continue;      // loose lepton pT cut
          if(passMuonLooseID(mu)) nLooseLep++;   // loose lepton selection
          if(nLooseLep>1) {  // extra lepton veto
            passSel=kFALSE;
            break;
          }
          
          if(fabs(mu->eta) > ETA_CUT)       continue;  // lepton |eta| cut
	  if(mu->pt < PT_CUT)               continue;  // lepton pT cut   
          if(!passMuonID(mu))               continue;  // lepton selection
          if(!(mu->hltMatchBits & trigObj)) continue;  // check trigger matching
  
	  passSel=kTRUE;
	  goodMuon = mu;
	}
	
	if(passSel) {
		  
	  /******** We have a W candidate! HURRAY! ********/
	    
	  nsel+=weight;
          nselvar+=weight*weight;
	  
	  LorentzVector vLep(goodMuon->pt, goodMuon->eta, goodMuon->phi, MUON_MASS);  
	  	  
	  //
	  // Fill tree
	  //
	  runNum   = info->runNum;
	  lumiSec  = info->lumiSec;
	  evtNum   = info->evtNum;
	  npv	   = pvArr->GetEntriesFast();
	  npu	   = info->nPU;
	  genVPt   = 0;
	  genVPhi  = 0;
	  genVy    = 0;
	  genVMass = 0;
	  genLepPt = 0;
	  genLepPhi= 0;
	  u1       = 0;
	  u2       = 0;
	  if(hasGen) {
	    genVPt   = gen->vpt;
            genVPhi  = gen->vphi;
	    genVy    = gen->vy;
	    genVMass = gen->vmass;
	    TVector2 vWPt((gen->vpt)*cos(gen->vphi),(gen->vpt)*sin(gen->vphi));
	    TVector2 vLepPt(vLep.Px(),vLep.Py());      
            TVector2 vMet((info->pfMET)*cos(info->pfMETphi), (info->pfMET)*sin(info->pfMETphi));        
            TVector2 vU = -1.0*(vMet+vLepPt);
            u1 = ((vWPt.Px())*(vU.Px()) + (vWPt.Py())*(vU.Py()))/(gen->vpt);  // u1 = (pT . u)/|pT|
            u2 = ((vWPt.Px())*(vU.Py()) - (vWPt.Py())*(vU.Px()))/(gen->vpt);  // u2 = (pT x u)/|pT|
	    
	    if(abs(gen->id_1)==EGenType::kMuon) { genLepPt = gen->vpt_1; genLepPhi = gen->vphi_1; }
	    if(abs(gen->id_2)==EGenType::kMuon) { genLepPt = gen->vpt_2; genLepPhi = gen->vphi_2; }
	  }
	  scale1fb = weight;
	  met	   = info->pfMET;
	  metPhi   = info->pfMETphi;
	  sumEt    = info->pfSumET;
	  mt       = sqrt( 2.0 * (vLep.Pt()) * (info->pfMET) * (1.0-cos(toolbox::deltaPhi(vLep.Phi(),info->pfMETphi))) );
	  q        = goodMuon->q;	  
	  lep      = &vLep;	  

	  ///// muon specific /////
	  trkIso     = goodMuon->trkIso03;
	  emIso      = goodMuon->emIso03;
	  hadIso     = goodMuon->hadIso03;
	  pfChIso    = goodMuon->pfChIso04;
	  pfGamIso   = goodMuon->pfGamIso04;
	  pfNeuIso   = goodMuon->pfNeuIso04;	  
	  pfCombIso  = goodMuon->pfChIso04 + TMath::Max(goodMuon->pfNeuIso04 + goodMuon->pfGamIso04 - 0.5*(goodMuon->puIso04),Double_t(0));
	  d0         = goodMuon->d0;
	  dz         = goodMuon->dz;
	  muNchi2    = goodMuon->muNchi2;
	  nPixHits   = goodMuon->nPixHits;
	  nTkLayers  = goodMuon->nTkLayers;
	  nMatch     = goodMuon->nMatch;
	  nValidHits = goodMuon->nValidHits;
	  typeBits   = goodMuon->typeBits;
	  
	  outTree->Fill();
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
  cout << " W -> mu nu" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectWm"); 
}
