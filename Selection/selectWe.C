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
#include "Math/LorentzVector.h"     // 4-vector class

#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples
#include "../Utils/MyTools.hh"      // various helper functions

// define structures to read in ntuple
#include "../Ntupler/interface/EWKAnaDefs.hh"
#include "../Ntupler/interface/TEventInfo.hh"
#include "../Ntupler/interface/TGenInfo.hh"
#include "../Ntupler/interface/TElectron.hh"
#include "../Ntupler/interface/TVertex.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// helper functions for lepton ID selection
#include "../Utils/LeptonIDCuts.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void selectWe(const TString conf,        // input file
              const TString outputDir,   // output directory
	      const Bool_t  doScaleCorr  // apply energy scale corrections?
) {
  gBenchmark->Start("selectWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT   = 20;
  const Double_t ETA_CUT  = 2.5;
  const Double_t ELE_MASS = 0.000511;
  
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  const Double_t escaleNbins  = 6;
  const Double_t escaleEta[]  = { 0.4,     0.8,     1.2,     1.4442,  2,        2.5 };
  const Double_t escaleCorr[] = { 1.00284, 1.00479, 1.00734, 1.00851, 1.00001,  0.982898 };
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
  ///// electron specific /////
  Float_t trkIso, emIso, hadIso;
  Float_t pfChIso, pfGamIso, pfNeuIso, pfCombIso;
  Float_t sigieie, hovere, eoverp, fbrem, ecalE;
  Float_t dphi, deta;
  Float_t d0, dz;
  UInt_t  isConv, nexphits, typeBits;
  LorentzVector *sc=0;
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo   *gen   = new mithep::TGenInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  
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
    if(isam==0 && !doScaleCorr) outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.raw.root");
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
    outTree->Branch("lep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep);  // lepton 4-vector
    ///// electron specific /////
    outTree->Branch("trkIso",    &trkIso,    "trkIso/F");     // track isolation of tag lepton
    outTree->Branch("emIso",     &emIso,     "emIso/F");      // ECAL isolation of tag lepton
    outTree->Branch("hadIso",    &hadIso,    "hadIso/F");     // HCAL isolation of tag lepton
    outTree->Branch("pfChIso",   &pfChIso,   "pfChIso/F");    // PF charged hadron isolation of lepton
    outTree->Branch("pfGamIso",  &pfGamIso,  "pfGamIso/F");   // PF photon isolation of lepton
    outTree->Branch("pfNeuIso",  &pfNeuIso,  "pfNeuIso/F");   // PF neutral hadron isolation of lepton
    outTree->Branch("pfCombIso", &pfCombIso, "pfCombIso/F");  // PF combined isolation of electron
    outTree->Branch("sigieie",   &sigieie,   "sigieie/F");    // sigma-ieta-ieta of electron
    outTree->Branch("hovere",    &hovere,    "hovere/F");     // H/E of electron
    outTree->Branch("eoverp",    &eoverp,    "eoverp/F");     // E/p of electron
    outTree->Branch("fbrem",     &fbrem,     "fbrem/F");      // brem fraction of electron
    outTree->Branch("dphi",      &dphi,	     "dphi/F");       // GSF track - ECAL dphi of electron
    outTree->Branch("deta",      &deta,      "deta/F");       // GSF track - ECAL deta of electron
    outTree->Branch("ecalE",     &ecalE,     "ecalE/F");      // ECAL energy of electron
    outTree->Branch("d0",        &d0,        "d0/F");         // transverse impact parameter of electron
    outTree->Branch("dz",        &dz,        "dz/F");         // longitudinal impact parameter of electron
    outTree->Branch("isConv",    &isConv,    "isConv/i");     // conversion filter flag of electron
    outTree->Branch("nexphits",  &nexphits,  "nexphits/i");   // number of missing expected inner hits of electron
    outTree->Branch("typeBits",  &typeBits,  "typeBits/i");   // electron type of electron
    outTree->Branch("sc",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc);   // electron Supercluster 4-vector
    
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
      eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("PV",       &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
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
        ULong64_t trigger = kHLT_Ele22_CaloIdL_CaloIsoVL;
	ULong64_t trigObj = kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj;   
        if(!(info->triggerBits & trigger)) continue;      
      
        // good vertex requirement
        if(!(info->hasGoodPV)) continue;
        pvArr->Clear();
        pvBr->GetEntry(ientry);
      
        //
	// SELECTION PROCEDURE:
	//  (1) Look for 1 good electron matched to trigger
	//  (2) Reject event if another electron is present passing looser cuts
	//
	electronArr->Clear();
        electronBr->GetEntry(ientry);
	Int_t nLooseLep=0;
	const mithep::TElectron *goodEle=0;
	Bool_t passSel=kFALSE;	
        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
	  
	  // check ECAL gap
	  if(fabs(ele->scEta)>=ECAL_GAP_LOW && fabs(ele->scEta)<=ECAL_GAP_HIGH) continue;
	  
	  Double_t escale=1;
	  if(doScaleCorr && isam==0) {
	    for(UInt_t ieta=0; ieta<escaleNbins; ieta++) {
	      if(fabs(ele->scEta)<escaleEta[ieta]) {
	        escale = escaleCorr[ieta];
		break;
	      }
	    }
	  }
	  
	  if(fabs(ele->scEta)   > 2.5) continue;                // loose lepton |eta| cut
          if(escale*(ele->scEt) < 20)  continue;                // loose lepton pT cut
          if(passEleLooseID(ele,info->rhoLowEta)) nLooseLep++;  // loose lepton selection
          if(nLooseLep>1) {  // extra lepton veto
            passSel=kFALSE;
            break;
          }
          
          if(fabs(ele->scEta)   > ETA_CUT)    continue;  // lepton |eta| cut
          if(escale*(ele->scEt) < PT_CUT)     continue;  // lepton pT cut
          if(!passEleID(ele,info->rhoLowEta)) continue;  // lepton selection
          if(!(ele->hltMatchBits & trigObj))  continue;  // check trigger matching
	  
	  passSel=kTRUE;
	  goodEle = ele;  
	}
	
	if(passSel) {
	  
	  /******** We have a W candidate! HURRAY! ********/
	    
	  nsel+=weight;
          nselvar+=weight*weight;
	  
	  Double_t escale=1;
	  if(doScaleCorr && isam==0) {
	    for(UInt_t ieta=0; ieta<escaleNbins; ieta++) {
	      if(fabs(goodEle->scEta)<escaleEta[ieta]) {
	        escale = escaleCorr[ieta];
		break;
	      }
	    }
	  }
	  
	  LorentzVector vLep(escale*(goodEle->pt), goodEle->eta, goodEle->phi, ELE_MASS);  
	  LorentzVector vSC(escale*(goodEle->scEt), goodEle->scEta, goodEle->scPhi, ELE_MASS); 	  
	  
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
	    
	    if(abs(gen->id_1)==EGenType::kElectron) { genLepPt = gen->vpt_1; genLepPhi = gen->vphi_1; }
	    if(abs(gen->id_2)==EGenType::kElectron) { genLepPt = gen->vpt_2; genLepPhi = gen->vphi_2; }
	  }
	  scale1fb = weight;
	  met	   = info->pfMET;
	  metPhi   = info->pfMETphi;
	  sumEt    = info->pfSumET;
	  mt       = sqrt( 2.0 * (vLep.Pt()) * (info->pfMET) * (1.0-cos(toolbox::deltaPhi(vLep.Phi(),info->pfMETphi))) );
	  q        = goodEle->q;	  
	  lep      = &vLep;	  
	  
	  ///// electron specific /////
	  sc	    = &vSC;
	  trkIso    = goodEle->trkIso03;
	  emIso     = goodEle->emIso03;
	  hadIso    = goodEle->hadIso03;
	  pfChIso   = goodEle->pfChIso03;
	  pfGamIso  = goodEle->pfGamIso03;
	  pfNeuIso  = goodEle->pfNeuIso03;	
	  pfCombIso = goodEle->pfChIso03 + TMath::Max(goodEle->pfNeuIso03 + goodEle->pfGamIso03 - (info->rhoLowEta)*getEffArea(goodEle->scEta), 0.);
	  sigieie   = goodEle->sigiEtaiEta;
	  hovere    = goodEle->HoverE;
	  eoverp    = goodEle->EoverP;
	  fbrem     = goodEle->fBrem;
	  dphi      = goodEle->deltaPhiIn;
	  deta      = goodEle->deltaEtaIn;
	  d0        = goodEle->d0;
	  dz        = goodEle->dz;
	  isConv    = goodEle->isConv;
	  nexphits  = goodEle->nExpHitsInner;
	  typeBits  = goodEle->typeBits;
	   
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
  delete electronArr;
  delete pvArr;
  
    
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
