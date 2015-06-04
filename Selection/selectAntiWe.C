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
#include "TLorentzVector.h"     // 4-vector class

#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
// lumi section selection with JSON files
//#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions
#endif

//=== MAIN MACRO ================================================================================================= 

void selectAntiWe(const TString conf="we.conf", // input file
                  const TString outputDir=".",   // output directory
	          const Bool_t  doScaleCorr=0   // apply energy scale corrections?
) {
  gBenchmark->Start("selectAntiWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT   = 20;
  const Double_t ETA_CUT  = 2.5;
  const Double_t ELE_MASS = 0.000511;

  const Double_t VETO_PT   = 20;
  const Double_t VETO_ETA  = 2.5;
  
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  const Double_t escaleNbins  = 6;
  const Double_t escaleEta[]  = { 0.4,     0.8,     1.2,     1.4442,  2,        2.5 };
  const Double_t escaleCorr[] = { 1.00284, 1.00479, 1.00734, 1.00851, 1.00001,  0.982898 };

  const Int_t BOSON_ID  = 24;
  const Int_t LEPTON_ID = 11;

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
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Float_t tkMet, tkMetPhi, tkSumEt, tkMt, tkU1, tkU2;
  Int_t   q;
  TLorentzVector *lep=0;
  ///// electron specific /////
  Float_t trkIso, emIso, hadIso;
  Float_t pfChIso, pfGamIso, pfNeuIso, pfCombIso;
  Float_t sigieie, hovere, eoverp, fbrem, ecalE;
  Float_t dphi, deta;
  Float_t d0, dz;
  UInt_t  isConv, nexphits, typeBits;
  TLorentzVector *sc=0;
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info     = new baconhep::TEventInfo();
  baconhep::TGenEventInfo   *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr  = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *electronArr = new TClonesArray("baconhep::TElectron");
  TClonesArray *scArr       = new TClonesArray("baconhep::TPhoton");
  TClonesArray *pvArr       = new TClonesArray("baconhep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;
  
  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    
    // Assume data sample is first sample in .conf file
    // If sample is empty (i.e. contains no ntuple files), skip to next sample
    if(isam==0 && !hasData) continue;

    // Assume signal sample is given name "we" -- flag to store GEN W kinematics
    Bool_t isSignal = (snamev[isam].CompareTo("we",TString::kIgnoreCase)==0);
    // flag to reject W->enu events for wrong flavor backgrounds
    Bool_t isWrongFlavor = (snamev[isam].CompareTo("wx",TString::kIgnoreCase)==0);
  
    CSample* samp = samplev[isam];
  
    //
    // Set up output ntuple
    //
    TString outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    if(isam==0 && !doScaleCorr) outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.raw.root");
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");

    outTree->Branch("runNum",   &runNum,   "runNum/i");      // event run number
    outTree->Branch("lumiSec",  &lumiSec,  "lumiSec/i");     // event lumi section
    outTree->Branch("evtNum",   &evtNum,   "evtNum/i");      // event number
    outTree->Branch("npv",      &npv,      "npv/i");         // number of primary vertices
    outTree->Branch("npu",      &npu,      "npu/i");         // number of in-time PU events (MC)
    outTree->Branch("id_1",     &id_1,     "id_1/i");        // PDF info -- parton ID for parton 1
    outTree->Branch("id_2",     &id_2,     "id_2/i");        // PDF info -- parton ID for parton 2
    outTree->Branch("x_1",      &x_1,      "x_1/d");         // PDF info -- x for parton 1
    outTree->Branch("x_2",      &x_2,      "x_2/d");         // PDF info -- x for parton 2
    outTree->Branch("xPDF_1",   &xPDF_1,   "xPDF_1/d");      // PDF info -- x*F for parton 1
    outTree->Branch("xPDF_2",   &xPDF_2,   "xPDF_2/d");      // PDF info -- x*F for parton 2
    outTree->Branch("scalePDF", &scalePDF, "scalePDF/d");    // PDF info -- energy scale of parton interaction
    outTree->Branch("weightPDF",&weightPDF,"weightPDF/d");   // PDF info -- PDF weight
    outTree->Branch("genV",     "TLorentzVector", &genV);    // GEN boson 4-vector (signal MC)
    outTree->Branch("genLep",   "TLorentzVector", &genLep);  // GEN lepton 4-vector (signal MC)
    outTree->Branch("genVPt",   &genVPt,   "genVPt/F");      // GEN boson pT (signal MC)
    outTree->Branch("genVPhi",  &genVPhi,  "genVPhi/F");     // GEN boson phi (signal MC)
    outTree->Branch("genVy",    &genVy,    "genVy/F");       // GEN boson rapidity (signal MC)
    outTree->Branch("genVMass", &genVMass, "genVMass/F");    // GEN boson mass (signal MC)
    outTree->Branch("genLepPt", &genLepPt, "genLepPt/F");    // GEN lepton pT (signal MC)
    outTree->Branch("genLepPhi",&genLepPhi,"genLepPhi/F");   // GEN lepton phi (signal MC)
    outTree->Branch("scale1fb", &scale1fb, "scale1fb/F");    // event weight per 1/fb (MC)
    outTree->Branch("met",      &met,      "met/F");         // MET
    outTree->Branch("metPhi",   &metPhi,   "metPhi/F");      // phi(MET)
    outTree->Branch("sumEt",    &sumEt,    "sumEt/F");       // Sum ET
    outTree->Branch("mt",       &mt,       "mt/F");          // transverse mass
    outTree->Branch("u1",       &u1,       "u1/F");          // parallel component of recoil
    outTree->Branch("u2",       &u2,       "u2/F");          // perpendicular component of recoil
    outTree->Branch("q",        &q,        "q/I");           // lepton charge
    outTree->Branch("lep",      "TLorentzVector", &lep);        // lepton 4-vector
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
    outTree->Branch("sc",        "TLorentzVector", &sc);      // supercluster 4-vector
    
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
      UInt_t trigger    = triggerMenu.getTriggerBit("HLT_Ele23_WP75_Gsf_v*");
      //need to clean this up
      UInt_t trigObjL1  = 4;//triggerMenu.getTriggerObjectBit("HLT_Ele22_WP75_Gsf_v*", "hltL1sL1SingleEG20");
      UInt_t trigObjHLT = 5;//triggerMenu.getTriggerObjectBit("HLT_Ele23_WP75_Gsf_v*", "hltEle23WP75GsfTrackIsoFilter");

      //Bool_t hasJSON = kFALSE;
      //baconhep::RunLumiRangeMap rlrm;
      //if(samp->jsonv[ifile].CompareTo("NONE")!=0) { 
      //hasJSON = kTRUE;
      //rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
      //}
  
      eventTree = (TTree*)infile->Get("Events");
      assert(eventTree);  
      eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
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
        eventTree->SetBranchAddress("Vertex", &pvArr); pvBr = eventTree->GetBranch("Vertex");
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
	//  (1) Look for 1 good electron matched to trigger
	//  (2) Reject event if another electron is present passing looser cuts
	//
	electronArr->Clear();
        electronBr->GetEntry(ientry);
	Int_t nLooseLep=0;
	const baconhep::TElectron *goodEle=0;
	Bool_t passSel=kFALSE;
        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
          const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[i]);
	  
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
	  
          if(fabs(ele->scEta)   > VETO_ETA) continue; // loose lepton |eta| cut
          if(escale*(ele->scEt) < VETO_PT)  continue; // loose lepton pT cut
//          if(passEleLooseID(ele,info->rhoIso)) nLooseLep++;  // loose lepton selection
          if(nLooseLep>1) {  // extra lepton veto
            passSel=kFALSE;
            break;
          }
          
          if(fabs(ele->scEta)   > ETA_CUT)     continue;  // lepton |eta| cut
          if(escale*(ele->scEt) < PT_CUT)      continue;  // lepton pT cut
          if(!passAntiEleID(ele,info->rhoIso)) continue;  // lepton anti-selection
          if(!(ele->hltMatchBits[trigObjHLT])) continue;  // check trigger matching
	  
	  passSel=kTRUE;
	  goodEle = ele;  
	}

	// veto w -> enu decay for wrong flavor background samples (needed for inclusive WToLNu sample)                                                
        if (isWrongFlavor) {
          TLorentzVector *vec=0, *lep1=0, *lep2=0;
          if (fabs(toolbox::flavor(genPartArr, BOSON_ID, vec, lep1, lep2))==LEPTON_ID) continue;
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
	  	  
	  TLorentzVector vLep; vLep.SetPtEtaPhiM(escale*(goodEle->pt), goodEle->eta, goodEle->phi, ELE_MASS);  
	  TLorentzVector vSC;  vSC.SetPtEtaPhiM(escale*(goodEle->scEt), goodEle->scEta, goodEle->scPhi, ELE_MASS); 	  
	  
	  //
	  // Fill tree
	  //
	  runNum   = info->runNum;
	  lumiSec  = info->lumiSec;
	  evtNum   = info->evtNum;
	  npv	   = hasVer ? pvArr->GetEntriesFast() : 0;
	  npu	   = info->nPU;
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
          id_1      = -999;
          id_2      = -999;
          x_1       = -999;
          x_2       = -999;
          xPDF_1    = -999;
          xPDF_2    = -999;
          scalePDF  = -999;
          weightPDF = -999;
	  if(isSignal) {
	    TLorentzVector *vec=0, *lep1=0, *lep2=0;
	    // veto wrong flavor events for signal sample
            if (fabs(toolbox::flavor(genPartArr, BOSON_ID, vec, lep1, lep2))!=LEPTON_ID) continue;
            if (vec && lep1) {
              genV     = vec;
              genLep   = lep1;
              genVPt   = vec->Pt();
              genVPhi  = vec->Phi();
              genVy    = vec->Rapidity();
              genVMass = vec->M();
              genLepPt = lep1->Pt();
              genLepPhi = lep1->Phi();

	      TVector2 vWPt((genVPt)*cos(genVPhi),(genVPt)*sin(genVPhi));
              TVector2 vLepPt(vLep.Px(),vLep.Py());
              TVector2 vMet((info->pfMET)*cos(info->pfMETphi), (info->pfMET)*sin(info->pfMETphi));
              TVector2 vU = -1.0*(vMet+vLepPt);
              u1 = ((vWPt.Px())*(vU.Px()) + (vWPt.Py())*(vU.Py()))/(genVPt);  // u1 = (pT . u)/|pT|
              u2 = ((vWPt.Px())*(vU.Py()) - (vWPt.Py())*(vU.Px()))/(genVPt);  // u2 = (pT x u)/|pT|
	      
              TVector2 vTkMet((info->trkMET)*cos(info->trkMETphi), (info->trkMET)*sin(info->trkMETphi));
              TVector2 vTkU = -1.0*(vTkMet+vLepPt);
              tkU1 = ((vWPt.Px())*(vTkU.Px()) + (vWPt.Py())*(vTkU.Py()))/(genVPt);  // u1 = (pT . u)/|pT|
              tkU2 = ((vWPt.Px())*(vTkU.Py()) - (vWPt.Py())*(vTkU.Px()))/(genVPt);  // u2 = (pT x u)/|pT|
	    }
	    id_1      = gen->id_1;
            id_2      = gen->id_2;
            x_1       = gen->x_1;
            x_2       = gen->x_2;
            xPDF_1    = gen->xPDF_1;
            xPDF_2    = gen->xPDF_2;
            scalePDF  = gen->scalePDF;
            weightPDF = gen->weight;
	  }
	  scale1fb = weight;
	  met	   = info->pfMET;
	  metPhi   = info->pfMETphi;
	  sumEt    = 0;
	  mt       = sqrt( 2.0 * (vLep.Pt()) * (info->pfMET) * (1.0-cos(toolbox::deltaPhi(vLep.Phi(),info->pfMETphi))) );
	  tkMet    = info->trkMET;
          tkMetPhi = info->trkMETphi;
          tkSumEt  = 0;
          tkMt     = sqrt( 2.0 * (vLep.Pt()) * (info->trkMET) * (1.0-cos(toolbox::deltaPhi(vLep.Phi(),info->trkMETphi))) );
	  q        = goodEle->q;
	  lep      = &vLep;
	  
	  ///// electron specific /////
	  sc	    = &vSC;
	  trkIso    = goodEle->trkIso;
          emIso     = goodEle->ecalIso;
          hadIso    = goodEle->hcalIso;
          pfChIso   = goodEle->chHadIso;
          pfGamIso  = goodEle->gammaIso;
          pfNeuIso  = goodEle->neuHadIso;
          pfCombIso = goodEle->chHadIso + TMath::Max(goodEle->neuHadIso + goodEle->gammaIso -
                                                     (info->rhoIso)*getEffArea(goodEle->scEta), 0.);
          sigieie   = goodEle->sieie;
          hovere    = goodEle->hovere;
	  eoverp    = goodEle->eoverp;
          fbrem     = goodEle->fbrem;
          dphi      = goodEle->dPhiIn;
          deta      = goodEle->dEtaIn;
          ecalE     = goodEle->ecalEnergy;
	  d0        = goodEle->d0;
	  dz        = goodEle->dz;
	  isConv    = goodEle->isConv;
          nexphits  = goodEle->nMissingHits;
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
      
  gBenchmark->Show("selectAntiWe"); 
}
