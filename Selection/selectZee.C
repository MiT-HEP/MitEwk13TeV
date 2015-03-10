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
#include "Math/LorentzVector.h"     // 4-vector class

#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples
#include "../Utils/MyTools.hh"      // various helper functions

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"

// lumi section selection with JSON files
//#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// helper functions for lepton ID selection
#include "../Utils/LeptonIDCuts.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void selectZee(const TString conf,        // input file
               const TString outputDir,   // output directory
	       const Bool_t  doScaleCorr  // apply energy scale corrections?
) {
  gBenchmark->Start("selectZee");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW  = 40;
  const Double_t MASS_HIGH = 200;
  const Double_t PT_CUT    = 20;
  const Double_t ETA_CUT   = 2.5;
  const Double_t ELE_MASS  = 0.000511;
  
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  const Double_t escaleNbins  = 6;
  const Double_t escaleEta[]  = { 0.4,     0.8,     1.2,     1.4442,  2,        2.5 };
  const Double_t escaleCorr[] = { 1.00284, 1.00479, 1.00734, 1.00851, 1.00001,  0.982898 };


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  enum { eEleEle2HLT=1, eEleEle1HLT, eEleEleNoSel, eEleSC };  // event category enum
  
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
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, u1, u2;
  Int_t   q1, q2;
  LorentzVector *dilep=0, *lep1=0, *lep2=0;
  ///// electron specific /////
  Float_t trkIso1, emIso1, hadIso1, trkIso2, emIso2, hadIso2;
  Float_t pfChIso1, pfGamIso1, pfNeuIso1, pfCombIso1, pfChIso2, pfGamIso2, pfNeuIso2, pfCombIso2;
  Float_t sigieie1, hovere1, eoverp1, fbrem1, ecalE1, sigieie2, hovere2, eoverp2, fbrem2, ecalE2;
  Float_t dphi1, deta1, dphi2, deta2;
  Float_t d01, dz1, d02, dz2;
  UInt_t  isConv1, nexphits1, typeBits1, isConv2, nexphits2, typeBits2; 
  LorentzVector *sc1=0, *sc2=0;
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info  = new baconhep::TEventInfo();
  baconhep::TGenEventInfo   *gen   = new baconhep::TGenEventInfo();
  TClonesArray *electronArr = new TClonesArray("baconhep::TElectron");
  TClonesArray *scArr       = new TClonesArray("baconhep::TPhoton");
  //TClonesArray *pvArr       = new TClonesArray("baconhep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;
  
  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    
    // Assume data sample is first sample in .conf file
    // If sample is empty (i.e. contains no ntuple files), skip to next sample
    if(isam==0 && !hasData) continue;
    
    // Assume signal sample is given name "zee"
    // If it's the signal sample, toggle flag to store GEN W kinematics
    Bool_t isSignal = (snamev[isam].CompareTo("zee",TString::kIgnoreCase)==0);  
    
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
    outTree->Branch("matchGen", &matchGen, "matchGen/i");   // event has both leptons matched to MC Z->ll
    outTree->Branch("category", &category, "category/i");   // dilepton category
    outTree->Branch("npv",      &npv,      "npv/i");        // number of primary vertices
    outTree->Branch("npu",      &npu,      "npu/i");        // number of in-time PU events (MC)
    outTree->Branch("genVPt",   &genVPt,   "genVPt/F");     // GEN boson pT (signal MC)
    outTree->Branch("genVPhi",  &genVPhi,  "genVPhi/F");    // GEN boson phi (signal MC)
    outTree->Branch("genVy",    &genVy,    "genVy/F");      // GEN boson rapidity (signal MC)
    outTree->Branch("genVMass", &genVMass, "genVMass/F");   // GEN boson mass (signal MC)
    outTree->Branch("scale1fb", &scale1fb, "scale1fb/F");   // event weight per 1/fb (MC)
    outTree->Branch("met",      &met,      "met/F");        // MET
    outTree->Branch("metPhi",   &metPhi,   "metPhi/F");     // phi(MET)
    outTree->Branch("sumEt",    &sumEt,    "sumEt/F");      // Sum ET
    outTree->Branch("u1",       &u1,       "u1/F");         // parallel component of recoil
    outTree->Branch("u2",       &u2,       "u2/F");         // perpendicular component of recoil
    outTree->Branch("q1",       &q1,       "q1/I");         // charge of tag lepton
    outTree->Branch("q2",       &q2,       "q2/I");         // charge of probe lepton
    outTree->Branch("dilep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &dilep);  // dilepton 4-vector
    outTree->Branch("lep1",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1);   // tag lepton 4-vector
    outTree->Branch("lep2",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep2);   // probe lepton 4-vector
    ///// electron specific /////
    outTree->Branch("trkIso1",    &trkIso1,    "trkIso1/F");     // track isolation of tag lepton
    outTree->Branch("trkIso2",    &trkIso2,    "trkIso2/F");     // track isolation of probe lepton
    outTree->Branch("emIso1",     &emIso1,     "emIso1/F");      // ECAL isolation of tag lepton
    outTree->Branch("emIso2",     &emIso2,     "emIso2/F");      // ECAL isolation of probe lepton
    outTree->Branch("hadIso1",    &hadIso1,    "hadIso1/F");     // HCAL isolation of tag lepton
    outTree->Branch("hadIso2",    &hadIso2,    "hadIso2/F");     // HCAL isolation of probe lepton
    outTree->Branch("pfChIso1",   &pfChIso1,   "pfChIso1/F");    // PF charged hadron isolation of tag lepton
    outTree->Branch("pfChIso2",   &pfChIso2,   "pfChIso2/F");    // PF charged hadron isolation of probe lepton
    outTree->Branch("pfGamIso1",  &pfGamIso1,  "pfGamIso1/F");   // PF photon isolation of tag lepton
    outTree->Branch("pfGamIso2",  &pfGamIso2,  "pfGamIso2/F");   // PF photon isolation of probe lepton
    outTree->Branch("pfNeuIso1",  &pfNeuIso1,  "pfNeuIso1/F");   // PF neutral hadron isolation of tag lepton
    outTree->Branch("pfNeuIso2",  &pfNeuIso2,  "pfNeuIso2/F");   // PF neutral hadron isolation of probe lepton
    outTree->Branch("pfCombIso1", &pfCombIso1, "pfCombIso1/F");  // PF combine isolation of tag lepton
    outTree->Branch("pfCombIso2", &pfCombIso2, "pfCombIso2/F");  // PF combined isolation of probe lepton    
    outTree->Branch("sigieie1",   &sigieie1,   "sigieie1/F");    // sigma-ieta-ieta of tag
    outTree->Branch("sigieie2",   &sigieie2,   "sigieie2/F");    // sigma-ieta-ieta of probe
    outTree->Branch("hovere1",    &hovere1,    "hovere1/F");     // H/E of tag
    outTree->Branch("hovere2",    &hovere2,    "hovere2/F");     // H/E of probe
    outTree->Branch("eoverp1",    &eoverp1,    "eoverp1/F");     // E/p of tag
    outTree->Branch("eoverp2",    &eoverp2,    "eoverp2/F");     // E/p of probe	 
    outTree->Branch("fbrem1",     &fbrem1,     "fbrem1/F");      // brem fraction of tag
    outTree->Branch("fbrem2",     &fbrem2,     "fbrem2/F");      // brem fraction of probe
    outTree->Branch("dphi1",      &dphi1,      "dphi1/F");       // GSF track - ECAL dphi of tag
    outTree->Branch("dphi2",      &dphi2,      "dphi2/F");       // GSF track - ECAL dphi of probe 	
    outTree->Branch("deta1",      &deta1,      "deta1/F");       // GSF track - ECAL deta of tag
    outTree->Branch("deta2",      &deta2,      "deta2/F");       // GSF track - ECAL deta of probe
    outTree->Branch("ecalE1",     &ecalE1,     "ecalE1/F");      // ECAL energy of tag
    outTree->Branch("ecalE2",     &ecalE2,     "ecalE2/F");      // ECAL energy of probe
    outTree->Branch("d01",        &d01,        "d01/F");	 // transverse impact parameter of tag
    outTree->Branch("d02",        &d02,        "d02/F");	 // transverse impact parameter of probe	  
    outTree->Branch("dz1",        &dz1,        "dz1/F");	 // longitudinal impact parameter of tag
    outTree->Branch("dz2",        &dz2,        "dz2/F");	 // longitudinal impact parameter of probe
    outTree->Branch("isConv1",    &isConv1,    "isConv1/i");     // conversion filter flag of tag lepton
    outTree->Branch("isConv2",    &isConv2,    "isConv2/i");     // conversion filter flag of probe lepton
    outTree->Branch("nexphits1",  &nexphits1,  "nexphits1/i");   // number of missing expected inner hits of tag lepton
    outTree->Branch("nexphits2",  &nexphits2,  "nexphits2/i");   // number of missing expected inner hits of probe lepton
    outTree->Branch("typeBits1",  &typeBits1,  "typeBits1/i");   // electron type of tag lepton
    outTree->Branch("typeBits2",  &typeBits2,  "typeBits2/i");   // electron type of probe lepton
    outTree->Branch("sc1",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc1);   // tag Supercluster 4-vector
    outTree->Branch("sc2",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc2);   // probe Supercluster 4-vector 
    
    //
    // loop through files
    //
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {  

      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... "; cout.flush();
      infile = new TFile(samp->fnamev[ifile]); 
      assert(infile);

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
      eventTree->SetBranchAddress("Photon",   &scArr);       TBranch *scBr       = eventTree->GetBranch("Photon");
      //eventTree->SetBranchAddress("PV",       &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
      Bool_t hasGen = eventTree->GetBranchStatus("Gen");
      TBranch *genBr=0;
      if(hasGen) {
        eventTree->SetBranchAddress("GenEvtInfo", &gen);
	genBr = eventTree->GetBranch("GenEvtInfo");
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
        //baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        //if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  

        // trigger requirement               
        //ULong64_t trigger = kHLT_Ele22_CaloIdL_CaloIsoVL;
	//ULong64_t trigObj = kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj;  
        //if(!(info->triggerBits & trigger)) continue;      
      
        // good vertex requirement
        //if(!(info->hasGoodPV)) continue;
        //pvArr->Clear();
        //pvBr->GetEntry(ientry);
      
        //
	// SELECTION PROCEDURE:
	//  (1) Find a good electron matched to trigger -> this will be the "tag"
	//  (2) Pair the tag with Supercluster probes which form a tag+probe mass inside 
	//      the Z window and divide candidates into exclusive categories as follows:
	//      (a) if probe SC is part of a good electron matched to trigger     -> EleEle2HLT category
	//      (b) if probe SC is part of a good electron not matched to trigger -> EleEle1HLT category
	//      (c) if probe SC is part of an electron failing selection cuts     -> EleEleNoSel category
	//      (d) if probe SC is not part of an ECAL driven electron            -> EleSC category
	//	
	electronArr->Clear();
        electronBr->GetEntry(ientry);
	scArr->Clear();
	scBr->GetEntry(ientry);
        for(Int_t i1=0; i1<electronArr->GetEntriesFast(); i1++) {
          const baconhep::TElectron *tag = (baconhep::TElectron*)((*electronArr)[i1]);
	  
	  // check ECAL gap
	  if(fabs(tag->scEta)>=ECAL_GAP_LOW && fabs(tag->scEta)<=ECAL_GAP_HIGH) continue;
	  
	  Double_t escale1=1;
	  if(doScaleCorr && isam==0) {
	    for(UInt_t ieta=0; ieta<escaleNbins; ieta++) {
	      if(fabs(tag->scEta)<escaleEta[ieta]) {
	        escale1 = escaleCorr[ieta];
		break;
	      }
	    }
	  }
	  
	  if(escale1*(tag->scEt) < PT_CUT)    continue;  // lepton pT cut
	  if(fabs(tag->scEta)    > ETA_CUT)   continue;  // lepton |eta| cut
	  if(!passEleID(tag,info->rhoIso)) continue;  // lepton selection
	  //if(!(tag->hltMatchBits & trigObj))  continue;  // check trigger matching
	  
	  LorentzVector vTag(escale1*(tag->pt), tag->eta, tag->phi, ELE_MASS);
	  LorentzVector vTagSC(escale1*(tag->scEt), tag->scEta, tag->scPhi, ELE_MASS);
	
	  for(Int_t j=0; j<scArr->GetEntriesFast(); j++) {
	    const baconhep::TPhoton *scProbe = (baconhep::TPhoton*)((*scArr)[j]);
	    if(scProbe->scID == tag->scID) continue;
	    
	    // check ECAL gap
	    if(fabs(scProbe->scEta)>=ECAL_GAP_LOW && fabs(scProbe->scEta)<=ECAL_GAP_HIGH) continue;
	    
	   Double_t escale2=1;
	    if(doScaleCorr && isam==0) {
	      for(UInt_t ieta=0; ieta<escaleNbins; ieta++) {
	        if(fabs(scProbe->scEta)<escaleEta[ieta]) {
	          escale2 = escaleCorr[ieta];
		  break;
	        }
	      }
	    }
	    
	    if(escale2*(scProbe->pt) < PT_CUT)  continue;  // Supercluster ET cut ("pt" = corrected by PV position)
	    if(fabs(scProbe->scEta)  > ETA_CUT) continue;  // Supercluster |eta| cuts
	    
	    const baconhep::TElectron *eleProbe=0;
	    Int_t iprobe=-1;
	    for(Int_t i2=0; i2<electronArr->GetEntriesFast(); i2++) {
	      if(i1==i2) continue;
	      const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[i2]);
	      if(!(ele->typeBits & baconhep::EEleType::kEcalDriven)) continue;
	      if(scProbe->scID==ele->scID) { 
	        eleProbe = ele; 
		iprobe   = i2;
		break; 
	      }
            }
	    
	    LorentzVector vProbe((eleProbe) ? escale2*(eleProbe->pt) : escale2*(scProbe->pt), 
	                         (eleProbe) ? eleProbe->eta : scProbe->eta, 
				 (eleProbe) ? eleProbe->phi : scProbe->phi, 
				 ELE_MASS);
	    LorentzVector vProbeSC((eleProbe) ? escale2*(eleProbe->scEt) : escale2*(scProbe->pt), 
	                           scProbe->scEta, scProbe->scPhi, ELE_MASS);
	    
	    // mass window
	    LorentzVector vDilep = vTag + vProbe;
	    if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
	    
	    // determine event category
	    UInt_t icat=0;
	    if(eleProbe) {
	      if(passEleID(eleProbe,info->rhoIso)) {
	        //if(eleProbe->hltMatchBits & trigObj) {
		//if(i1>iprobe) continue;  // make sure we don't double count EleEle2HLT category
		//icat=eEleEle2HLT;
		
		//} else {
		//icat=eEleEle1HLT;
		//}	      
	      } else { 
	        icat=eEleEleNoSel;
	      } 
	    } else { 
	      icat=eEleSC;
	    }
	    if(icat==0) continue;
	    
	    
	    /******** We have a Z candidate! HURRAY! ********/
	    
	    nsel+=weight;
            nselvar+=weight*weight;
	    
	    // Perform matching of dileptons to GEN leptons from Z decay
	    Bool_t hasGenMatch = kFALSE;
	    if(isSignal) {
	      Bool_t match1 = 0;//( (abs(gen->id_1)==11) && ((toolbox::deltaR(tag->eta, tag->phi, gen->eta_1, gen->phi_1) < 0.5)) )
	      //|| ( (abs(gen->id_2)==11) && ((toolbox::deltaR(tag->eta, tag->phi, gen->eta_2, gen->phi_2) < 0.5)) );

	      Bool_t match2 = 0;//( (abs(gen->id_1)==11) && ((toolbox::deltaR(vProbe.Eta(), vProbe.Phi(), gen->eta_1, gen->phi_1) < 0.5)) )
	      //|| ( (abs(gen->id_2)==11) && ((toolbox::deltaR(vProbe.Eta(), vProbe.Phi(), gen->eta_2, gen->phi_2) < 0.5)) );

	      if(match1 && match2) hasGenMatch = kTRUE;
	    };
	    
	    //
	    // Fill tree
	    //
	    runNum   = info->runNum;
	    lumiSec  = info->lumiSec;
	    evtNum   = info->evtNum;
	    matchGen = hasGenMatch ? 1 : 0;
	    category = icat;
	    npv      = 0; //pvArr->GetEntriesFast();
	    npu      = info->nPU;
	    //genVPt   = (hasGen) ? gen->vpt   : 0;
	    genVPt   = 0;
	    //genVPhi  = (hasGen) ? gen->vphi  : 0;
	    genVPhi  = 0;
	    //genVy    = (hasGen) ? gen->vy    : 0;
	    genVy    = 0;
	    //genVMass = (hasGen) ? gen->vmass : 0;
	    genVMass = 0;
	    scale1fb = weight;
	    met      = info->pfMET;
	    metPhi   = info->pfMETphi;
	    //sumEt    = info->pfSumET;
	    sumEt    = 0;
	    
	    lep1 = &vTag;
	    q1   = tag->q;
	    
	    lep2 = &vProbe;
	    q2   = (eleProbe) ? eleProbe->q : -(tag->q);
	    
	    dilep = &vDilep;
	    
	    TVector2 vZPt((vDilep.Pt())*cos(vDilep.Phi()),(vDilep.Pt())*sin(vDilep.Phi()));        
            TVector2 vMet((info->pfMET)*cos(info->pfMETphi), (info->pfMET)*sin(info->pfMETphi));        
            TVector2 vU = -1.0*(vMet+vZPt);
            u1 = ((vDilep.Px())*(vU.Px()) + (vDilep.Py())*(vU.Py()))/(vDilep.Pt());  // u1 = (pT . u)/|pT|
            u2 = ((vDilep.Px())*(vU.Py()) - (vDilep.Py())*(vU.Px()))/(vDilep.Pt());  // u2 = (pT x u)/|pT|
	  
	    ///// electron specific /////
	    sc1        = &vTagSC;
	    trkIso1    = tag->trkIso;
	    emIso1     = tag->ecalIso;
	    hadIso1    = tag->hcalIso;
	    pfChIso1   = tag->chHadIso;
	    pfGamIso1  = tag->gammaIso;	    
	    pfNeuIso1  = tag->neuHadIso;
	    pfCombIso1 = tag->chHadIso + TMath::Max(tag->neuHadIso + tag->gammaIso - (info->rhoIso)*getEffArea(tag->scEta), 0.);
	    sigieie1   = tag->sieie;
	    hovere1    = tag->hovere;
	    eoverp1    = tag->eoverp;
	    fbrem1     = tag->fbrem;
	    dphi1      = tag->dPhiIn;
	    deta1      = tag->dEtaIn;
	    ecalE1     = tag->ecalEnergy;
	    d01        = tag->d0;
	    dz1        = tag->dz;
	    isConv1    = tag->isConv;
	    nexphits1  = tag->nMissingHits;
	    typeBits1  = tag->typeBits;
	    
	    sc2        = &vProbeSC;
	    trkIso2    = (eleProbe) ? eleProbe->trkIso        : -1;
	    emIso2     = (eleProbe) ? eleProbe->ecalIso       : -1;
	    hadIso2    = (eleProbe) ? eleProbe->hcalIso       : -1;
	    pfChIso2   = (eleProbe) ? eleProbe->chHadIso      : -1;
	    pfGamIso2  = (eleProbe) ? eleProbe->gammaIso      : -1;
	    pfNeuIso2  = (eleProbe) ? eleProbe->neuHadIso     : -1;
	    
	    pfCombIso2 = (eleProbe) ? 
	                 eleProbe->chHadIso + TMath::Max(eleProbe->neuHadIso + eleProbe->gammaIso - (info->rhoIso)*getEffArea(eleProbe->scEta), 0.) : 
			 -1;
	    
	    sigieie2   = (eleProbe) ? eleProbe->sieie         : scProbe->sieie;
	    hovere2    = (eleProbe) ? eleProbe->hovere        : scProbe->hovere;
	    eoverp2    = (eleProbe) ? eleProbe->eoverp        : -1;
	    fbrem2     = (eleProbe) ? eleProbe->fbrem         : -1;
	    dphi2      = (eleProbe) ? eleProbe->dPhiIn        : -999;
	    deta2      = (eleProbe) ? eleProbe->dEtaIn        : -999;
	    ecalE2     = (eleProbe) ? eleProbe->ecalEnergy    : -999;
	    d02        = (eleProbe) ? eleProbe->d0            : -999;
	    dz2        = (eleProbe) ? eleProbe->dz            : -999;
	    isConv2    = (eleProbe) ? eleProbe->isConv        : 0;
	    nexphits2  = (eleProbe) ? eleProbe->nMissingHits  : 0;
	    typeBits2  = (eleProbe) ? eleProbe->typeBits      : 0; 
	    
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
  delete electronArr;
  delete scArr;
  //delete pvArr;
  
    
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
  if(doScaleCorr)
    cout << "  *** Scale corrections applied ***" << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectZee"); 
}
