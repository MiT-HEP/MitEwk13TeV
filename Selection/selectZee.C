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
#include "TRandom.h"

#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples
#include "../Utils/LeptonCorr.hh"   // electron scale and resolution corrections

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
#include "BaconAna/Utils/interface/RunLumiRangeMap.hh"

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions
#endif

//=== MAIN MACRO ================================================================================================= 

void selectZee(const TString conf="zee.conf", // input file
               const TString outputDir=".",   // output directory
	       const Bool_t  doScaleCorr=0    // apply energy scale corrections?
) {
  gBenchmark->Start("selectZee");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW  = 40;
  const Double_t MASS_HIGH = 200;
  const Double_t PT_CUT    = 22;
  const Double_t ETA_CUT   = 2.5;
  const Double_t ELE_MASS  = 0.000511;
  
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  const Double_t escaleNbins  = 2;
  const Double_t escaleEta[]  = { 1.4442, 2.5   };
  const Double_t escaleCorr[] = { 0.992,  1.009 };

  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 11;

  // load trigger menu
  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/pileup_rw_Golden.root", "read");

  // for systematics we need 3
  TH1D *h_rw = (TH1D*) f_rw->Get("h_rw_golden");
  TH1D *h_rw_up = (TH1D*) f_rw->Get("h_rw_up_golden");
  TH1D *h_rw_down = (TH1D*) f_rw->Get("h_rw_down_golden");

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
  TLorentzVector *genV=0;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t scale1fb,scale1fbUp,scale1fbDown;
  Float_t met, metPhi, sumEt, u1, u2;
  Float_t tkMet, tkMetPhi, tkSumEt, tkU1, tkU2;
  Float_t mvaMet, mvaMetPhi, mvaSumEt, mvaU1, mvaU2;
  Float_t puppiMet, puppiMetPhi, puppiSumEt, puppiU1, puppiU2;
  Int_t   q1, q2;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0;
  ///// electron specific /////
  Float_t trkIso1, emIso1, hadIso1, trkIso2, emIso2, hadIso2;
  Float_t pfChIso1, pfGamIso1, pfNeuIso1, pfCombIso1, pfChIso2, pfGamIso2, pfNeuIso2, pfCombIso2;
  Float_t sigieie1, hovere1, eoverp1, fbrem1, ecalE1, sigieie2, hovere2, eoverp2, fbrem2, ecalE2;
  Float_t dphi1, deta1, dphi2, deta2;
  Float_t d01, dz1, d02, dz2;
  UInt_t  isConv1, nexphits1, typeBits1, isConv2, nexphits2, typeBits2; 
  TLorentzVector *sc1=0, *sc2=0;
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr     = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *electronArr    = new TClonesArray("baconhep::TElectron");
  TClonesArray *scArr          = new TClonesArray("baconhep::TPhoton");
  TClonesArray *vertexArr      = new TClonesArray("baconhep::TVertex");

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
    
    // Assume signal sample is given name "zee" - flag to store GEN Z kinematics
    Bool_t isSignal = (snamev[isam].CompareTo("zee",TString::kIgnoreCase)==0);  
    // flag to reject Z->ee events when selecting at wrong-flavor background events
    Bool_t isWrongFlavor = (snamev[isam].CompareTo("zxx",TString::kIgnoreCase)==0);  
    
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
    outTree->Branch("id_1",       &id_1,       "id_1/i");        // PDF info -- parton ID for parton 1
    outTree->Branch("id_2",       &id_2,       "id_2/i");        // PDF info -- parton ID for parton 2
    outTree->Branch("x_1",        &x_1,        "x_1/d");         // PDF info -- x for parton 1
    outTree->Branch("x_2",        &x_2,        "x_2/d");         // PDF info -- x for parton 2
    outTree->Branch("xPDF_1",     &xPDF_1,     "xPDF_1/d");      // PDF info -- x*F for parton 1
    outTree->Branch("xPDF_2",     &xPDF_2,     "xPDF_2/d");      // PDF info -- x*F for parton 2
    outTree->Branch("scalePDF",   &scalePDF,   "scalePDF/d");    // PDF info -- energy scale of parton interaction
    outTree->Branch("weightPDF",  &weightPDF,  "weightPDF/d");   // PDF info -- PDF weight
    outTree->Branch("npv",        &npv,        "npv/i");         // number of primary vertices
    outTree->Branch("npu",        &npu,        "npu/i");         // number of in-time PU events (MC)
    outTree->Branch("genV",      "TLorentzVector",  &genV);      // GEN boson 4-vector
    outTree->Branch("genVPt",     &genVPt,     "genVPt/F");      // GEN boson pT (signal MC)
    outTree->Branch("genVPhi",    &genVPhi,    "genVPhi/F");     // GEN boson phi (signal MC)
    outTree->Branch("genVy",      &genVy,      "genVy/F");       // GEN boson rapidity (signal MC)
    outTree->Branch("genVMass",   &genVMass,   "genVMass/F");    // GEN boson mass (signal MC)
    outTree->Branch("scale1fb",   &scale1fb,   "scale1fb/F");    // event weight per 1/fb (MC)
    outTree->Branch("scale1fbUp",    &scale1fbUp,   "scale1fbUp/F");    // event weight per 1/fb (MC)
    outTree->Branch("scale1fbDown",    &scale1fbDown,   "scale1fbDown/F");    // event weight per 1/fb (MC)
    outTree->Branch("met",        &met,        "met/F");         // MET
    outTree->Branch("metPhi",     &metPhi,     "metPhi/F");      // phi(MET)
    outTree->Branch("sumEt",      &sumEt,      "sumEt/F");       // Sum ET
    outTree->Branch("u1",         &u1,         "u1/F");          // parallel component of recoil
    outTree->Branch("u2",         &u2,         "u2/F");          // perpendicular component of recoil
    outTree->Branch("tkMet",      &tkMet,      "tkMet/F");       // MET (track MET)
    outTree->Branch("tkMetPhi",   &tkMetPhi,   "tkMetPhi/F");    // phi(MET) (track MET)
    outTree->Branch("tkSumEt",    &tkSumEt,    "tkSumEt/F");     // Sum ET (track MET)
    outTree->Branch("tkU1",       &tkU1,       "tkU1/F");        // parallel component of recoil (track MET)
    outTree->Branch("tkU2",       &tkU2,       "tkU2/F");        // perpendicular component of recoil (track MET)
    outTree->Branch("mvaMet",     &mvaMet,     "mvaMet/F");      // MVA MET
    outTree->Branch("mvaMetPhi",  &mvaMetPhi,  "mvaMetPhi/F");   // phi(MVA MET)
    outTree->Branch("mvaSumEt",   &mvaSumEt,   "mvaSumEt/F");    // Sum ET (mva MET)
    outTree->Branch("mvaU1",      &mvaU1,      "mvaU1/F");       // parallel component of recoil (mva MET)
    outTree->Branch("mvaU2",      &mvaU2,      "mvaU2/F");       // perpendicular component of recoil (mva MET)
    outTree->Branch("puppiMet",    &puppiMet,   "puppiMet/F");      // Puppi MET
    outTree->Branch("puppiMetPhi", &puppiMetPhi,"puppiMetPhi/F");   // phi(Puppi MET)
    outTree->Branch("puppiSumEt",  &puppiSumEt, "puppiSumEt/F");    // Sum ET (Puppi MET)
    outTree->Branch("puppiU1",     &puppiU1,    "puppiU1/F");       // parallel component of recoil (Puppi MET)
    outTree->Branch("puppiU2",     &puppiU2,    "puppiU2/F");       // perpendicular component of recoil (Puppi MET)
    outTree->Branch("q1",         &q1,         "q1/I");          // charge of tag lepton
    outTree->Branch("q2",         &q2,         "q2/I");          // charge of probe lepton
    outTree->Branch("dilep",      "TLorentzVector",  &dilep);    // di-lepton 4-vector
    outTree->Branch("lep1",       "TLorentzVector",  &lep1);     // tag lepton 4-vector
    outTree->Branch("lep2",       "TLorentzVector",  &lep2);     // probe lepton 4-vector
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
    outTree->Branch("sc1",       "TLorentzVector",  &sc1);       // tag supercluster 4-vector
    outTree->Branch("sc2",       "TLorentzVector",  &sc2);       // probe supercluster 4-vector

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
      if(samp->jsonv[ifile].CompareTo("NONE")!=0) { 
	hasJSON = kTRUE;
	rlrm.addJSONFile(samp->jsonv[ifile].Data()); 
      }
  
      eventTree = (TTree*)infile->Get("Events");
      assert(eventTree);  
      eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Photon",   &scArr);       TBranch *scBr       = eventTree->GetBranch("Photon");
      eventTree->SetBranchAddress("PV",   &vertexArr);       TBranch *vertexBr = eventTree->GetBranch("PV");
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
	  totalWeight+=gen->weight*puWeight;
	  totalWeightUp+=gen->weight*puWeightUp;
	  totalWeightDown+=gen->weight*puWeightDown;
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
	
	// veto z -> xx decays for signal and z -> ee for bacground samples (needed for inclusive DYToLL sample)
	if (isWrongFlavor && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;
	else if (isSignal && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;

        // check for certified lumi (if applicable)
        baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        if(hasJSON && !rlrm.hasRunLumi(rl)) continue;  

        // trigger requirement
	if (!isEleTrigger(triggerMenu, info->triggerBits, isData)) continue;

        // good vertex requirement
        if(!(info->hasGoodPV)) continue;

	electronArr->Clear();
        electronBr->GetEntry(ientry);
	scArr->Clear();
	scBr->GetEntry(ientry);

	TLorentzVector vTag(0,0,0,0);
	TLorentzVector vTagSC(0,0,0,0);
	Double_t tagPt=0;
	Double_t Pt1=0;
	Double_t Pt2=0;
	Int_t itag=-1;
	Int_t tagscID=-1;
		
	for(Int_t i1=0; i1<electronArr->GetEntriesFast(); i1++) {
          const baconhep::TElectron *tag = (baconhep::TElectron*)((*electronArr)[i1]);
	  
	  // check ECAL gap
	  if(fabs(tag->scEta)>=ECAL_GAP_LOW && fabs(tag->scEta)<=ECAL_GAP_HIGH) continue;
	  
          // apply scale and resolution corrections to MC
          Double_t tagscEt_corr = tag->scEt;
          if(doScaleCorr && snamev[isam].CompareTo("data",TString::kIgnoreCase)!=0)
            tagscEt_corr = gRandom->Gaus(tag->scEt*getEleScaleCorr(tag->scEta,0),getEleResCorr(tag->scEta,0));
	  
	  if(tagscEt_corr        < PT_CUT)     continue;  // lepton pT cut
	  if(fabs(tag->scEta)    > ETA_CUT)    continue;  // lepton |eta| cut
	  if(!passEleID(tag,info->rhoIso))     continue;  // lepton selection
	
	  double El_Pt=0;
	  if(doScaleCorr) {
	    El_Pt=gRandom->Gaus(tag->pt*getEleScaleCorr(tag->scEta,0),getEleResCorr(tag->scEta,0));
	  }
	  else
	    {
	      El_Pt=tag->pt;
	    }

	  if(El_Pt>Pt1)
	    {
	      Pt2=Pt1;
	      Pt1=El_Pt;
	    }
	  else if(El_Pt>Pt2&&El_Pt<Pt1)
	    {
	      Pt2=El_Pt;
	    }

	  if(!isEleTriggerObj(triggerMenu, tag->hltMatchBits, kFALSE, isData)) continue;
	  
	  if(El_Pt<tagPt) continue;
	  
	  tagPt=El_Pt;
	  itag=i1;
	  tagscID=tag->scID;

	  // apply scale and resolution corrections to MC
          if(doScaleCorr && snamev[isam].CompareTo("data",TString::kIgnoreCase)!=0) {
            vTag.SetPtEtaPhiM(El_Pt, tag->eta, tag->phi, ELE_MASS);
            vTagSC.SetPtEtaPhiM(tagscEt_corr, tag->scEta, tag->scPhi, ELE_MASS);
          } else {
  	    vTag.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, ELE_MASS);
	    vTagSC.SetPtEtaPhiM(tag->scEt, tag->scEta, tag->scPhi, ELE_MASS);
          }

	  trkIso1    = tag->trkIso;
	  emIso1     = tag->ecalIso;
	  hadIso1    = tag->hcalIso;
	  pfChIso1   = tag->chHadIso;
	  pfGamIso1  = tag->gammaIso;	    
	  pfNeuIso1  = tag->neuHadIso;
	  pfCombIso1 = tag->chHadIso + TMath::Max(tag->neuHadIso + tag->gammaIso - (info->rhoIso)*getEffAreaEl(tag->scEta), 0.);
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
	  q1         = tag->q;

	}

	if(tagPt<Pt2) continue;

	TLorentzVector vProbe(0,0,0,0); TLorentzVector vProbeSC(0,0,0,0);
	Double_t probePt=0;
	Int_t iprobe=-1;
	Int_t passID=false;
	UInt_t icat=0;
	
	const baconhep::TElectron *eleProbe=0;

	for(Int_t j=0; j<scArr->GetEntriesFast(); j++) {
	  const baconhep::TPhoton *scProbe = (baconhep::TPhoton*)((*scArr)[j]);
	  
	  if(scProbe->scID == tagscID) continue;

	  // check ECAL gap
	  if(fabs(scProbe->eta)>=ECAL_GAP_LOW && fabs(scProbe->eta)<=ECAL_GAP_HIGH) continue;
	  
	  // apply scale and resolution corrections to MC
	  Double_t scProbept_corr = scProbe->pt;
	  if(doScaleCorr && snamev[isam].CompareTo("data",TString::kIgnoreCase)!=0)
	    scProbept_corr = gRandom->Gaus(scProbe->pt*getEleScaleCorr(scProbe->eta,0),getEleResCorr(scProbe->eta,0));
	  
	  if(scProbept_corr        < PT_CUT)  continue;  // Supercluster ET cut ("pt" = corrected by PV position)
	  if(fabs(scProbe->eta)  > ETA_CUT) continue;  // Supercluster |eta| cuts

	  for(Int_t i2=0; i2<electronArr->GetEntriesFast(); i2++) {
	    if(itag==i2) continue;
	    const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[i2]);
	    if(!(ele->typeBits & baconhep::EEleType::kEcalDriven)) continue;
	    if(scProbe->scID==ele->scID) { 
	      eleProbe = ele; 
	      iprobe   = i2;
	      break; 
	    }
	  }

	  double El_Pt=0;
	  if(doScaleCorr&&eleProbe) {
	    El_Pt=gRandom->Gaus(eleProbe->pt*getEleScaleCorr(scProbe->eta,0),getEleResCorr(scProbe->eta,0));
	  }
	  else if(!doScaleCorr&&eleProbe)
	    {
	      El_Pt=eleProbe->pt;
	    }
	  else
	    {
	      El_Pt=scProbept_corr;
	    }

	  if(passID&&eleProbe&&passEleID(eleProbe,info->rhoIso)&&El_Pt<probePt) continue;
	  if(passID&&eleProbe&&!passEleID(eleProbe,info->rhoIso)) continue;
	  if(passID&&!eleProbe) continue;
	  if(!passID&&eleProbe&&!passEleID(eleProbe,info->rhoIso)&&El_Pt<probePt) continue;
	  if(!passID&&!eleProbe&&El_Pt<probePt) continue;
	  if(!passID&&eleProbe&&passEleID(eleProbe,info->rhoIso)) passID=true;

	  probePt=El_Pt;

	  // apply scale and resolution corrections to MC
	  if(doScaleCorr && snamev[isam].CompareTo("data",TString::kIgnoreCase)!=0) {
	    vProbe.SetPtEtaPhiM((eleProbe) ? gRandom->Gaus(eleProbe->pt*getEleScaleCorr(scProbe->eta,0),getEleResCorr(scProbe->eta,0)) : scProbept_corr,
				(eleProbe) ? eleProbe->eta : scProbe->eta,
				(eleProbe) ? eleProbe->phi : scProbe->phi,
				ELE_MASS);
	    vProbeSC.SetPtEtaPhiM((eleProbe) ? gRandom->Gaus(eleProbe->scEt*getEleScaleCorr(scProbe->eta,0),getEleResCorr(scProbe->eta,0)) : gRandom->Gaus(scProbe->pt*getEleScaleCorr(scProbe->eta,0),getEleResCorr(scProbe->eta,0)),
				  scProbe->eta, scProbe->phi, ELE_MASS);
	  } else {
	    vProbe.SetPtEtaPhiM((eleProbe) ? eleProbe->pt : scProbe->pt,
				(eleProbe) ? eleProbe->eta : scProbe->eta,
				(eleProbe) ? eleProbe->phi : scProbe->phi,
				ELE_MASS);
	    vProbeSC.SetPtEtaPhiM((eleProbe) ? eleProbe->scEt : scProbe->pt,
				  scProbe->eta, scProbe->phi, ELE_MASS);
	  }

	  trkIso2    = (eleProbe) ? eleProbe->trkIso        : -1;
	  emIso2     = (eleProbe) ? eleProbe->ecalIso       : -1;
	  hadIso2    = (eleProbe) ? eleProbe->hcalIso       : -1;
	  pfChIso2   = (eleProbe) ? eleProbe->chHadIso      : -1;
	  pfGamIso2  = (eleProbe) ? eleProbe->gammaIso      : -1;
	  pfNeuIso2  = (eleProbe) ? eleProbe->neuHadIso     : -1;	    
	  pfCombIso2 = (eleProbe) ? 
	    eleProbe->chHadIso + TMath::Max(eleProbe->neuHadIso + eleProbe->gammaIso - 
					    (info->rhoIso)*getEffAreaEl(eleProbe->scEta), 0.) :  -1;
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
	  q2         = (eleProbe) ? eleProbe->q : -q1;

	  // determine event category
	  if(eleProbe) {
	    if(passEleID(eleProbe,info->rhoIso)) {
	      
	      if(isEleTriggerObj(triggerMenu, eleProbe->hltMatchBits, kFALSE, isData)) {
		icat=eEleEle2HLT;  
	      } 
	      else if(isEleTriggerObj(triggerMenu, eleProbe->hltMatchBits, kTRUE, isData)) {
		icat=eEleEle1HLT1L1; 
	      }
	      else { icat=eEleEle1HLT; }
	    }
	    else { icat=eEleEleNoSel; } 
	  } 
	  else { icat=eEleSC; }
	  
	}

	if(q1 == q2)         continue;  // opposite charge requirement

	// mass window
	TLorentzVector vDilep = vTag + vProbe;
	if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;

	if(icat==0) continue;

	//******** We have a Z candidate! HURRAY! ********
	nsel+=weight;
	nselvar+=weight*weight;

	// Perform matching of dileptons to GEN leptons from Z decay

	Int_t glepq1=-99;
	Int_t glepq2=-99;
	TLorentzVector *gvec=new TLorentzVector(0,0,0,0);
	TLorentzVector *glep1=new TLorentzVector(0,0,0,0);
	TLorentzVector *glep2=new TLorentzVector(0,0,0,0);
	TLorentzVector *gph=new TLorentzVector(0,0,0,0);
	Bool_t hasGenMatch = kFALSE;
	if(isSignal && hasGen) {
	  toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2,&glepq1,&glepq2,1);
	  
	  Bool_t match1 = ( ((glep1) && toolbox::deltaR(vTag.Eta(), vTag.Phi(), glep1->Eta(), glep1->Phi())<0.3) || 
			    ((glep2) && toolbox::deltaR(vTag.Eta(), vTag.Phi(), glep2->Eta(), glep2->Phi())<0.3) );
	  
	  Bool_t match2 = ( ((glep1) && toolbox::deltaR(vProbe.Eta(), vProbe.Phi(), glep1->Eta(), glep1->Phi())<0.3) || 
			    ((glep2) && toolbox::deltaR(vProbe.Eta(), vProbe.Phi(), glep2->Eta(), glep2->Phi())<0.3) );
	  
	  if(match1 && match2) {
	    hasGenMatch = kTRUE;
	    if (gvec!=0) {
	      genV=new TLorentzVector(0,0,0,0);
	      genV->SetPtEtaPhiM(gvec->Pt(), gvec->Eta(), gvec->Phi(), gvec->M());
	      genVPt   = gvec->Pt();
	      genVPhi  = gvec->Phi();
	      genVy    = gvec->Rapidity();
	      genVMass = gvec->M();
	    }
	    else {
	      TLorentzVector tvec=*glep1+*glep2;
	      genV=new TLorentzVector(0,0,0,0);
	      genV->SetPtEtaPhiM(tvec.Pt(), tvec.Eta(), tvec.Phi(), tvec.M());
	      genVPt   = tvec.Pt();
	      genVPhi  = tvec.Phi();
	      genVy    = tvec.Rapidity();
	      genVMass = tvec.M();
	    }
	    delete gvec;
	    delete glep1;
	    delete glep2;
	    glep1=0; glep2=0; gvec=0;
	  }
	  else {
	    genV     = new TLorentzVector(0,0,0,0);
	    genVPt   = -999;
	    genVPhi  = -999;
	    genVy    = -999;
	    genVMass = -999;
	  }
	}
	
	if (hasGen) {
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

	if (hasGenMatch) matchGen=1;
	else matchGen=0;

	category = icat;

	vertexArr->Clear();
	vertexBr->GetEntry(ientry);

	npv      = vertexArr->GetEntries();
	npu      = info->nPUmean;
	scale1fb = weight;
	scale1fbUp = weightUp;
	scale1fbDown = weightDown;
	met      = info->pfMETC;
	metPhi   = info->pfMETCphi;
	sumEt    = 0;
	tkMet    = info->trkMET;
	tkMetPhi = info->trkMETphi;
	tkSumEt  = 0;
	mvaMet   = info->mvaMET;
	mvaMetPhi = info->mvaMETphi; 
	mvaSumEt = 0;
	TVector2 vZPt((vDilep.Pt())*cos(vDilep.Phi()),(vDilep.Pt())*sin(vDilep.Phi()));

	puppiMet = info->puppET;
        puppiMetPhi = info->puppETphi;
	puppiSumEt = 0;
	lep1     = &vTag;
	lep2     = &vProbe;
	dilep    = &vDilep;
	sc1        = &vTagSC;
	sc2        = &vProbeSC;

	TVector2 vMet((info->pfMETC)*cos(info->pfMETCphi), (info->pfMETC)*sin(info->pfMETCphi));
	TVector2 vU = -1.0*(vMet+vZPt);
	u1 = ((vDilep.Px())*(vU.Px()) + (vDilep.Py())*(vU.Py()))/(vDilep.Pt());  // u1 = (pT . u)/|pT|
	u2 = ((vDilep.Px())*(vU.Py()) - (vDilep.Py())*(vU.Px()))/(vDilep.Pt());  // u2 = (pT x u)/|pT|
	
	TVector2 vTkMet((info->trkMET)*cos(info->trkMETphi), (info->trkMET)*sin(info->trkMETphi));        
	TVector2 vTkU = -1.0*(vTkMet+vZPt);
	tkU1 = ((vDilep.Px())*(vTkU.Px()) + (vDilep.Py())*(vTkU.Py()))/(vDilep.Pt());  // u1 = (pT . u)/|pT|
	tkU2 = ((vDilep.Px())*(vTkU.Py()) - (vDilep.Py())*(vTkU.Px()))/(vDilep.Pt());  // u2 = (pT x u)/|pT|
	
	TVector2 vMvaMet((info->mvaMET)*cos(info->mvaMETphi), (info->mvaMET)*sin(info->mvaMETphi));
	TVector2 vMvaU = -1.0*(vMvaMet+vZPt);
	mvaU1 = ((vDilep.Px())*(vMvaU.Px()) + (vDilep.Py())*(vMvaU.Py()))/(vDilep.Pt());  // u1 = (pT . u)/|pT|
	mvaU2 = ((vDilep.Px())*(vMvaU.Py()) - (vDilep.Py())*(vMvaU.Px()))/(vDilep.Pt());  // u2 = (pT x u)/|pT|
        
	TVector2 vPuppiMet((info->puppET)*cos(info->puppETphi), (info->puppET)*sin(info->puppETphi));
	TVector2 vPuppiU = -1.0*(vPuppiMet+vZPt);
	puppiU1 = ((vDilep.Px())*(vPuppiU.Px()) + (vDilep.Py())*(vPuppiU.Py()))/(vDilep.Pt());  // u1 = (pT . u)/|pT|
	puppiU2 = ((vDilep.Px())*(vPuppiU.Py()) - (vDilep.Py())*(vPuppiU.Px()))/(vDilep.Pt());  // u2 = (pT x u)/|pT|

	outTree->Fill();
	delete genV;
	genV=0, dilep=0, lep1=0, lep2=0, sc1=0, sc2=0;
      }
      delete infile;
      infile=0, eventTree=0;    
      
      cout << nsel  << " +/- " << sqrt(nselvar);
      if(!isData) cout << " per 1/fb";
      cout << endl;
    }
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
