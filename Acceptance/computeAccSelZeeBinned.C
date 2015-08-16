//================================================================================================
//
// Compute Z->ee acceptance at full selection level
//
//  * outputs results summary text file
//
//  [!!!] propagation of efficiency scale factor uncertainties no yet implemented
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TH1D.h>                   // histogram class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "TLorentzVector.h"         // 4-vector class
#include <TRandom3.h>

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions

// helper class to handle efficiency tables
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"
#endif

//=== MAIN MACRO ================================================================================================= 

void computeAccSelZeeBinned(const TString conf,            // input file
                            const TString outputDir        // output directory
) {
  gBenchmark->Start("computeAccSelZeeBinned");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.5;
  const Double_t ELE_MASS   = 0.000511;

  const Double_t ETA_BARREL = 1.4442;
  const Double_t ETA_ENDCAP = 1.566;

  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 11;
  
  // efficiency files
  const TString dataHLTEffName     = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/DataZee_EleHLTEff/eff.root";
  const TString dataHLTEffName_pos = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/DataZee_EleHLTEff/eff.root";
  const TString dataHLTEffName_neg = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/DataZee_EleHLTEff/eff.root";

  const TString zeeHLTEffName      = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/Zee_EleHLTEff/eff.root";
  const TString zeeHLTEffName_pos  = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/Zee_EleHLTEff/eff.root";
  const TString zeeHLTEffName_neg  = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/Zee_EleHLTEff/eff.root";
  
  const TString dataGsfSelEffName     = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/DataZee_EleGsfSelEff/eff.root";
  const TString dataGsfSelEffName_pos = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/DataZee_EleGsfSelEff/eff.root";
  const TString dataGsfSelEffName_neg = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/DataZee_EleGsfSelEff/eff.root";

  const TString zeeGsfSelEffName      = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/Zee_EleGsfSelEff/eff.root";
  const TString zeeGsfSelEffName_pos  = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/Zee_EleGsfSelEff/eff.root";
  const TString zeeGsfSelEffName_neg  = "/data/blue/cmedlock/wz-efficiency-results-coarsebinning/Zee_EleGsfSelEff/eff.root";

  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/pileup_weights_2015B.root", "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("npv_rw");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString> fnamev;  // file name per input file
  vector<TString> labelv;  // TLegend label per input file
  vector<Int_t>   colorv;  // plot color per input file
  vector<Int_t>   linev;   // plot line style per input file

  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    
    string fname;
    Int_t color, linesty;
    stringstream ss(line);
    ss >> fname >> color >> linesty;
    string label = line.substr(line.find('@')+1);
    fnamev.push_back(fname);
    labelv.push_back(label);
    colorv.push_back(color);
    linev.push_back(linesty);
  }
  ifs.close();

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  
  TH2D *h=0;
  
  //
  // HLT efficiency
  //
  cout << "Loading trigger efficiencies..." << endl;

  TFile *dataHLTEffFile_pos = new TFile(dataHLTEffName_pos);
  CEffUser2D dataHLTEff_pos;
  dataHLTEff_pos.loadEff((TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEffFile_neg = new TFile(dataHLTEffName_neg);
  CEffUser2D dataHLTEff_neg;
  dataHLTEff_neg.loadEff((TH2D*)dataHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zeeHLTEffFile_pos = new TFile(zeeHLTEffName_pos);
  CEffUser2D zeeHLTEff_pos;
  zeeHLTEff_pos.loadEff((TH2D*)zeeHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)zeeHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zeeHLTEffFile_neg = new TFile(zeeHLTEffName_neg);
  CEffUser2D zeeHLTEff_neg;
  zeeHLTEff_neg.loadEff((TH2D*)zeeHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)zeeHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeHLTEffFile_neg->Get("hErrhEtaPt"));
  
  h =(TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt");
  TH2D *hHLTErr_pos = new TH2D("hHLTErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                                 h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hHLTErr_neg = new TH2D("hHLTErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                                 h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  
  //
  // Selection efficiency
  //
  cout << "Loading GSF+selection efficiencies..." << endl;
  
  TFile *dataGsfSelEffFile_pos = new TFile(dataGsfSelEffName_pos);
  CEffUser2D dataGsfSelEff_pos;
  dataGsfSelEff_pos.loadEff((TH2D*)dataGsfSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataGsfSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataGsfSelEffFile_neg = new TFile(dataGsfSelEffName_neg);
  CEffUser2D dataGsfSelEff_neg;
  dataGsfSelEff_neg.loadEff((TH2D*)dataGsfSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataGsfSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEffFile_neg->Get("hErrhEtaPt"));

  TFile *zeeGsfSelEffFile_pos = new TFile(zeeGsfSelEffName_pos);
  CEffUser2D zeeGsfSelEff_pos;
  zeeGsfSelEff_pos.loadEff((TH2D*)zeeGsfSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEffFile_pos->Get("hErrhEtaPt"));

  TFile *zeeGsfSelEffFile_neg = new TFile(zeeGsfSelEffName_neg);
  CEffUser2D zeeGsfSelEff_neg;
  zeeGsfSelEff_neg.loadEff((TH2D*)zeeGsfSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEffFile_neg->Get("hErrhEtaPt"));
 
  h =(TH2D*)dataGsfSelEffFile_pos->Get("hEffEtaPt");
  TH2D *hGsfSelErr_pos = new TH2D("hGsfSelErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                                       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hGsfSelErr_neg = new TH2D("hGsfSelErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                                       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo   *info = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr     = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *electronArr    = new TClonesArray("baconhep::TElectron");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");

  TFile *infile=0;
  TTree *eventTree=0;
  
  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv;
  vector<Double_t> nSelCorrv, nSelCorrVarv;
  vector<Double_t> accv, accCorrv;
  vector<Double_t> accErrv, accErrCorrv;

  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = TFile::Open(fnamev[ifile]); 
    assert(infile);

    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info",              &info); TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("GenEvtInfo",         &gen); TBranch *genBr      = eventTree->GetBranch("GenEvtInfo");
    eventTree->SetBranchAddress("GenParticle", &genPartArr); TBranch *genPartBr  = eventTree->GetBranch("GenParticle");
    eventTree->SetBranchAddress("Electron",   &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");

    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelCorrv.push_back(0);
    nSelCorrVarv.push_back(0);

    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    //for(UInt_t ientry=0; ientry<100; ientry++) {
      genBr->GetEntry(ientry);
      genPartArr->Clear(); genPartBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);

      Double_t weight=gen->weight/23443.4;
      //Double_t weight=gen->weight;
      //std::cout << gen->weight << std::endl;
      //Double_t weight=1.0;

      if (fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      TLorentzVector *vec=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,1);
      if(vec->M()<MASS_LOW || vec->M()>MASS_HIGH) continue;      
      delete vec; delete lep1; delete lep2;

      vertexArr->Clear();
      vertexBr->GetEntry(ientry);
      double npv  = vertexArr->GetEntries();
      weight*=h_rw->GetBinContent(npv+1);
      nEvtsv[ifile]+=weight;        
      // trigger requirement               
      if (!isEleTrigger(triggerMenu, info->triggerBits, kFALSE)) continue;
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
      
      electronArr->Clear();
      electronBr->GetEntry(ientry);

      for(Int_t i1=0; i1<electronArr->GetEntriesFast(); i1++) {
  	const baconhep::TElectron *ele1 = (baconhep::TElectron*)((*electronArr)[i1]);
	
	// check ECAL gap
	if(fabs(ele1->scEta)>=ETA_BARREL && fabs(ele1->scEta)<=ETA_ENDCAP) continue;
	
	double ele1_pt = gRandom->Gaus(ele1->pt*getEleScaleCorr(ele1->scEta,0), getEleResCorr(ele1->scEta,0));

	if(ele1_pt	     < PT_CUT)	  continue;  // lepton pT cut
        if(fabs(ele1->scEta) > ETA_CUT)	  continue;  // lepton |eta| cut
        if(!passEleID(ele1,info->rhoIso)) continue;  // lepton selection
	//if(!isEleTriggerObj(triggerMenu, ele1->hltMatchBits, kFALSE, kFALSE)) continue;

        TLorentzVector vEle1(0,0,0,0);
	vEle1.SetPtEtaPhiM(ele1_pt, ele1->eta, ele1->phi, ELE_MASS);
	Bool_t isB1 = (fabs(ele1->scEta)<ETA_BARREL) ? kTRUE : kFALSE;

        for(Int_t i2=i1+1; i2<electronArr->GetEntriesFast(); i2++) {          
	  const baconhep::TElectron *ele2 = (baconhep::TElectron*)((*electronArr)[i2]);
	  
	  // check ECAL gap
	  if(fabs(ele2->scEta)>=ETA_BARREL && fabs(ele2->scEta)<=ETA_ENDCAP) continue;
	  double ele2_pt = gRandom->Gaus(ele2->scEt*getEleScaleCorr(ele2->scEta,0), getEleResCorr(ele2->scEta,0));
        
          if(ele2_pt        < PT_CUT)    continue;  // lepton pT cut
          if(fabs(ele2->scEta) > ETA_CUT)   continue;  // lepton |eta| cut
	  if(!passEleID(ele2,info->rhoIso)) continue;  // lepton selection

          TLorentzVector vEle2(0,0,0,0);
	  vEle2.SetPtEtaPhiM(ele2_pt, ele2->eta, ele2->phi, ELE_MASS);  
          Bool_t isB2 = (fabs(ele2->scEta)<ETA_BARREL) ? kTRUE : kFALSE;

	  if(!isEleTriggerObj(triggerMenu, ele1->hltMatchBits, kFALSE, kFALSE) && !isEleTriggerObj(triggerMenu, ele2->hltMatchBits, kFALSE, kFALSE)) continue;
	  
	  // mass window
          TLorentzVector vDilep = vEle1 + vEle2;
          if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
          
          /******** We have a Z candidate! HURRAY! ********/
          Double_t effdata, effmc;
          //Double_t sceta1 = (fabs(ele1->scEta)<2.5) ? ele1->scEta : 0.99*(ele1->scEta);
	  Double_t sceta1 = ele1->scEta;
          //Double_t sceta2 = (fabs(ele2->scEta)<2.5) ? ele2->scEta : 0.99*(ele2->scEta);
    	  Double_t sceta2 = ele2->scEta;

          Double_t corr=1;
	  
	  effdata=1; effmc=1;
          if(ele1->q>0) { 
            effdata *= (1.-dataHLTEff_pos.getEff(sceta1, ele1->scEt));
            effmc   *= (1.-zeeHLTEff_pos.getEff(sceta1, ele1->scEt));
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(sceta1, ele1->scEt)); 
            effmc   *= (1.-zeeHLTEff_neg.getEff(sceta1, ele1->scEt)); 
          }
          if(ele2->q>0) {
            effdata *= (1.-dataHLTEff_pos.getEff(sceta2, ele2->scEt)); 
            effmc   *= (1.-zeeHLTEff_pos.getEff(sceta2, ele2->scEt));
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(sceta2, ele2->scEt)); 
            effmc   *= (1.-zeeHLTEff_neg.getEff(sceta2, ele2->scEt));
          }
          effdata = 1.-effdata;
          effmc   = 1.-effmc;
          corr *= effdata/effmc;


          effdata=1; effmc=1;
          if(ele1->q>0) { 
            effdata *= dataGsfSelEff_pos.getEff(sceta1, ele1->scEt); 
            effmc   *= zeeGsfSelEff_pos.getEff(sceta1, ele1->scEt); 
          } else {
            effdata *= dataGsfSelEff_neg.getEff(sceta1, ele1->scEt); 
            effmc   *= zeeGsfSelEff_neg.getEff(sceta1, ele1->scEt); 
          }
          if(ele2->q>0) {
            effdata *= dataGsfSelEff_pos.getEff(sceta2, ele2->scEt); 
            effmc   *= zeeGsfSelEff_pos.getEff(sceta2, ele2->scEt);
          } else {
            effdata *= dataGsfSelEff_neg.getEff(sceta2, ele2->scEt); 
            effmc   *= zeeGsfSelEff_neg.getEff(sceta2, ele2->scEt);
          }
          corr *= effdata/effmc;
	  
	  // scale factor uncertainties
	  if(ele1->q>0) {	    
	    Double_t effdata = dataGsfSelEff_pos.getEff(sceta1, ele1->scEt);
	    Double_t errdata = TMath::Max(dataGsfSelEff_pos.getErrLow(sceta1, ele1->scEt), dataGsfSelEff_pos.getErrHigh(sceta1, ele1->scEt));
            Double_t effmc   = zeeGsfSelEff_pos.getEff(sceta1, ele1->scEt); 
	    Double_t errmc   = TMath::Max(zeeGsfSelEff_pos.getErrLow(sceta1, ele1->scEt), zeeGsfSelEff_pos.getErrHigh(sceta1, ele1->scEt));
	    Double_t errGsfSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    hGsfSelErr_pos->Fill(sceta1, ele1->scEt, errGsfSel);
	  } else {
	    Double_t effdata = dataGsfSelEff_neg.getEff(sceta1, ele1->scEt);
	    Double_t errdata = TMath::Max(dataGsfSelEff_neg.getErrLow(sceta1, ele1->scEt), dataGsfSelEff_neg.getErrHigh(sceta1, ele1->scEt));
            Double_t effmc   = zeeGsfSelEff_neg.getEff(sceta1, ele1->scEt); 
	    Double_t errmc   = TMath::Max(zeeGsfSelEff_neg.getErrLow(sceta1, ele1->scEt), zeeGsfSelEff_neg.getErrHigh(sceta1, ele1->scEt));
	    Double_t errGsfSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    hGsfSelErr_neg->Fill(sceta1, ele1->scEt, errGsfSel);
	  }

	  if(ele2->q>0) {	    
	    Double_t effdata = dataHLTEff_pos.getEff(sceta2, ele2->scEt);
	    Double_t errdata = TMath::Max(dataHLTEff_pos.getErrLow(sceta2, ele2->scEt), dataHLTEff_pos.getErrHigh(sceta2, ele2->scEt));
            Double_t effmc   = zeeHLTEff_pos.getEff(sceta2, ele2->scEt); 
	    Double_t errmc   = TMath::Max(zeeHLTEff_pos.getErrLow(sceta2, ele2->scEt), zeeHLTEff_pos.getErrHigh(sceta2, ele2->scEt));
	    Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    hHLTErr_pos->Fill(sceta2, ele2->scEt, errHLT);
	  } else {
	    Double_t effdata = dataHLTEff_neg.getEff(sceta2, ele2->scEt);
	    Double_t errdata = TMath::Max(dataHLTEff_neg.getErrLow(sceta2, ele2->scEt), dataHLTEff_neg.getErrHigh(sceta2, ele2->scEt));
            Double_t effmc   = zeeHLTEff_neg.getEff(sceta2, ele2->scEt); 
	    Double_t errmc   = TMath::Max(zeeHLTEff_neg.getErrLow(sceta2, ele2->scEt), zeeHLTEff_neg.getErrHigh(sceta2, ele2->scEt));
	    Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    hHLTErr_neg->Fill(sceta2, ele2->scEt, errHLT);
	  }

	  nSelv[ifile]+=weight;
	  nSelCorrv[ifile]+=weight*corr;
	  nSelCorrVarv[ifile]+=weight*weight*corr*corr;
	}
      }
    }
    
    Double_t var=0;
    for(Int_t iy=0; iy<=hHLTErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr_pos->GetNbinsX(); ix++) {
	Double_t err=hHLTErr_pos->GetBinContent(ix,iy);
	var+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hHLTErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr_neg->GetNbinsX(); ix++) {
	Double_t err=hHLTErr_neg->GetBinContent(ix,iy);
	var+=err*err;
      }
    }
    cout << "var1: " << var << endl;
    for(Int_t iy=0; iy<=hGsfSelErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hGsfSelErr_pos->GetNbinsX(); ix++) {
	Double_t err=hGsfSelErr_pos->GetBinContent(ix,iy);
	var+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hGsfSelErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hGsfSelErr_neg->GetNbinsX(); ix++) {
	Double_t err=hGsfSelErr_neg->GetBinContent(ix,iy);
	var+=err*err;
      }
    }
    cout << "var2: " << var << endl;

    nSelCorrVarv[ifile]+=var;

    // compute acceptances
    std::cout << nEvtsv[ifile] << " " << nSelv[ifile] << std::endl;
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);         accErrv.push_back(sqrt(accv[ifile]*(1. +accv[ifile])/nEvtsv[ifile]));
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]); accErrCorrv.push_back(accCorrv[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrv[ifile]*nSelCorrv[ifile]) + 1./nEvtsv[ifile]));
    
    delete infile;
    infile=0, eventTree=0;  
  }  
  delete info;
  delete gen;
  delete electronArr;  
    
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
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "   ================================================" << endl;
    cout << "    Label: " << labelv[ifile] << endl;
    cout << "     File: " << fnamev[ifile] << endl;
    cout << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "          nominal: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    cout << "     SF corrected: " << accCorrv[ifile] << " +/- " << accErrCorrv[ifile] << endl;
    cout << endl;
  }
  
  char txtfname[100];
  sprintf(txtfname,"%s/binned.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << " Z -> e e" << endl;
  txtfile << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  txtfile << "  pT > " << PT_CUT << endl;
  txtfile << "  |eta| < " << ETA_CUT << endl;
  txtfile << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "   ================================================" << endl;
    txtfile << "    Label: " << labelv[ifile] << endl;
    txtfile << "     File: " << fnamev[ifile] << endl;
    txtfile << endl;
    txtfile << "    *** Acceptance ***" << endl;
    txtfile << "          nominal: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    txtfile << "     SF corrected: " << accCorrv[ifile] << " +/- " << accErrCorrv[ifile] << endl;
    txtfile << endl;
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelZeeBinned"); 
}
