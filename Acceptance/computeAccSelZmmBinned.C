//================================================================================================
//
// Compute Z->mumu acceptance at full selection level
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
#include "Math/LorentzVector.h"     // 4-vector class

// define structures to read in ntuple
#include "../Ntupler/interface/EWKAnaDefs.hh"
#include "../Ntupler/interface/TEventInfo.hh"
#include "../Ntupler/interface/TGenInfo.hh"
#include "../Ntupler/interface/TMuon.hh"

// helper functions for lepton ID selection
#include "../Utils/LeptonIDCuts.hh"

// helper class to handle efficiency tables
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void computeAccSelZmmBinned(const TString conf,      // input file
                            const TString outputDir  // output directory
) {
  gBenchmark->Start("computeAccSelZmmBinned");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.1;
  const Double_t MUON_MASS  = 0.105658369;
  
  // efficiency files
  const TString dataHLTEffName_pos = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_MuHLTEff_pos/analysis/eff.root";
  const TString dataHLTEffName_neg = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_MuHLTEff_neg/analysis/eff.root";
  const TString zmmHLTEffName_pos  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zmm_MuHLTEff_pos/analysis/eff.root";
  const TString zmmHLTEffName_neg  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zmm_MuHLTEff_neg/analysis/eff.root";
  
  const TString dataSelEffName_pos = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_MuSelEff_pos/analysis/eff.root";
  const TString dataSelEffName_neg = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_MuSelEff_neg/analysis/eff.root";
  const TString zmmSelEffName_pos  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zmm_MuSelEff_pos/analysis/eff.root";
  const TString zmmSelEffName_neg  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zmm_MuSelEff_neg/analysis/eff.root";
  
  const TString dataTrkEffName_pos = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_MuTrkEff_pos/analysis/eff.root";
  const TString dataTrkEffName_neg = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_MuTrkEff_neg/analysis/eff.root";
  const TString zmmTrkEffName_pos  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zmm_MuTrkEff_pos/analysis/eff.root";
  const TString zmmTrkEffName_neg  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zmm_MuTrkEff_neg/analysis/eff.root";
  
  const TString dataStaEffName_pos = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_MuStaEff_iso_pos/analysis/eff.root";
  const TString dataStaEffName_neg = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_MuStaEff_iso_neg/analysis/eff.root";
  const TString zmmStaEffName_pos  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zmm_MuStaEff_iso_pos/analysis/eff.root";
  const TString zmmStaEffName_neg  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zmm_MuStaEff_iso_neg/analysis/eff.root";
  
  
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
    
  TFile *zmmHLTEffFile_pos = new TFile(zmmHLTEffName_pos);
  CEffUser2D zmmHLTEff_pos;
  zmmHLTEff_pos.loadEff((TH2D*)zmmHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmHLTEffFile_neg = new TFile(zmmHLTEffName_neg);
  CEffUser2D zmmHLTEff_neg;
  zmmHLTEff_neg.loadEff((TH2D*)zmmHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmHLTEffFile_neg->Get("hErrhEtaPt"));
  
  //
  // Selection efficiency
  //
  cout << "Loading selection efficiencies..." << endl;
  
  TFile *dataSelEffFile_pos = new TFile(dataSelEffName_pos);
  CEffUser2D dataSelEff_pos;
  dataSelEff_pos.loadEff((TH2D*)dataSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataSelEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataSelEffFile_neg = new TFile(dataSelEffName_neg);
  CEffUser2D dataSelEff_neg;
  dataSelEff_neg.loadEff((TH2D*)dataSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataSelEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmSelEffFile_pos = new TFile(zmmSelEffName_pos);
  CEffUser2D zmmSelEff_pos;
  zmmSelEff_pos.loadEff((TH2D*)zmmSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmSelEffFile_pos->Get("hErrhEtaPt"));

  TFile *zmmSelEffFile_neg = new TFile(zmmSelEffName_neg);
  CEffUser2D zmmSelEff_neg;
  zmmSelEff_neg.loadEff((TH2D*)zmmSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmSelEffFile_neg->Get("hErrhEtaPt"));
 
  //
  // Standalone efficiency
  //
  cout << "Loading standalone efficiencies..." << endl;
  
  TFile *dataStaEffFile_pos = new TFile(dataStaEffName_pos);
  CEffUser2D dataStaEff_pos;
  dataStaEff_pos.loadEff((TH2D*)dataStaEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataStaEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataStaEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataStaEffFile_neg = new TFile(dataStaEffName_neg);
  CEffUser2D dataStaEff_neg;
  dataStaEff_neg.loadEff((TH2D*)dataStaEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataStaEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataStaEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmStaEffFile_pos = new TFile(zmmStaEffName_pos);
  CEffUser2D zmmStaEff_pos;
  zmmStaEff_pos.loadEff((TH2D*)zmmStaEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmStaEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmStaEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmStaEffFile_neg = new TFile(zmmStaEffName_neg);
  CEffUser2D zmmStaEff_neg;
  zmmStaEff_neg.loadEff((TH2D*)zmmStaEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmStaEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmStaEffFile_neg->Get("hErrhEtaPt"));
 
  //
  // Tracker efficiency
  //
  cout << "Loading track efficiencies..." << endl;
  
  TFile *dataTrkEffFile_pos = new TFile(dataTrkEffName_pos);
  CEffUser2D dataTrkEff_pos;
  dataTrkEff_pos.loadEff((TH2D*)dataTrkEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataTrkEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataTrkEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataTrkEffFile_neg = new TFile(dataTrkEffName_neg);
  CEffUser2D dataTrkEff_neg;
  dataTrkEff_neg.loadEff((TH2D*)dataTrkEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataTrkEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataTrkEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmTrkEffFile_pos = new TFile(zmmTrkEffName_pos);
  CEffUser2D zmmTrkEff_pos;
  zmmTrkEff_pos.loadEff((TH2D*)zmmTrkEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmTrkEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmTrkEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmTrkEffFile_neg = new TFile(zmmTrkEffName_neg);
  CEffUser2D zmmTrkEff_neg;
  zmmTrkEff_neg.loadEff((TH2D*)zmmTrkEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmTrkEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmTrkEffFile_neg->Get("hErrhEtaPt"));


  // Data structures to store info from TTrees
  mithep::TEventInfo *info = new mithep::TEventInfo();
  mithep::TGenInfo   *gen  = new mithep::TGenInfo();
  TClonesArray *muonArr    = new TClonesArray("mithep::TMuon");
  
  TFile *infile=0;
  TTree *eventTree=0;
   
  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv, nSelCorrv;
  vector<Double_t> accv, accCorrv;
  vector<Double_t> accErrv, accCorrErrv;
  
  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);
  
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Gen",  &gen);     TBranch *genBr  = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");   

    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelCorrv.push_back(0);
    
    //
    // loop over events
    //      
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      genBr->GetEntry(ientry);
      if(gen->vmass<MASS_LOW || gen->vmass>MASS_HIGH) continue;
   
      infoBr->GetEntry(ientry);
    
      Double_t weight=1;
      nEvtsv[ifile]+=weight;
      
      // trigger requirement               
      ULong_t trigger = kHLT_Mu15_eta2p1;
      ULong_t trigObj = kHLT_Mu15_eta2p1_MuObj;   
      if(!(info->triggerBits & trigger)) continue;  
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
    
      muonArr->Clear();
      muonBr->GetEntry(ientry);
      for(Int_t i1=0; i1<muonArr->GetEntriesFast(); i1++) {
  	const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[i1]);

        if(mu1->pt	  < PT_CUT)  continue;  // lepton pT cut
        if(fabs(mu1->eta) > ETA_CUT) continue;  // lepton |eta| cut
        if(!passMuonID(mu1))	     continue;  // lepton selection
	
	LorentzVector vMu1(mu1->pt, mu1->eta, mu1->phi, MUON_MASS);

        for(Int_t i2=i1+1; i2<muonArr->GetEntriesFast(); i2++) {
          const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[i2]);
        
          if(mu1->q == mu2->q)	       continue;  // opposite charge requirement
          if(mu2->pt        < PT_CUT)  continue;  // lepton pT cut
          if(fabs(mu2->eta) > ETA_CUT) continue;  // lepton |eta| cut
	  if(!passMuonID(mu2))	       continue;  // lepton selection

          LorentzVector vMu2(mu2->pt, mu2->eta, mu2->phi, MUON_MASS);  

          // trigger match
	  if(!(mu1->hltMatchBits & trigObj) && !(mu2->hltMatchBits & trigObj)) continue;
	  
	  // mass window
          LorentzVector vDilep = vMu1 + vMu2;
          if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
          
          
          /******** We have a Z candidate! HURRAY! ********/
          
          Double_t effdata, effmc;
	  Double_t corr=1;
	  
	  effdata=1; effmc=1;    
          if(mu1->q>0) { 
            effdata *= (1.-dataHLTEff_pos.getEff(mu1->eta, mu1->pt)); 
            effmc   *= (1.-zmmHLTEff_pos.getEff(mu1->eta, mu1->pt)); 
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(mu1->eta, mu1->pt)); 
            effmc   *= (1.-zmmHLTEff_neg.getEff(mu1->eta, mu1->pt)); 
          }
          if(mu2->q>0) {
            effdata *= (1.-dataHLTEff_pos.getEff(mu2->eta, mu2->pt)); 
            effmc   *= (1.-zmmHLTEff_pos.getEff(mu2->eta, mu2->pt));
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(mu2->eta, mu2->pt)); 
            effmc   *= (1.-zmmHLTEff_neg.getEff(mu2->eta, mu2->pt));
          }
          effdata = 1.-effdata;
          effmc   = 1.-effmc;
          corr *= effdata/effmc;
    
          effdata=1; effmc=1;
          if(mu1->q>0) { 
            effdata *= dataSelEff_pos.getEff(mu1->eta, mu1->pt); 
            effmc   *= zmmSelEff_pos.getEff(mu1->eta, mu1->pt); 
          } else {
            effdata *= dataSelEff_neg.getEff(mu1->eta, mu1->pt); 
            effmc   *= zmmSelEff_neg.getEff(mu1->eta, mu1->pt); 
          }
          if(mu2->q>0) {
            effdata *= dataSelEff_pos.getEff(mu2->eta, mu2->pt); 
            effmc   *= zmmSelEff_pos.getEff(mu2->eta, mu2->pt);
          } else {
            effdata *= dataSelEff_neg.getEff(mu2->eta, mu2->pt); 
            effmc   *= zmmSelEff_neg.getEff(mu2->eta, mu2->pt);
          }
          corr *= effdata/effmc;
    
          effdata=1; effmc=1;
          if(mu1->q>0) { 
            effdata *= dataStaEff_pos.getEff(fabs(mu1->eta), mu1->pt); 
            effmc   *= zmmStaEff_pos.getEff(fabs(mu1->eta), mu1->pt); 
          } else {
            effdata *= dataStaEff_neg.getEff(fabs(mu1->eta), mu1->pt); 
            effmc   *= zmmStaEff_neg.getEff(fabs(mu1->eta), mu1->pt); 
          }
          if(mu2->q>0) {
            effdata *= dataStaEff_pos.getEff(fabs(mu2->eta), mu2->pt); 
            effmc   *= zmmStaEff_pos.getEff(fabs(mu2->eta), mu2->pt);
          } else {
            effdata *= dataStaEff_neg.getEff(fabs(mu2->eta), mu2->pt); 
            effmc   *= zmmStaEff_neg.getEff(fabs(mu2->eta), mu2->pt);
          }
          corr *= effdata/effmc;
    
          effdata=1; effmc=1;
          if(mu1->q>0) { 
            effdata *= dataTrkEff_pos.getEff(fabs(mu1->eta), mu1->pt); 
            effmc   *= zmmTrkEff_pos.getEff(fabs(mu1->eta), mu1->pt); 
          } else {
            effdata *= dataTrkEff_neg.getEff(fabs(mu1->eta), mu1->pt); 
            effmc   *= zmmTrkEff_neg.getEff(fabs(mu1->eta), mu1->pt); 
          }
          if(mu2->q>0) {
            effdata *= dataTrkEff_pos.getEff(fabs(mu2->eta), mu2->pt); 
            effmc   *= zmmTrkEff_pos.getEff(fabs(mu2->eta), mu2->pt);
          } else {
            effdata *= dataTrkEff_neg.getEff(fabs(mu2->eta), mu2->pt); 
            effmc   *= zmmTrkEff_neg.getEff(fabs(mu2->eta), mu2->pt);
          }
          corr *= effdata/effmc;
    
	  nSelv[ifile]    +=weight;
	  nSelCorrv[ifile]+=weight*corr;
        }
      }      
    }
   
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);     accErrv.push_back(accv[ifile]*sqrt((1.-accv[ifile])/nEvtsv[ifile]));
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]); accCorrErrv.push_back(accCorrv[ifile]*sqrt((1.-accCorrv[ifile])/nEvtsv[ifile]));
    
    delete infile;
    infile=0, eventTree=0;  
  }  
  delete info;
  delete gen;
  delete muonArr;

    
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
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "   ================================================" << endl;
    cout << "    Label: " << labelv[ifile] << endl;
    cout << "     File: " << fnamev[ifile] << endl;
    cout << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "          nominal: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    cout << "     SF corrected: " << accCorrv[ifile] << " +/- " << accCorrErrv[ifile] << endl;
    cout << endl;
  }
  
  char txtfname[100];
  sprintf(txtfname,"%s/binned.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << " Z -> mu mu" << endl;
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
    txtfile << "     SF corrected: " << accCorrv[ifile] << " +/- " << accCorrErrv[ifile] << endl;
    txtfile << endl;
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelZmmBinned"); 
}
