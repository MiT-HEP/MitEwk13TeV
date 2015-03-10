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
#include "Math/LorentzVector.h"     // 4-vector class

// define structures to read in ntuple
#include "../Ntupler/interface/EWKAnaDefs.hh"
#include "../Ntupler/interface/TEventInfo.hh"
#include "../Ntupler/interface/TGenInfo.hh"
#include "../Ntupler/interface/TElectron.hh"

// helper functions for lepton ID selection
#include "../Utils/LeptonIDCuts.hh"

// helper class to handle efficiency tables
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


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
  
  // efficiency files
  const TString dataHLTEffName     = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleHLTEff/analysis/eff.root";
  const TString dataHLTEffName_pos = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleHLTEff_pos/analysis/eff.root";
  const TString dataHLTEffName_neg = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleHLTEff_neg/analysis/eff.root";
  const TString zeeHLTEffName      = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleHLTEff/analysis/eff.root";
  const TString zeeHLTEffName_pos  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleHLTEff_pos/analysis/eff.root";
  const TString zeeHLTEffName_neg  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleHLTEff_neg/analysis/eff.root";
  
  const TString dataGsfSelEffName     = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleGsfSelEff/analysis/eff.root";
  const TString dataGsfSelEffName_pos = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleGsfSelEff_pos/analysis/eff.root";
  const TString dataGsfSelEffName_neg = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleGsfSelEff_neg/analysis/eff.root";
  const TString zeeGsfSelEffName      = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleGsfSelEff/analysis/eff.root";
  const TString zeeGsfSelEffName_pos  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleGsfSelEff_pos/analysis/eff.root";
  const TString zeeGsfSelEffName_neg  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleGsfSelEff_neg/analysis/eff.root";
    

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
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo   *gen   = new mithep::TGenInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  
  TFile *infile=0;
  TTree *eventTree=0;
  
  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv, nSelCorrv;
  vector<Double_t> statErr2v, effErr2v;
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
    eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Gen",      &gen);         TBranch *genBr      = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");

    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelCorrv.push_back(0);
    statErr2v.push_back(0);
    
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
      ULong_t trigger = kHLT_Ele22_CaloIdL_CaloIsoVL;
      ULong_t trigObj = kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj;   
      if(!(info->triggerBits & trigger)) continue;
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
    
      electronArr->Clear();
      electronBr->GetEntry(ientry);
      for(Int_t i1=0; i1<electronArr->GetEntriesFast(); i1++) {
  	const mithep::TElectron *ele1 = (mithep::TElectron*)((*electronArr)[i1]);
	
	// check ECAL gap
	if(fabs(ele1->scEta)>=ETA_BARREL && fabs(ele1->scEta)<=ETA_ENDCAP) continue;
        
	if(ele1->scEt	     < PT_CUT)	     continue;  // lepton pT cut
        if(fabs(ele1->scEta) > ETA_CUT)	     continue;  // lepton |eta| cut
        if(!passEleID(ele1,info->rhoLowEta)) continue;  // lepton selection

        LorentzVector vEle1(ele1->pt, ele1->eta, ele1->phi, ELE_MASS);
	Bool_t isB1 = (fabs(ele1->scEta)<ETA_BARREL) ? kTRUE : kFALSE;

        for(Int_t i2=i1+1; i2<electronArr->GetEntriesFast(); i2++) {          
	  const mithep::TElectron *ele2 = (mithep::TElectron*)((*electronArr)[i2]);
	  
	  // check ECAL gap
	  if(fabs(ele2->scEta)>=ETA_BARREL && fabs(ele2->scEta)<=ETA_ENDCAP) continue;
        
          if(ele2->scEt        < PT_CUT)       continue;  // lepton pT cut
          if(fabs(ele2->scEta) > ETA_CUT)      continue;  // lepton |eta| cut
	  if(!passEleID(ele2,info->rhoLowEta)) continue;  // lepton selection

          LorentzVector vEle2(ele2->pt, ele2->eta, ele2->phi, ELE_MASS);  
          Bool_t isB2 = (fabs(ele2->scEta)<ETA_BARREL) ? kTRUE : kFALSE;

          // trigger match
	  if(!(ele1->hltMatchBits & trigObj) && !(ele2->hltMatchBits & trigObj)) continue;
	  
	  // mass window
          LorentzVector vDilep = vEle1 + vEle2;
          if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
          
          
          /******** We have a Z candidate! HURRAY! ********/
       
          Double_t effdata, effmc;
          Double_t sceta1 = (fabs(ele1->scEta)<2.5) ? ele1->scEta : 0.99*(ele1->scEta);
          Double_t sceta2 = (fabs(ele2->scEta)<2.5) ? ele2->scEta : 0.99*(ele2->scEta);
    
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
            effdata *= dataGsfSelEff_pos.getEff(fabs(sceta1), ele1->scEt); 
            effmc   *= zeeGsfSelEff_pos.getEff(fabs(sceta1), ele1->scEt); 
          } else {
            effdata *= dataGsfSelEff_neg.getEff(fabs(sceta1), ele1->scEt); 
            effmc   *= zeeGsfSelEff_neg.getEff(fabs(sceta1), ele1->scEt); 
          }
          if(ele2->q>0) {
            effdata *= dataGsfSelEff_pos.getEff(fabs(sceta2), ele2->scEt); 
            effmc   *= zeeGsfSelEff_pos.getEff(fabs(sceta2), ele2->scEt);
          } else {
            effdata *= dataGsfSelEff_neg.getEff(fabs(sceta2), ele2->scEt); 
            effmc   *= zeeGsfSelEff_neg.getEff(fabs(sceta2), ele2->scEt);
          }
          corr *= effdata/effmc;
	  
	  nSelv[ifile]+=weight;
	  nSelCorrv[ifile]+=weight*corr;
	  statErr2v[ifile]+=weight*weight*corr*corr;
	  
	  // scale factor uncertainties
//	  Double_t dataerr=0, mcerr=0;
//	  if(ele1->q>0) {	    
//	    effdata = dataGsfSelEff_pos.getEff(fabs(sceta1), ele1->scEt);
//	    errdata = TMath::Max(dataGsfSelEff_pos.getErrl(fabs(sceta1), ele1->scEt), dataGsfSelEff_pos.getErrh(fabs(sceta1), ele1->scEt));
//            effmc   = zeeGsfSelEff_pos.getEff(fabs(sceta1), ele1->scEt); 
//	    errmc   = TMath::Max(zeeGsfSelEff_pos.getErrl(fabs(sceta1), ele1->scEt), zeeGsfSelEff_pos.getErrh(fabs(sceta1), ele1->scEt));
//	    Double_t errGsfSel = weight*corr*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
//	    hGsfSelErr_pos->Fill(fabs(sceta1), ele1->scEt, errGsfSel);
//
//	  } else {
//	    effdata = dataGsfSelEff_neg.getEff(fabs(sceta1), ele1->scEt);
//	    errdata = TMath::Max(dataGsfSelEff_neg.getErrl(fabs(sceta1), ele1->scEt), dataGsfSelEff_neg.getErrh(fabs(sceta1), ele1->scEt));
//            effmc   = zeeGsfSelEff_neg.getEff(fabs(sceta1), ele1->scEt); 
//	    errmc   = TMath::Max(zeeGsfSelEff_neg.getErrl(fabs(sceta1), ele1->scEt), zeeGsfSelEff_neg.getErrh(fabs(sceta1), ele1->scEt));
//	    Double_t errGsfSel = weight*corr*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
//	    hGsfSelErr_neg->Fill(fabs(sceta1), ele1->scEt, errGsfSel);
//	  }
//	  if(ele2->q>0) {	    
//	    effdata = dataGsfSelEff_pos.getEff(fabs(sceta2), ele2->scEt);
//	    errdata = TMath::Max(dataGsfSelEff_pos.getErrl(fabs(sceta2), ele2->scEt), dataGsfSelEff_pos.getErrh(fabs(sceta2), ele2->scEt));
//            effmc   = zeeGsfSelEff_pos.getEff(fabs(sceta2), ele2->scEt); 
//	    errmc   = TMath::Max(zeeGsfSelEff_pos.getErrl(fabs(sceta2), ele2->scEt), zeeGsfSelEff_pos.getErrh(fabs(sceta2), ele2->scEt));
//	    Double_t errGsfSel = weight*corr*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
//	    hGsfSelErr_pos->Fill(fabs(sceta2), ele2->scEt, errGsfSel);
//
//	  } else {
//	    effdata = dataGsfSelEff_neg.getEff(fabs(sceta2), ele2->scEt);
//	    errdata = TMath::Max(dataGsfSelEff_neg.getErrl(fabs(sceta2), ele2->scEt), dataGsfSelEff_neg.getErrh(fabs(sceta2), ele2->scEt));
//            effmc   = zeeGsfSelEff_neg.getEff(fabs(sceta2), ele2->scEt); 
//	    errmc   = TMath::Max(zeeGsfSelEff_neg.getErrl(fabs(sceta2), ele2->scEt), zeeGsfSelEff_neg.getErrh(fabs(sceta2), ele2->scEt));
//	    Double_t errGsfSel = weight*corr*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
//	    hGsfSelErr_neg->Fill(fabs(sceta2), ele2->scEt, errGsfSel);
//	  }
        }
      }      
    }
    
//    Double_t err2=0;
//    for(Int_t iy=0; iy<=hGsfSelErr_pos->GetNbinsY(); iy++) {
//      for(Int_t ix=0; ix<=hGsfSelErr_pos->GetNbinsX(); ix++) {
//	err=hSelErr_pos->GetBinContent(ix,iy);
//	err2+=err*err;
//      }
//    }
//    for(Int_t iy=0; iy<=hGsfSelErr_neg->GetNbinsY(); iy++) {
//      for(Int_t ix=0; ix<=hGsfSelErr_neg->GetNbinsX(); ix++) {
//	err=hSelErr_neg->GetBinContent(ix,iy);
//	err2+=err*err;
//      }
//    }
//    effErr2v.push_back(err2);
    
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);         accErrv.push_back(sqrt(accv[ifile]*(1.-accv[ifile])/nEvtsv[ifile]));
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]); accCorrErrv.push_back(sqrt(accCorrv[ifile]*(1.-accCorrv[ifile])/nEvtsv[ifile]));
    
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
    txtfile << "     SF corrected: " << accCorrv[ifile] << " +/- " << accCorrErrv[ifile] << endl;
    txtfile << endl;
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelZeeBinned"); 
}
