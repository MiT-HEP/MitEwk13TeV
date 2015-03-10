//================================================================================================
//
// Compute W->enu acceptance at full selection level
//
//  * outputs results summary text file
//
//  [!!!] propagation of efficiency scale factor uncertainties need to be checked
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TH1D.h>                   // histogram class
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

// define structures to read in ntuple
#include "../Ntupler/interface/EWKAnaDefs.hh"
#include "../Ntupler/interface/TEventInfo.hh"
#include "../Ntupler/interface/TGenInfo.hh"
#include "../Ntupler/interface/TElectron.hh"
#include "../Ntupler/interface/TVertex.hh"

// helper functions for lepton ID selection
#include "../Utils/LeptonIDCuts.hh"

// helper class to handle efficiency tables
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"
#endif


//=== MAIN MACRO ================================================================================================= 

void computeAccSelWe(const TString conf,       // input file
                     const TString outputDir,  // output directory
		     const Int_t   charge      // 0 = inclusive, +1 = W+, -1 = W-
) {
  gBenchmark->Start("computeAccSelWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.5;
  const Double_t ETA_BARREL = 1.4442;
  const Double_t ETA_ENDCAP = 1.566;
  
  // efficiency files
  TString dataHLTEffName("/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleHLTEff/analysis/eff.root");
  TString zeeHLTEffName("/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleHLTEff/analysis/eff.root");
  TString dataGsfSelEffName("/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleGsfSelEff/analysis/eff.root");
  TString zeeGsfSelEffName("/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleGsfSelEff/analysis/eff.root");
//  TString dataGsfEffName("../Efficiency/May23_EleGsfEff/analysis/eff.root");
//  TString zeeGsfEffName("../Efficiency/Zee_EleGsfEff/analysis/eff.root");
  if(charge==1) {
    dataHLTEffName = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleHLTEff_pos/analysis/eff.root";
    zeeHLTEffName  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleHLTEff_pos/analysis/eff.root";
    dataGsfSelEffName = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleGsfSelEff_pos/analysis/eff.root";
    zeeGsfSelEffName  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleGsfSelEff_pos/analysis/eff.root";
//    dataGsfEffName = "../Efficiency/May23_EleGsfEff_pos/analysis/eff.root";
//    zeeGsfEffName  = "../Efficiency/Zee_EleGsfEff_pos/analysis/eff.root";
  }
  if(charge==-1) {
    dataHLTEffName = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleHLTEff_neg/analysis/eff.root";
    zeeHLTEffName  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleHLTEff_neg/analysis/eff.root";
    dataGsfSelEffName = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/May23_EleGsfSelEff_neg/analysis/eff.root";
    zeeGsfSelEffName  = "/scratch/klawhorn/EWKAnaStore/8TeV/EfficiencyResults/Zee_EleGsfSelEff_neg/analysis/eff.root";
//    dataGsfEffName = "../Efficiency/May23_EleGsfEff_neg/analysis/eff.root";
//    zeeGsfEffName  = "../Efficiency/Zee_EleGsfEff_neg/analysis/eff.root";
  }

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
  // Get efficiency
  //
  TFile *dataHLTEffFile = new TFile(dataHLTEffName);
  CEffUser2D dataHLTEff;
  TH2D *hHLTErr=0, *hHLTErrB=0, *hHLTErrE=0;
  if(dataHLTEffName) {
    dataHLTEff.loadEff((TH2D*)dataHLTEffFile->Get("hEffEtaPt"),
                       (TH2D*)dataHLTEffFile->Get("hErrlEtaPt"),
                       (TH2D*)dataHLTEffFile->Get("hErrhEtaPt"));
    
    TH2D* h =(TH2D*)dataHLTEffFile->Get("hEffEtaPt");
    hHLTErr  = new TH2D("hHLTErr", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hHLTErrB = new TH2D("hHLTErrB","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hHLTErrE = new TH2D("hHLTErrE","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  }
  
  TFile *zeeHLTEffFile = new TFile(zeeHLTEffName);
  CEffUser2D zeeHLTEff;
  if(zeeHLTEffName) {
    zeeHLTEff.loadEff((TH2D*)zeeHLTEffFile->Get("hEffEtaPt"),
                      (TH2D*)zeeHLTEffFile->Get("hErrlEtaPt"),
                      (TH2D*)zeeHLTEffFile->Get("hErrhEtaPt"));
  }
  
  TFile *dataGsfSelEffFile = new TFile(dataGsfSelEffName);
  CEffUser2D dataGsfSelEff;
  TH2D *hGsfSelErr=0, *hGsfSelErrB=0, *hGsfSelErrE=0;
  if(dataGsfSelEffName) {
    dataGsfSelEff.loadEff((TH2D*)dataGsfSelEffFile->Get("hEffEtaPt"),
                       (TH2D*)dataGsfSelEffFile->Get("hErrlEtaPt"),
                       (TH2D*)dataGsfSelEffFile->Get("hErrhEtaPt"));
    
    TH2D* h =(TH2D*)dataGsfSelEffFile->Get("hEffEtaPt");
    hGsfSelErr  = new TH2D("hGsfSelErr", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hGsfSelErrB = new TH2D("hGsfSelErrB","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hGsfSelErrE = new TH2D("hGsfSelErrE","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  }
  
  TFile *zeeGsfSelEffFile = new TFile(zeeGsfSelEffName);
  CEffUser2D zeeGsfSelEff;
  if(zeeGsfSelEffName) {
    zeeGsfSelEff.loadEff((TH2D*)zeeGsfSelEffFile->Get("hEffEtaPt"),
                      (TH2D*)zeeGsfSelEffFile->Get("hErrlEtaPt"),
                      (TH2D*)zeeGsfSelEffFile->Get("hErrhEtaPt"));
  }
  
//  TFile *dataGsfEffFile = new TFile(dataGsfEffName);
//  CEffUser1D dataGsfEff;
//  TH1D *hGsfErr=0, *hGsfErrB=0, *hGsfErrE=0;
//  if(dataGsfEffName) {
//    dataGsfEff.loadEff((TGraphAsymmErrors*)dataGsfEffFile->Get("grEffEta"));
//    
//    TGraphAsymmErrors* gr =(TGraphAsymmErrors*)dataGsfEffFile->Get("grEffEta");
//    Double_t binning[gr->GetN()+1];
//    const Double_t *xval  = gr->GetX();
//    const Double_t *xerrl = gr->GetEXlow();
//    const Double_t *xerrh = gr->GetEXhigh();
//    binning[0] = xval[0]-xerrl[0];
//    for(Int_t i=0; i<gr->GetN(); i++) binning[i+1] = xval[i]+xerrh[i];
//    hGsfErr  = new TH1D("hGsfErr", "",gr->GetN(),binning);
//    hGsfErrB = new TH1D("hGsfErrB","",gr->GetN(),binning);
//    hGsfErrE = new TH1D("hGsfErrE","",gr->GetN(),binning);
//  }
//  
//  TFile *zeeGsfEffFile = new TFile(zeeGsfEffName);
//  CEffUser1D zeeGsfEff;
//  if(zeeGsfEffName) {
//    zeeGsfEff.loadEff((TGraphAsymmErrors*)zeeGsfEffFile->Get("grEffEta"));
//  }
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo   *gen   = new mithep::TGenInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  
  TFile *infile=0;
  TTree *eventTree=0;

  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv, nSelBv, nSelEv;
  vector<Double_t> accv, accBv, accEv;
  vector<Double_t> accErrv, accErrBv, accErrEv;
  vector<Double_t> nSelCorrv, nSelBCorrv, nSelECorrv;
  vector<Double_t> nSelCorrVarv, nSelBCorrVarv, nSelECorrVarv;
  vector<Double_t> accCorrv, accBCorrv, accECorrv;
  vector<Double_t> accErrCorrv, accErrBCorrv, accErrECorrv;
  
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
    nSelBv.push_back(0);
    nSelEv.push_back(0);
    nSelCorrv.push_back(0);
    nSelBCorrv.push_back(0);
    nSelECorrv.push_back(0);
    nSelCorrVarv.push_back(0);
    nSelBCorrVarv.push_back(0);
    nSelECorrVarv.push_back(0);
    
    for(Int_t iy=0; iy<=hHLTErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr->GetNbinsX(); ix++) {
        hHLTErr ->SetBinContent(ix,iy,0);
        hHLTErrB->SetBinContent(ix,iy,0);
        hHLTErrE->SetBinContent(ix,iy,0);
      }
    }
    for(Int_t iy=0; iy<=hGsfSelErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hGsfSelErr->GetNbinsX(); ix++) {
        hGsfSelErr ->SetBinContent(ix,iy,0);
        hGsfSelErrB->SetBinContent(ix,iy,0);
        hGsfSelErrE->SetBinContent(ix,iy,0);
      }
    }
//    for(Int_t ix=0; ix<=hGsfErr->GetNbinsX(); ix++) {
//      hGsfErr ->SetBinContent(ix,0);
//      hGsfErrB->SetBinContent(ix,0);
//      hGsfErrE->SetBinContent(ix,0);
//    }
    
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      genBr->GetEntry(ientry);
      if(charge==-1 && gen->id_1!= EGenType::kElectron) continue;  // check for W-
      if(charge== 1 && gen->id_2!=-EGenType::kElectron) continue;  // check for W+
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
      Int_t nLooseLep=0;
      const mithep::TElectron *goodEle=0;
      Bool_t passSel=kFALSE;
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
        
        // check ECAL gap
        if(fabs(ele->scEta)>=ETA_BARREL && fabs(ele->scEta)<=ETA_ENDCAP) continue;
        
        if(fabs(ele->scEta) > 2.5) continue;                  // loose lepton |eta| cut
        if(ele->scEt	    < 20)  continue;                  // loose lepton pT cut
        if(passEleLooseID(ele,info->rhoLowEta)) nLooseLep++;  // loose lepton selection
        if(nLooseLep>1) {  // extra lepton veto
          passSel=kFALSE;
          break;
        }
        
        if(fabs(ele->scEta) > ETA_CUT)      continue;  // lepton |eta| cut
        if(ele->scEt < PT_CUT)  	    continue;  // lepton pT cut
        if(!passEleID(ele,info->rhoLowEta)) continue;  // lepton selection
        if(!(ele->hltMatchBits & trigObj))  continue;  // check trigger matching
        
	if(charge!=0 && ele->q!=charge) continue;  // check charge (if necessary)
        
	passSel=kTRUE;
        goodEle = ele;  
      }
      
      if(passSel) {
        
	/******** We have a W candidate! HURRAY! ********/
        
	Bool_t isBarrel = (fabs(goodEle->scEta)<ETA_BARREL) ? kTRUE : kFALSE;
	
	// data/MC scale factor corrections
        Double_t corr=1;
        if(dataHLTEffFile && zeeHLTEffFile) {
	  Float_t sceta = goodEle->scEta;
	  if(fabs(goodEle->scEta)>=2.5) sceta *= 0.99;
          Double_t effdata = dataHLTEff.getEff(sceta, goodEle->scEt);
          Double_t effmc   = zeeHLTEff.getEff(sceta, goodEle->scEt);
          corr *= effdata/effmc;
        }
        if(dataGsfSelEffFile && zeeGsfSelEffFile) {
          Float_t sceta = TMath::Min(fabs(goodEle->scEta),Float_t(2.499));
	  Double_t effdata = dataGsfSelEff.getEff(sceta, goodEle->scEt);
          Double_t effmc   = zeeGsfSelEff.getEff(sceta, goodEle->scEt);
          corr *= effdata/effmc;
        }
//        if(dataGsfEffFile && zeeGsfEffFile) {
//          Float_t sceta = TMath::Min(fabs(goodEle->scEta),Float_t(2.499));
//	  Double_t effdata = dataGsfEff.getEff(sceta);
//          Double_t effmc   = zeeGsfEff.getEff(sceta);
//          corr *= effdata/effmc;
//        }
	
	// scale factor uncertainties
        if(dataHLTEffFile && zeeHLTEffFile) {
	  Float_t sceta = goodEle->scEta;
	  if(fabs(goodEle->scEta)>=2.5) sceta *= 0.99;
          Double_t effdata = dataHLTEff.getEff(sceta, goodEle->scEt);
          Double_t effmc   = zeeHLTEff.getEff(sceta, goodEle->scEt);
	  Double_t errdata = TMath::Max(dataHLTEff.getErrLow(sceta, goodEle->scEt),dataHLTEff.getErrHigh(sceta, goodEle->scEt));
	  Double_t errmc   = TMath::Max(zeeHLTEff.getErrLow(sceta, goodEle->scEt), zeeHLTEff.getErrHigh(sceta, goodEle->scEt));
	  Double_t err     = corr*sqrt(errdata*errdata/effdata/effdata+errmc*errmc/effmc/effmc);
	  hHLTErr->Fill(sceta,goodEle->scEt,err);
	  if(isBarrel) hHLTErrB->Fill(sceta,goodEle->scEt,err);
	  else         hHLTErrE->Fill(sceta,goodEle->scEt,err);
        }
        if(dataGsfSelEffFile && zeeGsfSelEffFile) {
          Float_t sceta = TMath::Min(fabs(goodEle->scEta),Float_t(2.499));
	  Double_t effdata = dataGsfSelEff.getEff(sceta, goodEle->scEt);
          Double_t effmc   = zeeGsfSelEff.getEff(sceta, goodEle->scEt);
          Double_t errdata = TMath::Max(dataGsfSelEff.getErrLow(sceta, goodEle->scEt),dataGsfSelEff.getErrHigh(sceta, goodEle->scEt));
	  Double_t errmc   = TMath::Max(zeeGsfSelEff.getErrLow(sceta, goodEle->scEt), zeeGsfSelEff.getErrHigh(sceta, goodEle->scEt));
	  Double_t err     = corr*sqrt(errdata*errdata/effdata/effdata+errmc*errmc/effmc/effmc);
	  hGsfSelErr->Fill(sceta,goodEle->scEt,err);
	  if(isBarrel) hGsfSelErrB->Fill(sceta,goodEle->scEt,err);
	  else         hGsfSelErrE->Fill(sceta,goodEle->scEt,err);
        }
//        if(dataGsfEffFile && zeeGsfEffFile) {
//          Float_t sceta = TMath::Min(fabs(goodEle->scEta),Float_t(2.499));
//	  Double_t effdata = dataGsfEff.getEff(sceta);
//          Double_t effmc   = zeeGsfEff.getEff(sceta);
//          Double_t errdata = TMath::Max(dataGsfEff.getErrLow(sceta),dataGsfEff.getErrHigh(sceta));
//	  Double_t errmc   = TMath::Max(zeeGsfEff.getErrLow(sceta), zeeGsfEff.getErrHigh(sceta));
//	  Double_t err     = corr*sqrt(errdata*errdata/effdata/effdata+errmc*errmc/effmc/effmc);  
//	  hGsfErr->Fill(sceta,err);
//	  if(isBarrel) hGsfErrB->Fill(sceta,err);
//	  else         hGsfErrE->Fill(sceta,err);
//        }	
	
	nSelv[ifile]+=weight;
	nSelCorrv[ifile]+=weight*corr;
	nSelCorrVarv[ifile]+=weight*weight*corr*corr;
  	if(isBarrel) { 
	  nSelBv[ifile]+=weight;
	  nSelBCorrv[ifile]+=weight*corr;
	  nSelBCorrVarv[ifile]+=weight*weight*corr*corr;
	  	
	} else { 
	  nSelEv[ifile]+=weight;
	  nSelECorrv[ifile]+=weight*corr;
	  nSelECorrVarv[ifile]+=weight*weight*corr*corr;
	}
      }
    }
    
    Double_t var=0, varB=0, varE=0;
    for(Int_t iy=0; iy<=hHLTErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr->GetNbinsX(); ix++) {
        Double_t err;
	err=hHLTErr->GetBinContent(ix,iy);  var+=err*err;
        err=hHLTErrB->GetBinContent(ix,iy); varB+=err*err;
        err=hHLTErrE->GetBinContent(ix,iy); varE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hGsfSelErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hGsfSelErr->GetNbinsX(); ix++) {
        Double_t err;
	err=hGsfSelErr->GetBinContent(ix,iy);  var+=err*err;
	err=hGsfSelErrB->GetBinContent(ix,iy); varB+=err*err;
	err=hGsfSelErrE->GetBinContent(ix,iy); varE+=err*err;
      }
    }
//    for(Int_t ix=0; ix<=hGsfErr->GetNbinsX(); ix++) {
//      Double_t err;
//      err=hGsfErr->GetBinContent(ix);  var+=err*err;
//      err=hGsfErrB->GetBinContent(ix); varB+=err*err;
//      err=hGsfErrE->GetBinContent(ix); varE+=err*err;
//    }
    nSelCorrVarv[ifile]+=var;
    nSelBCorrVarv[ifile]+=varB;
    nSelECorrVarv[ifile]+=varE;
    
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);   accErrv.push_back(sqrt(accv[ifile]*(1.-accv[ifile])/nEvtsv[ifile]));
    accBv.push_back(nSelBv[ifile]/nEvtsv[ifile]); accErrBv.push_back(sqrt(accBv[ifile]*(1.-accBv[ifile])/nEvtsv[ifile]));
    accEv.push_back(nSelEv[ifile]/nEvtsv[ifile]); accErrEv.push_back(sqrt(accEv[ifile]*(1.-accEv[ifile])/nEvtsv[ifile]));
    
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]);   accErrCorrv.push_back(accCorrv[ifile]*sqrt(nSelCorrVarv[ifile]/nSelCorrv[ifile]/nSelCorrv[ifile] + 1./nEvtsv[ifile]));
    accBCorrv.push_back(nSelBCorrv[ifile]/nEvtsv[ifile]); accErrBCorrv.push_back(accBCorrv[ifile]*sqrt(nSelBCorrVarv[ifile]/nSelBCorrv[ifile]/nSelBCorrv[ifile] + 1./nEvtsv[ifile]));
    accECorrv.push_back(nSelECorrv[ifile]/nEvtsv[ifile]); accErrECorrv.push_back(accECorrv[ifile]*sqrt(nSelECorrVarv[ifile]/nSelECorrv[ifile]/nSelECorrv[ifile] + 1./nEvtsv[ifile]));
    
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
  if(charge== 0) cout << " W -> e nu"  << endl;
  if(charge==-1) cout << " W- -> e nu" << endl;
  if(charge== 1) cout << " W+ -> e nu" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  cout << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  cout << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "   ================================================" << endl;
    cout << "    Label: " << labelv[ifile] << endl;
    cout << "     File: " << fnamev[ifile] << endl;
    cout << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "     barrel: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile];
    cout << "  ==eff corr==> " << accBCorrv[ifile] << " +/- " << accErrBCorrv[ifile] << endl;
    cout << "     endcap: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrEv[ifile];
    cout << "  ==eff corr==> " << accECorrv[ifile] << " +/- " << accErrECorrv[ifile] << endl;
    cout << "      total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile];
    cout << "  ==eff corr==> " << accCorrv[ifile]  << " +/- " << accErrCorrv[ifile] << endl;
    cout << endl;
  }
  
  char txtfname[100];
  sprintf(txtfname,"%s/sel.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  if(charge== 0) txtfile << " W -> e nu"  << endl;
  if(charge==-1) txtfile << " W- -> e nu" << endl;
  if(charge== 1) txtfile << " W+ -> e nu" << endl;
  txtfile << "  pT > " << PT_CUT << endl;
  txtfile << "  |eta| < " << ETA_CUT << endl;
  txtfile << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  txtfile << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  txtfile << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "   ================================================" << endl;
    txtfile << "    Label: " << labelv[ifile] << endl;
    txtfile << "     File: " << fnamev[ifile] << endl;
    txtfile << endl;
    txtfile << "    *** Acceptance ***" << endl;
    txtfile << "     barrel: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile];
    txtfile << "  ==eff corr==> " << accBCorrv[ifile] << " +/- " << accErrBCorrv[ifile] << endl;
    txtfile << "     endcap: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrBv[ifile];
    txtfile << "  ==eff corr==> " << accECorrv[ifile] << " +/- " << accErrECorrv[ifile] << endl;
    txtfile << "      total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile];
    txtfile << "  ==eff corr==> " << accCorrv[ifile]  << " +/- " << accErrCorrv[ifile] << endl;
    txtfile << endl;
  }
  txtfile.close();  
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelWe"); 
}
