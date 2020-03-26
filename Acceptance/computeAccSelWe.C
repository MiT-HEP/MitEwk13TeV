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
#include "TLorentzVector.h"
#include <TRandom3.h>
#include "TGraph.h"
#include "TF1.h"

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections
// #include "../EleScale/EnergyScaleCorrection_class.hh" //EGMSmear
#include "../EleScale/EnergyScaleCorrection.h" //EGMSmear

// helper class to handle efficiency tables
#include "../Utils/CEffUser1D.hh"
#include "../Utils/CEffUser2D.hh"
#include "../Utils/AppEffSF.cc"
#endif


//=== MAIN MACRO ================================================================================================= 

void computeAccSelWe(const TString conf,            // input file
			    const TString inputDir, // efficiency main directory
          const TString outputDir,        // output directory
		      const Int_t   charge,      // 0 = inclusive, +1 = W+, -1 = W-
			    const Int_t   doPU,
			    const Int_t   doScaleCorr,
			    const Int_t   sigma,
          const TString SysFileGSFSel="SysUnc_GSFSelEff.root",
          const bool is13TeV=1
) {
  gBenchmark->Start("computeAccSelWe");
  const int gainSeed = 12;
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.4;
  const Double_t ELE_MASS = 0.000511;

  const Double_t ETA_BARREL = 1.4442;
  const Double_t ETA_ENDCAP = 1.566;
  // const Double_t ETA_BARREL = 10.;
  // const Double_t ETA_ENDCAP = 10.;

  const Double_t VETO_PT   = 10;
  const Double_t VETO_ETA  = 2.4;

  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  const Int_t BOSON_ID  = 24;
  const Int_t LEPTON_ID = 11;
  
  const int muEtaNB = 12;
  const float muEtaRange[muEtaNB+1] = {-2.4,-2.0,-1.566,-1.4442,-1.0,-0.5,0,0.5,1.0,1.4442,1.566,2.0,2.4};
  const int muPtNB = 3;
  const float muPtRange[muPtNB+1] = {25,35,50,10000};
  // const int muPtNB = 11;
  // const float muPtRange[muPtNB+1] = {25,26,27,28,29,30,32,35,40,50,60,10000};
  
  const int NBptHLT = 12;
  const float ptrangeHLT[NBptHLT+1] = {25, 26.5, 28, 29.5, 31, 32.5, 35, 40, 45, 50, 60, 80, 10000};

  AppEffSF effs(inputDir);
  effs.loadHLT("EleHLTEff_aMCxPythia","Positive","Negative");
  effs.loadSel("EleGSFSelEff_aMCxPythia","Combined","Combined");
  // effs.loadSta("MuStaEff_aMCxPythia","Combined","Combined");
  effs.loadUncSel(SysFileGSFSel);

  const TString corrFiles = "../EleScale/Run2017_LowPU_v2";

  EnergyScaleCorrection eleCorr( corrFiles.Data(), EnergyScaleCorrection::ECALELF); 

  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/puWeights_76x.root", "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("puWeights");
 
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

  
  TH2D *hHLTErr_pos = new TH2D("hHLTErr_pos", "",muEtaNB,muEtaRange,NBptHLT,ptrangeHLT);
  TH2D *hHLTErr_neg = new TH2D("hHLTErr_neg", "",muEtaNB,muEtaRange,NBptHLT,ptrangeHLT);
  
  TH2D *hGsfSelErr_pos = new TH2D("hGsfSelErr_pos", "",muEtaNB,muEtaRange,muPtNB,muPtRange);
  TH2D *hGsfSelErr_neg = new TH2D("hGsfSelErr_neg", "",muEtaNB,muEtaRange,muPtNB,muPtRange);
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo    *info = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen  = new baconhep::TGenEventInfo();
  TClonesArray *electronArr = new TClonesArray("baconhep::TElectron");
  TClonesArray *genPartArr  = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  
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
  vector<Double_t> nSelCorrvFSR, nSelCorrvMC, nSelCorrvBkg, nSelCorrvTag;
  vector<Double_t> nSelCorrVarvFSR, nSelCorrVarvMC, nSelCorrVarvBkg, nSelCorrVarvTag;
  vector<Double_t> accCorrvFSR, accCorrvMC, accCorrvBkg, accCorrvTag;
  vector<Double_t> accErrCorrvFSR, accErrCorrvMC, accErrCorrvBkg, accErrCorrvTag;

  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");
 
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = TFile::Open(fnamev[ifile]); 
    assert(infile);
  
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info",             &info); TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("GenEvtInfo",        &gen); TBranch *genBr      = eventTree->GetBranch("GenEvtInfo");
    eventTree->SetBranchAddress("GenParticle",&genPartArr); TBranch *genPartBr  = eventTree->GetBranch("GenParticle");
    eventTree->SetBranchAddress("Electron",  &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");

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
    nSelCorrvFSR.push_back(0);
    nSelCorrVarvFSR.push_back(0);
    nSelCorrvMC.push_back(0);
    nSelCorrVarvMC.push_back(0);
    nSelCorrvBkg.push_back(0);
    nSelCorrVarvBkg.push_back(0);
    nSelCorrvTag.push_back(0);
    nSelCorrVarvTag.push_back(0);
    
    //
    // loop over events
    //
    // for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    for(UInt_t ientry=0; ientry<(uint)(0.1*eventTree->GetEntries()); ientry++) {
      if(ientry%100000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
      infoBr->GetEntry(ientry);
      genBr->GetEntry(ientry);      
      genPartArr->Clear(); genPartBr->GetEntry(ientry);
  
      if (charge==-1 && toolbox::flavor(genPartArr, -BOSON_ID)!=LEPTON_ID) continue;
      if (charge==1 && toolbox::flavor(genPartArr, BOSON_ID)!=-LEPTON_ID) continue;
      if (charge==0 && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      /*TLorentzVector *vec=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,1);*/
      Int_t glepq1=-99;
      Int_t glepq2=-99;
      TLorentzVector *gvec=new TLorentzVector(0,0,0,0);
      TLorentzVector *glep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *glep2=new TLorentzVector(0,0,0,0);
	    toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2,&glepq1,&glepq2,1);
   
      // TLorentzVector tvec=*glep1+*glep2;
      // TLorentzVector* genV=new TLorentzVector(0,0,0,0);
      // genV->SetPtEtaPhiM(tvec.Pt(), tvec.Eta(), tvec.Phi(), tvec.M());
      // genVPt   = tvec.Pt();
      // genVPhi  = tvec.Phi();
      // genVy    = tvec.Rapidity();
      // double genVMass = tvec.M();
      double mtgen  = sqrt( 2.0 * (glep1->Pt()) * (glep2->Pt()) * (1.0-cos(toolbox::deltaPhi(glep1->Phi(),glep2->Phi()))) );
      // if(mtgen < 40) continue;
      // if(mtgen < 40 || mtgen > 140) continue;
     // cout << "mass " << genVMass <<  "   mt " << mt << endl;
      
      vertexArr->Clear();
      vertexBr->GetEntry(ientry);
      double npv  = vertexArr->GetEntries();
      Double_t weight=gen->weight;
      if(doPU>0) weight*=h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));

      nEvtsv[ifile]+=weight;
      
      if (!isEleTrigger(triggerMenu, info->triggerBits, kFALSE, is13TeV)) continue;     
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
      
      electronArr->Clear();
      electronBr->GetEntry(ientry);
      Int_t nLooseLep=0;
      const baconhep::TElectron *goodEle=0;
      TLorentzVector vEle(0,0,0,0);
      TLorentzVector vElefinal(0,0,0,0);
      Bool_t passSel=kFALSE;
      int q = 0;
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
        const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[i]);
        vEle.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, ELE_MASS);
        // check ECAL gap
        // if(fabs(vEle.Eta())>=ECAL_GAP_LOW && fabs(vEle.Eta())<=ECAL_GAP_HIGH) continue;
        if(doScaleCorr && (ele->r9 < 1.)){
          float eleSmear = 0.;

          float eleAbsEta   = fabs(vEle.Eta());
          float eleEt       = vEle.E() / cosh(eleAbsEta);
          bool  eleisBarrel = eleAbsEta < 1.4442;

          float eleR9Prime = ele->r9; // r9 corrections MC only, none after 2016

          double eleRamdom = gRandom->Gaus(0,1);
          
          eleSmear = eleCorr.smearingSigma(info->runNum, eleEt, eleAbsEta, eleR9Prime, gainSeed, 0., 0.);
          float eleSmearEP = eleCorr.smearingSigma(info->runNum, eleEt, eleAbsEta, eleR9Prime, gainSeed, 1., 0.);
          float eleSmearEM = eleCorr.smearingSigma(info->runNum, eleEt, eleAbsEta, eleR9Prime, gainSeed, -1., 0.);

          if(sigma==0)(vEle) *= 1. + eleSmear * eleRamdom;
          else if(sigma == 1) (vEle) *= 1. + eleSmearEP * eleRamdom;
          else if(sigma == -1)(vEle) *= 1.  + eleSmearEM * eleRamdom;
        }

      
        if(fabs(vEle.Eta())    > VETO_ETA) continue;
        if(vEle.Pt()           < VETO_PT)  continue;
        if(passEleLooseID(ele, vEle, info->rhoIso)) nLooseLep++;

        if(nLooseLep>1) {  // extra lepton veto
          passSel=kFALSE;
          break;
        }

      if(vEle.Pt()           < PT_CUT)     continue;  // lepton pT cut
      if(fabs(vEle.Eta())    > ETA_CUT)    continue;  // lepton |eta| cut
      if(!passEleTightID(ele, vEle, info->rhoIso))     continue;  // lepton selection

      if(!isEleTriggerObj(triggerMenu, ele->hltMatchBits, kFALSE, kFALSE, is13TeV)) continue;
      q = ele->q;
      if(charge!=0 && ele->q!=charge) continue;  // check charge (if necessary)
    
	
      double mtreco  = sqrt( 2.0 * (ele->pt) * (info->pfMETC) * (1.0-cos(toolbox::deltaPhi(ele->phi,info->pfMETCphi))) );
      
      if(mtreco < 40) continue;
      // if(mtreco > 40 && mtreco < 140) continue;  
    
      passSel=kTRUE;
      goodEle = ele;  
      vElefinal = vEle;
    }
      
  if(passSel) {
        
    /******** We have a W candidate! HURRAY! ********/
    Bool_t isBarrel = (fabs(vElefinal.Eta())<ETA_BARREL) ? kTRUE : kFALSE;

    Double_t effdata, effmc;
    Double_t corr=1;
    Double_t effdataFSR, effdataMC, effdataBkg;

    effdata=1; effmc=1;
    Double_t corrFSR=1;
    Double_t corrMC=1;
    Double_t corrBkg=1;
    Double_t corrTag=1;
    effdataFSR=1; effdataMC=1; effdataBkg=1;  

          
    corr = effs.fullEfficiencies(&vElefinal,q);
    // corr = effs.dataOnly(&vElefinal,q);
    vector<double> uncs_gsf = effs.getUncSel(&vElefinal,q);
    
    corrFSR *= uncs_gsf[0]*effs.computeHLTSF(&vElefinal,q); // alternate fsr model
    corrMC  *= uncs_gsf[1]*effs.computeHLTSF(&vElefinal,q); // alternate mc gen model
    corrBkg *= uncs_gsf[2]*effs.computeHLTSF(&vElefinal,q); // alternate bkg model
    corrTag *= uncs_gsf[3]*effs.computeHLTSF(&vElefinal,q); // alternate bkg model

     
    double var=0.;        
    var += effs.statUncSel(&vElefinal, q, hGsfSelErr_pos, hGsfSelErr_neg, fabs(weight)*corr);
    var += effs.statUncHLT(&vElefinal, q, hHLTErr_pos, hHLTErr_neg, fabs(weight)*corr);
  
  
    nSelv[ifile]+=weight;
    nSelCorrv[ifile]+=weight*corr;
    // nSelCorrVarv[ifile]+=weight*weight*corr*corr;
    nSelCorrvFSR[ifile]+=weight*corrFSR;
	  nSelCorrvMC[ifile]+=weight*corrMC;
	  nSelCorrvBkg[ifile]+=weight*corrBkg;
	  nSelCorrvTag[ifile]+=weight*corrTag;
    nSelCorrVarvFSR[ifile]+=weight*weight*corrFSR*corrFSR;
	  nSelCorrVarvMC[ifile]+=weight*weight*corrMC*corrMC;
	  nSelCorrVarvBkg[ifile]+=weight*weight*corrBkg*corrBkg;
	  nSelCorrVarvTag[ifile]+=weight*weight*corrTag*corrTag;

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
    for(Int_t iy=1; iy<=hHLTErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=1; ix<=hHLTErr_pos->GetNbinsX(); ix++) {
        Double_t err=hHLTErr_pos->GetBinContent(ix,iy);
        var+=err*err;
        // cout << "hlt pos " << err*err << endl;
        err=hHLTErr_neg->GetBinContent(ix,iy);
        var+=err*err;
      }
    }
   cout << "var1: " << var << endl;
    for(Int_t iy=1; iy<=hGsfSelErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=1; ix<=hGsfSelErr_pos->GetNbinsX(); ix++) {
        Double_t err=hGsfSelErr_pos->GetBinContent(ix,iy);
        
        var+=err*err;
        err=hGsfSelErr_neg->GetBinContent(ix,iy);
        
        var+=err*err;
        // cout << "gsf pos " << err*err << endl;
      }
    }
    
    nSelCorrVarv[ifile]+=var;
    nSelBCorrVarv[ifile]+=varB;
    nSelECorrVarv[ifile]+=varE;
    nSelCorrVarvFSR[ifile]+=var;
    nSelCorrVarvMC[ifile]+=var;
    nSelCorrVarvBkg[ifile]+=var;
    nSelCorrVarvTag[ifile]+=var;
    
    
    
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);   accErrv.push_back(sqrt(accv[ifile]*(1.+accv[ifile])/nEvtsv[ifile]));
    accBv.push_back(nSelBv[ifile]/nEvtsv[ifile]); accErrBv.push_back(sqrt(accBv[ifile]*(1.+accBv[ifile])/nEvtsv[ifile]));
    accEv.push_back(nSelEv[ifile]/nEvtsv[ifile]); accErrEv.push_back(sqrt(accEv[ifile]*(1.+accEv[ifile])/nEvtsv[ifile]));
    
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]);   
    accErrCorrv.push_back(accCorrv[ifile]*sqrt(nSelCorrVarv[ifile]/nSelCorrv[ifile]/nSelCorrv[ifile] + 1./nEvtsv[ifile]));
    accBCorrv.push_back(nSelBCorrv[ifile]/nEvtsv[ifile]); 
    accErrBCorrv.push_back(accBCorrv[ifile]*sqrt(nSelBCorrVarv[ifile]/nSelBCorrv[ifile]/nSelBCorrv[ifile] + 1./nEvtsv[ifile]));
    accECorrv.push_back(nSelECorrv[ifile]/nEvtsv[ifile]); 
    accErrECorrv.push_back(accECorrv[ifile]*sqrt(nSelECorrVarv[ifile]/nSelECorrv[ifile]/nSelECorrv[ifile] + 1./nEvtsv[ifile]));
    
    accCorrvFSR.push_back(nSelCorrvFSR[ifile]/nEvtsv[ifile]); 
    accCorrvMC.push_back(nSelCorrvMC[ifile]/nEvtsv[ifile]); 
    accCorrvBkg.push_back(nSelCorrvBkg[ifile]/nEvtsv[ifile]); 
    accCorrvTag.push_back(nSelCorrvTag[ifile]/nEvtsv[ifile]); 
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]);  
    
    accErrCorrvFSR.push_back(accCorrvFSR[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrvFSR[ifile]*nSelCorrvFSR[ifile]) + 1./nEvtsv[ifile]));
    accErrCorrvMC.push_back(accCorrvMC[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrvMC[ifile]*nSelCorrvMC[ifile]) + 1./nEvtsv[ifile]));
    accErrCorrvBkg.push_back(accCorrvBkg[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrvBkg[ifile]*nSelCorrvBkg[ifile]) + 1./nEvtsv[ifile])); 
    accErrCorrvTag.push_back(accCorrvTag[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrvTag[ifile]*nSelCorrvTag[ifile]) + 1./nEvtsv[ifile])); 
    accErrCorrv.push_back(accCorrv[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrv[ifile]*nSelCorrv[ifile]) + 1./nEvtsv[ifile]));
    
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
    cout << "  ==total efficiency==> " <<  setw(4) << accCorrv[ifile]/accv[ifile] << endl;
    cout << "     SF corrected: " << accCorrv[ifile]    << " +/- " << accErrCorrv[ifile]    << endl;
    cout << "          pct: " << 100*accErrCorrv[ifile] /accCorrv[ifile] << endl;
    cout << "          FSR unc: " << accCorrvFSR[ifile] << " +/- " << accErrCorrvFSR[ifile] << endl;
    cout << "  ==pct diff (FSR) ==> " << 100*(accCorrvFSR[ifile]-accCorrv[ifile])/accCorrv[ifile] <<  endl;
    cout << "           MC unc: " << accCorrvMC[ifile]  << " +/- " << accErrCorrvMC[ifile]  << endl;
    cout << "  ==pct diff (MC) ==> " << 100*(accCorrvMC[ifile]-accCorrv[ifile])/accCorrv[ifile]  << endl;
    cout << "          Bkg unc: " << accCorrvBkg[ifile] << " +/- " << accErrCorrvBkg[ifile] << endl;
    cout << "  ==pct diff (bkg) ==> " << 100*(accCorrvBkg[ifile]-accCorrv[ifile])/accCorrv[ifile] << endl;
    cout << endl;
  }
  
  char txtfname[500];
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
    txtfile << "     endcap: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrEv[ifile];
    txtfile << "  ==eff corr==> " << accECorrv[ifile] << " +/- " << accErrECorrv[ifile] << endl;
    txtfile << "      total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile];
    txtfile << "  ==eff corr==> " << accCorrv[ifile]  << " +/- " << accErrCorrv[ifile] << endl;
    txtfile << "  ==total efficiency==> " <<  setw(4) << accCorrv[ifile]/accv[ifile] << endl;
    txtfile << "     SF corrected: " << accCorrv[ifile]    << " +/- " << accErrCorrv[ifile]    << endl;
    txtfile << "          FSR unc: " << accCorrvFSR[ifile] << " +/- " << accErrCorrvFSR[ifile] << endl;
    txtfile << "  ==pct diff (FSR) ==> " << 100*(accCorrvFSR[ifile]-accCorrv[ifile])/accCorrv[ifile]  << endl;
    txtfile << "           MC unc: " << accCorrvMC[ifile]  << " +/- " << accErrCorrvMC[ifile]  << endl;
    txtfile << "  ==pct diff (MC) ==> " << 100*(accCorrvMC[ifile]-accCorrv[ifile])/accCorrv[ifile]   <<endl;
    txtfile << "          Bkg unc: " << accCorrvBkg[ifile] << " +/- " << accErrCorrvBkg[ifile] << endl;
    txtfile << "  ==pct diff (bkg) ==> " << 100*(accCorrvBkg[ifile]-accCorrv[ifile])/accCorrv[ifile]  << endl;
    txtfile << endl;
  }
  txtfile.close();  
  
  // char txtfname[100];
  sprintf(txtfname,"%s/sel_nums_only.txt",outputDir.Data());
  ofstream txtfile2;
  txtfile2.open(txtfname);
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile2 << "uncorrected: " << accv[ifile]    << endl;
    txtfile2 << accCorrv[ifile]    << " " << accErrCorrv[ifile]    << endl;
    txtfile2 << accCorrvFSR[ifile] << endl;
    txtfile2 << accCorrvMC[ifile]  << endl;
    txtfile2 << accCorrvBkg[ifile] << endl;
    txtfile2 << accCorrvTag[ifile] << endl;
    // txtfile << accCorrvFSR[ifile]  << ", " << accCorrvMC[ifile] << ", " << accCorrvBkg[ifile] << ", " << accCorrvTag[ifile] << endl;

    txtfile2 << endl;
  }
  txtfile2.close();  
  
  // char txtfname[100];
  sprintf(txtfname,"%s/gsf_unc.txt",outputDir.Data());
  ofstream txtfile3;
  txtfile3.open(txtfname);
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile3 << accCorrv[ifile]      << endl;
    txtfile3 << accCorrvFSR[ifile] << endl;
    txtfile3 << accCorrvMC[ifile]  << endl;
    txtfile3 << accCorrvBkg[ifile] << endl;
    txtfile3 << accCorrvTag[ifile] << endl;
    // txtfile << accCorrvFSR[ifile]  << ", " << accCorrvMC[ifile] << ", " << accCorrvBkg[ifile] << ", " << accCorrvTag[ifile] << endl;

    txtfile3 << endl;
  }
  txtfile3.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelWe"); 
}
