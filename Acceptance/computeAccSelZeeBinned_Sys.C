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
#include "TGraph.h"
#include "TF1.h"

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections
#include "../EleScale/EnergyScaleCorrection.h" //EGMSmear

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions

// helper class to handle efficiency tables
#include "../Utils/CEffUser1D.hh"
#include "../Utils/CEffUser2D.hh"


#include "../Utils/AppEffSF.cc"
#endif

//=== MAIN MACRO ================================================================================================= 
void computeAccSelZeeBinned_Sys(const TString conf,            // input file
			    const TString inputDir,
                            const TString outputDir,        // output directory
          const TString  sqrts, 
			    const Int_t   doPU,
			    const Int_t   doScaleCorr,
			    const Int_t   sigma,
                const TString SysFileGSFSel="SysUnc_GSFSelEff.root",
                const bool is13TeV=1
) {
  gBenchmark->Start("computeAccSelZeeBinned");
  const int gainSeed = 12;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  // const Double_t ETA_CUT    = 2.5;
  const Double_t ETA_CUT    = 2.4; // to match muons
  // const Double_t ETA_CUT    = 1.444; // to match muons
  const Double_t ELE_MASS   = 0.000511;

  // const Double_t ETA_BARREL = 1.4442;
  // const Double_t ETA_ENDCAP = 1.566;
  const Double_t ETA_BARREL = 10.;
  const Double_t ETA_ENDCAP = 10.;

  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 11;


  const int muEtaNB = 12;
  const float muEtaRange[muEtaNB+1] = {-2.4,-2.0,-1.566,-1.4442,-1.0,-0.5,0,0.5,1.0,1.4442,1.566,2.0,2.4};
  const int muPtNB = 4;
  const float muPtRange[muPtNB+1] = {25,35,45,60,8000};

  AppEffSF effs(inputDir);
  effs.loadHLT("EleHLTEff_aMCxPythia","Positive","Negative");
  effs.loadSel("EleGSFSelEff_aMCxPythia","Positive","Negative");
  // effs.loadSta("MuStaEff_aMCxPythia","Combined","Combined");
  effs.loadUncSel(SysFileGSFSel);

  const TString corrFiles = "../EleScale/Run2017_LowPU_v2";
  EnergyScaleCorrection eleCorr( corrFiles.Data());// eleCorr.doScale= true; eleCorr.doSmearings =true;
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
  
  TH2D *h=0;

  TH2D *hHLTErr_pos = new TH2D("hHLTErr_pos", "",muEtaNB,muEtaRange,muPtNB,muPtRange);
  TH2D *hHLTErr_neg = new TH2D("hHLTErr_neg", "",muEtaNB,muEtaRange,muPtNB,muPtRange);
  
  TH2D *hGsfSelErr_pos = new TH2D("hGsfSelErr_pos", "",muEtaNB,muEtaRange,muPtNB,muPtRange);
  TH2D *hGsfSelErr_neg = new TH2D("hGsfSelErr_neg", "",muEtaNB,muEtaRange,muPtNB,muPtRange);

  
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
  vector<Double_t> nSelCorrvFSR, nSelCorrvMC, nSelCorrvBkg, nSelCorrvTag;
  vector<Double_t> nSelCorrVarvFSR, nSelCorrVarvMC, nSelCorrVarvBkg, nSelCorrVarvTag;
  vector<Double_t> accCorrvFSR, accCorrvMC, accCorrvBkg, accCorrvTag;
  vector<Double_t> accErrCorrvFSR, accErrCorrvMC, accErrCorrvBkg, accErrCorrvTag;

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
    nSelCorrv.push_back(0);      nSelCorrVarv.push_back(0);
    nSelCorrvFSR.push_back(0);   nSelCorrVarvFSR.push_back(0);
    nSelCorrvMC.push_back(0);    nSelCorrVarvMC.push_back(0);
    nSelCorrvBkg.push_back(0);   nSelCorrVarvBkg.push_back(0);
    nSelCorrvTag.push_back(0);   nSelCorrVarvTag.push_back(0);

    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    // for(UInt_t ientry=0.05*eventTree->GetEntries(); ientry<(uint)(0.1*eventTree->GetEntries()); ientry++) {
    // for(UInt_t ientry=0; ientry<(uint)(0.05*eventTree->GetEntries()); ientry++) {
    // for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry+=15) {
      if(ientry%100000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
      genBr->GetEntry(ientry);
      genPartArr->Clear(); genPartBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);

      Int_t glepq1=-99;
      Int_t glepq2=-99;

      if (fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      TLorentzVector *vec=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,&glepq1,&glepq2,1);
      if(vec->M()<MASS_LOW || vec->M()>MASS_HIGH) continue;      
      delete vec; delete lep1; delete lep2;

      vertexArr->Clear();
      vertexBr->GetEntry(ientry);
      double npv  = vertexArr->GetEntries();
      Double_t weight=gen->weight;
      if(doPU>0) weight*=h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));

      nEvtsv[ifile]+=weight;        

      // trigger requirement  
      if (!isEleTrigger(triggerMenu, info->triggerBits, kFALSE, is13TeV)) continue;      
      // if (!isEleTrigger(triggerMenu, info->triggerBits, kFALSE)) continue;
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
      
      electronArr->Clear();
      electronBr->GetEntry(ientry);

      for(Int_t i1=0; i1<electronArr->GetEntriesFast(); i1++) {
        const baconhep::TElectron *ele1 = (baconhep::TElectron*)((*electronArr)[i1]);

        TLorentzVector vEle1(0,0,0,0);
        vEle1.SetPtEtaPhiM(ele1->pt, ele1->eta, ele1->phi, ELE_MASS);
	
        // check ECAL gap
        if(fabs(vEle1.Eta())>=ETA_BARREL && fabs(vEle1.Eta())<=ETA_ENDCAP) continue;

        if(doScaleCorr && (ele1->r9 < 1.)){
          // set up variable and apply smear correction to ele1 
          float ele1Smear = 0.;
          float ele1Error = 0.;

          float ele1AbsEta   = fabs(vEle1.Eta());
          float ele1Et       = vEle1.E() / cosh(ele1AbsEta);
          
          double tagEcalE = ele1->ecalEnergy;
          double eTregress = tagEcalE/cosh(fabs(ele1->eta));
          bool  ele1isBarrel = ele1AbsEta < 1.4442;

          float ele1R9Prime = ele1->r9; // r9 corrections MC only
          double ele1Random = gRandom->Gaus(0,1);

        
          ele1Smear = eleCorr.smearingSigma(info->runNum, eTregress, ele1AbsEta, ele1R9Prime, gainSeed, 0., 0.);
          float ele1SmearEP = eleCorr.smearingSigma(info->runNum, eTregress, ele1AbsEta, ele1R9Prime, gainSeed, 1., 0.);
          float ele1SmearEM = eleCorr.smearingSigma(info->runNum, eTregress, ele1AbsEta, ele1R9Prime, gainSeed, -1., 0.);	

          if(sigma==0)  (vEle1) *= 1. + ele1Smear * ele1Random;
          else if(sigma==1) (vEle1) *= 1. + ele1SmearEP * ele1Random;
          else if(sigma==-1) (vEle1) *= 1. + ele1SmearEM * ele1Random;
          
        }

	
        if(vEle1.Pt()           < PT_CUT)     continue;  // lepton pT cut
        if(fabs(vEle1.Eta())    > ETA_CUT)    continue;  // lepton |eta| cut
        if(!passEleTightID(ele1, vEle1, info->rhoIso))     continue;  // lepton selection


        for(Int_t i2=i1+1; i2<electronArr->GetEntriesFast(); i2++) {         
          const baconhep::TElectron *ele2 = (baconhep::TElectron*)((*electronArr)[i2]);

          TLorentzVector vEle2(0,0,0,0);
          vEle2.SetPtEtaPhiM(ele2->pt, ele2->eta, ele2->phi, ELE_MASS); 

          if(fabs(vEle2.Eta())>=ETA_BARREL && fabs(vEle2.Eta())<=ETA_ENDCAP) continue;

          if(doScaleCorr && (ele2->r9 < 1.)){
            float ele2Smear = 0.;
            float ele2Error = 0.;

            float ele2AbsEta   = fabs(vEle2.Eta());
            float ele2Et       = vEle2.E() / cosh(ele2AbsEta);
            bool  ele2isBarrel = ele2AbsEta < 1.4442;

            float ele2R9Prime = ele2->r9; // r9 corrections MC only
              
            double ele2e = ele2->ecalEnergy;
            double eTregress2 = ele2e/cosh(fabs(ele2->eta));

	          double ele2Random = gRandom->Gaus(0,1);

          	 ele2Smear = eleCorr.smearingSigma(info->runNum, eTregress2, ele2AbsEta, ele2R9Prime, gainSeed, 0., 0.);
            float ele2SmearEP = eleCorr.smearingSigma(info->runNum, eTregress2, ele2AbsEta, ele2R9Prime, gainSeed, 1., 0.);
            float ele2SmearEM = eleCorr.smearingSigma(info->runNum, eTregress2, ele2AbsEta, ele2R9Prime, gainSeed, -1., 0.);	

            if(sigma==0) (vEle2) *= 1. + ele2Smear * ele2Random;
            else if(sigma==1) (vEle2) *= 1. + ele2SmearEP * ele2Random;
            else if(sigma==-1) (vEle2) *= 1. + ele2SmearEM * ele2Random;
          }
        
          if(ele1->q == ele2->q)	continue;
          if(vEle2.Pt()           < PT_CUT)     continue;  // lepton pT cut
          if(fabs(vEle2.Eta())    > ETA_CUT)    continue;  // lepton |eta| cut
          if(!passEleTightID(ele2, vEle2, info->rhoIso))     continue;  // lepton selection


          if(!isEleTriggerObj(triggerMenu, ele1->hltMatchBits, kFALSE, kFALSE, is13TeV)) continue;
          if(!isEleTriggerObj(triggerMenu, ele2->hltMatchBits, kFALSE, kFALSE, is13TeV)) continue;
          
          TLorentzVector vDilep = vEle1 + vEle2;
          if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
          
          /******** We have a Z candidate! HURRAY! ********/
          Double_t effdata, effmc;
          Double_t corr=1;
          Double_t effdataFSR, effdataMC, effdataBkg;
	  
          effdata=1; effmc=1;
          Double_t corrFSR=1;
          Double_t corrMC=1;
          Double_t corrBkg=1;
          Double_t corrTag=1;
          effdataFSR=1; effdataMC=1; effdataBkg=1;  
          int q1 = ele1->q;
          int q2 = ele2->q;

          
          corr = effs.fullEfficiencies(&vEle1,q1,&vEle2,q2);
          vector<double> uncs_gsf = effs.getUncSel(&vEle1,q1,&vEle2,q2);
          
          corrFSR *= uncs_gsf[0]*effs.computeHLTSF(&vEle1,q1,&vEle2,q2); // alternate fsr model
          corrMC  *= uncs_gsf[1]*effs.computeHLTSF(&vEle1,q1,&vEle2,q2); // alternate mc gen model
          corrBkg *= uncs_gsf[2]*effs.computeHLTSF(&vEle1,q1,&vEle2,q2); // alternate bkg model
          corrTag *= uncs_gsf[3]*effs.computeHLTSF(&vEle1,q1,&vEle2,q2); // alternate bkg model
          // corr *= effdata/effmc; // orig
	  
           
          double var=0.;        
          // var += effs.statUncSta(&l1, q1) + effs.statUncSta(&l2, q2);
          var += effs.statUncSel(&vEle1, q1, hGsfSelErr_pos, hGsfSelErr_neg, fabs(weight)*corr);
          var += effs.statUncSel(&vEle2, q2, hGsfSelErr_pos, hGsfSelErr_neg, fabs(weight)*corr);
          var += effs.statUncHLT(&vEle1, q1, hHLTErr_pos, hHLTErr_neg, fabs(weight)*corr);
          var += effs.statUncHLT(&vEle2, q2, hHLTErr_pos, hHLTErr_neg, fabs(weight)*corr);
          // cout << var1 << " " << var << endl;
          // std::cout << "event " << info->evtNum << " weight " << corr << std::endl;
          nSelv[ifile]    +=weight;
          nSelCorrvFSR[ifile]+=weight*corrFSR;
          nSelCorrvMC[ifile]+=weight*corrMC;
          nSelCorrvBkg[ifile]+=weight*corrBkg;
          nSelCorrvTag[ifile]+=weight*corrTag;
          nSelCorrv[ifile]+=weight*corr;
            // std::cout << "corr " << corr << " corr FSR " << corrFSR << "  corr MC " << corrMC << "  corr Bkg " << corrBkg << std::endl;
          nSelCorrVarvFSR[ifile]+=weight*weight*corrFSR*corrFSR;
          nSelCorrVarvMC[ifile]+=weight*weight*corrMC*corrMC;
          nSelCorrVarvBkg[ifile]+=weight*weight*corrBkg*corrBkg;
          nSelCorrVarvTag[ifile]+=weight*weight*corrTag*corrTag;
          nSelCorrVarv[ifile]+=weight*weight*corr*corr;
        }
      }
    }
    
    cout << "already there " << nSelCorrVarv[ifile]<< endl;
    
    Double_t var=0;
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
    std::cout << "blah  " << std::endl;
    nSelCorrVarvFSR[ifile]+=var;
    nSelCorrVarvMC[ifile]+=var;
    nSelCorrVarvBkg[ifile]+=var;
    nSelCorrVarv[ifile]+=var;
    cout << var << endl;
    std::cout << "comput acceptances" << std::endl;

    // compute acceptances
    std::cout << nEvtsv[ifile] << " " << nSelv[ifile] << std::endl;
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);   
    
    accCorrvFSR.push_back(nSelCorrvFSR[ifile]/nEvtsv[ifile]); 
    accCorrvMC.push_back(nSelCorrvMC[ifile]/nEvtsv[ifile]);   
    accCorrvBkg.push_back(nSelCorrvBkg[ifile]/nEvtsv[ifile]); 
    accCorrvTag.push_back(nSelCorrvTag[ifile]/nEvtsv[ifile]); 
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]); 
    
    
    accErrv.push_back(sqrt(accv[ifile]*(1.+accv[ifile])/nEvtsv[ifile]));
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
    cout << "     SF corrected: " << accCorrv[ifile]    << " +/- " << accErrCorrv[ifile]    << endl;
    cout << "          FSR unc: " << accCorrvFSR[ifile] << " +/- " << accErrCorrvFSR[ifile] << endl;
    cout << "           MC unc: " << accCorrvMC[ifile]  << " +/- " << accErrCorrvMC[ifile]  << endl;
    cout << "          Bkg unc: " << accCorrvBkg[ifile] << " +/- " << accErrCorrvBkg[ifile] << endl;
    cout << "          Tag unc: " << accCorrvTag[ifile] << " +/- " << accErrCorrvTag[ifile] << endl;
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
    txtfile << "     SF corrected: " << accCorrv[ifile]    << " +/- " << accErrCorrv[ifile]    << endl;
    txtfile << "          FSR unc: " << accCorrvFSR[ifile] << " +/- " << accErrCorrvFSR[ifile] << endl;
    txtfile << "           MC unc: " << accCorrvMC[ifile]  << " +/- " << accErrCorrvMC[ifile]  << endl;
    txtfile << "          Bkg unc: " << accCorrvBkg[ifile] << " +/- " << accErrCorrvBkg[ifile] << endl;
    txtfile << "          Bkg unc: " << accCorrvTag[ifile] << " +/- " << accErrCorrvTag[ifile] << endl;
    txtfile << endl;
  }
  txtfile.close();
  
 // char txtfname[100];
  sprintf(txtfname,"%s/sel_nums_only.txt",outputDir.Data());
  ofstream txtfile2;
  txtfile2.open(txtfname);
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
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
      
  gBenchmark->Show("computeAccSelZeeBinned"); 
}
