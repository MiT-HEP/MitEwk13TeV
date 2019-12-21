//================================================================================================
// 
// Perform fit to extract Z->ee signal and efficiency simultaneously
//
//  * outputs plots and fit results summary
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>                         // histogram class
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include <TGaxis.h>
#include "TLorentzVector.h"           // 4-vector class

#include "../Utils/MyTools.hh"	          // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections

// helper class to handle efficiency tables
#include "../Utils/CEffUser1D.hh"
#include "../Utils/CEffUser2D.hh"

#include "../EleScale/EnergyScaleCorrection.h"

#include <../Utils/AppEffSF.cc>

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);
// void drawWMetPlots(string plotname, TH1D *diff, RooRealVar &x, RooDataHist* dat, RooAddPdf* pdf, RooHistPdf* ewk, RooAbsPdf* qcd, RooHistPdf* wsigp, string lumitext, TH1D* hData);
//=== MAIN MACRO ================================================================================================= 

void fitZee(const TString  inputDir,    // input directory
	     const TString  outputDir,   // output directory
          const TString  sqrts, 
             const Double_t lumi,        // integrated luminosity (/fb)
	     const Bool_t   normToData=0 // draw MC normalized to data
) {
  std::cout<<"---------------- STARTING ZEE ------------------------"<<std::endl;
  gBenchmark->Start("plotZee");
  gStyle->SetTitleOffset(1.100,"Y");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  const Double_t ELE_MASS  = 0.000511;
  //
  // input ntuple file names
  //
  enum { eData, eZee, eEWK, eTop, eDib, eWx, eZxx};  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;

  fnamev.push_back(inputDir + TString("/") + TString("data_select.root")); typev.push_back(eData);
  fnamev.push_back(inputDir + TString("/") + TString("zee_select.root"));   typev.push_back(eZee);
  
  
    // fnamev.push_back(inputDir + TString("/") + TString("wx_select.root"));  typev.push_back(eWx);
  // fnamev.push_back(inputDir + TString("/") + TString("zxx_select.root"));  typev.push_back(eZxx);
  // fnamev.push_back(inputDir + TString("/") + TString("dib_select.root"));  typev.push_back(eDib);
  // fnamev.push_back(inputDir + TString("/") + TString("top_select.root"));  typev.push_back(eTop);
  
  fnamev.push_back(inputDir + TString("/") + TString("wx0_select.root"));  typev.push_back(eWx);
  fnamev.push_back(inputDir + TString("/") + TString("wx1_select.root"));  typev.push_back(eWx);
  fnamev.push_back(inputDir + TString("/") + TString("wx2_select.root"));  typev.push_back(eWx);
  fnamev.push_back(inputDir + TString("/") + TString("zxx_select.root"));  typev.push_back(eZxx);
  fnamev.push_back(inputDir + TString("/") + TString("ww_select.root"));  typev.push_back(eDib);
  fnamev.push_back(inputDir + TString("/") + TString("wz_select.root"));  typev.push_back(eDib);
  fnamev.push_back(inputDir + TString("/") + TString("zz_select.root"));  typev.push_back(eDib);
  fnamev.push_back(inputDir + TString("/") + TString("top1_select.root"));  typev.push_back(eTop);
  fnamev.push_back(inputDir + TString("/") + TString("top2_select.root"));  typev.push_back(eTop);
  fnamev.push_back(inputDir + TString("/") + TString("top3_select.root"));  typev.push_back(eTop);
 
  //
  // Fit options
  //
  // const Int_t    NBINS     = 120;
  const Int_t    NBINS    = 60;
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  const Double_t PT_CUT    = 25;
  // const Double_t PT_CUT    = 30;
  // const Double_t ETA_CUT   = 1.444;
  const Double_t ETA_CUT   = 2.4;//4;
    
  // const Double_t ECAL_GAP_LOW  = 1.4442;
  // const Double_t ECAL_GAP_HIGH = 1.566;
  const Double_t ECAL_GAP_LOW  = 10.;
  const Double_t ECAL_GAP_HIGH = 10.;
  
  // efficiency files

  TString baseDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_"+sqrts+"/results/Zee/";
  // TString baseDir = "/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV_v5_EleMedID2017/results/Zee/";
  AppEffSF effs(baseDir);
  // effs.loadHLT("EleHLTEff_aMCxPythia","Combined","Combined");
  effs.loadHLT("EleHLTEff_aMCxPythia","Positive","Negative");
  effs.loadSel("EleGSFSelEff_aMCxPythia","Combined","Combined");
  // effs.loadSel("EleGSFSelEff_aMCxPythia","Positive","Negative");
  string sysDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics/";
  string SysFileGSFSel = sysDir + "SysUnc_EleGSFSelEff.root";
  effs.loadUncSel(SysFileGSFSel);
  
  
  TH2D *hErr  = new TH2D("hErr", "",10,0,10,20,0,20);

  // //
  // // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zee_DataBkg.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);


  
  // plot output file format
  const TString format("all");


  Int_t yield = 0;
  Double_t yield_zee_up = 0, yield_zee_dn = 0;
  Double_t yield_zee_noPrefire = 0;
  Double_t yield_zee_pfJet=0, yield_zee_pfPhoton=0;
  Double_t yield_wm = 0;
  Double_t yield_zee = 0, yield_zee_unc=0;
  Double_t yield_ewk = 0, yield_ewk_unc=0;
  Double_t yield_top = 0, yield_top_unc=0;
 

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  // event category enumeration
  enum { eEleEle2HLT=1, eEleEle1HLT1L1, eEleEle1HLT, eEleEleNoSel, eEleSC };
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  

  // histograms for full selection
  double ZPtBins[]={0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,400,1000};
  double PhiStarBins[]={0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.01,0.012,0.014,0.016,0.018,0.021,0.024,0.027,0.030,0.034,0.038,0.044,0.050,0.058,0.066,0.076,0.088,0.10,0.12,0.14,0.16,0.18,0.20,0.24,0.28,0.34,0.42,0.52,0.64,0.8,1.0,1.5,2,3};
  double Lep1PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
  double Lep2PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150};
  double LepNegPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
  double LepPosPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};


  TH2D *hCompareElePtEcalE = new TH2D("hCompare","",200,0,100,200,0,100);
   
  enum{mcUp,mcDown,fsrUp,fsrDown,bkgUp,bkgDown,tagptUp,tagptDown,effsUp,effsDown,lepsfUp,lepsfDown,pfireUp,pfireDown};
  const string vWeight[]={"mcUp","mcDown","fsrUp","fsrDown","bkgUp","bkgDown","tagptUp","tagptDown","effstatUp","effstatDown","lepsfUp","lepsfDown","prefireUp","prefireDown"};
  int nWeight = sizeof(vWeight)/sizeof(vWeight[0]);

  TH1D *hData = new TH1D("hData","",NBINS,MASS_LOW,MASS_HIGH); hData->Sumw2();
  TH1D *hZee  = new TH1D("hZee", "",NBINS,MASS_LOW,MASS_HIGH); hZee->Sumw2();
  TH1D *hZxx  = new TH1D("hZxx", "",NBINS,MASS_LOW,MASS_HIGH); hZxx->Sumw2();
  TH1D *hWx   = new TH1D("hWx",  "",NBINS,MASS_LOW,MASS_HIGH); hWx->Sumw2();
  TH1D *hTtb  = new TH1D("hTtb", "",NBINS,MASS_LOW,MASS_HIGH); hTtb->Sumw2();
  TH1D *hDib  = new TH1D("hDib", "",NBINS,MASS_LOW,MASS_HIGH); hDib->Sumw2();
  TH1D *hEWK  = new TH1D("hEWK", "",NBINS,MASS_LOW,MASS_HIGH); hEWK->Sumw2();
  TH1D *hMC   = new TH1D("hMC",  "",NBINS,MASS_LOW,MASS_HIGH); hMC->Sumw2();
  
  // uncertainty shapes: 
  // vector of hists
  TH1D **hDataUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2();  
  TH1D **hZeeUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2();  
  TH1D **hZxxUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2(); 
  TH1D **hWxUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2(); 
  TH1D **hDibUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2(); 
  TH1D **hTtbUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2();   
  for(int j=0; j < nWeight; ++j){
    char hname[150];//char type[50];
    // sprintf(type,"%s",(vWeight[j]).c_str());
    sprintf(hname,"hData_%s",(vWeight[j]).c_str());   hDataUnc[j] = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
    cout << hDataUnc[j]->GetName() << endl;
    sprintf(hname,"hZee_%s",(vWeight[j]).c_str());    hZeeUnc[j]  = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
    sprintf(hname,"hZxx_%s",(vWeight[j]).c_str());    hZxxUnc[j]  = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
    sprintf(hname,"hWx_%s",(vWeight[j]).c_str());     hWxUnc[j]   = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
    sprintf(hname,"hTtb_%s",(vWeight[j]).c_str());    hTtbUnc[j]  = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
    sprintf(hname,"hDib_%s",(vWeight[j]).c_str());    hDibUnc[j]  = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
  }
  
  TH1D *hDataUp = new TH1D("hDataUp","",NBINS,MASS_LOW,MASS_HIGH); hDataUp->Sumw2();
  TH1D *hZeeUp  = new TH1D("hZeeUp", "",NBINS,MASS_LOW,MASS_HIGH); hZeeUp->Sumw2();
  
  TH1D *hDataDown = new TH1D("hDataDown","",NBINS,MASS_LOW,MASS_HIGH); hDataDown->Sumw2();
  TH1D *hZeeDown  = new TH1D("hZeeDown", "",NBINS,MASS_LOW,MASS_HIGH); hZeeDown->Sumw2();
  
  

  const int nBinsLep1Pt= sizeof(Lep1PtBins)/sizeof(double)-1;
  TH1D *hDataLep1Pt = new TH1D("hDataLep1Pt","",nBinsLep1Pt,Lep1PtBins); hDataLep1Pt->Sumw2();
  TH1D *hZeeLep1Pt  = new TH1D("hZeeLep1Pt", "",nBinsLep1Pt,Lep1PtBins); hZeeLep1Pt->Sumw2();
     
     
  TH1D *hGausRandNtuple = new TH1D("hGausRandNtuple","",100,-3,3); hGausRandNtuple->Sumw2();
  TH1D *hGausRandHere = new TH1D("hGausRandHere","",100,-3,3); hGausRandHere->Sumw2();
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t scale1fb, scale1fbUp, scale1fbDown, genVMass;
  Float_t prefireWeight, prefireUp=1, prefireDown=1;
  Float_t prefireJet, prefirePhoton;
  Float_t lep1error, lep2error;
  Int_t   q1, q2;
  TLorentzVector *lep1=0, *lep2=0;
  TLorentzVector *lep1_raw=0, *lep2_raw=0;
  TLorentzVector *dilep=0, *dilepSC = 0;
  TLorentzVector *sc1=0, *sc2=0;

 

  TFile *infile=0;
  TTree *intree=0;
  Float_t r91=0; 
  Float_t r92=0;
  Float_t random=0;

  Double_t nDib=0, nWx=0, nZxx=0;
  Double_t nDibUnc=0, nWxUnc=0, nZxxUnc=0;

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
  
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);	    assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);

    intree -> SetBranchStatus("*",0);
    intree -> SetBranchStatus("runNum",1);
    intree -> SetBranchStatus("lumiSec",1);
    intree -> SetBranchStatus("evtNum",1);
    intree -> SetBranchStatus("category",1);
    intree -> SetBranchStatus("npv",1);
    intree -> SetBranchStatus("npu",1);
    intree -> SetBranchStatus("prefirePhoton",1);
    intree -> SetBranchStatus("prefireJet",1);
    intree -> SetBranchStatus("prefireWeight",1);
    intree -> SetBranchStatus("prefireUp",1);
    intree -> SetBranchStatus("prefireDown",1);
    intree -> SetBranchStatus("scale1fb",1);
    intree -> SetBranchStatus("scale1fbUp",1);
    intree -> SetBranchStatus("scale1fbDown",1);
    intree -> SetBranchStatus("genVMass",1);
    intree -> SetBranchStatus("q1",1);
    intree -> SetBranchStatus("q2",1);
    intree -> SetBranchStatus("lep1",1);
    intree -> SetBranchStatus("lep2",1);
    intree -> SetBranchStatus("lep1_raw",1);
    intree -> SetBranchStatus("lep2_raw",1);
    intree -> SetBranchStatus("sc1",1);
    intree -> SetBranchStatus("sc2",1);
    intree -> SetBranchStatus("dilep",1);
    // intree -> SetBranchStatus("dilepSC",1);
    intree -> SetBranchStatus("r91",1);
    intree -> SetBranchStatus("r92",1);
    // intree -> SetBranchStatus("random",1);
    intree -> SetBranchStatus("lep1error",1);
    intree -> SetBranchStatus("lep2error",1);

    intree->SetBranchAddress("runNum",     &runNum);      // event run number
    intree->SetBranchAddress("lumiSec",    &lumiSec);     // event lumi section
    intree->SetBranchAddress("evtNum",     &evtNum);      // event number
    intree->SetBranchAddress("category",   &category);    // dilepton category
    intree->SetBranchAddress("npv",        &npv);	  // number of primary vertices
    intree->SetBranchAddress("npu",        &npu);	  // number of in-time PU events (MC)
    intree->SetBranchAddress("prefirePhoton", &prefirePhoton);    // prefire weight for 2017 conditions (MC)
    intree->SetBranchAddress("prefireJet", &prefireJet);    // prefire weight for 2017 conditions (MC)
    intree->SetBranchAddress("prefireWeight", &prefireWeight);    // prefire weight for 2017 conditions (MC)
    intree->SetBranchAddress("prefireUp",     &prefireUp);    // prefire weight for 2017 conditions (MC)
    intree->SetBranchAddress("prefireDown",   &prefireDown);    // prefire weight for 2017 conditions (MC)
    intree->SetBranchAddress("scale1fb",   &scale1fb);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp",   &scale1fbUp);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown",   &scale1fbDown);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("genVMass",   &genVMass);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("q1",         &q1);	  // charge of tag lepton
    intree->SetBranchAddress("q2",         &q2);	  // charge of probe lepton
    intree->SetBranchAddress("lep1",       &lep1);        // tag lepton 4-vector
    intree->SetBranchAddress("lep2",       &lep2);        // probe lepton 4-vector
    intree->SetBranchAddress("lep1_raw",       &lep1_raw);        // tag lepton 4-vector
    intree->SetBranchAddress("lep2_raw",       &lep2_raw);        // probe lepton 4-vector
    intree->SetBranchAddress("sc1",       &sc1);        // sc1 4-vector
    intree->SetBranchAddress("sc2",       &sc2);        // sc2 4-vector
    intree->SetBranchAddress("dilep",     &dilep);        // sc2 4-vector
    // intree->SetBranchAddress("dilepSC",   &dilepSC);      // sc2 4-vector
    intree->SetBranchAddress("r91",       &r91);        // sc2 4-vector
    intree->SetBranchAddress("r92",       &r92);        // sc2 4-vector
    // intree->SetBranchAddress("random",       &random);        // sc2 4-vector
    intree->SetBranchAddress("lep1error",     &lep1error);        // sc2 4-vector
    intree->SetBranchAddress("lep2error",    &lep2error);        // sc2 4-vector
  
    TH1D* hGenWeights;
    double totalNorm = 1.0;
    cout << "Hello " << endl;
    if(typev[ifile] != eData ){
      cout << "get gen weights" << endl;
      hGenWeights = (TH1D*)infile->Get("hGenWeights");
      totalNorm = hGenWeights->Integral();
      cout << totalNorm << endl;
    }
  
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    // for(UInt_t ientry=0; ientry<(uint)(0.1*intree->GetEntries()); ientry++) {
    // for(UInt_t ientry=0; ientry<1000; ientry++) {
      if(ientry%100000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " percent done with this file." << endl;
      intree->GetEntry(ientry);
      // std::cout << "r92 " << r92 << std::endl;
      if(fabs(lep1->Eta()) > ETA_CUT)   continue;      
      if(fabs(lep2->Eta()) > ETA_CUT)   continue;
      
      if(q1*q2>0) continue;
      if(lep1->Pt()        < PT_CUT)    continue;
      if(lep2->Pt()       < PT_CUT)    continue;
      
      // if(lep1->Pt()        > 30)    continue;
      // if(lep2->Pt()        > 30)    continue;
      
      hCompareElePtEcalE->Fill(lep1->Pt(),sc1->Pt());
      if(fabs(lep1->Eta())>=ECAL_GAP_LOW && fabs(lep1->Eta())<=ECAL_GAP_HIGH) continue;
      if(fabs(lep2->Eta())>=ECAL_GAP_LOW && fabs(lep2->Eta())<=ECAL_GAP_HIGH) continue;
      
      float mass = 0;
      float pt = 0;
      float rapidity = 0;
      float phiacop=0;
      float costhetastar=0;
      float phistar=0;
     
      Double_t weight=1;
      if(typev[ifile]!=eData) {
        weight *= scale1fb*prefireWeight*lumi/totalNorm;
      }  
      
      if(!(category==1) && !(category==2) && !(category==3)) continue;
      if(typev[ifile]==eData) {
    
        TLorentzVector el1;
        TLorentzVector el2;
        el1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),ELE_MASS);
        el2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),ELE_MASS);

        Double_t lp1 = el1.Pt();
        Double_t lp2 = el2.Pt();
        Double_t lq1 = q1;
        Double_t lq2 = q2;
        

        
        TLorentzVector l1, l2;
        if(lp1>lp2) {
          l1.SetPtEtaPhiM(lp1,lep1->Eta(),lep1->Phi(),ELE_MASS);
          l2.SetPtEtaPhiM(lp2,lep2->Eta(),lep2->Phi(),ELE_MASS);
	      } else {
          l1.SetPtEtaPhiM(lp2,lep2->Eta(),lep2->Phi(),ELE_MASS);
          l2.SetPtEtaPhiM(lp1,lep1->Eta(),lep1->Phi(),ELE_MASS);
        }

        TLorentzVector l1U, l2U;
        l1U.SetPtEtaPhiM(lp1*(1+lep1error),lep1->Eta(),lep1->Phi(),ELE_MASS);
        l2U.SetPtEtaPhiM(lp2*(1+lep2error),lep2->Eta(),lep2->Phi(),ELE_MASS);
        
        TLorentzVector l1D, l2D;
        l1D.SetPtEtaPhiM(lp1*(1-lep1error),lep1->Eta(),lep1->Phi(),ELE_MASS);
        l2D.SetPtEtaPhiM(lp2*(1-lep2error),lep2->Eta(),lep2->Phi(),ELE_MASS);
        
        double massU=(l1U+l2U).M();
        double massD=(l1D+l2D).M();
        
        mass=(l1+l2).M();
        // mass=dilepSC->M();
        pt =(l1+l2).Pt();
        rapidity = (l1+l2).Rapidity();

        phiacop=TMath::Pi()-fabs(l1.DeltaPhi(l2));
        if(lq1<0) costhetastar=tanh(float((l1.Rapidity()-l2.Rapidity())/2));
        else costhetastar=tanh(float((l2.Rapidity()-l1.Rapidity())/2));
        phistar=tan(phiacop/2)*sqrt(1-pow(costhetastar,2));
        
        if(mass        < MASS_LOW)  continue;
        if(mass        > MASS_HIGH) continue; 
        // if(l1.Pt()        < PT_CUT)    continue;
        // if(l2.Pt()        < PT_CUT)    continue;

    
        // hDataEG->Fill(dilepSC->M());

        hData->Fill(mass); 
        hDataUp->Fill(massD); 
        hDataDown->Fill(massU); 
        yield++;
	  
      } else {
        Double_t lp1 = lep1->Pt();
        Double_t lp2 = lep2->Pt();
        Double_t lq1 = q1;
        Double_t lq2 = q2;
        // double rand1;
        
        // // set the smearings here
        // for(int i = 0; i < 5; ++i){
          // double rand = gRandom->Gaus(0,1);
          // if(i==2) rand1=rand;
        // }
        // hGausRandHere->Fill(rand1);
        TLorentzVector l1, l2;
	      l1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),ELE_MASS);
	      l2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),ELE_MASS);
	  
        double mll=(l1+l2).M();
        Double_t effdata, effmc;
        Double_t corr=1;
        Double_t corrFSR=1, corrBkg=1, corrMC=1, corrTag=1;
        Double_t eff2Bindata, eff2Binmc;
        Double_t corr2Bin=1;
        Double_t corrUp=1;
        Double_t corrDown=1;
        Double_t effSigShapedata;
        Double_t corrSigShape=1;
        Double_t effBkgShapedata;
        Double_t corrBkgShape=1;
        
        if(mll       < MASS_LOW)  continue;
        if(mll       > MASS_HIGH) continue;
        // if(lp1        < PT_CUT)    continue;
        // if(lp2        < PT_CUT)    continue;
      
        // if(genVMass>MASS_LOW && genVMass<MASS_HIGH) continue;      
      
        corr = effs.fullEfficiencies(&l1,q1,&l2,q2);
        // cout << corr1 << " " << corr << endl;
        
        vector<double> uncs_gsf = effs.getUncSel(&l1,q1,&l2,q2);
        
        corrFSR *= uncs_gsf[0]*effs.computeHLTSF(&l1,q1,&l2,q2); // alternate fsr model
        corrMC  *= uncs_gsf[1]*effs.computeHLTSF(&l1,q1,&l2,q2); // alternate mc gen model
        corrBkg *= uncs_gsf[2]*effs.computeHLTSF(&l1,q1,&l2,q2); // alternate bkg model
        corrTag *= uncs_gsf[3]*effs.computeHLTSF(&l1,q1,&l2,q2); // alternate bkg model
        
        double var=0.;        
        // var += effs.statUncSta(&l1, q1) + effs.statUncSta(&l2, q2);
        var += effs.statUncSel(&l1, q1, hErr, hErr, fabs(weight)*corr);
        var += effs.statUncSel(&l2, q2, hErr, hErr, fabs(weight)*corr);
        var += effs.statUncHLTDilep(&l1, q1, &l2, q2);
        // cout << var1 << " " << var << endl;
        
        corrUp=corr+sqrt(var);
        corrDown=corr-sqrt(var);  
      
        // corr=1;
        mass = (l1+l2).M();
        
        TLorentzVector l1U, l2U;
        l1U.SetPtEtaPhiM(lp1*(1+lep1error),lep1->Eta(),lep1->Phi(),ELE_MASS);
        l2U.SetPtEtaPhiM(lp2*(1+lep2error),lep2->Eta(),lep2->Phi(),ELE_MASS);
        
        TLorentzVector l1D, l2D;
        l1D.SetPtEtaPhiM(lp1*(1-lep1error),lep1->Eta(),lep1->Phi(),ELE_MASS);
        l2D.SetPtEtaPhiM(lp2*(1-lep2error),lep2->Eta(),lep2->Phi(),ELE_MASS);
        
        double massU=(l1U+l2U).M();
        double massD=(l1D+l2D).M();

        if(typev[ifile]==eZee)  {
        yield_zee += weight*corr;
        yield_zee_up += weight*corrUp;
        yield_zee_dn += weight*corrDown;
        yield_zee_noPrefire += scale1fb*lumi*corr/totalNorm;
        yield_zee_pfPhoton += scale1fb*lumi*corr*prefirePhoton/totalNorm;
        yield_zee_pfJet += scale1fb*lumi*corr*prefireJet/totalNorm;
        yield_zee_unc += weight*weight*corr*corr;
        if(genVMass<MASS_LOW || genVMass>MASS_HIGH) yield_wm+= weight*corr;
        
        hZeeUnc[mcUp]->Fill(mass,weight*corrMC);
        hZeeUnc[mcDown]->Fill(mass,weight*(corr+(corr-corrMC)));
        hZeeUnc[fsrUp]->Fill(mass,weight*corrFSR);
        hZeeUnc[fsrDown]->Fill(mass,weight*(corr+(corr-corrFSR)));
        hZeeUnc[bkgUp]->Fill(mass,weight*corrBkg);
        hZeeUnc[bkgDown]->Fill(mass,weight*(corr+(corr-corrBkg)));
        hZeeUnc[tagptUp]->Fill(mass,weight*corrTag);
        hZeeUnc[tagptDown]->Fill(mass,weight*(corr+(corr-corrTag)));
        hZeeUnc[effsUp]->Fill(mass,weight*(corr+sqrt(var)));
        hZeeUnc[effsDown]->Fill(mass,weight*(corr-sqrt(var)));
        // roch up/down
        hZeeUnc[lepsfUp]->Fill((l1U+l2U).M(),weight*corr);
        hZeeUnc[lepsfDown]->Fill((l1D+l2D).M(),weight*corr);
        
        hZeeUnc[pfireUp]->Fill(mass,scale1fb*lumi*corr*prefireUp/totalNorm);
        hZeeUnc[pfireDown]->Fill(mass,scale1fb*lumi*corr*prefireDown/totalNorm);

        
        hZee->Fill(mass,weight*corr); 
        hZeeUp->Fill(massU,weight*corr); 
        hZeeDown->Fill(massD,weight*corr); 
        // hZeeEG->Fill(dilepSC->M(),weight*corr); 
        hMC->Fill(mass,weight*corr);
        hZeeLep1Pt->Fill(l1.Pt(),weight*corr); 
      }
      if(typev[ifile]==eZxx){
        nZxx+=weight*corr;
        nZxxUnc+=weight*weight*corr*corr;
        hZxx->Fill(mass,weight*corr); 
        //mcUp,mcDown,fsrUp,fsrDown,bkgUp,bkgDown,tagptUp,tagptDown,effsUp,effsDown,lepsfUp,lepsfDown,pfireUp,pfireDown
        hZxxUnc[mcUp]->Fill(mass,weight*corrMC);
        hZxxUnc[mcDown]->Fill(mass,weight*(corr+(corr-corrMC)));
        hZxxUnc[fsrUp]->Fill(mass,weight*corrFSR);
        hZxxUnc[fsrDown]->Fill(mass,weight*(corr+(corr-corrFSR)));
        hZxxUnc[bkgUp]->Fill(mass,weight*corrBkg);
        hZxxUnc[bkgDown]->Fill(mass,weight*(corr+(corr-corrBkg)));
        hZxxUnc[tagptUp]->Fill(mass,weight*corrTag);
        hZxxUnc[tagptDown]->Fill(mass,weight*(corr+(corr-corrTag)));
        hZxxUnc[effsUp]->Fill(mass,weight*(corr+sqrt(var)));
        hZxxUnc[effsUp]->Fill(mass,weight*(corr-sqrt(var)));
        
        // roch up/down
        hZxxUnc[lepsfUp]->Fill((l1U+l2U).M(),weight*corr);
        hZxxUnc[lepsfDown]->Fill((l1D+l2D).M(),weight*corr);
        
        hZxxUnc[pfireUp]->Fill(mass,scale1fb*lumi*corr*prefireUp/totalNorm);
        hZxxUnc[pfireDown]->Fill(mass,scale1fb*lumi*corr*prefireDown/totalNorm);
    } if(typev[ifile]==eWx){
      
        nWx+=weight*corr;
        nWxUnc+=weight*weight*corr*corr;
        hWx->Fill(mass,weight*corr); 
      
        hWxUnc[mcUp]->Fill(mass,weight*corrMC);
        hWxUnc[mcDown]->Fill(mass,weight*(corr+(corr-corrMC)));
        hWxUnc[fsrUp]->Fill(mass,weight*corrFSR);
        hWxUnc[fsrDown]->Fill(mass,weight*(corr+(corr-corrFSR)));
        hWxUnc[bkgUp]->Fill(mass,weight*corrBkg);
        hWxUnc[bkgDown]->Fill(mass,weight*(corr+(corr-corrBkg)));
        hWxUnc[tagptUp]->Fill(mass,weight*corrTag);
        hWxUnc[tagptDown]->Fill(mass,weight*(corr+(corr-corrTag)));
        hWxUnc[effsUp]->Fill(mass,weight*(corr+sqrt(var)));
        hWxUnc[effsDown]->Fill(mass,weight*(corr-sqrt(var)));
        
        // roch up/down
        hWxUnc[lepsfUp]->Fill((l1U+l2U).M(),weight*corr);
        hWxUnc[lepsfDown]->Fill((l1D+l2D).M(),weight*corr);
        

        hWxUnc[pfireUp]->Fill(mass,scale1fb*lumi*corr*prefireUp/totalNorm);
        hWxUnc[pfireDown]->Fill(mass,scale1fb*lumi*corr*prefireDown/totalNorm);
      
    } if(typev[ifile]==eDib){
        // cout << "blah " << endl;
        // cout << weight << " " << corr << " " << mass << endl;
        nDib+=weight*corr;
        nDibUnc+=weight*weight*corr*corr;
        hDib->Fill(mass,weight*corr); 
        
        // cout << "bleh << " << endl;
        hDibUnc[mcUp]->Fill(mass,weight*corrMC);
        hDibUnc[mcDown]->Fill(mass,weight*(corr+(corr-corrMC)));
        hDibUnc[fsrUp]->Fill(mass,weight*corrFSR);
        hDibUnc[fsrDown]->Fill(mass,weight*(corr+(corr-corrFSR)));
        hDibUnc[bkgUp]->Fill(mass,weight*corrBkg);
        hDibUnc[bkgDown]->Fill(mass,weight*(corr+(corr-corrBkg)));
        hDibUnc[tagptUp]->Fill(mass,weight*corrTag);
        hDibUnc[tagptDown]->Fill(mass,weight*(corr+(corr-corrTag)));
        hDibUnc[effsUp]->Fill(mass,weight*(corr+sqrt(var)));
        hDibUnc[effsDown]->Fill(mass,weight*(corr-sqrt(var)));
        
        // roch up/down
        hDibUnc[lepsfUp]->Fill((l1U+l2U).M(),weight*corr);
        hDibUnc[lepsfDown]->Fill((l1D+l2D).M(),weight*corr);
        
        hDibUnc[pfireUp]->Fill(mass,scale1fb*lumi*corr*prefireUp/totalNorm);
        hDibUnc[pfireDown]->Fill(mass,scale1fb*lumi*corr*prefireDown/totalNorm);

        // cout << "blah " << endl;
    }
      if(typev[ifile]==eEWK || typev[ifile]==eWx || typev[ifile]==eZxx || typev[ifile]==eDib)  {
	      yield_ewk += weight*corr;
	      yield_ewk_unc += weight*weight*corr*corr;
	      hEWK->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);
      }
      if(typev[ifile]==eTop)  {
	      yield_top += weight*corr;
	      yield_top_unc += weight*weight*corr*corr;
                
        hTtbUnc[mcUp]->Fill(mass,weight*corrMC);
        hTtbUnc[mcDown]->Fill(mass,weight*(corr+(corr-corrMC)));
        hTtbUnc[fsrUp]->Fill(mass,weight*corrFSR);
        hTtbUnc[fsrDown]->Fill(mass,weight*(corr+(corr-corrFSR)));
        hTtbUnc[bkgUp]->Fill(mass,weight*corrBkg);
        hTtbUnc[bkgDown]->Fill(mass,weight*(corr+(corr-corrBkg)));
        hTtbUnc[tagptUp]->Fill(mass,weight*corrTag);
        hTtbUnc[tagptDown]->Fill(mass,weight*(corr+(corr-corrTag)));
        hTtbUnc[effsUp]->Fill(mass,weight*(corr+sqrt(var)));
        hTtbUnc[effsDown]->Fill(mass,weight*(corr-sqrt(var)));
        
        // roch up/down
        hTtbUnc[lepsfUp]->Fill((l1U+l2U).M(),weight*corr);
        hTtbUnc[lepsfDown]->Fill((l1D+l2D).M(),weight*corr);
        

        hTtbUnc[pfireUp]->Fill(mass,scale1fb*lumi*corr*prefireUp/totalNorm);
        hTtbUnc[pfireDown]->Fill(mass,scale1fb*lumi*corr*prefireDown/totalNorm);


        
	      hTtb->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);
      }
    }
      }
    delete infile;
    infile=0, intree=0;
  }

  outFile->cd();
  
  TString histfname = outputDir + TString("/Zmumu_Histograms.root");
  TFile *histFile = new TFile(histfname,"RECREATE");
  histFile->cd();
  hData->Write();
  hZee->Write();
  hZxx->Write();
  hWx->Write();
  hTtb->Write();
  hDib->Write();
  
  for(int j = 0; j < nWeight; ++j){
    // hDataUnc[j]->Write();
    hZeeUnc[j]->Write();
    hZxxUnc[j]->Write();
    hWxUnc[j]->Write();
    hTtbUnc[j]->Write();
    hDibUnc[j]->Write();
  }
  
  histFile->Write();
  histFile->Close();
  
  
  TCanvas *ca = new TCanvas("ca","ca",800,800);

  cout << "Main values " << endl;
  cout << hData->GetMean()<< ", " << hZee->GetMean() << ", " <<  hData->GetMean()/hZee->GetMean() << endl;
  

  hDataLep1Pt->Write();
  


  outFile->Write();
  outFile->Close();
  

 // std::cout << "egamma ntuples size, MC: " << hZeeEG->Integral() << "  data: " << hDataEG->Integral() << std::endl;
 std::cout << "ours   ntuples size, MC: " << hZee->Integral()   << "  data: " << hData->Integral()   << std::endl;

  double MCscale=hData->Integral()/hMC->Integral();
  double GausScale = hGausRandHere->Integral()/hGausRandNtuple->Integral();
  hGausRandNtuple->Scale(GausScale);

  if(normToData)
    {

      cout<<"Normalized to data: "<<MCscale<<endl;
      
      hZee->Scale(MCscale);
      hZeeUp->Scale(MCscale);
      hZeeDown->Scale(MCscale);
      hMC->Scale(MCscale);
      hEWK->Scale(MCscale);
      hTtb->Scale(MCscale);

    }

  std::cout << hData->Integral() << std::endl;
  std::cout << hMC->Integral() << std::endl;
  std::cout << hZee->Integral() << std::endl;
  std::cout << hEWK->Integral() << std::endl;
  std::cout << hTtb->Integral() << std::endl;


  for(int j=0;j!=hDataLep1Pt->GetNbinsX();++j)
    {
      hDataLep1Pt->SetBinContent(j+1,hDataLep1Pt->GetBinContent(j+1)/hDataLep1Pt->GetBinWidth(j+1));
      hDataLep1Pt->SetBinError(j+1,hDataLep1Pt->GetBinError(j+1)/hDataLep1Pt->GetBinWidth(j+1));
    }
    
    
    TF1* g1 = new TF1("g1","gaus",0,1);
    TF1* g2 = new TF1("g2","gaus",0,1);
    hGausRandHere->Fit("gaus");
    hGausRandNtuple->Fit("gaus");
    std::cout << "here, mean,rms: " << hGausRandHere->GetMean() << " " << hGausRandHere->GetRMS() << std::endl;
    std::cout << "ntup, mean,rms: " << hGausRandNtuple->GetMean() << " " << hGausRandNtuple->GetRMS() << std::endl;
    std::cout << " params here:  mean " << g1->GetParameter(0) << " sigma " << g1->GetParameter(1) << std::endl;
    std::cout << " params ntup:  mean " << g2->GetParameter(0) << " sigma " << g2->GetParameter(1) << std::endl;

  TH1D *hGausDiff = makeDiffHist(hGausRandHere,hGausRandNtuple,"hGausDiff");
  hGausDiff->SetMarkerStyle(kFullCircle); 
  hGausDiff->SetMarkerSize(0.9);
  
  TH1D *hZeeDiff = makeDiffHist(hData,hMC,"hZeeDiff");
  hZeeDiff->SetMarkerStyle(kFullCircle); 
  hZeeDiff->SetMarkerSize(0.9);


  TH1D *massUnc = new TH1D("massUnc","massUnc",NBINS, MASS_LOW, MASS_HIGH);
  TH1D *massUnc2 = new TH1D("massUnc2","massUnc2",NBINS, MASS_LOW, MASS_HIGH);
  
	for(int i =1 ; i <= NBINS ; ++i){
		massUnc->SetBinContent(i,0);
		massUnc2->SetBinContent(i,0);
		massUnc2->SetBinError(i,0);
    double dataup = (hData->GetBinContent(i)-hDataUp->GetBinContent(i))/hData->GetBinContent(i);
    double datadown = (hData->GetBinContent(i)-hDataDown->GetBinContent(i))/hData->GetBinContent(i);
    double zeeup = (hZee->GetBinContent(i)-hZeeUp->GetBinContent(i))/hData->GetBinContent(i);
    // double zeedown = (hZee->GetBinContent(i)-hZeeDown->GetBinContent(i))/hData->GetBinContent(i);
		// massUnc->SetBinError(i,sqrt(dataup*dataup+datadown*datadown+zeeup*zeeup+zeedown*zeedown));
		massUnc->SetBinError(i,sqrt(dataup*dataup+zeeup*zeeup));
    cout << "Bin Error is " << massUnc->GetBinError(i) << endl;
	}
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  char ylabel[100];     // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000);  
  
  char normtext[100];
  sprintf(normtext,"MC normalized to data (#times %.2f)",MCscale);  

  string norm="";
  if(normToData)norm="_norm";
  
  // plot colors
  Int_t linecolorZ   = kOrange-3;
  Int_t fillcolorZ   = kOrange-2;
  Int_t linecolorEWK = kOrange+10;
  Int_t fillcolorEWK = kOrange+7;
  Int_t linecolorTop = kGreen+2;
  Int_t fillcolorTop = kGreen-5;
  Int_t ratioColor   = kGray+2;

  TCanvas *c = MakeCanvas("c","c",800,800);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.01);
  c->cd(1)->SetLeftMargin(0.15);  
  c->cd(1)->SetRightMargin(0.07);  
  c->cd(1)->SetTickx(1);
  c->cd(1)->SetTicky(1);  
  c->cd(2)->SetPad(0,0,1.0,0.3);
  c->cd(2)->SetTopMargin(0.05);
  c->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.15);
  c->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetTickx(1);
  c->cd(2)->SetTicky(1); 
  TGaxis::SetMaxDigits(3);
  
  
  c->cd(1);
  hCompareElePtEcalE->Draw("colz");
  c->SaveAs("lepPt.png");
  
  hGausRandHere->Draw();
  hGausRandNtuple->Draw("same");
  c->SaveAs("test.png");
  
  
  // char[150] plotname;
  // sprintf(nname,"isop%d",i); sprintf(plotname, "wep_fitmetp_bin%i",i);
  // drawWMetPlots(plotname, hMetpDiff, pfmet, dataMetp_[nname], pdfMetp_[i], pdfEWKp_[i], doTemplate?(RooAbsPdf*)pdfQCDp_[i]:(RooAbsPdf*)qcdp_[i]->model, pdfWep_[i], lumitext, hDataMetp2d[i]);
  
  std::cout << "peak dat " << hData->GetBinCenter(hData->GetMaximumBin()) << std::endl;
  std::cout << "peak zee " << hZee->GetBinCenter(hZee->GetMaximumBin()) << std::endl;
  

  sprintf(ylabel,"Events / %.1f GeV",hData->GetBinWidth(1));
  CPlot plotZee("zee"+norm,"","",ylabel);
  plotZee.AddHist1D(hData,"data","E");
  plotZee.AddToStack(hZee,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZee.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZee.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZee.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZee.SetYRange(0.01,1.2*(hData->GetMaximum() + sqrt(hData->GetMaximum())));
  plotZee.TransLegend(0.1,-0.05);
  plotZee.Draw(c,kFALSE,format,1);

  CPlot plotZeeDiff("zee"+norm,"","M(e^{+}e^{-}) [GeV]","#frac{Data-Pred}{Data}");
  plotZeeDiff.AddHist1D(massUnc,"E3",kGray,1,1);
  plotZeeDiff.AddHist1D(massUnc2,"E3",kBlack,1,1);
  plotZeeDiff.AddHist1D(hZeeDiff,"EX0",ratioColor);
  plotZeeDiff.SetYRange(-0.2,0.2);
  plotZeeDiff.AddLine(MASS_LOW, 0,MASS_HIGH, 0,kBlack,1);
  plotZeeDiff.AddLine(MASS_LOW, 0.1,MASS_HIGH, 0.1,kBlack,3);
  plotZeeDiff.AddLine(MASS_LOW,-0.1,MASS_HIGH,-0.1,kBlack,3);
  plotZeeDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZee2("zeelog"+norm,"","",ylabel);
  plotZee2.AddHist1D(hData,"data","E");
  plotZee2.AddToStack(hEWK,"EWK",fillcolorEWK,linecolorEWK);
  plotZee2.AddToStack(hTtb,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZee2.AddToStack(hZee,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZee2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZee2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZee2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZee2.SetLogy();
  plotZee2.SetYRange(1e-4*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plotZee2.TransLegend(0.1,-0.05);
  plotZee2.Draw(c,kTRUE,format,1);


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  


  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  cout << endl;

  cout << " The Zee event yield is " << yield << " +/-" << sqrt(yield) << "." << endl;
  cout << " The Zee expected event yield (incl. prefire) " << yield_zee << " +/-" << sqrt(yield_zee_unc) << "." << endl;
  cout << " The Zee expected event yield (no. prefire) " << yield_zee_noPrefire << " +/-" << sqrt(yield_zee_unc) << "." << endl;
  cout << " The Zee stat up event yield is " << yield_zee_up << " +/-" << sqrt(yield_zee_unc) << "." << endl;
  cout << " The Zee stat down event yield is " << yield_zee_dn << " +/-" << sqrt(yield_zee_unc) << "." << endl;
  cout << " The EWK event yield is " << yield_ewk << " +/-" << sqrt(yield_ewk_unc) << "." << endl;
  cout << "  -> Dib event yield is " << nDib      << " +/-" << sqrt(nDibUnc)       << "." << endl;
  cout << "  -> Zxx event yield is " << nZxx      << " +/-" << sqrt(nZxxUnc)       << "." << endl;
  cout << "  -> Wx  event yield is " << nWx       << " +/-" << sqrt(nWxUnc)        << "." << endl;
  cout << " The Top event yield is " << yield_top << " +/-" << sqrt(yield_top_unc) << "." << endl;
  cout << " The Zee Data Yield w/ ewk&top removed: " << yield - yield_ewk - yield_top << endl;
  cout << " Zee yield w/ all bkg removed: " << yield - yield_ewk - yield_top - yield_wm << endl;
  cout << " Prefire scale factor :" << yield_zee_noPrefire/yield_zee << endl;  
  cout << " Prefire Jets only : " << yield_zee_pfJet << "  scale fac: " << yield_zee_noPrefire/yield_zee_pfJet << endl;
  cout << " Prefire Photons only : " << yield_zee_pfPhoton << "  scale fac: " << yield_zee_noPrefire/yield_zee_pfPhoton << endl;
  

  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  
  ofstream txtfile;
  char txtfname[100];  
  sprintf(txtfname,"%s/zee_yields.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  

  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;  
  txtfile << endl;

  txtfile << " The Zee event yield is " << yield << " +/-" << sqrt(yield) << "." << endl;
  txtfile << " The Zee expected event yield (incl. prefire) " << yield_zee << " +/-" << sqrt(yield_zee_unc) << "." << endl;
  txtfile << " The Zee expected event yield (no. prefire) " << yield_zee_noPrefire << " +/-" << sqrt(yield_zee_unc) << "." << endl;
  txtfile << " The Zee stat up event yield is " << yield_zee_up << " +/-" << sqrt(yield_zee_unc) << "." << endl;
  txtfile << " The Zee stat down event yield is " << yield_zee_dn << " +/-" << sqrt(yield_zee_unc) << "." << endl;
  txtfile << " The EWK event yield is " << yield_ewk << " +/-" << sqrt(yield_ewk_unc) << "." << endl;
  txtfile << "  -> Dib event yield is " << nDib      << " +/-" << sqrt(nDibUnc)       << "." << endl;
  txtfile << "  -> Zxx event yield is " << nZxx      << " +/-" << sqrt(nZxxUnc)       << "." << endl;
  txtfile << "  -> Wx  event yield is " << nWx       << " +/-" << sqrt(nWxUnc)        << "." << endl;
  txtfile << " The Top event yield is " << yield_top << " +/-" << sqrt(yield_top_unc) << "." << endl;
  txtfile << " The Zee Data Yield w/ ewk&top removed: " << yield - yield_ewk - yield_top << endl;
  txtfile << " Zee yield w/ all bkg removed: " << yield - yield_ewk - yield_top - yield_wm << endl;
  txtfile << " Prefire scale factor :" << yield_zee_noPrefire/yield_zee << endl;
  txtfile << " Prefire Jets only : " << yield_zee_pfJet << "  scale fac: " << yield_zee_noPrefire/yield_zee_pfJet << endl;
  txtfile << " Prefire Photons only : " << yield_zee_pfPhoton << "  scale fac: " << yield_zee_noPrefire/yield_zee_pfPhoton << endl;
  txtfile << std::endl;
  txtfile.close();
  
  
  
  gBenchmark->Show("plotZee");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------

TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = (TH1D*)hData->Clone("hDiff");
  hDiff->SetName(name);
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    Double_t diff=0;
    Double_t err=0;
    if(hData->GetBinContent(ibin)!=0)
      {
	diff = diff0/hData->GetBinContent(ibin);
	err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((1.0/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
      }
    hDiff->SetBinContent(ibin,diff);
    hDiff->SetBinError(ibin,err);   
  }
  
  hDiff->GetYaxis()->SetTitleOffset(0.42);
  hDiff->GetYaxis()->SetTitleSize(0.13);
  hDiff->GetYaxis()->SetLabelSize(0.10);
  hDiff->GetYaxis()->SetNdivisions(104);
  hDiff->GetYaxis()->CenterTitle();
  hDiff->GetXaxis()->SetTitleOffset(1.2);
  hDiff->GetXaxis()->SetTitleSize(0.13);
  hDiff->GetXaxis()->SetLabelSize(0.12);
  hDiff->GetXaxis()->CenterTitle();
  
  return hDiff;
}

