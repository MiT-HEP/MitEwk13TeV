//================================================================================================
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
// #include "../Utils/CEffUser1D.hh"
// #include "../Utils/CEffUser2D.hh"

//helper class to handle rochester corrections
// // #include <rochcor2015r.h>
// #include <muresolution_run2r.h>
#include "../RochesterCorr/RoccoR.cc"
#include "../Utils/AppEffSF.cc"

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

void orderByLepPt(TLorentzVector &mu1, TLorentzVector &mu2, Int_t &q1, Int_t &q2, double mu_MASS);

//=== MAIN MACRO ================================================================================================= 

void fitZmm(const TString  inputDir,    // input directory
          const TString  outputDir,   // output directory
          const TString  sqrts, 
          const Double_t lumi,        // integrated luminosity (/fb)
          const Bool_t   normToData=0 //draw MC normalized to data
) {
  gBenchmark->Start("plotZmm");
  gStyle->SetTitleOffset(1.100,"Y");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  const Double_t mu_MASS  = 0.1057;
  //
  // input ntuple file names
  //
  enum { eData, eZmm, eEWK, eTop ,eDib, eZxx, eWx};  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;

  fnamev.push_back(inputDir + TString("/") + TString("data_select.root"));      typev.push_back(eData);
  fnamev.push_back(inputDir + TString("/") + TString("zmm_select.raw.root"));   typev.push_back(eZmm);
  fnamev.push_back(inputDir + TString("/") + TString("top_select.raw.root"));   typev.push_back(eTop);
  fnamev.push_back(inputDir + TString("/") + TString("wx_select.raw.root"));    typev.push_back(eWx);
  fnamev.push_back(inputDir + TString("/") + TString("zxx_select.raw.root"));   typev.push_back(eZxx);
  fnamev.push_back(inputDir + TString("/") + TString("dib_select.raw.root"));   typev.push_back(eDib);

  //
  // Fit options
  //
  const Int_t    NBINS     = 60;
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.4;

  const bool doRoch = true;
  // efficiency files
 
  // new eff SF helper code
  // constructor-> construct and intialize the main file path
  // add 
  TString baseDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_"+sqrts+"/results/Zmm/";
  AppEffSF effs(baseDir);
  effs.loadHLT("MuHLTEff_aMCxPythia","Positive","Negative");
  effs.loadSel("MuSITEff_aMCxPythia","Combined","Combined");
  effs.loadSta("MuStaEff_aMCxPythia","Combined","Combined");
  
  string sysDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics/";
  string sysFileSIT = sysDir + "SysUnc_MuSITEff.root";
  string sysFileSta = sysDir + "SysUnc_MuStaEff.root";
  effs.loadUncSel(sysFileSIT);
  effs.loadUncSta(sysFileSta);
  
  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_DataBkg.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);


  // plot output file format
  const TString format("all");

  Int_t yield = 0;  
  Double_t yield_zmm_up = 0, yield_zmm_dn = 0;
  Double_t yield_zmm_noPrefire = 0;
  Double_t yield_zmm_pfJet=0, yield_zmm_pfPhoton=0;
  Double_t yield_wm = 0;
  Double_t yield_zmm = 0, yield_zmm_unc=0;
  Double_t yield_ewk = 0, yield_ewk_unc=0;
  Double_t yield_top = 0, yield_top_unc=0;
   
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  // event category enumeration
  enum { eMuMu2HLT=1, eMuMu1HLT1mu1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk }; // event category enum

 
   // histograms for full selection
  double ZPtBins[]={0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,400,1000};
  double Lep1PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  enum{mcUp,mcDown,fsrUp,fsrDown,bkgUp,bkgDown,tagptUp,tagptDown,effsUp,effsDown,lepsfUp,lepsfDown,pfireUp,pfireDown};
  const string vWeight[]={"mcUp","mcDown","fsrUp","fsrDown","bkgUp","bkgDown","tagptUp","tagptDown","effsUp","effsDown","lepsfUp","lepsfDown","pfireUp","pfireDown"};
  int nWeight = sizeof(vWeight)/sizeof(vWeight[0]);

  TH1D *hData = new TH1D("hData","",NBINS,MASS_LOW,MASS_HIGH); hData->Sumw2();
  TH1D *hZmm  = new TH1D("hZmm", "",NBINS,MASS_LOW,MASS_HIGH); hZmm->Sumw2();
  TH1D *hZxx  = new TH1D("hZxx", "",NBINS,MASS_LOW,MASS_HIGH); hZxx->Sumw2();
  TH1D *hWx   = new TH1D("hWx",  "",NBINS,MASS_LOW,MASS_HIGH); hWx->Sumw2();
  TH1D *hTtb  = new TH1D("hTtb", "",NBINS,MASS_LOW,MASS_HIGH); hTtb->Sumw2();
  TH1D *hDib  = new TH1D("hDib", "",NBINS,MASS_LOW,MASS_HIGH); hDib->Sumw2();
  TH1D *hEWK  = new TH1D("hEWK", "",NBINS,MASS_LOW,MASS_HIGH); hEWK->Sumw2();
  TH1D *hMC   = new TH1D("hMC",  "",NBINS,MASS_LOW,MASS_HIGH); hMC->Sumw2();
  
  // uncertainty shapes: 
  // vector of hists
  TH1D **hDataUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2();  
  TH1D **hZmmUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2();  
  TH1D **hZxxUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2(); 
  TH1D **hWxUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2(); 
  TH1D **hDibUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2(); 
  TH1D **hTtbUnc  = new TH1D*[nWeight];// hAntiDataMetm->Sumw2();   
  for(int j=0; j < nWeight; ++j){
    char hname[150];//char type[50];
    // sprintf(type,"%s",(vWeight[j]).c_str());
    sprintf(hname,"hData_%s",(vWeight[j]).c_str());   hDataUnc[j] = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
    cout << hDataUnc[j]->GetName() << endl;
    sprintf(hname,"hZmm_%s",(vWeight[j]).c_str());    hZmmUnc[j]  = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
    sprintf(hname,"hZxx_%s",(vWeight[j]).c_str());    hZxxUnc[j]  = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
    sprintf(hname,"hWx_%s",(vWeight[j]).c_str());     hWxUnc[j]   = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
    sprintf(hname,"hTtb_%s",(vWeight[j]).c_str());    hTtbUnc[j]  = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
    sprintf(hname,"hDib_%s",(vWeight[j]).c_str());    hDibUnc[j]  = new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH);
  }

  TH1D *hDataNPV = new TH1D("hDataNPV","",50,0,50); hDataNPV->Sumw2();
  TH1D *hZmmNPV  = new TH1D("hZmmNPV", "",50,0,50); hZmmNPV->Sumw2();
  TH1D *hEWKNPV  = new TH1D("hEWKNPV", "",50,0,50); hEWKNPV->Sumw2();
  TH1D *hTtbNPV  = new TH1D("hTtbNPV", "",50,0,50); hTtbNPV->Sumw2();
  TH1D *hMCNPV   = new TH1D("hMCNPV",  "",50,0,50); hMCNPV->Sumw2();


  const int nBinsZPt= sizeof(ZPtBins)/sizeof(double)-1;
  TH1D *hDataZPt = new TH1D("hDataZPt","",nBinsZPt,ZPtBins); hDataZPt->Sumw2();
  TH1D *hZmmZPt  = new TH1D("hZmmZPt", "",nBinsZPt,ZPtBins); hZmmZPt->Sumw2();
  TH1D *hEWKZPt  = new TH1D("hEWKZPt", "",nBinsZPt,ZPtBins); hEWKZPt->Sumw2();
  TH1D *hTtbZPt  = new TH1D("hTtbZPt", "",nBinsZPt,ZPtBins); hTtbZPt->Sumw2();
  TH1D *hMCZPt   = new TH1D("hMCZPt",  "",nBinsZPt,ZPtBins); hMCZPt->Sumw2();

  TH1D *hEWKZPt_EffBin  = new TH1D("hEWKZPt_EffBin", "",nBinsZPt,ZPtBins); hEWKZPt_EffBin->Sumw2();
  TH1D *hTtbZPt_EffBin  = new TH1D("hTtbZPt_EffBin", "",nBinsZPt,ZPtBins); hTtbZPt_EffBin->Sumw2();
  TH1D *hEWKZPt_EffStatUp  = new TH1D("hEWKZPt_EffStatUp", "",nBinsZPt,ZPtBins); hEWKZPt_EffStatUp->Sumw2();
  TH1D *hTtbZPt_EffStatUp  = new TH1D("hTtbZPt_EffStatUp", "",nBinsZPt,ZPtBins); hTtbZPt_EffStatUp->Sumw2();
  TH1D *hEWKZPt_EffStatDown  = new TH1D("hEWKZPt_EffStatDown", "",nBinsZPt,ZPtBins); hEWKZPt_EffStatDown->Sumw2();
  TH1D *hTtbZPt_EffStatDown  = new TH1D("hTtbZPt_EffStatDown", "",nBinsZPt,ZPtBins); hTtbZPt_EffStatDown->Sumw2();
  TH1D *hEWKZPt_EffSigShape  = new TH1D("hEWKZPt_EffSigShape", "",nBinsZPt,ZPtBins); hEWKZPt_EffSigShape->Sumw2();
  TH1D *hTtbZPt_EffSigShape  = new TH1D("hTtbZPt_EffSigShape", "",nBinsZPt,ZPtBins); hTtbZPt_EffSigShape->Sumw2();
  TH1D *hEWKZPt_EffBkgShape  = new TH1D("hEWKZPt_EffBkgShape", "",nBinsZPt,ZPtBins); hEWKZPt_EffBkgShape->Sumw2();
  TH1D *hTtbZPt_EffBkgShape  = new TH1D("hTtbZPt_EffBkgShape", "",nBinsZPt,ZPtBins); hTtbZPt_EffBkgShape->Sumw2();

  const int nBinsLep1Pt= sizeof(Lep1PtBins)/sizeof(double)-1;
  TH1D *hDataLep1Pt = new TH1D("hDataLep1Pt","",nBinsLep1Pt,Lep1PtBins); hDataLep1Pt->Sumw2();
  TH1D *hZmmLep1Pt  = new TH1D("hZmmLep1Pt", "",nBinsLep1Pt,Lep1PtBins); hZmmLep1Pt->Sumw2();
  TH1D *hEWKLep1Pt  = new TH1D("hEWKLep1Pt", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt->Sumw2();
  TH1D *hTtbLep1Pt  = new TH1D("hTtbLep1Pt", "",nBinsLep1Pt,Lep1PtBins); hTtbLep1Pt->Sumw2();
  TH1D *hMCLep1Pt   = new TH1D("hMCLep1Pt",  "",nBinsLep1Pt,Lep1PtBins); hMCLep1Pt->Sumw2();

  TH1D *hEWKLep1Pt_EffBin  = new TH1D("hEWKLep1Pt_EffBin", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt_EffBin->Sumw2();
  TH1D *hTtbLep1Pt_EffBin  = new TH1D("hTtbLep1Pt_EffBin", "",nBinsLep1Pt,Lep1PtBins); hTtbLep1Pt_EffBin->Sumw2();
  TH1D *hEWKLep1Pt_EffStatUp  = new TH1D("hEWKLep1Pt_EffStatUp", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt_EffStatUp->Sumw2();
  TH1D *hTtbLep1Pt_EffStatUp  = new TH1D("hTtbLep1Pt_EffStatUp", "",nBinsLep1Pt,Lep1PtBins); hTtbLep1Pt_EffStatUp->Sumw2();
  TH1D *hEWKLep1Pt_EffStatDown  = new TH1D("hEWKLep1Pt_EffStatDown", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt_EffStatDown->Sumw2();
  TH1D *hTtbLep1Pt_EffStatDown  = new TH1D("hTtbLep1Pt_EffStatDown", "",nBinsLep1Pt,Lep1PtBins); hTtbLep1Pt_EffStatDown->Sumw2();
  TH1D *hEWKLep1Pt_EffSigShape  = new TH1D("hEWKLep1Pt_EffSigShape", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt_EffSigShape->Sumw2();
  TH1D *hTtbLep1Pt_EffSigShape  = new TH1D("hTtbLep1Pt_EffSigShape", "",nBinsLep1Pt,Lep1PtBins); hTtbLep1Pt_EffSigShape->Sumw2();
  TH1D *hEWKLep1Pt_EffBkgShape  = new TH1D("hEWKLep1Pt_EffBkgShape", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt_EffBkgShape->Sumw2();
  TH1D *hTtbLep1Pt_EffBkgShape  = new TH1D("hTtbLep1Pt_EffBkgShape", "",nBinsLep1Pt,Lep1PtBins); hTtbLep1Pt_EffBkgShape->Sumw2();
  
  
  TH2D *hErr  = new TH2D("hErr", "",10,0,10,20,0,20);
    
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t scale1fb, scale1fbUp, scale1fbDown, genVMass;
  Float_t prefireWeight, prefireUp, prefireDown;
  Float_t prefirePhoton, prefireJet;
  Int_t   q1, q2;
  TLorentzVector *lep1=0, *lep2=0;
  TLorentzVector *genlep1=0, *genlep2=0;
  Float_t genMuonPt1, genMuonPt2;

  Double_t nDib=0, nWx=0, nZxx=0;
  Double_t nDibUnc=0, nWxUnc=0, nZxxUnc=0;

  // Loading the Rochster Corrections
  RoccoR  rc("../RochesterCorr/RoccoR2017.txt");

  TFile *infile=0;
  TTree *intree=0;

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
    intree -> SetBranchStatus("prefireWeight",1);
    intree -> SetBranchStatus("prefirePhoton",1);
    intree -> SetBranchStatus("prefireJet",1);
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
    intree -> SetBranchStatus("genlep1",1);
    intree -> SetBranchStatus("genlep2",1);
    intree -> SetBranchStatus("genMuonPt1",1);
    intree -> SetBranchStatus("genMuonPt2",1);

    intree->SetBranchAddress("runNum",     &runNum);      // event run number
    intree->SetBranchAddress("lumiSec",    &lumiSec);     // event lumi section
    intree->SetBranchAddress("evtNum",     &evtNum);      // event number
    intree->SetBranchAddress("category",   &category);    // dilepton category
    intree->SetBranchAddress("npv",        &npv);	  // number of primary vertices
    intree->SetBranchAddress("npu",        &npu);	  // number of in-time PU events (MC)
    intree->SetBranchAddress("prefireWeight",   &prefireWeight);    // prefire weights for 2017 (MC)
    intree->SetBranchAddress("prefirePhoton",   &prefirePhoton);    // prefire weights for 2017 (MC)
    intree->SetBranchAddress("prefireJet",   &prefireJet);    // prefire weights for 2017 (MC)
    intree->SetBranchAddress("prefireUp",   &prefireUp);    // prefire weights for 2017 (MC)
    intree->SetBranchAddress("prefireDown",   &prefireDown);    // prefire weights for 2017 (MC)
    intree->SetBranchAddress("scale1fb",   &scale1fb);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp",   &scale1fbUp);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown",   &scale1fbDown);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("genVMass",   &genVMass);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("q1",         &q1);	  // charge of tag lepton
    intree->SetBranchAddress("q2",         &q2);	  // charge of probe lepton
    intree->SetBranchAddress("lep1",       &lep1);        // tag lepton 4-vector
    intree->SetBranchAddress("lep2",       &lep2);        // probe lepton 4-vector
    intree->SetBranchAddress("genlep1",       &genlep1);        // tag lepton 4-vector
    intree->SetBranchAddress("genlep2",       &genlep2);        // probe lepton 4-vector
    intree->SetBranchAddress("genMuonPt1",       &genMuonPt1);        // probe lepton 4-vector
    intree->SetBranchAddress("genMuonPt2",       &genMuonPt2);        // probe lepton 4-vector
    
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    // for(int ientry=0; ientry<0.1*(intree->GetEntries()); ientry++) {
    // for(int ientry=0; ientry<1000; ientry++) {
      if(ientry%100000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " percent done with this file." << endl;
      intree->GetEntry(ientry);
 
      // cout << "hello " << endl;
      if(fabs(lep1->Eta()) > ETA_CUT)   continue;      
      if(fabs(lep2->Eta()) > ETA_CUT)   continue;
      if(lep1->Pt() < PT_CUT) continue;
      if(lep2->Pt() < PT_CUT) continue;
      if(q1*q2>0) continue;
     
      // cout << "pass kinematics " << endl;
      float mass = 0;
      float pt = 0;
      // float rapidity = 0;
      // float phiacop=0;
      // float costhetastar=0;
      // float phistar=0;
     
      Double_t weight=1;      
      if(typev[ifile]!=eData) {
        weight *= scale1fb*prefireWeight*lumi;
      } 
      
      // fill Z events passing selection
        if(!(category==eMuMu2HLT) && !(category==eMuMu1HLT) && !(category==eMuMu1HLT1mu1)) continue;
        // cout << "pass trigger?" << endl;
        
        if(typev[ifile]==eData) {

          TLorentzVector mu1;
          TLorentzVector mu2;
          mu1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),mu_MASS);
          mu2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),mu_MASS);
          float qter1=1.0;
          float qter2=1.0;

          double dtSF1 = rc.kScaleDT(q1, mu1.Pt(), mu1.Eta(), mu1.Phi());//, s=0, m=0);
          double dtSF2 = rc.kScaleDT(q2, mu2.Pt(), mu2.Eta(), mu2.Phi());//s=0, m=0);
          if(doRoch){
            mu1*=dtSF1;
            mu2*=dtSF2;
          }
      

          Double_t lp1 = mu1.Pt();
          Double_t lp2 = mu2.Pt();
          Int_t q1 = q1;
          Int_t q2 = q2;

          // TLorentzVector mu1, mu2;
          orderByLepPt(mu1, mu2, q1, q2, mu_MASS);

          mass=(mu1+mu2).M();
          pt =(mu1+mu2).Pt();
          
          if(mass        < MASS_LOW)  continue;
          if(mass        > MASS_HIGH) continue;
          // if(mu1.Pt()        < PT_CUT)    continue;
          // if(mu2.Pt()        < PT_CUT)    continue;
            
          hData->Fill(mass); 
          hDataNPV->Fill(npv);
          hDataZPt->Fill(pt); 
          hDataLep1Pt->Fill(mu1.Pt()); 
          
          yield++;
	
	} else {

	  TLorentzVector mu1;
	  TLorentzVector mu2;
	  mu1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),mu_MASS);
	  mu2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),mu_MASS);
    TLorentzVector mu1u, mu1d, mu2u, mu2d;
    mu1u = mu1; mu1d = mu1;
    mu2u = mu2; mu2d = mu2;
    // std::cout << "MC pt " << mu1.Pt() << "  " << mu2.Pt() << std::endl;
     // cout << "gen leptons " << genMuonPt1 << " " << genMuonPt2 << endl;
     // cout << "other gen leptons " << genlep1->Pt() << " " << genlep2->Pt() << endl;
      
	  float qter1=1.0;
	  float qter2=1.0;
    double mcSF1 = rc.kSpreadMC(q1, mu1.Pt(), mu1.Eta(), mu1.Phi(), genMuonPt1);//, s=0, m=0);
    double mcSF2 = rc.kSpreadMC(q2, mu2.Pt(), mu2.Eta(), mu2.Phi(), genMuonPt2);//, s=0, m=0);
    if(typev[ifile]==eDib || genlep1->Pt()==0 || genlep2->Pt()==0) {
      mcSF1=1; 
      mcSF2=1;
    }// if gen stuff is messed up
    if(doRoch){
      mu1*=mcSF1;
      mu2*=mcSF2;
    }
    
    
    // std::cout << "corrected pt " << mu1.Pt() << "  " << mu2.Pt() << std::endl;

    double rand = gRandom->Uniform(1);
    double deltaMcSF1 = rc.kSpreadMCerror(q1, mu1.Pt(), mu1.Eta(), mu1.Phi(), genMuonPt1);
    double deltaMcSF2 = rc.kSpreadMCerror(q2, mu1.Pt(), mu1.Eta(), mu1.Phi(), genMuonPt2);
    mu1u*=(1+deltaMcSF1);
    mu1d*=(1-deltaMcSF1);
    mu2u*=(1+deltaMcSF2);
    mu2d*=(1-deltaMcSF2);


	  Double_t lp1 = mu1.Pt();
	  Double_t lp2 = mu2.Pt();

	  // TLorentzVector mu1, mu2;
    
    orderByLepPt(mu1, mu2, q1, q2, mu_MASS);

	  double mll=(mu1+mu2).M();
	  Double_t effdata, effmc;
	  Double_t corr=1;
	  Double_t corrFSR=1, corrMC=1, corrBkg=1, corrTag=1;
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

    // corr = effs.computeHLTSF(&mu1,q1,&mu2,q2)*effs.computeSelSF(&mu1,q1,&mu2,q2)*effs.computeStaSF(&mu1,q1,&mu2,q2);
    corr = effs.fullEfficiencies(&mu1,q1,&mu2,q2);
    // cout << corr1 << " " << corr << endl;
    
    vector<double> uncs_sta = effs.getUncSta(&mu1,q1,&mu2,q2);
    vector<double> uncs_sit = effs.getUncSel(&mu1,q1,&mu2,q2);
    
    corrFSR *= uncs_sta[0]*uncs_sit[0]*effs.computeHLTSF(&mu1,q1,&mu2,q2); // alternate fsr model
    corrMC  *= uncs_sta[1]*uncs_sit[1]*effs.computeHLTSF(&mu1,q1,&mu2,q2); // alternate mc gen model
    corrBkg *= uncs_sta[2]*uncs_sit[2]*effs.computeHLTSF(&mu1,q1,&mu2,q2); // alternate bkg model
    corrTag *= uncs_sta[3]*uncs_sit[3]*effs.computeHLTSF(&mu1,q1,&mu2,q2); // alternate bkg model
    // corr *= effdata/effmc; // orig
    
	  
	  double var=0.;        
    var += effs.statUncSta(&mu1, q1, hErr, hErr, fabs(weight)*corr);
    var += effs.statUncSta(&mu2, q2, hErr, hErr, fabs(weight)*corr);
    var += effs.statUncSel(&mu1, q1, hErr, hErr, fabs(weight)*corr);
    var += effs.statUncSel(&mu2, q2, hErr, hErr, fabs(weight)*corr);
    // // // // // // // var += effs.statUncHLT(&mu1, q1, hErr, hErr, fabs(weight)*corr); 
    // // // // // // // var += effs.statUncHLT(&mu2, q2, hErr, hErr, fabs(weight)*corr);
    var += effs.statUncHLTDilep(&mu1, q1, &mu2, q2);
    // cout << var1 << " " << var << endl;
	  
	  corrUp=corr+sqrt(var);
	  corrDown=corr-sqrt(var);  

    // corr=1;

	  mass = (mu1+mu2).M();
	  pt = (mu1+mu2).Pt();
	  // rapidity = (mu1+mu2).Rapidity();

    if(mass        < MASS_LOW)  continue;
    if(mass        > MASS_HIGH) continue;
    if(mu1.Pt()        < PT_CUT)    continue;
    if(mu2.Pt()        < PT_CUT)    continue;

	  if(typev[ifile]==eZmm) {//mcUp,mcDown,fsrUp,fsrDown,bkgUp,bkgDown,tagptUp,tagptDown,effsUp,effsDown,lepsfUp,lepsfDown,pfireUp,pfireDown
	      yield_zmm += weight*corr;
	      yield_zmm_unc += weight*weight*corr*corr;
        yield_zmm_up += weight*corrUp;
        yield_zmm_dn += weight*corrDown;
        yield_zmm_noPrefire += scale1fb*lumi*corr;
        yield_zmm_pfPhoton += scale1fb*lumi*corr*prefirePhoton;
        yield_zmm_pfJet += scale1fb*lumi*corr*prefireJet;
        if(genVMass<MASS_LOW || genVMass>MASS_HIGH) yield_wm+= weight*corr;
        
        hZmmUnc[mcUp]->Fill(mass,weight*corrMC);
        hZmmUnc[mcDown]->Fill(mass,weight*(corr+(corr-corrMC)));
        hZmmUnc[fsrUp]->Fill(mass,weight*corrFSR);
        hZmmUnc[fsrDown]->Fill(mass,weight*(corr+(corr-corrFSR)));
        hZmmUnc[bkgUp]->Fill(mass,weight*corrBkg);
        hZmmUnc[bkgDown]->Fill(mass,weight*(corr+(corr-corrBkg)));
        hZmmUnc[tagptUp]->Fill(mass,weight*corrTag);
        hZmmUnc[tagptDown]->Fill(mass,weight*(corr+(corr-corrTag)));
        hZmmUnc[effsUp]->Fill(mass,weight*(corr+sqrt(var)));
        hZmmUnc[effsDown]->Fill(mass,weight*(corr-sqrt(var)));
        // roch up/down
        hZmmUnc[lepsfUp]->Fill((mu1u+mu2u).M(),weight*corr);
        hZmmUnc[lepsfDown]->Fill((mu1d+mu2d).M(),weight*corr);
        
        hZmmUnc[pfireUp]->Fill(mass,prefireUp*corr);
        hZmmUnc[pfireDown]->Fill(mass,prefireDown*corr);
        
	      hZmm->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);
	      hZmmNPV->Fill(npv,weight*corr); 
	      hMCNPV->Fill(npv,weight*corr);
	      hZmmZPt->Fill(pt,weight*corr); 
	      hMCZPt->Fill(pt,weight*corr);
	      hZmmLep1Pt->Fill(mu1.Pt(),weight*corr); 
	      hMCLep1Pt->Fill(mu1.Pt(),weight*corr);
    } if(typev[ifile]==eZxx){
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
        hZxxUnc[effsDown]->Fill(mass,weight*(corr-sqrt(var)));
        
        // roch up/down
        hZxxUnc[lepsfUp]->Fill((mu1u+mu2u).M(),weight*corr);
        hZxxUnc[lepsfDown]->Fill((mu1d+mu2d).M(),weight*corr);
        
        hZxxUnc[pfireUp]->Fill(mass,prefireUp*corr);
        hZxxUnc[pfireDown]->Fill(mass,prefireDown*corr);
    } if(typev[ifile]==eWx){
        // cout << "wx? " << nWx << endl;
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
        hWxUnc[lepsfUp]->Fill((mu1u+mu2u).M(),weight*corr);
        hWxUnc[lepsfDown]->Fill((mu1d+mu2d).M(),weight*corr);
        
        hWxUnc[pfireUp]->Fill(mass,prefireUp*corr);
        hWxUnc[pfireDown]->Fill(mass,prefireDown*corr);
      
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
        hDibUnc[lepsfUp]->Fill((mu1u+mu2u).M(),weight*corr);
        hDibUnc[lepsfDown]->Fill((mu1d+mu2d).M(),weight*corr);
        
        hDibUnc[pfireUp]->Fill(mass,prefireUp*corr);
        hDibUnc[pfireDown]->Fill(mass,prefireDown*corr);
        // cout << "blah " << endl;
    } if(typev[ifile]==eEWK || typev[ifile]==eDib || typev[ifile]==eWx|| typev[ifile]==eZxx) {
	      yield_ewk += weight*corr;
	      yield_ewk_unc += weight*weight*corr*corr;
	      hEWK->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);

	      hEWKNPV->Fill(npv,weight*corr); 
	      hMCNPV->Fill(npv,weight*corr);

	      hEWKZPt->Fill(pt,weight*corr); 
	      hEWKZPt_EffBin->Fill(pt,weight*corr2Bin); 
	      hEWKZPt_EffStatUp->Fill(pt,weight*corrUp);
	      hEWKZPt_EffStatDown->Fill(pt,weight*corrDown);
	      hEWKZPt_EffSigShape->Fill(pt,weight*corrSigShape);
	      hEWKZPt_EffBkgShape->Fill(pt,weight*corrBkgShape);
	      hMCZPt->Fill(pt,weight*corr);

	      hEWKLep1Pt->Fill(mu1.Pt(),weight*corr);
	      hEWKLep1Pt_EffBin->Fill(mu1.Pt(),weight*corr2Bin);
	      hEWKLep1Pt_EffStatUp->Fill(mu1.Pt(),weight*corrUp);
	      hEWKLep1Pt_EffStatDown->Fill(mu1.Pt(),weight*corrDown);
	      hEWKLep1Pt_EffSigShape->Fill(mu1.Pt(),weight*corrSigShape);
	      hEWKLep1Pt_EffBkgShape->Fill(mu1.Pt(),weight*corrBkgShape);
	      hMCLep1Pt->Fill(mu1.Pt(),weight*corr);

	    } if(typev[ifile]==eTop) {
	      yield_top += weight*corr;
	      yield_top_unc += weight*weight*corr*corr;
	      hTtb->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);
        
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
        hTtbUnc[lepsfUp]->Fill((mu1u+mu2u).M(),weight*corr);
        hTtbUnc[lepsfDown]->Fill((mu1d+mu2d).M(),weight*corr);
        
        hTtbUnc[pfireUp]->Fill(mass,prefireUp*corr);
        hTtbUnc[pfireDown]->Fill(mass,prefireDown*corr);
      
        

	      hTtbNPV->Fill(npv,weight*corr); 
	      hMCNPV->Fill(npv,weight*corr);

	      hTtbZPt->Fill(pt,weight*corr); 
	      hTtbZPt_EffBin->Fill(pt,weight*corr2Bin); 
	      hTtbZPt_EffStatUp->Fill(pt,weight*corrUp);
	      hTtbZPt_EffStatDown->Fill(pt,weight*corrDown);
	      hTtbZPt_EffSigShape->Fill(pt,weight*corrSigShape); 
	      hTtbZPt_EffBkgShape->Fill(pt,weight*corrBkgShape); 
	      hMCZPt->Fill(pt,weight*corr);

	      hTtbLep1Pt->Fill(mu1.Pt(),weight*corr);
	      hTtbLep1Pt_EffBin->Fill(mu1.Pt(),weight*corr2Bin);
	      hTtbLep1Pt_EffStatUp->Fill(mu1.Pt(),weight*corrUp);
	      hTtbLep1Pt_EffStatDown->Fill(mu1.Pt(),weight*corrDown);
	      hTtbLep1Pt_EffSigShape->Fill(mu1.Pt(),weight*corrSigShape);
	      hTtbLep1Pt_EffBkgShape->Fill(mu1.Pt(),weight*corrBkgShape);
	      hMCLep1Pt->Fill(mu1.Pt(),weight*corr);

          }
        }
      // }
    }
    
    delete infile;
    infile=0, intree=0;
  } 

  outFile->cd();

  hDataZPt->Write();
  hEWKZPt->Write();
  hEWKZPt_EffBin->Write();
  hEWKZPt_EffStatUp->Write();
  hEWKZPt_EffStatDown->Write();
  hEWKZPt_EffSigShape->Write();
  hEWKZPt_EffBkgShape->Write();
  hTtbZPt->Write();
  hTtbZPt_EffBin->Write();
  hTtbZPt_EffStatUp->Write();
  hTtbZPt_EffStatDown->Write();
  hTtbZPt_EffSigShape->Write();
  hTtbZPt_EffBkgShape->Write();

  hDataLep1Pt->Write();
  hEWKLep1Pt->Write();
  hEWKLep1Pt_EffBin->Write();
  hEWKLep1Pt_EffStatUp->Write();
  hEWKLep1Pt_EffStatDown->Write();
  hEWKLep1Pt_EffSigShape->Write();
  hEWKLep1Pt_EffBkgShape->Write();
  hTtbLep1Pt->Write();
  hTtbLep1Pt_EffBin->Write();
  hTtbLep1Pt_EffStatUp->Write();
  hTtbLep1Pt_EffStatDown->Write();
  hTtbLep1Pt_EffSigShape->Write();
  hTtbLep1Pt_EffBkgShape->Write();

  outFile->Write();
  outFile->Close(); 
  
  TString histfname = outputDir + TString("/Zmumu_Histograms.root");
  TFile *histFile = new TFile(histfname,"RECREATE");
  histFile->cd();
  hData->Write();
  hZmm->Write();
  hZxx->Write();
  hWx->Write();
  hTtb->Write();
  hDib->Write();
  
  for(int j = 0; j < nWeight; ++j){
    // hDataUnc[j]->Write();
    hZmmUnc[j]->Write();
    hZxxUnc[j]->Write();
    hWxUnc[j]->Write();
    hTtbUnc[j]->Write();
    hDibUnc[j]->Write();
  }
  
  histFile->Write();
  histFile->Close();

  double MCscale=hData->Integral()/hMC->Integral();

  if(normToData)
    {

      cout<<"Normalized to data: "<<MCscale<<endl;
      
      hZmm->Scale(MCscale);
      hMC->Scale(MCscale);
      hEWK->Scale(MCscale);
      hTtb->Scale(MCscale);
      hZmmNPV->Scale(MCscale);
      hMCNPV->Scale(MCscale);
      hEWKNPV->Scale(MCscale);
      hTtbNPV->Scale(MCscale);
      hZmmZPt->Scale(MCscale);
      hMCZPt->Scale(MCscale);
      hEWKZPt->Scale(MCscale);
      hTtbZPt->Scale(MCscale);
      hZmmLep1Pt->Scale(MCscale);
      hMCLep1Pt->Scale(MCscale);
      hEWKLep1Pt->Scale(MCscale);
      hTtbLep1Pt->Scale(MCscale);
    }

  std::cout << hData->Integral() << std::endl;
  std::cout << hMC->Integral() << std::endl;
  std::cout << hZmm->Integral() << std::endl;
  std::cout << hEWK->Integral() << std::endl;
  std::cout << hTtb->Integral() << std::endl;

  for(int j=0;j!=hDataZPt->GetNbinsX();++j)
    {
      hDataZPt->SetBinContent(j+1,hDataZPt->GetBinContent(j+1)/hDataZPt->GetBinWidth(j+1));
      hDataZPt->SetBinError(j+1,hDataZPt->GetBinError(j+1)/hDataZPt->GetBinWidth(j+1));
      hMCZPt->SetBinContent(j+1,hMCZPt->GetBinContent(j+1)/hMCZPt->GetBinWidth(j+1));
      hMCZPt->SetBinError(j+1,hMCZPt->GetBinError(j+1)/hMCZPt->GetBinWidth(j+1));
      hZmmZPt->SetBinContent(j+1,hZmmZPt->GetBinContent(j+1)/hZmmZPt->GetBinWidth(j+1));
      hZmmZPt->SetBinError(j+1,hZmmZPt->GetBinError(j+1)/hZmmZPt->GetBinWidth(j+1));
      hEWKZPt->SetBinContent(j+1,hEWKZPt->GetBinContent(j+1)/hEWKZPt->GetBinWidth(j+1));
      hEWKZPt->SetBinError(j+1,hEWKZPt->GetBinError(j+1)/hEWKZPt->GetBinWidth(j+1));
      hTtbZPt->SetBinContent(j+1,hTtbZPt->GetBinContent(j+1)/hTtbZPt->GetBinWidth(j+1));
      hTtbZPt->SetBinError(j+1,hTtbZPt->GetBinError(j+1)/hTtbZPt->GetBinWidth(j+1));
    }

  for(int j=0;j!=hDataLep1Pt->GetNbinsX();++j)
    {
      hDataLep1Pt->SetBinContent(j+1,hDataLep1Pt->GetBinContent(j+1)/hDataLep1Pt->GetBinWidth(j+1));
      hDataLep1Pt->SetBinError(j+1,hDataLep1Pt->GetBinError(j+1)/hDataLep1Pt->GetBinWidth(j+1));
      hMCLep1Pt->SetBinContent(j+1,hMCLep1Pt->GetBinContent(j+1)/hMCLep1Pt->GetBinWidth(j+1));
      hMCLep1Pt->SetBinError(j+1,hMCLep1Pt->GetBinError(j+1)/hMCLep1Pt->GetBinWidth(j+1));
      hZmmLep1Pt->SetBinContent(j+1,hZmmLep1Pt->GetBinContent(j+1)/hZmmLep1Pt->GetBinWidth(j+1));
      hZmmLep1Pt->SetBinError(j+1,hZmmLep1Pt->GetBinError(j+1)/hZmmLep1Pt->GetBinWidth(j+1));
      hEWKLep1Pt->SetBinContent(j+1,hEWKLep1Pt->GetBinContent(j+1)/hEWKLep1Pt->GetBinWidth(j+1));
      hEWKLep1Pt->SetBinError(j+1,hEWKLep1Pt->GetBinError(j+1)/hEWKLep1Pt->GetBinWidth(j+1));
      hTtbLep1Pt->SetBinContent(j+1,hTtbLep1Pt->GetBinContent(j+1)/hTtbLep1Pt->GetBinWidth(j+1));
      hTtbLep1Pt->SetBinError(j+1,hTtbLep1Pt->GetBinError(j+1)/hTtbLep1Pt->GetBinWidth(j+1));
    }

  
  TH1D *hZmumuDiff = makeDiffHist(hData,hMC,"hZmumuDiff");
  hZmumuDiff->SetMarkerStyle(kFullCircle); 
  hZmumuDiff->SetMarkerSize(0.9);

  TH1D *hZmumuNPVDiff = makeDiffHist(hDataNPV,hMCNPV,"hZmumuNPVDiff");
  hZmumuNPVDiff->SetMarkerStyle(kFullCircle); 
  hZmumuNPVDiff->SetMarkerSize(0.9);

  TH1D *hZmumuZPtDiff = makeDiffHist(hDataZPt,hMCZPt,"hZmumuZPtDiff");
  hZmumuZPtDiff->SetMarkerStyle(kFullCircle); 
  hZmumuZPtDiff->SetMarkerSize(0.9);

  TH1D *hZmumuLep1PtDiff = makeDiffHist(hDataLep1Pt,hMCLep1Pt,"hZmumuLep1PtDiff");
  hZmumuLep1PtDiff->SetMarkerStyle(kFullCircle); 
  hZmumuLep1PtDiff->SetMarkerSize(0.9);
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  char ylabel[100];     // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000);  
  // sprintf(lumitext,"%.1f fb^{-1}  (5 TeV)",lumi/1000);  
  
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
  
  sprintf(ylabel,"Events / %.1f GeV",hData->GetBinWidth(1));
  CPlot plotZmumu("zmm"+norm,"","",ylabel);
  plotZmumu.AddHist1D(hData,"data","E");
  plotZmumu.AddToStack(hZmm,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumu.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumu.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumu.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumu.SetYRange(0.01,1.2*(hData->GetMaximum() + sqrt(hData->GetMaximum())));
  plotZmumu.TransLegend(0.1,-0.05);
  plotZmumu.Draw(c,kFALSE,format,1);

  CPlot plotZmumuDiff("zmm"+norm,"","M(#mu^{+}#mu^{-}) [GeV]","#frac{Data-Pred}{Data}");
  plotZmumuDiff.AddHist1D(hZmumuDiff,"EX0",ratioColor);
  plotZmumuDiff.SetYRange(-0.2,0.2);
  plotZmumuDiff.AddLine(MASS_LOW, 0,MASS_HIGH, 0,kBlack,1);
  plotZmumuDiff.AddLine(MASS_LOW, 0.1,MASS_HIGH, 0.1,kBlack,3);
  plotZmumuDiff.AddLine(MASS_LOW,-0.1,MASS_HIGH,-0.1,kBlack,3);
  plotZmumuDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumu2("zmmlog"+norm,"","",ylabel);
  plotZmumu2.AddHist1D(hData,"data","E");
  plotZmumu2.AddToStack(hEWK,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumu2.AddToStack(hTtb,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumu2.AddToStack(hZmm,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumu2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumu2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumu2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZmumu2.SetLogy();
  plotZmumu2.SetYRange(1e-4*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plotZmumu2.TransLegend(0.1,-0.05);
  plotZmumu2.Draw(c,kTRUE,format,1);

  //
  // NPV
  // 

  sprintf(ylabel,"Events");
  CPlot plotZmumuNPV("zmmNPV"+norm,"","",ylabel);
  plotZmumuNPV.AddHist1D(hDataNPV,"data","E");
  plotZmumuNPV.AddToStack(hZmmNPV,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuNPV.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuNPV.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuNPV.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuNPV.SetYRange(0.01,1.2*(hData->GetMaximum() + sqrt(hData->GetMaximum())));
  plotZmumuNPV.TransLegend(0.1,-0.05);
  plotZmumuNPV.Draw(c,kFALSE,format,1);

  CPlot plotZmumuNPVDiff("zmmNPV"+norm,"","npv","#frac{Data-Pred}{Data}");
  plotZmumuNPVDiff.AddHist1D(hZmumuNPVDiff,"EX0",ratioColor);
  plotZmumuNPVDiff.SetYRange(-0.2,0.2);
  plotZmumuNPVDiff.AddLine(0, 0,50, 0,kBlack,1);
  plotZmumuNPVDiff.AddLine(0, 0.1,50, 0.1,kBlack,3);
  plotZmumuNPVDiff.AddLine(0,-0.1,50,-0.1,kBlack,3);
  plotZmumuNPVDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuNPV2("zmmNPVlog"+norm,"","",ylabel);
  plotZmumuNPV2.AddHist1D(hDataNPV,"data","E");
  plotZmumuNPV2.AddToStack(hEWKNPV,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuNPV2.AddToStack(hTtbNPV,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuNPV2.AddToStack(hZmmNPV,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuNPV2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuNPV2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuNPV2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZmumuNPV2.SetLogy();
  plotZmumuNPV2.SetYRange(1e-4*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plotZmumuNPV2.TransLegend(0.1,-0.05);
  plotZmumuNPV2.Draw(c,kTRUE,format,1);

  //
  // Z Pt
  //   
  sprintf(ylabel,"Events / 1 GeV");
  CPlot plotZmumuZPt("zmmZPt"+norm,"","",ylabel);
  plotZmumuZPt.AddHist1D(hDataZPt,"data","E");
  plotZmumuZPt.AddToStack(hZmmZPt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuZPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuZPt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuZPt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuZPt.SetLogx();
  plotZmumuZPt.SetLogy(0);
  plotZmumuZPt.SetYRange(0.01,1.2*(hDataZPt->GetMaximum() + sqrt(hDataZPt->GetMaximum())));
  plotZmumuZPt.TransLegend(0.1,-0.05);
  plotZmumuZPt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuZPtDiff("zmmZPt"+norm,"","p_{T}^{#mu^{+}#mu^{-}} [GeV]","#frac{Data-Pred}{Data}");
  plotZmumuZPtDiff.AddHist1D(hZmumuZPtDiff,"EX0",ratioColor);
  plotZmumuZPtDiff.SetLogx();
  plotZmumuZPtDiff.SetYRange(-0.2,0.2);
  plotZmumuZPtDiff.AddLine(0, 0,1000, 0,kBlack,1);
  plotZmumuZPtDiff.AddLine(0, 0.1,1000, 0.1,kBlack,3);
  plotZmumuZPtDiff.AddLine(0,-0.1,1000,-0.1,kBlack,3);
  plotZmumuZPtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuZPt2("zmmZPtlog"+norm,"","",ylabel);
  plotZmumuZPt2.AddHist1D(hDataZPt,"data","E");
  plotZmumuZPt2.AddToStack(hEWKZPt,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuZPt2.AddToStack(hTtbZPt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuZPt2.AddToStack(hZmmZPt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);\
  plotZmumuZPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuZPt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuZPt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZmumuZPt2.SetLogx();
  plotZmumuZPt2.SetLogy();
  plotZmumuZPt2.SetYRange(1e-6*(hDataZPt->GetMaximum()),10*(hDataZPt->GetMaximum()));
  plotZmumuZPt2.TransLegend(0.1,-0.05);
  plotZmumuZPt2.Draw(c,kTRUE,format,1);

 //
  // Lep1 Pt
  //   
  sprintf(ylabel,"Events / 1 GeV");
  CPlot plotZmumuLep1Pt("zmmLep1Pt"+norm,"","",ylabel);
  plotZmumuLep1Pt.AddHist1D(hDataLep1Pt,"data","E");
  plotZmumuLep1Pt.AddToStack(hZmmLep1Pt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLep1Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep1Pt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLep1Pt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuLep1Pt.SetLogx();
  plotZmumuLep1Pt.SetLogy(0);
  plotZmumuLep1Pt.SetYRange(0.01,1.2*(hDataLep1Pt->GetMaximum() + sqrt(hDataLep1Pt->GetMaximum())));
  plotZmumuLep1Pt.TransLegend(0.1,-0.05);
  plotZmumuLep1Pt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLep1PtDiff("zmmLep1Pt"+norm,"","p_{T}(leading muon) [GeV]","#frac{Data-Pred}{Data}");
  plotZmumuLep1PtDiff.AddHist1D(hZmumuLep1PtDiff,"EX0",ratioColor);
  plotZmumuLep1PtDiff.SetLogx();
  plotZmumuLep1PtDiff.SetYRange(-0.2,0.2);
  plotZmumuLep1PtDiff.AddLine(25, 0,300, 0,kBlack,1);
  plotZmumuLep1PtDiff.AddLine(25, 0.1,300, 0.1,kBlack,3);
  plotZmumuLep1PtDiff.AddLine(25,-0.1,300,-0.1,kBlack,3);
  plotZmumuLep1PtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLep1Pt2("zmmLep1Ptlog"+norm,"","",ylabel);
  plotZmumuLep1Pt2.AddHist1D(hDataLep1Pt,"data","E");
  plotZmumuLep1Pt2.AddToStack(hEWKLep1Pt,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuLep1Pt2.AddToStack(hTtbLep1Pt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuLep1Pt2.AddToStack(hZmmLep1Pt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLep1Pt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep1Pt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLep1Pt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZmumuLep1Pt2.SetLogx();
  plotZmumuLep1Pt2.SetLogy();
  plotZmumuLep1Pt2.SetYRange(1e-5*(hDataLep1Pt->GetMaximum()),10*(hDataLep1Pt->GetMaximum()));
  plotZmumuLep1Pt2.TransLegend(0.1,-0.05);
  plotZmumuLep1Pt2.Draw(c,kTRUE,format,1);

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  


  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  cout << endl;

  cout << " The Zmm event yield is " << yield << " +/-" << sqrt(yield) << "." << endl;
  cout << " The Zmm expected event yield is " << yield_zmm << " +/-" << sqrt(yield_zmm_unc) << "." << endl;
  cout << " The Zmm expected event yield (no. prefire) " << yield_zmm_noPrefire << " +/-" << sqrt(yield_zmm_unc) << "." << endl;
  cout << " The Zmm stat up event yield is " << yield_zmm_up << " +/-" << sqrt(yield_zmm_unc) << "." << endl;
  cout << " The Zmm stat down event yield is " << yield_zmm_dn << " +/-" << sqrt(yield_zmm_unc) << "." << endl;
  cout << " The EWK event yield is " << yield_ewk << " +/-" << sqrt(yield_ewk_unc) << "." << endl;
  cout << "  -> Dib event yield is " << nDib      << " +/-" << sqrt(nDibUnc)       << "." << endl;
  cout << "  -> Zxx event yield is " << nZxx      << " +/-" << sqrt(nZxxUnc)       << "." << endl;
  cout << "  -> Wx  event yield is " << nWx       << " +/-" << sqrt(nWxUnc)        << "." << endl;
  cout << " The Top event yield is " << yield_top << " +/-" << sqrt(yield_top_unc) << "." << endl;
  cout << " The Zmm Data Yield w/ ewk&top removed: " << yield - yield_ewk - yield_top << endl;
  cout << " Zmm yield w/ all bkg removed: " << yield - yield_ewk - yield_top - yield_wm << endl;
  cout << " Prefire scale factor :" << yield_zmm_noPrefire/yield_zmm << endl;  
  cout << " Prefire Jets only : " << yield_zmm_pfJet << "  scale fac: " << yield_zmm_noPrefire/yield_zmm_pfJet << endl;
  cout << " Prefire Photons only : " << yield_zmm_pfPhoton << "  scale fac: " << yield_zmm_noPrefire/yield_zmm_pfPhoton << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  
  ofstream txtfile;
  char txtfname[100];  
  sprintf(txtfname,"%s/zmm_yields.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;  
  txtfile << endl;

  txtfile << " The Zmm event yield is " << yield << " +/-" << sqrt(yield) << "." << endl;
  txtfile << " The Zmm expected event yield is " << yield_zmm << " +/-" << sqrt(yield_zmm_unc) << "." << endl;
  txtfile << " The Zmm expected event yield (no. prefire) " << yield_zmm_noPrefire << " +/-" << sqrt(yield_zmm_unc) << "." << endl;
  txtfile << " The Zmm stat up event yield is " << yield_zmm_up << " +/-" << sqrt(yield_zmm_unc) << "." << endl;
  txtfile << " The Zmm stat down event yield is " << yield_zmm_dn << " +/-" << sqrt(yield_zmm_unc) << "." << endl;
  txtfile << " The EWK event yield is " << yield_ewk << " +/-" << sqrt(yield_ewk_unc) << "." << endl;
  txtfile << "  -> Dib event yield is " << nDib      << " +/-" << sqrt(nDibUnc)       << "." << endl;
  txtfile << "  -> Zxx event yield is " << nZxx      << " +/-" << sqrt(nZxxUnc)       << "." << endl;
  txtfile << "  -> Wx  event yield is " << nWx       << " +/-" << sqrt(nWxUnc)        << "." << endl;
  txtfile << " The Top event yield is " << yield_top << " +/-" << sqrt(yield_top_unc) << "." << endl;
  txtfile << " The Zmm Data Yield w/ ewk&top removed: " << yield - yield_ewk - yield_top << endl;
  txtfile << " Zmm yield w/ all bkg removed: " << yield - yield_ewk - yield_top - yield_wm << endl;
  txtfile << " Prefire scale factor :" << yield_zmm_noPrefire/yield_zmm << endl;
  txtfile << " Prefire Jets only : " << yield_zmm_pfJet << "  scale fac: " << yield_zmm_noPrefire/yield_zmm_pfJet << endl;
  txtfile << " Prefire Photons only : " << yield_zmm_pfPhoton << "  scale fac: " << yield_zmm_noPrefire/yield_zmm_pfPhoton << endl;
  txtfile << std::endl;
  txtfile.close();
  
  

  gBenchmark->Show("plotZmm");
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

void orderByLepPt(TLorentzVector &mu1, TLorentzVector &mu2,Int_t &q1, Int_t &q2, double mu_MASS){
  if(mu2.Pt() > mu1.Pt()){
    TLorentzVector tmp;
    tmp.SetPtEtaPhiM(mu1.Pt(), mu1.Eta(), mu1.Phi(), mu_MASS);
    mu1 = mu2;
    mu2 = tmp;

    q1*= -1;
    q2*= -1;
  }
  return;
}