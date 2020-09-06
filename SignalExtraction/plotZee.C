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

//=== MAIN MACRO ================================================================================================= 

void plotZee(const TString  inputDir,    // input directory
	     const TString  outputDir,   // output directory
             const TString sqrts,
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
  enum { eData, eZee, eEWK, eTop};  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;

  fnamev.push_back(inputDir + TString("/") + TString("data_select.root")); typev.push_back(eData);
  fnamev.push_back(inputDir + TString("/") + TString("zee_select.root"));   typev.push_back(eZee);
  // fnamev.push_back(inputDir + TString("/") + TString("ewk_select.root"));  typev.push_back(eEWK);
  // fnamev.push_back(inputDir + TString("/") + TString("top_select.root"));  typev.push_back(eTop);
  
  if(sqrts == "5TeV" ){
    fnamev.push_back(inputDir + TString("/") + TString("wx_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("zxx_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("wz_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("ww_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("zz_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("top_select.root"));  typev.push_back(eTop);
  } else {
    fnamev.push_back(inputDir + TString("/") + TString("wx0_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("wx1_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("wx2_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("zxx_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("ww_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("wz_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("zz_select.root"));  typev.push_back(eEWK);
    fnamev.push_back(inputDir + TString("/") + TString("top1_select.root"));  typev.push_back(eTop);
    fnamev.push_back(inputDir + TString("/") + TString("top2_select.root"));  typev.push_back(eTop);
    fnamev.push_back(inputDir + TString("/") + TString("top3_select.root"));  typev.push_back(eTop);
  }
 
  //
  // Fit options
  //
  // const Int_t    NBINS     = 120;
  const Int_t    NBINS     = 60;
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  const Double_t PT_CUT    = 25;
  // const Double_t PT_CUT    = 32;
  // const Double_t ETA_CUT   = 1.444;
  const Double_t ETA_CUT   = 2.4;
    
  // const Double_t ECAL_GAP_LOW  = 1.4442;
  // const Double_t ECAL_GAP_HIGH = 1.566;
  const Double_t ECAL_GAP_LOW  = 10.;
  const Double_t ECAL_GAP_HIGH = 10.;
  
  // efficiency files

  TString baseDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_"+sqrts+"/results/Zee/";
  // std::string baseDirs = baseDir.Data();
  AppEffSF effs(baseDir);
  effs.loadHLT("EleHLTEff_aMCxPythia","Positive","Negative");
  effs.loadSel("EleGSFSelEff_aMCxPythia","Combined","Combined");
  // effs.loadSta("MuStaEff_aMCxPythia","Combined","Combined");

  const TString corrFiles = "../EleScale/Run2017_LowPU_v2";
  EnergyScaleCorrection ec( corrFiles.Data(),EnergyScaleCorrection::ECALELF);

  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zee_DataBkg.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);

  TFile *rajdeepData = new TFile("root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/MINIAODNTUPLE/94X_dataRun2_ReReco_EOY17_v6/HighEGJet-Run2017H-noSkim-17Nov2017-ReReco-miniAOD/306896-307082/306896-307082_13TeV_LowPU/Moriond18MC1042s/HighEGJet-Run2017H-noSkim-17Nov2017-ReReco-miniAOD-allRange.root","OPEN");
  TFile *rajdeepMC = new TFile("/eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/MINIAODNTUPLE/94X_mc2017_realistic_v10/DYJetsToLL_fixECALGT_LowPU/allRange/MCfixECAL/DYJetsToLL_fixECALGT_LowPU-allRange.root","OPEN");
  
  TTree *treeData = (TTree*)rajdeepData->Get("selected");// assert(intree);
  TTree *treeMC = (TTree*)rajdeepMC->Get("selected");// assert(intree);

  TH1D *hDataEG = new TH1D("hDataEG","",NBINS,MASS_LOW,MASS_HIGH); //hDataEG->Sumw2();
  TH1D *hZeeEG  = new TH1D("hZeeEG", "",NBINS,MASS_LOW,MASS_HIGH); //hZeeEG->Sumw2();
  // treeData->Draw("invMass_ECAL_ele>>hDataEGa");
  Float_t invMass=0;
  UInt_t runNumber=0;
  UShort_t lumiBlock=0;
  ULong64_t eventNumber=0;
  Float_t invMass_ECAL_ele=0;
  Float_t etaEle[3];
  Float_t phiEle[3];
  Float_t R9Ele[3];
  Float_t etaSCEle[3];
  Float_t energy_ECAL_ele[3];
  Float_t ele1E;
  Float_t ele2E;
  TLorentzVector clep1, clep2;
  TLorentzVector eclep1, eclep2;
  TLorentzVector dl;
  // for(int i =0; i < treeData->GetEntriesFast();i++){
  for(int i =0; i < 1000;i++){
    if(i%100000==0) cout << "Processing event " << i << ". " << (double)i/(double)treeData->GetEntries()*100 << " percent done with this file." << endl;
    treeData -> SetBranchStatus("runNumber",1);
    treeData -> SetBranchStatus("lumiBlock",1);
    treeData -> SetBranchStatus("eventNumber",1);
    treeData -> SetBranchStatus("invMass",1);
    treeData -> SetBranchStatus("invMass_ECAL_ele",1);
    treeData -> SetBranchStatus("etaEle",1);
    treeData -> SetBranchStatus("phiEle",1);
    treeData -> SetBranchStatus("R9Ele",1);
    treeData -> SetBranchStatus("etaSCEle",1);
    treeData -> SetBranchStatus("energy_ECAL_ele",1);
    treeData -> SetBranchStatus("ele1E",1);
    treeData -> SetBranchStatus("ele2E",1);
    treeData->SetBranchAddress("invMass",       &invMass);        // sc2 4-vector
    treeData->SetBranchAddress("runNumber",       &runNumber);        // sc2 4-vector
    treeData->SetBranchAddress("lumiBlock",       &lumiBlock);        // sc2 4-vector
    treeData->SetBranchAddress("eventNumber",       &eventNumber);        // sc2 4-vector
    treeData->SetBranchAddress("invMass_ECAL_ele",       &invMass_ECAL_ele);        // sc2 4-vector
    treeData->SetBranchAddress("etaEle",       &etaEle);        // sc2 4-vector
    treeData->SetBranchAddress("phiEle",       &phiEle);        // sc2 4-vector
    treeData->SetBranchAddress("R9Ele",        &R9Ele);        // sc2 4-vector
    treeData->SetBranchAddress("etaSCEle",       &etaSCEle);        // sc2 4-vector
    treeData->SetBranchAddress("energy_ECAL_ele",       &energy_ECAL_ele);        // sc2 4-vector
    treeData->SetBranchAddress("ele1E",       &ele1E);// sc2 4-vector
    treeData->SetBranchAddress("ele2E",       &ele2E);        // sc2 4-vector
    treeData->GetEntry(i);
    // let's look at run 30717
    if(runNumber!=307017)continue;// || lumiBlock != 81) continue;
    // if(fabs(etaEle[0])>ETA_CUT||fabs(etaEle[1])>ETA_CUT) continue;
    if(fabs(etaEle[0])>ETA_CUT||fabs(etaEle[1])>ETA_CUT) continue;
    if(fabs(etaEle[0])>=ECAL_GAP_LOW && fabs(etaEle[0])<=ECAL_GAP_HIGH) continue;
    if(fabs(etaEle[1])>=ECAL_GAP_LOW && fabs(etaEle[1])<=ECAL_GAP_HIGH) continue;
    // std::cout << "eta ele " << fabs(etaEle[0]) << " " << fabs(etaEle[1]) << std::endl;
    // if(energy_ECAL_ele[0]/cosh(fabs(etaEle[0]))<PT_CUT||energy_ECAL_ele[1]/cosh(fabs(etaEle[1]))<PT_CUT) continue;
    if(energy_ECAL_ele[0]/cosh(fabs(etaEle[0]))<PT_CUT||energy_ECAL_ele[1]/cosh(fabs(etaEle[1]))<PT_CUT) continue;
    clep1.SetPtEtaPhiM(energy_ECAL_ele[0]/cosh(fabs(etaEle[0])), etaEle[0], phiEle[0],ELE_MASS);
    clep2.SetPtEtaPhiM(energy_ECAL_ele[1]/cosh(fabs(etaEle[1])), etaEle[1], phiEle[1],ELE_MASS);
    eclep1=clep1;
    eclep2=clep2;
    
    double scale = ec.scaleCorr(runNumber, clep1.Pt(), fabs(etaEle[0]), R9Ele[0]);
    clep1*=scale;
    scale = ec.scaleCorr(runNumber, clep2.Pt(), fabs(etaEle[1]), R9Ele[1]);
    clep2*=scale;
    
    if(clep1.Pt() < PT_CUT || clep2.Pt() < PT_CUT) continue;
    dl = clep1+clep2;
    
    double mass = dl.M();
    // hDataEG->Fill(invMass_ECAL_ele);
    // hDataEG->Fill(mass);
  }
  
  // for(int i =0; i < treeMC->GetEntriesFast();i++){
  for(int i =0; i < 1000;i++){
    if(i%100000==0) cout << "Processing event " << i << ". " << (double)i/(double)treeMC->GetEntries()*100 << " percent done with this file." << endl;
    treeMC -> SetBranchStatus("invMass",1);
    treeMC -> SetBranchStatus("runNumber",1);
    treeMC -> SetBranchStatus("invMass_ECAL_ele",1);
    treeMC -> SetBranchStatus("etaEle",1);
    treeMC -> SetBranchStatus("phiEle",1);
    treeMC -> SetBranchStatus("R9Ele",1);
    treeMC -> SetBranchStatus("etaSCEle",1);
    treeMC -> SetBranchStatus("ele1E",1);
    treeMC -> SetBranchStatus("ele2E",1);
    treeMC -> SetBranchStatus("energy_ECAL_ele",1);
    treeMC->SetBranchAddress("invMass",       &invMass);        // sc2 4-vector
    treeMC->SetBranchAddress("runNumber",       &runNumber);        // sc2 4-vector
    treeMC->SetBranchAddress("invMass_ECAL_ele",       &invMass_ECAL_ele);        // sc2 4-vector
    treeMC->SetBranchAddress("etaEle",       &etaEle);        // sc2 4-vector
    treeMC->SetBranchAddress("phiEle",       &phiEle);        // sc2 4-vector
    treeMC->SetBranchAddress("R9Ele",        &R9Ele);        // sc2 4-vector
    treeMC->SetBranchAddress("etaSCEle",       &etaSCEle);        // sc2 4-vector
    treeMC->SetBranchAddress("energy_ECAL_ele",       &energy_ECAL_ele);        // sc2 4-vector
    treeMC->SetBranchAddress("ele1E",       &ele1E);        // sc2 4-vector
    treeMC->SetBranchAddress("ele2E",       &ele2E);        // sc2 4-vector
    treeMC->GetEntry(i);
    // std::cout << "eta ele " << fabs(etaEle[0]) << " " << fabs(etaEle[1]) << std::endl;
    // if(fabs(etaEle[0])>ETA_CUT||fabs(etaEle[1])>ETA_CUT) continue;
    if(fabs(etaEle[0])>ETA_CUT||fabs(etaEle[1])>ETA_CUT) continue;
    if(fabs(etaEle[0])>=ECAL_GAP_LOW && fabs(etaEle[0])<=ECAL_GAP_HIGH) continue;
    if(fabs(etaEle[1])>=ECAL_GAP_LOW && fabs(etaEle[1])<=ECAL_GAP_HIGH) continue;
    // if(energy_ECAL_ele[0]/cosh(fabs(etaEle[0]))<PT_CUT||energy_ECAL_ele[1]/cosh(fabs(etaEle[1]))<PT_CUT) continue;
    if(energy_ECAL_ele[0]/cosh(fabs(etaEle[0]))<PT_CUT||energy_ECAL_ele[1]/cosh(fabs(etaEle[1]))<PT_CUT) continue;
    
    clep1.SetPtEtaPhiM(energy_ECAL_ele[0]/cosh(fabs(etaEle[0])), etaEle[0], phiEle[0],ELE_MASS);
    clep2.SetPtEtaPhiM(energy_ECAL_ele[1]/cosh(fabs(etaEle[1])), etaEle[1], phiEle[1],ELE_MASS);
    // double rand1 = gRandom->Gaus(0,1);
    // double rand2 = gRandom->Gaus(0,1);
    // double tagSmear = ec.smearingSigma(runNumber, clep1.Pt(), fabs(etaEle[0]), R9Ele[0], 12, 0., 0.);
    // clep1*= 1+ rand1*tagSmear;
    // tagSmear = ec.smearingSigma(runNumber, clep2.Pt(), fabs(etaEle[1]), R9Ele[1], 12, 0., 0.);
    // clep2*= 1+rand2*tagSmear;
    dl = clep1+clep2;
    
    if(clep1.Pt() < PT_CUT || clep2.Pt() < PT_CUT) continue;
    double mass = dl.M();
    // hZeeEG->Fill(invMass_ECAL_ele);
    // hZeeEG->Fill(mass);
  }
  
  // // TH1D *hDataEG = (TH1D*)gDirectory->Get("hDataEGa");
  // // treeMC->Draw("invMass_ECAL_ele>>hZeeEGa");
  // // TH1D *hZeeEG = (TH1D*)gDirectory->Get("hZeeEGa");
  
  std::cout << "integrals " <<  hZeeEG->Integral() << " " << hDataEG->Integral() << std::endl;
  
  // plot output file format
  const TString format("all");

  // setup efficiency shape systematics
  // TFile *GsfSelSigSysFile = new TFile(GsfSelEffSignalShapeSys);
  // TH2D *hGsfSelSigSys = (TH2D*)GsfSelSigSysFile->Get("h");
  // TFile *GsfSelBkgSysFile = new TFile(GsfSelEffBackgroundShapeSys);
  // TH2D *hGsfSelBkgSys = (TH2D*)GsfSelBkgSysFile->Get("h");

  Int_t yield = 0;
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
  double ZPtBins[]={0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,200};
  double PhiStarBins[]={0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.01,0.012,0.014,0.016,0.018,0.021,0.024,0.027,0.030,0.034,0.038,0.044,0.050,0.058,0.066,0.076,0.088,0.10,0.12,0.14,0.16,0.18,0.20,0.24,0.28,0.34,0.42,0.52,0.64,0.8,1.0,1.5,2,3};
  double Lep1PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150};
  double Lep2PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150};
  double LepNegPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
  double LepPosPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};


  TH2D *hCompareElePtEcalE = new TH2D("hCompare","",200,0,100,200,0,100);

  // TH1D *hData = new TH1D("hDataEG","",NBINS,MASS_LOW,MASS_HIGH); hData->Sumw2();
  // TH1D *hZee  = new TH1D("hZeeEG", "",NBINS,MASS_LOW,MASS_HIGH); hZee->Sumw2();

  TH1D *hData = new TH1D("hData","",NBINS,MASS_LOW,MASS_HIGH); hData->Sumw2();
  TH1D *hZee  = new TH1D("hZee", "",NBINS,MASS_LOW,MASS_HIGH); hZee->Sumw2();
  TH1D *hEWK  = new TH1D("hEWK", "",NBINS,MASS_LOW,MASS_HIGH); hEWK->Sumw2();
  TH1D *hTop  = new TH1D("hTop", "",NBINS,MASS_LOW,MASS_HIGH); hTop->Sumw2();
  TH1D *hMC   = new TH1D("hMC",  "",NBINS,MASS_LOW,MASS_HIGH); hMC->Sumw2();

  TH1D *hDataNPV = new TH1D("hDataNPV","",50,0,50); hDataNPV->Sumw2();
  TH1D *hZeeNPV  = new TH1D("hZeeNPV", "",50,0,50); hZeeNPV->Sumw2();
  TH1D *hEWKNPV  = new TH1D("hEWKNPV", "",50,0,50); hEWKNPV->Sumw2();
  TH1D *hTopNPV  = new TH1D("hTopNPV", "",50,0,50); hTopNPV->Sumw2();
  TH1D *hMCNPV   = new TH1D("hMCNPV",  "",50,0,50); hMCNPV->Sumw2();


  const int nBinsZPt= sizeof(ZPtBins)/sizeof(double)-1;
  TH1D *hDataZPt = new TH1D("hDataZPt","",nBinsZPt,ZPtBins); hDataZPt->Sumw2();
  TH1D *hZeeZPt  = new TH1D("hZeeZPt", "",nBinsZPt,ZPtBins); hZeeZPt->Sumw2();
  TH1D *hEWKZPt  = new TH1D("hEWKZPt", "",nBinsZPt,ZPtBins); hEWKZPt->Sumw2();
  TH1D *hTopZPt  = new TH1D("hTopZPt", "",nBinsZPt,ZPtBins); hTopZPt->Sumw2();
  TH1D *hMCZPt   = new TH1D("hMCZPt",  "",nBinsZPt,ZPtBins); hMCZPt->Sumw2();

  TH1D *hEWKZPt_EffBin  = new TH1D("hEWKZPt_EffBin", "",nBinsZPt,ZPtBins); hEWKZPt_EffBin->Sumw2();
  TH1D *hTopZPt_EffBin  = new TH1D("hTopZPt_EffBin", "",nBinsZPt,ZPtBins); hTopZPt_EffBin->Sumw2();
  TH1D *hEWKZPt_EffStatUp  = new TH1D("hEWKZPt_EffStatUp", "",nBinsZPt,ZPtBins); hEWKZPt_EffStatUp->Sumw2();
  TH1D *hTopZPt_EffStatUp  = new TH1D("hTopZPt_EffStatUp", "",nBinsZPt,ZPtBins); hTopZPt_EffStatUp->Sumw2();
  TH1D *hEWKZPt_EffStatDown  = new TH1D("hEWKZPt_EffStatDown", "",nBinsZPt,ZPtBins); hEWKZPt_EffStatDown->Sumw2();
  TH1D *hTopZPt_EffStatDown  = new TH1D("hTopZPt_EffStatDown", "",nBinsZPt,ZPtBins); hTopZPt_EffStatDown->Sumw2();
  TH1D *hEWKZPt_EffSigShape  = new TH1D("hEWKZPt_EffSigShape", "",nBinsZPt,ZPtBins); hEWKZPt_EffSigShape->Sumw2();
  TH1D *hTopZPt_EffSigShape  = new TH1D("hTopZPt_EffSigShape", "",nBinsZPt,ZPtBins); hTopZPt_EffSigShape->Sumw2();
  TH1D *hEWKZPt_EffBkgShape  = new TH1D("hEWKZPt_EffBkgShape", "",nBinsZPt,ZPtBins); hEWKZPt_EffBkgShape->Sumw2();
  TH1D *hTopZPt_EffBkgShape  = new TH1D("hTopZPt_EffBkgShape", "",nBinsZPt,ZPtBins); hTopZPt_EffBkgShape->Sumw2();

  const int nBinsPhiStar= sizeof(PhiStarBins)/sizeof(double)-1;
  TH1D *hDataPhiStar = new TH1D("hDataPhiStar","",nBinsPhiStar,PhiStarBins); hDataPhiStar->Sumw2();
  TH1D *hZeePhiStar  = new TH1D("hZeePhiStar", "",nBinsPhiStar,PhiStarBins); hZeePhiStar->Sumw2();
  TH1D *hEWKPhiStar  = new TH1D("hEWKPhiStar", "",nBinsPhiStar,PhiStarBins); hEWKPhiStar->Sumw2();
  TH1D *hTopPhiStar  = new TH1D("hTopPhiStar", "",nBinsPhiStar,PhiStarBins); hTopPhiStar->Sumw2();
  TH1D *hMCPhiStar   = new TH1D("hMCPhiStar",  "",nBinsPhiStar,PhiStarBins); hMCPhiStar->Sumw2();
  TH1D *hEWKPhiStar_EffBin  = new TH1D("hEWKPhiStar_EffBin", "",nBinsPhiStar,PhiStarBins); hEWKPhiStar_EffBin->Sumw2();
  TH1D *hTopPhiStar_EffBin  = new TH1D("hTopPhiStar_EffBin", "",nBinsPhiStar,PhiStarBins); hTopPhiStar_EffBin->Sumw2();
  TH1D *hEWKPhiStar_EffStatUp  = new TH1D("hEWKPhiStar_EffStatUp", "",nBinsPhiStar,PhiStarBins); hEWKPhiStar_EffStatUp->Sumw2();
  TH1D *hTopPhiStar_EffStatUp  = new TH1D("hTopPhiStar_EffStatUp", "",nBinsPhiStar,PhiStarBins); hTopPhiStar_EffStatUp->Sumw2();
  TH1D *hEWKPhiStar_EffStatDown  = new TH1D("hEWKPhiStar_EffStatDown", "",nBinsPhiStar,PhiStarBins); hEWKPhiStar_EffStatDown->Sumw2();
  TH1D *hTopPhiStar_EffStatDown  = new TH1D("hTopPhiStar_EffStatDown", "",nBinsPhiStar,PhiStarBins); hTopPhiStar_EffStatDown->Sumw2();
  TH1D *hEWKPhiStar_EffSigShape  = new TH1D("hEWKPhiStar_EffSigShape", "",nBinsPhiStar,PhiStarBins); hEWKPhiStar_EffSigShape->Sumw2();
  TH1D *hTopPhiStar_EffSigShape  = new TH1D("hTopPhiStar_EffSigShape", "",nBinsPhiStar,PhiStarBins); hTopPhiStar_EffSigShape->Sumw2();
  TH1D *hEWKPhiStar_EffBkgShape  = new TH1D("hEWKPhiStar_EffBkgShape", "",nBinsPhiStar,PhiStarBins); hEWKPhiStar_EffBkgShape->Sumw2();
  TH1D *hTopPhiStar_EffBkgShape  = new TH1D("hTopPhiStar_EffBkgShape", "",nBinsPhiStar,PhiStarBins); hTopPhiStar_EffBkgShape->Sumw2();

  
  TH1D *hDataZRap = new TH1D("hDataZRap","",24,0,2.4); hDataZRap->Sumw2();
  TH1D *hZeeZRap  = new TH1D("hZeeZRap", "",24,0,2.4); hZeeZRap->Sumw2();
  TH1D *hEWKZRap  = new TH1D("hEWKZRap", "",24,0,2.4); hEWKZRap->Sumw2();
  TH1D *hTopZRap  = new TH1D("hTopZRap", "",24,0,2.4); hTopZRap->Sumw2();
  TH1D *hMCZRap   = new TH1D("hMCZRap",  "",24,0,2.4); hMCZRap->Sumw2();
  
  
  TH1D *hZeeZRapUp  = new TH1D("hZeeZRapUp", "",24,0,2.4); hZeeZRapUp->Sumw2();
  TH1D *hZeeZRapDown  = new TH1D("hZeeZRapDown", "",24,0,2.4); hZeeZRapDown->Sumw2();

  TH1D *hEWKZRap_EffBin  = new TH1D("hEWKZRap_EffBin", "",24,0,2.4); hEWKZRap_EffBin->Sumw2();
  TH1D *hTopZRap_EffBin  = new TH1D("hTopZRap_EffBin", "",24,0,2.4); hTopZRap_EffBin->Sumw2();
  TH1D *hEWKZRap_EffStatUp  = new TH1D("hEWKZRap_EffStatUp", "",24,0,2.4); hEWKZRap_EffStatUp->Sumw2();
  TH1D *hTopZRap_EffStatUp  = new TH1D("hTopZRap_EffStatUp", "",24,0,2.4); hTopZRap_EffStatUp->Sumw2();
  TH1D *hEWKZRap_EffStatDown  = new TH1D("hEWKZRap_EffStatDown", "",24,0,2.4); hEWKZRap_EffStatDown->Sumw2();
  TH1D *hTopZRap_EffStatDown  = new TH1D("hTopZRap_EffStatDown", "",24,0,2.4); hTopZRap_EffStatDown->Sumw2();
  TH1D *hEWKZRap_EffSigShape  = new TH1D("hEWKZRap_EffSigShape", "",24,0,2.4); hEWKZRap_EffSigShape->Sumw2();
  TH1D *hTopZRap_EffSigShape  = new TH1D("hTopZRap_EffSigShape", "",24,0,2.4); hTopZRap_EffSigShape->Sumw2();
  TH1D *hEWKZRap_EffBkgShape  = new TH1D("hEWKZRap_EffBkgShape", "",24,0,2.4); hEWKZRap_EffBkgShape->Sumw2();
  TH1D *hTopZRap_EffBkgShape  = new TH1D("hTopZRap_EffBkgShape", "",24,0,2.4); hTopZRap_EffBkgShape->Sumw2();


  const int nBinsLep1Pt= sizeof(Lep1PtBins)/sizeof(double)-1;
  TH1D *hDataLep1Pt = new TH1D("hDataLep1Pt","",nBinsLep1Pt,Lep1PtBins); hDataLep1Pt->Sumw2();
  TH1D *hZeeLep1Pt  = new TH1D("hZeeLep1Pt", "",nBinsLep1Pt,Lep1PtBins); hZeeLep1Pt->Sumw2();
  TH1D *hEWKLep1Pt  = new TH1D("hEWKLep1Pt", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt->Sumw2();
  TH1D *hTopLep1Pt  = new TH1D("hTopLep1Pt", "",nBinsLep1Pt,Lep1PtBins); hTopLep1Pt->Sumw2();
  TH1D *hMCLep1Pt   = new TH1D("hMCLep1Pt",  "",nBinsLep1Pt,Lep1PtBins); hMCLep1Pt->Sumw2();

  TH1D *hEWKLep1Pt_EffBin  = new TH1D("hEWKLep1Pt_EffBin", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt_EffBin->Sumw2();
  TH1D *hTopLep1Pt_EffBin  = new TH1D("hTopLep1Pt_EffBin", "",nBinsLep1Pt,Lep1PtBins); hTopLep1Pt_EffBin->Sumw2();
  TH1D *hEWKLep1Pt_EffStatUp  = new TH1D("hEWKLep1Pt_EffStatUp", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt_EffStatUp->Sumw2();
  TH1D *hTopLep1Pt_EffStatUp  = new TH1D("hTopLep1Pt_EffStatUp", "",nBinsLep1Pt,Lep1PtBins); hTopLep1Pt_EffStatUp->Sumw2();
  TH1D *hEWKLep1Pt_EffStatDown  = new TH1D("hEWKLep1Pt_EffStatDown", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt_EffStatDown->Sumw2();
  TH1D *hTopLep1Pt_EffStatDown  = new TH1D("hTopLep1Pt_EffStatDown", "",nBinsLep1Pt,Lep1PtBins); hTopLep1Pt_EffStatDown->Sumw2();
  TH1D *hEWKLep1Pt_EffSigShape  = new TH1D("hEWKLep1Pt_EffSigShape", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt_EffSigShape->Sumw2();
  TH1D *hTopLep1Pt_EffSigShape  = new TH1D("hTopLep1Pt_EffSigShape", "",nBinsLep1Pt,Lep1PtBins); hTopLep1Pt_EffSigShape->Sumw2();
  TH1D *hEWKLep1Pt_EffBkgShape  = new TH1D("hEWKLep1Pt_EffBkgShape", "",nBinsLep1Pt,Lep1PtBins); hEWKLep1Pt_EffBkgShape->Sumw2();
  TH1D *hTopLep1Pt_EffBkgShape  = new TH1D("hTopLep1Pt_EffBkgShape", "",nBinsLep1Pt,Lep1PtBins); hTopLep1Pt_EffBkgShape->Sumw2();

  const int nBinsLep2Pt= sizeof(Lep2PtBins)/sizeof(double)-1;
  TH1D *hDataLep2Pt = new TH1D("hDataLep2Pt","",nBinsLep2Pt,Lep2PtBins); hDataLep2Pt->Sumw2();
  TH1D *hZeeLep2Pt  = new TH1D("hZeeLep2Pt", "",nBinsLep2Pt,Lep2PtBins); hZeeLep2Pt->Sumw2();
  TH1D *hEWKLep2Pt  = new TH1D("hEWKLep2Pt", "",nBinsLep2Pt,Lep2PtBins); hEWKLep2Pt->Sumw2();
  TH1D *hTopLep2Pt  = new TH1D("hTopLep2Pt", "",nBinsLep2Pt,Lep2PtBins); hTopLep2Pt->Sumw2();
  TH1D *hMCLep2Pt   = new TH1D("hMCLep2Pt",  "",nBinsLep2Pt,Lep2PtBins); hMCLep2Pt->Sumw2();

  TH1D *hEWKLep2Pt_EffBin  = new TH1D("hEWKLep2Pt_EffBin", "",nBinsLep2Pt,Lep2PtBins); hEWKLep2Pt_EffBin->Sumw2();
  TH1D *hTopLep2Pt_EffBin  = new TH1D("hTopLep2Pt_EffBin", "",nBinsLep2Pt,Lep2PtBins); hTopLep2Pt_EffBin->Sumw2();
  TH1D *hEWKLep2Pt_EffStatUp  = new TH1D("hEWKLep2Pt_EffStatUp", "",nBinsLep2Pt,Lep2PtBins); hEWKLep2Pt_EffStatUp->Sumw2();
  TH1D *hTopLep2Pt_EffStatUp  = new TH1D("hTopLep2Pt_EffStatUp", "",nBinsLep2Pt,Lep2PtBins); hTopLep2Pt_EffStatUp->Sumw2();
  TH1D *hEWKLep2Pt_EffStatDown  = new TH1D("hEWKLep2Pt_EffStatDown", "",nBinsLep2Pt,Lep2PtBins); hEWKLep2Pt_EffStatDown->Sumw2();
  TH1D *hTopLep2Pt_EffStatDown  = new TH1D("hTopLep2Pt_EffStatDown", "",nBinsLep2Pt,Lep2PtBins); hTopLep2Pt_EffStatDown->Sumw2();
  TH1D *hEWKLep2Pt_EffSigShape  = new TH1D("hEWKLep2Pt_EffSigShape", "",nBinsLep2Pt,Lep2PtBins); hEWKLep2Pt_EffSigShape->Sumw2();
  TH1D *hTopLep2Pt_EffSigShape  = new TH1D("hTopLep2Pt_EffSigShape", "",nBinsLep2Pt,Lep2PtBins); hTopLep2Pt_EffSigShape->Sumw2();
  TH1D *hEWKLep2Pt_EffBkgShape  = new TH1D("hEWKLep2Pt_EffBkgShape", "",nBinsLep2Pt,Lep2PtBins); hEWKLep2Pt_EffBkgShape->Sumw2();
  TH1D *hTopLep2Pt_EffBkgShape  = new TH1D("hTopLep2Pt_EffBkgShape", "",nBinsLep2Pt,Lep2PtBins); hTopLep2Pt_EffBkgShape->Sumw2();


  const int nBinsLepNegPt= sizeof(LepNegPtBins)/sizeof(double)-1;
  TH1D *hDataLepNegPt = new TH1D("hDataLepNegPt","",nBinsLepNegPt,LepNegPtBins); hDataLepNegPt->Sumw2();
  TH1D *hZeeLepNegPt  = new TH1D("hZeeLepNegPt", "",nBinsLepNegPt,LepNegPtBins); hZeeLepNegPt->Sumw2();
  TH1D *hEWKLepNegPt  = new TH1D("hEWKLepNegPt", "",nBinsLepNegPt,LepNegPtBins); hEWKLepNegPt->Sumw2();
  TH1D *hTopLepNegPt  = new TH1D("hTopLepNegPt", "",nBinsLepNegPt,LepNegPtBins); hTopLepNegPt->Sumw2();
  TH1D *hMCLepNegPt   = new TH1D("hMCLepNegPt",  "",nBinsLepNegPt,LepNegPtBins); hMCLepNegPt->Sumw2();

  TH1D *hEWKLepNegPt_EffBin  = new TH1D("hEWKLepNegPt_EffBin", "",nBinsLepNegPt,LepNegPtBins); hEWKLepNegPt_EffBin->Sumw2();
  TH1D *hTopLepNegPt_EffBin  = new TH1D("hTopLepNegPt_EffBin", "",nBinsLepNegPt,LepNegPtBins); hTopLepNegPt_EffBin->Sumw2();
  TH1D *hEWKLepNegPt_EffStatUp  = new TH1D("hEWKLepNegPt_EffStatUp", "",nBinsLepNegPt,LepNegPtBins); hEWKLepNegPt_EffStatUp->Sumw2();
  TH1D *hTopLepNegPt_EffStatUp  = new TH1D("hTopLepNegPt_EffStatUp", "",nBinsLepNegPt,LepNegPtBins); hTopLepNegPt_EffStatUp->Sumw2();
  TH1D *hEWKLepNegPt_EffStatDown  = new TH1D("hEWKLepNegPt_EffStatDown", "",nBinsLepNegPt,LepNegPtBins); hEWKLepNegPt_EffStatDown->Sumw2();
  TH1D *hTopLepNegPt_EffStatDown  = new TH1D("hTopLepNegPt_EffStatDown", "",nBinsLepNegPt,LepNegPtBins); hTopLepNegPt_EffStatDown->Sumw2();
  TH1D *hEWKLepNegPt_EffSigShape  = new TH1D("hEWKLepNegPt_EffSigShape", "",nBinsLepNegPt,LepNegPtBins); hEWKLepNegPt_EffSigShape->Sumw2();
  TH1D *hTopLepNegPt_EffSigShape  = new TH1D("hTopLepNegPt_EffSigShape", "",nBinsLepNegPt,LepNegPtBins); hTopLepNegPt_EffSigShape->Sumw2();
  TH1D *hEWKLepNegPt_EffBkgShape  = new TH1D("hEWKLepNegPt_EffBkgShape", "",nBinsLepNegPt,LepNegPtBins); hEWKLepNegPt_EffBkgShape->Sumw2();
  TH1D *hTopLepNegPt_EffBkgShape  = new TH1D("hTopLepNegPt_EffBkgShape", "",nBinsLepNegPt,LepNegPtBins); hTopLepNegPt_EffBkgShape->Sumw2();


  const int nBinsLepPosPt= sizeof(LepPosPtBins)/sizeof(double)-1;
  TH1D *hDataLepPosPt = new TH1D("hDataLepPosPt","",nBinsLepPosPt,LepPosPtBins); hDataLepPosPt->Sumw2();
  TH1D *hZeeLepPosPt  = new TH1D("hZeeLepPosPt", "",nBinsLepPosPt,LepPosPtBins); hZeeLepPosPt->Sumw2();
  TH1D *hEWKLepPosPt  = new TH1D("hEWKLepPosPt", "",nBinsLepPosPt,LepPosPtBins); hEWKLepPosPt->Sumw2();
  TH1D *hTopLepPosPt  = new TH1D("hTopLepPosPt", "",nBinsLepPosPt,LepPosPtBins); hTopLepPosPt->Sumw2();
  TH1D *hMCLepPosPt   = new TH1D("hMCLepPosPt",  "",nBinsLepPosPt,LepPosPtBins); hMCLepPosPt->Sumw2();

  TH1D *hEWKLepPosPt_EffBin  = new TH1D("hEWKLepPosPt_EffBin", "",nBinsLepPosPt,LepPosPtBins); hEWKLepPosPt_EffBin->Sumw2();
  TH1D *hTopLepPosPt_EffBin  = new TH1D("hTopLepPosPt_EffBin", "",nBinsLepPosPt,LepPosPtBins); hTopLepPosPt_EffBin->Sumw2();
  TH1D *hEWKLepPosPt_EffStatUp  = new TH1D("hEWKLepPosPt_EffStatUp", "",nBinsLepPosPt,LepPosPtBins); hEWKLepPosPt_EffStatUp->Sumw2();
  TH1D *hTopLepPosPt_EffStatUp  = new TH1D("hTopLepPosPt_EffStatUp", "",nBinsLepPosPt,LepPosPtBins); hTopLepPosPt_EffStatUp->Sumw2();
  TH1D *hEWKLepPosPt_EffStatDown  = new TH1D("hEWKLepPosPt_EffStatDown", "",nBinsLepPosPt,LepPosPtBins); hEWKLepPosPt_EffStatDown->Sumw2();
  TH1D *hTopLepPosPt_EffStatDown  = new TH1D("hTopLepPosPt_EffStatDown", "",nBinsLepPosPt,LepPosPtBins); hTopLepPosPt_EffStatDown->Sumw2();
  TH1D *hEWKLepPosPt_EffSigShape  = new TH1D("hEWKLepPosPt_EffSigShape", "",nBinsLepPosPt,LepPosPtBins); hEWKLepPosPt_EffSigShape->Sumw2();
  TH1D *hTopLepPosPt_EffSigShape  = new TH1D("hTopLepPosPt_EffSigShape", "",nBinsLepPosPt,LepPosPtBins); hTopLepPosPt_EffSigShape->Sumw2();
  TH1D *hEWKLepPosPt_EffBkgShape  = new TH1D("hEWKLepPosPt_EffBkgShape", "",nBinsLepPosPt,LepPosPtBins); hEWKLepPosPt_EffBkgShape->Sumw2();
  TH1D *hTopLepPosPt_EffBkgShape  = new TH1D("hTopLepPosPt_EffBkgShape", "",nBinsLepPosPt,LepPosPtBins); hTopLepPosPt_EffBkgShape->Sumw2();



  TH1D *hDataLep1Eta = new TH1D("hDataLep1Eta","",24,0,2.4); hDataLep1Eta->Sumw2();
  TH1D *hZeeLep1Eta  = new TH1D("hZeeLep1Eta", "",24,0,2.4); hZeeLep1Eta->Sumw2();
  TH1D *hEWKLep1Eta  = new TH1D("hEWKLep1Eta", "",24,0,2.4); hEWKLep1Eta->Sumw2();
  TH1D *hTopLep1Eta  = new TH1D("hTopLep1Eta", "",24,0,2.4); hTopLep1Eta->Sumw2();
  TH1D *hMCLep1Eta   = new TH1D("hMCLep1Eta",  "",24,0,2.4); hMCLep1Eta->Sumw2();

  TH1D *hEWKLep1Eta_EffBin  = new TH1D("hEWKLep1Eta_EffBin", "",24,0,2.4); hEWKLep1Eta_EffBin->Sumw2();
  TH1D *hTopLep1Eta_EffBin  = new TH1D("hTopLep1Eta_EffBin", "",24,0,2.4); hTopLep1Eta_EffBin->Sumw2();
  TH1D *hEWKLep1Eta_EffStatUp  = new TH1D("hEWKLep1Eta_EffStatUp", "",24,0,2.4); hEWKLep1Eta_EffStatUp->Sumw2();
  TH1D *hTopLep1Eta_EffStatUp  = new TH1D("hTopLep1Eta_EffStatUp", "",24,0,2.4); hTopLep1Eta_EffStatUp->Sumw2();
  TH1D *hEWKLep1Eta_EffStatDown  = new TH1D("hEWKLep1Eta_EffStatDown", "",24,0,2.4); hEWKLep1Eta_EffStatDown->Sumw2();
  TH1D *hTopLep1Eta_EffStatDown  = new TH1D("hTopLep1Eta_EffStatDown", "",24,0,2.4); hTopLep1Eta_EffStatDown->Sumw2();
  TH1D *hEWKLep1Eta_EffSigShape  = new TH1D("hEWKLep1Eta_EffSigShape", "",24,0,2.4); hEWKLep1Eta_EffSigShape->Sumw2();
  TH1D *hTopLep1Eta_EffSigShape  = new TH1D("hTopLep1Eta_EffSigShape", "",24,0,2.4); hTopLep1Eta_EffSigShape->Sumw2();
  TH1D *hEWKLep1Eta_EffBkgShape  = new TH1D("hEWKLep1Eta_EffBkgShape", "",24,0,2.4); hEWKLep1Eta_EffBkgShape->Sumw2();
  TH1D *hTopLep1Eta_EffBkgShape  = new TH1D("hTopLep1Eta_EffBkgShape", "",24,0,2.4); hTopLep1Eta_EffBkgShape->Sumw2();


  TH1D *hDataLep2Eta = new TH1D("hDataLep2Eta","",24,0,2.4); hDataLep2Eta->Sumw2();
  TH1D *hZeeLep2Eta  = new TH1D("hZeeLep2Eta", "",24,0,2.4); hZeeLep2Eta->Sumw2();
  TH1D *hEWKLep2Eta  = new TH1D("hEWKLep2Eta", "",24,0,2.4); hEWKLep2Eta->Sumw2();
  TH1D *hTopLep2Eta  = new TH1D("hTopLep2Eta", "",24,0,2.4); hTopLep2Eta->Sumw2();
  TH1D *hMCLep2Eta   = new TH1D("hMCLep2Eta",  "",24,0,2.4); hMCLep2Eta->Sumw2();

  TH1D *hEWKLep2Eta_EffBin  = new TH1D("hEWKLep2Eta_EffBin", "",24,0,2.4); hEWKLep2Eta_EffBin->Sumw2();
  TH1D *hTopLep2Eta_EffBin  = new TH1D("hTopLep2Eta_EffBin", "",24,0,2.4); hTopLep2Eta_EffBin->Sumw2();
  TH1D *hEWKLep2Eta_EffStatUp  = new TH1D("hEWKLep2Eta_EffStatUp", "",24,0,2.4); hEWKLep2Eta_EffStatUp->Sumw2();
  TH1D *hTopLep2Eta_EffStatUp  = new TH1D("hTopLep2Eta_EffStatUp", "",24,0,2.4); hTopLep2Eta_EffStatUp->Sumw2();
  TH1D *hEWKLep2Eta_EffStatDown  = new TH1D("hEWKLep2Eta_EffStatDown", "",24,0,2.4); hEWKLep2Eta_EffStatDown->Sumw2();
  TH1D *hTopLep2Eta_EffStatDown  = new TH1D("hTopLep2Eta_EffStatDown", "",24,0,2.4); hTopLep2Eta_EffStatDown->Sumw2();
  TH1D *hEWKLep2Eta_EffSigShape  = new TH1D("hEWKLep2Eta_EffSigShape", "",24,0,2.4); hEWKLep2Eta_EffSigShape->Sumw2();
  TH1D *hTopLep2Eta_EffSigShape  = new TH1D("hTopLep2Eta_EffSigShape", "",24,0,2.4); hTopLep2Eta_EffSigShape->Sumw2();
  TH1D *hEWKLep2Eta_EffBkgShape  = new TH1D("hEWKLep2Eta_EffBkgShape", "",24,0,2.4); hEWKLep2Eta_EffBkgShape->Sumw2();
  TH1D *hTopLep2Eta_EffBkgShape  = new TH1D("hTopLep2Eta_EffBkgShape", "",24,0,2.4); hTopLep2Eta_EffBkgShape->Sumw2();
     
     
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
  Float_t prefireWeight, prefireUp, prefireDown;
  Int_t   q1, q2;
  TLorentzVector *lep1=0, *lep2=0;
  TLorentzVector *lep1_raw=0, *lep2_raw=0;
  TLorentzVector *dilep=0, *dilepSC = 0;
  TLorentzVector *sc1=0, *sc2=0;
  
  TFile *infile=0;
  TTree *intree=0;
  Float_t r91=0; 
  Float_t r92=0;
  Float_t random;

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
    intree -> SetBranchStatus("dilepSC",1);
    intree -> SetBranchStatus("r91",1);
    intree -> SetBranchStatus("r92",1);
    intree -> SetBranchStatus("random",1);

    intree->SetBranchAddress("runNum",     &runNum);      // event run number
    intree->SetBranchAddress("lumiSec",    &lumiSec);     // event lumi section
    intree->SetBranchAddress("evtNum",     &evtNum);      // event number
    intree->SetBranchAddress("category",   &category);    // dilepton category
    intree->SetBranchAddress("npv",        &npv);	  // number of primary vertices
    intree->SetBranchAddress("npu",        &npu);	  // number of in-time PU events (MC)
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
    intree->SetBranchAddress("dilepSC",   &dilepSC);      // sc2 4-vector
    intree->SetBranchAddress("r91",       &r91);        // sc2 4-vector
    intree->SetBranchAddress("r92",       &r92);        // sc2 4-vector
    intree->SetBranchAddress("random",       &random);        // sc2 4-vector
      
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
    // for(UInt_t ientry=0; ientry<(UInt_t)(0.1*intree->GetEntries()); ientry++) {
        if(ientry%100000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " percent done with this file." << endl;
      intree->GetEntry(ientry);
      // std::cout << "r92 " << r92 << std::endl;
      if(fabs(lep1->Eta()) > ETA_CUT)   continue;      
      if(fabs(lep2->Eta()) > ETA_CUT)   continue;
      if(q1*q2>0) continue;
      if(lep1->Pt()        < PT_CUT)    continue;
      if(lep2->Pt()       < PT_CUT)    continue;
      
      
      hCompareElePtEcalE->Fill(lep1->Pt(),sc1->Pt());
      // hCompareElePtEcalE->Fill(lep2->Pt(),sc2->Pt());
      
      // if(fabs(lep1->Eta()) < ETA_CUT || fabs(lep2->Eta()) < ETA_CUT) continue;
      float ETA_LOW = 0; float ETA_HIGH = 1.0;
      
      float mass = 0;
      float pt = 0;
      float rapidity = 0;
      float phiacop=0;
      float costhetastar=0;
      float phistar=0;
     
      Double_t weight=1, weightUp=1, weightDown=1;
      if(typev[ifile]!=eData) {
	    weight *= scale1fb*prefireWeight*lumi/totalNorm;
	    weightUp *= scale1fb*prefireUp*lumi/totalNorm;
	    weightDown *= scale1fb*prefireDown*lumi/totalNorm;
	    // weight *= scale1fb*lumi/totalNorm;
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



        mass=(l1+l2).M();
        pt =(l1+l2).Pt();
        rapidity = (l1+l2).Rapidity();

        phiacop=TMath::Pi()-fabs(l1.DeltaPhi(l2));
        if(lq1<0) costhetastar=tanh(float((l1.Rapidity()-l2.Rapidity())/2));
        else costhetastar=tanh(float((l2.Rapidity()-l1.Rapidity())/2));
        phistar=tan(phiacop/2)*sqrt(1-pow(costhetastar,2));
        
        if(mass        < MASS_LOW)  continue;
        if(mass        > MASS_HIGH) continue;
        if(l1.Pt()        < PT_CUT)    continue;
        if(l2.Pt()        < PT_CUT)    continue;
        hDataEG->Fill(dilepSC->M());

        hData->Fill(mass); 
        hDataNPV->Fill(npv);
        hDataZPt->Fill(pt); 
        hDataPhiStar->Fill(phistar); 
        hDataLep1Pt->Fill(l1.Pt()); 
        hDataLep2Pt->Fill(l2.Pt()); 
        if(lq1<0) {
          hDataLepNegPt->Fill(l1.Pt()); 
          hDataLepPosPt->Fill(l2.Pt());
        } else  {
          hDataLepNegPt->Fill(l2.Pt()); 
          hDataLepPosPt->Fill(l1.Pt());
        }
        hDataLep1Eta->Fill(fabs(l1.Eta())); 
        hDataLep2Eta->Fill(fabs(l2.Eta())); 
        hDataZRap->Fill(fabs(rapidity));
        
        yield++;
	  
      } else {

        Double_t lp1 = lep1->Pt();
        Double_t lp2 = lep2->Pt();
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

        double tagSmear1 = ec.smearingSigma(runNumber, sc1->Pt(), fabs(sc1->Eta()), r91, 12, 0., 0.);
        double tagSmear2 = ec.smearingSigma(runNumber, sc2->Pt(), fabs(sc2->Eta()), r92, 12, 0., 0.);
        hGausRandNtuple->Fill(random);

        
        double mll=(l1+l2).M();
        Double_t effdata, effmc;
        Double_t corr=1;
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
        if(lp1        < PT_CUT)    continue;
        if(lp2        < PT_CUT)    continue;
        
        double var = 0;
    
        corr = effs.fullEfficiencies(&l1,q1,&l2,q2);
        mass = (l1+l2).M();
        pt = (l1+l2).Pt();
        rapidity = (l1+l2).Rapidity();

        phiacop=TMath::Pi()-fabs(l1.DeltaPhi(l2));
        if(lq1<0) costhetastar=tanh(float((l1.Rapidity()-l2.Rapidity())/2));
        else costhetastar=tanh(float((l2.Rapidity()-l1.Rapidity())/2));
        phistar=tan(phiacop/2)*sqrt(1-pow(costhetastar,2));
	  
      if(typev[ifile]==eZee)  { 
        yield_zee += weight*corr;
        yield_zee_unc += weight*weight*corr*corr;
        hZee->Fill(mass,weight*corr); 
        hZeeEG->Fill(dilepSC->M(),weight*corr); 
        hMC->Fill(mass,weight*corr);
        hZeeNPV->Fill(npv,weight*corr); 
        hMCNPV->Fill(npv,weight*corr);
        hZeeZPt->Fill(pt,weight*corr); 
        hMCZPt->Fill(pt,weight*corr);
        hZeePhiStar->Fill(phistar,weight*corr); 
        hMCPhiStar->Fill(phistar,weight*corr);
        hZeeZRap->Fill(fabs(rapidity),weight*corr); 
        hZeeZRapUp->Fill(fabs(rapidity),weightUp*corr); 
        hZeeZRapDown->Fill(fabs(rapidity),weightDown*corr); 
        hMCZRap->Fill(fabs(rapidity),weight*corr);
        hZeeLep1Pt->Fill(l1.Pt(),weight*corr); 
        hMCLep1Pt->Fill(l1.Pt(),weight*corr);
        hZeeLep2Pt->Fill(l2.Pt(),weight*corr); 
        hMCLep2Pt->Fill(l2.Pt(),weight*corr);
        if(lq1<0)	{
          hZeeLepNegPt->Fill(l1.Pt(),weight*corr); 
          hMCLepNegPt->Fill(l1.Pt(),weight*corr);
          hZeeLepPosPt->Fill(l2.Pt(),weight*corr); 
          hMCLepPosPt->Fill(l2.Pt(),weight*corr);
        } else {
          hZeeLepNegPt->Fill(l2.Pt(),weight*corr); 
          hMCLepNegPt->Fill(l2.Pt(),weight*corr);
          hZeeLepPosPt->Fill(l1.Pt(),weight*corr); 
          hMCLepPosPt->Fill(l1.Pt(),weight*corr);
        }
        hZeeLep1Eta->Fill(fabs(l1.Eta()),weight*corr); 
        hMCLep1Eta->Fill(fabs(l1.Eta()),weight*corr);
        hZeeLep2Eta->Fill(fabs(l2.Eta()),weight*corr); 
        hMCLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	    }
      if(typev[ifile]==eEWK) { 
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
	      
	      hEWKPhiStar->Fill(phistar,weight*corr);
	      hEWKPhiStar_EffBin->Fill(phistar,weight*corr2Bin); 
	      hEWKPhiStar_EffStatUp->Fill(phistar,weight*corrUp);
	      hEWKPhiStar_EffStatDown->Fill(phistar,weight*corrDown);
	      hEWKPhiStar_EffSigShape->Fill(phistar,weight*corrSigShape);
	      hEWKPhiStar_EffBkgShape->Fill(phistar,weight*corrBkgShape);
	      hMCPhiStar->Fill(phistar,weight*corr);
	      
	      hEWKZRap->Fill(fabs(rapidity),weight*corr); 
	      hEWKZRap_EffBin->Fill(fabs(rapidity),weight*corr2Bin);
	      hEWKZRap_EffStatUp->Fill(fabs(rapidity),weight*corrUp);
	      hEWKZRap_EffStatDown->Fill(fabs(rapidity),weight*corrDown);
	      hEWKZRap_EffSigShape->Fill(fabs(rapidity),weight*corrSigShape);
	      hEWKZRap_EffBkgShape->Fill(fabs(rapidity),weight*corrBkgShape);
	      hMCZRap->Fill(fabs(rapidity),weight*corr);
	      
	      hEWKLep1Pt->Fill(l1.Pt(),weight*corr);
	      hEWKLep1Pt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
	      hEWKLep1Pt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
	      hEWKLep1Pt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
	      hEWKLep1Pt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
	      hEWKLep1Pt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
	      hMCLep1Pt->Fill(l1.Pt(),weight*corr);
	      
	      hEWKLep2Pt->Fill(l2.Pt(),weight*corr);
	      hEWKLep2Pt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
	      hEWKLep2Pt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
	      hEWKLep2Pt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
	      hEWKLep2Pt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
	      hEWKLep2Pt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
	      hMCLep2Pt->Fill(l2.Pt(),weight*corr);
	      
	      if(lq1<0) {
          hEWKLepNegPt->Fill(l1.Pt(),weight*corr);
          hEWKLepNegPt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
          hEWKLepNegPt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
          hEWKLepNegPt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
          hEWKLepNegPt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
          hEWKLepNegPt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
          hMCLepNegPt->Fill(l1.Pt(),weight*corr);
          
          hEWKLepPosPt->Fill(l2.Pt(),weight*corr);
          hEWKLepPosPt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
          hEWKLepPosPt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
          hEWKLepPosPt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
          hEWKLepPosPt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
          hEWKLepPosPt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
          hMCLepPosPt->Fill(l2.Pt(),weight*corr);
        } else {
          hEWKLepNegPt->Fill(l2.Pt(),weight*corr);
          hEWKLepNegPt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
          hEWKLepNegPt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
          hEWKLepNegPt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
          hEWKLepNegPt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
          hEWKLepNegPt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
          hMCLepNegPt->Fill(l2.Pt(),weight*corr);
          
          hEWKLepPosPt->Fill(l1.Pt(),weight*corr);
          hEWKLepPosPt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
          hEWKLepPosPt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
          hEWKLepPosPt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
          hEWKLepPosPt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
          hEWKLepPosPt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
          hMCLepPosPt->Fill(l1.Pt(),weight*corr);
        }
	      
	      hEWKLep1Eta->Fill(fabs(l1.Eta()),weight*corr); 
	      hEWKLep1Eta_EffBin->Fill(fabs(l1.Eta()),weight*corr2Bin);
	      hEWKLep1Eta_EffStatUp->Fill(fabs(l1.Eta()),weight*corrUp);
	      hEWKLep1Eta_EffStatDown->Fill(fabs(l1.Eta()),weight*corrDown);
	      hEWKLep1Eta_EffSigShape->Fill(fabs(l1.Eta()),weight*corrSigShape);
	      hEWKLep1Eta_EffBkgShape->Fill(fabs(l1.Eta()),weight*corrBkgShape);
	      hMCLep1Eta->Fill(fabs(l1.Eta()),weight*corr);
	      
	      hEWKLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	      hEWKLep2Eta_EffBin->Fill(fabs(l2.Eta()),weight*corr2Bin);
	      hEWKLep2Eta_EffStatUp->Fill(fabs(l2.Eta()),weight*corrUp);
	      hEWKLep2Eta_EffStatDown->Fill(fabs(l2.Eta()),weight*corrDown);
	      hEWKLep2Eta_EffSigShape->Fill(fabs(l2.Eta()),weight*corrSigShape);
	      hEWKLep2Eta_EffBkgShape->Fill(fabs(l2.Eta()),weight*corrBkgShape);
	      hMCLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	    }
      if(typev[ifile]==eTop) {
	      yield_top += weight*corr;
	      yield_top_unc += weight*weight*corr*corr;
	      hTop->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);

	      hTopNPV->Fill(npv,weight*corr); 
	      hMCNPV->Fill(npv,weight*corr);

	      hTopZPt->Fill(pt,weight*corr); 
	      hTopZPt_EffBin->Fill(pt,weight*corr2Bin); 
	      hTopZPt_EffStatUp->Fill(pt,weight*corrUp);
	      hTopZPt_EffStatDown->Fill(pt,weight*corrDown);
	      hTopZPt_EffSigShape->Fill(pt,weight*corrSigShape); 
	      hTopZPt_EffBkgShape->Fill(pt,weight*corrBkgShape);
	      hMCZPt->Fill(pt,weight*corr);

	      hTopPhiStar->Fill(phistar,weight*corr);
	      hTopPhiStar_EffBin->Fill(phistar,weight*corr2Bin); 
	      hTopPhiStar_EffStatUp->Fill(phistar,weight*corrUp);
	      hTopPhiStar_EffStatDown->Fill(phistar,weight*corrDown);
	      hTopPhiStar_EffSigShape->Fill(phistar,weight*corrSigShape); 
	      hTopPhiStar_EffBkgShape->Fill(phistar,weight*corrBkgShape); 
	      hMCPhiStar->Fill(phistar,weight*corr);

	      hTopZRap->Fill(fabs(rapidity),weight*corr); 
	      hTopZRap_EffBin->Fill(fabs(rapidity),weight*corr2Bin);
	      hTopZRap_EffStatUp->Fill(fabs(rapidity),weight*corrUp);
	      hTopZRap_EffStatDown->Fill(fabs(rapidity),weight*corrDown);
	      hTopZRap_EffSigShape->Fill(fabs(rapidity),weight*corrSigShape); 
	      hTopZRap_EffBkgShape->Fill(fabs(rapidity),weight*corrBkgShape); 
	      hMCZRap->Fill(fabs(rapidity),weight*corr);
	      
	      hTopLep1Pt->Fill(l1.Pt(),weight*corr);
	      hTopLep1Pt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
	      hTopLep1Pt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
	      hTopLep1Pt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
	      hTopLep1Pt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
	      hTopLep1Pt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
	      hMCLep1Pt->Fill(l1.Pt(),weight*corr);
	      
	      hTopLep2Pt->Fill(l2.Pt(),weight*corr);
	      hTopLep2Pt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
	      hTopLep2Pt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
	      hTopLep2Pt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
	      hTopLep2Pt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
	      hTopLep2Pt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
	      hMCLep2Pt->Fill(l2.Pt(),weight*corr);

	      if(lq1<0) {
          hTopLepNegPt->Fill(l1.Pt(),weight*corr);
          hTopLepNegPt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
          hTopLepNegPt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
          hTopLepNegPt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
          hTopLepNegPt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
          hTopLepNegPt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
          hMCLepNegPt->Fill(l1.Pt(),weight*corr);
          
          hTopLepPosPt->Fill(l2.Pt(),weight*corr);
          hTopLepPosPt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
          hTopLepPosPt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
          hTopLepPosPt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
          hTopLepPosPt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
          hTopLepPosPt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
          hMCLepPosPt->Fill(l2.Pt(),weight*corr);
        } else  {
          hTopLepNegPt->Fill(l2.Pt(),weight*corr);
          hTopLepNegPt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
          hTopLepNegPt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
          hTopLepNegPt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
          hTopLepNegPt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
          hTopLepNegPt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
          hMCLepNegPt->Fill(l2.Pt(),weight*corr);
          
          hTopLepPosPt->Fill(l1.Pt(),weight*corr);
          hTopLepPosPt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
          hTopLepPosPt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
          hTopLepPosPt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
          hTopLepPosPt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
          hTopLepPosPt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
          hMCLepPosPt->Fill(l1.Pt(),weight*corr);
        }

	      hTopLep1Eta->Fill(fabs(l1.Eta()),weight*corr);
	      hTopLep1Eta_EffBin->Fill(fabs(l1.Eta()),weight*corr2Bin);
	      hTopLep1Eta_EffStatUp->Fill(fabs(l1.Eta()),weight*corrUp);
	      hTopLep1Eta_EffStatDown->Fill(fabs(l1.Eta()),weight*corrDown);
	      hTopLep1Eta_EffSigShape->Fill(fabs(l1.Eta()),weight*corrSigShape);
	      hTopLep1Eta_EffBkgShape->Fill(fabs(l1.Eta()),weight*corrBkgShape);
	      hMCLep1Eta->Fill(fabs(l1.Eta()),weight*corr);

	      hTopLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	      hTopLep2Eta_EffBin->Fill(fabs(l2.Eta()),weight*corr2Bin);
	      hTopLep2Eta_EffStatUp->Fill(fabs(l2.Eta()),weight*corrUp);
	      hTopLep2Eta_EffStatDown->Fill(fabs(l2.Eta()),weight*corrDown);
	      hTopLep2Eta_EffSigShape->Fill(fabs(l2.Eta()),weight*corrSigShape);
	      hTopLep2Eta_EffBkgShape->Fill(fabs(l2.Eta()),weight*corrBkgShape);
	      hMCLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	    }
      }
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
  hTopZPt->Write();
  hTopZPt_EffBin->Write();
  hTopZPt_EffStatUp->Write();
  hTopZPt_EffStatDown->Write();
  hTopZPt_EffSigShape->Write();
  hTopZPt_EffBkgShape->Write();

  hDataPhiStar->Write();
  hEWKPhiStar->Write();
  hEWKPhiStar_EffBin->Write();
  hEWKPhiStar_EffStatUp->Write();
  hEWKPhiStar_EffStatDown->Write();
  hEWKPhiStar_EffSigShape->Write();
  hEWKPhiStar_EffBkgShape->Write();
  hTopPhiStar->Write();
  hTopPhiStar_EffBin->Write();
  hTopPhiStar_EffStatUp->Write();
  hTopPhiStar_EffStatDown->Write();
  hTopPhiStar_EffSigShape->Write();
  hTopPhiStar_EffBkgShape->Write();

 
  hDataZRap->Write();
  hEWKZRap->Write();
  hEWKZRap_EffBin->Write();
  hEWKZRap_EffStatUp->Write();
  hEWKZRap_EffStatDown->Write();
  hEWKZRap_EffSigShape->Write();
  hEWKZRap_EffBkgShape->Write();
  hTopZRap->Write();
  hTopZRap_EffBin->Write();
  hTopZRap_EffStatUp->Write();
  hTopZRap_EffStatDown->Write();
  hTopZRap_EffSigShape->Write();
  hTopZRap_EffBkgShape->Write();

  hDataLep1Pt->Write();
  hEWKLep1Pt->Write();
  hEWKLep1Pt_EffBin->Write();
  hEWKLep1Pt_EffStatUp->Write();
  hEWKLep1Pt_EffStatDown->Write();
  hEWKLep1Pt_EffSigShape->Write();
  hEWKLep1Pt_EffBkgShape->Write();
  hTopLep1Pt->Write();
  hTopLep1Pt_EffBin->Write();
  hTopLep1Pt_EffStatUp->Write();
  hTopLep1Pt_EffStatDown->Write();
  hTopLep1Pt_EffSigShape->Write();
  hTopLep1Pt_EffBkgShape->Write();

  hDataLep2Pt->Write();
  hEWKLep2Pt->Write();
  hEWKLep2Pt_EffBin->Write();
  hEWKLep2Pt_EffStatUp->Write();
  hEWKLep2Pt_EffStatDown->Write();
  hEWKLep2Pt_EffSigShape->Write();
  hEWKLep2Pt_EffBkgShape->Write();
  hTopLep2Pt->Write();
  hTopLep2Pt_EffBin->Write();
  hTopLep2Pt_EffStatUp->Write();
  hTopLep2Pt_EffStatDown->Write();
  hTopLep2Pt_EffSigShape->Write();
  hTopLep2Pt_EffBkgShape->Write();

  hDataLepNegPt->Write();
  hEWKLepNegPt->Write();
  hEWKLepNegPt_EffBin->Write();
  hEWKLepNegPt_EffStatUp->Write();
  hEWKLepNegPt_EffStatDown->Write();
  hEWKLepNegPt_EffSigShape->Write();
  hEWKLepNegPt_EffBkgShape->Write();
  hTopLepNegPt->Write();
  hTopLepNegPt_EffBin->Write();
  hTopLepNegPt_EffStatUp->Write();
  hTopLepNegPt_EffStatDown->Write();
  hTopLepNegPt_EffSigShape->Write();
  hTopLepNegPt_EffBkgShape->Write();

  hDataLepPosPt->Write();
  hEWKLepPosPt->Write();
  hEWKLepPosPt_EffBin->Write();
  hEWKLepPosPt_EffStatUp->Write();
  hEWKLepPosPt_EffStatDown->Write();
  hEWKLepPosPt_EffSigShape->Write();
  hEWKLepPosPt_EffBkgShape->Write();
  hTopLepPosPt->Write();
  hTopLepPosPt_EffBin->Write();
  hTopLepPosPt_EffStatUp->Write();
  hTopLepPosPt_EffStatDown->Write();
  hTopLepPosPt_EffSigShape->Write();
  hTopLepPosPt_EffBkgShape->Write();

  hDataLep1Eta->Write();
  hEWKLep1Eta->Write();
  hEWKLep1Eta_EffBin->Write();
  hEWKLep1Eta_EffStatUp->Write();
  hEWKLep1Eta_EffStatDown->Write();
  hEWKLep1Eta_EffSigShape->Write();
  hEWKLep1Eta_EffBkgShape->Write();
  hTopLep1Eta->Write();
  hTopLep1Eta_EffBin->Write();
  hTopLep1Eta_EffStatUp->Write();
  hTopLep1Eta_EffStatDown->Write();
  hTopLep1Eta_EffSigShape->Write();
  hTopLep1Eta_EffBkgShape->Write();

  hDataLep2Eta->Write();
  hEWKLep2Eta->Write();
  hEWKLep2Eta_EffBin->Write();
  hEWKLep2Eta_EffStatUp->Write();
  hEWKLep2Eta_EffStatDown->Write();
  hEWKLep2Eta_EffSigShape->Write();
  hEWKLep2Eta_EffBkgShape->Write();
  hTopLep2Eta->Write();
  hTopLep2Eta_EffBin->Write();
  hTopLep2Eta_EffStatUp->Write();
  hTopLep2Eta_EffStatDown->Write();
  hTopLep2Eta_EffSigShape->Write();
  hTopLep2Eta_EffBkgShape->Write();

  outFile->Write();
  outFile->Close();
  

 std::cout << "egamma ntuples size, MC: " << hZeeEG->Integral() << "  data: " << hDataEG->Integral() << std::endl;
 std::cout << "ours   ntuples size, MC: " << hZee->Integral()   << "  data: " << hData->Integral()   << std::endl;

  double MCscale=hData->Integral()/hMC->Integral();
  double scaleDat=hData->Integral()/hDataEG->Integral();
  // double scaleMC=hZee->Integral()/hZeeEG->Integral();
  // double scaleme=hData->Integral()/hZee->Integral();
  double scaleyou=hDataEG->Integral()/hZeeEG->Integral();

  hDataEG->Scale(scaleDat);
  hZeeEG->Scale(MCscale);
  double GausScale = hGausRandHere->Integral()/hGausRandNtuple->Integral();
  hGausRandNtuple->Scale(GausScale);

  if(normToData)
    {

      cout<<"Normalized to data: "<<MCscale<<endl;
      
      hZee->Scale(MCscale);
      hMC->Scale(MCscale);
      hEWK->Scale(MCscale);
      hTop->Scale(MCscale);
      hZeeNPV->Scale(MCscale);
      hMCNPV->Scale(MCscale);
      hEWKNPV->Scale(MCscale);
      hTopNPV->Scale(MCscale);
      hZeeZPt->Scale(MCscale);
      hMCZPt->Scale(MCscale);
      hEWKZPt->Scale(MCscale);
      hTopZPt->Scale(MCscale);
      hZeePhiStar->Scale(MCscale);
      hMCPhiStar->Scale(MCscale);
      hEWKPhiStar->Scale(MCscale);
      hTopPhiStar->Scale(MCscale);
      hZeeZRap->Scale(MCscale);
      hZeeZRapUp->Scale(MCscale);
      hZeeZRapDown->Scale(MCscale);
      hMCZRap->Scale(MCscale);
      hEWKZRap->Scale(MCscale);
      hTopZRap->Scale(MCscale);
      hZeeLep1Pt->Scale(MCscale);
      hMCLep1Pt->Scale(MCscale);
      hEWKLep1Pt->Scale(MCscale);
      hTopLep1Pt->Scale(MCscale);
      hZeeLep2Pt->Scale(MCscale);
      hMCLep2Pt->Scale(MCscale);
      hEWKLep2Pt->Scale(MCscale);
      hTopLep2Pt->Scale(MCscale);
      hZeeLepNegPt->Scale(MCscale);
      hMCLepNegPt->Scale(MCscale);
      hEWKLepNegPt->Scale(MCscale);
      hTopLepNegPt->Scale(MCscale);
      hZeeLepPosPt->Scale(MCscale);
      hMCLepPosPt->Scale(MCscale);
      hEWKLepPosPt->Scale(MCscale);
      hTopLepPosPt->Scale(MCscale);
      hZeeLep1Eta->Scale(MCscale);
      hMCLep1Eta->Scale(MCscale);
      hEWKLep1Eta->Scale(MCscale);
      hTopLep1Eta->Scale(MCscale);
      hZeeLep2Eta->Scale(MCscale);
      hMCLep2Eta->Scale(MCscale);
      hEWKLep2Eta->Scale(MCscale);
      hTopLep2Eta->Scale(MCscale);
    }

  std::cout << hData->Integral() << std::endl;
  std::cout << hMC->Integral() << std::endl;
  std::cout << hZee->Integral() << std::endl;
  std::cout << hEWK->Integral() << std::endl;
  std::cout << hTop->Integral() << std::endl;

  for(int j=0;j!=hDataZPt->GetNbinsX();++j)
    {
      hDataZPt->SetBinContent(j+1,hDataZPt->GetBinContent(j+1)/hDataZPt->GetBinWidth(j+1));
      hDataZPt->SetBinError(j+1,hDataZPt->GetBinError(j+1)/hDataZPt->GetBinWidth(j+1));
      hMCZPt->SetBinContent(j+1,hMCZPt->GetBinContent(j+1)/hMCZPt->GetBinWidth(j+1));
      hMCZPt->SetBinError(j+1,hMCZPt->GetBinError(j+1)/hMCZPt->GetBinWidth(j+1));
      hZeeZPt->SetBinContent(j+1,hZeeZPt->GetBinContent(j+1)/hZeeZPt->GetBinWidth(j+1));
      hZeeZPt->SetBinError(j+1,hZeeZPt->GetBinError(j+1)/hZeeZPt->GetBinWidth(j+1));
      hEWKZPt->SetBinContent(j+1,hEWKZPt->GetBinContent(j+1)/hEWKZPt->GetBinWidth(j+1));
      hEWKZPt->SetBinError(j+1,hEWKZPt->GetBinError(j+1)/hEWKZPt->GetBinWidth(j+1));
      hTopZPt->SetBinContent(j+1,hTopZPt->GetBinContent(j+1)/hTopZPt->GetBinWidth(j+1));
      hTopZPt->SetBinError(j+1,hTopZPt->GetBinError(j+1)/hTopZPt->GetBinWidth(j+1));
    }

  for(int j=0;j!=hDataPhiStar->GetNbinsX();++j)
    {
      hDataPhiStar->SetBinContent(j+1,hDataPhiStar->GetBinContent(j+1)/hDataPhiStar->GetBinWidth(j+1));
      hDataPhiStar->SetBinError(j+1,hDataPhiStar->GetBinError(j+1)/hDataPhiStar->GetBinWidth(j+1));
      hMCPhiStar->SetBinContent(j+1,hMCPhiStar->GetBinContent(j+1)/hMCPhiStar->GetBinWidth(j+1));
      hMCPhiStar->SetBinError(j+1,hMCPhiStar->GetBinError(j+1)/hMCPhiStar->GetBinWidth(j+1));
      hZeePhiStar->SetBinContent(j+1,hZeePhiStar->GetBinContent(j+1)/hZeePhiStar->GetBinWidth(j+1));
      hZeePhiStar->SetBinError(j+1,hZeePhiStar->GetBinError(j+1)/hZeePhiStar->GetBinWidth(j+1));
      hEWKPhiStar->SetBinContent(j+1,hEWKPhiStar->GetBinContent(j+1)/hEWKPhiStar->GetBinWidth(j+1));
      hEWKPhiStar->SetBinError(j+1,hEWKPhiStar->GetBinError(j+1)/hEWKPhiStar->GetBinWidth(j+1));
      hTopPhiStar->SetBinContent(j+1,hTopPhiStar->GetBinContent(j+1)/hTopPhiStar->GetBinWidth(j+1));
      hTopPhiStar->SetBinError(j+1,hTopPhiStar->GetBinError(j+1)/hTopPhiStar->GetBinWidth(j+1));
    }

  for(int j=0;j!=hDataLep1Pt->GetNbinsX();++j)
    {
      hDataLep1Pt->SetBinContent(j+1,hDataLep1Pt->GetBinContent(j+1)/hDataLep1Pt->GetBinWidth(j+1));
      hDataLep1Pt->SetBinError(j+1,hDataLep1Pt->GetBinError(j+1)/hDataLep1Pt->GetBinWidth(j+1));
      hMCLep1Pt->SetBinContent(j+1,hMCLep1Pt->GetBinContent(j+1)/hMCLep1Pt->GetBinWidth(j+1));
      hMCLep1Pt->SetBinError(j+1,hMCLep1Pt->GetBinError(j+1)/hMCLep1Pt->GetBinWidth(j+1));
      hZeeLep1Pt->SetBinContent(j+1,hZeeLep1Pt->GetBinContent(j+1)/hZeeLep1Pt->GetBinWidth(j+1));
      hZeeLep1Pt->SetBinError(j+1,hZeeLep1Pt->GetBinError(j+1)/hZeeLep1Pt->GetBinWidth(j+1));
      hEWKLep1Pt->SetBinContent(j+1,hEWKLep1Pt->GetBinContent(j+1)/hEWKLep1Pt->GetBinWidth(j+1));
      hEWKLep1Pt->SetBinError(j+1,hEWKLep1Pt->GetBinError(j+1)/hEWKLep1Pt->GetBinWidth(j+1));
      hTopLep1Pt->SetBinContent(j+1,hTopLep1Pt->GetBinContent(j+1)/hTopLep1Pt->GetBinWidth(j+1));
      hTopLep1Pt->SetBinError(j+1,hTopLep1Pt->GetBinError(j+1)/hTopLep1Pt->GetBinWidth(j+1));
    }

  for(int j=0;j!=hDataLep2Pt->GetNbinsX();++j)
    {
      hDataLep2Pt->SetBinContent(j+1,hDataLep2Pt->GetBinContent(j+1)/hDataLep2Pt->GetBinWidth(j+1));
      hDataLep2Pt->SetBinError(j+1,hDataLep2Pt->GetBinError(j+1)/hDataLep2Pt->GetBinWidth(j+1));
      hMCLep2Pt->SetBinContent(j+1,hMCLep2Pt->GetBinContent(j+1)/hMCLep2Pt->GetBinWidth(j+1));
      hMCLep2Pt->SetBinError(j+1,hMCLep2Pt->GetBinError(j+1)/hMCLep2Pt->GetBinWidth(j+1));
      hZeeLep2Pt->SetBinContent(j+1,hZeeLep2Pt->GetBinContent(j+1)/hZeeLep2Pt->GetBinWidth(j+1));
      hZeeLep2Pt->SetBinError(j+1,hZeeLep2Pt->GetBinError(j+1)/hZeeLep2Pt->GetBinWidth(j+1));
      hEWKLep2Pt->SetBinContent(j+1,hEWKLep2Pt->GetBinContent(j+1)/hEWKLep2Pt->GetBinWidth(j+1));
      hEWKLep2Pt->SetBinError(j+1,hEWKLep2Pt->GetBinError(j+1)/hEWKLep2Pt->GetBinWidth(j+1));
      hTopLep2Pt->SetBinContent(j+1,hTopLep2Pt->GetBinContent(j+1)/hTopLep2Pt->GetBinWidth(j+1));
      hTopLep2Pt->SetBinError(j+1,hTopLep2Pt->GetBinError(j+1)/hTopLep2Pt->GetBinWidth(j+1));
    }

  for(int j=0;j!=hDataLepNegPt->GetNbinsX();++j)
    {
      hDataLepNegPt->SetBinContent(j+1,hDataLepNegPt->GetBinContent(j+1)/hDataLepNegPt->GetBinWidth(j+1));
      hDataLepNegPt->SetBinError(j+1,hDataLepNegPt->GetBinError(j+1)/hDataLepNegPt->GetBinWidth(j+1));
      hMCLepNegPt->SetBinContent(j+1,hMCLepNegPt->GetBinContent(j+1)/hMCLepNegPt->GetBinWidth(j+1));
      hMCLepNegPt->SetBinError(j+1,hMCLepNegPt->GetBinError(j+1)/hMCLepNegPt->GetBinWidth(j+1));
      hZeeLepNegPt->SetBinContent(j+1,hZeeLepNegPt->GetBinContent(j+1)/hZeeLepNegPt->GetBinWidth(j+1));
      hZeeLepNegPt->SetBinError(j+1,hZeeLepNegPt->GetBinError(j+1)/hZeeLepNegPt->GetBinWidth(j+1));
      hEWKLepNegPt->SetBinContent(j+1,hEWKLepNegPt->GetBinContent(j+1)/hEWKLepNegPt->GetBinWidth(j+1));
      hEWKLepNegPt->SetBinError(j+1,hEWKLepNegPt->GetBinError(j+1)/hEWKLepNegPt->GetBinWidth(j+1));
      hTopLepNegPt->SetBinContent(j+1,hTopLepNegPt->GetBinContent(j+1)/hTopLepNegPt->GetBinWidth(j+1));
      hTopLepNegPt->SetBinError(j+1,hTopLepNegPt->GetBinError(j+1)/hTopLepNegPt->GetBinWidth(j+1));
    }

  for(int j=0;j!=hDataLepPosPt->GetNbinsX();++j)
    {
      hDataLepPosPt->SetBinContent(j+1,hDataLepPosPt->GetBinContent(j+1)/hDataLepPosPt->GetBinWidth(j+1));
      hDataLepPosPt->SetBinError(j+1,hDataLepPosPt->GetBinError(j+1)/hDataLepPosPt->GetBinWidth(j+1));
      hMCLepPosPt->SetBinContent(j+1,hMCLepPosPt->GetBinContent(j+1)/hMCLepPosPt->GetBinWidth(j+1));
      hMCLepPosPt->SetBinError(j+1,hMCLepPosPt->GetBinError(j+1)/hMCLepPosPt->GetBinWidth(j+1));
      hZeeLepPosPt->SetBinContent(j+1,hZeeLepPosPt->GetBinContent(j+1)/hZeeLepPosPt->GetBinWidth(j+1));
      hZeeLepPosPt->SetBinError(j+1,hZeeLepPosPt->GetBinError(j+1)/hZeeLepPosPt->GetBinWidth(j+1));
      hEWKLepPosPt->SetBinContent(j+1,hEWKLepPosPt->GetBinContent(j+1)/hEWKLepPosPt->GetBinWidth(j+1));
      hEWKLepPosPt->SetBinError(j+1,hEWKLepPosPt->GetBinError(j+1)/hEWKLepPosPt->GetBinWidth(j+1));
      hTopLepPosPt->SetBinContent(j+1,hTopLepPosPt->GetBinContent(j+1)/hTopLepPosPt->GetBinWidth(j+1));
      hTopLepPosPt->SetBinError(j+1,hTopLepPosPt->GetBinError(j+1)/hTopLepPosPt->GetBinWidth(j+1));
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

  TH1D *hZeeDiff_EGcomp = makeDiffHist(hDataEG,hZeeEG,"hZeeDiff_EGcomp");
  hZeeDiff_EGcomp->SetMarkerStyle(kFullCircle); 
  hZeeDiff_EGcomp->SetMarkerSize(0.9);

  TH1D *hZeeNPVDiff = makeDiffHist(hDataNPV,hMCNPV,"hZeeNPVDiff");
  hZeeNPVDiff->SetMarkerStyle(kFullCircle); 
  hZeeNPVDiff->SetMarkerSize(0.9);

  TH1D *hZeeZPtDiff = makeDiffHist(hDataZPt,hMCZPt,"hZeeZPtDiff");
  hZeeZPtDiff->SetMarkerStyle(kFullCircle); 
  hZeeZPtDiff->SetMarkerSize(0.9);

  TH1D *hZeePhiStarDiff = makeDiffHist(hDataPhiStar,hMCPhiStar,"hZeePhiStarDiff");
  hZeePhiStarDiff->SetMarkerStyle(kFullCircle); 
  hZeePhiStarDiff->SetMarkerSize(0.9);

  TH1D *hZeeZRapDiff = makeDiffHist(hDataZRap,hMCZRap,"hZeeZRapDiff");
  hZeeZRapDiff->SetMarkerStyle(kFullCircle); 
  hZeeZRapDiff->SetMarkerSize(0.9);

  TH1D *hZeeLep1PtDiff = makeDiffHist(hDataLep1Pt,hMCLep1Pt,"hZeeLep1PtDiff");
  hZeeLep1PtDiff->SetMarkerStyle(kFullCircle); 
  hZeeLep1PtDiff->SetMarkerSize(0.9);

  TH1D *hZeeLep2PtDiff = makeDiffHist(hDataLep2Pt,hMCLep2Pt,"hZeeLep2PtDiff");
  hZeeLep2PtDiff->SetMarkerStyle(kFullCircle); 
  hZeeLep2PtDiff->SetMarkerSize(0.9);

  TH1D *hZeeLepNegPtDiff = makeDiffHist(hDataLepNegPt,hMCLepNegPt,"hZeeLepNegPtDiff");
  hZeeLepNegPtDiff->SetMarkerStyle(kFullCircle); 
  hZeeLepNegPtDiff->SetMarkerSize(0.9);

  TH1D *hZeeLepPosPtDiff = makeDiffHist(hDataLepPosPt,hMCLepPosPt,"hZeeLepPosPtDiff");
  hZeeLepPosPtDiff->SetMarkerStyle(kFullCircle); 
  hZeeLepPosPtDiff->SetMarkerSize(0.9);

  TH1D *hZeeLep1EtaDiff = makeDiffHist(hDataLep1Eta,hMCLep1Eta,"hZeeLep1EtaDiff");
  hZeeLep1EtaDiff->SetMarkerStyle(kFullCircle); 
  hZeeLep1EtaDiff->SetMarkerSize(0.9);

  TH1D *hZeeLep2EtaDiff = makeDiffHist(hDataLep2Eta,hMCLep2Eta,"hZeeLep2EtaDiff");
  hZeeLep2EtaDiff->SetMarkerStyle(kFullCircle); 
  hZeeLep2EtaDiff->SetMarkerSize(0.9);
  
  
  TH1D *hZeeZRapDiffUnc = new TH1D("hZeeZRapDiffUnc","hZeeZRapDiffUnc",24, 0.0, 2.4);
  TH1D *hLine = new TH1D("hLine","hLine", 24, 0.0, 2.4);
  
	for(int i =1 ; i <= NBINS ; ++i){
		hZeeZRapDiffUnc->SetBinContent(i,0);
		hLine->SetBinContent(i,0);
		hLine->SetBinError(i,0);
    // double dataup = (hData->GetBinContent(i)-hZeeZRapDiffUp->GetBinContent(i))/hData->GetBinContent(i);
    // double datadown = (hData->GetBinContent(i)-hZeeZRapDiffDown->GetBinContent(i))/hData->GetBinContent(i);
    double zeeup = (hZeeZRap->GetBinContent(i)-hZeeZRapUp->GetBinContent(i))/hZeeZRap->GetBinContent(i);
    double zeedown = (hZeeZRap->GetBinContent(i)-hZeeZRapDown->GetBinContent(i))/hZeeZRap->GetBinContent(i);
		// massUnc->SetBinError(i,sqrt(dataup*dataup+datadown*datadown+zeeup*zeeup+zeedown*zeedown));
		hZeeZRapDiffUnc->SetBinError(i,sqrt(zeeup*zeeup+zeedown*zeedown));
    cout << "Bin Error is " << hZeeZRapDiffUnc->GetBinError(i) << endl;
	}
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  char ylabel[100];     // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000);  
  
  char normtext[100];
  // sprintf(normtext,"MC normalized to data (#times %.2f)",MCscale);  
  sprintf(normtext,"",MCscale);  

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
  
  
  
  std::cout << "peak dat " << hData->GetBinCenter(hData->GetMaximumBin()) << std::endl;
  std::cout << "peak zee " << hZee->GetBinCenter(hZee->GetMaximumBin()) << std::endl;
  
    sprintf(ylabel,"Events / %.1f GeV",hGausRandHere->GetBinWidth(1));
  CPlot plotGaus("gaus_rand"+norm,"","",ylabel);
  plotGaus.AddHist1D(hGausRandHere,"this code","E");
  // plotZee.AddToStack(hZee,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotGaus.AddToStack(hGausRandNtuple,"ntple",fillcolorZ,linecolorZ);
  plotGaus.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotGaus.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotGaus.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotGaus.SetYRange(0.01,1.2*(hGausRandHere->GetMaximum() + sqrt(hGausRandHere->GetMaximum())));
  plotGaus.TransLegend(0.1,-0.05);
  plotGaus.Draw(c,kFALSE,format,1);

  CPlot plotGausDiff("gaus_rand"+norm,"","M(e^{+}e^{-}) [GeV]","#frac{Data-Pred}{Data}");
  plotGausDiff.AddHist1D(hGausDiff,"EX0",ratioColor);
  plotGausDiff.SetYRange(-0.2,0.2);
  plotGausDiff.AddLine(-5, 0,5, 0,kBlack,1);
  plotGausDiff.AddLine(-5, 0.1,5, 0.1,kBlack,3);
  plotGausDiff.AddLine(-5,-0.1,5,-0.1,kBlack,3);
  plotGausDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotGaus2("gaus_rand_log"+norm,"","",ylabel);
  plotGaus2.AddHist1D(hGausRandHere,"mine","E");
  // plotZee2.AddToStack(hEWK,"EWK",fillcolorEWK,linecolorEWK);
  // plotZee2.AddToStack(hTop,"t#bar{t}",fillcolorTop,linecolorTop);
  plotGaus2.AddToStack(hGausRandNtuple,"yours",fillcolorZ,linecolorZ);
  plotGaus2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotGaus2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotGaus2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotGaus2.SetLogy();
  plotGaus2.SetYRange(1e-4*(hGausRandHere->GetMaximum()),10*(hGausRandHere->GetMaximum()));
  plotGaus2.TransLegend(0.1,-0.05);
  plotGaus2.Draw(c,kTRUE,format,1);
  
  
  
  TH1D *hDatDiff = makeDiffHist(hData,hDataEG,"hDatDiff");
  hDatDiff->SetMarkerStyle(kFullCircle); 
  hDatDiff->SetMarkerSize(0.9);
  
  // Shitty plots: Compare data
  sprintf(ylabel,"Events / %.1f GeV",hData->GetBinWidth(1));
  CPlot plotZeea("DATA_zee"+norm,"","",ylabel);
  plotZeea.AddHist1D(hData,"mine","E");
  // plotZee.AddToStack(hZee,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeea.AddToStack(hDataEG,"yours",fillcolorZ,linecolorZ);
  plotZeea.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeea.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeea.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZeea.SetYRange(0.01,1.2*(hData->GetMaximum() + sqrt(hData->GetMaximum())));
  plotZeea.TransLegend(0.1,-0.05);
  plotZeea.Draw(c,kFALSE,format,1);

  CPlot plotZeeDiffa("DATA_zee"+norm,"","M(e^{+}e^{-}) [GeV]","#frac{Data-Pred}{Data}");
  plotZeeDiffa.AddHist1D(hDatDiff,"EX0",ratioColor);
  plotZeeDiffa.SetYRange(-0.2,0.2);
  plotZeeDiffa.AddLine(MASS_LOW, 0,MASS_HIGH, 0,kBlack,1);
  plotZeeDiffa.AddLine(MASS_LOW, 0.1,MASS_HIGH, 0.1,kBlack,3);
  plotZeeDiffa.AddLine(MASS_LOW,-0.1,MASS_HIGH,-0.1,kBlack,3);
  plotZeeDiffa.Draw(c,kTRUE,format,2);
  
  CPlot plotZee2a("DATA_zeelog"+norm,"","",ylabel);
  plotZee2a.AddHist1D(hData,"mine","E");
  // plotZee2.AddToStack(hEWK,"EWK",fillcolorEWK,linecolorEWK);
  // plotZee2.AddToStack(hTop,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZee2a.AddToStack(hDataEG,"yours",fillcolorZ,linecolorZ);
  plotZee2a.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZee2a.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZee2a.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZee2a.SetLogy();
  plotZee2a.SetYRange(1e-4*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plotZee2a.TransLegend(0.1,-0.05);
  plotZee2a.Draw(c,kTRUE,format,1);
  
  
  //-------------------------------
  
  TH1D *hEGMCDiff = makeDiffHist(hZee,hZeeEG,"hEGMCDiff");
  hEGMCDiff->SetMarkerStyle(kFullCircle); 
  hEGMCDiff->SetMarkerSize(0.9);
  
    // Shitty plots: Compare MC
  sprintf(ylabel,"Events / %.1f GeV",hZee->GetBinWidth(1));
  CPlot plotZeeb("MC_zee"+norm,"","",ylabel);
  plotZeeb.AddHist1D(hZee,"mine","E");
  // plotZee.AddToStack(hZee,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeb.AddToStack(hZeeEG,"yours",fillcolorZ,linecolorZ);
  plotZeeb.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeb.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeb.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZeeb.SetYRange(0.01,1.2*(hZeeEG->GetMaximum() + sqrt(hZeeEG->GetMaximum())));
  plotZeeb.TransLegend(0.1,-0.05);
  plotZeeb.Draw(c,kFALSE,format,1);

  CPlot plotZeeDiffb("MC_zee"+norm,"","M(e^{+}e^{-}) [GeV]","#frac{Data-Pred}{Data}");
  plotZeeDiffb.AddHist1D(hEGMCDiff,"EX0",ratioColor);
  plotZeeDiffb.SetYRange(-0.2,0.2);
  plotZeeDiffb.AddLine(MASS_LOW, 0,MASS_HIGH, 0,kBlack,1);
  plotZeeDiffb.AddLine(MASS_LOW, 0.1,MASS_HIGH, 0.1,kBlack,3);
  plotZeeDiffb.AddLine(MASS_LOW,-0.1,MASS_HIGH,-0.1,kBlack,3);
  plotZeeDiffb.Draw(c,kTRUE,format,2);
  
  CPlot plotZee2b("MC_zeelog"+norm,"","",ylabel);
  plotZee2b.AddHist1D(hZee,"mine","E");
  // plotZee2.AddToStack(hEWK,"EWK",fillcolorEWK,linecolorEWK);
  // plotZee2.AddToStack(hTop,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZee2b.AddToStack(hZeeEG,"yours",fillcolorZ,linecolorZ);
  plotZee2b.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZee2b.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZee2b.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZee2b.SetLogy();
  plotZee2b.SetYRange(1e-4*(hZeeEG->GetMaximum()),10*(hZeeEG->GetMaximum()));
  plotZee2b.TransLegend(0.1,-0.05);
  plotZee2b.Draw(c,kTRUE,format,1);
  
  
  // TH1D *hRawEGDiff = makeDiffHist(hDataEG,hZeeEG,"hRawEGDiff");
  // hRawEGDiff->SetMarkerStyle(kFullCircle); 
  // hRawEGDiff->SetMarkerSize(0.9);
  
  // // Shitty plots: Compare data
  // sprintf(ylabel,"Events / %.1f GeV",hData->GetBinWidth(1));
  // CPlot plotZeea("EG_zee"+norm,"","",ylabel);
  // plotZeea.AddHist1D(hDataEG,"dat","E");
  // // plotZee.AddToStack(hZee,"Z#rightarrowee",fillcolorZ,linecolorZ);
  // plotZeea.AddToStack(hZeeEG,"mc",fillcolorZ,linecolorZ);
  // plotZeea.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  // plotZeea.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  // if(normToData)plotZeea.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  // plotZeea.SetYRange(0.01,1.2*(hData->GetMaximum() + sqrt(hData->GetMaximum())));
  // plotZeea.TransLegend(0.1,-0.05);
  // plotZeea.Draw(c,kFALSE,format,1);

  // CPlot plotZeeDiffa("EG_zee"+norm,"","M(e^{+}e^{-}) [GeV]","#frac{Data-Pred}{Data}");
  // plotZeeDiffa.AddHist1D(hRawEGDiff,"EX0",ratioColor);
  // plotZeeDiffa.SetYRange(-0.2,0.2);
  // plotZeeDiffa.AddLine(MASS_LOW, 0,MASS_HIGH, 0,kBlack,1);
  // plotZeeDiffa.AddLine(MASS_LOW, 0.1,MASS_HIGH, 0.1,kBlack,3);
  // plotZeeDiffa.AddLine(MASS_LOW,-0.1,MASS_HIGH,-0.1,kBlack,3);
  // plotZeeDiffa.Draw(c,kTRUE,format,2);
  
  // CPlot plotZee2a("EG_zeelog"+norm,"","",ylabel);
  // plotZee2a.AddHist1D(hDataEG,"dat","E");
  // // plotZee2.AddToStack(hEWK,"EWK",fillcolorEWK,linecolorEWK);
  // // plotZee2.AddToStack(hTop,"t#bar{t}",fillcolorTop,linecolorTop);
  // plotZee2a.AddToStack(hZeeEG,"mc",fillcolorZ,linecolorZ);
  // plotZee2a.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  // plotZee2a.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  // if(normToData)plotZee2a.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  // plotZee2a.SetLogy();
  // plotZee2a.SetYRange(1e-4*(hDataEG->GetMaximum()),10*(hDataEG->GetMaximum()));
  // plotZee2a.TransLegend(0.1,-0.05);
  // plotZee2a.Draw(c,kTRUE,format,1);
  
  
  // //-------------------------------
  
  // TH1D *hRawDiff = makeDiffHist(hData,hZee,"hRawDiff");
  // hRawDiff->SetMarkerStyle(kFullCircle); 
  // hRawDiff->SetMarkerSize(0.9);
  
    // // Shitty plots: Compare MC
  // sprintf(ylabel,"Events / %.1f GeV",hZee->GetBinWidth(1));
  // CPlot plotZeeb("me_zee"+norm,"","",ylabel);
  // plotZeeb.AddHist1D(hData,"dat","E");
  // // plotZee.AddToStack(hZee,"Z#rightarrowee",fillcolorZ,linecolorZ);
  // plotZeeb.AddToStack(hZee,"mc",fillcolorZ,linecolorZ);
  // plotZeeb.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  // plotZeeb.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  // if(normToData)plotZeeb.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  // plotZeeb.SetYRange(0.01,1.2*(hZeeEG->GetMaximum() + sqrt(hZeeEG->GetMaximum())));
  // plotZeeb.TransLegend(0.1,-0.05);
  // plotZeeb.Draw(c,kFALSE,format,1);

  // CPlot plotZeeDiffb("me_zee"+norm,"","M(e^{+}e^{-}) [GeV]","#frac{Data-Pred}{Data}");
  // plotZeeDiffb.AddHist1D(hRawDiff,"EX0",ratioColor);
  // plotZeeDiffb.SetYRange(-0.2,0.2);
  // plotZeeDiffb.AddLine(MASS_LOW, 0,MASS_HIGH, 0,kBlack,1);
  // plotZeeDiffb.AddLine(MASS_LOW, 0.1,MASS_HIGH, 0.1,kBlack,3);
  // plotZeeDiffb.AddLine(MASS_LOW,-0.1,MASS_HIGH,-0.1,kBlack,3);
  // plotZeeDiffb.Draw(c,kTRUE,format,2);
  
  // CPlot plotZee2b("me_zeelog"+norm,"","",ylabel);
  // plotZee2b.AddHist1D(hData,"dat","E");
  // // plotZee2.AddToStack(hEWK,"EWK",fillcolorEWK,linecolorEWK);
  // // plotZee2.AddToStack(hTop,"t#bar{t}",fillcolorTop,linecolorTop);
  // plotZee2b.AddToStack(hZee,"mc",fillcolorZ,linecolorZ);
  // plotZee2b.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  // plotZee2b.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  // if(normToData)plotZee2b.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  // plotZee2b.SetLogy();
  // plotZee2b.SetYRange(1e-4*(hData->GetMaximum()),10*(hData->GetMaximum()));
  // plotZee2b.TransLegend(0.1,-0.05);
  // plotZee2b.Draw(c,kTRUE,format,1);
  
  
  /////////////////////////////////////////////////////////////
  
  // EG
  sprintf(ylabel,"Events / %.1f GeV",hData->GetBinWidth(1));
  CPlot plotZee_EGcomp("zee_EGcomp"+norm,"","",ylabel);
  plotZee_EGcomp.AddHist1D(hDataEG,"data","E");
  plotZee_EGcomp.AddToStack(hZeeEG,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZee_EGcomp.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZee_EGcomp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZee_EGcomp.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZee_EGcomp.SetYRange(0.01,1.2*(hDataEG->GetMaximum() + sqrt(hDataEG->GetMaximum())));
  plotZee_EGcomp.TransLegend(0.1,-0.05);
  plotZee_EGcomp.Draw(c,kFALSE,format,1);

  CPlot plotZeeDiff_EGcomp("zee_EGcomp"+norm,"","M(e^{+}e^{-}) [GeV]","#frac{Data-Pred}{Data}");
  plotZeeDiff_EGcomp.AddHist1D(hZeeDiff_EGcomp,"EX0",ratioColor);
  plotZeeDiff_EGcomp.SetYRange(-0.2,0.2);
  plotZeeDiff_EGcomp.AddLine(MASS_LOW, 0,MASS_HIGH, 0,kBlack,1);
  plotZeeDiff_EGcomp.AddLine(MASS_LOW, 0.1,MASS_HIGH, 0.1,kBlack,3);
  plotZeeDiff_EGcomp.AddLine(MASS_LOW,-0.1,MASS_HIGH,-0.1,kBlack,3);
  plotZeeDiff_EGcomp.Draw(c,kTRUE,format,2);
  
  CPlot plotZee2_EGcomp("zeelog_EGcomp"+norm,"","",ylabel);
  plotZee2_EGcomp.AddHist1D(hDataEG,"data","E");
  // plotZee2.AddToStack(hEWK,"EWK",fillcolorEWK,linecolorEWK);
  // plotZee2.AddToStack(hTop,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZee2_EGcomp.AddToStack(hZeeEG,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZee2_EGcomp.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZee2_EGcomp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZee2_EGcomp.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZee2_EGcomp.SetLogy();
  plotZee2_EGcomp.SetYRange(1e-4*(hDataEG->GetMaximum()),10*(hDataEG->GetMaximum()));
  plotZee2_EGcomp.TransLegend(0.1,-0.05);
  plotZee2_EGcomp.Draw(c,kTRUE,format,1);
  // //
  
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
  plotZeeDiff.AddHist1D(hZeeDiff,"EX0",ratioColor);
  plotZeeDiff.SetYRange(-0.2,0.2);
  plotZeeDiff.AddLine(MASS_LOW, 0,MASS_HIGH, 0,kBlack,1);
  plotZeeDiff.AddLine(MASS_LOW, 0.1,MASS_HIGH, 0.1,kBlack,3);
  plotZeeDiff.AddLine(MASS_LOW,-0.1,MASS_HIGH,-0.1,kBlack,3);
  plotZeeDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZee2("zeelog"+norm,"","",ylabel);
  plotZee2.AddHist1D(hData,"data","E");
  plotZee2.AddToStack(hEWK,"EWK",fillcolorEWK,linecolorEWK);
  plotZee2.AddToStack(hTop,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZee2.AddToStack(hZee,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZee2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZee2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZee2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZee2.SetLogy();
  plotZee2.SetYRange(1e-4*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plotZee2.TransLegend(0.1,-0.05);
  plotZee2.Draw(c,kTRUE,format,1);

  //
  // NPV
  // 

  sprintf(ylabel,"Events");
  CPlot plotZeeNPV("zeeNPV"+norm,"","",ylabel);
  plotZeeNPV.AddHist1D(hDataNPV,"data","E");
  plotZeeNPV.AddToStack(hZeeNPV,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeNPV.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeNPV.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeNPV.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZeeNPV.SetYRange(0.01,1.2*(hData->GetMaximum() + sqrt(hData->GetMaximum())));
  plotZeeNPV.TransLegend(0.1,-0.05);
  plotZeeNPV.Draw(c,kFALSE,format,1);

  CPlot plotZeeNPVDiff("zeeNPV"+norm,"","npv","#frac{Data-Pred}{Data}");
  plotZeeNPVDiff.AddHist1D(hZeeNPVDiff,"EX0",ratioColor);
  plotZeeNPVDiff.SetYRange(-0.2,0.2);
  plotZeeNPVDiff.AddLine(0, 0,50, 0,kBlack,1);
  plotZeeNPVDiff.AddLine(0, 0.1,50, 0.1,kBlack,3);
  plotZeeNPVDiff.AddLine(0,-0.1,50,-0.1,kBlack,3);
  plotZeeNPVDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZeeNPV2("zeeNPVlog"+norm,"","",ylabel);
  plotZeeNPV2.AddHist1D(hDataNPV,"data","E");
  plotZeeNPV2.AddToStack(hEWKNPV,"EWK",fillcolorEWK,linecolorEWK);
  plotZeeNPV2.AddToStack(hTopNPV,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZeeNPV2.AddToStack(hZeeNPV,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeNPV2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeNPV2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeNPV2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZeeNPV2.SetLogy();
  plotZeeNPV2.SetYRange(1e-4*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plotZeeNPV2.TransLegend(0.1,-0.05);
  plotZeeNPV2.Draw(c,kTRUE,format,1);

  //
  // Z Pt
  //   
  sprintf(ylabel,"Events / 1 GeV");
  CPlot plotZeeZPt("zeeZPt"+norm,"","",ylabel);
  plotZeeZPt.AddHist1D(hDataZPt,"data","E");
  plotZeeZPt.AddToStack(hZeeZPt,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeZPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeZPt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeZPt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  // plotZeeZPt.SetLogx();
  plotZeeZPt.SetLogy(0);
  plotZeeZPt.SetYRange(0.01,1.2*(hDataZPt->GetMaximum() + sqrt(hDataZPt->GetMaximum())));
  plotZeeZPt.TransLegend(0.1,-0.05);
  plotZeeZPt.Draw(c,kFALSE,format,1);

  CPlot plotZeeZPtDiff("zeeZPt"+norm,"","p_{T}^{e^{+}e^{-}} [GeV]","#frac{Data-Pred}{Data}");
  plotZeeZPtDiff.AddHist1D(hZeeZPtDiff,"EX0",ratioColor);
  // plotZeeZPtDiff.SetLogx();
  plotZeeZPtDiff.SetYRange(-0.2,0.2);
  plotZeeZPtDiff.AddLine(0, 0,1000, 0,kBlack,1);
  plotZeeZPtDiff.AddLine(0, 0.1,1000, 0.1,kBlack,3);
  plotZeeZPtDiff.AddLine(0,-0.1,1000,-0.1,kBlack,3);
  plotZeeZPtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZeeZPt2("zeeZPtlog"+norm,"","",ylabel);
  plotZeeZPt2.AddHist1D(hDataZPt,"data","E");
  plotZeeZPt2.AddToStack(hEWKZPt,"EWK",fillcolorEWK,linecolorEWK);
  plotZeeZPt2.AddToStack(hTopZPt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZeeZPt2.AddToStack(hZeeZPt,"Z#rightarrowee",fillcolorZ,linecolorZ);\
  plotZeeZPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeZPt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeZPt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  // plotZeeZPt2.SetLogx();
  plotZeeZPt2.SetLogy();
  plotZeeZPt2.SetYRange(1e-6*(hDataZPt->GetMaximum()),10*(hDataZPt->GetMaximum()));
  plotZeeZPt2.TransLegend(0.1,-0.05);
  plotZeeZPt2.Draw(c,kTRUE,format,1);

  //
  // Phi*
  //   
  sprintf(ylabel,"Events / 1.0");
  CPlot plotZeePhiStar("zeePhiStar"+norm,"","",ylabel);
  plotZeePhiStar.AddHist1D(hDataPhiStar,"data","E");
  plotZeePhiStar.AddToStack(hZeePhiStar,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeePhiStar.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeePhiStar.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeePhiStar.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  // plotZeePhiStar.SetLogx();
  plotZeePhiStar.SetLogy(0);
  plotZeePhiStar.SetYRange(0.01,1.2*(hDataPhiStar->GetMaximum() + sqrt(hDataPhiStar->GetMaximum())));
  plotZeePhiStar.TransLegend(0.1,-0.05);
  plotZeePhiStar.Draw(c,kFALSE,format,1);

  CPlot plotZeePhiStarDiff("zeePhiStar"+norm,"","#phi_{#eta}*","#frac{Data-Pred}{Data}");
  plotZeePhiStarDiff.AddHist1D(hZeePhiStarDiff,"EX0",ratioColor);
  // plotZeePhiStarDiff.SetLogx();
  plotZeePhiStarDiff.SetYRange(-0.2,0.2);
  plotZeePhiStarDiff.AddLine(0, 0,3, 0,kBlack,1);
  plotZeePhiStarDiff.AddLine(0, 0.1,3, 0.1,kBlack,3);
  plotZeePhiStarDiff.AddLine(0,-0.1,3,-0.1,kBlack,3);
  plotZeePhiStarDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZeePhiStar2("zeePhiStarlog"+norm,"","",ylabel);
  plotZeePhiStar2.AddHist1D(hDataPhiStar,"data","E");
  plotZeePhiStar2.AddToStack(hEWKPhiStar,"EWK",fillcolorEWK,linecolorEWK);
  plotZeePhiStar2.AddToStack(hTopPhiStar,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZeePhiStar2.AddToStack(hZeePhiStar,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeePhiStar2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeePhiStar2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeePhiStar2.AddTextBox(normtext,0.21,0.75,0.45,0.80,0,13,0.03,-1);
  // plotZeePhiStar2.SetLogx();
  plotZeePhiStar2.SetLogy();
  plotZeePhiStar2.SetYRange(1e-5*(hDataPhiStar->GetMaximum()),10*(hDataPhiStar->GetMaximum()));
  plotZeePhiStar2.TransLegend(0.1,-0.05);
  plotZeePhiStar2.Draw(c,kTRUE,format,1);

  //
  // Z Rapidity
  //   
  sprintf(ylabel,"Events / %.1f ",hDataZRap->GetBinWidth(1));
  CPlot plotZeeZRap("zeeZRap"+norm,"","",ylabel);
  plotZeeZRap.AddHist1D(hDataZRap,"data","E");
  plotZeeZRap.AddToStack(hZeeZRap,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeZRap.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeZRap.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeZRap.AddTextBox(normtext,0.21,0.75,0.45,0.80,0,13,0.03,-1);
  plotZeeZRap.SetLogy(0);
  plotZeeZRap.SetYRange(0.01,1.25*(hDataZRap->GetMaximum() + sqrt(hDataZRap->GetMaximum())));
  plotZeeZRap.TransLegend(0.1,-0.05);
  plotZeeZRap.Draw(c,kFALSE,format,1);

  CPlot plotZeeZRapDiff("zeeZRap"+norm,"","|y^{e^{+}e^{-}}|","#frac{Data-Pred}{Data}");
  plotZeeZRapDiff.AddHist1D(hZeeZRapDiffUnc,"E3",kGray,1,1);
  plotZeeZRapDiff.AddHist1D(hLine,"E3",kBlack,1,1);
  plotZeeZRapDiff.AddHist1D(hZeeZRapDiff,"EX0",ratioColor);
  plotZeeZRapDiff.SetYRange(-0.2,0.2);
  plotZeeZRapDiff.AddLine(0, 0,2.4, 0,kBlack,1);
  plotZeeZRapDiff.AddLine(0, 0.1,2.4, 0.1,kBlack,3);
  plotZeeZRapDiff.AddLine(0,-0.1,2.4,-0.1,kBlack,3);
  plotZeeZRapDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZeeZRap2("zeeZRaplog"+norm,"","",ylabel);
  plotZeeZRap2.AddHist1D(hDataZRap,"data","E");
  plotZeeZRap2.AddToStack(hEWKZRap,"EWK",fillcolorEWK,linecolorEWK);
  plotZeeZRap2.AddToStack(hTopZRap,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZeeZRap2.AddToStack(hZeeZRap,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeZRap2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeZRap2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeZRap2.AddTextBox(normtext,0.21,0.75,0.45,0.80,0,13,0.03,-1);
  plotZeeZRap2.SetLogy();
  plotZeeZRap2.SetYRange(1e-4*(hDataZRap->GetMaximum()),500*(hDataZRap->GetMaximum()));
  plotZeeZRap2.TransLegend(0.1,-0.05);
  plotZeeZRap2.Draw(c,kTRUE,format,1);

 //
  // Lep1 Pt
  //   
  sprintf(ylabel,"Events / 1 GeV");
  CPlot plotZeeLep1Pt("zeeLep1Pt"+norm,"","",ylabel);
  plotZeeLep1Pt.AddHist1D(hDataLep1Pt,"data","E");
  plotZeeLep1Pt.AddToStack(hZeeLep1Pt,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeLep1Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLep1Pt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLep1Pt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  // plotZeeLep1Pt.SetLogx();
  plotZeeLep1Pt.SetLogy(0);
  plotZeeLep1Pt.SetYRange(0.01,1.2*(hDataLep1Pt->GetMaximum() + sqrt(hDataLep1Pt->GetMaximum())));
  plotZeeLep1Pt.TransLegend(0.1,-0.05);
  plotZeeLep1Pt.Draw(c,kFALSE,format,1);

  CPlot plotZeeLep1PtDiff("zeeLep1Pt"+norm,"","p_{T}(leading electron) [GeV]","#frac{Data-Pred}{Data}");
  plotZeeLep1PtDiff.AddHist1D(hZeeLep1PtDiff,"EX0",ratioColor);
  // plotZeeLep1PtDiff.SetLogx();
  plotZeeLep1PtDiff.SetYRange(-0.2,0.2);
  plotZeeLep1PtDiff.AddLine(25, 0,300, 0,kBlack,1);
  plotZeeLep1PtDiff.AddLine(25, 0.1,300, 0.1,kBlack,3);
  plotZeeLep1PtDiff.AddLine(25,-0.1,300,-0.1,kBlack,3);
  plotZeeLep1PtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZeeLep1Pt2("zeeLep1Ptlog"+norm,"","",ylabel);
  plotZeeLep1Pt2.AddHist1D(hDataLep1Pt,"data","E");
  plotZeeLep1Pt2.AddToStack(hEWKLep1Pt,"EWK",fillcolorEWK,linecolorEWK);
  plotZeeLep1Pt2.AddToStack(hTopLep1Pt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZeeLep1Pt2.AddToStack(hZeeLep1Pt,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeLep1Pt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLep1Pt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLep1Pt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  // plotZeeLep1Pt2.SetLogx();
  plotZeeLep1Pt2.SetLogy();
  plotZeeLep1Pt2.SetYRange(1e-5*(hDataLep1Pt->GetMaximum()),10*(hDataLep1Pt->GetMaximum()));
  plotZeeLep1Pt2.TransLegend(0.1,-0.05);
  plotZeeLep1Pt2.Draw(c,kTRUE,format,1);

  //
  // Lep2 Pt
  //   
  sprintf(ylabel,"Events / 1 GeV");
  CPlot plotZeeLep2Pt("zeeLep2Pt"+norm,"","",ylabel);
  plotZeeLep2Pt.AddHist1D(hDataLep2Pt,"data","E");
  plotZeeLep2Pt.AddToStack(hZeeLep2Pt,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeLep2Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLep2Pt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLep2Pt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  // plotZeeLep2Pt.SetLogx();
  plotZeeLep2Pt.SetLogy(0);
  plotZeeLep2Pt.SetYRange(0.01,1.2*(hDataLep2Pt->GetMaximum() + sqrt(hDataLep2Pt->GetMaximum())));
  plotZeeLep2Pt.TransLegend(0.1,-0.05);
  plotZeeLep2Pt.Draw(c,kFALSE,format,1);

  CPlot plotZeeLep2PtDiff("zeeLep2Pt"+norm,"","p_{T}(2nd leading electron) [GeV]","#frac{Data-Pred}{Data}");
  plotZeeLep2PtDiff.AddHist1D(hZeeLep2PtDiff,"EX0",ratioColor);
  // plotZeeLep2PtDiff.SetLogx();
  plotZeeLep2PtDiff.SetYRange(-0.2,0.2);
  plotZeeLep2PtDiff.AddLine(25, 0,150, 0,kBlack,1);
  plotZeeLep2PtDiff.AddLine(25, 0.1,150, 0.1,kBlack,3);
  plotZeeLep2PtDiff.AddLine(25,-0.1,150,-0.1,kBlack,3);
  plotZeeLep2PtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZeeLep2Pt2("zeeLep2Ptlog"+norm,"","",ylabel);
  plotZeeLep2Pt2.AddHist1D(hDataLep2Pt,"data","E");
  plotZeeLep2Pt2.AddToStack(hEWKLep2Pt,"EWK",fillcolorEWK,linecolorEWK);
  plotZeeLep2Pt2.AddToStack(hTopLep2Pt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZeeLep2Pt2.AddToStack(hZeeLep2Pt,"Z#rightarrowee",fillcolorZ,linecolorZ);\
  plotZeeLep2Pt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLep2Pt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLep2Pt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  // plotZeeLep2Pt2.SetLogx();
  plotZeeLep2Pt2.SetLogy();
  plotZeeLep2Pt2.SetYRange(1e-5*(hDataLep2Pt->GetMaximum()),10*(hDataLep2Pt->GetMaximum()));
  plotZeeLep2Pt2.TransLegend(0.1,-0.05);
  plotZeeLep2Pt2.Draw(c,kTRUE,format,1);

  //
  // LepNeg Pt
  //   
  sprintf(ylabel,"Events / 1 GeV");
  CPlot plotZeeLepNegPt("zeeLepNegPt"+norm,"","",ylabel);
  plotZeeLepNegPt.AddHist1D(hDataLepNegPt,"data","E");
  plotZeeLepNegPt.AddToStack(hZeeLepNegPt,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeLepNegPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLepNegPt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLepNegPt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  // plotZeeLepNegPt.SetLogx();
  plotZeeLepNegPt.SetLogy(0);
  plotZeeLepNegPt.SetYRange(0.01,1.2*(hDataLepNegPt->GetMaximum() + sqrt(hDataLepNegPt->GetMaximum())));
  plotZeeLepNegPt.TransLegend(0.1,-0.05);
  plotZeeLepNegPt.Draw(c,kFALSE,format,1);

  CPlot plotZeeLepNegPtDiff("zeeLepNegPt"+norm,"","p_{T}^{e^{-}} [GeV]","#frac{Data-Pred}{Data}");
  plotZeeLepNegPtDiff.AddHist1D(hZeeLepNegPtDiff,"EX0",ratioColor);
  // plotZeeLepNegPtDiff.SetLogx();
  plotZeeLepNegPtDiff.SetYRange(-0.2,0.2);
  plotZeeLepNegPtDiff.AddLine(25, 0,300, 0,kBlack,1);
  plotZeeLepNegPtDiff.AddLine(25, 0.1,300, 0.1,kBlack,3);
  plotZeeLepNegPtDiff.AddLine(25,-0.1,300,-0.1,kBlack,3);
  plotZeeLepNegPtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZeeLepNegPt2("zeeLepNegPtlog"+norm,"","",ylabel);
  plotZeeLepNegPt2.AddHist1D(hDataLepNegPt,"data","E");
  plotZeeLepNegPt2.AddToStack(hEWKLepNegPt,"EWK",fillcolorEWK,linecolorEWK);
  plotZeeLepNegPt2.AddToStack(hTopLepNegPt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZeeLepNegPt2.AddToStack(hZeeLepNegPt,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeLepNegPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLepNegPt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLepNegPt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  // plotZeeLepNegPt2.SetLogx();
  plotZeeLepNegPt2.SetLogy();
  plotZeeLepNegPt2.SetYRange(1e-5*(hDataLepNegPt->GetMaximum()),10*(hDataLepNegPt->GetMaximum()));
  plotZeeLepNegPt2.TransLegend(0.1,-0.05);
  plotZeeLepNegPt2.Draw(c,kTRUE,format,1);

  //
  // LepPos Pt
  //   
  sprintf(ylabel,"Events / 1 GeV");
  CPlot plotZeeLepPosPt("zeeLepPosPt"+norm,"","",ylabel);
  plotZeeLepPosPt.AddHist1D(hDataLepPosPt,"data","E");
  plotZeeLepPosPt.AddToStack(hZeeLepPosPt,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeLepPosPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLepPosPt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLepPosPt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  // plotZeeLepPosPt.SetLogx();
  plotZeeLepPosPt.SetLogy(0);
  plotZeeLepPosPt.SetYRange(0.01,1.2*(hDataLepPosPt->GetMaximum() + sqrt(hDataLepPosPt->GetMaximum())));
  plotZeeLepPosPt.TransLegend(0.1,-0.05);
  plotZeeLepPosPt.Draw(c,kFALSE,format,1);

  CPlot plotZeeLepPosPtDiff("zeeLepPosPt"+norm,"","p_{T}^{e^{+}} [GeV]","#frac{Data-Pred}{Data}");
  plotZeeLepPosPtDiff.AddHist1D(hZeeLepPosPtDiff,"EX0",ratioColor);
  // plotZeeLepPosPtDiff.SetLogx();
  plotZeeLepPosPtDiff.SetYRange(-0.2,0.2);
  plotZeeLepPosPtDiff.AddLine(25, 0,300, 0,kBlack,1);
  plotZeeLepPosPtDiff.AddLine(25, 0.1,300, 0.1,kBlack,3);
  plotZeeLepPosPtDiff.AddLine(25,-0.1,300,-0.1,kBlack,3);
  plotZeeLepPosPtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZeeLepPosPt2("zeeLepPosPtlog"+norm,"","",ylabel);
  plotZeeLepPosPt2.AddHist1D(hDataLepPosPt,"data","E");
  plotZeeLepPosPt2.AddToStack(hEWKLepPosPt,"EWK",fillcolorEWK,linecolorEWK);
  plotZeeLepPosPt2.AddToStack(hTopLepPosPt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZeeLepPosPt2.AddToStack(hZeeLepPosPt,"Z#rightarrowee",fillcolorZ,linecolorZ);\
  plotZeeLepPosPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLepPosPt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLepPosPt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  // plotZeeLepPosPt2.SetLogx();
  plotZeeLepPosPt2.SetLogy();
  plotZeeLepPosPt2.SetYRange(1e-5*(hDataLepPosPt->GetMaximum()),10*(hDataLepPosPt->GetMaximum()));
  plotZeeLepPosPt2.TransLegend(0.1,-0.05);
  plotZeeLepPosPt2.Draw(c,kTRUE,format,1);

  //
  // Lep1 Eta
  //   
  sprintf(ylabel,"Events / %.1f ",hDataLep1Eta->GetBinWidth(1));
  CPlot plotZeeLep1Eta("zeeLep1Eta"+norm,"","",ylabel);
  plotZeeLep1Eta.AddHist1D(hDataLep1Eta,"data","E");
  plotZeeLep1Eta.AddToStack(hZeeLep1Eta,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeLep1Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLep1Eta.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLep1Eta.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZeeLep1Eta.SetLogy(0);
  plotZeeLep1Eta.SetYRange(0.01,1.2*(hDataLep1Eta->GetMaximum() + sqrt(hDataLep1Eta->GetMaximum())));
  plotZeeLep1Eta.TransLegend(0.1,-0.05);
  plotZeeLep1Eta.Draw(c,kFALSE,format,1);

  CPlot plotZeeLep1EtaDiff("zeeLep1Eta"+norm,"","|#eta| (leading electron)","#frac{Data-Pred}{Data}");
  plotZeeLep1EtaDiff.AddHist1D(hZeeLep1EtaDiff,"EX0",ratioColor);
  plotZeeLep1EtaDiff.SetYRange(-0.2,0.2);
  plotZeeLep1EtaDiff.AddLine(0, 0,2.4, 0,kBlack,1);
  plotZeeLep1EtaDiff.AddLine(0, 0.1,2.4, 0.1,kBlack,3);
  plotZeeLep1EtaDiff.AddLine(0,-0.1,2.4,-0.1,kBlack,3);
  plotZeeLep1EtaDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZeeLep1Eta2("zeeLep1Etalog"+norm,"","",ylabel);
  plotZeeLep1Eta2.AddHist1D(hDataLep1Eta,"data","E");
  plotZeeLep1Eta2.AddToStack(hEWKLep1Eta,"EWK",fillcolorEWK,linecolorEWK);
  plotZeeLep1Eta2.AddToStack(hTopLep1Eta,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZeeLep1Eta2.AddToStack(hZeeLep1Eta,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeLep1Eta2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLep1Eta2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLep1Eta2.AddTextBox(normtext,0.21,0.75,0.45,0.80,0,13,0.03,-1);
  plotZeeLep1Eta2.SetLogy();
  plotZeeLep1Eta2.SetYRange(1e-4*(hDataLep1Eta->GetMaximum()),500*(hDataLep1Eta->GetMaximum()));
  plotZeeLep1Eta2.TransLegend(0.1,-0.05);
  plotZeeLep1Eta2.Draw(c,kTRUE,format,1);

  //
  // Lep2 Eta
  //   
  sprintf(ylabel,"Events / %.1f ",hDataLep2Eta->GetBinWidth(1));
  CPlot plotZeeLep2Eta("zeeLep2Eta"+norm,"","",ylabel);
  plotZeeLep2Eta.AddHist1D(hDataLep2Eta,"data","E");
  plotZeeLep2Eta.AddToStack(hZeeLep2Eta,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeLep2Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLep2Eta.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLep2Eta.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZeeLep2Eta.SetLogy(0);
  plotZeeLep2Eta.SetYRange(0.01,1.2*(hDataLep2Eta->GetMaximum() + sqrt(hDataLep2Eta->GetMaximum())));
  plotZeeLep2Eta.TransLegend(0.1,-0.05);
  plotZeeLep2Eta.Draw(c,kFALSE,format,1);

  CPlot plotZeeLep2EtaDiff("zeeLep2Eta"+norm,"","|#eta| (2nd leading electron)","#frac{Data-Pred}{Data}");
  plotZeeLep2EtaDiff.AddHist1D(hZeeLep2EtaDiff,"EX0",ratioColor);
  plotZeeLep2EtaDiff.SetYRange(-0.2,0.2);
  plotZeeLep2EtaDiff.AddLine(0, 0,2.4, 0,kBlack,1);
  plotZeeLep2EtaDiff.AddLine(0, 0.1,2.4, 0.1,kBlack,3);
  plotZeeLep2EtaDiff.AddLine(0,-0.1,2.4,-0.1,kBlack,3);
  plotZeeLep2EtaDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZeeLep2Eta2("zeeLep2Etalog"+norm,"","",ylabel);
  plotZeeLep2Eta2.AddHist1D(hDataLep2Eta,"data","E");
  plotZeeLep2Eta2.AddToStack(hEWKLep2Eta,"EWK",fillcolorEWK,linecolorEWK);
  plotZeeLep2Eta2.AddToStack(hTopLep2Eta,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZeeLep2Eta2.AddToStack(hZeeLep2Eta,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZeeLep2Eta2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZeeLep2Eta2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZeeLep2Eta2.AddTextBox(normtext,0.21,0.75,0.45,0.80,0,13,0.03,-1);
  plotZeeLep2Eta2.SetLogy();
  plotZeeLep2Eta2.SetYRange(1e-4*(hDataLep2Eta->GetMaximum()),500*(hDataLep2Eta->GetMaximum()));
  plotZeeLep2Eta2.TransLegend(0.1,-0.05);
  plotZeeLep2Eta2.Draw(c,kTRUE,format,1);

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  


  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  cout << endl;

  cout << " The Zee event yield is " << yield << " +/-" << sqrt(yield) << "." << endl;
  cout << " The Zee expected event yield is " << yield_zee << " +/-" << sqrt(yield_zee_unc) << "." << endl;
  cout << " The EWK event yield is " << yield_ewk << " +/-" << sqrt(yield_ewk_unc) << "." << endl;
  cout << " The Top event yield is " << yield_top << " +/-" << sqrt(yield_top_unc) << "." << endl;

  
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
  txtfile << " The Zee expected event yield is " << yield_zee << " +/-" << sqrt(yield_zee_unc) << "." << endl;
  txtfile << " The EWK event yield is " << yield_ewk << " +/-" << sqrt(yield_ewk_unc) << "." << endl;
  txtfile << " The Top event yield is " << yield_top << " +/-" << sqrt(yield_top_unc) << "." << endl;
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

