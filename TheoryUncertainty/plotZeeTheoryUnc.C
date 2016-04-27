//================================================================================================
//
// Make plots of various distributions after Zee selection
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>   
#include <TGraphAsymmErrors.h> 
#include <TLatex.h>
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include <TGaxis.h>
#include "TLorentzVector.h"               // 4-vector class

#include "ConfParse.hh"                   // input conf file parser
#include "../Utils/CSample.hh"            // helper class to handle samples
#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/LeptonCorr.hh"


#endif

string int2string(int i) {
  stringstream ss;
  string ret;
  ss << i;
  ss >> ret;
  return ret;
}

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);
TGraphAsymmErrors* TH1TOTGraphAsymmErrors(TH1 *h1);
TGraphAsymmErrors* myMakeBandSymmetric(TGraphAsymmErrors* g0, TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);
void myAddtoBand(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2);
TGraphAsymmErrors* myTGraphErrorsDivide_noErrGraph2(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);

//=== MAIN MACRO ================================================================================================= 

void plotZeeTheoryUnc(const TString  conf,            // input file
		      const TString  inputDir,        // input directory
		      const TString  outputDir,       // output directory
		      const Double_t lumi             // integrated luminosity (/fb)
		      ) {
  gBenchmark->Start("plotZeeTheoryUnc");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  const Double_t ELE_MASS  = 0.000511;
  
  const TString format("png");
  
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.5;
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
 
  vector<TString>  snamev;      // sample name 
  vector<CSample*> samplev;     // data/MC samples

  //
  // parse .conf file
  //
  confParse(conf, snamev, samplev);
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);

  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  npv, npu;
  UInt_t  triggerDec;
  UInt_t  goodPV;
  UInt_t  matchTrigger;
  UInt_t  ngenlep;
  TLorentzVector *genlep1=0, *genlep2=0;
  Int_t   genq1, genq2;
  UInt_t nlep;
  TLorentzVector *lep1=0, *lep2=0;
  Int_t   q1, q2;
  Float_t scale1fbGen,scale1fb, scale1fbUp, scale1fbDown;
  vector<float> *lheweight = new vector<float>();

  
  TFile *infile=0;
  TTree *intree=0;

  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    
    // Read input file and get the TTrees
    TString infilename=inputDir + TString("/") + snamev[isam] + TString("_select.raw.root");
    cout << "Processing " << infilename << "..." << endl;
    infile = TFile::Open(infilename);assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);
    
    intree->SetBranchAddress("runNum",   &runNum);     // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);    // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);     // event number
    intree->SetBranchAddress("matchGen", &matchGen);   // event has both leptons matched to MC Z->ll
    intree->SetBranchAddress("npv",      &npv);        // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);        // number of in-time PU events (MC)
    intree->SetBranchAddress("triggerDec",   &triggerDec);    // event pass the trigger
    intree->SetBranchAddress("goodPV",   &goodPV);    // event has a good PV
    intree->SetBranchAddress("matchTrigger",   &matchTrigger);    // event has at least one lepton matched to the trigger
    intree->SetBranchAddress("ngenlep",     &ngenlep);      // number of gen leptons
    intree->SetBranchAddress("genlep1",   &genlep1);     // gen lepton1 4-vector
    intree->SetBranchAddress("genlep2",   &genlep2);     // gen lepton2 4-vector
    intree->SetBranchAddress("genq1",     &genq1);     // charge gen lepton1
    intree->SetBranchAddress("genq2",     &genq2);     // charge gen lepton2
    intree->SetBranchAddress("nlep",     &nlep);      // number of leptons
    intree->SetBranchAddress("lep1",       &lep1);     // lepton1 4-vector
    intree->SetBranchAddress("lep2",       &lep2);     // lepton2 4-vector
    intree->SetBranchAddress("q1",       &q1);     // charge lepton1
    intree->SetBranchAddress("q2",       &q2);     // charge lepton2
    intree->SetBranchAddress("scale1fbGen",   &scale1fbGen);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fb",   &scale1fb);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp",   &scale1fbUp);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown",   &scale1fbDown);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("lheweight",  &lheweight);
    
    //
    // Set up output file
    //
    TString outfilename = outputDir + TString("/")+ snamev[isam] + TString("_PDFUnc.root");
    TFile *outFile = new TFile(outfilename,"RECREATE");
    
    //
    // Create histograms
    //
    double ZPtBins[]={0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,400,1000};
    double PhiStarBins[]={0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.01,0.012,0.014,0.016,0.018,0.021,0.024,0.027,0.030,0.034,0.038,0.044,0.050,0.058,0.066,0.076,0.088,0.10,0.12,0.14,0.16,0.18,0.20,0.24,0.28,0.34,0.42,0.52,0.64,0.8,1.0,1.5,2,3};
    double Lep1PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
    double Lep2PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150};
    double LepNegPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
    double LepPosPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
    
    const int nBinsZPt= sizeof(ZPtBins)/sizeof(double)-1;
    TH1D *hZPtTruthNominal  = new TH1D("hZPtTruthNominal","",nBinsZPt,ZPtBins); hZPtTruthNominal->Sumw2();
    TH1D *hZPtTruthScaleUp  = new TH1D("hZPtTruthScaleUp","",nBinsZPt,ZPtBins); hZPtTruthScaleUp->Sumw2();
    TH1D *hZPtTruthScaleDown  = new TH1D("hZPtTruthScaleDown","",nBinsZPt,ZPtBins); hZPtTruthScaleDown->Sumw2();
    TH1D *hZPtTruthPDFUp  = new TH1D("hZPtTruthPDFUp","",nBinsZPt,ZPtBins); hZPtTruthPDFUp->Sumw2();
    TH1D *hZPtTruthPDFDown  = new TH1D("hZPtTruthPDFDown","",nBinsZPt,ZPtBins); hZPtTruthPDFDown->Sumw2();
    TH1D *hZPtTruthAlphasUp  = new TH1D("hZPtTruthAlphasUp","",nBinsZPt,ZPtBins); hZPtTruthAlphasUp->Sumw2();
    TH1D *hZPtTruthAlphasDown  = new TH1D("hZPtTruthAlphasDown","",nBinsZPt,ZPtBins); hZPtTruthAlphasDown->Sumw2();
    
    const int nBinsPhiStar= sizeof(PhiStarBins)/sizeof(double)-1;
    TH1D *hPhiStarTruthNominal  = new TH1D("hPhiStarTruthNominal","",nBinsPhiStar,PhiStarBins); hPhiStarTruthNominal->Sumw2();
    TH1D *hPhiStarTruthScaleUp  = new TH1D("hPhiStarTruthScaleUp","",nBinsPhiStar,PhiStarBins); hPhiStarTruthScaleUp->Sumw2();
    TH1D *hPhiStarTruthScaleDown  = new TH1D("hPhiStarTruthScaleDown","",nBinsPhiStar,PhiStarBins); hPhiStarTruthScaleDown->Sumw2();
    TH1D *hPhiStarTruthPDFUp  = new TH1D("hPhiStarTruthPDFUp","",nBinsPhiStar,PhiStarBins); hPhiStarTruthPDFUp->Sumw2();
    TH1D *hPhiStarTruthPDFDown  = new TH1D("hPhiStarTruthPDFDown","",nBinsPhiStar,PhiStarBins); hPhiStarTruthPDFDown->Sumw2();
    TH1D *hPhiStarTruthAlphasUp  = new TH1D("hPhiStarTruthAlphasUp","",nBinsPhiStar,PhiStarBins); hPhiStarTruthAlphasUp->Sumw2();
    TH1D *hPhiStarTruthAlphasDown  = new TH1D("hPhiStarTruthAlphasDown","",nBinsPhiStar,PhiStarBins); hPhiStarTruthAlphasDown->Sumw2();
    
    
    TH1D *hZRapTruthNominal  = new TH1D("hZRapTruthNominal","",24,0,2.4); hZRapTruthNominal->Sumw2();
    TH1D *hZRapTruthScaleUp  = new TH1D("hZRapTruthScaleUp","",24,0,2.4); hZRapTruthScaleUp->Sumw2();
    TH1D *hZRapTruthScaleDown  = new TH1D("hZRapTruthScaleDown","",24,0,2.4); hZRapTruthScaleDown->Sumw2();
    TH1D *hZRapTruthPDFUp  = new TH1D("hZRapTruthPDFUp","",24,0,2.4); hZRapTruthPDFUp->Sumw2();
    TH1D *hZRapTruthPDFDown  = new TH1D("hZRapTruthPDFDown","",24,0,2.4); hZRapTruthPDFDown->Sumw2();
    TH1D *hZRapTruthAlphasUp  = new TH1D("hZRapTruthAlphasUp","",24,0,2.4); hZRapTruthAlphasUp->Sumw2();
    TH1D *hZRapTruthAlphasDown  = new TH1D("hZRapTruthAlphasDown","",24,0,2.4); hZRapTruthAlphasDown->Sumw2();
    

    const int nBinsLep1Pt= sizeof(Lep1PtBins)/sizeof(double)-1;
    TH1D *hLep1PtTruthNominal  = new TH1D("hLep1PtTruthNominal","",nBinsLep1Pt,Lep1PtBins); hLep1PtTruthNominal->Sumw2();
    TH1D *hLep1PtTruthScaleUp  = new TH1D("hLep1PtTruthScaleUp","",nBinsLep1Pt,Lep1PtBins); hLep1PtTruthScaleUp->Sumw2();
    TH1D *hLep1PtTruthScaleDown  = new TH1D("hLep1PtTruthScaleDown","",nBinsLep1Pt,Lep1PtBins); hLep1PtTruthScaleDown->Sumw2();
    TH1D *hLep1PtTruthPDFUp  = new TH1D("hLep1PtTruthPDFUp","",nBinsLep1Pt,Lep1PtBins); hLep1PtTruthPDFUp->Sumw2();
    TH1D *hLep1PtTruthPDFDown  = new TH1D("hLep1PtTruthPDFDown","",nBinsLep1Pt,Lep1PtBins); hLep1PtTruthPDFDown->Sumw2();
    TH1D *hLep1PtTruthAlphasUp  = new TH1D("hLep1PtTruthAlphasUp","",nBinsLep1Pt,Lep1PtBins); hLep1PtTruthAlphasUp->Sumw2();
    TH1D *hLep1PtTruthAlphasDown  = new TH1D("hLep1PtTruthAlphasDown","",nBinsLep1Pt,Lep1PtBins); hLep1PtTruthAlphasDown->Sumw2();
    
    
    const int nBinsLep2Pt= sizeof(Lep2PtBins)/sizeof(double)-1;
    TH1D *hLep2PtTruthNominal  = new TH1D("hLep2PtTruthNominal","",nBinsLep2Pt,Lep2PtBins); hLep2PtTruthNominal->Sumw2();
    TH1D *hLep2PtTruthScaleUp  = new TH1D("hLep2PtTruthScaleUp","",nBinsLep2Pt,Lep2PtBins); hLep2PtTruthScaleUp->Sumw2();
    TH1D *hLep2PtTruthScaleDown  = new TH1D("hLep2PtTruthScaleDown","",nBinsLep2Pt,Lep2PtBins); hLep2PtTruthScaleDown->Sumw2();
    TH1D *hLep2PtTruthPDFUp  = new TH1D("hLep2PtTruthPDFUp","",nBinsLep2Pt,Lep2PtBins); hLep2PtTruthPDFUp->Sumw2();
    TH1D *hLep2PtTruthPDFDown  = new TH1D("hLep2PtTruthPDFDown","",nBinsLep2Pt,Lep2PtBins); hLep2PtTruthPDFDown->Sumw2();
    TH1D *hLep2PtTruthAlphasUp  = new TH1D("hLep2PtTruthAlphasUp","",nBinsLep2Pt,Lep2PtBins); hLep2PtTruthAlphasUp->Sumw2();
    TH1D *hLep2PtTruthAlphasDown  = new TH1D("hLep2PtTruthAlphasDown","",nBinsLep2Pt,Lep2PtBins); hLep2PtTruthAlphasDown->Sumw2();
    
    
    const int nBinsLepNegPt= sizeof(LepNegPtBins)/sizeof(double)-1;
    TH1D *hLepNegPtTruthNominal  = new TH1D("hLepNegPtTruthNominal","",nBinsLepNegPt,LepNegPtBins); hLepNegPtTruthNominal->Sumw2();
    TH1D *hLepNegPtTruthScaleUp  = new TH1D("hLepNegPtTruthScaleUp","",nBinsLepNegPt,LepNegPtBins); hLepNegPtTruthScaleUp->Sumw2();
    TH1D *hLepNegPtTruthScaleDown  = new TH1D("hLepNegPtTruthScaleDown","",nBinsLepNegPt,LepNegPtBins); hLepNegPtTruthScaleDown->Sumw2();
    TH1D *hLepNegPtTruthPDFUp  = new TH1D("hLepNegPtTruthPDFUp","",nBinsLepNegPt,LepNegPtBins); hLepNegPtTruthPDFUp->Sumw2();
    TH1D *hLepNegPtTruthPDFDown  = new TH1D("hLepNegPtTruthPDFDown","",nBinsLepNegPt,LepNegPtBins); hLepNegPtTruthPDFDown->Sumw2();
    TH1D *hLepNegPtTruthAlphasUp  = new TH1D("hLepNegPtTruthAlphasUp","",nBinsLepNegPt,LepNegPtBins); hLepNegPtTruthAlphasUp->Sumw2();
    TH1D *hLepNegPtTruthAlphasDown  = new TH1D("hLepNegPtTruthAlphasDown","",nBinsLepNegPt,LepNegPtBins); hLepNegPtTruthAlphasDown->Sumw2();
    
    
    const int nBinsLepPosPt= sizeof(LepPosPtBins)/sizeof(double)-1;
    TH1D *hLepPosPtTruthNominal  = new TH1D("hLepPosPtTruthNominal","",nBinsLepPosPt,LepPosPtBins); hLepPosPtTruthNominal->Sumw2();
    TH1D *hLepPosPtTruthScaleUp  = new TH1D("hLepPosPtTruthScaleUp","",nBinsLepPosPt,LepPosPtBins); hLepPosPtTruthScaleUp->Sumw2();
    TH1D *hLepPosPtTruthScaleDown  = new TH1D("hLepPosPtTruthScaleDown","",nBinsLepPosPt,LepPosPtBins); hLepPosPtTruthScaleDown->Sumw2();
    TH1D *hLepPosPtTruthPDFUp  = new TH1D("hLepPosPtTruthPDFUp","",nBinsLepPosPt,LepPosPtBins); hLepPosPtTruthPDFUp->Sumw2();
    TH1D *hLepPosPtTruthPDFDown  = new TH1D("hLepPosPtTruthPDFDown","",nBinsLepPosPt,LepPosPtBins); hLepPosPtTruthPDFDown->Sumw2();
    TH1D *hLepPosPtTruthAlphasUp  = new TH1D("hLepPosPtTruthAlphasUp","",nBinsLepPosPt,LepPosPtBins); hLepPosPtTruthAlphasUp->Sumw2();
    TH1D *hLepPosPtTruthAlphasDown  = new TH1D("hLepPosPtTruthAlphasDown","",nBinsLepPosPt,LepPosPtBins); hLepPosPtTruthAlphasDown->Sumw2();
    
    TH1D *hLep1EtaTruthNominal  = new TH1D("hLep1EtaTruthNominal","",24,0,2.4); hLep1EtaTruthNominal->Sumw2();
    TH1D *hLep1EtaTruthScaleUp  = new TH1D("hLep1EtaTruthScaleUp","",24,0,2.4); hLep1EtaTruthScaleUp->Sumw2();
    TH1D *hLep1EtaTruthScaleDown  = new TH1D("hLep1EtaTruthScaleDown","",24,0,2.4); hLep1EtaTruthScaleDown->Sumw2();
    TH1D *hLep1EtaTruthPDFUp  = new TH1D("hLep1EtaTruthPDFUp","",24,0,2.4); hLep1EtaTruthPDFUp->Sumw2();
    TH1D *hLep1EtaTruthPDFDown  = new TH1D("hLep1EtaTruthPDFDown","",24,0,2.4); hLep1EtaTruthPDFDown->Sumw2();
    TH1D *hLep1EtaTruthAlphasUp  = new TH1D("hLep1EtaTruthAlphasUp","",24,0,2.4); hLep1EtaTruthAlphasUp->Sumw2();
    TH1D *hLep1EtaTruthAlphasDown  = new TH1D("hLep1EtaTruthAlphasDown","",24,0,2.4); hLep1EtaTruthAlphasDown->Sumw2();
    
    TH1D *hLep2EtaTruthNominal  = new TH1D("hLep2EtaTruthNominal","",24,0,2.4); hLep2EtaTruthNominal->Sumw2();
    TH1D *hLep2EtaTruthScaleUp  = new TH1D("hLep2EtaTruthScaleUp","",24,0,2.4); hLep2EtaTruthScaleUp->Sumw2();
    TH1D *hLep2EtaTruthScaleDown  = new TH1D("hLep2EtaTruthScaleDown","",24,0,2.4); hLep2EtaTruthScaleDown->Sumw2();
    TH1D *hLep2EtaTruthPDFUp  = new TH1D("hLep2EtaTruthPDFUp","",24,0,2.4); hLep2EtaTruthPDFUp->Sumw2();
    TH1D *hLep2EtaTruthPDFDown  = new TH1D("hLep2EtaTruthPDFDown","",24,0,2.4); hLep2EtaTruthPDFDown->Sumw2();
    TH1D *hLep2EtaTruthAlphasUp  = new TH1D("hLep2EtaTruthAlphasUp","",24,0,2.4); hLep2EtaTruthAlphasUp->Sumw2();
    TH1D *hLep2EtaTruthAlphasDown  = new TH1D("hLep2EtaTruthAlphasDown","",24,0,2.4); hLep2EtaTruthAlphasDown->Sumw2();
    
    TH1D *hZPtTruth[111];
    TH1D *hPhiStarTruth[111];
    TH1D *hZRapTruth[111];
    TH1D *hLep1PtTruth[111];
    TH1D *hLep2PtTruth[111];
    TH1D *hLepNegPtTruth[111];
    TH1D *hLepPosPtTruth[111];
    TH1D *hLep1EtaTruth[111];
    TH1D *hLep2EtaTruth[111];
    for(int i=0;i!=111;i++)
      {
	string si=int2string(i);
	hZPtTruth[i]  = new TH1D((string("hZPtTruth")+si).c_str(),"",nBinsZPt,ZPtBins); hZPtTruth[i]->Sumw2();
	hPhiStarTruth[i]  = new TH1D((string("hPhiStarTruth")+si).c_str(),"",nBinsPhiStar,PhiStarBins); hPhiStarTruth[i]->Sumw2();
	hZRapTruth[i]  = new TH1D((string("hZRapTruth")+si).c_str(),"",24,0,2.4); hZRapTruth[i]->Sumw2();
	hLep1PtTruth[i]  = new TH1D((string("hLep1PtTruth")+si).c_str(),"",nBinsLep1Pt,Lep1PtBins); hLep1PtTruth[i]->Sumw2();
	hLep2PtTruth[i]  = new TH1D((string("hLep2PtTruth")+si).c_str(),"",nBinsLep2Pt,Lep2PtBins); hLep2PtTruth[i]->Sumw2();
	hLepNegPtTruth[i]  = new TH1D((string("hLepNegPtTruth")+si).c_str(),"",nBinsLepNegPt,LepNegPtBins); hLepNegPtTruth[i]->Sumw2();
	hLepPosPtTruth[i]  = new TH1D((string("hLepPosPtTruth")+si).c_str(),"",nBinsLepPosPt,LepPosPtBins); hLepPosPtTruth[i]->Sumw2();
	hLep1EtaTruth[i]  = new TH1D((string("hLep1EtaTruth")+si).c_str(),"",24,0,2.4); hLep1EtaTruth[i]->Sumw2();
	hLep2EtaTruth[i]  = new TH1D((string("hLep2EtaTruth")+si).c_str(),"",24,0,2.4); hLep2EtaTruth[i]->Sumw2();
      }
    
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      
      TLorentzVector *gendilep=new TLorentzVector(0,0,0,0);
      gendilep->operator+=(*genlep1);
      gendilep->operator+=(*genlep2);
      
    double genphiacop=0;
    double gencosthetastar=0;
    double genphistar=0;

    genphiacop=TMath::Pi()-fabs(genlep1->DeltaPhi(*genlep2));
    if(genq1<0) gencosthetastar=TMath::TanH((genlep1->Rapidity()-genlep2->Rapidity())/2);
    else gencosthetastar=TMath::TanH((genlep2->Rapidity()-genlep1->Rapidity())/2);
    genphistar=TMath::Tan(genphiacop/2)*sqrt(1-pow(gencosthetastar,2));
    
    bool isGen=false;
   
    if(ngenlep>=2&&gendilep->M()>MASS_LOW&&gendilep->M()<MASS_HIGH&&genlep1->Pt()>=PT_CUT&&genlep2->Pt()>=PT_CUT&&fabs(genlep1->Eta())<=ETA_CUT&&fabs(genlep2->Eta())<=ETA_CUT)
      {
	isGen=true;
      }

    Double_t genweight[111];
    for(int i=0;i!=111;i++)
      {
	genweight[i]=1;
	genweight[i]*= scale1fbGen*lumi;
	if(i>0)
	  {
	    genweight[i]*=(*lheweight)[i];
	  }
      }
  
    if(isGen)
      {
	for(int i=0;i!=111;i++)
	  {
	    hZPtTruth[i] ->Fill(gendilep->Pt(),genweight[i]);
	    hPhiStarTruth[i] ->Fill(genphistar,genweight[i]);
	    hZRapTruth[i] ->Fill(fabs(gendilep->Rapidity()),genweight[i]);
	    hLep1PtTruth[i] ->Fill(genlep1->Pt(),genweight[i]);
	    hLep2PtTruth[i] ->Fill(genlep2->Pt(),genweight[i]);
	    if(genq1<0)
	      {
		hLepNegPtTruth[i] ->Fill(genlep1->Pt(),genweight[i]);
		hLepPosPtTruth[i] ->Fill(genlep2->Pt(),genweight[i]);
	      }
	    else
	      {
		hLepNegPtTruth[i] ->Fill(genlep2->Pt(),genweight[i]);
		hLepPosPtTruth[i] ->Fill(genlep1->Pt(),genweight[i]);
	      }
	    hLep1EtaTruth[i] ->Fill(fabs(genlep1->Eta()),genweight[i]);
	    hLep2EtaTruth[i] ->Fill(fabs(genlep2->Eta()),genweight[i]);
	  }
      }
  }
  delete infile;
  infile=0, intree=0; 

  //---------------------------------------------------------------------------
  //                             Z Pt
  //---------------------------------------------------------------------------

  double pdfUncZPtUp[34];
  double pdfUncZPtDown[34];
  double alphasUncZPtUp[34];
  double alphasUncZPtDown[34];
  double scaleUncZPtUp[34];
  double scaleUncZPtDown[34];

  int npdfZPtP=0;
  int npdfZPtM=0;

  for(int i=0;i!=hZPtTruth[0]->GetNbinsX();i++)
    {
      scaleUncZPtUp[i]=0;
      scaleUncZPtDown[i]=0;
      pdfUncZPtUp[i]=0;
      pdfUncZPtDown[i]=0;
      for(int j=1;j!=109;j++)
	{
	  if(j<=8)
	    {
	      if(hZPtTruth[0]->GetBinContent(i+1)<hZPtTruth[j]->GetBinContent(i+1))
		{
		  if((hZPtTruth[j]->GetBinContent(i+1)-hZPtTruth[0]->GetBinContent(i+1))>scaleUncZPtUp[i])
		    {
		      scaleUncZPtUp[i]=hZPtTruth[j]->GetBinContent(i+1)-hZPtTruth[0]->GetBinContent(i+1);
		    }
		}
	      else
		{
		  if((hZPtTruth[0]->GetBinContent(i+1)-hZPtTruth[j]->GetBinContent(i+1))>scaleUncZPtDown[i])
		    {
		      scaleUncZPtDown[i]=hZPtTruth[0]->GetBinContent(i+1)-hZPtTruth[j]->GetBinContent(i+1);
		    }
		}
	    }
	  else
	    {
	      if(hZPtTruth[0]->GetBinContent(i+1)<hZPtTruth[j]->GetBinContent(i+1))
		{
		  pdfUncZPtUp[i]+=pow(hZPtTruth[0]->GetBinContent(i+1)-hZPtTruth[j]->GetBinContent(i+1),2);
		  npdfZPtP++;
		}
	      else
		{
		  pdfUncZPtDown[i]+=pow(hZPtTruth[j]->GetBinContent(i+1)-hZPtTruth[0]->GetBinContent(i+1),2);
		  npdfZPtM++;
		}
	    }
	}
      alphasUncZPtUp[i]=fabs(hZPtTruth[110]->GetBinContent(i+1)-hZPtTruth[109]->GetBinContent(i+1))/2;
      alphasUncZPtDown[i]=fabs(hZPtTruth[110]->GetBinContent(i+1)-hZPtTruth[109]->GetBinContent(i+1))/2;

      hZPtTruthNominal->SetBinContent(i+1,hZPtTruth[0]->GetBinContent(i+1));
      hZPtTruthNominal->SetBinError(i+1,hZPtTruth[0]->GetBinError(i+1));

      hZPtTruthScaleUp->SetBinContent(i+1,hZPtTruth[0]->GetBinContent(i+1)+scaleUncZPtUp[i]);
      hZPtTruthScaleDown->SetBinContent(i+1,hZPtTruth[0]->GetBinContent(i+1)-scaleUncZPtDown[i]);

      hZPtTruthPDFUp->SetBinContent(i+1,hZPtTruth[0]->GetBinContent(i+1)+sqrt(pdfUncZPtUp[i]/(npdfZPtP-1)));
      hZPtTruthPDFDown->SetBinContent(i+1,hZPtTruth[0]->GetBinContent(i+1)-sqrt(pdfUncZPtDown[i]/(npdfZPtM-1)));

      hZPtTruthAlphasUp->SetBinContent(i+1,hZPtTruth[0]->GetBinContent(i+1)+alphasUncZPtUp[i]);
      hZPtTruthAlphasDown->SetBinContent(i+1,hZPtTruth[0]->GetBinContent(i+1)-alphasUncZPtDown[i]);      
    }
  for(int j=0;j!=hZPtTruth[0]->GetNbinsX();++j)
    {
      hZPtTruthNominal->SetBinContent(j+1,hZPtTruthNominal->GetBinContent(j+1)/hZPtTruthNominal->GetBinWidth(j+1));
      hZPtTruthNominal->SetBinError(j+1,hZPtTruthNominal->GetBinError(j+1)/hZPtTruthNominal->GetBinWidth(j+1));

      hZPtTruthScaleUp->SetBinContent(j+1,hZPtTruthScaleUp->GetBinContent(j+1)/hZPtTruthScaleUp->GetBinWidth(j+1));
      hZPtTruthScaleUp->SetBinError(j+1,hZPtTruthScaleUp->GetBinError(j+1)/hZPtTruthScaleUp->GetBinWidth(j+1));
      hZPtTruthScaleDown->SetBinContent(j+1,hZPtTruthScaleDown->GetBinContent(j+1)/hZPtTruthScaleDown->GetBinWidth(j+1));
      hZPtTruthScaleDown->SetBinError(j+1,hZPtTruthScaleDown->GetBinError(j+1)/hZPtTruthScaleDown->GetBinWidth(j+1));

      hZPtTruthPDFUp->SetBinContent(j+1,hZPtTruthPDFUp->GetBinContent(j+1)/hZPtTruthPDFUp->GetBinWidth(j+1));
      hZPtTruthPDFUp->SetBinError(j+1,hZPtTruthPDFUp->GetBinError(j+1)/hZPtTruthPDFUp->GetBinWidth(j+1));
      hZPtTruthPDFDown->SetBinContent(j+1,hZPtTruthPDFDown->GetBinContent(j+1)/hZPtTruthPDFDown->GetBinWidth(j+1));
      hZPtTruthPDFDown->SetBinError(j+1,hZPtTruthPDFDown->GetBinError(j+1)/hZPtTruthPDFDown->GetBinWidth(j+1));

      hZPtTruthAlphasUp->SetBinContent(j+1,hZPtTruthAlphasUp->GetBinContent(j+1)/hZPtTruthAlphasUp->GetBinWidth(j+1));
      hZPtTruthAlphasUp->SetBinError(j+1,hZPtTruthAlphasUp->GetBinError(j+1)/hZPtTruthAlphasUp->GetBinWidth(j+1));
      hZPtTruthAlphasDown->SetBinContent(j+1,hZPtTruthAlphasDown->GetBinContent(j+1)/hZPtTruthAlphasDown->GetBinWidth(j+1));
      hZPtTruthAlphasDown->SetBinError(j+1,hZPtTruthAlphasDown->GetBinError(j+1)/hZPtTruthAlphasDown->GetBinWidth(j+1));
    }
  hZPtTruthNominal->Scale(1./lumi);
  hZPtTruthScaleUp->Scale(1./lumi);
  hZPtTruthScaleDown->Scale(1./lumi);
  hZPtTruthPDFUp->Scale(1./lumi);
  hZPtTruthPDFDown->Scale(1./lumi);
  hZPtTruthAlphasUp->Scale(1./lumi);
  hZPtTruthAlphasDown->Scale(1./lumi);
  

  //---------------------------------------------------------------------------
  //                             PhiStar
  //---------------------------------------------------------------------------

  double pdfUncPhiStarUp[27];
  double pdfUncPhiStarDown[27];
  double alphasUncPhiStarUp[27];
  double alphasUncPhiStarDown[27];
  double scaleUncPhiStarUp[27];
  double scaleUncPhiStarDown[27];

  int npdfPhiStarP=0;
  int npdfPhiStarM=0;

  
  for(int i=0;i!=hPhiStarTruth[0]->GetNbinsX();i++)
    {
      scaleUncPhiStarUp[i]=0;
      scaleUncPhiStarDown[i]=0;
      pdfUncPhiStarUp[i]=0;
      pdfUncPhiStarDown[i]=0;
      for(int j=1;j!=109;j++)
	{
	  if(j<=8)
	    {
	      if(hPhiStarTruth[0]->GetBinContent(i+1)<hPhiStarTruth[j]->GetBinContent(i+1))
		{
		  if((hPhiStarTruth[j]->GetBinContent(i+1)-hPhiStarTruth[0]->GetBinContent(i+1))>scaleUncPhiStarUp[i])
		    {
		      scaleUncPhiStarUp[i]=hPhiStarTruth[j]->GetBinContent(i+1)-hPhiStarTruth[0]->GetBinContent(i+1);
		    }
		}
	      else
		{
		  if((hPhiStarTruth[0]->GetBinContent(i+1)-hPhiStarTruth[j]->GetBinContent(i+1))>scaleUncPhiStarDown[i])
		    {
		      scaleUncPhiStarDown[i]=hPhiStarTruth[0]->GetBinContent(i+1)-hPhiStarTruth[j]->GetBinContent(i+1);
		    }
		}
	    }
	  else
	    {
	      if(hPhiStarTruth[0]->GetBinContent(i+1)<hPhiStarTruth[j]->GetBinContent(i+1))
		{
		  pdfUncPhiStarUp[i]+=pow(hPhiStarTruth[0]->GetBinContent(i+1)-hPhiStarTruth[j]->GetBinContent(i+1),2);
		  npdfPhiStarP++;
		}
	      else
		{
		  pdfUncPhiStarDown[i]+=pow(hPhiStarTruth[j]->GetBinContent(i+1)-hPhiStarTruth[0]->GetBinContent(i+1),2);
		  npdfPhiStarM++;
		}
	    }
	}
      alphasUncPhiStarUp[i]=fabs(hPhiStarTruth[110]->GetBinContent(i+1)-hPhiStarTruth[109]->GetBinContent(i+1))/2;
      alphasUncPhiStarDown[i]=fabs(hPhiStarTruth[110]->GetBinContent(i+1)-hPhiStarTruth[109]->GetBinContent(i+1))/2;

      hPhiStarTruthNominal->SetBinContent(i+1,hPhiStarTruth[0]->GetBinContent(i+1));
      hPhiStarTruthNominal->SetBinError(i+1,hPhiStarTruth[0]->GetBinError(i+1));

      hPhiStarTruthScaleUp->SetBinContent(i+1,hPhiStarTruth[0]->GetBinContent(i+1)+scaleUncPhiStarUp[i]);
      hPhiStarTruthScaleDown->SetBinContent(i+1,hPhiStarTruth[0]->GetBinContent(i+1)-scaleUncPhiStarDown[i]);

      hPhiStarTruthPDFUp->SetBinContent(i+1,hPhiStarTruth[0]->GetBinContent(i+1)+sqrt(pdfUncPhiStarUp[i]/(npdfPhiStarP-1)));
      hPhiStarTruthPDFDown->SetBinContent(i+1,hPhiStarTruth[0]->GetBinContent(i+1)-sqrt(pdfUncPhiStarDown[i]/(npdfPhiStarM-1)));

      hPhiStarTruthAlphasUp->SetBinContent(i+1,hPhiStarTruth[0]->GetBinContent(i+1)+alphasUncPhiStarUp[i]);
      hPhiStarTruthAlphasDown->SetBinContent(i+1,hPhiStarTruth[0]->GetBinContent(i+1)-alphasUncPhiStarDown[i]);      
    }
  for(int j=0;j!=hPhiStarTruth[0]->GetNbinsX();++j)
    {
      hPhiStarTruthNominal->SetBinContent(j+1,hPhiStarTruthNominal->GetBinContent(j+1)/hPhiStarTruthNominal->GetBinWidth(j+1));
      hPhiStarTruthNominal->SetBinError(j+1,hPhiStarTruthNominal->GetBinError(j+1)/hPhiStarTruthNominal->GetBinWidth(j+1));

      hPhiStarTruthScaleUp->SetBinContent(j+1,hPhiStarTruthScaleUp->GetBinContent(j+1)/hPhiStarTruthScaleUp->GetBinWidth(j+1));
      hPhiStarTruthScaleUp->SetBinError(j+1,hPhiStarTruthScaleUp->GetBinError(j+1)/hPhiStarTruthScaleUp->GetBinWidth(j+1));
      hPhiStarTruthScaleDown->SetBinContent(j+1,hPhiStarTruthScaleDown->GetBinContent(j+1)/hPhiStarTruthScaleDown->GetBinWidth(j+1));
      hPhiStarTruthScaleDown->SetBinError(j+1,hPhiStarTruthScaleDown->GetBinError(j+1)/hPhiStarTruthScaleDown->GetBinWidth(j+1));

      hPhiStarTruthPDFUp->SetBinContent(j+1,hPhiStarTruthPDFUp->GetBinContent(j+1)/hPhiStarTruthPDFUp->GetBinWidth(j+1));
      hPhiStarTruthPDFUp->SetBinError(j+1,hPhiStarTruthPDFUp->GetBinError(j+1)/hPhiStarTruthPDFUp->GetBinWidth(j+1));
      hPhiStarTruthPDFDown->SetBinContent(j+1,hPhiStarTruthPDFDown->GetBinContent(j+1)/hPhiStarTruthPDFDown->GetBinWidth(j+1));
      hPhiStarTruthPDFDown->SetBinError(j+1,hPhiStarTruthPDFDown->GetBinError(j+1)/hPhiStarTruthPDFDown->GetBinWidth(j+1));

      hPhiStarTruthAlphasUp->SetBinContent(j+1,hPhiStarTruthAlphasUp->GetBinContent(j+1)/hPhiStarTruthAlphasUp->GetBinWidth(j+1));
      hPhiStarTruthAlphasUp->SetBinError(j+1,hPhiStarTruthAlphasUp->GetBinError(j+1)/hPhiStarTruthAlphasUp->GetBinWidth(j+1));
      hPhiStarTruthAlphasDown->SetBinContent(j+1,hPhiStarTruthAlphasDown->GetBinContent(j+1)/hPhiStarTruthAlphasDown->GetBinWidth(j+1));
      hPhiStarTruthAlphasDown->SetBinError(j+1,hPhiStarTruthAlphasDown->GetBinError(j+1)/hPhiStarTruthAlphasDown->GetBinWidth(j+1));
    }
  hPhiStarTruthNominal->Scale(1./lumi);
  hPhiStarTruthScaleUp->Scale(1./lumi);
  hPhiStarTruthScaleDown->Scale(1./lumi);
  hPhiStarTruthPDFUp->Scale(1./lumi);
  hPhiStarTruthPDFDown->Scale(1./lumi);
  hPhiStarTruthAlphasUp->Scale(1./lumi);
  hPhiStarTruthAlphasDown->Scale(1./lumi);

  //---------------------------------------------------------------------------
  //                             Z Rapidity
  //---------------------------------------------------------------------------

  double pdfUncZRapUp[24];
  double pdfUncZRapDown[24];
  double alphasUncZRapUp[24];
  double alphasUncZRapDown[24];
  double scaleUncZRapUp[24];
  double scaleUncZRapDown[24];

  int npdfZRapP=0;
  int npdfZRapM=0;

  for(int i=0;i!=hZRapTruth[0]->GetNbinsX();i++)
    {
      scaleUncZRapUp[i]=0;
      scaleUncZRapDown[i]=0;
      pdfUncZRapUp[i]=0;
      pdfUncZRapDown[i]=0;
      for(int j=1;j!=109;j++)
	{
	  if(j<=8)
	    {
	      if(hZRapTruth[0]->GetBinContent(i+1)<hZRapTruth[j]->GetBinContent(i+1))
		{
		  if((hZRapTruth[j]->GetBinContent(i+1)-hZRapTruth[0]->GetBinContent(i+1))>scaleUncZRapUp[i])
		    {
		      scaleUncZRapUp[i]=hZRapTruth[j]->GetBinContent(i+1)-hZRapTruth[0]->GetBinContent(i+1);
		    }
		}
	      else
		{
		  if((hZRapTruth[0]->GetBinContent(i+1)-hZRapTruth[j]->GetBinContent(i+1))>scaleUncZRapDown[i])
		    {
		      scaleUncZRapDown[i]=hZRapTruth[0]->GetBinContent(i+1)-hZRapTruth[j]->GetBinContent(i+1);
		    }
		}
	    }
	  else
	    {
	      if(hZRapTruth[0]->GetBinContent(i+1)<hZRapTruth[j]->GetBinContent(i+1))
		{
		  pdfUncZRapUp[i]+=pow(hZRapTruth[0]->GetBinContent(i+1)-hZRapTruth[j]->GetBinContent(i+1),2);
		  npdfZRapP++;
		}
	      else
		{
		  pdfUncZRapDown[i]+=pow(hZRapTruth[j]->GetBinContent(i+1)-hZRapTruth[0]->GetBinContent(i+1),2);
		  npdfZRapM++;
		}
	    }
	}
      alphasUncZRapUp[i]=fabs(hZRapTruth[110]->GetBinContent(i+1)-hZRapTruth[109]->GetBinContent(i+1))/2;
      alphasUncZRapDown[i]=fabs(hZRapTruth[110]->GetBinContent(i+1)-hZRapTruth[109]->GetBinContent(i+1))/2;

      hZRapTruthNominal->SetBinContent(i+1,hZRapTruth[0]->GetBinContent(i+1));
      hZRapTruthNominal->SetBinError(i+1,hZRapTruth[0]->GetBinError(i+1));

      hZRapTruthScaleUp->SetBinContent(i+1,hZRapTruth[0]->GetBinContent(i+1)+scaleUncZRapUp[i]);
      hZRapTruthScaleDown->SetBinContent(i+1,hZRapTruth[0]->GetBinContent(i+1)-scaleUncZRapDown[i]);

      hZRapTruthPDFUp->SetBinContent(i+1,hZRapTruth[0]->GetBinContent(i+1)+sqrt(pdfUncZRapUp[i]/(npdfZRapP-1)));
      hZRapTruthPDFDown->SetBinContent(i+1,hZRapTruth[0]->GetBinContent(i+1)-sqrt(pdfUncZRapDown[i]/(npdfZRapM-1)));

      hZRapTruthAlphasUp->SetBinContent(i+1,hZRapTruth[0]->GetBinContent(i+1)+alphasUncZRapUp[i]);
      hZRapTruthAlphasDown->SetBinContent(i+1,hZRapTruth[0]->GetBinContent(i+1)-alphasUncZRapDown[i]);      
    }
  for(int j=0;j!=hZRapTruth[0]->GetNbinsX();++j)
    {
      hZRapTruthNominal->SetBinContent(j+1,hZRapTruthNominal->GetBinContent(j+1)/hZRapTruthNominal->GetBinWidth(j+1));
      hZRapTruthNominal->SetBinError(j+1,hZRapTruthNominal->GetBinError(j+1)/hZRapTruthNominal->GetBinWidth(j+1));

      hZRapTruthScaleUp->SetBinContent(j+1,hZRapTruthScaleUp->GetBinContent(j+1)/hZRapTruthScaleUp->GetBinWidth(j+1));
      hZRapTruthScaleUp->SetBinError(j+1,hZRapTruthScaleUp->GetBinError(j+1)/hZRapTruthScaleUp->GetBinWidth(j+1));
      hZRapTruthScaleDown->SetBinContent(j+1,hZRapTruthScaleDown->GetBinContent(j+1)/hZRapTruthScaleDown->GetBinWidth(j+1));
      hZRapTruthScaleDown->SetBinError(j+1,hZRapTruthScaleDown->GetBinError(j+1)/hZRapTruthScaleDown->GetBinWidth(j+1));

      hZRapTruthPDFUp->SetBinContent(j+1,hZRapTruthPDFUp->GetBinContent(j+1)/hZRapTruthPDFUp->GetBinWidth(j+1));
      hZRapTruthPDFUp->SetBinError(j+1,hZRapTruthPDFUp->GetBinError(j+1)/hZRapTruthPDFUp->GetBinWidth(j+1));
      hZRapTruthPDFDown->SetBinContent(j+1,hZRapTruthPDFDown->GetBinContent(j+1)/hZRapTruthPDFDown->GetBinWidth(j+1));
      hZRapTruthPDFDown->SetBinError(j+1,hZRapTruthPDFDown->GetBinError(j+1)/hZRapTruthPDFDown->GetBinWidth(j+1));

      hZRapTruthAlphasUp->SetBinContent(j+1,hZRapTruthAlphasUp->GetBinContent(j+1)/hZRapTruthAlphasUp->GetBinWidth(j+1));
      hZRapTruthAlphasUp->SetBinError(j+1,hZRapTruthAlphasUp->GetBinError(j+1)/hZRapTruthAlphasUp->GetBinWidth(j+1));
      hZRapTruthAlphasDown->SetBinContent(j+1,hZRapTruthAlphasDown->GetBinContent(j+1)/hZRapTruthAlphasDown->GetBinWidth(j+1));
      hZRapTruthAlphasDown->SetBinError(j+1,hZRapTruthAlphasDown->GetBinError(j+1)/hZRapTruthAlphasDown->GetBinWidth(j+1));
    }
  hZRapTruthNominal->Scale(1./lumi);
  hZRapTruthScaleUp->Scale(1./lumi);
  hZRapTruthScaleDown->Scale(1./lumi);
  hZRapTruthPDFUp->Scale(1./lumi);
  hZRapTruthPDFDown->Scale(1./lumi);
  hZRapTruthAlphasUp->Scale(1./lumi);
  hZRapTruthAlphasDown->Scale(1./lumi);

  //---------------------------------------------------------------------------
  //                             Lep1 Pt
  //---------------------------------------------------------------------------

  double pdfUncLep1PtUp[25];
  double pdfUncLep1PtDown[25];
  double alphasUncLep1PtUp[25];
  double alphasUncLep1PtDown[25];
  double scaleUncLep1PtUp[25];
  double scaleUncLep1PtDown[25];

  int npdfLep1PtP=0;
  int npdfLep1PtM=0;

  for(int i=0;i!=hLep1PtTruth[0]->GetNbinsX();i++)
    {
      scaleUncLep1PtUp[i]=0;
      scaleUncLep1PtDown[i]=0;
      pdfUncLep1PtUp[i]=0;
      pdfUncLep1PtDown[i]=0;
      for(int j=1;j!=109;j++)
	{
	  if(j<=8)
	    {
	      if(hLep1PtTruth[0]->GetBinContent(i+1)<hLep1PtTruth[j]->GetBinContent(i+1))
		{
		  if((hLep1PtTruth[j]->GetBinContent(i+1)-hLep1PtTruth[0]->GetBinContent(i+1))>scaleUncLep1PtUp[i])
		    {
		      scaleUncLep1PtUp[i]=hLep1PtTruth[j]->GetBinContent(i+1)-hLep1PtTruth[0]->GetBinContent(i+1);
		    }
		}
	      else
		{
		  if((hLep1PtTruth[0]->GetBinContent(i+1)-hLep1PtTruth[j]->GetBinContent(i+1))>scaleUncLep1PtDown[i])
		    {
		      scaleUncLep1PtDown[i]=hLep1PtTruth[0]->GetBinContent(i+1)-hLep1PtTruth[j]->GetBinContent(i+1);
		    }
		}
	    }
	  else
	    {
	      if(hLep1PtTruth[0]->GetBinContent(i+1)<hLep1PtTruth[j]->GetBinContent(i+1))
		{
		  pdfUncLep1PtUp[i]+=pow(hLep1PtTruth[0]->GetBinContent(i+1)-hLep1PtTruth[j]->GetBinContent(i+1),2);
		  npdfLep1PtP++;
		}
	      else
		{
		  pdfUncLep1PtDown[i]+=pow(hLep1PtTruth[j]->GetBinContent(i+1)-hLep1PtTruth[0]->GetBinContent(i+1),2);
		  npdfLep1PtM++;
		}
	    }
	}
      alphasUncLep1PtUp[i]=fabs(hLep1PtTruth[110]->GetBinContent(i+1)-hLep1PtTruth[109]->GetBinContent(i+1))/2;
      alphasUncLep1PtDown[i]=fabs(hLep1PtTruth[110]->GetBinContent(i+1)-hLep1PtTruth[109]->GetBinContent(i+1))/2;

      hLep1PtTruthNominal->SetBinContent(i+1,hLep1PtTruth[0]->GetBinContent(i+1));
      hLep1PtTruthNominal->SetBinError(i+1,hLep1PtTruth[0]->GetBinError(i+1));

      hLep1PtTruthScaleUp->SetBinContent(i+1,hLep1PtTruth[0]->GetBinContent(i+1)+scaleUncLep1PtUp[i]);
      hLep1PtTruthScaleDown->SetBinContent(i+1,hLep1PtTruth[0]->GetBinContent(i+1)-scaleUncLep1PtDown[i]);

      hLep1PtTruthPDFUp->SetBinContent(i+1,hLep1PtTruth[0]->GetBinContent(i+1)+sqrt(pdfUncLep1PtUp[i]/(npdfLep1PtP-1)));
      hLep1PtTruthPDFDown->SetBinContent(i+1,hLep1PtTruth[0]->GetBinContent(i+1)-sqrt(pdfUncLep1PtDown[i]/(npdfLep1PtM-1)));

      hLep1PtTruthAlphasUp->SetBinContent(i+1,hLep1PtTruth[0]->GetBinContent(i+1)+alphasUncLep1PtUp[i]);
      hLep1PtTruthAlphasDown->SetBinContent(i+1,hLep1PtTruth[0]->GetBinContent(i+1)-alphasUncLep1PtDown[i]);      
    }
  for(int j=0;j!=hLep1PtTruth[0]->GetNbinsX();++j)
    {
      hLep1PtTruthNominal->SetBinContent(j+1,hLep1PtTruthNominal->GetBinContent(j+1)/hLep1PtTruthNominal->GetBinWidth(j+1));
      hLep1PtTruthNominal->SetBinError(j+1,hLep1PtTruthNominal->GetBinError(j+1)/hLep1PtTruthNominal->GetBinWidth(j+1));

      hLep1PtTruthScaleUp->SetBinContent(j+1,hLep1PtTruthScaleUp->GetBinContent(j+1)/hLep1PtTruthScaleUp->GetBinWidth(j+1));
      hLep1PtTruthScaleUp->SetBinError(j+1,hLep1PtTruthScaleUp->GetBinError(j+1)/hLep1PtTruthScaleUp->GetBinWidth(j+1));
      hLep1PtTruthScaleDown->SetBinContent(j+1,hLep1PtTruthScaleDown->GetBinContent(j+1)/hLep1PtTruthScaleDown->GetBinWidth(j+1));
      hLep1PtTruthScaleDown->SetBinError(j+1,hLep1PtTruthScaleDown->GetBinError(j+1)/hLep1PtTruthScaleDown->GetBinWidth(j+1));

      hLep1PtTruthPDFUp->SetBinContent(j+1,hLep1PtTruthPDFUp->GetBinContent(j+1)/hLep1PtTruthPDFUp->GetBinWidth(j+1));
      hLep1PtTruthPDFUp->SetBinError(j+1,hLep1PtTruthPDFUp->GetBinError(j+1)/hLep1PtTruthPDFUp->GetBinWidth(j+1));
      hLep1PtTruthPDFDown->SetBinContent(j+1,hLep1PtTruthPDFDown->GetBinContent(j+1)/hLep1PtTruthPDFDown->GetBinWidth(j+1));
      hLep1PtTruthPDFDown->SetBinError(j+1,hLep1PtTruthPDFDown->GetBinError(j+1)/hLep1PtTruthPDFDown->GetBinWidth(j+1));

      hLep1PtTruthAlphasUp->SetBinContent(j+1,hLep1PtTruthAlphasUp->GetBinContent(j+1)/hLep1PtTruthAlphasUp->GetBinWidth(j+1));
      hLep1PtTruthAlphasUp->SetBinError(j+1,hLep1PtTruthAlphasUp->GetBinError(j+1)/hLep1PtTruthAlphasUp->GetBinWidth(j+1));
      hLep1PtTruthAlphasDown->SetBinContent(j+1,hLep1PtTruthAlphasDown->GetBinContent(j+1)/hLep1PtTruthAlphasDown->GetBinWidth(j+1));
      hLep1PtTruthAlphasDown->SetBinError(j+1,hLep1PtTruthAlphasDown->GetBinError(j+1)/hLep1PtTruthAlphasDown->GetBinWidth(j+1));
    }
  hLep1PtTruthNominal->Scale(1./lumi);
  hLep1PtTruthScaleUp->Scale(1./lumi);
  hLep1PtTruthScaleDown->Scale(1./lumi);
  hLep1PtTruthPDFUp->Scale(1./lumi);
  hLep1PtTruthPDFDown->Scale(1./lumi);
  hLep1PtTruthAlphasUp->Scale(1./lumi);
  hLep1PtTruthAlphasDown->Scale(1./lumi);

  

  //---------------------------------------------------------------------------
  //                             Lep2 Pt
  //---------------------------------------------------------------------------

  double pdfUncLep2PtUp[20];
  double pdfUncLep2PtDown[20];
  double alphasUncLep2PtUp[20];
  double alphasUncLep2PtDown[20];
  double scaleUncLep2PtUp[20];
  double scaleUncLep2PtDown[20];

  int npdfLep2PtP=0;
  int npdfLep2PtM=0;

  for(int i=0;i!=hLep2PtTruth[0]->GetNbinsX();i++)
    {
      scaleUncLep2PtUp[i]=0;
      scaleUncLep2PtDown[i]=0;
      pdfUncLep2PtUp[i]=0;
      pdfUncLep2PtDown[i]=0;
      for(int j=1;j!=109;j++)
	{
	  if(j<=8)
	    {
	      if(hLep2PtTruth[0]->GetBinContent(i+1)<hLep2PtTruth[j]->GetBinContent(i+1))
		{
		  if((hLep2PtTruth[j]->GetBinContent(i+1)-hLep2PtTruth[0]->GetBinContent(i+1))>scaleUncLep2PtUp[i])
		    {
		      scaleUncLep2PtUp[i]=hLep2PtTruth[j]->GetBinContent(i+1)-hLep2PtTruth[0]->GetBinContent(i+1);
		    }
		}
	      else
		{
		  if((hLep2PtTruth[0]->GetBinContent(i+1)-hLep2PtTruth[j]->GetBinContent(i+1))>scaleUncLep2PtDown[i])
		    {
		      scaleUncLep2PtDown[i]=hLep2PtTruth[0]->GetBinContent(i+1)-hLep2PtTruth[j]->GetBinContent(i+1);
		    }
		}
	    }
	  else
	    {
	      if(hLep2PtTruth[0]->GetBinContent(i+1)<hLep2PtTruth[j]->GetBinContent(i+1))
		{
		  pdfUncLep2PtUp[i]+=pow(hLep2PtTruth[0]->GetBinContent(i+1)-hLep2PtTruth[j]->GetBinContent(i+1),2);
		  npdfLep2PtP++;
		}
	      else
		{
		  pdfUncLep2PtDown[i]+=pow(hLep2PtTruth[j]->GetBinContent(i+1)-hLep2PtTruth[0]->GetBinContent(i+1),2);
		  npdfLep2PtM++;
		}
	    }
	}
      alphasUncLep2PtUp[i]=fabs(hLep2PtTruth[110]->GetBinContent(i+1)-hLep2PtTruth[109]->GetBinContent(i+1))/2;
      alphasUncLep2PtDown[i]=fabs(hLep2PtTruth[110]->GetBinContent(i+1)-hLep2PtTruth[109]->GetBinContent(i+1))/2;

      hLep2PtTruthNominal->SetBinContent(i+1,hLep2PtTruth[0]->GetBinContent(i+1));
      hLep2PtTruthNominal->SetBinError(i+1,hLep2PtTruth[0]->GetBinError(i+1));

      hLep2PtTruthScaleUp->SetBinContent(i+1,hLep2PtTruth[0]->GetBinContent(i+1)+scaleUncLep2PtUp[i]);
      hLep2PtTruthScaleDown->SetBinContent(i+1,hLep2PtTruth[0]->GetBinContent(i+1)-scaleUncLep2PtDown[i]);

      hLep2PtTruthPDFUp->SetBinContent(i+1,hLep2PtTruth[0]->GetBinContent(i+1)+sqrt(pdfUncLep2PtUp[i]/(npdfLep2PtP-1)));
      hLep2PtTruthPDFDown->SetBinContent(i+1,hLep2PtTruth[0]->GetBinContent(i+1)-sqrt(pdfUncLep2PtDown[i]/(npdfLep2PtM-1)));

      hLep2PtTruthAlphasUp->SetBinContent(i+1,hLep2PtTruth[0]->GetBinContent(i+1)+alphasUncLep2PtUp[i]);
      hLep2PtTruthAlphasDown->SetBinContent(i+1,hLep2PtTruth[0]->GetBinContent(i+1)-alphasUncLep2PtDown[i]);      
    }
  for(int j=0;j!=hLep2PtTruth[0]->GetNbinsX();++j)
    {
      hLep2PtTruthNominal->SetBinContent(j+1,hLep2PtTruthNominal->GetBinContent(j+1)/hLep2PtTruthNominal->GetBinWidth(j+1));
      hLep2PtTruthNominal->SetBinError(j+1,hLep2PtTruthNominal->GetBinError(j+1)/hLep2PtTruthNominal->GetBinWidth(j+1));

      hLep2PtTruthScaleUp->SetBinContent(j+1,hLep2PtTruthScaleUp->GetBinContent(j+1)/hLep2PtTruthScaleUp->GetBinWidth(j+1));
      hLep2PtTruthScaleUp->SetBinError(j+1,hLep2PtTruthScaleUp->GetBinError(j+1)/hLep2PtTruthScaleUp->GetBinWidth(j+1));
      hLep2PtTruthScaleDown->SetBinContent(j+1,hLep2PtTruthScaleDown->GetBinContent(j+1)/hLep2PtTruthScaleDown->GetBinWidth(j+1));
      hLep2PtTruthScaleDown->SetBinError(j+1,hLep2PtTruthScaleDown->GetBinError(j+1)/hLep2PtTruthScaleDown->GetBinWidth(j+1));

      hLep2PtTruthPDFUp->SetBinContent(j+1,hLep2PtTruthPDFUp->GetBinContent(j+1)/hLep2PtTruthPDFUp->GetBinWidth(j+1));
      hLep2PtTruthPDFUp->SetBinError(j+1,hLep2PtTruthPDFUp->GetBinError(j+1)/hLep2PtTruthPDFUp->GetBinWidth(j+1));
      hLep2PtTruthPDFDown->SetBinContent(j+1,hLep2PtTruthPDFDown->GetBinContent(j+1)/hLep2PtTruthPDFDown->GetBinWidth(j+1));
      hLep2PtTruthPDFDown->SetBinError(j+1,hLep2PtTruthPDFDown->GetBinError(j+1)/hLep2PtTruthPDFDown->GetBinWidth(j+1));

      hLep2PtTruthAlphasUp->SetBinContent(j+1,hLep2PtTruthAlphasUp->GetBinContent(j+1)/hLep2PtTruthAlphasUp->GetBinWidth(j+1));
      hLep2PtTruthAlphasUp->SetBinError(j+1,hLep2PtTruthAlphasUp->GetBinError(j+1)/hLep2PtTruthAlphasUp->GetBinWidth(j+1));
      hLep2PtTruthAlphasDown->SetBinContent(j+1,hLep2PtTruthAlphasDown->GetBinContent(j+1)/hLep2PtTruthAlphasDown->GetBinWidth(j+1));
      hLep2PtTruthAlphasDown->SetBinError(j+1,hLep2PtTruthAlphasDown->GetBinError(j+1)/hLep2PtTruthAlphasDown->GetBinWidth(j+1));
    }
  hLep2PtTruthNominal->Scale(1./lumi);
  hLep2PtTruthScaleUp->Scale(1./lumi);
  hLep2PtTruthScaleDown->Scale(1./lumi);
  hLep2PtTruthPDFUp->Scale(1./lumi);
  hLep2PtTruthPDFDown->Scale(1./lumi);
  hLep2PtTruthAlphasUp->Scale(1./lumi);
  hLep2PtTruthAlphasDown->Scale(1./lumi);

  //---------------------------------------------------------------------------
  //                             LepNeg Pt
  //---------------------------------------------------------------------------

  double pdfUncLepNegPtUp[25];
  double pdfUncLepNegPtDown[25];
  double alphasUncLepNegPtUp[25];
  double alphasUncLepNegPtDown[25];
  double scaleUncLepNegPtUp[25];
  double scaleUncLepNegPtDown[25];

  int npdfLepNegPtP=0;
  int npdfLepNegPtM=0;

  for(int i=0;i!=hLepNegPtTruth[0]->GetNbinsX();i++)
    {
      scaleUncLepNegPtUp[i]=0;
      scaleUncLepNegPtDown[i]=0;
      pdfUncLepNegPtUp[i]=0;
      pdfUncLepNegPtDown[i]=0;
      for(int j=1;j!=109;j++)
	{
	  if(j<=8)
	    {
	      if(hLepNegPtTruth[0]->GetBinContent(i+1)<hLepNegPtTruth[j]->GetBinContent(i+1))
		{
		  if((hLepNegPtTruth[j]->GetBinContent(i+1)-hLepNegPtTruth[0]->GetBinContent(i+1))>scaleUncLepNegPtUp[i])
		    {
		      scaleUncLepNegPtUp[i]=hLepNegPtTruth[j]->GetBinContent(i+1)-hLepNegPtTruth[0]->GetBinContent(i+1);
		    }
		}
	      else
		{
		  if((hLepNegPtTruth[0]->GetBinContent(i+1)-hLepNegPtTruth[j]->GetBinContent(i+1))>scaleUncLepNegPtDown[i])
		    {
		      scaleUncLepNegPtDown[i]=hLepNegPtTruth[0]->GetBinContent(i+1)-hLepNegPtTruth[j]->GetBinContent(i+1);
		    }
		}
	    }
	  else
	    {
	      if(hLepNegPtTruth[0]->GetBinContent(i+1)<hLepNegPtTruth[j]->GetBinContent(i+1))
		{
		  pdfUncLepNegPtUp[i]+=pow(hLepNegPtTruth[0]->GetBinContent(i+1)-hLepNegPtTruth[j]->GetBinContent(i+1),2);
		  npdfLepNegPtP++;
		}
	      else
		{
		  pdfUncLepNegPtDown[i]+=pow(hLepNegPtTruth[j]->GetBinContent(i+1)-hLepNegPtTruth[0]->GetBinContent(i+1),2);
		  npdfLepNegPtM++;
		}
	    }
	}
      alphasUncLepNegPtUp[i]=fabs(hLepNegPtTruth[110]->GetBinContent(i+1)-hLepNegPtTruth[109]->GetBinContent(i+1))/2;
      alphasUncLepNegPtDown[i]=fabs(hLepNegPtTruth[110]->GetBinContent(i+1)-hLepNegPtTruth[109]->GetBinContent(i+1))/2;

      hLepNegPtTruthNominal->SetBinContent(i+1,hLepNegPtTruth[0]->GetBinContent(i+1));
      hLepNegPtTruthNominal->SetBinError(i+1,hLepNegPtTruth[0]->GetBinError(i+1));

      hLepNegPtTruthScaleUp->SetBinContent(i+1,hLepNegPtTruth[0]->GetBinContent(i+1)+scaleUncLepNegPtUp[i]);
      hLepNegPtTruthScaleDown->SetBinContent(i+1,hLepNegPtTruth[0]->GetBinContent(i+1)-scaleUncLepNegPtDown[i]);

      hLepNegPtTruthPDFUp->SetBinContent(i+1,hLepNegPtTruth[0]->GetBinContent(i+1)+sqrt(pdfUncLepNegPtUp[i]/(npdfLepNegPtP-1)));
      hLepNegPtTruthPDFDown->SetBinContent(i+1,hLepNegPtTruth[0]->GetBinContent(i+1)-sqrt(pdfUncLepNegPtDown[i]/(npdfLepNegPtM-1)));

      hLepNegPtTruthAlphasUp->SetBinContent(i+1,hLepNegPtTruth[0]->GetBinContent(i+1)+alphasUncLepNegPtUp[i]);
      hLepNegPtTruthAlphasDown->SetBinContent(i+1,hLepNegPtTruth[0]->GetBinContent(i+1)-alphasUncLepNegPtDown[i]);      
    }
  for(int j=0;j!=hLepNegPtTruth[0]->GetNbinsX();++j)
    {
      hLepNegPtTruthNominal->SetBinContent(j+1,hLepNegPtTruthNominal->GetBinContent(j+1)/hLepNegPtTruthNominal->GetBinWidth(j+1));
      hLepNegPtTruthNominal->SetBinError(j+1,hLepNegPtTruthNominal->GetBinError(j+1)/hLepNegPtTruthNominal->GetBinWidth(j+1));

      hLepNegPtTruthScaleUp->SetBinContent(j+1,hLepNegPtTruthScaleUp->GetBinContent(j+1)/hLepNegPtTruthScaleUp->GetBinWidth(j+1));
      hLepNegPtTruthScaleUp->SetBinError(j+1,hLepNegPtTruthScaleUp->GetBinError(j+1)/hLepNegPtTruthScaleUp->GetBinWidth(j+1));
      hLepNegPtTruthScaleDown->SetBinContent(j+1,hLepNegPtTruthScaleDown->GetBinContent(j+1)/hLepNegPtTruthScaleDown->GetBinWidth(j+1));
      hLepNegPtTruthScaleDown->SetBinError(j+1,hLepNegPtTruthScaleDown->GetBinError(j+1)/hLepNegPtTruthScaleDown->GetBinWidth(j+1));

      hLepNegPtTruthPDFUp->SetBinContent(j+1,hLepNegPtTruthPDFUp->GetBinContent(j+1)/hLepNegPtTruthPDFUp->GetBinWidth(j+1));
      hLepNegPtTruthPDFUp->SetBinError(j+1,hLepNegPtTruthPDFUp->GetBinError(j+1)/hLepNegPtTruthPDFUp->GetBinWidth(j+1));
      hLepNegPtTruthPDFDown->SetBinContent(j+1,hLepNegPtTruthPDFDown->GetBinContent(j+1)/hLepNegPtTruthPDFDown->GetBinWidth(j+1));
      hLepNegPtTruthPDFDown->SetBinError(j+1,hLepNegPtTruthPDFDown->GetBinError(j+1)/hLepNegPtTruthPDFDown->GetBinWidth(j+1));

      hLepNegPtTruthAlphasUp->SetBinContent(j+1,hLepNegPtTruthAlphasUp->GetBinContent(j+1)/hLepNegPtTruthAlphasUp->GetBinWidth(j+1));
      hLepNegPtTruthAlphasUp->SetBinError(j+1,hLepNegPtTruthAlphasUp->GetBinError(j+1)/hLepNegPtTruthAlphasUp->GetBinWidth(j+1));
      hLepNegPtTruthAlphasDown->SetBinContent(j+1,hLepNegPtTruthAlphasDown->GetBinContent(j+1)/hLepNegPtTruthAlphasDown->GetBinWidth(j+1));
      hLepNegPtTruthAlphasDown->SetBinError(j+1,hLepNegPtTruthAlphasDown->GetBinError(j+1)/hLepNegPtTruthAlphasDown->GetBinWidth(j+1));
    }
  hLepNegPtTruthNominal->Scale(1./lumi);
  hLepNegPtTruthScaleUp->Scale(1./lumi);
  hLepNegPtTruthScaleDown->Scale(1./lumi);
  hLepNegPtTruthPDFUp->Scale(1./lumi);
  hLepNegPtTruthPDFDown->Scale(1./lumi);
  hLepNegPtTruthAlphasUp->Scale(1./lumi);
  hLepNegPtTruthAlphasDown->Scale(1./lumi);

  //---------------------------------------------------------------------------
  //                             LepPos Pt
  //---------------------------------------------------------------------------

  double pdfUncLepPosPtUp[25];
  double pdfUncLepPosPtDown[25];
  double alphasUncLepPosPtUp[25];
  double alphasUncLepPosPtDown[25];
  double scaleUncLepPosPtUp[25];
  double scaleUncLepPosPtDown[25];

  int npdfLepPosPtP=0;
  int npdfLepPosPtM=0;

  for(int i=0;i!=hLepPosPtTruth[0]->GetNbinsX();i++)
    {
      scaleUncLepPosPtUp[i]=0;
      scaleUncLepPosPtDown[i]=0;
      pdfUncLepPosPtUp[i]=0;
      pdfUncLepPosPtDown[i]=0;
      for(int j=1;j!=109;j++)
	{
	  if(j<=8)
	    {
	      if(hLepPosPtTruth[0]->GetBinContent(i+1)<hLepPosPtTruth[j]->GetBinContent(i+1))
		{
		  if((hLepPosPtTruth[j]->GetBinContent(i+1)-hLepPosPtTruth[0]->GetBinContent(i+1))>scaleUncLepPosPtUp[i])
		    {
		      scaleUncLepPosPtUp[i]=hLepPosPtTruth[j]->GetBinContent(i+1)-hLepPosPtTruth[0]->GetBinContent(i+1);
		    }
		}
	      else
		{
		  if((hLepPosPtTruth[0]->GetBinContent(i+1)-hLepPosPtTruth[j]->GetBinContent(i+1))>scaleUncLepPosPtDown[i])
		    {
		      scaleUncLepPosPtDown[i]=hLepPosPtTruth[0]->GetBinContent(i+1)-hLepPosPtTruth[j]->GetBinContent(i+1);
		    }
		}
	    }
	  else
	    {
	      if(hLepPosPtTruth[0]->GetBinContent(i+1)<hLepPosPtTruth[j]->GetBinContent(i+1))
		{
		  pdfUncLepPosPtUp[i]+=pow(hLepPosPtTruth[0]->GetBinContent(i+1)-hLepPosPtTruth[j]->GetBinContent(i+1),2);
		  npdfLepPosPtP++;
		}
	      else
		{
		  pdfUncLepPosPtDown[i]+=pow(hLepPosPtTruth[j]->GetBinContent(i+1)-hLepPosPtTruth[0]->GetBinContent(i+1),2);
		  npdfLepPosPtM++;
		}
	    }
	}
      alphasUncLepPosPtUp[i]=fabs(hLepPosPtTruth[110]->GetBinContent(i+1)-hLepPosPtTruth[109]->GetBinContent(i+1))/2;
      alphasUncLepPosPtDown[i]=fabs(hLepPosPtTruth[110]->GetBinContent(i+1)-hLepPosPtTruth[109]->GetBinContent(i+1))/2;

      hLepPosPtTruthNominal->SetBinContent(i+1,hLepPosPtTruth[0]->GetBinContent(i+1));
      hLepPosPtTruthNominal->SetBinError(i+1,hLepPosPtTruth[0]->GetBinError(i+1));

      hLepPosPtTruthScaleUp->SetBinContent(i+1,hLepPosPtTruth[0]->GetBinContent(i+1)+scaleUncLepPosPtUp[i]);
      hLepPosPtTruthScaleDown->SetBinContent(i+1,hLepPosPtTruth[0]->GetBinContent(i+1)-scaleUncLepPosPtDown[i]);

      hLepPosPtTruthPDFUp->SetBinContent(i+1,hLepPosPtTruth[0]->GetBinContent(i+1)+sqrt(pdfUncLepPosPtUp[i]/(npdfLepPosPtP-1)));
      hLepPosPtTruthPDFDown->SetBinContent(i+1,hLepPosPtTruth[0]->GetBinContent(i+1)-sqrt(pdfUncLepPosPtDown[i]/(npdfLepPosPtM-1)));

      hLepPosPtTruthAlphasUp->SetBinContent(i+1,hLepPosPtTruth[0]->GetBinContent(i+1)+alphasUncLepPosPtUp[i]);
      hLepPosPtTruthAlphasDown->SetBinContent(i+1,hLepPosPtTruth[0]->GetBinContent(i+1)-alphasUncLepPosPtDown[i]);      
    }
  for(int j=0;j!=hLepPosPtTruth[0]->GetNbinsX();++j)
    {
      hLepPosPtTruthNominal->SetBinContent(j+1,hLepPosPtTruthNominal->GetBinContent(j+1)/hLepPosPtTruthNominal->GetBinWidth(j+1));
      hLepPosPtTruthNominal->SetBinError(j+1,hLepPosPtTruthNominal->GetBinError(j+1)/hLepPosPtTruthNominal->GetBinWidth(j+1));

      hLepPosPtTruthScaleUp->SetBinContent(j+1,hLepPosPtTruthScaleUp->GetBinContent(j+1)/hLepPosPtTruthScaleUp->GetBinWidth(j+1));
      hLepPosPtTruthScaleUp->SetBinError(j+1,hLepPosPtTruthScaleUp->GetBinError(j+1)/hLepPosPtTruthScaleUp->GetBinWidth(j+1));
      hLepPosPtTruthScaleDown->SetBinContent(j+1,hLepPosPtTruthScaleDown->GetBinContent(j+1)/hLepPosPtTruthScaleDown->GetBinWidth(j+1));
      hLepPosPtTruthScaleDown->SetBinError(j+1,hLepPosPtTruthScaleDown->GetBinError(j+1)/hLepPosPtTruthScaleDown->GetBinWidth(j+1));

      hLepPosPtTruthPDFUp->SetBinContent(j+1,hLepPosPtTruthPDFUp->GetBinContent(j+1)/hLepPosPtTruthPDFUp->GetBinWidth(j+1));
      hLepPosPtTruthPDFUp->SetBinError(j+1,hLepPosPtTruthPDFUp->GetBinError(j+1)/hLepPosPtTruthPDFUp->GetBinWidth(j+1));
      hLepPosPtTruthPDFDown->SetBinContent(j+1,hLepPosPtTruthPDFDown->GetBinContent(j+1)/hLepPosPtTruthPDFDown->GetBinWidth(j+1));
      hLepPosPtTruthPDFDown->SetBinError(j+1,hLepPosPtTruthPDFDown->GetBinError(j+1)/hLepPosPtTruthPDFDown->GetBinWidth(j+1));

      hLepPosPtTruthAlphasUp->SetBinContent(j+1,hLepPosPtTruthAlphasUp->GetBinContent(j+1)/hLepPosPtTruthAlphasUp->GetBinWidth(j+1));
      hLepPosPtTruthAlphasUp->SetBinError(j+1,hLepPosPtTruthAlphasUp->GetBinError(j+1)/hLepPosPtTruthAlphasUp->GetBinWidth(j+1));
      hLepPosPtTruthAlphasDown->SetBinContent(j+1,hLepPosPtTruthAlphasDown->GetBinContent(j+1)/hLepPosPtTruthAlphasDown->GetBinWidth(j+1));
      hLepPosPtTruthAlphasDown->SetBinError(j+1,hLepPosPtTruthAlphasDown->GetBinError(j+1)/hLepPosPtTruthAlphasDown->GetBinWidth(j+1));
    }
  hLepPosPtTruthNominal->Scale(1./lumi);
  hLepPosPtTruthScaleUp->Scale(1./lumi);
  hLepPosPtTruthScaleDown->Scale(1./lumi);
  hLepPosPtTruthPDFUp->Scale(1./lumi);
  hLepPosPtTruthPDFDown->Scale(1./lumi);
  hLepPosPtTruthAlphasUp->Scale(1./lumi);
  hLepPosPtTruthAlphasDown->Scale(1./lumi);

  //---------------------------------------------------------------------------
  //                             Lep1 Eta
  //---------------------------------------------------------------------------

  double pdfUncLep1EtaUp[24];
  double pdfUncLep1EtaDown[24];
  double alphasUncLep1EtaUp[24];
  double alphasUncLep1EtaDown[24];
  double scaleUncLep1EtaUp[24];
  double scaleUncLep1EtaDown[24];

  int npdfLep1EtaP=0;
  int npdfLep1EtaM=0;

  for(int i=0;i!=hLep1EtaTruth[0]->GetNbinsX();i++)
    {
      scaleUncLep1EtaUp[i]=0;
      scaleUncLep1EtaDown[i]=0;
      pdfUncLep1EtaUp[i]=0;
      pdfUncLep1EtaDown[i]=0;
      for(int j=1;j!=109;j++)
	{
	  if(j<=8)
	    {
	      if(hLep1EtaTruth[0]->GetBinContent(i+1)<hLep1EtaTruth[j]->GetBinContent(i+1))
		{
		  if((hLep1EtaTruth[j]->GetBinContent(i+1)-hLep1EtaTruth[0]->GetBinContent(i+1))>scaleUncLep1EtaUp[i])
		    {
		      scaleUncLep1EtaUp[i]=hLep1EtaTruth[j]->GetBinContent(i+1)-hLep1EtaTruth[0]->GetBinContent(i+1);
		    }
		}
	      else
		{
		  if((hLep1EtaTruth[0]->GetBinContent(i+1)-hLep1EtaTruth[j]->GetBinContent(i+1))>scaleUncLep1EtaDown[i])
		    {
		      scaleUncLep1EtaDown[i]=hLep1EtaTruth[0]->GetBinContent(i+1)-hLep1EtaTruth[j]->GetBinContent(i+1);
		    }
		}
	    }
	  else
	    {
	      if(hLep1EtaTruth[0]->GetBinContent(i+1)<hLep1EtaTruth[j]->GetBinContent(i+1))
		{
		  pdfUncLep1EtaUp[i]+=pow(hLep1EtaTruth[0]->GetBinContent(i+1)-hLep1EtaTruth[j]->GetBinContent(i+1),2);
		  npdfLep1EtaP++;
		}
	      else
		{
		  pdfUncLep1EtaDown[i]+=pow(hLep1EtaTruth[j]->GetBinContent(i+1)-hLep1EtaTruth[0]->GetBinContent(i+1),2);
		  npdfLep1EtaM++;
		}
	    }
	}
      alphasUncLep1EtaUp[i]=fabs(hLep1EtaTruth[110]->GetBinContent(i+1)-hLep1EtaTruth[109]->GetBinContent(i+1))/2;
      alphasUncLep1EtaDown[i]=fabs(hLep1EtaTruth[110]->GetBinContent(i+1)-hLep1EtaTruth[109]->GetBinContent(i+1))/2;

      hLep1EtaTruthNominal->SetBinContent(i+1,hLep1EtaTruth[0]->GetBinContent(i+1));
      hLep1EtaTruthNominal->SetBinError(i+1,hLep1EtaTruth[0]->GetBinError(i+1));

      hLep1EtaTruthScaleUp->SetBinContent(i+1,hLep1EtaTruth[0]->GetBinContent(i+1)+scaleUncLep1EtaUp[i]);
      hLep1EtaTruthScaleDown->SetBinContent(i+1,hLep1EtaTruth[0]->GetBinContent(i+1)-scaleUncLep1EtaDown[i]);

      hLep1EtaTruthPDFUp->SetBinContent(i+1,hLep1EtaTruth[0]->GetBinContent(i+1)+sqrt(pdfUncLep1EtaUp[i]/(npdfLep1EtaP-1)));
      hLep1EtaTruthPDFDown->SetBinContent(i+1,hLep1EtaTruth[0]->GetBinContent(i+1)-sqrt(pdfUncLep1EtaDown[i]/(npdfLep1EtaM-1)));

      hLep1EtaTruthAlphasUp->SetBinContent(i+1,hLep1EtaTruth[0]->GetBinContent(i+1)+alphasUncLep1EtaUp[i]);
      hLep1EtaTruthAlphasDown->SetBinContent(i+1,hLep1EtaTruth[0]->GetBinContent(i+1)-alphasUncLep1EtaDown[i]);      
    }
  for(int j=0;j!=hLep1EtaTruth[0]->GetNbinsX();++j)
    {
      hLep1EtaTruthNominal->SetBinContent(j+1,hLep1EtaTruthNominal->GetBinContent(j+1)/hLep1EtaTruthNominal->GetBinWidth(j+1));
      hLep1EtaTruthNominal->SetBinError(j+1,hLep1EtaTruthNominal->GetBinError(j+1)/hLep1EtaTruthNominal->GetBinWidth(j+1));

      hLep1EtaTruthScaleUp->SetBinContent(j+1,hLep1EtaTruthScaleUp->GetBinContent(j+1)/hLep1EtaTruthScaleUp->GetBinWidth(j+1));
      hLep1EtaTruthScaleUp->SetBinError(j+1,hLep1EtaTruthScaleUp->GetBinError(j+1)/hLep1EtaTruthScaleUp->GetBinWidth(j+1));
      hLep1EtaTruthScaleDown->SetBinContent(j+1,hLep1EtaTruthScaleDown->GetBinContent(j+1)/hLep1EtaTruthScaleDown->GetBinWidth(j+1));
      hLep1EtaTruthScaleDown->SetBinError(j+1,hLep1EtaTruthScaleDown->GetBinError(j+1)/hLep1EtaTruthScaleDown->GetBinWidth(j+1));

      hLep1EtaTruthPDFUp->SetBinContent(j+1,hLep1EtaTruthPDFUp->GetBinContent(j+1)/hLep1EtaTruthPDFUp->GetBinWidth(j+1));
      hLep1EtaTruthPDFUp->SetBinError(j+1,hLep1EtaTruthPDFUp->GetBinError(j+1)/hLep1EtaTruthPDFUp->GetBinWidth(j+1));
      hLep1EtaTruthPDFDown->SetBinContent(j+1,hLep1EtaTruthPDFDown->GetBinContent(j+1)/hLep1EtaTruthPDFDown->GetBinWidth(j+1));
      hLep1EtaTruthPDFDown->SetBinError(j+1,hLep1EtaTruthPDFDown->GetBinError(j+1)/hLep1EtaTruthPDFDown->GetBinWidth(j+1));

      hLep1EtaTruthAlphasUp->SetBinContent(j+1,hLep1EtaTruthAlphasUp->GetBinContent(j+1)/hLep1EtaTruthAlphasUp->GetBinWidth(j+1));
      hLep1EtaTruthAlphasUp->SetBinError(j+1,hLep1EtaTruthAlphasUp->GetBinError(j+1)/hLep1EtaTruthAlphasUp->GetBinWidth(j+1));
      hLep1EtaTruthAlphasDown->SetBinContent(j+1,hLep1EtaTruthAlphasDown->GetBinContent(j+1)/hLep1EtaTruthAlphasDown->GetBinWidth(j+1));
      hLep1EtaTruthAlphasDown->SetBinError(j+1,hLep1EtaTruthAlphasDown->GetBinError(j+1)/hLep1EtaTruthAlphasDown->GetBinWidth(j+1));
    }
  hLep1EtaTruthNominal->Scale(1./lumi);
  hLep1EtaTruthScaleUp->Scale(1./lumi);
  hLep1EtaTruthScaleDown->Scale(1./lumi);
  hLep1EtaTruthPDFUp->Scale(1./lumi);
  hLep1EtaTruthPDFDown->Scale(1./lumi);
  hLep1EtaTruthAlphasUp->Scale(1./lumi);
  hLep1EtaTruthAlphasDown->Scale(1./lumi);

  

  //---------------------------------------------------------------------------
  //                             Lep2 Eta
  //---------------------------------------------------------------------------

  double pdfUncLep2EtaUp[24];
  double pdfUncLep2EtaDown[24];
  double alphasUncLep2EtaUp[24];
  double alphasUncLep2EtaDown[24];
  double scaleUncLep2EtaUp[24];
  double scaleUncLep2EtaDown[24];

  int npdfLep2EtaP=0;
  int npdfLep2EtaM=0;

  for(int i=0;i!=hLep2EtaTruth[0]->GetNbinsX();i++)
    {
      scaleUncLep2EtaUp[i]=0;
      scaleUncLep2EtaDown[i]=0;
      pdfUncLep2EtaUp[i]=0;
      pdfUncLep2EtaDown[i]=0;
      for(int j=1;j!=109;j++)
	{
	  if(j<=8)
	    {
	      if(hLep2EtaTruth[0]->GetBinContent(i+1)<hLep2EtaTruth[j]->GetBinContent(i+1))
		{
		  if((hLep2EtaTruth[j]->GetBinContent(i+1)-hLep2EtaTruth[0]->GetBinContent(i+1))>scaleUncLep2EtaUp[i])
		    {
		      scaleUncLep2EtaUp[i]=hLep2EtaTruth[j]->GetBinContent(i+1)-hLep2EtaTruth[0]->GetBinContent(i+1);
		    }
		}
	      else
		{
		  if((hLep2EtaTruth[0]->GetBinContent(i+1)-hLep2EtaTruth[j]->GetBinContent(i+1))>scaleUncLep2EtaDown[i])
		    {
		      scaleUncLep2EtaDown[i]=hLep2EtaTruth[0]->GetBinContent(i+1)-hLep2EtaTruth[j]->GetBinContent(i+1);
		    }
		}
	    }
	  else
	    {
	      if(hLep2EtaTruth[0]->GetBinContent(i+1)<hLep2EtaTruth[j]->GetBinContent(i+1))
		{
		  pdfUncLep2EtaUp[i]+=pow(hLep2EtaTruth[0]->GetBinContent(i+1)-hLep2EtaTruth[j]->GetBinContent(i+1),2);
		  npdfLep2EtaP++;
		}
	      else
		{
		  pdfUncLep2EtaDown[i]+=pow(hLep2EtaTruth[j]->GetBinContent(i+1)-hLep2EtaTruth[0]->GetBinContent(i+1),2);
		  npdfLep2EtaM++;
		}
	    }
	}
      alphasUncLep2EtaUp[i]=fabs(hLep2EtaTruth[110]->GetBinContent(i+1)-hLep2EtaTruth[109]->GetBinContent(i+1))/2;
      alphasUncLep2EtaDown[i]=fabs(hLep2EtaTruth[110]->GetBinContent(i+1)-hLep2EtaTruth[109]->GetBinContent(i+1))/2;

      hLep2EtaTruthNominal->SetBinContent(i+1,hLep2EtaTruth[0]->GetBinContent(i+1));
      hLep2EtaTruthNominal->SetBinError(i+1,hLep2EtaTruth[0]->GetBinError(i+1));

      hLep2EtaTruthScaleUp->SetBinContent(i+1,hLep2EtaTruth[0]->GetBinContent(i+1)+scaleUncLep2EtaUp[i]);
      hLep2EtaTruthScaleDown->SetBinContent(i+1,hLep2EtaTruth[0]->GetBinContent(i+1)-scaleUncLep2EtaDown[i]);

      hLep2EtaTruthPDFUp->SetBinContent(i+1,hLep2EtaTruth[0]->GetBinContent(i+1)+sqrt(pdfUncLep2EtaUp[i]/(npdfLep2EtaP-1)));
      hLep2EtaTruthPDFDown->SetBinContent(i+1,hLep2EtaTruth[0]->GetBinContent(i+1)-sqrt(pdfUncLep2EtaDown[i]/(npdfLep2EtaM-1)));

      hLep2EtaTruthAlphasUp->SetBinContent(i+1,hLep2EtaTruth[0]->GetBinContent(i+1)+alphasUncLep2EtaUp[i]);
      hLep2EtaTruthAlphasDown->SetBinContent(i+1,hLep2EtaTruth[0]->GetBinContent(i+1)-alphasUncLep2EtaDown[i]);      
    }
  for(int j=0;j!=hLep2EtaTruth[0]->GetNbinsX();++j)
    {
      hLep2EtaTruthNominal->SetBinContent(j+1,hLep2EtaTruthNominal->GetBinContent(j+1)/hLep2EtaTruthNominal->GetBinWidth(j+1));
      hLep2EtaTruthNominal->SetBinError(j+1,hLep2EtaTruthNominal->GetBinError(j+1)/hLep2EtaTruthNominal->GetBinWidth(j+1));

      hLep2EtaTruthScaleUp->SetBinContent(j+1,hLep2EtaTruthScaleUp->GetBinContent(j+1)/hLep2EtaTruthScaleUp->GetBinWidth(j+1));
      hLep2EtaTruthScaleUp->SetBinError(j+1,hLep2EtaTruthScaleUp->GetBinError(j+1)/hLep2EtaTruthScaleUp->GetBinWidth(j+1));
      hLep2EtaTruthScaleDown->SetBinContent(j+1,hLep2EtaTruthScaleDown->GetBinContent(j+1)/hLep2EtaTruthScaleDown->GetBinWidth(j+1));
      hLep2EtaTruthScaleDown->SetBinError(j+1,hLep2EtaTruthScaleDown->GetBinError(j+1)/hLep2EtaTruthScaleDown->GetBinWidth(j+1));

      hLep2EtaTruthPDFUp->SetBinContent(j+1,hLep2EtaTruthPDFUp->GetBinContent(j+1)/hLep2EtaTruthPDFUp->GetBinWidth(j+1));
      hLep2EtaTruthPDFUp->SetBinError(j+1,hLep2EtaTruthPDFUp->GetBinError(j+1)/hLep2EtaTruthPDFUp->GetBinWidth(j+1));
      hLep2EtaTruthPDFDown->SetBinContent(j+1,hLep2EtaTruthPDFDown->GetBinContent(j+1)/hLep2EtaTruthPDFDown->GetBinWidth(j+1));
      hLep2EtaTruthPDFDown->SetBinError(j+1,hLep2EtaTruthPDFDown->GetBinError(j+1)/hLep2EtaTruthPDFDown->GetBinWidth(j+1));

      hLep2EtaTruthAlphasUp->SetBinContent(j+1,hLep2EtaTruthAlphasUp->GetBinContent(j+1)/hLep2EtaTruthAlphasUp->GetBinWidth(j+1));
      hLep2EtaTruthAlphasUp->SetBinError(j+1,hLep2EtaTruthAlphasUp->GetBinError(j+1)/hLep2EtaTruthAlphasUp->GetBinWidth(j+1));
      hLep2EtaTruthAlphasDown->SetBinContent(j+1,hLep2EtaTruthAlphasDown->GetBinContent(j+1)/hLep2EtaTruthAlphasDown->GetBinWidth(j+1));
      hLep2EtaTruthAlphasDown->SetBinError(j+1,hLep2EtaTruthAlphasDown->GetBinError(j+1)/hLep2EtaTruthAlphasDown->GetBinWidth(j+1));
    }
  hLep2EtaTruthNominal->Scale(1./lumi);
  hLep2EtaTruthScaleUp->Scale(1./lumi);
  hLep2EtaTruthScaleDown->Scale(1./lumi);
  hLep2EtaTruthPDFUp->Scale(1./lumi);
  hLep2EtaTruthPDFDown->Scale(1./lumi);
  hLep2EtaTruthAlphasUp->Scale(1./lumi);
  hLep2EtaTruthAlphasDown->Scale(1./lumi);
 
  outFile->cd();
  outFile->Write();
  outFile->Close();
  }
  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl; 

  gBenchmark->Show("plotZeeTheoryUnc"); 
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------

TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = (TH1D*)hData->Clone("hDiff");
  hDiff->SetName(name);
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff=0;
    Double_t err=0;
    if(hData->GetBinContent(ibin)!=0)
      {
	diff = hFit->GetBinContent(ibin)/hData->GetBinContent(ibin);
	err = hFit->GetBinError(ibin)/hData->GetBinContent(ibin);
      }
    hDiff->SetBinContent(ibin,diff);
    hDiff->SetBinError(ibin,err);   
  }
  
  hDiff->GetYaxis()->SetTitleOffset(0.55);
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

TGraphAsymmErrors* TH1TOTGraphAsymmErrors(TH1 *h1){


 if (!h1) cout << "TH1TOTGraph: histogram not found !" << endl;

 TGraphAsymmErrors* g1= new TGraphAsymmErrors();

 Double_t x, y, exh,exl, eyh,eyl;
 for (Int_t i=0; i<h1->GetNbinsX(); i++) {
   y=h1->GetBinContent(i+1);
   eyl=h1->GetBinError(i+1);
   eyh=h1->GetBinError(i+1);
   x=h1->GetBinCenter(i+1);
   exl=h1->GetBinWidth(i+1)/2;
   exh=h1->GetBinWidth(i+1)/2;

   g1->SetPoint(i,x,y);
   g1->SetPointError(i,exl,exh,eyl,eyh);

 }
 return g1;
}

TGraphAsymmErrors* myMakeBandSymmetric(TGraphAsymmErrors* g0, TGraphAsymmErrors* g1,TGraphAsymmErrors* g2) {
   TGraphAsymmErrors* g3= new TGraphAsymmErrors();

  Double_t  x1=0., y1=0., x2=0., y2=0., y0=0, x3=0.;

  Double_t dum;
  for (Int_t i=0; i<g1->GetN(); i++) {
    g0->GetPoint(i, x1,y0);
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    if (i==g1->GetN()-1) x2=x1;
    else                 g2->GetPoint(i+1,x2,dum);

    if (i==0)            x3=x1;
    else                 g2->GetPoint(i-1,x3,dum);

    Double_t tmp=y2;
    if (y1<y2) {y2=y1; y1=tmp;}
    Double_t y3=y0;
    g3->SetPoint(i,x1,y3);

    Double_t binwl=(x1-x3)/2.;
    Double_t binwh=(x2-x1)/2.;
    if (binwl==0.)  binwl= binwh;
    if (binwh==0.)  binwh= binwl;
    g3->SetPointError(i,g0->GetErrorXlow(i),g0->GetErrorXhigh(i),(fabs(y3-y2)+fabs(y1-y3))/2,(fabs(y3-y2)+fabs(y1-y3))/2);

  }
  return g3;
}

void myAddtoBand(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2) {

  Double_t  x1=0., y1=0.,  y2=0., y0=0;
  
  if (g1->GetN()!=g2->GetN())
    cout << " graphs have not the same # of elements " <<g1->GetN()<<" "<<g2->GetN()<< endl;
  Double_t* EYhigh1 = g1-> GetEYhigh();
  Double_t* EYlow1  = g1-> GetEYlow();
  Double_t* EYhigh2 = g2-> GetEYhigh();
  Double_t* EYlow2  = g2-> GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    if (y1==0) y1=1;
    if (y2==0) y2=1;

    Double_t eyh1=0., eyl1=0.,eyh2=0., eyl2=0.;
  
    eyh1=EYhigh1[i];
    eyh2=EYhigh2[i];
    eyh2=sqrt(eyh1*eyh1+eyh2*eyh2);
    g2->SetPointEYhigh(i,eyh2);
    eyl1=EYlow1[i];
    eyl2=EYlow2[i];
    eyl2=sqrt(eyl1*eyl1+eyl2*eyl2);
    g2->SetPointEYlow (i,eyl2);
  }
  return;

}

TGraphAsymmErrors* myTGraphErrorsDivide_noErrGraph2(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2) {

  const Int_t debug=0; 

  TGraphAsymmErrors* g3= new TGraphAsymmErrors();
  Int_t n1=g1->GetN();
  Int_t n2=g2->GetN();

  if (n1!=n2) {
    printf(" vectors do not have same number of entries !  \n");
   return g3;
  }

  Double_t   x1=0.,   y1=0., x2=0., y2=0.;
  Double_t dx1h=0., dx1l=0.;
  Double_t dy1h=0., dy1l=0.;
  Double_t dy2h=0., dy2l=0.;

  Double_t* X1 = g1->GetX();
  Double_t* Y1 = g1->GetY();
  Double_t* EXhigh1 = g1->GetEXhigh();
  Double_t* EXlow1 =  g1->GetEXlow();
  Double_t* EYhigh1 = g1->GetEYhigh();
  Double_t* EYlow1 =  g1->GetEYlow();

  Double_t* X2 = g2->GetX();
  Double_t* Y2 = g2->GetY();
  Double_t* EXhigh2 = g2->GetEXhigh();
  Double_t* EXlow2 =  g2->GetEXlow();
  Double_t* EYhigh2 = g2->GetEYhigh();
  Double_t* EYlow2 =  g2->GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i,x1,y1);
    g2->GetPoint(i,x2,y2);
    dx1h  = EXhigh1[i];
    dx1l  = EXlow1[i];
    if (y1!=0.) dy1h  = EYhigh1[i]/y1;
    else        dy1h  = 0.;
    if (y2!=0.) dy2h  = EYhigh2[i]/y2;
    else        dy2h  = 0.;
    if (y1!=0.) dy1l  = EYlow1 [i]/y1;
    else        dy1l  = 0.;
    if (y2!=0.) dy2l  = EYlow2 [i]/y2;
    else        dy2l  = 0.;
    
    if (y2!=0.) g3->SetPoint(i, x1,y1/y2);
    else       g3->SetPoint(i, x1,y2);
    Double_t el=0.; Double_t eh=0.;

    if (y2!=0.) el=sqrt(dy1l*dy1l)*(y1/y2);
    if (y2!=0.) eh=sqrt(dy1h*dy1h)*(y1/y2);

    g3->SetPointError(i,dx1h,dx1l,fabs(el),fabs(eh));

  }  
  return g3;

}
