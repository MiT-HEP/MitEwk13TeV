//================================================================================================
//
// Perform fit to extract W->munu signal
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

#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
// RooFit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "RooMinuit.h"

// Bacon
#include "BaconAna/DataFormats/interface/TGenParticle.hh"  

// Custom tools from this package
#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"              // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../Utils/AppEffSF.cc"

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

void fillMETs( TH1D** h, vector<double> met, int nMET, double wgtLum );
void fillWeights( TH1D** h, double met, int nWeight, vector<double> wgtLum );

void makeUncMT(vector<Double_t> &metVars, vector<Double_t> &metVarsPhi, TLorentzVector* lep);

// check a flat Iso cut
bool pass2015Iso(double relIso, double eta){
  if(abs(eta)<1.4442){
    if(relIso >= 0.0354) return false;
    // if(relIso >= 0.02) return false;
  } else {
    if(relIso >= 0.0646) return false;
    // if(relIso >= 0.02) return false;
  }
  return true;
}


// global variables yolo
const  Int_t linecolorW   = kOrange-3;
const  Int_t fillcolorW   = kOrange-2;
const  Int_t linecolorEWK = kOrange+10;
const  Int_t fillcolorEWK = kOrange+7;
const  Int_t linecolorQCD = kViolet+2;
const  Int_t fillcolorQCD = kViolet-5;
const  Int_t ratioColor   = kGray+2;

const Int_t    NBINS   = 50;
const Double_t METMIN  = 40;
const Double_t METMAX  = 140;

//=== MAIN MACRO ================================================================================================= 

void fitWlnu(   const TString  outputDir,   // output directory 
                const TString  ntupleDir,
                const TString  flav,
                const Double_t lumi        // integrated luminosity (/fb)'
) {
  gBenchmark->Start("fitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================     

  // file format for output plots
  const TString format("png"); 
  
  bool doMTCut = false;
  // use this to change if we use mT or MET
  // bool overwriteMET = true;
  bool doTemplate = true;
  // don't change this part
  bool doTransversMass = true;

  double yscale=0.5;

  // main set
  enum{no,cent,eta,keys,ru,rd,stat0,stat1,stat2,stat3,stat4,stat5,stat6,stat7,stat8,stat9};
  const string vaMET[]={"no","main","eta","keys","ru","rd","stat0","stat1","stat2","stat3","stat4","stat5","stat6","stat7","stat8","stat9"};
  int nMET = sizeof(vaMET)/sizeof(vaMET[0]);
  std::vector<string> vMET;
  for(int i = 0; i < nMET; i++){vMET.push_back(vaMET[i]);} // this is weird
  
  enum{main,mc,fsr,bkg,tagpt,effstat,pfireu,pfired};
  const string vaWeight[]={"eff","mc","fsr","bkg","tagpt","effstat","pfireu","pfired"};
  int nWeight = sizeof(vaWeight)/sizeof(vaWeight[0]);
  std::vector<string> vWeight;
  for(int i = 0; i < nWeight; i++){vWeight.push_back(vaWeight[i]);}
  
  std::cout << "size of weight array is " << nWeight << std::endl;
  std::cout << "size of met array is " << nMET << std::endl;
  
  Double_t vIsoBins[] = {0.0,0.15,0.30,0.45,0.55,0.65};
  int nIsoBins = sizeof(vIsoBins)/sizeof(vIsoBins[0])-1;
  std::cout << "size of isobin array is " << nIsoBins << std::endl;
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  
  const Double_t MT_CUT = 40.0;
  const Double_t mu_MASS = 0.1057;
  
  TH2D *hErr  = new TH2D("hErr", "",10,0,10,20,0,20);
  
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  // -----------------------------------------------------

   
  //
  // input ntuple file names
  //
  enum {eData, eWlnu, eZxx, eWx, eTtb, eDib, eQCD, eAntiData, eAntiWlnu, eAntiQCD, eAntiTtb, eAntiDib, eAntiWx, eAntiZxx };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  
  
  if( lumi < 250 ) {
    if(flav.CompareTo("Wenu") == 0 ){
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/we0_select.root"  ));  typev.push_back(eWlnu);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/we1_select.root"  ));  typev.push_back(eWlnu);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/we2_select.root"  ));  typev.push_back(eWlnu);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx0_select.root"  ));  typev.push_back(eWx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx1_select.root"  ));  typev.push_back(eWx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx2_select.root"  ));  typev.push_back(eWx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zxx_select.root"  ));  typev.push_back(eZxx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zz_select.root"   ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/ww_select.root"   ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wz_select.root"   ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top1_select.root" ));  typev.push_back(eTtb);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top2_select.root" ));  typev.push_back(eTtb);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top3_select.root" ));  typev.push_back(eTtb);

      fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/we0_select.root")); typev.push_back(eAntiWlnu);
      fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/we1_select.root")); typev.push_back(eAntiWlnu);
      fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/we2_select.root")); typev.push_back(eAntiWlnu);
      
    } else if(flav.CompareTo("Wmunu") == 0 ){
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wm0_select.raw.root" ));  typev.push_back(eWlnu);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wm1_select.raw.root" ));  typev.push_back(eWlnu);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wm2_select.raw.root" ));  typev.push_back(eWlnu);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx0_select.raw.root" ));  typev.push_back(eWx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx1_select.raw.root" ));  typev.push_back(eWx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx2_select.raw.root" ));  typev.push_back(eWx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zxx_select.raw.root" ));  typev.push_back(eZxx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zz_select.raw.root"  ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/ww_select.raw.root"  ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wz_select.raw.root"  ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top1_select.raw.root"));  typev.push_back(eTtb);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top2_select.raw.root"));  typev.push_back(eTtb);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top3_select.raw.root"));  typev.push_back(eTtb);

      fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wm0_select.root")); typev.push_back(eAntiWlnu);
      fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wm1_select.root")); typev.push_back(eAntiWlnu);
      fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wm2_select.root")); typev.push_back(eAntiWlnu);
    }
    

    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/data_select.root"    ));  typev.push_back(eData);
      
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/data_select.root"));  typev.push_back(eAntiData);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx0_select.root" ));  typev.push_back(eAntiWx);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx1_select.root" ));  typev.push_back(eAntiWx);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx2_select.root" ));  typev.push_back(eAntiWx);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zxx_select.root" ));  typev.push_back(eAntiZxx);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/ww_select.root"  ));  typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wz_select.root"  ));  typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zz_select.root"  ));  typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top1_select.root"));  typev.push_back(eAntiTtb);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top2_select.root"));  typev.push_back(eAntiTtb);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top3_select.root"));  typev.push_back(eAntiTtb);
  }
 
  
  if( lumi > 250 ){ // The 5 TeV
    if(flav.CompareTo("Wmunu") == 0){
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wm_select.raw.root"  )); typev.push_back(eWlnu);    
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/data_select.root"    ));  typev.push_back(eData);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx_select.raw.root"  ));  typev.push_back(eWx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zxx_select.raw.root" ));  typev.push_back(eZxx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zz_select.raw.root"  ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/ww_select.raw.root"  ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wz_select.raw.root"  ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top_select.raw.root" ));  typev.push_back(eTtb);
      
      fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wm_select.root"    )); typev.push_back(eAntiWlnu);
    }else if(flav.CompareTo("Wenu") == 0){
      fnamev.push_back(ntupleDir+TString("/") +flav+TString("/ntuples/we_select.root"  ));  typev.push_back(eWlnu);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/data_select.root" ));  typev.push_back(eData);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx_select.root"   ));  typev.push_back(eWx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zxx_select.root"  ));  typev.push_back(eZxx);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zz_select.root"   ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/ww_select.root"   ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wz_select.root"   ));  typev.push_back(eDib);
      fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top_select.root"  ));  typev.push_back(eTtb);
      
      fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/we_select.root" ));  typev.push_back(eAntiWlnu);
    }
  
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/data_select.root")); typev.push_back(eAntiData);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx_select.root"  )); typev.push_back(eAntiWx);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zxx_select.root" )); typev.push_back(eAntiZxx);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/ww_select.root"  )); typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wz_select.root"  )); typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zz_select.root"  )); typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top_select.root" )); typev.push_back(eAntiTtb);
  }

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  
  // All teh uncertainties? 
  TH1D ***hWlnupMETU  = new TH1D**[nIsoBins];   TH1D ***hWlnumMETU  = new TH1D**[nIsoBins]; 
  TH1D ***hWlnupMETD  = new TH1D**[nIsoBins];   TH1D ***hWlnumMETD  = new TH1D**[nIsoBins];
  
  TH1D ***hEWKpMETU  = new TH1D**[nIsoBins];  TH1D ***hEWKmMETU  = new TH1D**[nIsoBins];
  TH1D ***hEWKpMETD  = new TH1D**[nIsoBins];  TH1D ***hEWKmMETD  = new TH1D**[nIsoBins];
  
  TH1D ***hWxpMETU  = new TH1D**[nIsoBins];  TH1D ***hWxmMETU  = new TH1D**[nIsoBins];
  TH1D ***hWxpMETD  = new TH1D**[nIsoBins];  TH1D ***hWxmMETD  = new TH1D**[nIsoBins];
  
  TH1D ***hZxxpMETU  = new TH1D**[nIsoBins];  TH1D ***hZxxmMETU  = new TH1D**[nIsoBins];
  TH1D ***hZxxpMETD  = new TH1D**[nIsoBins];  TH1D ***hZxxmMETD  = new TH1D**[nIsoBins];
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D ***hWlnupWeightU  = new TH1D**[nIsoBins];  TH1D ***hWlnumWeightU  = new TH1D**[nIsoBins];
  TH1D ***hWlnupWeightD  = new TH1D**[nIsoBins];  TH1D ***hWlnumWeightD  = new TH1D**[nIsoBins];
  
  TH1D ***hEWKpWeightU   = new TH1D**[nIsoBins];  TH1D ***hEWKmWeightU   = new TH1D**[nIsoBins];
  TH1D ***hEWKpWeightD   = new TH1D**[nIsoBins];  TH1D ***hEWKmWeightD   = new TH1D**[nIsoBins];
  
  TH1D ***hWxpWeightU    = new TH1D**[nIsoBins];  TH1D ***hWxmWeightU    = new TH1D**[nIsoBins];
  TH1D ***hWxpWeightD    = new TH1D**[nIsoBins];  TH1D ***hWxmWeightD    = new TH1D**[nIsoBins];
  
  TH1D ***hZxxpWeightU   = new TH1D**[nIsoBins];  TH1D ***hZxxmWeightU   = new TH1D**[nIsoBins];
  TH1D ***hZxxpWeightD   = new TH1D**[nIsoBins];  TH1D ***hZxxmWeightD   = new TH1D**[nIsoBins];
  
  TH1D ***hDibpWeightU   = new TH1D**[nIsoBins];  TH1D ***hDibmWeightU   = new TH1D**[nIsoBins];
  TH1D ***hDibpWeightD   = new TH1D**[nIsoBins];  TH1D ***hDibmWeightD   = new TH1D**[nIsoBins];
  
  TH1D ***hTtbpWeightU   = new TH1D**[nIsoBins];  TH1D ***hTtbmWeightU   = new TH1D**[nIsoBins];
  TH1D ***hTtbpWeightD   = new TH1D**[nIsoBins];  TH1D ***hTtbmWeightD   = new TH1D**[nIsoBins];
  
  // do a loop to create the second array
  for(int i=0; i < nIsoBins; ++i){
    // w signal recoil
    hWlnupMETU[i] = new TH1D*[nMET];    hWlnumMETU[i] = new TH1D*[nMET];
    hWlnupMETD[i] = new TH1D*[nMET];    hWlnumMETD[i] = new TH1D*[nMET];
    
    // ewk total recoil
    hEWKpMETU[i] = new TH1D*[nMET];    hEWKmMETU[i] = new TH1D*[nMET];
    hEWKpMETD[i] = new TH1D*[nMET];    hEWKmMETD[i] = new TH1D*[nMET];
    
    // wx  recoil
    hWxpMETU[i] = new TH1D*[nMET];    hWxmMETU[i] = new TH1D*[nMET];
    hWxpMETD[i] = new TH1D*[nMET];    hWxmMETD[i] = new TH1D*[nMET];
    
    // zxx recoil
    hZxxpMETU[i] = new TH1D*[nMET];    hZxxmMETU[i] = new TH1D*[nMET];
    hZxxpMETD[i] = new TH1D*[nMET];    hZxxmMETD[i] = new TH1D*[nMET];
    
    // w signal efficiency 
    hWlnupWeightU[i] = new TH1D*[nWeight];   hWlnumWeightU[i] = new TH1D*[nWeight];
    hWlnupWeightD[i] = new TH1D*[nWeight];   hWlnumWeightD[i] = new TH1D*[nWeight];
    
    // ewk total efficiency
    hEWKpWeightU[i]  = new TH1D*[nWeight];   hEWKmWeightU[i] = new TH1D*[nWeight];
    hEWKpWeightD[i]  = new TH1D*[nWeight];   hEWKmWeightD[i] = new TH1D*[nWeight];
    
    // wx  efficiency
    hWxpWeightU[i]   = new TH1D*[nWeight];   hWxmWeightU[i] = new TH1D*[nWeight];
    hWxpWeightD[i]   = new TH1D*[nWeight];   hWxmWeightD[i] = new TH1D*[nWeight];
    
    // zxx efficiency
    hZxxpWeightU[i]  = new TH1D*[nWeight];   hZxxmWeightU[i] = new TH1D*[nWeight];
    hZxxpWeightD[i]  = new TH1D*[nWeight];   hZxxmWeightD[i] = new TH1D*[nWeight];

    
    // diboson efficiency
    hDibpWeightU[i] = new TH1D*[nWeight];    hDibmWeightU[i] = new TH1D*[nWeight];
    hDibpWeightD[i] = new TH1D*[nWeight];    hDibmWeightD[i] = new TH1D*[nWeight];
    
    // ttbar efficiency
    hTtbpWeightU[i] = new TH1D*[nWeight];    hTtbmWeightU[i] = new TH1D*[nWeight];
    hTtbpWeightD[i] = new TH1D*[nWeight];    hTtbmWeightD[i] = new TH1D*[nWeight];
    
  }
  
  // -----------------------------------------------------------------------
  //           the main shapes  
  // -----------------------------------------------------------------------
  TH1D **hDataMetm2d  = new TH1D*[nIsoBins]; 
  TH1D **hDataMetp2d  = new TH1D*[nIsoBins];
  
  TH1D **hWlnuMetp2d  = new TH1D*[nIsoBins];
  TH1D **hWlnuMetm2d  = new TH1D*[nIsoBins];
  
  TH1D **hEWKMetp2d   = new TH1D*[nIsoBins];
  TH1D **hEWKMetm2d   = new TH1D*[nIsoBins];
  
  TH1D **hDibMetp2d   = new TH1D*[nIsoBins];
  TH1D **hDibMetm2d   = new TH1D*[nIsoBins];
  
  TH1D **hTtbMetp2d   = new TH1D*[nIsoBins];
  TH1D **hTtbMetm2d   = new TH1D*[nIsoBins];
  
  TH1D **hWxMetp2d   = new TH1D*[nIsoBins];
  TH1D **hWxMetm2d   = new TH1D*[nIsoBins];
  
  TH1D **hZxxMetp2d   = new TH1D*[nIsoBins];
  TH1D **hZxxMetm2d   = new TH1D*[nIsoBins];

  TH1D **hQCDMetp2d   = new TH1D*[nIsoBins];
  TH1D **hQCDMetm2d   = new TH1D*[nIsoBins];
  
  
  TH1D **hMetpIsoValues = new TH1D*[nIsoBins];
  TH1D **hMetmIsoValues = new TH1D*[nIsoBins];
  // Create a histogram pointer in each space in the array
  for(int i = 0; i < nIsoBins; i++){
    hDataMetp2d[i]  = new TH1D(("hDataMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hDataMetm2d[i]  = new TH1D(("hDataMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    hWlnuMetp2d[i]  = new TH1D(("hWlnuMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hWlnuMetm2d[i]  = new TH1D(("hWlnuMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    hEWKMetp2d[i]   = new TH1D(("hEwkMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hEWKMetm2d[i]   = new TH1D(("hEwkMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    hDibMetp2d[i]   = new TH1D(("hDibMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hDibMetm2d[i]   = new TH1D(("hDibMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    hTtbMetp2d[i]   = new TH1D(("hTtbMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hTtbMetm2d[i]   = new TH1D(("hTtbMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
	
    hWxMetp2d[i]    = new TH1D(("hWxMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hWxMetm2d[i]    = new TH1D(("hWxMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
   
    hZxxMetp2d[i]   = new TH1D(("hZxxMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hZxxMetm2d[i]   = new TH1D(("hZxxMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    hQCDMetp2d[i]   = new TH1D(("hQCDMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hQCDMetm2d[i]   = new TH1D(("hQCDMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    // Create a loop over the # of uncertainty shapes to produce the uncertainty histograms for the "up" shapes
    // Do the recoil ones here
    for(int j = 0; j < nMET; ++j){
      char hname[150]; char type[50];
      if      (j==ru) sprintf(type,"%i_rochUp",i);
      else if (j==rd) sprintf(type,"%i_rochDown",i);
      else            sprintf(type,"%i_%sUp",i,(vMET[j]).c_str());
      

      // Wlnu
      sprintf(hname,"hWlnuMetpBin%s",type);   hWlnupMETU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hWlnuMetmBin%s",type);   hWlnumMETU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      // ewk sum
      sprintf(hname,"hEWKMetpBin%s", type);   hEWKpMETU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hEWKMetmBin%s", type);   hEWKmMETU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      // wx
      sprintf(hname,"hWxMetpBin%s",  type);   hWxpMETU[i][j]   = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hWxMetmBin%s",  type);   hWxmMETU[i][j]   = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      // zxx
      sprintf(hname,"hZxxMetpBin%s", type);   hZxxpMETU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hZxxMetmBin%s", type);   hZxxmMETU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
      // }
    }
    
    for(int j=0; j < nWeight; ++j){
      char hname[150]; char type[50];
      //sdu,sdd,smu,smd,pfireu,pfired
      if      (j==pfireu) sprintf(type,"%i_prefireUp",i);
      else if (j==pfired) sprintf(type,"%i_prefireDown",i);
      else                sprintf(type,"%i_%sUp",i,(vWeight[j]).c_str());
      
      
      sprintf(hname,"hWlnuMetpBin%s",type);  hWlnupWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hWlnuMetmBin%s",type);  hWlnumWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      sprintf(hname,"hEwkMetpBin%s", type);  hEWKpWeightU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hEwkMetmBin%s", type);  hEWKmWeightU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      sprintf(hname,"hWxMetpBin%s",  type);  hWxpWeightU[i][j]   = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hWxMetmBin%s",  type);  hWxmWeightU[i][j]   = new TH1D(hname,"",NBINS,METMIN,METMAX);
  
      sprintf(hname,"hZxxMetpBin%s", type);  hZxxpWeightU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hZxxMetmBin%s", type);  hZxxmWeightU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      sprintf(hname,"hDibMetpBin%s", type);  hDibpWeightU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hDibMetmBin%s", type);  hDibmWeightU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
  
      sprintf(hname,"hTtbMetpBin%s", type);  hTtbpWeightU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hTtbMetmBin%s", type);  hTtbmWeightU[i][j]  = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
    }
    
    hMetpIsoValues[i] = new TH1D(("hMetpIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],vIsoBins[i+1]);
    hMetmIsoValues[i] = new TH1D(("hMetmIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],vIsoBins[i+1]);
  }
  

  // Some test plots: 
  // pfiso/pt vs rel comb iso in barrel
  // in endcap
  TH2D *hBarrelIsoPosData = new TH2D("hBarrelIsoPosData","hBarrelIsoPosData",100,25,100,100,0,0.50);
  TH2D *hEndcapIsoPosData = new TH2D("hEndcapIsoPosData","hEndcapIsoPosData",100,25,100,100,0,0.50);
  
  TH2D *hBarrelIsoPosWsig = new TH2D("hBarrelIsoPosWsig","hBarrelIsoPosWsig",100,25,100,100,0,0.50);
  TH2D *hEndcapIsoPosWsig = new TH2D("hEndcapIsoPosWsig","hEndcapIsoPosWsig",100,25,100,100,0,0.50);
  
  TH2D *hBarrelIsoNegData = new TH2D("hBarrelIsoNegData","hBarrelIsoNegData",100,25,100,100,0,0.50);
  TH2D *hEndcapIsoNegData = new TH2D("hEndcapIsoNegData","hEndcapIsoNegData",100,25,100,100,0,0.50);
  
  TH2D *hBarrelIsoNegWsig = new TH2D("hBarrelIsoNegWsig","hBarrelIsoNegWsig",100,25,100,100,0,0.50);
  TH2D *hEndcapIsoNegWsig = new TH2D("hEndcapIsoNegWsig","hEndcapIsoNegWsig",100,25,100,100,0,0.50);

  double tolerance = ROOT::Math::MinimizerOptions::DefaultTolerance();
  string algo = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();
  string type = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
  int strategy= ROOT::Math::MinimizerOptions::DefaultStrategy();

  int precision= ROOT::Math::MinimizerOptions::DefaultPrecision();
  int MaxFunctionCalls= ROOT::Math::MinimizerOptions::DefaultMaxFunctionCalls();
  int MaxIterations= ROOT::Math::MinimizerOptions::DefaultMaxIterations();

  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum, npv, npu;
  Float_t genVPt, genVPhi, genVy, genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown, prefireWeight, prefirePhoton, prefireJet;
  Float_t prefireUp, prefireDown;
  Float_t met, metPhi, sumEt, mt, u1, u2, deta;
  Int_t   q;
  TLorentzVector *lep=0, *lep_raw=0, *genV=0, *genLep=0;
  Float_t pfChIso, pfGamIso, pfNeuIso, pfCombIso, trkIso;
  Double_t mtCorr;
  Double_t effSFweight=1, relIso;


  vector<Double_t>  *metVars=0, *metVarsPhi=0;
  vector<Double_t>  *evtWeight=0;
    
  
  TFile *infile=0;
  TTree *intree=0;

  double noPrefire_Wp = 0, prefire_Wp=0, prefireJet_Wp=0, prefirePhoton_Wp=0;
  double noPrefire_Wm = 0, prefire_Wm=0, prefireJet_Wm=0, prefirePhoton_Wm=0;

  //
  // Loop over files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);    assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);
    std::cout << "type = " << typev[ifile] << std::endl;

    intree->SetBranchAddress("runNum",   &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    intree->SetBranchAddress("npv",      &npv);       // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);       // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genVy",    &genVy);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genLepPt",   &genLepPt);    // GEN lepton pT (signal MC)
    intree->SetBranchAddress("genLepPhi",  &genLepPhi);   // GEN lepton phi (signal MC)
    intree->SetBranchAddress("prefireWeight",  &prefireWeight);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("prefireUp",      &prefireUp);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("prefireDown",    &prefireDown);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("prefirePhoton",  &prefirePhoton);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("prefireJet",    &prefireJet);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("scale1fb",      &scale1fb);  // MCwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp",    &scale1fbUp);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown",  &scale1fbDown);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("sumEt",         &sumEt);     // Sum ET
    intree->SetBranchAddress("mt",            &mt);        // transverse mass
    intree->SetBranchAddress("mtCorr",        &mtCorr);        // transverse mass
    intree->SetBranchAddress("q",             &q);         // lepton charge
    intree->SetBranchAddress("lep",           &lep);       // lepton 4-vector
    intree->SetBranchAddress("lep_raw",       &lep_raw);       // lepton 4-vector
    intree->SetBranchAddress("genLep",        &genLep);       // lepton 4-vector
    intree->SetBranchAddress("genV",          &genV);       // lepton 4-vector
    intree->SetBranchAddress("trkIso",        &trkIso);
    intree->SetBranchAddress("pfChIso",       &pfChIso);
    intree->SetBranchAddress("pfGamIso",      &pfGamIso);
    intree->SetBranchAddress("pfNeuIso",      &pfNeuIso);
    intree->SetBranchAddress("pfCombIso",     &pfCombIso);       // lepton 4-vector
    intree->SetBranchAddress("relIso",        &relIso);       // relative isolation for the lepton
    intree->SetBranchAddress("evtWeight",     &evtWeight); // eventwgtLum[main] vector
    intree->SetBranchAddress("metVars",       &metVars);            // contains the different met variations
    intree->SetBranchAddress("metVarsPhi",    &metVarsPhi);         // met phi for variations
  
    TH1D* hGenWeights; double totalNorm = 1.0;
    // cout << "Hello " << endl;
    // if(typev[ifile] != eData && typev[ifile] != eAntiData && !(flav.CompareTo("Wmunu") == 0)){
    if(typev[ifile] != eData && typev[ifile] != eAntiData){
      // cout << "get gen weights" << endl;
      hGenWeights = (TH1D*)infile->Get("hGenWeights");
      totalNorm = hGenWeights->Integral();
      // cout << totalNorm << endl;
    }
  
    UInt_t iterator=20;
    std::cout << "Number of Events = " << intree->GetEntries() << std::endl;
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    // for(UInt_t ientry=0; ientry<((UInt_t)intree->GetEntries()); ientry+=iterator) {
      intree->GetEntry(ientry);
      if(ientry%1000000==0)  cout << "Event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " % done with this file." << endl;

      // Lepton Kinematic cut for the Nth time...
      if(lep->Pt()        < PT_CUT) continue;
      if(fabs(lep->Eta()) > ETA_CUT) continue;

      vector<double> wgtLum;

      // accidentally set up prefire uncertainties in a way that a very rare event gets nan for the up/down variation
      for(int jt=0; jt < nWeight; jt++) {
        wgtLum.push_back(lumi*(isnan((*evtWeight)[jt]) ? (*evtWeight)[main] : (*evtWeight)[jt])/totalNorm);
      }
      
      if(doTransversMass && typev[ifile] != eData && typev[ifile] != eAntiData){
        makeUncMT(*metVars, *metVarsPhi, lep);
      } else if (doTransversMass) (*metVars)[no] = mtCorr;

      if(doMTCut&&(mtCorr<MT_CUT)) continue;//std::cout << " pass mt " << std::endl;
      if(typev[ifile]==eData) {
        if(q>0) {
          hDataMetp2d[0]->Fill((*metVars)[no]);
          hQCDMetp2d[0] ->Fill((*metVars)[no]);
          hMetpIsoValues[0]->Fill(relIso);
        } else {
          hDataMetm2d[0]->Fill((*metVars)[no]);
          hQCDMetm2d[0] ->Fill((*metVars)[no]); 
          hMetmIsoValues[0]->Fill(relIso);
        }
      } else if(typev[ifile]==eAntiData) {
        if(abs(lep->Eta())<1.4442){
          if(q > 0)hBarrelIsoPosData->Fill(lep->Pt(),pfCombIso/lep->Pt());
          else hBarrelIsoNegData->Fill(lep->Pt(),relIso);
        } else {
          if(q > 0)hEndcapIsoPosData->Fill(lep->Pt(),pfCombIso/lep->Pt());
          else hEndcapIsoNegData->Fill(lep->Pt(),relIso);
        }
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0) { 
              hDataMetp2d[it]->Fill((*metVars)[no]);
              hQCDMetp2d[it] ->Fill((*metVars)[no]);
              hMetpIsoValues[it]->Fill(relIso);
            } else { 
              hDataMetm2d[it]->Fill((*metVars)[no]);
              hQCDMetm2d[it] ->Fill((*metVars)[no]);
              hMetmIsoValues[it]->Fill(relIso);
              break;
            }
          }
        }
      } else if(typev[ifile]==eWlnu) {
        int bin=0;
        if(typev[ifile]==eWlnu){
          if(q > 0){
            noPrefire_Wp+=wgtLum[main]/prefireWeight;
            prefire_Wp+=wgtLum[main];
            prefireJet_Wp+=wgtLum[main]*prefireJet/prefireWeight;
            prefirePhoton_Wp+=wgtLum[main]*prefirePhoton/prefireWeight;
          } else {
            noPrefire_Wm+=wgtLum[main]/prefireWeight;
            prefire_Wm+=wgtLum[main];
            prefireJet_Wm+=wgtLum[main]*prefireJet/prefireWeight;
            prefirePhoton_Wm+=wgtLum[main]*prefirePhoton/prefireWeight;
          }
          
        }
        
        if(q>0){
          hWlnuMetp2d[0]->Fill((*metVars)[cent],wgtLum[main]);
          fillMETs(hWlnupMETU[0],(*metVars),nMET,wgtLum[main]);
          fillWeights(hWlnupWeightU[0],(*metVars)[cent],nWeight,wgtLum);

        } else {
          hWlnuMetm2d[0]->Fill((*metVars)[cent],wgtLum[main]);
          fillMETs(hWlnumMETU[0],(*metVars),nMET,wgtLum[main]);
          fillWeights(hWlnumWeightU[0],(*metVars)[cent],nWeight,wgtLum);
        }
      } else if(typev[ifile]==eWx) {
        if(q>0){
          hWxMetp2d[0]->Fill((*metVars)[cent],wgtLum[main]);
          fillMETs(hWxpMETU[0],(*metVars),nMET,wgtLum[main]);
          fillWeights(hWxpWeightU[0],(*metVars)[cent],nWeight,wgtLum);
        } else {
          hWxMetm2d[0]->Fill((*metVars)[cent],wgtLum[main]);
          fillMETs(hWxmMETU[0],(*metVars),nMET,wgtLum[main]);
          fillWeights(hWxmWeightU[0],(*metVars)[cent],nWeight,wgtLum);
        }
      } else if(typev[ifile]==eZxx){
        if(q>0){
          hZxxMetp2d[0]->Fill((*metVars)[cent],wgtLum[main]);
          fillMETs(hZxxpMETU[0],(*metVars),nMET,wgtLum[main]);
          fillWeights(hZxxpWeightU[0],(*metVars)[cent],nWeight,wgtLum);
        } else {
          hZxxMetm2d[0]->Fill((*metVars)[cent],wgtLum[main]);
          fillMETs(hZxxmMETU[0],(*metVars),nMET,wgtLum[0]);
          fillWeights(hZxxmWeightU[0],(*metVars)[cent],nWeight,wgtLum);
        }
      } else if(typev[ifile]==eDib) {
        if(q>0){
          hDibMetp2d[0]->Fill((*metVars)[no],wgtLum[main]);
          fillWeights(hDibpWeightU[0],(*metVars)[no],nWeight,wgtLum);
        } else {
          hDibMetm2d[0]->Fill((*metVars)[no],wgtLum[main]);
          fillWeights(hDibmWeightU[0],(*metVars)[no],nWeight,wgtLum);
        }
      } else if(typev[ifile]==eTtb) {
        if(q>0){
          hTtbMetp2d[0]->Fill((*metVars)[no],wgtLum[main]);
          fillWeights(hTtbpWeightU[0],(*metVars)[no],nWeight,wgtLum);
        } else {
          hTtbMetm2d[0]->Fill((*metVars)[no],wgtLum[main]);
          fillWeights(hTtbmWeightU[0],(*metVars)[no],nWeight,wgtLum);
        }
      } else if(typev[ifile]==eAntiWlnu){
        if(abs(lep->Eta())<1.4442){
          if(q > 0)hBarrelIsoPosWsig->Fill(lep->Pt(),pfCombIso/lep->Pt());
          else hBarrelIsoNegWsig->Fill(lep->Pt(),relIso);
        } else {
          if(q > 0)hEndcapIsoPosWsig->Fill(lep->Pt(),pfCombIso/lep->Pt());
          else hEndcapIsoNegWsig->Fill(lep->Pt(),relIso);
        }
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0) {              
              hWlnuMetp2d[it]->Fill((*metVars)[no],wgtLum[main]);
              fillMETs(hWlnupMETU[it],(*metVars),nMET,wgtLum[main]);
              fillWeights(hWlnupWeightU[it],(*metVars)[no],nWeight,wgtLum);
            } else {
              hWlnuMetm2d[it]->Fill((*metVars)[no],wgtLum[main]);
              fillMETs(hWlnumMETU[it],(*metVars),nMET,wgtLum[main]);
              fillWeights(hWlnumWeightU[it],(*metVars)[no],nWeight,wgtLum);
            }
          }
        }
      } else if(typev[ifile]==eAntiWx){
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hWxMetp2d[it]->Fill((*metVars)[no],wgtLum[main]);
              fillMETs(hWxpMETU[it],(*metVars),nMET,wgtLum[main]);
              fillWeights(hWxpWeightU[it],(*metVars)[no],nWeight,wgtLum);
            } else {
              hWxMetm2d[it]->Fill((*metVars)[no],wgtLum[main]);
              fillMETs(hWxmMETU[it],(*metVars),nMET,wgtLum[main]);
              fillWeights(hWxmWeightU[it],(*metVars)[no],nWeight,wgtLum);
            }
            break;
          }
        }
      } else if(typev[ifile]==eAntiZxx){
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hZxxMetp2d[it]->Fill((*metVars)[no],wgtLum[main]);
              fillMETs(hZxxpMETU[it],(*metVars),nMET,wgtLum[main]);
              fillWeights(hZxxpWeightU[it],(*metVars)[no],nWeight,wgtLum);
            } else {
              hZxxMetm2d[it]->Fill((*metVars)[no],wgtLum[main]);
              fillMETs(hZxxmMETU[it],(*metVars),nMET,wgtLum[main]);
              fillWeights(hZxxmWeightU[it],(*metVars)[no],nWeight,wgtLum);
            }
          }
        }
      } else if(typev[ifile]==eAntiDib){
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hDibMetp2d[it]->Fill((*metVars)[no],wgtLum[main]);
              fillWeights(hDibpWeightU[it],(*metVars)[no],nWeight,wgtLum);
            } else {
              hDibMetm2d[it]->Fill((*metVars)[no],wgtLum[main]);
              fillWeights(hDibmWeightU[it],(*metVars)[no],nWeight,wgtLum);
            }
          }
        }
      } else if(typev[ifile]==eAntiTtb){
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hTtbMetp2d[it]->Fill((*metVars)[no],wgtLum[main]);
              fillWeights(hTtbpWeightU[it],(*metVars)[no],nWeight,wgtLum);
            } else {
              hTtbMetm2d[it]->Fill((*metVars)[no],wgtLum[main]);
              fillWeights(hTtbmWeightU[it],(*metVars)[no],nWeight,wgtLum);
            }
          }
        }
      } 
    }
  }
  
  
  // Print the yields in each of the isolation bins
  ofstream txtIso;
  char txtnameIso[100];
  std::cout << "Printing Yields by Isolation Bin: " << std::endl;
  sprintf(txtnameIso,"%s/isoBins_Yields.txt",CPlot::sOutDir.Data());
  txtIso.open(txtnameIso);
  assert(txtIso.is_open());
	txtIso << "W+ " << std::endl << "data" << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtIso << "bin " << i << "  Data# "  << hDataMetp2d[i]->Integral() << std::endl;
		txtIso << "bin " << i << "  Wsig# "  << hWlnuMetp2d[i]->Integral() << std::endl;
		// txtIso << "bin " << i << "  ewk # "  << hEWKMetp2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "   wx # "  << hWxMetp2d[i]->Integral()   << std::endl;
		txtIso << "bin " << i << "  zxx # "  << hZxxMetp2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "  ttb # "  << hTtbMetp2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "  dib # "  << hDibMetp2d[i]->Integral()  << std::endl;
	}
	
	txtIso << "W- " << std::endl << "data" << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtIso << "bin " << i << "  Data# "  << hDataMetm2d[i]->Integral() << std::endl;
		txtIso << "bin " << i << "  Wsig# "  << hWlnuMetm2d[i]->Integral() << std::endl;
		// txtIso << "bin " << i << "  ewk # "  << hEWKMetm2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "   wx # "  << hWxMetm2d[i]->Integral()   << std::endl;
		txtIso << "bin " << i << "  zxx # "  << hZxxMetm2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "  ttb # "  << hTtbMetm2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "  dib # "  << hDibMetm2d[i]->Integral()  << std::endl;
	}
	txtIso.close();
  
  
  char sname[50];
  for(int it = 0; it < nIsoBins; it++){
    
    hEWKMetp2d[it]->Add(hTtbMetp2d[it],1);
    hEWKMetp2d[it]->Add(hDibMetp2d[it],1);
    hEWKMetp2d[it]->Add(hZxxMetp2d[it],1);
    hEWKMetp2d[it]->Add(hWxMetp2d[it],1);
    
    hEWKMetm2d[it]->Add(hTtbMetm2d[it],1);
    hEWKMetm2d[it]->Add(hDibMetm2d[it],1);
    hEWKMetm2d[it]->Add(hZxxMetm2d[it],1);
    hEWKMetm2d[it]->Add(hWxMetm2d[it],1);

    hQCDMetp2d[it]->Add(hEWKMetp2d[it],-1);
    hQCDMetp2d[it]->Add(hWlnuMetp2d[it],-1);
	
    hQCDMetm2d[it]->Add(hEWKMetm2d[it],-1);
    hQCDMetm2d[it]->Add(hWlnuMetm2d[it],-1);

    for(Int_t ibin=1; ibin<=hQCDMetp2d[it]->GetNbinsX(); ibin++) {
      if(hQCDMetp2d[it]->GetBinContent(ibin)<0) {
        hQCDMetp2d[it]->SetBinContent(ibin,0);
        hQCDMetp2d[it]->SetBinError(ibin,1.8);
      }
      if(hQCDMetm2d[it]->GetBinContent(ibin)<0) {
        hQCDMetm2d[it]->SetBinContent(ibin,0);
        hQCDMetm2d[it]->SetBinError(ibin,1.8);
      }
    }
   
    
    for(int j=0; j< nWeight; j++){
      hEWKpWeightU[it][j]->Add(hTtbpWeightU[it][j],1);
      hEWKpWeightU[it][j]->Add(hDibpWeightU[it][j],1);
      hEWKpWeightU[it][j]->Add(hZxxpWeightU[it][j],1);
      hEWKpWeightU[it][j]->Add(hWxpWeightU[it][j],1);
      
      hEWKmWeightU[it][j]->Add(hTtbmWeightU[it][j],1);
      hEWKmWeightU[it][j]->Add(hDibmWeightU[it][j],1);
      hEWKmWeightU[it][j]->Add(hZxxmWeightU[it][j],1);
      hEWKmWeightU[it][j]->Add(hWxmWeightU[it][j],1);
    }
    for(int j=0; j < nMET; j++){
      hEWKpMETU[it][j]->Add(hZxxpMETU[it][j],1);
      hEWKpMETU[it][j]->Add(hWxpMETU[it][j],1);
      
      hEWKmMETU[it][j]->Add(hZxxmMETU[it][j],1);
      hEWKmMETU[it][j]->Add(hWxmMETU[it][j],1);
    }
    
  }
  // std::cout << "blahasfadsfa " << std::endl;
  for(int it=0; it<1; it++){
    
    for(int i=0; i < nWeight; i++){
      std::cout << hWlnupWeightU[it][i]->Integral() << std::endl;
    }
    for(int i=0; i < nMET; i++){
      std::cout << hWlnupMETU[it][i]->Integral() << std::endl;
    }
  }
  
  delete infile;
  infile=0, intree=0;   
   
  ofstream txtfilePF;
  char txtfname2[150];
  std::cout << "Printing prefiring  values" << std::endl;
  sprintf(txtfname2,"%s/prefire_yields.txt",CPlot::sOutDir.Data());
  txtfilePF.open(txtfname2);
  assert(txtfilePF.is_open());
  txtfilePF << "---- W+ ----- " << endl;
  txtfilePF << "No PF: " << noPrefire_Wp << endl;
  txtfilePF << "PF: " << prefire_Wp << "  scale factor: " << noPrefire_Wp/prefire_Wp << endl;
  txtfilePF << "PF (jet only): " << prefireJet_Wp << "  scale factor: " << noPrefire_Wp/prefireJet_Wp << endl;
  txtfilePF << "PF (photon only): " << prefirePhoton_Wp << "  scale factor: " << noPrefire_Wp/prefirePhoton_Wp << endl << endl;
  txtfilePF << "---- W- ----- " << endl;
  txtfilePF << "No PF: " << noPrefire_Wm << endl;
  txtfilePF << "PF: " << prefire_Wm << "  scale factor: " << noPrefire_Wm/prefire_Wm << endl;
  txtfilePF << "PF (jet only): " << prefireJet_Wm << "  scale factor: " << noPrefire_Wm/prefireJet_Wm << endl;
  txtfilePF << "PF (photon only): " << prefirePhoton_Wm << "  scale factor: " << noPrefire_Wm/prefirePhoton_Wm << endl << endl;
  txtfilePF.close();
   
  ofstream txtfile2;
  // char txtfname2[100];
  std::cout << "Printing We+. " << std::endl;
  sprintf(txtfname2,"%s/isoAvg.txt",CPlot::sOutDir.Data());
  txtfile2.open(txtfname2);
  assert(txtfile2.is_open());
 
  char plotname[100];
  TCanvas *c2 = MakeCanvas("c2","c2",800,600);

  char nname[50];
  char formula[50];
  //  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                     up and down shapes
  //  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int j=0; j < nIsoBins; ++j){
    // W+ signal shape "down" uncertainty shapes
    TH1D *hh_diff;
    for(int k=0; k < nMET; ++k){
      if(k==ru||k==rd) continue;
      
      sprintf(nname,"hWlnuMetpBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hWlnupMETU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
      hWlnupMETD[j][k] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnupMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWlnuMetmBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hWlnumMETU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
      hWlnumMETD[j][k] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnumMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWEwkMetpBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hEWKpMETU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
      hEWKpMETD[j][k] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKpMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hEwkMetmBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hEWKmMETU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
      hEWKmMETD[j][k] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKmMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetpBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hWxpMETU[j][k]->Clone("diff"); hh_diff->Add(hWxMetp2d[j],-1);
      hWxpMETD[j][k] = (TH1D*) hWxMetp2d[j]->Clone(nname); hWxpMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetmBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hWxmMETU[j][k]->Clone("diff"); hh_diff->Add(hWxMetm2d[j],-1);
      hWxmMETD[j][k] = (TH1D*) hWxMetm2d[j]->Clone(nname); hWxmMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetpBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hZxxpMETU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetp2d[j],-1);
      hZxxpMETD[j][k] = (TH1D*) hZxxMetp2d[j]->Clone(nname); hZxxpMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetmBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hZxxmMETU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetm2d[j],-1);
      hZxxmMETD[j][k] = (TH1D*) hZxxMetm2d[j]->Clone(nname); hZxxmMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
    }
    
    for(int k=0; k < nWeight; ++k){
      if(k==pfireu||k==pfired)continue;
        // sprintf(type,"%i_prefireDown",i);

      sprintf(nname,"hWlnuMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hWlnupWeightU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
      hWlnupWeightD[j][k] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnupWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWlnuMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hWlnumWeightU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
      hWlnumWeightD[j][k] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnumWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hEwkMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hEWKpWeightU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
      hEWKpWeightD[j][k] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKpWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hEwkMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hEWKmWeightU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
      hEWKmWeightD[j][k] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKmWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hWxpWeightU[j][k]->Clone("diff"); hh_diff->Add(hWxMetp2d[j],-1);
      hWxpWeightD[j][k] = (TH1D*) hWxMetp2d[j]->Clone(nname); hWxpWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hWxmWeightU[j][k]->Clone("diff"); hh_diff->Add(hWxMetm2d[j],-1);
      hWxmWeightD[j][k] = (TH1D*) hWxMetm2d[j]->Clone(nname); hWxmWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hZxxpWeightU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetp2d[j],-1);
      hZxxpWeightD[j][k] = (TH1D*) hZxxMetp2d[j]->Clone(nname); hZxxpWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hZxxmWeightU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetm2d[j],-1);
      hZxxmWeightD[j][k] = (TH1D*) hZxxMetm2d[j]->Clone(nname); hZxxmWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hDibMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hDibpWeightU[j][k]->Clone("diff"); hh_diff->Add(hDibMetp2d[j],-1);
      hDibpWeightD[j][k] = (TH1D*) hDibMetp2d[j]->Clone(nname); hDibpWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hDibMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hDibmWeightU[j][k]->Clone("diff"); hh_diff->Add(hDibMetm2d[j],-1);
      hDibmWeightD[j][k] = (TH1D*) hDibMetm2d[j]->Clone(nname); hDibmWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hTtbMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hTtbpWeightU[j][k]->Clone("diff"); hh_diff->Add(hTtbMetp2d[j],-1);
      hTtbpWeightD[j][k] = (TH1D*) hTtbMetp2d[j]->Clone(nname); hTtbpWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hTtbMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hTtbmWeightU[j][k]->Clone("diff"); hh_diff->Add(hTtbMetm2d[j],-1);
      hTtbmWeightD[j][k] = (TH1D*) hTtbMetm2d[j]->Clone(nname); hTtbmWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
    }
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Fill the histograms
  //////////////////////////////////////////////////////////////////////////
  TH1 *qcdp2 = (TH1*) hQCDMetp2d[2]->Clone();   TH1 *qcdp3=(TH1*) hQCDMetp2d[3]->Clone();
  TH1 *qcdm2 = (TH1*) hQCDMetm2d[2]->Clone();   TH1 *qcdm3=(TH1*) hQCDMetm2d[3]->Clone();
  qcdp3->Scale(qcdp2->Integral()/qcdp3->Integral());
  qcdm3->Scale(qcdm2->Integral()/qcdm3->Integral());
  TH1D *hQCDMetpBin2_shapeUp   = (TH1D*) hQCDMetp2d[2]->Clone();  
  hQCDMetpBin2_shapeUp->SetTitle("hQCDMetpBin2_shapeUp");
  hQCDMetpBin2_shapeUp->SetName("hQCDMetpBin2_shapeUp");
  TH1D *hQCDMetpBin2_shapeDown = (TH1D*) hQCDMetp2d[2]->Clone();
  hQCDMetpBin2_shapeDown->SetTitle("hQCDMetpBin2_shapeDown");
  hQCDMetpBin2_shapeDown->SetName("hQCDMetpBin2_shapeDown");
  TH1D *hQCDMetmBin2_shapeUp   = (TH1D*) hQCDMetm2d[2]->Clone();
  hQCDMetmBin2_shapeUp->SetTitle("hQCDMetmBin2_shapeUp");
  hQCDMetmBin2_shapeUp->SetName("hQCDMetmBin2_shapeUp");
  TH1D *hQCDMetmBin2_shapeDown  = (TH1D*) hQCDMetm2d[2]->Clone();
  hQCDMetmBin2_shapeDown->SetTitle("hQCDMetmBin2_shapeDown");
  hQCDMetmBin2_shapeDown->SetName("hQCDMetmBin2_shapeDown");

  for(Int_t ibin=1; ibin<=qcdp2->GetNbinsX(); ibin++) 
    {
      double targetp = qcdp3->GetBinContent(ibin);
      double medianp = qcdp2->GetBinContent(ibin);
      double targetm = qcdm3->GetBinContent(ibin);
      double medianm = qcdm2->GetBinContent(ibin);
      hQCDMetpBin2_shapeUp->SetBinContent(ibin,targetp);
      hQCDMetmBin2_shapeUp->SetBinContent(ibin,targetm);
      hQCDMetpBin2_shapeDown->SetBinContent(ibin,2*medianp-targetp);
      hQCDMetmBin2_shapeDown->SetBinContent(ibin,2*medianm-targetm);
    }
  //create qcd up/down shapes

  TString histfname = outputDir + TString("/Wlnu_Histograms.root");
  TFile *histFile = new TFile(histfname,"RECREATE");
  histFile->cd();
  hQCDMetpBin2_shapeUp->Write();
  hQCDMetmBin2_shapeUp->Write();
  hQCDMetpBin2_shapeDown->Write();
  hQCDMetmBin2_shapeDown->Write();
  for(int j = 0; j < nIsoBins; ++j){
    hDataMetp2d[j]->Write();    hDataMetm2d[j]->Write();
    hWlnuMetp2d[j]->Write();    hWlnuMetm2d[j]->Write();
    hEWKMetp2d[j] ->Write();    hEWKMetm2d[j] ->Write();
    hZxxMetp2d[j] ->Write();    hZxxMetm2d[j] ->Write();
    hWxMetp2d[j]  ->Write();    hWxMetm2d[j]  ->Write();
    hTtbMetp2d[j] ->Write();    hTtbMetm2d[j] ->Write();
    hDibMetp2d[j] ->Write();    hDibMetm2d[j] ->Write();
    hQCDMetp2d[j] ->Write();    hQCDMetm2d[j] ->Write();

    for(int k=0; k < nMET; k++){
    // std::cout << "writing MET" << std::endl;
      hWlnupMETU[j][k]->Write();      hWlnumMETU[j][k]->Write();
      hEWKpMETU[j][k] ->Write();      hEWKmMETU[j][k] ->Write();
      hWxpMETU[j][k]  ->Write();      hWxmMETU[j][k]  ->Write();
      hZxxpMETU[j][k] ->Write();      hZxxmMETU[j][k] ->Write();
      
      // std::cout << "writing MET DOWN" << std::endl;
      if(k==rd||k==ru)continue;
      hWlnupMETD[j][k]->Write();     hWlnumMETD[j][k]->Write();
      hWxpMETD  [j][k]->Write();     hWxmMETD  [j][k]->Write();
      hEWKpMETD [j][k]->Write();     hEWKmMETD [j][k]->Write();
      hZxxpMETD [j][k]->Write();     hZxxmMETD [j][k]->Write();
      
    }
    
    for(int k=0; k < nWeight; k++){
      hWlnupWeightU[j][k]->Write();      hWlnumWeightU[j][k]->Write();
      hEWKpWeightU [j][k]->Write();      hEWKmWeightU [j][k]->Write();
      hWxpWeightU  [j][k]->Write();      hWxmWeightU  [j][k]->Write();
      hZxxpWeightU [j][k]->Write();      hZxxmWeightU [j][k]->Write();
      hDibpWeightU [j][k]->Write();      hDibmWeightU [j][k]->Write();
      hTtbpWeightU [j][k]->Write();      hTtbmWeightU [j][k]->Write();
      if(k==pfireu||k==pfired)continue;
            
      hWlnupWeightD[j][k]->Write();      hWlnumWeightD[j][k]->Write();
      hEWKpWeightD [j][k]->Write();      hEWKmWeightD [j][k]->Write();
      hWxpWeightD  [j][k]->Write();      hWxmWeightD  [j][k]->Write();
      hZxxpWeightD [j][k]->Write();      hZxxmWeightD [j][k]->Write();
      hDibpWeightD [j][k]->Write();      hDibmWeightD [j][k]->Write();
      hTtbpWeightD [j][k]->Write();      hTtbmWeightD [j][k]->Write();
    }
  }
  histFile->Write();
  histFile->Close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  gBenchmark->Show("fitWm");
}


//=== FUNCTION DEFINITIONS ======================================================================================


void fillMETs( TH1D** h, vector<double> met, int nMET, double wgt){
  for(int k =2; k < nMET; k++) h[k]->Fill(met[k] , wgt);
  return;
}

void fillWeights( TH1D** h, double met, int nWeight, vector<double> wgt ){
  for(int k =1; k < nWeight; k++) h[k] ->Fill(met , wgt[k]);
  return;
}

void makeUncMT(vector<Double_t> &metVars, vector<Double_t> &metVarsPhi, TLorentzVector* lep){
  for(int i = 0; i < metVars.size(); ++i){
    double mtCorr  = sqrt(2.0*(lep->Pt()) * (metVars[i]) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metVarsPhi[i]))));
    metVars[i] = mtCorr;
  }
  return;
}
