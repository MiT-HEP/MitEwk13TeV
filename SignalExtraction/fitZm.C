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

// helper class to handle efficiency tables
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"

#include "BaconAna/DataFormats/interface/TGenParticle.hh"  

#include "../Utils/MyTools.hh"	          // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
//#include "../Utils/RecoilCorrector.hh"    // class to handle recoil corrections for MET
//#include "../Utils/RecoilCorrector_htautau_hist.hh"
// #include "../Utils/RecoilCorrector_addJets.hh"
// #include "../Utils/RecoilCorrector_asymtest.hh"
#include "../Utils/RecoilCorrector_asym2.hh"
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections

// #include "ZBackgrounds.hh"

#include <rochcor2015r.h>
#include <muresolution_run2r.h>

// RooFit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);
//  TCanvas *b = MakeCanvas("b","b",800,600);
// print correlations of fitted parameters
void printCorrelations(ostream& os, RooFitResult *res);

// print chi2 test and KS test results
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
			   const Double_t ksprob, const Double_t ksprobpe);
               
// make webpage
void makeHTML(const TString outDir);
string int2string(int i) {
  stringstream ss;
  string ret;
  ss << i;
  ss >> ret;
  return ret;
}

std::string outputdir;

void fitWithGaus(TH2F *inputHist, std::string name, std::vector<double> &mean, std::vector<double> &meanErr, std::vector<double> &sigma, std::vector<double> &sigmaErr){
  
  std::vector<TH1F*> vHistos;
  TF1* gaus;
//  TCanvas *b = MakeCanvas("b","b",800,600);
  gStyle->SetOptStat(1111);
  //std::cout << "before loop" << std::endl;
  for(int bin = 1; bin < inputHist->GetNbinsX()+1; bin++){
     // add gaus fit
     //std::cout << "first" << std::endl;
    vHistos.push_back((TH1F*)inputHist->ProjectionY(name.c_str(),bin, bin));
    //std::cout << "projected" << std::endl;
    if(vHistos.back()->GetEntries() > 0)
    {
       // std::cout << "fit it" << std::endl;
        gaus =new TF1("gaus","gaus",-100,100);
        gaus->SetParameters(500.,1.,0.2);
        vHistos.back()->Fit("gaus","Q");
        gaus = vHistos.back()->GetFunction("gaus");
       
//        vHistos.back()->Draw();
//        b->SetLogy(kTRUE);
//        b->Print((outputdir+"/"+name+"bin_"+int2string(bin)+".png").c_str());
        
//         mean.push_back(-gaus->GetParameter(1));///((bins[bin-1]+bins[bin])*0.5));
//         sigma.push_back(gaus->GetParameter(2));
//         meanErr.push_back(gaus->GetParError(1));///((bins[bin]-bins[bin-1])*0.5));
//         sigmaErr.push_back(gaus->GetParError(2));

        mean.push_back(-vHistos.back()->GetMean());///((bins[bin-1]+bins[bin])*0.5));
        sigma.push_back(vHistos.back()->GetRMS());
        meanErr.push_back(vHistos.back()->GetMeanError());///((bins[bin]-bins[bin-1])*0.5));
        sigmaErr.push_back(vHistos.back()->GetRMSError());

    }
  }
//  delete b;
  return;
}


//=== MAIN MACRO ================================================================================================= 

void fitZm(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)
	       const Double_t nsigma=0,     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
           const TString input_section = "1"
) {
  gBenchmark->Start("fitZm");

  bool doFootprint=false;
  // bool doRecoilplot=true;
  bool doRecoilplot=false;
  // bool doDiago = false;
  bool doDiago = false;
  bool doKeys=false;
  bool doEta=false;
  bool doInclusive=true;
  bool doToys = false;
  bool doShittyRecoil=true;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  outputdir = outputDir;
  const Double_t mu_MASS  = 0.1057;

  // MET histogram binning and range
  const Int_t    NBINS   = 75;
  const Double_t METMAX  = 150;
  
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  

  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;

  enum { eMuMu2HLT=1, eMuMu1HLT1L1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk }; // event category enum

//   TString pufname = "../Tools/pileup_weights_2015B.root";

  // file format for output plots
  const TString format("png"); 
  rochcor2015 *rmcor = new rochcor2015();

  const TString directory1("../Recoil/");
  const TString directory("/afs/cern.ch/user/d/dalfonso/public/WZ/JULY5");
  // for Puppi, inclusive
  int rec_sig = 1;
  RecoilCorrector *recoilCorr = new  RecoilCorrector("","");
  RecoilCorrector *recoilCorrm = new  RecoilCorrector("","");
      RecoilCorrector *recoilCorr05 = new  RecoilCorrector("","");
      RecoilCorrector *recoilCorrm05 = new  RecoilCorrector("","");
      RecoilCorrector *recoilCorr051 = new  RecoilCorrector("","");
        RecoilCorrector *recoilCorr1 = new  RecoilCorrector("","");
   RecoilCorrector *recoilCorrShit = new  RecoilCorrector("","");
   char recName[150];
   // sprintf(recName,"../Recoil/PoorMansRecoilCorr_sec%s.root",input_section.Data());
   sprintf(recName,"../Recoil/PoorMansRecoilCorr_blah.root");
   recoilCorrShit->loadFileRatio(recName);
  if(doDiago && !doEta){
  std::cout << "loading diagonal recoil workspaces W+" << std::endl;
  recoilCorr->loadRooWorkspacesDiagMCtoCorrect(Form("%s/ZmmMCPuppi_v2/",directory1.Data()), rec_sig);
  recoilCorr->loadRooWorkspacesDiagData(Form("%s/ZmmDataPuppi_bkg_v2/",directory1.Data()), rec_sig);
  recoilCorr->loadRooWorkspacesDiagMC(Form("%s/ZmmMCPuppi_v2/",directory1.Data()), rec_sig);
  
  
  
  } else if (doDiago && doEta){
      recoilCorr05->loadRooWorkspacesDiagMCtoCorrect(Form("%s/ZmmMCPuppi_rap05/",directory1.Data()), rec_sig);
      // recoilCorr05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
      recoilCorr05->loadRooWorkspacesDiagData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap05/",directory1.Data()), rec_sig);
      recoilCorr05->loadRooWorkspacesDiagMC(Form("%s/ZmmMCPuppi_rap05/",directory1.Data()), rec_sig);
      
      // recoilCorrm05->loadRooWorkspacesDiagMCtoCorrect(Form("%s/ZmmMCPuppi_rap05/",directory1.Data()), rec_sig);
      // // recoilCorrm05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
      // recoilCorrm05->loadRooWorkspacesDiagData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap05/",directory1.Data()), rec_sig);
      // recoilCorrm05->loadRooWorkspacesDiagMC(Form("%s/ZmmMCPuppi_rap05/",directory1.Data()), rec_sig);

      recoilCorr051->loadRooWorkspacesDiagMCtoCorrect(Form("%s/ZmmMCPuppi_rap05-1/",directory1.Data()), rec_sig);
      // recoilCorr051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05-1/",directory.Data()));
      recoilCorr051->loadRooWorkspacesDiagData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap05-1/",directory1.Data()), rec_sig);
      recoilCorr051->loadRooWorkspacesDiagMC(Form("%s/ZmmMCPuppi_rap05-1/",directory1.Data()), rec_sig);
      
      recoilCorr1->loadRooWorkspacesDiagMCtoCorrect(Form("%s/ZmmMCPuppi_rap1/",directory1.Data()), rec_sig);
      // recoilCorr1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap1/",directory.Data()));
      recoilCorr1->loadRooWorkspacesDiagData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap1/",directory1.Data()), rec_sig);
      recoilCorr1->loadRooWorkspacesDiagMC(Form("%s/ZmmMCPuppi_rap1/",directory1.Data()), rec_sig);
  
  } else {
  // // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/JUNE25_ZmmMCPuppi/",directory1.Data()));
  // // recoilCorr->loadRooWorkspacesData(Form("%s/JUNE25_ZmmDataPuppi_bkg/",directory1.Data()));
  // // recoilCorr->loadRooWorkspacesMC(Form("%s/JUNE25_ZmmMCPuppi/",directory1.Data()));
  // // // // // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/APR4_ZmmMCPuppi/",directory1.Data()));

  // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMCPuppi_puFix/",directory1.Data()));
  recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMCPuppi_genBin/",directory1.Data()));
     recoilCorr->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg/",directory1.Data()));
    recoilCorr->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_genBin/",directory1.Data()));
  
  // PF MET
    // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMCPF_lowPU13/",directory1.Data()));
  // recoilCorr->loadRooWorkspacesData(Form("%s/ZmmDataPF_bkgTopEWK_lowPU13/",directory1.Data()));
  // recoilCorr->loadRooWorkspacesMC(Form("%s/ZmmMCPF_lowPU13/",directory1.Data()));
  
  // // make the diagonalized and central consistent
  // // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMCPuppi_puUp/",directory1.Data()));
  // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMCPuppi_puFix/",directory1.Data()));
  // recoilCorr->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_v3/",directory1.Data()));
  // // recoilCorr->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_puUp/",directory1.Data()));
  // recoilCorr->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_puFix/",directory1.Data()));
  
    // // // // // recoilCorr->loadRooWorkspacesData(Form("%s/ZeeDataPuppi_bkgTopEWK/",directory.Data()));
  // // // the ones i'm actually using
    // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPuppi_genFix/",directory.Data()));
   // recoilCorr->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK/",directory.Data()));
    // recoilCorr->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi/",directory.Data()));
  // // // // // // recoilCorr->loadRooWorkspacesMC(Form("%s/ZeeMCPuppi/",directory.Data()));
  // just test these
      // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/ZeeMCPuppi/",directory.Data()));
  // recoilCorr->loadRooWorkspacesData(Form("%s/ZeeDataPuppi_bkgTopEWK//",directory.Data()));
    // recoilCorr->loadRooWorkspacesMC(Form("%s/ZeeMCPuppi/",directory.Data()));
  // // // // recoilCorr->loadRooWorkspacesMC(Form("%s/ZeeMCPuppi/",directory.Data()));
  
  }
  if(doEta && !doDiago){
  // RecoilCorrector *recoilCorr05 = new  RecoilCorrector("","");
  recoilCorr05->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMCPuppi_rap05/",directory1.Data()));
  // recoilCorr05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
  recoilCorr05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap05/",directory1.Data()));
  recoilCorr05->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05/",directory1.Data()));
  
  // RecoilCorrector *recoilCorrm05 = new  RecoilCorrector("","");
  recoilCorrm05->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMCPuppi_rap05/",directory1.Data()));
  // recoilCorrm05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
  recoilCorrm05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap05/",directory1.Data()));
  recoilCorrm05->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05/",directory1.Data()));

  // RecoilCorrector *recoilCorr051 = new  RecoilCorrector("","");
  recoilCorr051->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMCPuppi_rap05-1/",directory1.Data()));
  // recoilCorr051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05-1/",directory.Data()));
  recoilCorr051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap05-1/",directory1.Data()));
  recoilCorr051->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05-1/",directory1.Data()));
  
    // RecoilCorrector *recoilCorr1 = new  RecoilCorrector("","");
  recoilCorr1->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMCPuppi_rap1/",directory1.Data()));
  // recoilCorr1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap1/",directory.Data()));
  recoilCorr1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap1/",directory1.Data()));
  recoilCorr1->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap1/",directory1.Data()));
  }
  
  // //   // ---------------------- KEYS - inclusive-------------
std::cout << "load inclusive RooKeys files" << std::endl;
  RecoilCorrector *recoilCorrKeys = new  RecoilCorrector("","");
  // recoilCorrKeys->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmpMCPuppi_keys/",directory.Data()));
  recoilCorrKeys->loadRooWorkspacesMCtoCorrectKeys(Form("%s/ZmmMCPuppi_keys/",directory.Data()));
  recoilCorrKeys->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK/",directory.Data()));
  recoilCorrKeys->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi/",directory.Data()));
  
  RecoilCorrector *recoilCorrKeysm = new  RecoilCorrector("","");
  // recoilCorrKeysm->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmmMCPuppi_keys/",directory.Data()));
  recoilCorrKeysm->loadRooWorkspacesMCtoCorrectKeys(Form("%s/ZmmMCPuppi_keys/",directory.Data()));
  recoilCorrKeysm->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK/",directory.Data()));
  recoilCorrKeysm->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi/",directory.Data()));
  
  // Load the Z data and Z MC Pt spectra
  // TFile *_rat1 = new TFile("zptratio.root");
  // TH1D *hh_mc;// = new TH1D("hh_diff","hh_diff",75,0,150);
  // hh_mc = (TH1D*)_rat1->Get("hZPtTruthNominal"); 
  // hh_mc->Scale(1/hh_mc->Integral()); // normalize
  
  TFile *_rat2 = new TFile("zPtDataMCresid_noEWK.root");
  TH1D *hh_diff;// = new TH1D("hh_diff","hh_diff",75,0,150);
  hh_diff = (TH1D*)_rat2->Get("hZptRatio");
  // hh_diff->Scale(1/hh_diff->Integral()); // normalize
  // hh_diff->Divide(hh_mc);

      // // efficiency files
      // const TString baseDir = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/10Sections/"; 
  // const TString dataHLTEffName_pos = baseDir + "MuHLTEff/MGpositive_part"+input_section+"/eff.root";
  // const TString dataHLTEffName_neg = baseDir + "MuHLTEff/MGnegative_part"+input_section+"/eff.root";
  // const TString zmmHLTEffName_pos  = baseDir + "MuHLTEff/CTpositive_part"+input_section+"/eff.root";
  // const TString zmmHLTEffName_neg  = baseDir + "MuHLTEff/CTnegative_part"+input_section+"/eff.root";

  // const TString dataSelEffName_pos = baseDir + "MuSITEff/MGpositive_part"+input_section+"/eff.root";
  // const TString dataSelEffName_neg = baseDir + "MuSITEff/MGnegative_part"+input_section+"/eff.root";
  // const TString zmmSelEffName_pos  = baseDir + "MuSITEff/CTpositive_part"+input_section+"/eff.root";
  // const TString zmmSelEffName_neg  = baseDir + "MuSITEff/CTnegative_part"+input_section+"/eff.root";

  // const TString dataTrkEffName_pos = baseDir + "MuSITEff/MGpositive_part"+input_section+"/eff.root";
  // const TString dataTrkEffName_neg = baseDir + "MuSITEff/MGnegative_part"+input_section+"/eff.root";
  // const TString zmmTrkEffName_pos  = baseDir + "MuSITEff/CTpositive_part"+input_section+"/eff.root";
  // const TString zmmTrkEffName_neg  = baseDir + "MuSITEff/CTnegative_part"+input_section+"/eff.root";

  // const TString dataStaEffName_pos = baseDir + "MuStaEff/MGpositive_part"+input_section+"/eff.root";
  // const TString dataStaEffName_neg = baseDir + "MuStaEff/MGnegative_part"+input_section+"/eff.root";
  // const TString zmmStaEffName_pos  = baseDir + "MuStaEff/CTpositive_part"+input_section+"/eff.root";
  // const TString zmmStaEffName_neg  = baseDir + "MuStaEff/CTnegative_part"+input_section+"/eff.root";

  // // efficiency files 2Bins

  // const TString dataHLTEff2BinName_pos = baseDir + "MuHLTEff/MGpositive_part"+input_section+"/eff.root";
  // const TString dataHLTEff2BinName_neg = baseDir + "MuHLTEff/MGnegative_part"+input_section+"/eff.root";
  // const TString zmmHLTEff2BinName_pos  = baseDir + "MuHLTEff/CTpositive_part"+input_section+"/eff.root";
  // const TString zmmHLTEff2BinName_neg  = baseDir + "MuHLTEff/CTnegative_part"+input_section+"/eff.root";

  // const TString dataSelEff2BinName_pos = baseDir + "MuSITEff/MGpositive_part"+input_section+"/eff.root";
  // const TString dataSelEff2BinName_neg = baseDir + "MuSITEff/MGnegative_part"+input_section+"/eff.root";
  // const TString zmmSelEff2BinName_pos  = baseDir + "MuSITEff/CTpositive_part"+input_section+"/eff.root";
  // const TString zmmSelEff2BinName_neg  = baseDir + "MuSITEff/CTnegative_part"+input_section+"/eff.root";

  // const TString dataTrkEff2BinName_pos = baseDir + "MuSITEff/MGpositive_part"+input_section+"/eff.root";
  // const TString dataTrkEff2BinName_neg = baseDir + "MuSITEff/MGnegative_part"+input_section+"/eff.root";
  // const TString zmmTrkEff2BinName_pos  = baseDir + "MuSITEff/CTpositive_part"+input_section+"/eff.root";
  // const TString zmmTrkEff2BinName_neg  = baseDir + "MuSITEff/CTnegative_part"+input_section+"/eff.root";

  // const TString dataStaEff2BinName_pos = baseDir + "MuStaEff/MGpositive_part"+input_section+"/eff.root";
  // const TString dataStaEff2BinName_neg = baseDir + "MuStaEff/MGnegative_part"+input_section+"/eff.root";
  // const TString zmmStaEff2BinName_pos  = baseDir + "MuStaEff/CTpositive_part"+input_section+"/eff.root";
  // const TString zmmStaEff2BinName_neg  = baseDir + "MuStaEff/CTnegative_part"+input_section+"/eff.root";

  // TString StaEffSignalShapeSys     = baseDir + "Results/MuStaSigSys.root";
  // TString StaEffBackgroundShapeSys = baseDir + "Results/MuStaBkgSys.root";
  // TString SelEffSignalShapeSys     = baseDir + "Results/MuSITSigSys.root";
  // TString SelEffBackgroundShapeSys = baseDir + "Results/MuSITBkgSys.root";
  
  const TString baseDir = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/";
  const TString dataHLTEffName_pos = baseDir + "MuHLTEff/MGpositive/eff.root";
  const TString dataHLTEffName_neg = baseDir + "MuHLTEff/MGnegative/eff.root";
  const TString zmmHLTEffName_pos  = baseDir + "MuHLTEff/CTpositive/eff.root";
  const TString zmmHLTEffName_neg  = baseDir + "MuHLTEff/CTnegative/eff.root";

  const TString dataSelEffName_pos = baseDir + "MuSITEff/MGpositive_FineBin/eff.root";
  const TString dataSelEffName_neg = baseDir + "MuSITEff/MGnegative_FineBin/eff.root";
  const TString zmmSelEffName_pos  = baseDir + "MuSITEff/CTpositive/eff.root";
  const TString zmmSelEffName_neg  = baseDir + "MuSITEff/CTnegative/eff.root";

  const TString dataTrkEffName_pos = baseDir + "MuSITEff/MGpositive_FineBin/eff.root";
  const TString dataTrkEffName_neg = baseDir + "MuSITEff/MGnegative_FineBin/eff.root";
  const TString zmmTrkEffName_pos  = baseDir + "MuSITEff/CTpositive/eff.root";
  const TString zmmTrkEffName_neg  = baseDir + "MuSITEff/CTnegative/eff.root";

  const TString dataStaEffName_pos = baseDir + "MuStaEff/MGpositive/eff.root";
  const TString dataStaEffName_neg = baseDir + "MuStaEff/MGnegative/eff.root";
  const TString zmmStaEffName_pos  = baseDir + "MuStaEff/CTpositive/eff.root";
  const TString zmmStaEffName_neg  = baseDir + "MuStaEff/CTnegative/eff.root";

  // efficiency files 2Bins

  const TString dataHLTEff2BinName_pos = baseDir + "MuHLTEff/MGpositive/eff.root";
  const TString dataHLTEff2BinName_neg = baseDir + "MuHLTEff/MGnegative/eff.root";
  const TString zmmHLTEff2BinName_pos  = baseDir + "MuHLTEff/CTpositive/eff.root";
  const TString zmmHLTEff2BinName_neg  = baseDir + "MuHLTEff/CTnegative/eff.root";

  const TString dataSelEff2BinName_pos = baseDir + "MuSITEff/MGpositive_FineBin/eff.root";
  const TString dataSelEff2BinName_neg = baseDir + "MuSITEff/MGnegative_FineBin/eff.root";
  const TString zmmSelEff2BinName_pos  = baseDir + "MuSITEff/CTpositive/eff.root";
  const TString zmmSelEff2BinName_neg  = baseDir + "MuSITEff/CTnegative/eff.root";

  const TString dataTrkEff2BinName_pos = baseDir + "MuSITEff/MGpositive_FineBin/eff.root";
  const TString dataTrkEff2BinName_neg = baseDir + "MuSITEff/MGnegative_FineBin/eff.root";
  const TString zmmTrkEff2BinName_pos  = baseDir + "MuSITEff/CTpositive/eff.root";
  const TString zmmTrkEff2BinName_neg  = baseDir + "MuSITEff/CTnegative/eff.root";

  const TString dataStaEff2BinName_pos = baseDir + "MuStaEff/MGpositive/eff.root";
  const TString dataStaEff2BinName_neg = baseDir + "MuStaEff/MGnegative/eff.root";
  const TString zmmStaEff2BinName_pos  = baseDir + "MuStaEff/CTpositive/eff.root";
  const TString zmmStaEff2BinName_neg  = baseDir + "MuStaEff/CTnegative/eff.root";

  TString StaEffSignalShapeSys     = baseDir + "Results/MuStaSigSys.root";
  TString StaEffBackgroundShapeSys = baseDir + "Results/MuStaBkgSys.root";
  TString SelEffSignalShapeSys     = baseDir + "Results/MuSITSigSys.root";
  TString SelEffBackgroundShapeSys = baseDir + "Results/MuSITBkgSys.root";

  //
  // input ntuple file names
  //
  enum { eData, eZmumu, eEWK, eBKG, eQCD, eAntiData, eAntiWmunu, eAntiEWK };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;

    // fnamev.push_back("/data/t3home000/sabrandt/2018_11_11_LowPUFlat_Aram/Zmumu/ntuples/data_select.root"); typev.push_back(eData);
  // fnamev.push_back("/data/t3home000/sabrandt/2018_11_11_LowPUFlat_Aram/Zmumu/ntuples/zmm_select.raw.root"); typev.push_back(eZmumu);
  
  // // fnamev.push_back("/data/t3home000/sabrandt/2018_09_07_Masters_10sec/Lumi10Parts_Part"+input_section+"/Zmumu/ntuples/data_select.root"); typev.push_back(eData);
  fnamev.push_back("/data/t3home000/sabrandt/2018_09_07_Masters_Incl/Zmumu/ntuples/data_select.root"); typev.push_back(eData);
  fnamev.push_back("/data/t3home000/sabrandt/2018_09_07_Masters_Incl_2/Zmumu/ntuples/zmm_select.raw.root"); typev.push_back(eZmumu);
  fnamev.push_back("/data/t3home000/sabrandt/2018_09_07_Masters_Incl/Zmumu/ntuples/zxx_select.raw.root"); typev.push_back(eBKG);
  fnamev.push_back("/data/t3home000/sabrandt/2018_09_07_Masters_Incl/Zmumu/ntuples/wx_select.raw.root"); typev.push_back(eBKG);
  fnamev.push_back("/data/t3home000/sabrandt/2018_09_07_Masters_Incl/Zmumu/ntuples/zz_select.raw.root"); typev.push_back(eEWK);
  fnamev.push_back("/data/t3home000/sabrandt/2018_09_07_Masters_Incl/Zmumu/ntuples/wz_select.raw.root"); typev.push_back(eEWK);
  fnamev.push_back("/data/t3home000/sabrandt/2018_09_07_Masters_Incl/Zmumu/ntuples/ww_select.raw.root"); typev.push_back(eEWK);
  fnamev.push_back("/data/t3home000/sabrandt/2018_09_07_Masters_Incl/Zmumu/ntuples/top_select.raw.root"); typev.push_back(eEWK);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  //
  // Declare MET histograms
  //
  TH1D *hDataMet   = new TH1D("hDataMet","",  NBINS,0,METMAX); hDataMet->Sumw2();
  TH1D *hDataMetm  = new TH1D("hDataMetm","", NBINS,0,METMAX); hDataMetm->Sumw2();  
  TH1D *hDataMetp  = new TH1D("hDataMetp","", NBINS,0,METMAX); hDataMetp->Sumw2();
  TH1D *hWmunuMet  = new TH1D("hWmunuMet","", NBINS,0,METMAX); hWmunuMet->Sumw2();
  TH1D *hWmunuMetp = new TH1D("hWmunuMetp","",NBINS,0,METMAX); hWmunuMetp->Sumw2();
  TH1D *hWmunuMetm = new TH1D("hWmunuMetm","",NBINS,0,METMAX); hWmunuMetm->Sumw2();
  TH1D *hEWKMet    = new TH1D("hEWKMet", "",  NBINS,0,METMAX); hEWKMet->Sumw2();
  TH1D *hEWKMetp   = new TH1D("hEWKMetp", "", NBINS,0,METMAX); hEWKMetp->Sumw2();
  TH1D *hEWKMetm   = new TH1D("hEWKMetm", "", NBINS,0,METMAX); hEWKMetm->Sumw2();
  
  TH1D *hDataRecoilp  = new TH1D("hDataRecoilp","", NBINS,0,METMAX); hDataRecoilp->Sumw2();
  TH1D *hWmunuRecoilp = new TH1D("hWmunuRecoilp","",NBINS,0,METMAX); hWmunuRecoilp->Sumw2();
  TH1D *hEWKRecoilp   = new TH1D("hEWKRecoilp", "", NBINS,0,METMAX); hEWKRecoilp->Sumw2();
  
  TH1D *hDataZpt  = new TH1D("hDataZpt","", NBINS,0,METMAX); hDataZpt->Sumw2();
  TH1D *hWmunuZpt = new TH1D("hWmunuZpt","",NBINS,0,METMAX); hWmunuZpt->Sumw2();
  TH1D *hEWKZpt   = new TH1D("hEWKZpt", "", NBINS,0,METMAX); hEWKZpt->Sumw2();
  
  TH1D *hWmunuMet_RecoilUp  = new TH1D("hWmunuMet_RecoilUp", "",NBINS,0,METMAX); hWmunuMet_RecoilUp->Sumw2();
  TH1D *hWmunuMetp_RecoilUp = new TH1D("hWmunuMetp_RecoilUp","",NBINS,0,METMAX); hWmunuMetp_RecoilUp->Sumw2();
  TH1D *hWmunuMetm_RecoilUp = new TH1D("hWmunuMetm_RecoilUp","",NBINS,0,METMAX); hWmunuMetm_RecoilUp->Sumw2();
  TH1D *hWmunuMet_RecoilDown  = new TH1D("hWmunuMet_RecoilDown", "",NBINS,0,METMAX); hWmunuMet_RecoilDown->Sumw2();
  TH1D *hWmunuMetp_RecoilDown = new TH1D("hWmunuMetp_RecoilDown","",NBINS,0,METMAX); hWmunuMetp_RecoilDown->Sumw2();
  TH1D *hWmunuMetm_RecoilDown = new TH1D("hWmunuMetm_RecoilDown","",NBINS,0,METMAX); hWmunuMetm_RecoilDown->Sumw2();

  TH1D *hWmunuMet_RecoilCUp  = new TH1D("hWmunuMet_RecoilCUp", "",NBINS,0,METMAX); hWmunuMet_RecoilCUp->Sumw2();
  TH1D *hWmunuMetp_RecoilCUp = new TH1D("hWmunuMetp_RecoilCUp","",NBINS,0,METMAX); hWmunuMetp_RecoilCUp->Sumw2();
  TH1D *hWmunuMetm_RecoilCUp = new TH1D("hWmunuMetm_RecoilCUp","",NBINS,0,METMAX); hWmunuMetm_RecoilCUp->Sumw2();
  TH1D *hWmunuMet_RecoilCDown  = new TH1D("hWmunuMet_RecoilCDown", "",NBINS,0,METMAX); hWmunuMet_RecoilCDown->Sumw2();
  TH1D *hWmunuMetp_RecoilCDown = new TH1D("hWmunuMetp_RecoilCDown","",NBINS,0,METMAX); hWmunuMetp_RecoilCDown->Sumw2();
  TH1D *hWmunuMetm_RecoilCDown = new TH1D("hWmunuMetm_RecoilCDown","",NBINS,0,METMAX); hWmunuMetm_RecoilCDown->Sumw2();

  TH1D *hWmunuMet_ScaleUp  = new TH1D("hWmunuMet_ScaleUp", "",NBINS,0,METMAX); hWmunuMet_ScaleUp->Sumw2();
  TH1D *hWmunuMetp_ScaleUp = new TH1D("hWmunuMetp_ScaleUp","",NBINS,0,METMAX); hWmunuMetp_ScaleUp->Sumw2();
  TH1D *hWmunuMetm_ScaleUp = new TH1D("hWmunuMetm_ScaleUp","",NBINS,0,METMAX); hWmunuMetm_ScaleUp->Sumw2();
  TH1D *hWmunuMet_ScaleDown  = new TH1D("hWmunuMet_ScaleDown", "",NBINS,0,METMAX); hWmunuMet_ScaleDown->Sumw2();
  TH1D *hWmunuMetp_ScaleDown = new TH1D("hWmunuMetp_ScaleDown","",NBINS,0,METMAX); hWmunuMetp_ScaleDown->Sumw2();
  TH1D *hWmunuMetm_ScaleDown = new TH1D("hWmunuMetm_ScaleDown","",NBINS,0,METMAX); hWmunuMetm_ScaleDown->Sumw2();

  TH1D *hAntiDataMet   = new TH1D("hAntiDataMet","",  NBINS,0,METMAX); hAntiDataMet->Sumw2();
  TH1D *hAntiDataMetm  = new TH1D("hAntiDataMetm","", NBINS,0,METMAX); hAntiDataMetm->Sumw2();  
  TH1D *hAntiDataMetp  = new TH1D("hAntiDataMetp","", NBINS,0,METMAX); hAntiDataMetp->Sumw2();
  TH1D *hAntiWmunuMet  = new TH1D("hAntiWmunuMet","", NBINS,0,METMAX); hAntiWmunuMet->Sumw2();
  TH1D *hAntiWmunuMetp = new TH1D("hAntiWmunuMetp","",NBINS,0,METMAX); hAntiWmunuMetp->Sumw2();
  TH1D *hAntiWmunuMetm = new TH1D("hAntiWmunuMetm","",NBINS,0,METMAX); hAntiWmunuMetm->Sumw2();
  TH1D *hAntiEWKMet    = new TH1D("hAntiEWKMet", "",  NBINS,0,METMAX); hAntiEWKMet->Sumw2();
  TH1D *hAntiEWKMetp   = new TH1D("hAntiEWKMetp", "", NBINS,0,METMAX); hAntiEWKMetp->Sumw2();
  TH1D *hAntiEWKMetm   = new TH1D("hAntiEWKMetm", "", NBINS,0,METMAX); hAntiEWKMetm->Sumw2();

  Double_t vEtaBinsFP[] = {2.5,2.4,2.3,2.2,2.1,2.0,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1.0,-1.1,-1.2,-1.3,-1.4,-1.5,-1.6,-1.7,-1.8,-1.9,-2.0,-2.1,-2.2,-2.3,-2.4,-2.5};
  Int_t nBinsFP = sizeof(vEtaBinsFP)/sizeof(vEtaBinsFP[0])-1;
  vector<Float_t> upper(nBinsFP);
  vector<Float_t> lower(nBinsFP);
  std::cout << "nEtaBins = " << nBinsFP << std::endl;
  
  // vector<Float_t> upper;
  for(int i = 0; i < nBinsFP; ++i){
	  std::cout << "bin val = " << vEtaBinsFP[i] << std::endl;
	  upper[i] = vEtaBinsFP[i];
	  lower[i] = vEtaBinsFP[i+1];
	  // lower.push_back(vEtaBinsFP[i+1]);
  }
   
  Double_t vListOfFPs[] = {0, -0.00109938, -0.00267664, -0.00935168, -0.0116393, -0.00877814, -0.00073066, -0.00168777, -0.0114005, -0.0102761, -0.010612, -0.00345421, -0.0092375, -0.00850389, -0.00988786, -0.00506354, -0.0064118, -0.00408221, -0.00398614, -0.0172608, -0.0104789, -0.007598, -0.0177603, -0.0129501, -0.00846278, -0.010801, -0.0179897, -0.014265, -0.004627, -0.00198961, -0.00940792, -0.0137091, -0.00842374, -0.00451072, -0.0134227, -0.013647, -0.000902762, -0.00728768, -0.00516646, -0.0112141, -0.000557311, -0.0124703, -0.0169047, -0.00473276, -0.00623816, -0.00367185, -0.0133091, -0.00340216, -0.0278066, 0};
  vector<Float_t> fpcorr;
  // std::cout << "size of list " << 
  for(int i = 0; i < nBinsFP; ++i){
	  std::cout << "loop bin " << i << std::endl;
	  std::cout << vListOfFPs[i] << std::endl;
	  fpcorr.push_back(vListOfFPs[i]);
  }
  
  double bins[] = {0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,275,300};
  // double bins[] = {0,2,5,10,20,30,40,50,60,70,80,90,100};
// double bins[] ={0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,80,85,90,95,100/*,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,275,300*/};
//     double bins[] = {0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,80,85,90,95,100};// ={0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100};// {0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100/*,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200*/};
  int nbins = sizeof(bins)/sizeof(double)-1; // calculate the number of bins
  
  TH2F *hU1vsZpt_rsp = new TH2F("hU1vsZpt_rsp","",nbins, bins,400,-12.0,12.0); // Get resp. & res. for Parallel component from this
  TH2F *hU1vsZpt = new TH2F("hU1vsZpt","",nbins, bins,400, -500,300);
  TH2F *hU2vsZpt = new TH2F("hU2vsZpt","",nbins, bins,300,-500,300); // Get resp. & res. for perpendicular component from this
  
  TH2F *hU1vsZpt_rsp_MC = new TH2F("hU1vsZpt_rsp_MC","",nbins, bins,400,-12.0,12.0); // Get resp. & res. for Parallel component from this
  TH2F *hU1vsZpt_MC = new TH2F("hU1vsZpt_MC","",nbins, bins,400, -500,300);
  TH2F *hU2vsZpt_MC = new TH2F("hU2vsZpt_MC","",nbins, bins,300,-500,300); // Get resp. & res. for perpendicular component from this
  
  TH2F *hU1vsZpt_rsp_MC_raw = new TH2F("hU1vsZpt_rsp_MC_raw","",nbins, bins,400,-12.0,12.0); // Get resp. & res. for Parallel component from this
  TH2F *hU1vsZpt_MC_raw = new TH2F("hU1vsZpt_MC_raw","",nbins, bins,400, -500,300);
  TH2F *hU2vsZpt_MC_raw = new TH2F("hU2vsZpt_MC_raw","",nbins, bins,300,-500,300); // Get resp. & res. for perpendicular component from this
  
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
// 
  
  
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi;
  Float_t scale1fb,scale1fbUp,scale1fbDown;
  Float_t met, metPhi, sumEt,  u1, u2;
  Int_t   q1, q2;
  UInt_t  category;
  //TLorentzVector *lep=0;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0;
    
  TFile *infile=0;
  TTree *intree=0;

  //
  // Loop over files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);	  assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);

    intree->SetBranchAddress("runNum",      &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",     &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",      &evtNum);    // event number
    intree->SetBranchAddress("npv",         &npv);       // number of primary vertices
    intree->SetBranchAddress("npu",         &npu);       // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",      &genVPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",     &genVPhi);   // GEN W boson phi (signal MC)   
    intree->SetBranchAddress("scale1fb",    &scale1fb);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp",    &scale1fbUp);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown",    &scale1fbDown);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("puppiMet",    &met);       // MET
    intree->SetBranchAddress("puppiMetPhi", &metPhi);    // phi(MET)
    // intree->SetBranchAddress("met",    &met);       // MET
    // intree->SetBranchAddress("metPhi", &metPhi);    // phi(MET)
    intree->SetBranchAddress("sumEt",       &sumEt);     // Sum ET
    intree->SetBranchAddress("puppiU1",     &u1);        // parallel component of recoil
    intree->SetBranchAddress("puppiU2",     &u2);        // perpendicular component of recoil
    // intree->SetBranchAddress("u1",     &u1);        // parallel component of recoil
    // intree->SetBranchAddress("u2",     &u2);        // perpendicular component of recoil
    intree->SetBranchAddress("q1",          &q1);	  // charge of tag lepton
    intree->SetBranchAddress("q2",          &q2);	  // charge of probe lepton
    //intree->SetBranchAddress("lep",       &lep);       // lepton 4-vector
    intree->SetBranchAddress("lep1",        &lep1);        // tag lepton 4-vector
    intree->SetBranchAddress("lep2",        &lep2);        // probe lepton 4-vector
    intree->SetBranchAddress("dilep",       &dilep);       // dilepton 4-vector
    intree->SetBranchAddress("category",    &category);    // dilepton category
    
    
    Double_t mt=-999;

    //
    // loop over events
    //
    // for(UInt_t ientry=0; ientry<500000; ientry++) {
    // for(UInt_t ientry=0; ientry<(int)(intree->GetEntries()*0.1); ientry++) {
    // for(UInt_t ientry=(int)(intree->GetEntries()*0.1); ientry<(int)(intree->GetEntries()*0.2); ientry++) {
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);

//       if(!((category==eMuMu2HLT) || (category==eMuMu1HLT) || (category==eMuMu1HLT1L1))) continue;
      if(ientry%100000==0) std::cout << "On Entry.... " << ientry << std::endl;
      if(!((category==1) || (category==2) || (category==3))) continue;

      double pU1         = 0;  //--
      double pU2         = 0;  //--
      
      Double_t effdata, effmc;
      Double_t corr=1;

      if(fabs(lep1->Eta()) > ETA_CUT) continue;
      if(fabs(lep2->Eta()) > ETA_CUT) continue;
	  
	  TVector2 vLepRaw1((lep1->Pt())*cos(lep1->Phi()),(lep1->Pt())*sin(lep1->Phi()));
      TVector2 vLepRaw2((lep2->Pt())*cos(lep2->Phi()),(lep2->Pt())*sin(lep2->Phi()));
      TVector2 vDilepRaw = vLepRaw1 + vLepRaw2;
	  // if(vDilepRaw.Mod() < 25 || vDilepRaw.Mod() > 37.5) continue;
    
// //       if(nJets30 <2) continue;
      // if(dilep->Pt() < 25) continue;
//   
      double tl1Pt = 0;
      double tl1Phi = 0;
      
      effdata=1; effmc=1;    
          if(q1>0) { 
            effdata *= (1.-dataHLTEff_pos.getEff(lep1->Eta(), lep1->Pt())); 
            effmc   *= (1.-zmmHLTEff_pos.getEff(lep1->Eta(), lep1->Pt())); 
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(lep1->Eta(), lep1->Pt())); 
            effmc   *= (1.-zmmHLTEff_neg.getEff(lep1->Eta(), lep1->Pt())); 
          }
          if(q2>0) {
            effdata *= (1.-dataHLTEff_pos.getEff(lep2->Eta(), lep2->Pt())); 
            effmc   *= (1.-zmmHLTEff_pos.getEff(lep2->Eta(), lep2->Pt()));
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(lep2->Eta(), lep2->Pt())); 
            effmc   *= (1.-zmmHLTEff_neg.getEff(lep2->Eta(), lep2->Pt()));
          }
          effdata = 1.-effdata;
          effmc   = 1.-effmc;
          corr *= effdata/effmc;
    
          effdata=1; effmc=1;
          if(q1>0) { 
            effdata *= dataSelEff_pos.getEff(lep1->Eta(), lep1->Pt()); 
            effmc   *= zmmSelEff_pos.getEff(lep1->Eta(), lep1->Pt()); 
          } else {
            effdata *= dataSelEff_neg.getEff(lep1->Eta(), lep1->Pt()); 
            effmc   *= zmmSelEff_neg.getEff(lep1->Eta(), lep1->Pt()); 
          }
          if(q2>0) {
            effdata *= dataSelEff_pos.getEff(lep2->Eta(), lep2->Pt()); 
            effmc   *= zmmSelEff_pos.getEff(lep2->Eta(), lep2->Pt());
          } else {
            effdata *= dataSelEff_neg.getEff(lep2->Eta(), lep2->Pt()); 
            effmc   *= zmmSelEff_neg.getEff(lep2->Eta(), lep2->Pt());
          }
          corr *= effdata/effmc;
    
          effdata=1; effmc=1;
          if(q1>0) { 
            effdata *= dataStaEff_pos.getEff(lep1->Eta(), lep1->Pt()); 
            effmc   *= zmmStaEff_pos.getEff(lep1->Eta(), lep1->Pt()); 
          } else {
            effdata *= dataStaEff_neg.getEff(lep1->Eta(), lep1->Pt()); 
            effmc   *= zmmStaEff_neg.getEff(lep1->Eta(), lep1->Pt()); 
          }
          if(q2>0) {
            effdata *= dataStaEff_pos.getEff(lep2->Eta(), lep2->Pt()); 
            effmc   *= zmmStaEff_pos.getEff(lep2->Eta(), lep2->Pt());
          } else {
            effdata *= dataStaEff_neg.getEff(lep2->Eta(), lep2->Pt()); 
            effmc   *= zmmStaEff_neg.getEff(lep2->Eta(), lep2->Pt());
          }
          corr *= effdata/effmc;
    
	  /*  effdata=1; effmc=1;
          if(mu1->q>0) { 
            effdata *= dataTrkEff_pos.getEff(mu1->eta, mu1->pt); 
            effmc   *= zmmTrkEff_pos.getEff(mu1->eta, mu1->pt); 
          } else {
            effdata *= dataTrkEff_neg.getEff(mu1->eta, mu1->pt); 
            effmc   *= zmmTrkEff_neg.getEff(mu1->eta, mu1->pt); 
          }
          if(mu2->q>0) {
            effdata *= dataTrkEff_pos.getEff(mu2->eta, mu2->pt); 
            effmc   *= zmmTrkEff_pos.getEff(mu2->eta, mu2->pt);
          } else {
            effdata *= dataTrkEff_neg.getEff(mu2->eta, mu2->pt); 
            effmc   *= zmmTrkEff_neg.getEff(mu2->eta, mu2->pt);
          }
          //corr *= effdata/effmc;*/
		  // set corr to 1 for low PU
	  // corr = 1;
      
      if(typev[ifile]==eData) {
        
        TLorentzVector mu1;
        TLorentzVector mu2;
        mu1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),mu_MASS);
        mu2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),mu_MASS);
        float qter1=1.0;
        float qter2=1.0;
		// temporarily remove rochester corr for low PU
        rmcor->momcor_data(mu1,q1,0,qter1);
        rmcor->momcor_data(mu2,q2,0,qter2);
        if(mu1.Pt()        < PT_CUT)  continue;
        if(mu2.Pt()        < PT_CUT)  continue;
//         if(lep1->Pt()        < PT_CUT)  continue;
//         if(lep2->Pt()        < PT_CUT)  continue;
        //if(dilep->M()        < MASS_LOW)  continue;
        //if(dilep->M()        > MASS_HIGH) continue;
	if((mu1+mu2).M()        < MASS_LOW)  continue;
	if((mu1+mu2).M()        > MASS_HIGH) continue;
        
        
		 TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
         TVector2 vLepCor1((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
         TVector2 vLepCor2((mu2.Pt())*cos(mu2.Phi()),(mu2.Pt())*sin(mu2.Phi()));
         TVector2 vDilepCor = vLepCor1 + vLepCor2;
           // std::cout << "vDilepCor " << vDilepCor.Mod() <<  "  orig " << dilep->Pt() << "  raw " << vDilepRaw.Mod() << std::endl;
         Double_t corrMetWithLepton = (vMetCorr + vDilepRaw - vDilepCor).Mod();
         Double_t corrMetPhiLepton = (vMetCorr + vDilepRaw - vDilepCor).Phi();
        // For making W-like met spectrum
//       if(q1 > 0) {tl1Pt = mu1.Pt(); tl1Phi = mu1.Phi();}
//       if(q2 > 0) {tl1Pt = mu2.Pt(); tl1Phi = mu2.Phi();}
//       double lMX = -tl1Pt*cos(tl1Phi) - u1*cos(tl1Phi) + u2*sin(tl1Phi);
//       double lMY = -tl1Pt*sin(tl1Phi) - u1*sin(tl1Phi) - u2*cos(tl1Phi);
//       met= sqrt(lMX*lMX + lMY*lMY);
      
        hDataMet->Fill(met);
        if(q1*q2<0) {
	      double pUX  = corrMetWithLepton*cos(corrMetPhiLepton) + vDilepCor.Mod()*cos(vDilepCor.Phi());
          double pUY  = corrMetWithLepton*sin(corrMetPhiLepton) + vDilepCor.Mod()*sin(vDilepCor.Phi());
          double pU   = sqrt(pUX*pUX+pUY*pUY);
          double pCos = - (pUX*cos(vDilepCor.Phi()) + pUY*sin(vDilepCor.Phi()))/pU;
          double pSin =   (pUX*sin(vDilepCor.Phi()) - pUY*cos(vDilepCor.Phi()))/pU;
          double pU1  = pU*pCos; // U1 in sample to Correct (WMC or ZMC)
          double pU2  = pU*pSin; // U2 in sample to Correct (WMC or ZMC)
          // hDataRecoilp->Fill(abs(pU1));
          // hDataMetp->Fill(abs(pU2));
          hDataRecoilp->Fill(pU);
          hDataMetp->Fill(corrMetWithLepton);
          hU1vsZpt_rsp->Fill(vDilepCor.Mod(),pU1/vDilepCor.Mod());
          hU1vsZpt->Fill(vDilepCor.Mod(),pU1);
          hU2vsZpt->Fill(vDilepCor.Mod(),pU2);
          hDataZpt->Fill(vDilepCor.Mod());
        }
        else    { hDataMetm->Fill(corrMetPhiLepton); }
      } else if(typev[ifile]==eAntiData) {
        hAntiDataMet->Fill(met);
        if(q1*q2<0) { hAntiDataMetp->Fill(met); } 
        else    { hAntiDataMetm->Fill(met); }
      } else {
        Double_t weight = 1;
        weight *= scale1fb*lumi*corr;
        // weight *= scale1fbUp*lumi*corr;
        // weight *= scale1fbDown*lumi*corr;
        
        TLorentzVector mu1;TLorentzVector mu2;
        mu1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),mu_MASS);
        mu2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),mu_MASS);
        float qter1=1.0;
        float qter2=1.0;


		// remove roch cor for low PU
        rmcor->momcor_mc(mu1,q1,0,qter1);
        rmcor->momcor_mc(mu2,q2,0,qter2);
//         Double_t lepPt = mu1.Pt();
          Double_t lp1 = mu1.Pt();
          Double_t lp2 = mu2.Pt();
          TLorentzVector l1, l2, dl;
          l1.SetPtEtaPhiM(lp1,lep1->Eta(),lep1->Phi(),mu_MASS);
          l2.SetPtEtaPhiM(lp2,lep2->Eta(),lep2->Phi(),mu_MASS);
          dl=l1+l2;
          double mll=dl.M();
		  

		 Double_t eleCorr1=0.0;
		 Double_t eleCorr2=0.0;
         if(doFootprint){
            for (Int_t i=0; i<nBinsFP; i++) {
              if( l1.Eta()>lower[i] && l1.Eta()<upper[i] ) {
                eleCorr1=fpcorr[i];
              }
            }//footprint?
            for (Int_t i=0; i<nBinsFP; i++) {
              if( l2.Eta()>lower[i] && l2.Eta()<upper[i] ) {
                eleCorr2=fpcorr[i];
              }
            }
        }
		 TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
         TVector2 vLepCor1((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
         TVector2 vLepCor2((mu2.Pt())*cos(mu2.Phi()),(mu2.Pt())*sin(mu2.Phi()));
         // TVector2 vDilepCor = vLepCor1 + vLepCor2;
		 vLepCor1=vLepCor1*(1-eleCorr1);
         vLepCor2=vLepCor2*(1-eleCorr2);
           
         TVector2 vDilepCor = vLepCor1 + vLepCor2;
		 // if(vDilepCor.)
         Double_t corrMetWithLepton = (vMetCorr + vDilepRaw - vDilepCor).Mod();
         Double_t corrMetPhiLepton = (vMetCorr + vDilepRaw - vDilepCor).Phi();
		          
          if(l1.Pt()        < PT_CUT)  continue;
          if(l2.Pt()        < PT_CUT)  continue;
          if(mll        < MASS_LOW)  continue;
          if(mll        > MASS_HIGH) continue;

/*        // For making W-like MET spectrum
        if(q1 > 0) {tl1Pt = mu1.Pt(); tl1Phi = mu1.Phi();}
        else if(q2 > 0) {tl1Pt = mu2.Pt(); tl1Phi = mu2.Phi();}
        double lMX = -tl1Pt*cos(tl1Phi) - u1*cos(tl1Phi) + u2*sin(tl1Phi);
        double lMY = -tl1Pt*sin(tl1Phi) - u1*sin(tl1Phi) - u2*cos(tl1Phi);
        met= sqrt(lMX*lMX + lMY*lMY);*/
        
        if(typev[ifile]==eZmumu || typev[ifile]==eBKG) {
          Double_t corrMet=met, corrMetPhi=metPhi;
	  
         
        double u1r = u1;
        double u2r = u2;

	  //
          hWmunuMet->Fill(corrMet,weight);
          if(q1*q2<0) {
              pU1 = 0; pU2 = 0;
	          double bin = 0;
	          for(int i = 0; i <= hh_diff->GetNbinsX();++i){
        		if(vDilepCor.Mod() > hh_diff->GetBinLowEdge(i) && vDilepCor.Mod() < hh_diff->GetBinLowEdge(i+1)){ bin = i; break; }
	          }
	          double w2 = hh_diff->GetBinContent(bin);
             double recoilWeight = 1;
	      if(doKeys) {
                  recoilCorrKeys->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,vDilepCor.Mod(),vDilepCor.Phi(),vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
	      }  else if(doEta) {
                if(fabs(dilep->Eta())<0.5)
                  recoilCorr05->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,vDilepCor.Mod(),vDilepCor.Phi(),vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
                else if (fabs(dilep->Eta())>=0.5 && fabs(dilep->Eta())<1.0)
                  recoilCorr051->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,vDilepCor.Mod(),vDilepCor.Phi(),vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
                else
                  recoilCorr1->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,vDilepCor.Mod(),vDilepCor.Phi(),vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doDiago); 
          } else if(doToys){
            recoilCorr->CorrectFromToys(corrMetWithLepton,corrMetPhiLepton,vDilepCor.Mod(),vDilepCor.Phi(),vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0);
          
          }else if(doInclusive){

	      // recoilCorr->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,vDilepCor.Mod(),vDilepCor.Phi(),vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
	      recoilCorr->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,genVPt,genVPhi,vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
	      }

		  
	      //recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,tl1Pt,tl1Pt,pU1,pU2,0,0,0);
//               recoilCorr->CorrectFromToys(corrMet,corrMetPhi,genVPt,genVPhi,dl.Pt(),dl.Phi(),pU1,pU2,0,0,0);
           double pUX  = corrMetWithLepton*cos(corrMetPhiLepton) + vDilepCor.Mod()*cos(vDilepCor.Phi());
           double pUY  = corrMetWithLepton*sin(corrMetPhiLepton) + vDilepCor.Mod()*sin(vDilepCor.Phi());
           double pU   = sqrt(pUX*pUX+pUY*pUY);
           double pCos = - (pUX*cos(vDilepCor.Phi()) + pUY*sin(vDilepCor.Phi()))/pU;
           double pSin =   (pUX*sin(vDilepCor.Phi()) - pUY*cos(vDilepCor.Phi()))/pU;
           pU1   = pU*pCos; // U1 in data
           pU2   = pU*pSin; // U2 in data
           // std::cout << "u1 vs u1: " << pU1 << " vs " << u1r << std::endl;

          if(doShittyRecoil){
               recoilCorrShit->CorrectShitty(corrMetWithLepton,corrMetPhiLepton,genVPt,genVPhi,vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,recoilWeight);
          }

              // hU1vsZpt_rsp_MC_raw->Fill(dilep->Pt(),u1/dilep->Pt(),weight);
              // hU1vsZpt_MC_raw->Fill(dilep->Pt(),u1,weight);
              // hU2vsZpt_MC_raw->Fill(dilep->Pt(),u2,weight);

              // // hU1vsZpt_rsp_MC->Fill(dilep->Pt(),pU1/dilep->Pt(),weight);
              // // hU1vsZpt_MC->Fill(dilep->Pt(),pU1,weight);
              // // hU2vsZpt_MC->Fill(dilep->Pt(),pU2,weight);
			 // hU1vsZpt_rsp_MC->Fill(vDilepCor.Mod(),pU1/vDilepCor.Mod(),weight);
            // hU1vsZpt_MC->Fill(vDilepCor.Mod(),pU1,weight);
            // hU2vsZpt_MC->Fill(vDilepCor.Mod(),pU2,weight);
			  
			  // TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
            // calculate the corrected MET
            // TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
            // Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
            // Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();
              
	      if(typev[ifile]==eZmumu)
		{
			  hU1vsZpt_rsp_MC_raw->Fill(dilep->Pt(),u1/dilep->Pt(),weight);
              hU1vsZpt_MC_raw->Fill(vDilepCor.Mod(),u1,weight);
              hU2vsZpt_MC_raw->Fill(dilep->Pt(),u2,weight);

              // hU1vsZpt_rsp_MC->Fill(dilep->Pt(),pU1/dilep->Pt(),weight);
              // hU1vsZpt_MC->Fill(dilep->Pt(),pU1,weight);
              // hU2vsZpt_MC->Fill(dilep->Pt(),pU2,weight);
			 hU1vsZpt_rsp_MC->Fill(vDilepCor.Mod(),pU1/vDilepCor.Mod(),weight*recoilWeight);
            hU1vsZpt_MC->Fill(vDilepCor.Mod(),pU1,weight*recoilWeight);
            hU2vsZpt_MC->Fill(vDilepCor.Mod(),pU2,weight*recoilWeight);
			
		  // hWmunuMetp->Fill(abs(pU2),weight*w2);
		  // hWmunuRecoilp->Fill(abs(pU1),weight*w2);			
		  hWmunuMetp->Fill(corrMetWithLepton,weight*w2*recoilWeight);
		  hWmunuRecoilp->Fill(pU,weight*w2*recoilWeight);
          hWmunuZpt->Fill(genVPt,weight*w2*recoilWeight);
		  corrMet=met, corrMetPhi=metPhi;
		}
	      else
		{
		  hEWKMetp->Fill(corrMetWithLepton,weight);
		  hEWKRecoilp->Fill(pU,weight);
          hEWKZpt->Fill(vDilepCor.Mod(),weight);
		  corrMet=met, corrMetPhi=metPhi;
		}
              //hWmunuMetp->Fill(corrMet,weight*w2); 
              //hWmunuMetp->Fill(pU,weight*w2); 
              corrMet=met, corrMetPhi=metPhi;
            }
	  
          hWmunuMet_RecoilUp->Fill(corrMetWithLepton,weight);
          if(q1*q2 < 0) 
          {
            pU1 = 0; pU2 = 0; 
//             recoilCorr->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,dl.Pt(),dl.Phi(),pU1,pU2,4,4,0);
//             recoilCorr->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,dl.Pt(),dl.Phi(),pU1,pU2,10,10,0);
            hWmunuMetp_RecoilUp->Fill(corrMet,weight); 
            corrMet=met, corrMetPhi=metPhi;
          }
	    
          hWmunuMet_RecoilDown->Fill(corrMet,weight);
          if(q1*q2<0) 
          {
            pU1 = 0; pU2 = 0; 
//             recoilCorr->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,dl.Pt(),dl.Phi(),pU1,pU2,-4,-4,0);
//             recoilCorr->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,dl.Pt(),dl.Phi(),pU1,pU2,-10,-10,0);
            hWmunuMetp_RecoilDown->Fill(corrMet,weight);
            corrMet=met, corrMetPhi=metPhi;
          }
			
		}
        if(typev[ifile]==eAntiWmunu) {
          Double_t corrMet=met, corrMetPhi=metPhi;
          hAntiWmunuMet->Fill(corrMet,weight);
          if(q1*q2 < 0) { hAntiWmunuMetp->Fill(corrMet,weight); } 
          else    { hAntiWmunuMetm->Fill(corrMet,weight); }
        }
        if(typev[ifile]==eEWK) {
	  // met and dilepton do not have muon/scale correction
		   // TVector2 vDilepCor = vLepCor1 + vLepCor2;
           // Double_t corrMetWithLepton = (vMetCorr + vDilepRaw - vDilepCor).Mod();
           // Double_t corrMetPhiLepton = (vMetCorr + vDilepRaw - vDilepCor).Phi();
	      double pUX  = corrMetWithLepton*cos(corrMetPhiLepton) + vDilepCor.Mod()*cos(vDilepCor.Phi());
	      double pUY  = corrMetWithLepton*sin(corrMetPhiLepton) + vDilepCor.Mod()*sin(vDilepCor.Phi());
	      double pU   = sqrt(pUX*pUX+pUY*pUY);
          // if(mu1.Pt()        < PT_CUT)  continue;
          // if(mu2.Pt()        < PT_CUT)  continue;
          // if(dilep->M()        < MASS_LOW)  continue;
          // if(dilep->M()        > MASS_HIGH) continue;
          hEWKMet->Fill(corrMetWithLepton,weight);
          if(q1*q2<0) { 
			  hEWKRecoilp->Fill(pU,weight); 
			  hEWKMetp->Fill(corrMetWithLepton,weight); 
			  hEWKZpt->Fill(vDilepCor.Mod(),weight); 
		  }
          else    { hEWKMetm->Fill(met,weight); }
        }
        if(typev[ifile]==eAntiEWK) {
          hAntiEWKMet->Fill(met,weight);
          if(q1*q2 < 0) { hAntiEWKMetp->Fill(met,weight); }
          else    { hAntiEWKMetm->Fill(met,weight); }
        }
      }
    }
  }
  delete infile;
  infile=0, intree=0;   
  
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  //
  RooRealVar nSig("nSig","nSig",1.0*(hDataMet->Integral()),0.7*hDataMet->Integral(),1.5*hDataMet->Integral());
  RooRealVar nQCD("nQCD","nQCD",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar cewk("cewk","cewk",0.1,0,5) ;
  cewk.setVal(hEWKMet->Integral()/hWmunuMet->Integral());
  cewk.setConstant(kTRUE);
  RooFormulaVar nEWK("nEWK","nEWK","cewk*nSig",RooArgList(nSig,cewk));
  RooRealVar nAntiSig("nAntiSig","nAntiSig",0.05*(hAntiDataMet->Integral()),0,hAntiDataMet->Integral());
  RooRealVar nAntiQCD("nAntiQCD","nAntiQCD",0.9*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar dewk("dewk","dewk",0.1,0,5) ;
  dewk.setVal(hAntiEWKMet->Integral()/hAntiWmunuMet->Integral());
  dewk.setConstant(kTRUE);
  RooFormulaVar nAntiEWK("nAntiEWK","nAntiEWK","dewk*nAntiSig",RooArgList(nAntiSig,dewk));
  
  // RooRealVar nSigp("nSigp","nSigp",1.0*(hWmunuMetp->Integral()),0.7*hDataMetp->Integral(),1.5*hDataMetp->Integral());
  RooRealVar nSigp("nSigp","nSigp",hDataMetp->Integral(),0.7*hDataMetp->Integral(),1.5*hDataMetp->Integral());
  //RooRealVar nQCDp("nQCDp","nQCDp",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nQCDp("nQCDp","nQCDp",25000.,0,hDataMetp->Integral());
  RooRealVar cewkp("cewkp","cewkp",0.01,0,5) ;
  cewkp.setVal(hEWKMetp->Integral()/hWmunuMetp->Integral());
  // cewkp.setVal(0);
  cewkp.setConstant(kTRUE);
  RooFormulaVar nEWKp("nEWKp","nEWKp","cewkp*nSigp",RooArgList(nSigp,cewkp));
  RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",0.05*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar dewkp("dewkp","dewkp",0.1,0,5) ;
  dewkp.setVal(hAntiEWKMetp->Integral()/hAntiWmunuMetp->Integral());
  dewkp.setConstant(kTRUE);
  RooFormulaVar nAntiEWKp("nAntiEWKp","nAntiEWKp","dewkp*nAntiSigp",RooArgList(nAntiSigp,dewkp));
  
  RooRealVar nSigRp("nSigRp","nSigRp",1.0*(hWmunuRecoilp->Integral()),0.7*hDataRecoilp->Integral(),1.5*hDataRecoilp->Integral());
  //RooRealVar nQCDp("nQCDp","nQCDp",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nQCDRp("nQCDRp","nQCDRp",0.01*hDataRecoilp->Integral(),0,hDataRecoilp->Integral());
  RooRealVar cewkRp("cewkRp","cewkRp",0.01,0,5) ;
  cewkRp.setVal(hEWKRecoilp->Integral()/hWmunuRecoilp->Integral());
  // cewkRp.setVal(0);
  cewkRp.setConstant(kTRUE);
  RooFormulaVar nEWKRp("nEWKRp","nEWKRp","cewkRp*nSigRp",RooArgList(nSigRp,cewkRp));
  
    RooRealVar nSigZp("nSigZp","nSigZp",1.0*(hWmunuZpt->Integral()),0.7*hDataZpt->Integral(),1.5*hDataZpt->Integral());
  //RooRealVar nQCDp("nQCDp","nQCDp",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nQCDZp("nQCDZp","nQCDZp",0.01*hDataZpt->Integral(),0,hDataZpt->Integral());
  RooRealVar cewkZp("cewkZp","cewkZp",0.01,0,5) ;
  // cewkZp.setVal(hEWKZpt->Integral()/hWmunuZpt->Integral());
  cewkZp.setVal(0);
  cewkZp.setConstant(kTRUE);
  RooFormulaVar nEWKZp("nEWKZp","nEWKZp","cewkZp*nSigZp",RooArgList(nSigZp,cewkZp));
  
  
  RooRealVar nSigm("nSigm","nSigm",1.0*(hDataMetm->Integral()),0,hDataMetm->Integral());
  //RooRealVar nQCDm("nQCDm","nQCDm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  RooRealVar nQCDm("nQCDm","nQCDm",25000,0,hDataMetm->Integral());
  RooRealVar cewkm("cewkm","cewkm",0.1,0,5) ;
  cewkm.setVal(hEWKMetm->Integral()/hWmunuMetm->Integral());
  cewkm.setConstant(kTRUE);
  RooFormulaVar nEWKm("nEWKm","nEWKm","cewkm*nSigm",RooArgList(nSigm,cewkm));  
  RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",0.05*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  RooRealVar dewkm("dewkm","dewkm",0.1,0,5) ;
  dewkm.setVal(hAntiEWKMetm->Integral()/hAntiWmunuMetm->Integral());
  dewkm.setConstant(kTRUE);
  RooFormulaVar nAntiEWKm("nAntiEWKm","nAntiEWKm","dewkm*nAntiSigm",RooArgList(nAntiSigm,dewkm));

  //
  // Construct PDFs for fitting
  //
  RooRealVar pfmet("pfmet","pfmet",0,METMAX);
  pfmet.setBins(NBINS);
   
  // Signal PDFs
  RooDataHist wmunuMet ("wmunuMET", "wmunuMET", RooArgSet(pfmet),hWmunuMet);  RooHistPdf pdfWm ("wm", "wm", pfmet,wmunuMet, 1);
  RooDataHist wmunuMetp("wmunuMETp","wmunuMETp",RooArgSet(pfmet),hWmunuMetp); RooHistPdf pdfWmp("wmp","wmp",pfmet,wmunuMetp,1);
  RooDataHist wmunuMetm("wmunuMETm","wmunuMETm",RooArgSet(pfmet),hWmunuMetm); RooHistPdf pdfWmm("wmm","wmm",pfmet,wmunuMetm,1); 
  RooDataHist wmunuMet_RecoilUp("wmunuMET_RecoilUp", "wmunuMET_RecoilUp", RooArgSet(pfmet),hWmunuMet_RecoilUp);  RooHistPdf pdfWm_RecoilUp("wm_RecoilUp", "wm_RecoilUp", pfmet,wmunuMet_RecoilUp, 1);
  RooDataHist wmunuMetp_RecoilUp("wmunuMETp_RecoilUp","wmunuMETp_RecoilUp",RooArgSet(pfmet),hWmunuMetp_RecoilUp); RooHistPdf pdfWmp_RecoilUp("wmp_RecoilUp","wmp_RecoilUp",pfmet,wmunuMetp_RecoilUp,1);
  RooDataHist wmunuMetm_RecoilUp("wmunuMETm_RecoilUp","wmunuMETm_RecoilUp",RooArgSet(pfmet),hWmunuMetm_RecoilUp); RooHistPdf pdfWmm_RecoilUp("wmm_RecoilUp","wmm_RecoilUp",pfmet,wmunuMetm_RecoilUp,1); 
  RooDataHist wmunuMet_RecoilDown("wmunuMET_RecoilDown", "wmunuMET_RecoilDown", RooArgSet(pfmet),hWmunuMet_RecoilDown);  RooHistPdf pdfWm_RecoilDown("wm_RecoilDown", "wm_RecoilDown", pfmet,wmunuMet_RecoilDown, 1);
  RooDataHist wmunuMetp_RecoilDown("wmunuMETp_RecoilDown","wmunuMETp_RecoilDown",RooArgSet(pfmet),hWmunuMetp_RecoilDown); RooHistPdf pdfWmp_RecoilDown("wmp_RecoilDown","wmp_RecoilDown",pfmet,wmunuMetp_RecoilDown,1);
  RooDataHist wmunuMetm_RecoilDown("wmunuMETm_RecoilDown","wmunuMETm_RecoilDown",RooArgSet(pfmet),hWmunuMetm_RecoilDown); RooHistPdf pdfWmm_RecoilDown("wmm_RecoilDown","wmm_RecoilDown",pfmet,wmunuMetm_RecoilDown,1); 
   RooDataHist wmunuMet_ScaleUp("wmunuMET_ScaleUp", "wmunuMET_ScaleUp", RooArgSet(pfmet),hWmunuMet_ScaleUp);  RooHistPdf pdfWm_ScaleUp("wm_ScaleUp", "wm_ScaleUp", pfmet,wmunuMet_ScaleUp, 1);
  RooDataHist wmunuMetp_ScaleUp("wmunuMETp_ScaleUp","wmunuMETp_ScaleUp",RooArgSet(pfmet),hWmunuMetp_ScaleUp); RooHistPdf pdfWmp_ScaleUp("wmp_ScaleUp","wmp_ScaleUp",pfmet,wmunuMetp_ScaleUp,1);
  RooDataHist wmunuMetm_ScaleUp("wmunuMETm_ScaleUp","wmunuMETm_ScaleUp",RooArgSet(pfmet),hWmunuMetm_ScaleUp); RooHistPdf pdfWmm_ScaleUp("wmm_ScaleUp","wmm_ScaleUp",pfmet,wmunuMetm_ScaleUp,1); 
  RooDataHist wmunuMet_ScaleDown("wmunuMET_ScaleDown", "wmunuMET_ScaleDown", RooArgSet(pfmet),hWmunuMet_ScaleDown);  RooHistPdf pdfWm_ScaleDown("wm_ScaleDown", "wm_ScaleDown", pfmet,wmunuMet_ScaleDown, 1);
  RooDataHist wmunuMetp_ScaleDown("wmunuMETp_ScaleDown","wmunuMETp_ScaleDown",RooArgSet(pfmet),hWmunuMetp_ScaleDown); RooHistPdf pdfWmp_ScaleDown("wmp_ScaleDown","wmp_ScaleDown",pfmet,wmunuMetp_ScaleDown,1);
  RooDataHist wmunuMetm_ScaleDown("wmunuMETm_ScaleDown","wmunuMETm_ScaleDown",RooArgSet(pfmet),hWmunuMetm_ScaleDown); RooHistPdf pdfWmm_ScaleDown("wmm_ScaleDown","wmm_ScaleDown",pfmet,wmunuMetm_ScaleDown,1); 
  
  
  RooDataHist wmunuRecoilp("wmunuRecoilp","wmunuRecoilp",RooArgSet(pfmet),hWmunuRecoilp); RooHistPdf pdfWmpRecoil("wmprec","wmprec",pfmet,wmunuRecoilp,1);
  RooDataHist wmunuZpt("wmunuZpt","wmunuZpt",RooArgSet(pfmet),hWmunuZpt); RooHistPdf pdfWmZpt("wmzpt","wmzpt",pfmet,wmunuZpt,1);
  
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", pfmet,ewkMet, 1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",pfmet,ewkMetp,1); 
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",pfmet,ewkMetm,1); 

  RooDataHist ewkRecoilp("ewkRecoilp","ewkRecoilp",RooArgSet(pfmet),hEWKRecoilp); RooHistPdf pdfEWKpRecoil("ewkprec","ewkprec",pfmet,ewkRecoilp,1); 
  RooDataHist ewkZpt("ewkZpt","ewkZpt",RooArgSet(pfmet),hEWKZpt); RooHistPdf pdfEWKZpt("pdfEWKZpt","pdfEWKZpt",pfmet,ewkZpt,1); 
  // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK),   RooArgList(nSig,nEWK));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp),RooArgList(nSigp,nEWKp));
  // RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp),RooArgList(nSigp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm),RooArgList(nSigm,nEWKm));
  
  
  RooAddPdf pdfRecoilp("pdRecoilp","pdfRecoilp",RooArgList(pdfWmpRecoil,pdfEWKpRecoil),RooArgList(nSigRp,nEWKRp));
  // RooAddPdf pdfZpt("pdfZpt","pdfZpt",RooArgList(pdfWmZpt,pdfEWKZpt),RooArgList(nSigZp,nEWKZp)); 
  // RooAddPdf pdfRecoilp("pdRecoilp","pdfRecoilp",RooArgList(pdfWmpRecoil),RooArgList(nSigRp));
  RooAddPdf pdfZpt("pdfZpt","pdfZpt",RooArgList(pdfWmZpt),RooArgList(nSigZp));
 
  //
  // Perform fits
  //

  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet), hDataMet);
  RooDataHist dataMetp("dataMetp", "dataMetp", RooArgSet(pfmet), hDataMetp);
  RooDataHist dataMetm("dataMetm", "dataMetm", RooArgSet(pfmet), hDataMetm);
  
  RooDataHist dataRecoilp("dataRecoilp", "dataRecoilp", RooArgSet(pfmet), hDataRecoilp);
  RooDataHist dataZpt("dataZpt", "dataZpt", RooArgSet(pfmet), hDataZpt);

  cout << "Starting values for Wmunu yields: " << endl;
  cout << "   sig: " << hWmunuMet->Integral() << endl;
  cout << "   EWK: " << hEWKMet->Integral() << endl;
  cout << "   qcd: " << hDataMet->Integral()-hWmunuMet->Integral()-hEWKMet->Integral() << endl;

  cout << "Starting values for Wmunu_p yields: " << endl;
  cout << "   sig: " << hWmunuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;
  cout << "   qcd: " << hDataMetp->Integral()-hWmunuMetp->Integral()-hEWKMetp->Integral() << endl;

  cout << "Starting values for Wmunu_m yields: " << endl;
  cout << "   sig: " << hWmunuMetm->Integral() << endl;
  cout << "   EWK: " << hEWKMetm->Integral() << endl;
  cout << "   qcd: " << hDataMetm->Integral()-hWmunuMetm->Integral()-hEWKMetm->Integral() << endl;

  RooWorkspace combine_workspace("combine_workspace");
  combine_workspace.import(dataMet);
  combine_workspace.import(dataMetp);
  combine_workspace.import(dataRecoilp);
  combine_workspace.import(dataMetm);

  combine_workspace.import(pdfWm);
  combine_workspace.import(pdfWmpRecoil);
  combine_workspace.import(pdfWmp);
  combine_workspace.import(pdfWmm);
  combine_workspace.import(pdfWm_RecoilUp);
  combine_workspace.import(pdfWmp_RecoilUp);
  combine_workspace.import(pdfWmm_RecoilUp);
  combine_workspace.import(pdfWm_RecoilDown);
  combine_workspace.import(pdfWmp_RecoilDown);
  combine_workspace.import(pdfWmm_RecoilDown);
  combine_workspace.import(pdfWm_ScaleUp);
  combine_workspace.import(pdfWmp_ScaleUp);
  combine_workspace.import(pdfWmm_ScaleUp);
  combine_workspace.import(pdfWm_ScaleDown);
  combine_workspace.import(pdfWmp_ScaleDown);
  combine_workspace.import(pdfWmm_ScaleDown);
  combine_workspace.import(pdfEWK);
  combine_workspace.import(pdfEWKp);
  combine_workspace.import(pdfEWKpRecoil);
  combine_workspace.import(pdfEWKm);
  
  char nname[100];
  sprintf(nname, "%s/Zmumu_pdfTemplates.root",CPlot::sOutDir.Data());
  combine_workspace.writeToFile(nname);

  //  RooFitResult *fitRes = pdfMet.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE));//dataTotal.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE));
  RooFitResult *fitResp = pdfMetp.fitTo(dataMetp,Extended(),Minos(kTRUE),Save(kTRUE)); 
  RooFitResult *fitRecoilp = pdfRecoilp.fitTo(dataRecoilp,Extended(),Minos(kTRUE),Save(kTRUE)); 
  RooFitResult *fitZpt = pdfZpt.fitTo(dataZpt,Extended(),Minos(kTRUE),Save(kTRUE)); 
  //  RooFitResult *fitResm = pdfMetm.fitTo(dataMetm,Extended(),Minos(kTRUE),Save(kTRUE));
 
  // Use histogram version of fitted PDFs to make ratio plots
  // (Will also use PDF histograms later for Chi^2 and KS tests)
  //
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", pfmet));
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  hMetDiff->SetMarkerStyle(kFullCircle);
  hMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", pfmet));
  std::cout << "nsig p " << nSigp.getVal()<< "  newp " << nEWKp.getVal()<< "  pdf integral  "  << hPdfMetp->Integral() << std::endl;
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal())/hPdfMetp->Integral());
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);
  
  TH1D *hPdfRecoilp = (TH1D*)(pdfRecoilp.createHistogram("hPdfRecoilp", pfmet));
  hPdfRecoilp->Scale((nSigRp.getVal()+nEWKRp.getVal())/hPdfRecoilp->Integral());
  TH1D *hRecoilpDiff = makeDiffHist(hDataRecoilp,hPdfRecoilp,"hRecoilpDiff");
  hRecoilpDiff->SetMarkerStyle(kFullCircle);
  hRecoilpDiff->SetMarkerSize(0.9);
  
  TH1D *hPdfZpt = (TH1D*)(pdfZpt.createHistogram("hPdfZpt", pfmet));
  hPdfZpt->Scale((nSigRp.getVal()+nEWKRp.getVal())/hPdfZpt->Integral());
  TH1D *hZptDiff = makeDiffHist(hDataZpt,hPdfZpt,"hZptDiff");
  hZptDiff->SetMarkerStyle(kFullCircle);
  hZptDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet));
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle); 
  hMetmDiff->SetMarkerSize(0.9);

 
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  
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
  gStyle->SetTitleOffset(1.100,"Y");
  TGaxis::SetMaxDigits(3);
  
  TCanvas *c2 = MakeCanvas("c2","c2",800,800);
  c2->Divide(1,0,0,0);
  c2->cd(1)->SetPad(0,0.3,1.0,1.0);
  c2->cd(1)->SetTopMargin(0.1);
  c2->cd(1)->SetBottomMargin(0.01);
  c2->cd(1)->SetLeftMargin(0.10);  
  c2->cd(1)->SetRightMargin(0.07);  
  c2->cd(1)->SetTickx(1);
  c2->cd(1)->SetTicky(1);  
//   c2->cd(2)->SetPad(0,0,1.0,0.3);
//   c2->cd(2)->SetTopMargin(0.05);
//   c2->cd(2)->SetBottomMargin(0.45);
//   c2->cd(2)->SetLeftMargin(0.15);
//   c2->cd(2)->SetRightMargin(0.07);
//   c2->cd(2)->SetTickx(1);
//   c2->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.400,"Y");
//   TGaxis::SetMaxDigits(3);
  TFile f("zPtDataMCresid_noEWK_2.root","RECREATE");
  c2->cd(1);
  TH1D* hZptRatio = (TH1D*) hDataZpt->Clone("hZptRatio");
  hZptRatio->Divide(hWmunuZpt);
  hZptRatio->Draw();
  hZptRatio->Write();
  f.Write();
  // f.
  // c2->SaveAs("zPtDataMCresid.root");
  
  
//   TProfile *hRespDataU1Prof = hResponseDataU1->ProfileX("hResponseDataU1_pfx",1,-1);
//   TProfile *hRespMCU1Prof = hResponseMCU1->ProfileX("hResponseMCU1_pfx",1,-1);
  
  std::vector<double> Z_pT, Z_pT_Err;
  
  std::vector<double> U1_Mean,U1_Sigma, U1_Mean_rsp, U1_Sigma_rsp;
  std::vector<double> U1_Mean_Err, U1_Sigma_Err, U1_Mean_Err_rsp, U1_Sigma_Err_rsp;
  std::vector<double> U2_Mean, U2_Sigma;
  std::vector<double> U2_Mean_Err, U2_Sigma_Err;

  std::vector<double> U1_Mean_MC,U1_Sigma_MC, U1_Mean_rsp_MC, U1_Sigma_rsp_MC;
  std::vector<double> U1_Mean_Err_MC, U1_Sigma_Err_MC, U1_Mean_Err_rsp_MC, U1_Sigma_Err_rsp_MC;
  std::vector<double> U2_Mean_MC, U2_Sigma_MC;
  std::vector<double> U2_Mean_Err_MC, U2_Sigma_Err_MC;
  
  std::vector<double> U1_Mean_MC_raw,U1_Sigma_MC_raw, U1_Mean_rsp_MC_raw, U1_Sigma_rsp_MC_raw;
  std::vector<double> U1_Mean_Err_MC_raw, U1_Sigma_Err_MC_raw, U1_Mean_Err_rsp_MC_raw, U1_Sigma_Err_rsp_MC_raw;
  std::vector<double> U2_Mean_MC_raw, U2_Sigma_MC_raw;
  std::vector<double> U2_Mean_Err_MC_raw, U2_Sigma_Err_MC_raw;
  
  fitWithGaus(hU1vsZpt_rsp, "data_puppiU1_rsp_", U1_Mean_rsp, U1_Mean_Err_rsp, U1_Sigma_rsp, U1_Sigma_Err_rsp);
  fitWithGaus(hU1vsZpt, "data_puppiU1_", U1_Mean, U1_Mean_Err, U1_Sigma, U1_Sigma_Err);
  fitWithGaus(hU2vsZpt, "data_puppiU2_", U2_Mean, U2_Mean_Err, U2_Sigma, U2_Sigma_Err);
    
  fitWithGaus(hU1vsZpt_rsp_MC, "mc_puppiU1_rsp_", U1_Mean_rsp_MC, U1_Mean_Err_rsp_MC, U1_Sigma_rsp_MC, U1_Sigma_Err_rsp_MC);
  fitWithGaus(hU1vsZpt_MC,     "mc_puppiU1_", U1_Mean_MC,     U1_Mean_Err_MC,     U1_Sigma_MC,     U1_Sigma_Err_MC);
  fitWithGaus(hU2vsZpt_MC,     "mc_puppiU2_", U2_Mean_MC,     U2_Mean_Err_MC,     U2_Sigma_MC,     U2_Sigma_Err_MC);
  
  fitWithGaus(hU1vsZpt_rsp_MC_raw, "mc_puppiU1_rsp_raw_", U1_Mean_rsp_MC_raw, U1_Mean_Err_rsp_MC_raw, U1_Sigma_rsp_MC_raw, U1_Sigma_Err_rsp_MC_raw);
  fitWithGaus(hU1vsZpt_MC_raw,     "mc_puppiU1_raw_", U1_Mean_MC_raw,     U1_Mean_Err_MC_raw,     U1_Sigma_MC_raw,     U1_Sigma_Err_MC_raw);
  fitWithGaus(hU2vsZpt_MC_raw,     "mc_puppiU2_raw_", U2_Mean_MC_raw,     U2_Mean_Err_MC_raw,     U2_Sigma_MC_raw,     U2_Sigma_Err_MC_raw);
  
  c->cd(1);
  c2->cd(1);
  TH2D *DivideU1 = (TH2D*) hU1vsZpt_MC->Clone("DivideU1");
  TH2D *DivideU2 = (TH2D*) hU2vsZpt_MC->Clone("DivideU2");
  DivideU1->Divide(hU1vsZpt_MC_raw);
  DivideU2->Divide(hU2vsZpt_MC_raw);

  DivideU1->Draw("colz");
  DivideU1->GetZaxis()->SetRangeUser(0,3);
  c->SaveAs((outputdir+"/MCRawCorrRatioU1.png").c_str());
  DivideU2->Draw("colz");
  DivideU2->GetZaxis()->SetRangeUser(0,3);
  c->SaveAs((outputdir+"/MCRawCorrRatioU2.png").c_str());
  
  for(int bin = 1; bin < hU1vsZpt->GetNbinsX()+1; bin++){
    Z_pT.push_back((bins[bin-1]+bins[bin])*0.5);
    Z_pT_Err.push_back((bins[bin]-bins[bin-1])*0.5);
  }
  
  TGraphErrors *response_U1      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U1_Mean_rsp[0], &Z_pT_Err[0], &U1_Mean_Err_rsp[0]);
  TGraphErrors *resolution_U1    = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U1_Sigma[0],&Z_pT_Err[0], &U1_Sigma_Err[0]);
  
  TGraphErrors *response_U2      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U2_Mean[0], &Z_pT_Err[0], &U2_Mean_Err[0]);
  TGraphErrors *resolution_U2    = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U2_Sigma[0],&Z_pT_Err[0], &U2_Sigma_Err[0]);
  
  TGraphErrors *response_U1_MC      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U1_Mean_rsp_MC[0], &Z_pT_Err[0], &U1_Mean_Err_rsp_MC[0]);
  TGraphErrors *resolution_U1_MC    = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U1_Sigma_MC[0],&Z_pT_Err[0], &U1_Sigma_Err_MC[0]);
  
  TGraphErrors *response_U2_MC      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U2_Mean_MC[0], &Z_pT_Err[0], &U2_Mean_Err_MC[0]);
  TGraphErrors *resolution_U2_MC    = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U2_Sigma_MC[0],&Z_pT_Err[0], &U2_Sigma_Err_MC[0]);
  
  TGraphErrors *response_U1_MC_raw      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U1_Mean_rsp_MC_raw[0], &Z_pT_Err[0], &U1_Mean_Err_rsp_MC_raw[0]);
  TGraphErrors *resolution_U1_MC_raw    = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U1_Sigma_MC_raw[0],    &Z_pT_Err[0], &U1_Sigma_Err_MC_raw[0]);
  
  TGraphErrors *response_U2_MC_raw      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U2_Mean_MC_raw[0], &Z_pT_Err[0], &U2_Mean_Err_MC_raw[0]);
  TGraphErrors *resolution_U2_MC_raw    = new TGraphErrors(Z_pT.size(), &Z_pT[0], &U2_Sigma_MC_raw[0],&Z_pT_Err[0], &U2_Sigma_Err_MC_raw[0]);
  
  
  std::vector<double> CorrU1Sigma; std::vector<double> CorrU1SigmaErr;
  std::vector<double> CorrU1Sigma_MC; std::vector<double> CorrU1SigmaErr_MC;
  std::vector<double> CorrU1Sigma_MC_raw; std::vector<double> CorrU1SigmaErr_MC_raw;
  
  std::vector<double> CorrU2Sigma; std::vector<double> CorrU2SigmaErr;
  std::vector<double> CorrU2Sigma_MC; std::vector<double> CorrU2SigmaErr_MC;
  std::vector<double> CorrU2Sigma_MC_raw; std::vector<double> CorrU2SigmaErr_MC_raw;
  
  for(int i = 0; i < Z_pT.size(); ++i){
    CorrU1Sigma.push_back(U1_Sigma[i]/U1_Mean_rsp[i]);
    CorrU1SigmaErr.push_back(U1_Sigma_Err[i]/U1_Mean_rsp[i]);
    
    CorrU1Sigma_MC.push_back(U1_Sigma_MC[i]/U1_Mean_rsp_MC[i]);
    CorrU1SigmaErr_MC.push_back(U1_Sigma_Err_MC[i]/U1_Mean_rsp_MC[i]);
    
    CorrU1Sigma_MC_raw.push_back(U1_Sigma_MC_raw[i]/U1_Mean_rsp_MC_raw[i]);
    CorrU1SigmaErr_MC_raw.push_back(U1_Sigma_Err_MC_raw[i]/U1_Mean_rsp_MC_raw[i]);
    
    CorrU2Sigma.push_back(U2_Sigma[i]/U1_Mean_rsp[i]);
    CorrU2SigmaErr.push_back(U2_Sigma_Err[i]/U1_Mean_rsp[i]);
    
    CorrU2Sigma_MC.push_back(U2_Sigma_MC[i]/U1_Mean_rsp_MC[i]);
    CorrU2SigmaErr_MC.push_back(U2_Sigma_Err_MC[i]/U1_Mean_rsp_MC[i]);
    
    CorrU2Sigma_MC_raw.push_back(U2_Sigma_MC_raw[i]/U1_Mean_rsp_MC_raw[i]);
    CorrU2SigmaErr_MC_raw.push_back(U2_Sigma_Err_MC_raw[i]/U1_Mean_rsp_MC_raw[i]);
  }
  
  TGraphErrors *corrResp_U1      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &CorrU1Sigma[0], &Z_pT_Err[0], &CorrU1SigmaErr[0]);
  TGraphErrors *corrResp_U2      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &CorrU2Sigma[0], &Z_pT_Err[0], &CorrU2SigmaErr[0]);
  
  TGraphErrors *corrResp_U1_MC      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &CorrU1Sigma_MC[0], &Z_pT_Err[0], &CorrU1SigmaErr_MC[0]);
  TGraphErrors *corrResp_U2_MC      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &CorrU2Sigma_MC[0], &Z_pT_Err[0], &CorrU2SigmaErr_MC[0]);
  
  TGraphErrors *corrResp_U1_MC_raw      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &CorrU1Sigma_MC_raw[0], &Z_pT_Err[0], &CorrU1SigmaErr_MC_raw[0]);
  TGraphErrors *corrResp_U2_MC_raw      = new TGraphErrors(Z_pT.size(), &Z_pT[0], &CorrU2Sigma_MC_raw[0], &Z_pT_Err[0], &CorrU2SigmaErr_MC_raw[0]);
  
  
  // Plot resolution of parallel resolution component
  CPlot plot_resolution_u1("resolution_par","","p_{T}(Z) [GeV]","#sigma(U_{#parallel}) [GeV]");
  plot_resolution_u1.AddGraph(resolution_U1        ,"Data"    ,"", kBlack );
  plot_resolution_u1.AddGraph(resolution_U1_MC_raw ,"MC w/o Recoil Corr"   ,"", kBlue   ,20);
  plot_resolution_u1.AddGraph(resolution_U1_MC     , "MC w/ Recoil Corr" ,"", kRed ,34);
//   plot_resolution_u1.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  plot_resolution_u1.SetYRange(0.0,40);
  plot_resolution_u1.Draw(c2,kTRUE,format);

  // Plot resolution of perpendicular resolution component
  CPlot plot_resolution_u2("resolution_prp","","p_{T}(Z) [GeV]","#sigma(U_{#perp}  ) [GeV]");
  plot_resolution_u2.AddGraph(resolution_U2        ,"Data"    ,"", kBlack ); // Draw data with closed circles
  plot_resolution_u2.AddGraph(resolution_U2_MC_raw ,"MC w/o Recoil Corr"   ,"", kBlue   ,20);
  plot_resolution_u2.AddGraph(resolution_U2_MC     ,"MC w/ Recoil Corr"  ,"", kRed ,34); // Draw MC with open circles
//   plot_resolution_u2.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  plot_resolution_u2.SetYRange(0.0,30.0);
  plot_resolution_u2.Draw(c2,kTRUE,format);
  
  // Plot response of parallel resolution component
  CPlot plot_response_u1("response_par","","p_{T}(Z) [GeV]","<U_{#parallel} / Z p_{T} >");
  plot_response_u1.AddGraph(response_U1        ,"Data"  ,"", kBlack ); // Draw data with closed circles
  plot_response_u1.AddGraph(response_U1_MC_raw ,"MC w/o Recoil Corr" ,"", kBlue   ,20);
  plot_response_u1.AddGraph(response_U1_MC     ,"MC w/ Recoil Corr","", kRed ,34); // Draw MC with open circles
//   plot_response_u1.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  plot_response_u1.SetYRange(0.0,1.1);
  plot_response_u1.Draw(c2,kTRUE,format);
  
  // Plot response of perpendicular resolution component
  CPlot plot_response_u2("response_prp","","p_{T}(Z) [GeV]","<U_{#perp}  >");
  plot_response_u2.AddGraph(response_U2        ,"Data"  ,"", kBlack ); // Draw data with closed circles
  plot_response_u2.AddGraph(response_U2_MC_raw ,"MC w/o Recoil Corr" ,"", kBlue   ,20);
  plot_response_u2.AddGraph(response_U2_MC     ,"MC w/ Recoil Corr","", kRed ,34); // Draw MC with open circles
//   plot_response_u2.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  plot_response_u2.SetYRange(-1.0, 3.0);
  plot_response_u2.Draw(c2,kTRUE,format);
  
    // Plot response of perpendicular resolution component
  CPlot corrResp_par("corrResp_par","","p_{T}(Z) [GeV]","#sigma(U_{#parallel})/<U_{#parallel} / Z p_{T} > [GeV]");
  corrResp_par.AddGraph(corrResp_U1        ,"Data"  ,"", kBlack ); // Draw data with closed circles
  corrResp_par.AddGraph(corrResp_U1_MC_raw ,"MC w/o Recoil Corr" ,"", kBlue   ,20);
  corrResp_par.AddGraph(corrResp_U1_MC     ,"MC w/ Recoil Corr","", kRed ,34); // Draw MC with open circles
//   corrResp_par.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  corrResp_par.SetYRange(0.0, 50.0);
  corrResp_par.Draw(c2,kTRUE,format);
  
    // Plot response of perpendicular resolution component
  CPlot corrResp_prp("corrResp_prp","","p_{T}(Z) [GeV]","#sigma(U_{#perp}  )/<U_{#parallel} / Z p_{T}  > [GeV]");
  corrResp_prp.AddGraph(corrResp_U2        ,"Data"  ,"", kBlack ); // Draw data with closed circles
  corrResp_prp.AddGraph(corrResp_U2_MC_raw ,"MC w/o Recoil Corr" ,"", kBlue   ,20);
  corrResp_prp.AddGraph(corrResp_U2_MC     ,"MC w/ Recoil Corr","", kRed ,34); // Draw MC with open circles
//   corrResp_prp.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  corrResp_prp.SetYRange(0.0, 50.0);
  corrResp_prp.Draw(c2,kTRUE,format);
  
  
  char name1[100]; char name2[100];
  // sprintf(name1,"PoorMansRecoilCorr_sec%s.root",input_section.Data());
  sprintf(name1,"PoorMansRecoilCorr_blah.root");
   TFile fcorr(name1,"RECREATE");
  // loop through all the bins of the 2-d plot, project out all of the individual bins, and plot them on top of each other
  for (int i = 1; i <= hU2vsZpt->GetNbinsX(); i++){
	  double mean=0; double rms=0;
	  sprintf(name2,"%f < Z pT < %f",hU2vsZpt->GetXaxis()->GetBinLowEdge(i),hU2vsZpt->GetXaxis()->GetBinLowEdge(i+1));
	  
	  sprintf(name1,"resolution_U1_bin%d",i);
	  TH1D* hU1vsZpt_pfy = hU1vsZpt->ProjectionY("hU1vsZpt_pfy",i,i);
	  TH1D* hU1vsZpt_MC_pfy = hU1vsZpt_MC->ProjectionY("hU1vsZpt_MC_pfy",i,i);
	  TH1D* hU1vsZpt_MC_raw_pfy = hU1vsZpt_MC_raw->ProjectionY("hU1vsZpt_MC_raw_pfy",i,i);
	  mean=hU1vsZpt_pfy->GetMean(); rms=hU1vsZpt_pfy->GetRMS();
	  hU1vsZpt_pfy->GetXaxis()->SetRangeUser(mean-4*rms,mean+4*rms);
	  hU1vsZpt_MC_pfy->GetXaxis()->SetRangeUser(mean-4*rms,mean+4*rms);
	  hU1vsZpt_MC_raw_pfy->GetXaxis()->SetRangeUser(mean-4*rms,mean+4*rms);
	  hU1vsZpt_pfy->Scale(1/hU1vsZpt_pfy->Integral());
	  hU1vsZpt_MC_pfy->Scale(1/hU1vsZpt_MC_pfy->Integral());
	  hU1vsZpt_MC_raw_pfy->Scale(1/hU1vsZpt_MC_raw_pfy->Integral());
	  CPlot corrResp_par(name1,name2,"U_{#parallel}","Nevents");
      corrResp_par.AddHist1D(hU1vsZpt_pfy        ,"Data"  ,"", kBlack ); // Draw data with closed circles
      corrResp_par.AddHist1D(hU1vsZpt_MC_raw_pfy ,"MC w/o Recoil Corr" ,"", kBlue   ,20);
      corrResp_par.AddHist1D(hU1vsZpt_MC_pfy     ,"MC w/ Recoil Corr","", kRed ,34); // Draw MC with open circles
      // corrResp_prp.SetYRange(0.0, 50.0);
      corrResp_par.Draw(c2,kTRUE,format);
      sprintf(name1,"hRecoilU1_bin%d",i);
      TH1D *hU1Diff = (TH1D*) hU1vsZpt_pfy->Clone(name1);
      hU1Diff->Divide(hU1vsZpt_MC_raw_pfy);
      hU1Diff->Write();
      std::cout << "bin" << i << "  raw mean " << hU1vsZpt_MC_raw_pfy->GetMean() << "  raw rms " << hU1vsZpt_MC_raw_pfy->GetRMS() << std::endl;
      std::cout << "bin" << i << "  par mean " << hU1vsZpt_MC_pfy->GetMean() << "  par rms " << hU1vsZpt_MC_pfy->GetRMS() << std::endl;
      std::cout << "bin" << i << "  dat mean " << hU1vsZpt_pfy->GetMean() << "  dat rms " << hU1vsZpt_pfy->GetRMS() << std::endl;
 
	  
	  delete hU1vsZpt_pfy; delete hU1vsZpt_MC_pfy; delete hU1vsZpt_MC_raw_pfy;
      delete hU1Diff;
	  
      sprintf(name1,"resolution_U2_bin%d",i);
	  TH1D* hU2vsZpt_pfy = hU2vsZpt->ProjectionY("hU2vsZpt_pfy",i,i);
	  TH1D* hU2vsZpt_MC_pfy = hU2vsZpt_MC->ProjectionY("hU2vsZpt_MC_pfy",i,i);
	  TH1D* hU2vsZpt_MC_raw_pfy = hU2vsZpt_MC_raw->ProjectionY("hU2vsZpt_MC_raw_pfy",i,i);
	  mean=hU2vsZpt_pfy->GetMean(); rms=hU2vsZpt_pfy->GetRMS();
	  hU2vsZpt_pfy->GetXaxis()->SetRangeUser(mean-4*rms,mean+4*rms);
	  hU2vsZpt_MC_pfy->GetXaxis()->SetRangeUser(mean-4*rms,mean+4*rms);
	  hU2vsZpt_MC_raw_pfy->GetXaxis()->SetRangeUser(mean-4*rms,mean+4*rms);
	  CPlot corrResp_prp(name1,name2,"U_{#perp}","Nevents");
	  hU2vsZpt_pfy->Scale(1/hU2vsZpt_pfy->Integral());
	  hU2vsZpt_MC_pfy->Scale(1/hU2vsZpt_MC_pfy->Integral());
	  hU2vsZpt_MC_raw_pfy->Scale(1/hU2vsZpt_MC_raw_pfy->Integral());
      corrResp_prp.AddHist1D(hU2vsZpt_pfy        ,"Data"  ,"", kBlack ); // Draw data with closed circles
      corrResp_prp.AddHist1D(hU2vsZpt_MC_raw_pfy ,"MC w/o Recoil Corr" ,"", kBlue   ,20);
      corrResp_prp.AddHist1D(hU2vsZpt_MC_pfy     ,"MC w/ Recoil Corr","", kRed ,34); // Draw MC with open circles
      // corrResp_prp.SetYRange(0.0, 50.0);
      corrResp_prp.Draw(c2,kTRUE,format);
      
      sprintf(name1,"hRecoilU2_bin%d",i);
      TH1D *hU2Diff = (TH1D*) hU2vsZpt_pfy->Clone(name1);
      hU2Diff->Divide(hU2vsZpt_MC_raw_pfy);
      hU2Diff->Write();
	  
	  delete hU2vsZpt_pfy; delete hU2vsZpt_MC_pfy; delete hU2vsZpt_MC_raw_pfy;
      delete hU2Diff;
      
      
	  
  }
  
  fcorr.Write();
//   return;
  char ylabel[100];  // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 8 TeV",lumi*1000.);
  else         sprintf(lumitext,"%.0f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);
  
  // plot colors
  Int_t linecolorW   = kOrange-3;
  Int_t fillcolorW   = kOrange-2;
  Int_t linecolorEWK = kOrange+10;
  Int_t fillcolorEWK = kOrange+7;
  Int_t linecolorQCD = kViolet+2;
  Int_t fillcolorQCD = kViolet-5;
  Int_t ratioColor   = kGray+2;
  
  //
  // Dummy histograms for TLegend
  // (I can't figure out how to properly pass RooFit objects...)
  //
  TH1D *hDummyData = new TH1D("hDummyData","",0,0,10);
  hDummyData->SetMarkerStyle(kFullCircle);
  hDummyData->SetMarkerSize(0.9);
  
  TH1D *hDummyW = new TH1D("hDummyW","",0,0,10);
  hDummyW->SetLineColor(linecolorW);
  hDummyW->SetFillColor(fillcolorW);
  hDummyW->SetFillStyle(1001);
  
  TH1D *hDummyEWK = new TH1D("hDummyEWK","",0,0,10);
  hDummyEWK->SetLineColor(linecolorEWK);
  hDummyEWK->SetFillColor(fillcolorEWK);
  hDummyEWK->SetFillStyle(1001);
   
  //
  // W MET plot
  //
  RooPlot *wmframe = pfmet.frame(Bins(NBINS)); 
  wmframe->GetYaxis()->SetNdivisions(505);
  dataMet.plotOn(wmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMet.plotOn(wmframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMet.plotOn(wmframe,LineColor(linecolorW));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfEWK)),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfEWK)),LineColor(linecolorEWK));
  //pdfMet.plotOn(wmframe,Components(RooArgSet(*(qcd.model))),FillColor(fillcolorQCD),DrawOption("F"));
  //pdfMet.plotOn(wmframe,Components(RooArgSet(*(qcd.model))),LineColor(linecolorQCD));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfWm)),LineColor(linecolorW),LineStyle(2));
  dataMet.plotOn(wmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMet->GetBinWidth(1));
  CPlot plotMet("fitmet",wmframe,"","",ylabel);
  plotMet.SetLegend(0.68,0.57,0.93,0.77);
  plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMet.GetLegend()->AddEntry(hDummyW,"Z#rightarrow#mu#mu","F");
  plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  //plotMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMet.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotMet.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  plotMet.SetYRange(0.1,1.1*(hDataMet->GetMaximum()));
  plotMet.Draw(c,kFALSE,format,1);

  CPlot plotMetDiff("fitmet","","#slash{E}_{T} [GeV]","#chi");
  //CPlot plotMetDiff("fitmet","","mT [GeV]","#chi");
  plotMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
  plotMetDiff.SetYRange(-8,8);
  plotMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotMetDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotMetDiff.Draw(c,kTRUE,format,2);
  
  plotMet.SetName("fitmetlog");
  plotMet.SetLogy();
  plotMet.SetYRange(1e-4*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
  //  plotMet.SetYRange(1e-3*(hDataMet->GetMaximum()),10*(hDataMet->GetMaximum()));
  plotMet.Draw(c,kTRUE,format,1);
   

   //
  // W+ MET plot
  //
  RooPlot *wepframe = pfmet.frame(Bins(NBINS));    
  wepframe->GetYaxis()->SetNdivisions(505);
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetp.plotOn(wepframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetp.plotOn(wepframe,LineColor(linecolorW));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp)),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp)),LineColor(linecolorEWK));
  
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfWmp)),LineColor(linecolorW),LineStyle(2));
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
  CPlot plotMetp("fitmetp",wepframe,"","",ylabel);
  plotMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetp.GetLegend()->AddEntry(hDummyW,"Z#rightarrow#mu^{+}#mu^{-}","F");
  plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  //plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetp.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotMetp.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
//  plotMetp.SetYRange(0.1,1.1*(hDataMetp->GetMaximum()));
  //plotMetp.SetYRange(0.1,5000);
  plotMetp.Draw(c,kFALSE,format,1);


  TString titleX="#slash{E}_{T} [GeV]";
  if(doRecoilplot) titleX="Recoil  [GeV]";
  CPlot plotMetpDiff("fitmetp","",titleX.Data(),"#frac{Data-Pred}{Data}");
  plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
  plotMetpDiff.SetYRange(-0.2,0.2);
  plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotMetpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotMetpDiff.Draw(c,kTRUE,format,2);
  
  plotMetp.SetName("fitmetplog");
  plotMetp.SetLogy();
  plotMetp.SetYRange(1e-4*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
  plotMetp.Draw(c,kTRUE,format,1);
  
  
  RooPlot *weprframe = pfmet.frame(Bins(NBINS));    
  weprframe->GetYaxis()->SetNdivisions(505);
  dataRecoilp.plotOn(weprframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfRecoilp.plotOn(weprframe,FillColor(fillcolorW),DrawOption("F"));
  pdfRecoilp.plotOn(weprframe,LineColor(linecolorW));
  pdfRecoilp.plotOn(weprframe,Components(RooArgSet(pdfEWKpRecoil)),FillColor(fillcolorEWK),DrawOption("F"));
  pdfRecoilp.plotOn(weprframe,Components(RooArgSet(pdfEWKpRecoil)),LineColor(linecolorEWK));
  
  pdfRecoilp.plotOn(weprframe,Components(RooArgSet(pdfWmpRecoil)),LineColor(linecolorW),LineStyle(2));
  dataRecoilp.plotOn(weprframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataRecoilp->GetBinWidth(1));
  CPlot plotRecoilp("fitrecoilp",weprframe,"","",ylabel);
  plotRecoilp.SetLegend(0.68,0.57,0.93,0.77);
  plotRecoilp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotRecoilp.GetLegend()->AddEntry(hDummyW,"Z#rightarrow#mu^{+}#mu^{-}","F");
  plotRecoilp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  //plotRecoilp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotRecoilp.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotRecoilp.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
//  plotRecoilp.SetYRange(0.1,1.1*(hDataRecoilp->GetMaximum()));
  //plotRecoilp.SetYRange(0.1,5000);
  plotRecoilp.Draw(c,kFALSE,format,1);


  TString titleX2="Recoil [GeV]";
  if(doRecoilplot) titleX="Recoil  [GeV]";
  CPlot plotRecoilpDiff("fitrecoilp","",titleX.Data(),"#frac{Data-Pred}{Data}");
  plotRecoilpDiff.AddHist1D(hRecoilpDiff,"EX0",ratioColor);
  plotRecoilpDiff.SetYRange(-0.2,0.2);
  plotRecoilpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotRecoilpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotRecoilpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotRecoilpDiff.Draw(c,kTRUE,format,2);
  
  plotRecoilp.SetName("fitrecoilplog");
  plotRecoilp.SetLogy();
  plotRecoilp.SetYRange(1e-4*(hDataRecoilp->GetMaximum()),10*(hDataRecoilp->GetMaximum()));
  plotRecoilp.Draw(c,kTRUE,format,1);
  
  
  
  RooPlot *wepzframe = pfmet.frame(Bins(NBINS));    
  wepzframe->GetYaxis()->SetNdivisions(505);
  dataZpt.plotOn(wepzframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfZpt.plotOn(wepzframe,FillColor(fillcolorW),DrawOption("F"));
  pdfZpt.plotOn(wepzframe,LineColor(linecolorW));
  pdfZpt.plotOn(wepzframe,Components(RooArgSet(pdfEWKZpt)),FillColor(fillcolorEWK),DrawOption("F"));
  pdfZpt.plotOn(wepzframe,Components(RooArgSet(pdfEWKZpt)),LineColor(linecolorEWK));
  
  pdfZpt.plotOn(wepzframe,Components(RooArgSet(pdfWmZpt)),LineColor(linecolorW),LineStyle(2));
  dataZpt.plotOn(wepzframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataZpt->GetBinWidth(1));
  CPlot plotZpt("fitzpt",wepzframe,"","",ylabel);
  plotZpt.SetLegend(0.68,0.57,0.93,0.77);
  plotZpt.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotZpt.GetLegend()->AddEntry(hDummyW,"Z#rightarrow#mu^{+}#mu^{-}","F");
  plotZpt.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  //plotRecoilp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotZpt.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotZpt.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
//  plotRecoilp.SetYRange(0.1,1.1*(hDataRecoilp->GetMaximum()));
  //plotRecoilp.SetYRange(0.1,5000);
  plotZpt.Draw(c,kFALSE,format,1);


  titleX2="Z pT [GeV]";
  // if(doRecoilplot) titleX="Recoil  [GeV]";
  CPlot plotZptDiff("fitzpt","",titleX.Data(),"#frac{Data-Pred}{Data}");
  plotZptDiff.AddHist1D(hZptDiff,"EX0",ratioColor);
  plotZptDiff.SetYRange(-0.2,0.2);
  plotZptDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotZptDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotZptDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotZptDiff.Draw(c,kTRUE,format,2);
  
  plotZpt.SetName("fitzptlog");
  plotZpt.SetLogy();
  plotZpt.SetYRange(1e-4*(hDataZpt->GetMaximum()),10*(hDataZpt->GetMaximum()));
  plotZpt.Draw(c,kTRUE,format,1);
  
  
  
  //
  // W- MET plot
  //
  RooPlot *wmmframe = pfmet.frame(Bins(NBINS)); 
  wmmframe->GetYaxis()->SetNdivisions(505);
  dataMetm.plotOn(wmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetm.plotOn(wmmframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetm.plotOn(wmmframe,LineColor(linecolorW));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfEWKm)),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfEWKm)),LineColor(linecolorEWK));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfWmm)),LineColor(linecolorW),LineStyle(2));
  dataMetm.plotOn(wmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
  CPlot plotMetm("fitmetm",wmmframe,"","",ylabel);
  plotMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetm.GetLegend()->AddEntry(hDummyW,"Z#rightarrow#mu#mu","F");
  plotMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetm.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotMetm.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
//  plotMetm.SetYRange(0.1,1.1*(hDataMetm->GetMaximum()));
//plotMetm.SetYRange(0.1,4100);
  plotMetm.Draw(c,kFALSE,format,1);

  CPlot plotMetmDiff("fitmetm","","#slash{E}_{T} [GeV]","#chi");
  //CPlot plotMetmDiff("fitmetm","","mT [GeV]","#chi");
  plotMetmDiff.AddHist1D(hMetmDiff,"EX0",ratioColor);
  plotMetmDiff.SetYRange(-8,8);
  plotMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetmDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotMetmDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotMetmDiff.Draw(c,kTRUE,format,2);
  
  plotMetm.SetName("fitmetmlog");
  plotMetm.SetLogy();
  plotMetm.SetYRange(1e-3*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
  plotMetm.Draw(c,kTRUE,format,1);
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  
  //
  // Write fit results
  //
  ofstream txtfile;
  char txtfname[100];    
  
  ios_base::fmtflags flags;
  
  Double_t chi2prob, chi2ndf;
  Double_t ksprob, ksprobpe;
  
  chi2prob = hDataMet->Chi2Test(hPdfMet,"PUW");
  chi2ndf  = hDataMet->Chi2Test(hPdfMet,"CHI2/NDFUW");
  ksprob   = hDataMet->KolmogorovTest(hPdfMet);
  ksprobpe = hDataMet->KolmogorovTest(hPdfMet,"DX");
  sprintf(txtfname,"%s/fitresWm.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << hDataMet->Integral() << endl;
  //  txtfile << "  Signal: " << nSig.getVal() << " +/- " << nSig.getPropagatedError(*fitRes) << endl;
  //txtfile << "     QCD: " << nQCD.getVal() << " +/- " << nQCD.getPropagatedError(*fitRes) << endl;
  //  txtfile << "   Other: " << nEWK.getVal() << " +/- " << nEWK.getPropagatedError(*fitRes) << endl;
  
  chi2prob = hDataMetp->Chi2Test(hPdfMetp,"PUW");
  chi2ndf  = hDataMetp->Chi2Test(hPdfMetp,"CHI2/NDFUW");
  ksprob   = hDataMetp->KolmogorovTest(hPdfMetp);
  ksprobpe = hDataMetp->KolmogorovTest(hPdfMetp,"DX");  
  // sprintf(txtfname,"%s/fitresWmp.txt",CPlot::sOutDir.Data());
  // txtfile.open(txtfname);
  // assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << hDataMetp->Integral() << endl;
  txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResp) << endl;
  //txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResp) << endl;
  txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResp) << endl;
  txtfile << endl; 
  
    txtfile << "signal: " << hWmunuMetp->Integral() << endl;
  txtfile << "ewk: " << hEWKMetp->Integral() << endl;
  txtfile << "sum: " << hWmunuMetp->Integral()+hEWKMetp->Integral() << endl;
  txtfile << "data: " << hDataMetp->Integral() << endl;
  
  txtfile.flags(flags);
    fitResp->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResp);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  
  txtfile << endl;
  txtfile.flags(flags);

  /*  
  fitRes->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitRes);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();
  */
  
  chi2prob = hDataMetp->Chi2Test(hPdfMetp,"PUW");
  chi2ndf  = hDataMetp->Chi2Test(hPdfMetp,"CHI2/NDFUW");
  ksprob   = hDataMetp->KolmogorovTest(hPdfMetp);
  ksprobpe = hDataMetp->KolmogorovTest(hPdfMetp,"DX");  
  sprintf(txtfname,"%s/fitresWmp.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << hDataMetp->Integral() << endl;
  txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResp) << endl;
  //txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResp) << endl;
  txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResp) << endl;
  txtfile << endl; 
  txtfile.flags(flags);
  
  fitResp->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResp);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();

  chi2prob = hDataMetm->Chi2Test(hPdfMetm,"PUW");
  chi2ndf  = hDataMetm->Chi2Test(hPdfMetm,"CHI2/NDFUW");
  ksprob   = hDataMetm->KolmogorovTest(hPdfMetm);
  ksprobpe = hDataMetm->KolmogorovTest(hPdfMetm,"DX");  
  sprintf(txtfname,"%s/fitresWmm.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << hDataMetm->Integral() << endl;
  //txtfile << "  Signal: " << nSigm.getVal() << " +/- " << nSigm.getPropagatedError(*fitResm) << endl;
  //txtfile << "     QCD: " << nQCDm.getVal() << " +/- " << nQCDm.getPropagatedError(*fitResm) << endl;
  //txtfile << "   Other: " << nEWKm.getVal() << " +/- " << nEWKm.getPropagatedError(*fitResm) << endl;
  //txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResm) << endl;
  //txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResm) << endl;
  //txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResm) << endl;
  txtfile << endl;
  txtfile.flags(flags);
  
  //fitResm->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  //printCorrelations(txtfile, fitResm);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();

  makeHTML(outputDir);
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  gBenchmark->Show("fitZm");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    Double_t diff = diff0/hData->GetBinContent(ibin);
    if(hData->GetBinContent(ibin) == 0) diff = 0;
    std::cout << "data " << hData->GetBinContent(ibin) << std::endl;
    std::cout << "fits " << hFit->GetBinContent(ibin) << std::endl;
//     Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((1.0/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))*(hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    if(hData->GetBinContent(ibin) == 0) err = 0;
    //Double_t err = sqrt(hData->GetBinContent(ibin));
    //if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    //if(err>0) hDiff->SetBinContent(ibin,diff/err);
    //else      hDiff->SetBinContent(ibin,0);
    std::cout << "err = " << err << std::endl;
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
//--------------------------------------------------------------------------------------------------
void printCorrelations(ostream& os, RooFitResult *res)
{
  ios_base::fmtflags flags = os.flags();
  const RooArgList parlist = res->floatParsFinal();
  
  os << "  Correlation Matrix" << endl;
  os << " --------------------" << endl;
  for(Int_t i=0; i<parlist.getSize(); i++) {
    for(Int_t j=0; j<parlist.getSize(); j++) 
      os << "  " << setw(7) << setprecision(4) << fixed << res->correlationMatrix()(i,j);    
    os << endl;
  }
  os.flags(flags);
}

//--------------------------------------------------------------------------------------------------
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
			   const Double_t ksprob, const Double_t ksprobpe)
{
  ios_base::fmtflags flags = os.flags();
  
  os << "  Chi2 Test" << endl;
  os << " -----------" << endl;
  os << "       prob = " << chi2prob << endl;
  os << "   chi2/ndf = " << chi2ndf << endl;
  os << endl;
  os << "  KS Test" << endl;
  os << " ---------" << endl;
  os << "   prob = " << ksprob << endl;
  os << "   prob = " << ksprobpe << " with 1000 pseudo-experiments" << endl;
  os << endl;
 
  os.flags(flags);
}
			   
//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/WmunuFitPlots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Wmunu</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmet.png\"><img src=\"fitmet.png\" alt=\"fitmet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetp.png\"><img src=\"fitmetp.png\" alt=\"fitmetp.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetm.png\"><img src=\"fitmetm.png\" alt=\"fitmetm.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetlog.png\"><img src=\"fitmetlog.png\" alt=\"fitmetlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetplog.png\"><img src=\"fitmetplog.png\" alt=\"fitmetplog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetmlog.png\"><img src=\"fitmetmlog.png\" alt=\"fitmetmlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimet.png\"><img src=\"fitantimet.png\" alt=\"fitantimet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetp.png\"><img src=\"fitantimetp.png\" alt=\"fitantimetp.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetm.png\"><img src=\"fitantimetm.png\" alt=\"fitantimetm.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetlog.png\"><img src=\"fitantimetlog.png\" alt=\"fitantimetlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetplog.png\"><img src=\"fitantimetplog.png\" alt=\"fitantimetplog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetmlog.png\"><img src=\"fitantimetmlog.png\" alt=\"fitantimetmlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();  
}
