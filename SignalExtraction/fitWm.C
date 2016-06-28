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

#include "BaconAna/DataFormats/interface/TGenParticle.hh"  

#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"              // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../Utils/RecoilCorrector_asym2.hh"
// #include "../Utils/RecoilCorrector_addJets.hh"
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections

// #include "ZBackgrounds.hh"

//helper class to handle rochester corrections
#include <rochcor2015r.h>
#include <muresolution_run2r.h>

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
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

// print correlations of fitted parameters
void printCorrelations(ostream& os, RooFitResult *res);

// print chi2 test and KS test results
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
               const Double_t ksprob, const Double_t ksprobpe);

// make webpage
void makeHTML(const TString outDir);


//=== MAIN MACRO ================================================================================================= 

void fitWm(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)'
           const Double_t lumi2, // lumi for the anti-isolation trigger
       const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
) {
  gBenchmark->Start("fitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  // MET histogram binning and range
  const Int_t    NBINS   = 75;
  const Double_t METMAX  = 150;
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  
  const Double_t mu_MASS = 0.1057;
  const int NTOYS = 100;

  TString pufname = "../Tools/pileup_weights_2015B.root";
  
//   TString ZmmShapes = "./ZmmShapes/Zmumu_pdfTemplates_BkgSub.root";

  // file format for output plots
  const TString format("png"); 


//   // for Puppi
  RecoilCorrector *recoilCorr = new  RecoilCorrector("../Recoil/ZmmMCPuppi_newBacon/fits_puppi.root","fcnPF"); // get tgraph from here? what i guess its mean so it doesn't matter?
  recoilCorr->loadRooWorkspacesData("../Recoil/ZmmDataPuppi_newBacon/");
  recoilCorr->loadRooWorkspacesMC("../Recoil/WmpMC_newBacon/");
  
  RecoilCorrector *recoilCorrm = new  RecoilCorrector("../Recoil/ZmmMCPuppi_newBacon/fits_puppi.root","fcnPF"); // get tgraph from here? what i guess its mean so it doesn't matter?
  recoilCorrm->loadRooWorkspacesData("../Recoil/ZmmDataPuppi_newBacon/");
  recoilCorrm->loadRooWorkspacesMC("../Recoil/WmmMC_newBacon/");
  
  
//     // for PF
//   RecoilCorrector *recoilCorr = new  RecoilCorrector("../Recoil/ZmmMCPF/fits_pf.root","fcnPF"); // get tgraph from here? what i guess its mean so it doesn't matter?
//   recoilCorr->loadRooWorkspacesData("../Recoil/ZmmDataPF/");
//   recoilCorr->loadRooWorkspacesMC("../Recoil/WmpMCPF/");
//   
//   RecoilCorrector *recoilCorrm = new  RecoilCorrector("../Recoil/ZmmMCPF/fits_pf.root","fcnPF"); // get tgraph from here? what i guess its mean so it doesn't matter?
//   recoilCorrm->loadRooWorkspacesData("../Recoil/ZmmDataPF/");
//   recoilCorrm->loadRooWorkspacesMC("../Recoil/WmmMCPF/");
  
  
  TFile *_rdWmp = new TFile("shapeDiff/WmpPF_relDiff.root");
  TFile *_rdWmm = new TFile("shapeDiff/WmmPF_relDiff.root");
  
  TH1D *hh_diffm = new TH1D("hh_diffm","hh_diffm",75,0,150);
  TH1D *hh_diffp = new TH1D("hh_diffp","hh_diffp",75,0,150);
  
  hh_diffm = (TH1D*)_rdWmm->Get("hh_diff");
  hh_diffp = (TH1D*)_rdWmp->Get("hh_diff");
   
  
 /*Fit to data
  Par  0                    p0 = -0.23882
  Par  1                    p1 = 0.0250329
  Par  2                    p2 = 0.000162729
  Par  3                    p3 = -5.98736e-05
  Par  4                    p4 = 1.89858e-06
  Par  5                    p5 = -2.48375e-08
  Par  6                    p6 = 1.46665e-10
  Par  7                    p7 = -3.22598e-13
 */
//     TF1 *fit = new TF1("fit","pol7");
//   fit->SetParameter(0,-0.23882);
//   fit->SetParameter(1,0.0250329);
//   fit->SetParameter(2, 0.000162729);
//   fit->SetParameter(3,-5.98736e-05);
//   fit->SetParameter(4,1.89858e-06);
//   fit->SetParameter(5,-2.48375e-08);
//   fit->SetParameter(6, 1.46665e-10);
//   fit->SetParameter(7,-3.22598e-13);
 
/* fit to nominal MC model  
Par  0                    p0 = -0.249167                         
 Par  1                    p1 = 0.0363433                         
 Par  2                    p2 = -0.00129434                       
 Par  3                    p3 = 8.57774e-06                       
 Par  4                    p4 = 3.60219e-07                       
 Par  5                    p5 = -7.52258e-09                      
 Par  6                    p6 = 5.3105e-11                        
 Par  7                    p7 = -1.28952e-13  */

/*calculate it correctly...
 Par  0                    p0 = -0.22887
 Par  1                    p1 = 0.0348117
 Par  2                    p2 = -0.00139386
 Par  3                    p3 = 2.03053e-05
 Par  4                    p4 = -2.71221e-08
 Par  5                    p5 = -1.69981e-09
 Par  6                    p6 = 1.19046e-11
 Par  7                    p7 = -1.7821e-14 */

 
//   TF1 *fit = new TF1("fit","pol7");
//   fit->SetParameter(0,-0.22887);
//   fit->SetParameter(1,0.0348117);
//   fit->SetParameter(2, -0.00139386);
//   fit->SetParameter(3,2.03053e-05);
//   fit->SetParameter(4,-2.71221e-08);
//   fit->SetParameter(5,-1.69981e-09);
//   fit->SetParameter(6, 1.19046e-11 );
//   fit->SetParameter(7,-1.7821e-14);
//   fit->Draw("same");
  //   return;
  // TH1D *hh_diff = (TH1D*) hh_diff->Clone("z_corr2");
//   for(int i = 0; i < 75; i++){hh_diff->SetBinContent(i,fit->Eval(2*i+1));}
  //   hh_diff->Draw();
  

  //
  // input ntuple file names
  //
  enum { eData, eWmunu, eEWK, eQCD, eAntiData, eAntiWmunu, eAntiEWK };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  
  
  fnamev.push_back("/afs/cern.ch/user/s/sabrandt/work/public/SM/newBacon/Wmunu/ntuples/data_select.root");     typev.push_back(eData);
  fnamev.push_back("/afs/cern.ch/user/s/sabrandt/work/public/SM/newBacon/Wmunu/ntuples/wm_select.raw.root");   typev.push_back(eWmunu);
  fnamev.push_back("/afs/cern.ch/user/s/sabrandt/work/public/SM/newBacon/Wmunu/ntuples/ewk_select.raw.root");  typev.push_back(eEWK);
  fnamev.push_back("/afs/cern.ch/user/s/sabrandt/work/public/SM/newBacon/Wmunu/ntuples/top_select.raw.root");  typev.push_back(eEWK);
  
// // //   // The anti-selection stuff
//   fnamev.push_back("/afs/cern.ch/user/s/sabrandt/work/public/SM/newBacon/AntiWmunu/ntuples/data_select.root"); typev.push_back(eData);
//   fnamev.push_back("/afs/cern.ch/user/s/sabrandt/work/public/SM/newBacon/AntiWmunu/ntuples/wm_select.root");   typev.push_back(eWmunu);
//   fnamev.push_back("/afs/cern.ch/user/s/sabrandt/work/public/SM/newBacon/AntiWmunu/ntuples/ewk_select.root");  typev.push_back(eEWK);
//   fnamev.push_back("/afs/cern.ch/user/s/sabrandt/work/public/SM/newBacon/AntiWmunu/ntuples/top_select.root");  typev.push_back(eEWK);


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
  
  

  
/*  
  TH1D *hEWKMet    = new TH1D("hEWKMet", "",  NBINS,0,METMAX); hEWKMet->Sumw2();
  TH1D *hEWKMetp   = new TH1D("hEWKMetp", "", NBINS,0,METMAX); hEWKMetp->Sumw2();
  TH1D *hEWKMetm   = new TH1D("hEWKMetm", "", NBINS,0,METMAX); hEWKMetm->Sumw2();*/
  
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  TLorentzVector *lep=0, *genV=0;
  Float_t pfChIso, pfGamIso, pfNeuIso;
    
  //Setting up rochester corrections
//   vector<rochcor2015> vRocToys;
//   for (int i=0 ; i<NTOYS; ++i) vRocToys.push_back(rochcor2015(1234+i*1000));
  
  rochcor2015 *rmcor = new rochcor2015();
  
  TFile *infile=0;
  TTree *intree=0;

  //
  // Loop over files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);    assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);

    intree->SetBranchAddress("runNum",   &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    intree->SetBranchAddress("npv",      &npv);       // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);       // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)   
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
//     intree->SetBranchAddress("mvaMet",      &met);       // MET
//     intree->SetBranchAddress("mvaMetPhi",   &metPhi);    // phi(MET)
//     intree->SetBranchAddress("met",      &met);       // MET
//     intree->SetBranchAddress("metPhi",   &metPhi);    // phi(MET)
    intree->SetBranchAddress("puppiMet",      &met);       // MET
    intree->SetBranchAddress("puppiMetPhi",   &metPhi);    // phi(MET)
    intree->SetBranchAddress("sumEt",    &sumEt);     // Sum ET
    intree->SetBranchAddress("mt",       &mt);        // transverse mass
    intree->SetBranchAddress("puppiU1",       &u1);        // parallel component of recoil
    intree->SetBranchAddress("puppiU2",       &u2);        // perpendicular component of recoil
//     intree->SetBranchAddress("u1",       &u1);        // parallel component of recoil
//     intree->SetBranchAddress("u2",       &u2);        // perpendicular component of recoil
    intree->SetBranchAddress("q",        &q);         // lepton charge
    intree->SetBranchAddress("lep",      &lep);       // lepton 4-vector
    intree->SetBranchAddress("genV",     &genV);       // lepton 4-vector
    intree->SetBranchAddress("pfChIso",  &pfChIso);
    intree->SetBranchAddress("pfGamIso", &pfGamIso);
    intree->SetBranchAddress("pfNeuIso", &pfNeuIso);
  
    Double_t mt=-999;

    //
    // loop over events
    //
    std::cout << "Number of Events = " << intree->GetEntries() << std::endl;
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
//       std::cout << ientry << std::endl;
      intree->GetEntry(ientry);
      if(ientry%100000==0) std::cout << "On Entry.... " << ientry << std::endl;

      double pU1         = 0;  //--
      double pU2         = 0;  //--
      
//       if(fabs(lep->Eta()) > 1.0) continue;

    
      mt     = sqrt( 2.0 * (lep->Pt()) * (met) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metPhi))) );
      
      
      if(typev[ifile]==eData) {
          
        // Apply the Rochester Corrections to data
        TLorentzVector mu1;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        float qter1=1.0;
        rmcor->momcor_data(mu1,q,0,qter1);
        if(mu1.Pt()        < PT_CUT)  continue;
        hDataMet->Fill(met);
        if(q>0) { hDataMetp->Fill(met); } 
        else    { hDataMetm->Fill(met); }
      } else if(typev[ifile]==eAntiData) {
        // Apply the Rochester Corrections to data (anti-isolation)
        TLorentzVector mu1;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        float qter1=1.0;
        rmcor->momcor_data(mu1,q,0,qter1);
        if(mu1.Pt()        < PT_CUT)  continue;
        
        hAntiDataMet->Fill(met);
        if(q>0) { hAntiDataMetp->Fill(met); } 
        else    { hAntiDataMetm->Fill(met); }   
      } else {
        Double_t weight = 1;
        weight *= scale1fb*lumi;
        
        Double_t weight2 = 1;
        weight2 *= scale1fb*lumi2;
        
        // Do some Rochester corrections for MC
        TLorentzVector mu1;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        float qter1=1.0;

        rmcor->momcor_mc(mu1,q,0,qter1);
        Double_t lepPt = mu1.Pt();
        
        if(typev[ifile]==eWmunu) {
          Double_t corrMet=met, corrMetPhi=metPhi;
          
          if(lepPt        > PT_CUT) {
              //recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),nsigma,q); 
            hWmunuMet->Fill(corrMet,weight);
            if(q>0) {
              pU1 = 0; pU2 = 0; 
//           recoilCorr->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
//               if(genV->Eta() < 1.0) recoilCorr1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
//               if(genV->Eta() > 1.0) continue;
//               if(genV->Eta() > 1.0) recoilCorr2->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
              recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
//                recoilCorr->CorrectFromToys(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
//           recoilCorr->CorrectType0(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
              hWmunuMetp->Fill(corrMet,weight); 
              corrMet=met, corrMetPhi=metPhi;
            } else { 
              pU1 = 0; pU2 = 0; 
//           recoilCorrm->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
// //               if(genV->Eta() > 1.0) continue;
//               if(genV->Eta() > 1.0) recoilCorrm2->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
              recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
//               recoilCorrm->CorrectFromToys(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
//           recoilCorrm->CorrectType0(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
              hWmunuMetm->Fill(corrMet,weight); 
              corrMet=met, corrMetPhi=metPhi;
            }
            //recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),1,q);
            hWmunuMet_RecoilUp->Fill(corrMet,weight);
            if(q>0) 
            {
            pU1 = 0; pU2 = 0; 
//           recoilCorr->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,4,4);
//           recoilCorr->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,4,4);
            hWmunuMetp_RecoilUp->Fill(corrMet,weight); 
            corrMet=met, corrMetPhi=metPhi;
            }
            else    
            { 
            pU1 = 0; pU2 = 0; 
//           recoilCorrm->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,4,4);
//           recoilCorrm->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,4,4);
            hWmunuMetm_RecoilUp->Fill(corrMet,weight);
            corrMet=met, corrMetPhi=metPhi;
            }
            //recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),-1,q);
            hWmunuMet_RecoilDown->Fill(corrMet,weight);
            if(q>0) 
            {
            pU1 = 0; pU2 = 0; 
  //           recoilCorr->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,-4,-4);
  //           recoilCorr->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,-4,-4);
            hWmunuMetp_RecoilDown->Fill(corrMet,weight);
            corrMet=met, corrMetPhi=metPhi;
            } 
            else    
            { 
            pU1 = 0; pU2 = 0; 
  //           recoilCorrm->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,-4,-4);
  //           recoilCorrm->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,-4,-4);
            hWmunuMetm_RecoilDown->Fill(corrMet,weight);
            corrMet=met, corrMetPhi=metPhi;
            }
          }
          Double_t lepPtup = (gRandom->Gaus((lep->Pt())*getEleScaleCorr(lep->Eta(),1),getEleResCorr(lep->Eta(),1)));  // (!) uncomment to apply scale/res corrections to MC
          if(lepPtup        > PT_CUT) {
            corrMet=met, corrMetPhi=metPhi;
            //recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPtup,lep->Phi(),0,q);
            hWmunuMet_ScaleUp->Fill(corrMet,weight);
            if(q>0){
            pU1 = 0; pU2 = 0; 
  //           recoilCorr->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,lepPtup,lep->Phi(),pU1,pU2,0);
  //           recoilCorr->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,lepPtup,lep->Phi(),pU1,pU2,0);
            hWmunuMetp_ScaleUp->Fill(corrMet,weight); 
            corrMet=met, corrMetPhi=metPhi;
            } else {
              pU1 = 0; pU2 = 0; 
    //           recoilCorrm->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,lepPtup,lep->Phi(),pU1,pU2,0);
    //           recoilCorrm->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,lepPtup,lep->Phi(),pU1,pU2,0);
              hWmunuMetm_ScaleUp->Fill(corrMet,weight);
              corrMet=met, corrMetPhi=metPhi;
            }
          }
          Double_t lepPtdown = (gRandom->Gaus((lep->Pt())*getEleScaleCorr(lep->Eta(),-1),getEleResCorr(lep->Eta(),-1)));  // (!) uncomment to apply scale/res corrections to MC
          if(lepPtdown        > PT_CUT) {
            corrMet=met, corrMetPhi=metPhi;
            //recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPtdown,lep->Phi(),0,q);
            hWmunuMet_ScaleDown->Fill(corrMet,weight);
            if(q>0) {
              pU1 = 0; pU2 = 0; 
  //           recoilCorr->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,lepPtup,lep->Phi(),pU1,pU2,0);
  //           recoilCorr->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,lepPtup,lep->Phi(),pU1,pU2,0);
              hWmunuMetp_ScaleDown->Fill(corrMet,weight);
              corrMet=met, corrMetPhi=metPhi;
            } else { 
              pU1 = 0; pU2 = 0; 
  //           recoilCorrm->CorrectType2FunG(corrMet,corrMetPhi,genVPt,genVPhi,lepPtdown,lep->Phi(),pU1,pU2,0);
  //           recoilCorrm->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,lepPtdown,lep->Phi(),pU1,pU2,0);
              hWmunuMetm_ScaleDown->Fill(corrMet,weight); 
            }
          }
        }
        if(typev[ifile]==eAntiWmunu) {
          Double_t corrMet=met, corrMetPhi=metPhi;
          // apply recoil corrections to W MC
          // old leptoon scale thing
//           Double_t lepPt = (gRandom->Gaus((lep->Pt())*getMuScaleCorr(lep->Eta(),0),getMuResCorr(lep->Eta(),0)));  // (!) uncomment to apply scale/res corrections to MC

          hAntiWmunuMet->Fill(corrMet,weight2);
          if(q>0) {               
            pU1 = 0; pU2 = 0; 
//             recoilCorr->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
            recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
            hAntiWmunuMetp->Fill(corrMet,weight2); 
            corrMet = met; corrMetPhi = metPhi;
          } 
          else { 
            pU1 = 0; pU2 = 0; 
//             recoilCorrm->CorrectType2FromGraph(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
            recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
            hAntiWmunuMetm->Fill(corrMet,weight2);
            corrMet = met; corrMetPhi = metPhi; 
          }
        }
        if(typev[ifile]==eEWK) {
          if(lep->Pt()        < PT_CUT)  continue;
          hEWKMet->Fill(met,weight);
          if(q>0) { hEWKMetp->Fill(met,weight); }
          else    { hEWKMetm->Fill(met,weight); }
        }
        if(typev[ifile]==eAntiEWK) {
          hAntiEWKMet->Fill(met,weight);
          if(q>0) { hAntiEWKMetp->Fill(met,weight2); }
          else    { hAntiEWKMetm->Fill(met,weight2); }
        }
      }
    }
  }
  delete infile;
  infile=0, intree=0;   
 
  // Here rescale the up/down whatever
//   Calculate the shapes for W+
  TH1D *corrP = (TH1D*) hWmunuMetp->Clone("up");
  *corrP = (*hh_diffp)*(*hWmunuMetp);
  hWmunuMetp_RecoilUp->Add(corrP,hWmunuMetp,-1);
  hWmunuMetp_RecoilDown->Add(corrP,hWmunuMetp,1);
  hWmunuMetp_ScaleUp->Add(corrP,hWmunuMetp,-1);
  hWmunuMetp_ScaleDown->Add(corrP,hWmunuMetp,1);
  // Calculate the shapes for W-
  TH1D *corrM = (TH1D*) hWmunuMetm->Clone("up");
  *corrM = (*hh_diffm)*(*hWmunuMetm);
  hWmunuMetm_RecoilUp->Add(corrM,hWmunuMetm,-1);
  hWmunuMetm_RecoilDown->Add(corrM,hWmunuMetm,1);
  hWmunuMetm_ScaleUp->Add(corrM,hWmunuMetm,-1);
  hWmunuMetm_ScaleDown->Add(corrM,hWmunuMetm,1);
  
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  //
  RooRealVar nSig("nSig","nSig",0.7*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar nQCD("nQCD","nQCD",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar cewk("cewk","cewk",0.1,0,5) ;
  cewk.setVal(hEWKMet->Integral()/hWmunuMet->Integral());
  cewk.setConstant(kTRUE);
  RooFormulaVar nEWK("nEWK","nEWK","cewk*nSig",RooArgList(nSig,cewk));
  RooRealVar nAntiSig("nAntiSig","nAntiSig",hAntiWmunuMet->Integral()*0.9,0,hAntiDataMet->Integral());
  RooRealVar nAntiQCD("nAntiQCD","nAntiQCD",0.9*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar dewk("dewk","dewk",0.1,0,5) ;
  dewk.setVal(hAntiEWKMet->Integral()/hAntiWmunuMet->Integral());
  dewk.setConstant(kTRUE);
//   nAntiSig.setConstant(kTRUE);
  RooFormulaVar nAntiEWK("nAntiEWK","nAntiEWK","dewk*nAntiSig",RooArgList(nAntiSig,dewk));
  
  RooRealVar nSigp("nSigp","nSigp",0.7*(hDataMetp->Integral()),0,hDataMetp->Integral());
//   RooRealVar nSigp("nSigp","nSigp",90000,0,hDataMetp->Integral());
  //RooRealVar nQCDp("nQCDp","nQCDp",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nQCDp("nQCDp","nQCDp",hDataMetp->Integral()*0.3,0,hDataMetp->Integral());
  RooRealVar cewkp("cewkp","cewkp",0.1,0,5) ;
  cewkp.setVal(hEWKMetp->Integral()/hWmunuMetp->Integral());
  cewkp.setConstant(kTRUE);
  RooFormulaVar nEWKp("nEWKp","nEWKp","cewkp*nSigp",RooArgList(nSigp,cewkp));
  RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWmunuMetp->Integral()*0.9,0,hAntiDataMetp->Integral());
  RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar dewkp("dewkp","dewkp",0.1,0,5) ;
  dewkp.setVal(hAntiEWKMetp->Integral()/hAntiWmunuMetp->Integral());
  dewkp.setConstant(kTRUE);
//   nAntiSigp.setConstant(kTRUE);
  RooFormulaVar nAntiEWKp("nAntiEWKp","nAntiEWKp","dewkp*nAntiSigp",RooArgList(nAntiSigp,dewkp));
  
  RooRealVar nSigm("nSigm","nSigm",0.7*(hDataMetm->Integral()),0,hDataMetm->Integral());
//   RooRealVar nSigm("nSigm","nSigm",75000,0,hDataMetm->Integral());
  //RooRealVar nQCDm("nQCDm","nQCDm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  RooRealVar nQCDm("nQCDm","nQCDm",hDataMetm->Integral()*0.3,0,hDataMetm->Integral());
  RooRealVar cewkm("cewkm","cewkm",0.1,0,5) ;
  cewkm.setVal(hEWKMetm->Integral()/hWmunuMetm->Integral());
  cewkm.setConstant(kTRUE);
  RooFormulaVar nEWKm("nEWKm","nEWKm","cewkm*nSigm",RooArgList(nSigm,cewkm));  
  RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",hAntiWmunuMetm->Integral()*0.9,0,hAntiDataMetm->Integral());
  RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  RooRealVar dewkm("dewkm","dewkm",0.1,0,5) ;
  dewkm.setVal(hAntiEWKMetm->Integral()/hAntiWmunuMetm->Integral());
  dewkm.setConstant(kTRUE);
//   nAntiSigm.setConstant(kTRUE);
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
  
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", pfmet,ewkMet, 1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",pfmet,ewkMetp,1); 
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",pfmet,ewkMetm,1); 
  
//   // test using the reversed dEta, dPhi cuts as background
//   RooDataHist qcdMet ("qcdMET", "qcdMET", RooArgSet(pfmet),hQCDMet);  RooHistPdf pdfQCD ("qcd", "qcd", pfmet,qcdMet, 1);
//   RooDataHist qcdMetp("qcdMetp","qcdMetp",RooArgSet(pfmet),hQCDMetp); RooHistPdf pdfQCDp("qcdp","qcdp",pfmet,qcdMetp,1); 
//   RooDataHist qcdMetm("qcdMetm","qcdMetm",RooArgSet(pfmet),hQCDMetm); RooHistPdf pdfQCDm("qcdm","qcdm",pfmet,qcdMetm,1); 
  
  // QCD Pdfs
  

  
  //CExponential qcd(pfmet,kTRUE);
  //CExponential qcdp(pfmet,kTRUE);
  //CExponential qcdm(pfmet,kTRUE);
  // comment back in for qcd functional form
  CPepeModel1 qcd("qcd",pfmet);
  CPepeModel1 qcdp("qcdp",pfmet);
  CPepeModel1 qcdm("qcdm",pfmet);
  //CPepeModel2 qcdm("qcdm",pfmet);
  //CPepeModel1 qcdm("qcdm",pfmet);
//  qcdp.a1->setConstant(kTRUE);
//   qcdp.a1->setVal(0.24);
//    qcdm.a1->setVal(0.24);

//   RooRealVar a1ConstMeanP("a1ConstMeanP","a1ConstMeanP",0.18);
//   RooRealVar a1ConstSigmaP("a1ConstSigmaP","a1ConstSigmaP",0.003);
// //   RooRealVar a1ConstMeanM("a1ConstMeanM","a1ConstMeanM",0.215);
// //   RooRealVar a1ConstSigmaM("a1ConstSigmaM","a1ConstSigmaM",0.007);
// 
//   RooRealVar sigConstMeanP("sigConstMeanP","sigConstMeanP",17.6);
//   RooRealVar sigConstSigmaP("sigConstSigmaP","sigConstSigmaP",0.2);
//   RooRealVar sigConstMeanM("sigConstMeanM","sigConstMeanM",16.9);
//   RooRealVar sigConstSigmaM("sigConstSigmaM","sigConstSigmaM",0.5);
//   
// // RooRealVar f("f","f",0.5,0.,1.) ;
//   RooGaussian fconsta1p("fconsta1p","fconsta1p",*(qcdp.a1),a1ConstMeanP,a1ConstSigmaP);
// //   RooGaussian fconsta1m("fconsta1m","fconsta1m",*(qcdm.a1),a1ConstMeanM,a1ConstSigmaM);
//   
//   RooGaussian fconstsigp("fconstsigp","fconstsigp",*(qcdp.sigma),sigConstMeanP,sigConstSigmaP);
//   RooGaussian fconstsigm("fconstsigm","fconstsigm",*(qcdm.sigma),sigConstMeanM,sigConstSigmaM);
//  
//   RooProdPdf qcdpc("qcdpc","qcdpc",RooArgSet(*(qcdp.model),fconsta1p,fconstsigp));
//   RooProdPdf qcdmc("qcdmc","qcdmc",RooArgSet(*(qcdm.model),fconstsigm));
//   RooProdPdf qcdmc("qcdmc","qcdmc",RooArgSet(qcdp,fconstraint));
 
  // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));
  
//   // constrained background?
//   RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
//   RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,qcdpc),RooArgList(nSigp,nEWKp,nQCDp));
//   RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,qcdmc),RooArgList(nSigm,nEWKm,nQCDm));
  
  
    
  
  // Anti-Signal PDFs
  RooDataHist awmunuMet ("awmunuMET", "awmunuMET", RooArgSet(pfmet),hAntiWmunuMet);  RooHistPdf apdfWm ("awm", "awm", pfmet,awmunuMet, 1);
  RooDataHist awmunuMetp("awmunuMETp","awmunuMETp",RooArgSet(pfmet),hAntiWmunuMetp); RooHistPdf apdfWmp("awmp","awmp",pfmet,awmunuMetp,1);
  RooDataHist awmunuMetm("awmunuMETm","awmunuMETm",RooArgSet(pfmet),hAntiWmunuMetm); RooHistPdf apdfWmm("awmm","awmm",pfmet,awmunuMetm,1); 
  
  // Anti-EWK+top PDFs
  RooDataHist aewkMet ("aewkMET", "aewkMET", RooArgSet(pfmet),hAntiEWKMet);  RooHistPdf apdfEWK ("aewk", "aewk", pfmet,aewkMet, 1);
  RooDataHist aewkMetp("aewkMETp","aewkMETp",RooArgSet(pfmet),hAntiEWKMetp); RooHistPdf apdfEWKp("aewkp","aewkp",pfmet,aewkMetp,1); 
  RooDataHist aewkMetm("aewkMETm","aewkMETm",RooArgSet(pfmet),hAntiEWKMetm); RooHistPdf apdfEWKm("aewkm","aewkm",pfmet,aewkMetm,1); 
  
//   // Anti-QCD Pdfs
//   CPepeModel1 aqcd("aqcd",pfmet, qcd.a1);
//   CPepeModel1 aqcdp("aqcdp",pfmet, qcdp.a1);
//   CPepeModel1 aqcdm("aqcdm",pfmet, qcdm.a1);

  CPepeModel1 aqcd("aqcd",pfmet);
  CPepeModel1 aqcdp("aqcdp",pfmet);
  CPepeModel1 aqcdm("aqcdm",pfmet);
  
  // Anti-selection PDFs
  RooAddPdf apdfMet ("apdfMet", "apdfMet", RooArgList(apdfWm,apdfEWK,*(aqcd.model)),   RooArgList(nAntiSig,nAntiEWK,nAntiQCD));  
  RooAddPdf apdfMetp("apdfMetp","apdfMetp",RooArgList(apdfWmp,apdfEWKp,*(aqcdp.model)),RooArgList(nAntiSigp,nAntiEWKp,nAntiQCDp));
  RooAddPdf apdfMetm("apdfMetm","apdfMetm",RooArgList(apdfWmm,apdfEWKm,*(aqcdm.model)),RooArgList(nAntiSigm,nAntiEWKm,nAntiQCDm));
  
  // PDF for simultaneous fit
  RooCategory rooCat("rooCat","rooCat");
  rooCat.defineType("Selectp");
  rooCat.defineType("Selectm");
  
  RooSimultaneous pdfTotal("pdfTotal","pdfTotal",rooCat);
  pdfTotal.addPdf(pdfMet, "Select");
  //pdfTotal.addPdf(apdfMet,"Anti");
  
  RooSimultaneous pdfTotalp("pdfTotalp","pdfTotalp",rooCat);
  pdfTotalp.addPdf(pdfMetp, "Selectp");
//   pdfTotalp.addPdf(pdfMetm,"Selectm");
  //pdfTotalp.addPdf(apdfMetp,"Anti");
  
  RooSimultaneous pdfTotalm("pdfTotalm","pdfTotalm",rooCat);
  pdfTotalm.addPdf(pdfMetm, "Select");
  //pdfTotalm.addPdf(apdfMetm,"Anti");

  
  //
  // Perform fits
  //

  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet), hDataMet);
  RooDataHist dataMetp("dataMetp", "dataMetp", RooArgSet(pfmet), hDataMetp);
  RooDataHist dataMetm("dataMetm", "dataMetm", RooArgSet(pfmet), hDataMetm);

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

//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",1000000,0,hDataMet->Integral()*0.5);
  RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",1000000,0,hDataMet->Integral()*0.5);

  RooWorkspace combine_workspace("combine_workspace");
  combine_workspace.import(dataMet);
  combine_workspace.import(dataMetp);
  combine_workspace.import(dataMetm);
  combine_workspace.import(pepe1Pdf_qcdp_norm);
  combine_workspace.import(pepe1Pdf_qcdm_norm);

  combine_workspace.import(pdfWm);
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
  combine_workspace.import(pdfEWKm);
  combine_workspace.import(*(qcd.model));
  //combine_workspace.import(qcdpn);
  combine_workspace.import(*(qcdp.model));
  combine_workspace.import(*(qcdm.model));

  combine_workspace.writeToFile("Wmunu_pdfTemplates.root");

  RooDataHist antiMet("antiMet", "antiMet", RooArgSet(pfmet), hAntiDataMet);
  RooDataHist dataTotal("dataTotal","dataTotal", RooArgList(pfmet), Index(rooCat),
            Import("Select", dataMet),
            Import("Anti",   antiMet));
  RooFitResult *fitRes = 0;//pdfMet.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE));//dataTotal.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE));
  
 
  RooDataHist antiMetp("antiMetp", "antiMetp", RooArgSet(pfmet), hAntiDataMetp);
  RooDataHist dataTotalp("dataTotalp","dataTotalp", RooArgList(pfmet), Index(rooCat),
             Import("Select", dataMetp),
             Import("Anti",   antiMetp));
             
  RooFitResult *fitResp = pdfMetp.fitTo(dataMetp,Extended(),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResAntip = apdfMetp.fitTo(antiMetp,Extended(),Minos(kTRUE),Save(kTRUE));
// RooFitResult *fitResp = pdfMetp.fitTo(dataTotalp,Extended(),Minos(kTRUE),Save(kTRUE));
  
  RooDataHist antiMetm("antiMetm", "antiMetm", RooArgSet(pfmet), hAntiDataMetm);
//   RooDataHist dataMetp2("dataMetp2", "dataMetp2", RooArgSet(pfmet), hDataMetp);
  RooDataHist dataTotalm("dataTotalm","dataTotalm", RooArgList(pfmet), Index(rooCat),
             Import("Select", dataMetm),
             Import("Anti", antiMetm));
  RooFitResult *fitResm = pdfMetm.fitTo(dataMetm,Extended(),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResm = pdfMetm.fitTo(dataTotalm,Extended(),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResAntim = apdfMetm.fitTo(antiMetm,Extended(),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResm = pdfTotalp.fitTo(dataTotalm,Extended(),Minos(kTRUE),Save(kTRUE));
  
  //
  // Use histogram version of fitted PDFs to make ratio plots
  // (Will also use PDF histograms later for Chi^2 and KS tests)
  //
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", pfmet));
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  hMetDiff->SetMarkerStyle(kFullCircle);
  hMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", pfmet));
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet));
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle); 
  hMetmDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfAntiMet = (TH1D*)(apdfMet.createHistogram("hPdfAntiMet", pfmet));
  hPdfAntiMet->Scale((nAntiSig.getVal()+nAntiEWK.getVal()+nAntiQCD.getVal())/hPdfAntiMet->Integral());
  TH1D *hAntiMetDiff = makeDiffHist(hAntiDataMet,hPdfAntiMet,"hAntiMetDiff");
  hAntiMetDiff->SetMarkerStyle(kFullCircle);
  hAntiMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfAntiMetp = (TH1D*)(apdfMetp.createHistogram("hPdfAntiMetp", pfmet));
  hPdfAntiMetp->Scale((nAntiSigp.getVal()+nAntiEWKp.getVal()+nAntiQCDp.getVal())/hPdfAntiMetp->Integral());
  TH1D *hAntiMetpDiff = makeDiffHist(hAntiDataMetp,hPdfAntiMetp,"hAntiMetpDiff");
  hAntiMetpDiff->SetMarkerStyle(kFullCircle);
  hAntiMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfAntiMetm = (TH1D*)(apdfMetm.createHistogram("hPdfAntiMetm", pfmet));
  hPdfAntiMetm->Scale((nAntiSigm.getVal()+nAntiEWKm.getVal()+nAntiQCDm.getVal())/hPdfAntiMetm->Integral());
  TH1D *hAntiMetmDiff = makeDiffHist(hAntiDataMetm,hPdfAntiMetm,"hAntiMetmDiff");
  hAntiMetmDiff->SetMarkerStyle(kFullCircle); 
  hAntiMetmDiff->SetMarkerSize(0.9);
   
  
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
  
  TH1D *hDummyQCD = new TH1D("hDummyQCD","",0,0,10);
  hDummyQCD->SetLineColor(linecolorQCD);
  hDummyQCD->SetFillColor(fillcolorQCD);
  hDummyQCD->SetFillStyle(1001);
   
  //
  // W MET plot
  //
  RooPlot *wmframe = pfmet.frame(Bins(NBINS)); 
  wmframe->GetYaxis()->SetNdivisions(505);
  dataMet.plotOn(wmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMet.plotOn(wmframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMet.plotOn(wmframe,LineColor(linecolorW));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfEWK,*(qcd.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfEWK,*(qcd.model))),LineColor(linecolorEWK));
  pdfMet.plotOn(wmframe,Components(RooArgSet(*(qcd.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMet.plotOn(wmframe,Components(RooArgSet(*(qcd.model))),LineColor(linecolorQCD));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfWm)),LineColor(linecolorW),LineStyle(2));
  dataMet.plotOn(wmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMet->GetBinWidth(1));
  CPlot plotMet("fitmet",wmframe,"","",ylabel);
  plotMet.SetLegend(0.68,0.57,0.93,0.77);
  plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMet.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#nu","F");
  plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMet.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotMet.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  plotMet.SetYRange(0.1,1.1*(hDataMet->GetMaximum()));
  plotMet.Draw(c,kFALSE,format,1);

  CPlot plotMetDiff("fitmet","","#slash{E}_{T} [GeV]","MC/Data");
  //CPlot plotMetDiff("fitmet","","mT [GeV]","#chi");
  plotMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
  plotMetDiff.SetYRange(-1,1);
  plotMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotMetDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotMetDiff.Draw(c,kTRUE,format,2);
  
  plotMet.SetName("fitmetlog");
  plotMet.SetLogy();
  plotMet.SetYRange(1e-3*(hDataMet->GetMaximum()),10*(hDataMet->GetMaximum()));
  plotMet.Draw(c,kTRUE,format,1);
    
  RooPlot *awmframe = pfmet.frame(Bins(NBINS));    
  antiMet.plotOn(awmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  apdfMet.plotOn(awmframe,FillColor(fillcolorW),DrawOption("F"));
  apdfMet.plotOn(awmframe,LineColor(linecolorW));
  apdfMet.plotOn(awmframe,Components(RooArgSet(apdfEWK,*(aqcd.model))),FillColor(fillcolorEWK),DrawOption("F"));
  apdfMet.plotOn(awmframe,Components(RooArgSet(apdfEWK,*(aqcd.model))),LineColor(linecolorEWK));
  apdfMet.plotOn(awmframe,Components(RooArgSet(*(aqcd.model))),FillColor(fillcolorQCD),DrawOption("F"));
  apdfMet.plotOn(awmframe,Components(RooArgSet(*(aqcd.model))),LineColor(linecolorQCD));
  apdfMet.plotOn(awmframe,Components(RooArgSet(apdfWm)),LineColor(linecolorW),LineStyle(2));
  antiMet.plotOn(awmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hAntiDataMet->GetBinWidth(1));
  CPlot plotAntiMet("fitantimet",awmframe,"","",ylabel);
  plotAntiMet.SetLegend(0.68,0.57,0.93,0.77);
  plotAntiMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotAntiMet.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#nu","F");
  plotAntiMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotAntiMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotAntiMet.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotAntiMet.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  plotAntiMet.SetYRange(0.1,1.1*(hAntiDataMet->GetMaximum())); 
  plotAntiMet.Draw(c,kFALSE,format,1);

  CPlot plotAntiMetDiff("fitantimet","","#slash{E}_{T} [GeV]","#chi");
  plotAntiMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
  plotAntiMetDiff.SetYRange(-0.2,0.2);
  plotAntiMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotAntiMetDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotAntiMetDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotAntiMetDiff.Draw(c,kTRUE,format,2);
  
  plotAntiMet.SetName("fitantimetlog");
  plotAntiMet.SetLogy();
  plotAntiMet.SetYRange(1e-3*(hAntiDataMet->GetMaximum()),10*(hAntiDataMet->GetMaximum()));
  plotAntiMet.Draw(c,kTRUE,format,1);
    
  //
  // W+ MET plot
  //
  RooPlot *wmpframe = pfmet.frame(Bins(NBINS));
  wmpframe->GetYaxis()->SetNdivisions(505);
  dataMetp.plotOn(wmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetp.plotOn(wmpframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetp.plotOn(wmpframe,LineColor(linecolorW));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),LineColor(linecolorEWK));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(*(qcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(*(qcdp.model))),LineColor(linecolorQCD));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfWmp)),LineColor(linecolorW),LineStyle(2));
  dataMetp.plotOn(wmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
  CPlot plotMetp("fitmetp",wmpframe,"","",ylabel);
  plotMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrow#mu^{+}#nu","F");
  plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetp.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotMetp.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
//  plotMetp.SetYRange(0.1,1.1*(hDataMetp->GetMaximum()));
  //plotMetp.SetYRange(0.1,4100);
  plotMetp.Draw(c,kFALSE,format,1);

  CPlot plotMetpDiff("fitmetp","","#slash{E}_{T} [GeV]","#chi");
  //CPlot plotMetpDiff("fitmetp","","mT [GeV]","#chi");
  plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
  plotMetpDiff.SetYRange(-0.2,0.2);
  plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetpDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotMetpDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotMetpDiff.Draw(c,kTRUE,format,2);
  
  plotMetp.SetName("fitmetplog");
  plotMetp.SetLogy();
  plotMetp.SetYRange(1e-3*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
  plotMetp.Draw(c,kTRUE,format,1);

  RooPlot *awmpframe = pfmet.frame(Bins(NBINS));    
  antiMetp.plotOn(awmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  apdfMetp.plotOn(awmpframe,FillColor(fillcolorW),DrawOption("F"));
  apdfMetp.plotOn(awmpframe,LineColor(linecolorW));
  apdfMetp.plotOn(awmpframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
  apdfMetp.plotOn(awmpframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),LineColor(linecolorEWK));
  apdfMetp.plotOn(awmpframe,Components(RooArgSet(*(aqcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
  apdfMetp.plotOn(awmpframe,Components(RooArgSet(*(aqcdp.model))),LineColor(linecolorQCD));
  apdfMetp.plotOn(awmpframe,Components(RooArgSet(apdfWmp)),LineColor(linecolorW),LineStyle(2));
  antiMetp.plotOn(awmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetp->GetBinWidth(1));
  CPlot plotAntiMetp("fitantimetp",awmpframe,"","",ylabel);
  plotAntiMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotAntiMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotAntiMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrow#mu^{+}#nu","F");
  plotAntiMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotAntiMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotAntiMetp.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotAntiMetp.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
//  plotAntiMetp.SetYRange(0.1,1.1*(hAntiDataMetp->GetMaximum()));
plotAntiMetp.SetYRange(0.1,1500);
  plotAntiMetp.Draw(c,kFALSE,format,1);

  CPlot plotAntiMetpDiff("fitantimetp","","#slash{E}_{T} [GeV]","#chi");
  plotAntiMetpDiff.AddHist1D(hAntiMetpDiff,"EX0",ratioColor);
  plotAntiMetpDiff.SetYRange(-0.2,0.2);
  plotAntiMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotAntiMetpDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotAntiMetpDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotAntiMetpDiff.Draw(c,kTRUE,format,2);
  
  plotAntiMetp.SetName("fitantimetplog");
  plotAntiMetp.SetLogy();
  plotAntiMetp.SetYRange(1e-3*(hAntiDataMetp->GetMaximum()),10*(hAntiDataMetp->GetMaximum()));
  plotAntiMetp.Draw(c,kTRUE,format,1);
  
  //
  // W- MET plot
  //
  RooPlot *wmmframe = pfmet.frame(Bins(NBINS)); 
  wmmframe->GetYaxis()->SetNdivisions(505);
  dataMetm.plotOn(wmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetm.plotOn(wmmframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetm.plotOn(wmmframe,LineColor(linecolorW));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),LineColor(linecolorEWK));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(*(qcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(*(qcdm.model))),LineColor(linecolorQCD));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfWmm)),LineColor(linecolorW),LineStyle(2));
  dataMetm.plotOn(wmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
  CPlot plotMetm("fitmetm",wmmframe,"","",ylabel);
  plotMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
  plotMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetm.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotMetm.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
//  plotMetm.SetYRange(0.1,1.1*(hDataMetm->GetMaximum()));
//plotMetm.SetYRange(0.1,4100);
  plotMetm.Draw(c,kFALSE,format,1);

  CPlot plotMetmDiff("fitmetm","","#slash{E}_{T} [GeV]","#chi");
  //CPlot plotMetmDiff("fitmetm","","mT [GeV]","#chi");
  plotMetmDiff.AddHist1D(hMetmDiff,"EX0",ratioColor);
  plotMetmDiff.SetYRange(-0.2,0.2);
  plotMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetmDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotMetmDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotMetmDiff.Draw(c,kTRUE,format,2);
  
  plotMetm.SetName("fitmetmlog");
  plotMetm.SetLogy();
  plotMetm.SetYRange(1e-3*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
  plotMetm.Draw(c,kTRUE,format,1);

  RooPlot *awmmframe = pfmet.frame(Bins(NBINS)); 
  antiMetm.plotOn(awmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  apdfMetm.plotOn(awmmframe,FillColor(fillcolorW),DrawOption("F"));
  apdfMetm.plotOn(awmmframe,LineColor(linecolorW));
  apdfMetm.plotOn(awmmframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
  apdfMetm.plotOn(awmmframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),LineColor(linecolorEWK));
  apdfMetm.plotOn(awmmframe,Components(RooArgSet(*(aqcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
  apdfMetm.plotOn(awmmframe,Components(RooArgSet(*(aqcdm.model))),LineColor(linecolorQCD));
  apdfMetm.plotOn(awmmframe,Components(RooArgSet(apdfWmm)),LineColor(linecolorW),LineStyle(2));
  antiMetm.plotOn(awmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
  CPlot plotAntiMetm("fitantimetm",awmmframe,"","",ylabel);
  plotAntiMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotAntiMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotAntiMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
  plotAntiMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotAntiMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotAntiMetm.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotAntiMetm.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
//  plotAntiMetm.SetYRange(0.1,1.1*(hAntiDataMetm->GetMaximum()));
plotAntiMetm.SetYRange(0.1,1500);
  plotAntiMetm.Draw(c,kFALSE,format,1);

  CPlot plotAntiMetmDiff("fitantimetm","","#slash{E}_{T} [GeV]","#chi");
  plotAntiMetmDiff.AddHist1D(hAntiMetmDiff,"EX0",ratioColor);
  plotAntiMetmDiff.SetYRange(-0.2,0.2);
  plotAntiMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotAntiMetmDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotAntiMetmDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotAntiMetmDiff.Draw(c,kTRUE,format,2);
  
  plotAntiMetm.SetName("fitantimetmlog");
  plotAntiMetm.SetLogy();
  plotAntiMetm.SetYRange(1e-3*(hAntiDataMetm->GetMaximum()),10*(hAntiDataMetm->GetMaximum()));
  plotAntiMetm.Draw(c,kTRUE,format,1);

    
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
  //txtfile << "  Signal: " << nSig.getVal() << " +/- " << nSig.getPropagatedError(*fitRes) << endl;
  //txtfile << "     QCD: " << nQCD.getVal() << " +/- " << nQCD.getPropagatedError(*fitRes) << endl;
  //txtfile << "   Other: " << nEWK.getVal() << " +/- " << nEWK.getPropagatedError(*fitRes) << endl;
  txtfile << endl;
  txtfile.flags(flags);
  
  //fitRes->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  //printCorrelations(txtfile, fitRes);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();
  
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
  //txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResp) << endl;
  //txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResp) << endl;
  //txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResp) << endl;
  txtfile << endl; 
  txtfile.flags(flags);
  
  //fitResp->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  //printCorrelations(txtfile, fitResp);
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
  txtfile << "  Signal: " << nSigm.getVal() << " +/- " << nSigm.getPropagatedError(*fitResm) << endl;
  txtfile << "     QCD: " << nQCDm.getVal() << " +/- " << nQCDm.getPropagatedError(*fitResm) << endl;
  txtfile << "   Other: " << nEWKm.getVal() << " +/- " << nEWKm.getPropagatedError(*fitResm) << endl;
  txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResm) << endl;
  txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResm) << endl;
  txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResm) << endl;
  
  txtfile << "reversed isolation sample" << endl;
  txtfile << "Selected: " << hDataMetm->Integral() << endl;
  txtfile << "  Signal: " << nAntiSigm.getVal() << " +/- " << nAntiSigm.getPropagatedError(*fitResm) << endl;
  txtfile << "     QCD: " << nAntiQCDm.getVal() << " +/- " << nAntiQCDm.getPropagatedError(*fitResm) << endl;
  txtfile << "   Other: " << nAntiEWKm.getVal() << " +/- " << nAntiEWKm.getPropagatedError(*fitResm) << endl;
  txtfile << "  Signal: " << nAntiSigp.getVal() << " +/- " << nAntiSigp.getPropagatedError(*fitResm) << endl;
  txtfile << "     QCD: " << nAntiQCDp.getVal() << " +/- " << nAntiQCDp.getPropagatedError(*fitResm) << endl;
  txtfile << "   Other: " << nAntiEWKp.getVal() << " +/- " << nAntiEWKp.getPropagatedError(*fitResm) << endl;
  txtfile << endl;
  txtfile.flags(flags);
  
  fitResm->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  fitResp->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  fitResAntim->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  fitResAntip->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResm);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();

  makeHTML(outputDir);
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  gBenchmark->Show("fitWm");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    
    Double_t diff = diff0/hData->GetBinContent(ibin);
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))*(hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
//     if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    
    hDiff->SetBinContent(ibin,diff);
//     else      hDiff->SetBinContent(ibin,0);
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