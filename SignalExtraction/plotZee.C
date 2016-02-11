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
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

//=== MAIN MACRO ================================================================================================= 

void plotZee(const TString  outputDir,   // output directory
             const Double_t lumi         // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotZee");
  gStyle->SetTitleOffset(1.100,"Y");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  const Double_t ELE_MASS  = 0.000511;
  //
  // input ntuple file names
  //
  enum { eData, eZee, eEWK };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;

  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_app4/Zee/ntuples/data_select.root"); typev.push_back(eData);
  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_app4/Zee/ntuples/zee_select.raw.root");   typev.push_back(eZee);
  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_app4/Zee/ntuples/top_select.raw.root");   typev.push_back(eEWK);
  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_app4/Zee/ntuples/ewk_select.raw.root");  typev.push_back(eEWK);
 
  //
  // Fit options
  //
  const Int_t    NBINS     = 60;
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  //const Int_t    NBINS     = 20;
  //const Double_t MASS_LOW  = 80;
  //const Double_t MASS_HIGH = 100;  
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.5;
  
  // plot output file format
  const TString format("pdf");

   // efficiency files
  const TString dataHLTEffName     = "/data/blue/cmedlock/wz-efficiency-results/DataZee_EleHLTEff/eff.root";
  const TString dataHLTEffName_pos = "/data/blue/cmedlock/wz-efficiency-results/DataZee_EleHLTEff/eff.root";
  const TString dataHLTEffName_neg = "/data/blue/cmedlock/wz-efficiency-results/DataZee_EleHLTEff/eff.root";
  const TString zeeHLTEffName      = "/data/blue/cmedlock/wz-efficiency-results/Zee_EleHLTEff/eff.root";
  const TString zeeHLTEffName_pos  = "/data/blue/cmedlock/wz-efficiency-results/Zee_EleHLTEff/eff.root";
  const TString zeeHLTEffName_neg  = "/data/blue/cmedlock/wz-efficiency-results/Zee_EleHLTEff/eff.root";
  
  const TString dataGsfSelEffName     = "/data/blue/cmedlock/wz-efficiency-results/DataZee_EleGsfSelEff/eff.root";
  const TString dataGsfSelEffName_pos = "/data/blue/cmedlock/wz-efficiency-results/DataZee_EleGsfSelEff/eff.root";
  const TString dataGsfSelEffName_neg = "/data/blue/cmedlock/wz-efficiency-results/DataZee_EleGsfSelEff/eff.root";
  const TString zeeGsfSelEffName      = "/data/blue/cmedlock/wz-efficiency-results/Zee_EleGsfSelEff/eff.root";
  const TString zeeGsfSelEffName_pos  = "/data/blue/cmedlock/wz-efficiency-results/Zee_EleGsfSelEff/eff.root";
  const TString zeeGsfSelEffName_neg  = "/data/blue/cmedlock/wz-efficiency-results/Zee_EleGsfSelEff/eff.root";

  const TString corrFiles = "../EleScale/76X_16DecRereco_2015";
  const bool doScaleAndSmear = true;

  //data
  EnergyScaleCorrection_class eleCorr( scaleFile); eleCorr.doScale= true; eleCorr.doSmear = true;
 

  Int_t yield = 0;
  Double_t yield_zee = 0, yield_zee_unc=0;
  Double_t yield_ewk = 0, yield_ewk_unc=0; 

  TString pufname = "../Tools/pileup_weights_2015B.root";
   
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  // event category enumeration
  enum { eEleEle2HLT=1, eEleEle1HLT1L1, eEleEle1HLT, eEleEleNoSel, eEleSC };
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/pileup_weights_2015B.root", "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("npv_rw");


  // histograms for full selection (EleEle2HLT + EleEle1HLT)
  TH1D *hData = new TH1D("hData","",NBINS,MASS_LOW,MASS_HIGH); hData->Sumw2();
  TH1D *hZee  = new TH1D("hZee", "",NBINS,MASS_LOW,MASS_HIGH); hZee->Sumw2();
  TH1D *hEWK  = new TH1D("hEWK", "",NBINS,MASS_LOW,MASS_HIGH); hEWK->Sumw2();
  TH1D *hMC   = new TH1D("hMC",  "",NBINS,MASS_LOW,MASS_HIGH); hMC->Sumw2();
    
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, u1, u2;
  Float_t puWeight;
  Int_t   q1, q2;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0;
  TLorentzVector *sc1=0, *sc2=0;
  Float_t r91,r92;
  
  TH2D *h=0;
   //
  // HLT efficiency
  //
  cout << "Loading trigger efficiencies..." << endl;

  TFile *dataHLTEffFile_pos = new TFile(dataHLTEffName_pos);
  CEffUser2D dataHLTEff_pos;
  dataHLTEff_pos.loadEff((TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEffFile_neg = new TFile(dataHLTEffName_neg);
  CEffUser2D dataHLTEff_neg;
  dataHLTEff_neg.loadEff((TH2D*)dataHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zeeHLTEffFile_pos = new TFile(zeeHLTEffName_pos);
  CEffUser2D zeeHLTEff_pos;
  zeeHLTEff_pos.loadEff((TH2D*)zeeHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)zeeHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zeeHLTEffFile_neg = new TFile(zeeHLTEffName_neg);
  CEffUser2D zeeHLTEff_neg;
  zeeHLTEff_neg.loadEff((TH2D*)zeeHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)zeeHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeHLTEffFile_neg->Get("hErrhEtaPt"));
  
  h =(TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt");
  TH2D *hHLTErr_pos = new TH2D("hHLTErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                                 h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hHLTErr_neg = new TH2D("hHLTErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                                 h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  
  //
  // Selection efficiency
  //
  cout << "Loading GSF+selection efficiencies..." << endl;
  
  TFile *dataGsfSelEffFile_pos = new TFile(dataGsfSelEffName_pos);
  CEffUser2D dataGsfSelEff_pos;
  dataGsfSelEff_pos.loadEff((TH2D*)dataGsfSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataGsfSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataGsfSelEffFile_neg = new TFile(dataGsfSelEffName_neg);
  CEffUser2D dataGsfSelEff_neg;
  dataGsfSelEff_neg.loadEff((TH2D*)dataGsfSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataGsfSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEffFile_neg->Get("hErrhEtaPt"));

  TFile *zeeGsfSelEffFile_pos = new TFile(zeeGsfSelEffName_pos);
  CEffUser2D zeeGsfSelEff_pos;
  zeeGsfSelEff_pos.loadEff((TH2D*)zeeGsfSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEffFile_pos->Get("hErrhEtaPt"));

  TFile *zeeGsfSelEffFile_neg = new TFile(zeeGsfSelEffName_neg);
  CEffUser2D zeeGsfSelEff_neg;
  zeeGsfSelEff_neg.loadEff((TH2D*)zeeGsfSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEffFile_neg->Get("hErrhEtaPt"));
 
  h =(TH2D*)dataGsfSelEffFile_pos->Get("hEffEtaPt");
  TH2D *hGsfSelErr_pos = new TH2D("hGsfSelErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                                       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hGsfSelErr_neg = new TH2D("hGsfSelErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                                       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *markus = new TH2D("plot","plot",60,60,120,60,0.8,1.2);


  TFile *infile=0;
  TTree *intree=0;

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
  
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);	    assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);

    intree->SetBranchAddress("runNum",   &runNum);     // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);    // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);     // event number
    intree->SetBranchAddress("matchGen", &matchGen);   // event has both leptons matched to MC Z->ll
    intree->SetBranchAddress("category", &category);   // dilepton category
    intree->SetBranchAddress("npv",      &npv);	       // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);	       // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);     // GEN Z boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);    // GEN Z boson phi (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);   // event weight per 1/fb (MC)
    intree->SetBranchAddress("met",      &met);	       // MET
    intree->SetBranchAddress("metPhi",   &metPhi);     // phi(MET)
    intree->SetBranchAddress("sumEt",    &sumEt);      // Sum ET
    intree->SetBranchAddress("u1",       &u1);	       // parallel component of recoil
    intree->SetBranchAddress("u2",       &u2);	       // perpendicular component of recoil
    intree->SetBranchAddress("q1",       &q1);	       // charge of tag lepton
    intree->SetBranchAddress("q2",       &q2);	       // charge of probe lepton
    intree->SetBranchAddress("dilep",    &dilep);      // dilepton 4-vector
    intree->SetBranchAddress("lep1",     &lep1);       // tag lepton 4-vector
    intree->SetBranchAddress("lep2",     &lep2);       // probe lepton 4-vector
    intree->SetBranchAddress("sc1",      &sc1);        // tag Supercluster 4-vector
    intree->SetBranchAddress("sc2",      &sc2);        // probe Supercluster 4-vector 
    intree->SetBranchAddress("r91",       &r91);	       // r9
    intree->SetBranchAddress("r92",       &r92);	       // r9
    intree->SetBranchAddress("puWeight",      &puWeight);        // pu weight
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);

      if(fabs(lep1->Eta()) > ETA_CUT)   continue;      
      if(fabs(lep2->Eta()) > ETA_CUT)   continue;


      if ( doScaleAndSmear )
      {

      	float smear1 = 0.0, scale1 = 1.0;
      	float aeta1= fabs(lep1->Eta());
      	float et1 = lep1.E() / cosh(aeta1);
      	bool eb1 = aeta1 < 1.4442;
      	float smear2 = 0.0, scale2 = 1.0;
      	float aeta2= fabs(lep2->Eta());
      	float et2 = lep2.E() / cosh(aeta2);
      	bool eb2 = aeta2 < 1.4442;

       	if( typev[ifile]!=eData)
       	{// DATA
       		scale1 = eleCorr.ScaleCorrection( runNum, eb1, r91, aeta1, et1);
       		scale2 = eleCorr.ScaleCorrection( runNum, eb2, r92, aeta2, et2);
		//newEcalEnergy      = electron.getNewEnergy() * scale;
		//newEcalEnergyError = std::hypot(electron.getNewEnergyError() * scale, smear * newEcalEnergy);
		lep1 *= scale;
		lep2 *= scale;
       	}
       	else
       	{ // MC
       		smear1 = eleCorr.getSmearingSigma( runNum, lep1->E(), eb1, r91, aeta1); 
       		smear2 = eleCorr.getSmearingSigma( runNum, lep2->E(), eb2, r92, aeta2); 
		double corr1 = 1.0 + smear1 * gRandom->Gaus(0,1);
		double corr2 = 1.0 + smear2 * gRandom->Gaus(0,1);
		//float newEnergy = lep1->E() * corr1;
		//float newEcalEnergyError = std::hypot(electron.getNewEnergyError() * corr, smear * newEcalEnergy);
		lep1 * = corr;
		lep2 * = corr;
       	}
      } // do scale and smear

      //Double_t lp1 = gRandom->Gaus(lep1->Pt()*getEleScaleCorr(lep1->Eta(),1), getEleResCorr(lep1->Eta(),-1));
      //Double_t lp2 = gRandom->Gaus(lep2->Pt()*getEleScaleCorr(lep2->Eta(),1), getEleResCorr(lep2->Eta(),-1));
      //TLorentzVector l1, l2;
      //l1.SetPtEtaPhiM(lp1,lep1->Eta(),lep1->Phi(),ELE_MASS);
      //l2.SetPtEtaPhiM(lp2,lep2->Eta(),lep2->Phi(),ELE_MASS);
      //double mass=(l1+l2).M();
      Float_t mass = dilep->M();
      
      Double_t weight=1;
      if(typev[ifile]!=eData) {
	weight *= scale1fb*lumi;
	weight*= h_rw->GetBinContent(npv+1);
      }
      
      // fill Z events passing selection (EleEle2HLT + EleEle1HLT)
      if((category==eEleEle2HLT) || (category==eEleEle1HLT) || (category==eEleEle1HLT1L1)) {
      if(typev[ifile]==eData) { 
	 
	if(dilep->M()       < MASS_LOW)  continue;
	if(dilep->M()       > MASS_HIGH) continue;
	if(lep1->Pt()        < PT_CUT)    continue;
	if(lep2->Pt()        < PT_CUT)    continue;
	
	//if(mass       < MASS_LOW)  continue;
	//if(mass       > MASS_HIGH) continue;
	//if(lp1        < PT_CUT)    continue;
	//if(lp2        < PT_CUT)    continue;
	  //if(q1*q2>0)  continue;
	hData->Fill(mass); 
	
	yield++;
	  
      } else {
	Double_t lp1 = gRandom->Gaus(lep1->Pt()*getEleScaleCorr(lep1->Eta(),0), getEleResCorr(lep1->Eta(),0));
	Double_t lp2 = gRandom->Gaus(lep2->Pt()*getEleScaleCorr(lep2->Eta(),0), getEleResCorr(lep2->Eta(),0));
	  TLorentzVector l1, l2;
	  l1.SetPtEtaPhiM(lp1,lep1->Eta(),lep1->Phi(),ELE_MASS);
	  l2.SetPtEtaPhiM(lp2,lep2->Eta(),lep2->Phi(),ELE_MASS);
	  double mll=(l1+l2).M();
	  if(mll       < MASS_LOW)  continue;
	  if(mll       > MASS_HIGH) continue;
	  if(lp1        < PT_CUT)    continue;
	  if(lp2        < PT_CUT)    continue;
	  Double_t effdata, effmc;
	  Double_t corr=1;
	  effdata=1; effmc=1;
          if(q1>0) { 
            effdata *= (1.-dataHLTEff_pos.getEff(sc1->Eta(), sc1->Pt())); 
            effmc   *= (1.-zeeHLTEff_pos.getEff(sc1->Eta(), sc1->Pt())); 
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(sc1->Eta(), sc1->Pt())); 
            effmc   *= (1.-zeeHLTEff_neg.getEff(sc1->Eta(), sc1->Pt())); 
          }
          if(q2>0) {
            effdata *= (1.-dataHLTEff_pos.getEff(sc2->Eta(), sc2->Pt())); 
            effmc   *= (1.-zeeHLTEff_pos.getEff(sc2->Eta(), sc2->Pt()));
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(sc2->Eta(), sc2->Pt())); 
            effmc   *= (1.-zeeHLTEff_neg.getEff(sc2->Eta(), sc2->Pt()));
          }
          effdata = 1.-effdata;
          effmc   = 1.-effmc;
          corr *= effdata/effmc;
    
          effdata=1; effmc=1;
          if(q1>0) { 
            effdata *= dataGsfSelEff_pos.getEff(sc1->Eta(), sc1->Pt()); 
            effmc   *= zeeGsfSelEff_pos.getEff(sc1->Eta(), sc1->Pt()); 
          } else {
            effdata *= dataGsfSelEff_neg.getEff(sc1->Eta(), sc1->Pt()); 
            effmc   *= zeeGsfSelEff_neg.getEff(sc1->Eta(), sc1->Pt()); 
          }
          if(q2>0) {
            effdata *= dataGsfSelEff_pos.getEff(sc2->Eta(), sc2->Pt()); 
            effmc   *= zeeGsfSelEff_pos.getEff(sc2->Eta(), sc2->Pt());
          } else {
            effdata *= dataGsfSelEff_neg.getEff(sc2->Eta(), sc2->Pt()); 
            effmc   *= zeeGsfSelEff_neg.getEff(sc2->Eta(), sc2->Pt());
          }
          corr *= effdata/effmc;
	  corr=1;
	  //TLorentzVector slep1 = (*lep1);
	  //slep1 *= gRandom->Gaus(slep1.Pt(), getEleResCorr(sc1->Eta(),0))/slep1.Pt();
	  
	  //TLorentzVector slep2 = (*lep2);
	  //slep2 *= gRandom->Gaus(slep2.Pt(), getEleResCorr(sc2->Eta(),0))/slep2.Pt();
	  mass = (l1+l2).M();	
	  
	  if(typev[ifile]==eZee) 
	    { 
	      yield_zee += weight*corr;
	      yield_zee_unc += weight*weight*corr*corr;
	      hZee->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);
	      markus->Fill(mass, corr);
	    }
	  if(typev[ifile]==eEWK) 
	    { 
	      yield_ewk += weight*corr;
	      yield_ewk_unc += weight*weight*corr*corr;
	      hEWK->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);
	    }
	}
      }
    }
    
    delete infile;
    infile=0, intree=0;
  }

  hZee->Scale((hData->Integral()-hEWK->Integral())/hZee->Integral());
  hMC->Scale(hData->Integral()/hMC->Integral());

  TH1D *hZeeDiff = makeDiffHist(hData,hMC,"hZeeDiff");
  hZeeDiff->SetMarkerStyle(kFullCircle); 
  hZeeDiff->SetMarkerSize(0.9);
  


  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  char ylabel[100];     // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  sprintf(lumitext,"%.0f pb^{-1}  (13 TeV)",lumi);  
  
  // plot colors
  Int_t linecolorZ   = kOrange-3;
  Int_t fillcolorZ   = kOrange-2;
  Int_t linecolorEWK = kOrange+10;
  Int_t fillcolorEWK = kOrange+7;
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
  
  //
  // EleEle2HLT + EleEle1HLT categories
  //   
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hData->GetBinWidth(1));
  CPlot plotZee("zee","","",ylabel);
  plotZee.AddHist1D(hData,"data","E");
  plotZee.AddToStack(hZee,"Z#rightarrowee",fillcolorZ,linecolorZ);
  //plotZee.AddTextBox(lumitext,0.60,0.92,0.92,0.99,0);
  //plotZee.AddTextBox("CMS Preliminary",0.555,0.80,0.905,0.86,0);
  plotZee.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZee.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZee.SetYRange(0.01,1.2*(hData->GetMaximum() + sqrt(hData->GetMaximum())));
  plotZee.TransLegend(-0.35,-0.15);
  plotZee.Draw(c,kFALSE,format,1);

  CPlot plotZeeDiff("zee","","M(e^{+}e^{-}) [GeV/c^{2}]","#frac{Data-Pred}{Data}");
  plotZeeDiff.AddHist1D(hZeeDiff,"EX0",ratioColor);
  plotZeeDiff.SetYRange(-0.2,0.2);
  plotZeeDiff.AddLine(MASS_LOW, 0,MASS_HIGH, 0,kBlack,1);
  plotZeeDiff.AddLine(MASS_LOW, 0.1,MASS_HIGH, 0.1,kBlack,3);
  plotZeeDiff.AddLine(MASS_LOW,-0.1,MASS_HIGH,-0.1,kBlack,3);
  plotZeeDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZee2("zeelog","","",ylabel);
  plotZee2.AddHist1D(hData,"data","E");
  plotZee2.AddToStack(hEWK,"EWK+t#bar{t}",fillcolorEWK,linecolorEWK);
  plotZee2.AddToStack(hZee,"Z#rightarrowee",fillcolorZ,linecolorZ);
  plotZee2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZee2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  //plotZee2.AddTextBox(lumitext,0.60,0.91,0.92,0.98,0);plotZee2.SetName("zeelog");
  //plotZee2.AddTextBox("CMS Preliminary",0.555,0.80,0.905,0.86,0);
  plotZee2.SetLogy();
  plotZee2.SetYRange(1e-4*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plotZee2.TransLegend(-0.395,-0.15);
  plotZee2.Draw(c,kTRUE,format,1);

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  cout << endl;

  cout << " Total Zee event yield is " << yield << "." << endl;
  cout << " The Zee expected event yield is " << yield_zee << " +/-" << sqrt(yield_zee_unc) << "." << endl;
  cout << " The EWK event yield is " << yield_ewk << " +/-" << sqrt(yield_ewk_unc) << "." << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  gBenchmark->Show("plotZee");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
/*TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    
    Double_t err = sqrt(hData->GetBinContent(ibin));
    if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    
    if(err>0) hDiff->SetBinContent(ibin,diff/err);
    else      hDiff->SetBinContent(ibin,0);
    hDiff->SetBinError(ibin,1);   
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
  }*/
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    Double_t diff = diff0/hData->GetBinContent(ibin);
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((1.0/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    //Double_t err = sqrt(hData->GetBinContent(ibin));
    //if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    //if(err>0) hDiff->SetBinContent(ibin,diff/err);
    //else      hDiff->SetBinContent(ibin,0);
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
