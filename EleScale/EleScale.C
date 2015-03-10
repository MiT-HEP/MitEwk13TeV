#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                    // access to gROOT, entry point to ROOT system
#include <TSystem.h>                  // interface to OS
#include <TStyle.h>                   // class to handle ROOT plotting style
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TCanvas.h>                  // class for drawing
#include <TH1D.h>                     // 1D histograms
#include <TGraphErrors.h>             // graphs
#include <vector>                     // STL vector class
#include <utility>                    // For STL pair class
#include <map>                        // STL map class
#include <iostream>                   // standard I/O
#include <iomanip>                    // functions to format standard I/O
#include <fstream>                    // functions for file I/O
#include <string>                     // C++ string class
#include <sstream>                    // class for parsing strings
#include "Math/LorentzVector.h"       // 4-vector class

#include "../Utils/CPlot.hh"          // helper class for plots
#include "../Utils/MitStyleRemix.hh"  // style settings for drawing
#endif

// RooFit headers
#include "RooRealVar.h"
#include "RooLinearVar.h"
#include "RooFormulaVar.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooFitResult.h"


typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void EleScale() {
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  // event category enumeration
  enum { eEleEle2HLT=1, eEleEle1HLT, eEleEleNoSel, eEleSC };
  
  TString outputDir = "Results";
  TString pufname = ""; 
  
  vector<TString> infilenamev;
  infilenamev.push_back("/scratch/ksung/EWKAna/8TeV/Selection/Zee/ntuples/data_select.raw.root");  // data
  infilenamev.push_back("/scratch/ksung/EWKAna/8TeV/Selection/Zee/ntuples/zee_select.root");       // MC
  
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.5;
  const Double_t ELE_MASS  = 0.000511;  
  
  vector<pair<Double_t,Double_t> > scEta_limits;
  scEta_limits.push_back(make_pair(0.0,0.4));
  scEta_limits.push_back(make_pair(0.4,0.8));
  scEta_limits.push_back(make_pair(0.8,1.2));
  scEta_limits.push_back(make_pair(1.2,1.4442));
  scEta_limits.push_back(make_pair(1.566,2.0));
  scEta_limits.push_back(make_pair(2.0,2.5));

  CPlot::sOutDir = outputDir;
  
  const TString format("png");
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
   
  enum { eData=0, eMC };
  
//  TFile *pufile = new TFile(pufname); assert(pufile);
//  TH1D  *puWeights = (TH1D*)pufile->Get("puWeights");
  
  char hname[100];
  vector<TH1D*> hMCv, hDatav;  
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    for(UInt_t jbin=ibin; jbin<scEta_limits.size(); jbin++) {
      sprintf(hname,"mc_%i_%i",ibin,jbin);
      hMCv.push_back(new TH1D(hname,"",100,MASS_LOW,MASS_HIGH));
      hMCv.back()->Sumw2();
      
      sprintf(hname,"data_%i_%i",ibin,jbin);
      hDatav.push_back(new TH1D(hname,"",100,MASS_LOW,MASS_HIGH));
      hDatav.back()->Sumw2();
    }
  }
  
  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Int_t   q1, q2;
  LorentzVector *dilep=0, *lep1=0, *lep2=0;
  ///// electron specific /////
  LorentzVector *sc1=0, *sc2=0;
  
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    cout << "Processing " << infilenamev[ifile] << "..." << endl;
    TFile *infile = new TFile(infilenamev[ifile]); assert(infile);
    TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
  
    intree->SetBranchAddress("runNum",   &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    intree->SetBranchAddress("matchGen", &matchGen);  // event has both leptons matched to MC Z->ll
    intree->SetBranchAddress("category", &category);  // dilepton category
    intree->SetBranchAddress("npv",      &npv);	      // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);	      // number of in-time PU events (MC)
    intree->SetBranchAddress("q1",       &q1);	      // charge of lead lepton
    intree->SetBranchAddress("q2",       &q2);	      // charge of trail lepton
    intree->SetBranchAddress("dilep",    &dilep);     // dilepton 4-vector
    intree->SetBranchAddress("lep1",     &lep1);      // lead lepton 4-vector
    intree->SetBranchAddress("lep2",     &lep2);      // trail lepton 4-vector
    intree->SetBranchAddress("sc1",      &sc1);	      // lead Supercluster 4-vector
    intree->SetBranchAddress("sc2",      &sc2);	      // trail Supercluster 4-vector 
  
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      
      Double_t weight = 1;
      if(ifile==eMC) {
        if(!matchGen) continue;
//        weight *= puWeights->GetBinContent(npu+1);
      }
      
      if((category!=eEleEle2HLT) && (category!=eEleEle1HLT)) continue;
      if(q1 == q2) continue;
      if(dilep->M()	  < MASS_LOW)  continue;
      if(dilep->M()	  > MASS_HIGH) continue;
      if(sc1->Pt()	  < PT_CUT)    continue;
      if(sc2->Pt()	  < PT_CUT)    continue;
      if(fabs(sc1->Eta()) > ETA_CUT)   continue;      
      if(fabs(sc2->Eta()) > ETA_CUT)   continue;
    
      LorentzVector vLep1(sc1->Pt(), lep1->Eta(), lep1->Phi(), ELE_MASS);
      LorentzVector vLep2(sc2->Pt(), lep2->Eta(), lep2->Phi(), ELE_MASS);
      //LorentzVector vLep1 = (*lep1);
      //LorentzVector vLep2 = (*lep2);
      LorentzVector vDilep = vLep1 + vLep2;
    
      Int_t bin1=-1, bin2=-1;
      for(UInt_t i=0; i<scEta_limits.size(); i++) {
        Double_t etalow  = scEta_limits.at(i).first;
        Double_t etahigh = scEta_limits.at(i).second;
        if(fabs(sc1->Eta())>=etalow && fabs(sc1->Eta())<=etahigh) bin1=i;
        if(fabs(sc2->Eta())>=etalow && fabs(sc2->Eta())<=etahigh) bin2=i;
      }
      assert(bin1>=0);
      assert(bin2>=0);
      Int_t ibin= (bin1<=bin2) ? bin1 : bin2;
      Int_t jbin= (bin1<=bin2) ? bin2 : bin1;
      
      UInt_t n=jbin-ibin;
      for(Int_t k=0; k<ibin; k++)
        n+=(scEta_limits.size()-k);
      
      if(ifile==eData) hDatav[n]->Fill(vDilep.M(),weight);
      if(ifile==eMC)   hMCv[n]->Fill(vDilep.M(),weight);
    }  
    delete infile;
    infile=0, intree=0;
  }
  
  //
  // Fit for energy scale and resolution corrections
  //
  char vname[100];  // buffer for RooFit object names
  
  char pname[100];
  char str1[100];
  char str2[100];
  TCanvas *c = MakeCanvas("c","c",800,600);
  
  // Dummy histograms for TLegend (I can't figure out how to properly pass RooFit objects...)
  TH1D *hDummyData = new TH1D("hDummyData","",0,0,10);
  hDummyData->SetMarkerStyle(kFullCircle);
  hDummyData->SetMarkerSize(0.9);
  TH1D *hDummyMC = new TH1D("hDummyMC","",0,0,10);
  hDummyMC->SetLineColor(kBlue);
  hDummyMC->SetFillColor(kBlue);
  hDummyMC->SetFillStyle(3002);  
  TH1D *hDummyFit = new TH1D("hDummyFit","",0,0,10);
  hDummyFit->SetLineColor(kGreen+2);

  RooRealVar mass("mass","M_{ee}",60.0,120.0,"GeV") ;
  mass.setBins(1600,"cache");
  
  RooRealVar massmc("massmc","massmc",0.0,150.0,"GeV");  // mass variable for building MC template

  RooCategory zscEta_cat("zscEta_cat","zscEta_cat");  
  RooSimultaneous combscalefit("combscalefit","combscalefit",zscEta_cat);
  
  map<string,TH1*> hmap;  // Mapping of category labels and data histograms
  
  RooArgList scalebins;   // List of RooRealVars storing per bin energy scale corrections
  RooArgList sigmabins;   // List of RooRealVars storing per bin energy resolution corrections
  Int_t intOrder = 1;     // Interpolation order for       
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    sprintf(vname,"scale_%i",ibin);
    RooRealVar *scalebinned = new RooRealVar(vname,vname,1.0,0.8,1.2);      
    scalebins.add(*scalebinned);
    
    sprintf(vname,"sigma_%i",ibin);
    RooRealVar *sigmabinned = new RooRealVar(vname,vname,0.8,0.0,5.0);      
    sigmabins.add(*sigmabinned);
  }
    
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    for(UInt_t jbin=ibin; jbin<scEta_limits.size(); jbin++) {
      UInt_t n=jbin-ibin;
      for(UInt_t k=0; k<ibin; k++)
        n+=(scEta_limits.size()-k);
      
      sprintf(vname,"masslinearshifted_%i_%i",ibin,jbin);
      RooFormulaVar *masslinearshifted = new RooFormulaVar(vname,vname,"sqrt(@0*@1)",RooArgList(*scalebins.at(ibin),*scalebins.at(jbin)));

      sprintf(vname,"massshiftedscEta_%i_%i",ibin,jbin);
      RooLinearVar *massshiftedscEta = new RooLinearVar(vname,vname,mass,*masslinearshifted,RooConst(0.0));        
      
      // MC-based template
      sprintf(vname,"zmassmcscEta_%i_%i",ibin,jbin);
      RooDataHist *zmassmcscEta = new RooDataHist(vname,vname,RooArgList(massmc),hMCv[n]);      
      sprintf(vname,"masstemplatescEta_%i_%i",ibin,jbin);
      RooHistPdf *masstemplatescEta = new RooHistPdf(vname,vname,RooArgList(*massshiftedscEta),RooArgList(massmc),*zmassmcscEta,intOrder);             

      // Gaussian smearing function
      sprintf(vname,"sigmascEta_%i_%i",ibin,jbin);
      RooFormulaVar *sigmascEta = new RooFormulaVar(vname,vname,"TMath::Max(0.01,sqrt(@0*@0+@1*@1))",RooArgList(*sigmabins.at(ibin),*sigmabins.at(jbin)));        
      sprintf(vname,"resscEta_%i_%i",ibin,jbin);
      RooGaussian *resscEta = new RooGaussian(vname,vname,mass,RooConst(0.),*sigmascEta);

      // Fit model: MC-template convoluted with Gaussian
      sprintf(vname,"fftscEta_%i_%i",ibin,jbin);
      RooFFTConvPdf *fftscEta = new RooFFTConvPdf(vname,vname,mass,*masstemplatescEta,*resscEta);
      fftscEta->setBufferStrategy(RooFFTConvPdf::Flat);
      
      // Add bin as a category
      char zscEta_catname[100];
      sprintf(zscEta_catname,"zscEta_cat_%i_%i",ibin,jbin);
      zscEta_cat.defineType(zscEta_catname); 
      zscEta_cat.setLabel(zscEta_catname);
      hmap.insert(pair<string,TH1*>(zscEta_catname,hDatav[n]));      
      combscalefit.addPdf(*fftscEta,zscEta_catname);
    }
  }
  
  // perform fit
  RooDataHist zdatascEta_comb("zdatascEta_comb","zdatascEta_comb",RooArgList(mass),zscEta_cat,hmap,1.0/(hDatav.size()));
  combscalefit.fitTo(zdatascEta_comb,PrintEvalErrors(kFALSE),Minos(kFALSE),Strategy(0),Minimizer("Minuit2",""));

  Double_t xval[scEta_limits.size()];
  Double_t xerr[scEta_limits.size()];
  Double_t scaleDatatoMC[scEta_limits.size()];
  Double_t scaleDatatoMCerr[scEta_limits.size()];
  Double_t scaleMCtoData[scEta_limits.size()];
  Double_t scaleMCtoDataerr[scEta_limits.size()];
  Double_t sigmaMCtoData[scEta_limits.size()];
  Double_t sigmaMCtoDataerr[scEta_limits.size()];
  
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    Double_t etalow  = scEta_limits.at(ibin).first;
    Double_t etahigh = scEta_limits.at(ibin).second;
    
    xval[ibin] = 0.5*(etahigh+etalow);
    xerr[ibin] = 0.5*(etahigh-etalow);
    
    scaleDatatoMC[ibin]    = ((RooRealVar*)scalebins.at(ibin))->getVal();
    scaleDatatoMCerr[ibin] = ((RooRealVar*)scalebins.at(ibin))->getError();

    scaleMCtoData[ibin]    = 1.0/scaleDatatoMC[ibin];
    scaleMCtoDataerr[ibin] = scaleDatatoMCerr[ibin]/scaleDatatoMC[ibin]/scaleDatatoMC[ibin];
    
    sigmaMCtoData[ibin]    = ((RooRealVar*)sigmabins.at(ibin))->getVal();
    sigmaMCtoDataerr[ibin] = ((RooRealVar*)sigmabins.at(ibin))->getError();
  }
  TGraphErrors *grScaleDatatoMC = new TGraphErrors(scEta_limits.size(),xval,scaleDatatoMC,xerr,scaleDatatoMCerr);
  TGraphErrors *grScaleMCtoData = new TGraphErrors(scEta_limits.size(),xval,scaleMCtoData,xerr,scaleMCtoDataerr);
  TGraphErrors *grSigmaMCtoData = new TGraphErrors(scEta_limits.size(),xval,sigmaMCtoData,xerr,sigmaMCtoDataerr);
  
  CPlot plotScale1("ele_scale_datatomc","","Supercluster |#eta|","Data scale correction");
  plotScale1.AddGraph(grScaleDatatoMC,"",kBlue);
  plotScale1.SetYRange(0.99,1.01);
  plotScale1.AddLine(0,1,2.75,1,kBlack,7);
  plotScale1.Draw(c,kTRUE,format);
  
  CPlot plotScale2("ele_scale_mctodata","","Supercluster |#eta|","MC#rightarrowData scale correction");
  plotScale2.AddGraph(grScaleMCtoData,"",kBlue);
  plotScale2.SetYRange(0.99,1.01);
  plotScale2.AddLine(0,1,2.75,1,kBlack,7);
  plotScale2.Draw(c,kTRUE,format);

  CPlot plotRes("ele_res_mctodata","","Supercluster |#eta|","MC#rightarrowData additional smear [GeV]");
  plotRes.AddGraph(grSigmaMCtoData,"",kBlue);
  plotRes.SetYRange(0,1.6);
  plotRes.Draw(c,kTRUE,format);

  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    for(UInt_t jbin=ibin; jbin<scEta_limits.size(); jbin++) {
      UInt_t n=jbin-ibin;
      for(UInt_t k=0; k<ibin; k++)
        n+=(scEta_limits.size()-k);
      
      // Post-fit plot
      RooPlot *frame = mass.frame();
      char catname[100]; sprintf(catname,"zscEta_cat_%i_%i",ibin,jbin);
      char cutstr[100];  sprintf(cutstr,"zscEta_cat==zscEta_cat::%s",catname); 
      RooDataHist zmc(catname,catname,RooArgList(mass),hMCv[n]);
      RooHistPdf mctemplate(catname,catname,RooArgList(mass),zmc,intOrder);          
      mctemplate.plotOn(frame,LineColor(kBlue),LineWidth(1),Normalization(hDatav[n]->GetEntries()));            
      mctemplate.plotOn(frame,LineColor(kBlue),FillColor(kBlue),FillStyle(3002),DrawOption("F"),Normalization(hDatav[n]->GetEntries()));      
      zdatascEta_comb.plotOn(frame,Cut(cutstr),MarkerStyle(kFullCircle),MarkerSize(1.0),DrawOption("ZP"));
      combscalefit.plotOn(frame,Slice(zscEta_cat,catname),ProjWData(RooArgSet(mass,catname),zdatascEta_comb),
                          LineColor(kGreen+2));           
      sprintf(pname,"postfit_%i_%i",ibin,jbin);
      sprintf(str1,"[%.1f, %.1f]",scEta_limits.at(ibin).first,scEta_limits.at(ibin).second);
      sprintf(str2,"[%.1f, %.1f]",scEta_limits.at(jbin).first,scEta_limits.at(jbin).second);
      CPlot plot(pname,frame,"","m(e^{+}e^{-}) [GeV/c^{2}]","Events / 0.6 GeV/c^{2}");
      plot.AddTextBox(str1,0.21,0.80,0.45,0.87,0,kBlack,-1);
      plot.AddTextBox(str2,0.21,0.73,0.45,0.80,0,kBlack,-1);
      plot.SetLegend(0.75,0.64,0.93,0.88);
      plot.GetLegend()->AddEntry(hDummyData,"Data","PL");
      plot.GetLegend()->AddEntry(hDummyMC,"Sim","FL");
      plot.GetLegend()->AddEntry(hDummyFit,"Fit","L");
      plot.Draw(c,kTRUE,format);
    }
  }

      
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
   
  ofstream txtfile;
  char txtfname[100];    
  sprintf(txtfname,"%s/summary.txt",outputDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  txtfile << "  Data->MC scale correction" << endl;
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    Double_t etalow  = scEta_limits.at(ibin).first;
    Double_t etahigh = scEta_limits.at(ibin).second;
    txtfile << "$" << etalow << " < |\\eta| < " << etahigh << "$ & ";
    txtfile << "$" << ((RooRealVar*)scalebins.at(ibin))->getVal() << "$ \\pm $" << ((RooRealVar*)scalebins.at(ibin))->getError() << "$ \\\\" << endl;
  }
  txtfile << endl;
  txtfile << "  MC->Data resolution correction [GeV]" << endl;
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    Double_t etalow  = scEta_limits.at(ibin).first;
    Double_t etahigh = scEta_limits.at(ibin).second;
    txtfile << etalow << " < |\\eta| < " << etahigh << " & ";
    txtfile << "$" << ((RooRealVar*)sigmabins.at(ibin))->getVal() << "$ \\pm $" << ((RooRealVar*)sigmabins.at(ibin))->getError() << "$ \\\\" << endl;
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
}
