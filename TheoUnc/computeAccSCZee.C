//================================================================================================
//
// Compute Z->ee acceptance at supercluster level
//
//  * outputs results summary text file
//
//________________________________________________________________________________________________

#if !defined(__CINT__) //|| defined(__MAKECINT__)
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
#include "LHAPDF/LHAPDF.h"

// various helper functions
#include "../Utils/MyTools.hh"

// define structures to read in ntuple
//#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

Double_t ecalEta(const Double_t eta, const Double_t z, const Double_t rho);


//=== MAIN MACRO ================================================================================================= 

void computeAccSCZee(const TString conf="zee.conf",       // input file
                     const TString outputDir="idfk"   // output directory
) {
  gBenchmark->Start("computeAccSCZee");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.5;
  const Double_t ELE_MASS   = 0.000511;  
  const Double_t ETA_BARREL = 1.4442;
  const Double_t ETA_ENDCAP = 1.566;
  
  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 11;

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
 
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info  = new baconhep::TEventInfo();
  baconhep::TGenEventInfo   *gen   = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr  = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *photonArr   = new TClonesArray("baconhep::TPhoton");
  
  TFile *infile=0;
  TTree *eventTree=0;
  
  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv, nSelBBv, nSelBEv, nSelEEv;
  vector<Double_t> accv, accBBv, accBEv, accEEv;
  vector<Double_t> accErrv, accErrBBv, accErrBEv, accErrEEv;

  LHAPDF::PDF* nomPdf = LHAPDF::mkPDF(292200);
  //LHAPDF::PDF* testPdf = LHAPDF::mkPDF(13100);
  LHAPDF::PDF* testPdf = LHAPDF::mkPDF(25100);
    
  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = TFile::Open(fnamev[ifile]); 
    assert(infile);
  
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info",       &info);  TBranch *infoBr   = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("GenEvtInfo", &gen);   TBranch *genBr	 = eventTree->GetBranch("GenEvtInfo");
    eventTree->SetBranchAddress("GenParticle", &genPartArr); TBranch *genPartBr	 = eventTree->GetBranch("GenParticle");
    eventTree->SetBranchAddress("Photon", &photonArr); TBranch *photonBr = eventTree->GetBranch("Photon");   

    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelBBv.push_back(0);
    nSelBEv.push_back(0);
    nSelEEv.push_back(0);
    
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);
      genPartArr->Clear(); genPartBr->GetEntry(ientry);


      if (fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      TLorentzVector *vec=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,1);
      if(vec->M()<MASS_LOW || vec->M()>MASS_HIGH) continue;

      Int_t fid_1 = gen->id_1;
      Int_t fid_2 = gen->id_2;
      Double_t weight=gen->weight*testPdf->xfxQ(fid_1, gen->x_1, gen->scalePDF)*testPdf->xfxQ(fid_2, gen->x_2, gen->scalePDF)/(nomPdf->xfxQ(fid_1, gen->x_1, gen->scalePDF)*nomPdf->xfxQ(fid_2, gen->x_2, gen->scalePDF));
      nEvtsv[ifile]+=weight;

      photonArr->Clear(); photonBr->GetEntry(ientry);
      const baconhep::TPhoton *ele1=0, *ele2=0;
      Double_t drMin1=0.5, drMin2=0.5;
      for(Int_t i=0; i<photonArr->GetEntriesFast(); i++) {
  	const baconhep::TPhoton *sc = (baconhep::TPhoton*)((*photonArr)[i]);
	
	// check ECAL gap
	if(fabs(sc->scEta)>=ETA_BARREL && fabs(sc->scEta)<=ETA_ENDCAP) continue;
        
	if(sc->pt	   < PT_CUT)  continue;  // Supercluster ET cut corrected for PV
        if(fabs(sc->scEta) > ETA_CUT) continue;  // Supercluster |eta| cut
	  
	Double_t decRho = sqrt((info->pvx)*(info->pvx) + (info->pvy)*(info->pvy));
	Double_t ecaleta1 = ecalEta(lep1->Eta(),info->pvz,decRho);
	Double_t ecaleta2 = ecalEta(lep2->Eta(),info->pvz,decRho);
	Double_t dr1 = toolbox::deltaR(sc->scEta,sc->scPhi,ecaleta1,lep1->Phi());
	Double_t dr2 = toolbox::deltaR(sc->scEta,sc->scPhi,ecaleta2,lep2->Phi());
	
	if(dr1<drMin1) {
	  drMin1 = dr1;
	  ele1 = sc;
	}
	
	if(dr2<drMin2) {
	  drMin2 = dr2;
	  ele2 = sc;
	}
      }
      
      if(!ele1 || !ele2) continue;
      
      // mass window
      TLorentzVector vEle1; vEle1.SetPtEtaPhiM(ele1->pt, ele1->scEta, ele1->scPhi, ELE_MASS);
      TLorentzVector vEle2; vEle2.SetPtEtaPhiM(ele2->pt, ele2->scEta, ele2->scPhi, ELE_MASS);
      TLorentzVector vDilep = vEle1 + vEle2;
      if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue; 
      
      /******** We have a Z candidate! HURRAY! ********/
      //cout << ele1->scEta << ", " << ele2->scEta << endl;
      Bool_t isB1 = (fabs(ele1->scEta)<ETA_BARREL) ? kTRUE : kFALSE;
      Bool_t isB2 = (fabs(ele2->scEta)<ETA_BARREL) ? kTRUE : kFALSE;
      
      nSelv[ifile]+=weight;

      if(isB1 && isB2)        nSelBBv[ifile]+=weight;
      else if(!isB1 && !isB2) nSelEEv[ifile]+=weight;
      else		      nSelBEv[ifile]+=weight;	  

    }
    
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);     accErrv.push_back(sqrt(accv[ifile]*(1.-accv[ifile])/nEvtsv[ifile]));
    accBBv.push_back(nSelBBv[ifile]/nEvtsv[ifile]); accErrBBv.push_back(sqrt(accBBv[ifile]*(1.-accBBv[ifile])/nEvtsv[ifile]));
    accBEv.push_back(nSelBEv[ifile]/nEvtsv[ifile]); accErrBEv.push_back(sqrt(accBEv[ifile]*(1.-accBEv[ifile])/nEvtsv[ifile]));
    accEEv.push_back(nSelEEv[ifile]/nEvtsv[ifile]); accErrEEv.push_back(sqrt(accEEv[ifile]*(1.-accEEv[ifile])/nEvtsv[ifile]));
    
    delete infile;
    infile=0, eventTree=0;  
  }  
  delete info;
  delete gen;
  delete photonArr;
  
    
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
  cout << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  cout << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  cout << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "   ================================================" << endl;
    cout << "    Label: " << labelv[ifile] << endl;
    cout << "     File: " << fnamev[ifile] << endl;
    cout << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "     barrel-barrel: " << setw(12) << nSelBBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBBv[ifile] << " +/- " << accErrBBv[ifile] << endl;
    cout << "     barrel-endcap: " << setw(12) << nSelBEv[ifile] << " / " << nEvtsv[ifile] << " = " << accBEv[ifile] << " +/- " << accErrBEv[ifile] << endl;
    cout << "     endcap-endcap: " << setw(12) << nSelEEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEEv[ifile] << " +/- " << accErrEEv[ifile] << endl;
    cout << "             total: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    cout << endl;
  }
  
  char txtfname[100];
  sprintf(txtfname,"%s/sc.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << " Z -> e e" << endl;
  txtfile << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  txtfile << "  pT > " << PT_CUT << endl;
  txtfile << "  |eta| < " << ETA_CUT << endl;
  txtfile << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  txtfile << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  txtfile << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "   ================================================" << endl;
    txtfile << "    Label: " << labelv[ifile] << endl;
    txtfile << "     File: " << fnamev[ifile] << endl;
    txtfile << endl;
    txtfile << "    *** Acceptance ***" << endl;
    txtfile << "     barrel-barrel: " << setw(12) << nSelBBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBBv[ifile] << " +/- " << accErrBBv[ifile] << endl;
    txtfile << "     barrel-endcap: " << setw(12) << nSelBEv[ifile] << " / " << nEvtsv[ifile] << " = " << accBEv[ifile] << " +/- " << accErrBEv[ifile] << endl;
    txtfile << "     endcap-endcap: " << setw(12) << nSelEEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEEv[ifile] << " +/- " << accErrEEv[ifile] << endl;
    txtfile << "             total: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    txtfile << endl;
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSCZee"); 
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
Double_t ecalEta(const Double_t eta, const Double_t z, const Double_t rho) {
  // From https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideEcalRecoClustering
  
  const Double_t R_ECAL            = 136.5;
  const Double_t Z_ENDCAP          = 328.0;
  const Double_t ETA_BARREL_ENDCAP = 1.479;
  const Double_t pi                = 3.14159265358979;

  if(eta!= 0.) {
    Double_t theta = 0.0;
    Double_t zecal = (R_ECAL - rho)*sinh(eta) + z;
      
    if(zecal != 0.0) theta = atan(R_ECAL/zecal);
    if(theta<0.0) theta = theta + pi;

    Double_t detEta = -log(tan(0.5*theta));
      
    if(fabs(detEta) > ETA_BARREL_ENDCAP) {
      Double_t zend = Z_ENDCAP;
      if(eta<0.0) zend = -zend;
      Double_t zlen = zend - z;
      Double_t rr = zlen/sinh(eta);
      theta = atan((rr+rho)/zend);
      if(theta<0.0) theta = theta + pi;
      detEta = - log(tan(0.5*theta));
    }
    return detEta;
      
  } else {
    return eta;
  }
}
