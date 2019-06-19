#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <fstream>

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMinuit.h"

// const int NBpt  = 2;
// const float ptrange[NBpt+1]   = {25., 40., 8000.};
// const int NBeta = 12;
// const float etarange[NBeta+1] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,0,0.3,0.9,1.2,1.6,2.1,2.4};


const int NBeta = 12;
const float etarange[NBeta+1] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,0,0.3,0.9,1.2,1.6,2.1,2.4};
const int NBpt = 8;
const float ptrange[NBpt+1] = {25,30,35,40,45,50,60,80,8000};

// TString masterDir="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/testReweights_v2_2/results/";
// TString masterDir="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/LowPU2017ID_13TeV_v0/results/TOYS/";
TString masterDir="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/LowPU2017ID_13TeV_v1/results/Zmm/Data/";
TString amcnlo="_aMCxPythia_v2/";
TString minlo="_minloxPythia_v2/";
// TString bkgDir="_aMCxPythia_v1_POWBKG_v1/"; // should be exp vs Powerlaw
TString pythia="_POWxPythia_v2/"; // should be exp vs Powerlaw
TString photos="_POWxPhotos_v2/"; // should be exp vs Powerlaw
TString powerlaw="_POWBKG_v2/"; // should be exp vs Powerlaw

TString subfolder="";

// TString sigDirMC="_aMCxPythia_v0_POWBKG_v0/";
// TString sigDirFSR="_aMCxPythia_v0_POWBKG_v0/";
// TString bkgDir="_aMCxPythia_v0_POWBKG_v0/"; // should be exp vs Powerlaw

// For the closure tests
// TString sigDirMC="_aMCxPythia_v0_aMCxPythia_v0/";
// TString sigDirFSR="_aMCxPythia_v0_aMCxPythia_v0/";
// TString bkgDir="_aMCxPythia_v0_aMCxPythia_v0/"; 

const string charges[2] = {"Negative","Positive"};
// const string charges[2] = {"Combined","Combined"};

// void makeTH2DMuDirect(TString effTypeSig = "MuStaEff"){
void makeTH2DMuDirect(TString effTypeSig = "MuSITEff"){
    // alternate input option is MuStaEff, HLT should not have fits?MuSITEff
    bool doAbs=false; 
    vector<double> vSigMCNeg;
    vector<double> vSigMCPos;
    vector<double> vSigFSRNeg;
    vector<double> vSigFSRPos;
    vector<double> vBkgNeg; 
    vector<double> vBkgPos;
    
    char infilename[200];
    sprintf(infilename,"%s%s%s%s%s/eff.root",masterDir.Data(),effTypeSig.Data(),amcnlo.Data(),subfolder.Data(),charges[0].c_str());
    TFile *amcP = new TFile(infilename); TH2D *amcnloPos = (TH2D*) amcP->Get("hEffEtaPt");
    sprintf(infilename,"%s%s%s%s%s/eff.root",masterDir.Data(),effTypeSig.Data(),amcnlo.Data(),subfolder.Data(),charges[1].c_str());
    TFile *amcN = new TFile(infilename); TH2D *amcnloNeg = (TH2D*) amcN->Get("hEffEtaPt");
    
    sprintf(infilename,"%s%s%s%s%s/eff.root",masterDir.Data(),effTypeSig.Data(),minlo.Data(),subfolder.Data(),charges[0].c_str());
    TFile *minP = new TFile(infilename);TH2D *minloPos = (TH2D*) minP->Get("hEffEtaPt");
    sprintf(infilename,"%s%s%s%s%s/eff.root",masterDir.Data(),effTypeSig.Data(),minlo.Data(),subfolder.Data(),charges[1].c_str());
    TFile *minN = new TFile(infilename);TH2D *minloNeg = (TH2D*) minN->Get("hEffEtaPt");
    
    sprintf(infilename,"%s%s%s%s%s/eff.root",masterDir.Data(),effTypeSig.Data(),pythia.Data(),subfolder.Data(),charges[0].c_str());
    TFile *pytP = new TFile(infilename);TH2D *pythiaPos = (TH2D*) pytP->Get("hEffEtaPt");
    sprintf(infilename,"%s%s%s%s%s/eff.root",masterDir.Data(),effTypeSig.Data(),pythia.Data(),subfolder.Data(),charges[1].c_str());
    TFile *pytN = new TFile(infilename);TH2D *pythiaNeg = (TH2D*) pytN->Get("hEffEtaPt");
    
    sprintf(infilename,"%s%s%s%s%s/eff.root",masterDir.Data(),effTypeSig.Data(),photos.Data(),subfolder.Data(),charges[0].c_str());
    TFile *photP = new TFile(infilename);TH2D *photosPos = (TH2D*) photP->Get("hEffEtaPt");
    sprintf(infilename,"%s%s%s%s%s/eff.root",masterDir.Data(),effTypeSig.Data(),photos.Data(),subfolder.Data(),charges[1].c_str());
    TFile *photN = new TFile(infilename);TH2D *photosNeg = (TH2D*) photN->Get("hEffEtaPt");
    
    sprintf(infilename,"%s%s%s%s%s/eff.root",masterDir.Data(),effTypeSig.Data(),powerlaw.Data(),subfolder.Data(),charges[0].c_str());
    TFile *plP = new TFile(infilename);TH2D *plPos = (TH2D*) plP->Get("hEffEtaPt");
    sprintf(infilename,"%s%s%s%s%s/eff.root",masterDir.Data(),effTypeSig.Data(),powerlaw.Data(),subfolder.Data(),charges[1].c_str());
    TFile *plN = new TFile(infilename);TH2D *plNeg = (TH2D*) plN->Get("hEffEtaPt");
    
    char histname[100];
    sprintf(histname,"h%sSigMCNeg",effTypeSig.Data());
	TH2D *hSigMCNeg = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"h%sSigMCPos",effTypeSig.Data());
	TH2D *hSigMCPos = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"h%sSigFSRNeg",effTypeSig.Data());
    TH2D *hSigFSRNeg = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"h%sSigFSRPos",effTypeSig.Data());
    TH2D *hSigFSRPos = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"h%sBkgNeg",effTypeSig.Data());
    TH2D *hBkgNeg = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"h%sBkgPos",effTypeSig.Data());
    TH2D *hBkgPos = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    

	for(int ipt=0; ipt<NBpt; ipt++){

		for(int ieta=0; ieta<NBeta; ieta++){
            
			hSigFSRNeg->SetBinContent(ieta+1, ipt+1, photosNeg->GetBinContent(ieta+1,ipt+1)/pythiaNeg->GetBinContent(ieta+1,ipt+1));
			hSigFSRPos->SetBinContent(ieta+1, ipt+1, photosPos->GetBinContent(ieta+1,ipt+1)/pythiaPos->GetBinContent(ieta+1,ipt+1));
			hSigMCNeg->SetBinContent(ieta+1, ipt+1, minloNeg->GetBinContent(ieta+1,ipt+1)/amcnloNeg->GetBinContent(ieta+1,ipt+1));
			hSigMCPos->SetBinContent(ieta+1, ipt+1, minloPos->GetBinContent(ieta+1,ipt+1)/amcnloPos->GetBinContent(ieta+1,ipt+1));
			hBkgNeg->SetBinContent(ieta+1, ipt+1, plNeg->GetBinContent(ieta+1,ipt+1)/amcnloNeg->GetBinContent(ieta+1,ipt+1));
			hBkgPos->SetBinContent(ieta+1, ipt+1, plPos->GetBinContent(ieta+1,ipt+1)/amcnloNeg->GetBinContent(ieta+1,ipt+1));
		}
	}
    char outfilename[100];
    sprintf(outfilename,"SysUnc_%s.root",effTypeSig.Data());
	TFile *f = new TFile(outfilename, "RECREATE");
	hSigFSRNeg->Write();
	hSigFSRPos->Write();
	hSigMCNeg->Write();
	hSigMCPos->Write();
	hBkgNeg->Write();
	hBkgPos->Write();
	f->Close();
}
