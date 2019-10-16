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

const int NBpt  = 2;
const float ptrange[NBpt+1]   = {25., 40., 8000.};
const int NBeta = 12;
const float etarange[NBeta+1] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,0,0.3,0.9,1.2,1.6,2.1,2.4};


// const int NBeta = 12;
// const float etarange[NBeta+1] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,0,0.3,0.9,1.2,1.6,2.1,2.4};
// const int NBpt = 8;
// const float ptrange[NBpt+1] = {25,30,35,40,45,50,60,80,8000};

// TString mainDir="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/testReweights_v2_2/results/";
// TString mainDir="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/LowPU2017ID_13TeV_v0/results/TOYS/";
TString mainDir="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results/TOYS_v2/";
TString sigDirFSR="_POWxPythia_POWxPhotos/";
TString sigDirMC ="_aMCxPythia_minloxPythia/";
TString bkgDir   ="_aMCxPythia_POWBKG/"; // should be exp vs Powerlaw
// TString tagDir   ="_aMCxPythia_aMCxPythia_tagPt/"; // should be exp vs Powerlaw
TString effDirD = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results/Zmm/Data/";
TString effDirM = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results/Zmm/MC/";
TString amc="_aMCxPythia/";
TString tagpt="_aMCxPythia_tagPt/"; // should be exp vs Powerlaw

TString outDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics";
TString subf="";

// TString sigDirMC="_aMCxPythia_v0_POWBKG_v0/";
// TString sigDirFSR="_aMCxPythia_v0_POWBKG_v0/";
// TString bkgDir="_aMCxPythia_v0_POWBKG_v0/"; // should be exp vs Powerlaw

// For the closure tests
// TString sigDirMC="_aMCxPythia_v0_aMCxPythia_v0/";
// TString sigDirFSR="_aMCxPythia_v0_aMCxPythia_v0/";
// TString bkgDir="_aMCxPythia_v0_aMCxPythia_v0/"; 

// const string charges[2] = {"Negative","Positive"};
// const string charges[2] = {"Combined","Combined"};

// const vector<TString> charges{"Negative","Positive"};
const vector<TString> charges{"Combined", "Combined"};

// the alternative values for the Standalone 
const float fsrReplace[NBpt*NBeta]={ 0.999918, 1.003591, 1.000001, 1.000711, 1.000439, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000};
const float mcReplace[NBpt*NBeta]={0.999980, 0.992656, 0.999727, 0.996864, 1.001382, 1.003100, 1.004811, 1.000910, 1.002059, 0.998573, 0.999180, 0.995128, 1.000510, 0.999844, 0.999514, 1.000040, 1.000328, 1.001896, 1.001451, 1.001363, 1.000182, 1.000526, 0.999576, 1.001489};
const float bkgReplace[NBpt*NBeta]={0.997122, 0.990444, 0.998531, 0.993576, 0.996215, 0.993380, 0.994157, 0.976469, 0.998377, 0.998529, 0.986246, 0.993661, 0.997157, 0.998264, 0.998694, 0.998115, 0.998386, 0.999483, 0.999141, 0.999100, 0.998073, 0.999231, 0.999083, 0.999723};

void makeTH2DMu(TString eType = "MuStaEff"){
// void makeTH2DMu(TString eType = "MuSITEff"){
    // alternate input option is MuStaEff, HLT should not have fits?MuSITEff
    bool doAbs=false; 
    vector<double> vMCNeg;
    vector<double> vMCPos;
    vector<double> vFSRNeg;
    vector<double> vFSRPos;
    vector<double> vBkgNeg; 
    vector<double> vBkgPos;
    vector<double> vTagNeg; 
    vector<double> vTagPos;
    
    
    double value=0; // read in value is the % diff from the central value
    char infilename[200];
    sprintf(infilename,"%s%s%s%s%s/eff.root",effDirD.Data(),eType.Data(),amc.Data(),subf.Data(),charges[0].Data());
    TFile *amcPd = new TFile(infilename); TH2D *amcPosDat = (TH2D*) amcPd->Get("hEffEtaPt");
    sprintf(infilename,"%s%s%s%s%s/eff.root",effDirD.Data(),eType.Data(),amc.Data(),subf.Data(),charges[1].Data());
    TFile *amcNd = new TFile(infilename); TH2D *amcNegDat = (TH2D*) amcNd->Get("hEffEtaPt");
    
    sprintf(infilename,"%s%s%s%s%s/eff.root",effDirD.Data(),eType.Data(),tagpt.Data(),subf.Data(),charges[0].Data());
    TFile *tagPd = new TFile(infilename);TH2D *tagPosDat = (TH2D*) tagPd->Get("hEffEtaPt");
    sprintf(infilename,"%s%s%s%s%s/eff.root",effDirD.Data(),eType.Data(),tagpt.Data(),subf.Data(),charges[1].Data());
    TFile *tagNd = new TFile(infilename);TH2D *tagNegDat = (TH2D*) tagNd->Get("hEffEtaPt");
    
    sprintf(infilename,"%s%s%s%s%s/eff.root",effDirM.Data(),eType.Data(),amc.Data(),subf.Data(),charges[0].Data());
    TFile *amcPm = new TFile(infilename); TH2D *amcPosMC = (TH2D*) amcPm->Get("hEffEtaPt");
    sprintf(infilename,"%s%s%s%s%s/eff.root",effDirM.Data(),eType.Data(),amc.Data(),subf.Data(),charges[1].Data());
    TFile *amcNm = new TFile(infilename); TH2D *amcNegMC = (TH2D*) amcNm->Get("hEffEtaPt");
    
    sprintf(infilename,"%s%s%s%s%s/eff.root",effDirM.Data(),eType.Data(),tagpt.Data(),subf.Data(),charges[0].Data());
    TFile *tagPm = new TFile(infilename);TH2D *tagPosMC = (TH2D*) tagPm->Get("hEffEtaPt");
    sprintf(infilename,"%s%s%s%s%s/eff.root",effDirM.Data(),eType.Data(),tagpt.Data(),subf.Data(),charges[1].Data());
    TFile *tagNm = new TFile(infilename);TH2D *tagNegMC = (TH2D*) tagNm->Get("hEffEtaPt");
    
  for(int i = 0; i < NBeta*NBpt; ++i){
    
    
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",mainDir.Data(),eType.Data(),sigDirFSR.Data(),subf.Data(),charges[0].Data(),i);
    ifstream infile1(infilename); assert (infile1); infile1>>value; value=(doAbs?fabs(value):value);vFSRNeg.push_back(value); value=0;
    std::cout << infilename << vFSRNeg.back() << std::endl;
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",mainDir.Data(),eType.Data(),sigDirFSR.Data(),subf.Data(),charges[1].Data(),i);
    ifstream infile2(infilename); assert (infile2);infile2>>value; value=(doAbs?fabs(value):value);vFSRPos.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",mainDir.Data(),eType.Data(),sigDirMC.Data(),subf.Data(),charges[0].Data(),i);
    ifstream infile3(infilename); assert (infile3);infile3>>value; value=(doAbs?fabs(value):value);vMCNeg.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",mainDir.Data(),eType.Data(),sigDirMC.Data(),subf.Data(),charges[1].Data(),i);
    ifstream infile4(infilename); assert (infile4);infile4>>value;value=(doAbs?fabs(value):value);vMCPos.push_back(value); value=0;    
    std::cout << infilename << vMCPos.back() << std::endl;
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",mainDir.Data(),eType.Data(),bkgDir.Data(),subf.Data(),charges[0].Data(),i);
    ifstream infile5(infilename); assert (infile5);infile5>>value; value=(doAbs?fabs(value):value);vBkgNeg.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",mainDir.Data(),eType.Data(),bkgDir.Data(),subf.Data(),charges[1].Data(),i);
    ifstream infile6(infilename); infile6>>value; value=(doAbs?fabs(value):value); vBkgPos.push_back(value); value=0; 
    
    if(i==2||i==10||i==14) {vFSRNeg[i]=fsrReplace[i]-1;vFSRPos[i]=fsrReplace[i]-1;}
    // if(i==2||i==10||i==14) {vMCNeg[i]=mcReplace[i]-1;vMCPos[i]=mcReplace[i]-1;}
    // if(i==10||i==14) {vBkgNeg[i]=bkgReplace[i]-1;vBkgPos[i]=bkgReplace[i]-1;}
  }
    
    char histname[100];
    sprintf(histname,"hMCNeg");
    TH2D *hMCNeg = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"hMCPos");
    TH2D *hMCPos = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"hFSRNeg");
    TH2D *hFSRNeg = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"hFSRPos");
    TH2D *hFSRPos = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"hBkgNeg");
    TH2D *hBkgNeg = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"hBkgPos");
    TH2D *hBkgPos = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"hTagNeg");
    TH2D *hTagNeg = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"hTagPos");
    TH2D *hTagPos = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);

	for(int ipt=0; ipt<NBpt; ipt++){

		for(int ieta=0; ieta<NBeta; ieta++){
      
      std::cout << ipt*NBeta+ieta << "  " << 1.+vFSRNeg[ipt*NBeta+ieta] << " " << 1.+vMCNeg[ipt*NBeta+ieta] << std::endl;
			hFSRNeg->SetBinContent(ieta+1, ipt+1, 1.+vFSRNeg[ipt*NBeta+ieta]);
			hFSRPos->SetBinContent(ieta+1, ipt+1, 1.+vFSRPos[ipt*NBeta+ieta]);
			hMCNeg ->SetBinContent(ieta+1, ipt+1, 1.+vMCNeg[ipt*NBeta+ieta]);
			hMCPos ->SetBinContent(ieta+1, ipt+1, 1.+vMCPos[ipt*NBeta+ieta]);
			hBkgNeg->SetBinContent(ieta+1, ipt+1, 1.+vBkgNeg[ipt*NBeta+ieta]);
			hBkgPos->SetBinContent(ieta+1, ipt+1, 1.+vBkgPos[ipt*NBeta+ieta]);
      
      // hSigFSRNeg->SetBinContent(ieta+1, ipt+1, 1.-vFSRNeg[ipt*NBeta+ieta]);
			// hSigFSRPos->SetBinContent(ieta+1, ipt+1, 1.-vFSRPos[ipt*NBeta+ieta]);
			// hSigMCNeg->SetBinContent(ieta+1, ipt+1, 1.-vMCNeg[ipt*NBeta+ieta]);
			// hSigMCPos->SetBinContent(ieta+1, ipt+1, 1.-vMCPos[ipt*NBeta+ieta]);
			// hBkgNeg->SetBinContent(ieta+1, ipt+1, 1.-vBkgNeg[ipt*NBeta+ieta]);
			// hBkgPos->SetBinContent(ieta+1, ipt+1, 1.-vBkgPos[ipt*NBeta+ieta]);
      
      double num = tagPosDat->GetBinContent(ieta+1, ipt+1) / amcPosDat->GetBinContent(ieta+1, ipt+1);
      double dnm = tagPosMC ->GetBinContent(ieta+1, ipt+1) / amcPosMC ->GetBinContent(ieta+1, ipt+1);
      hTagNeg->SetBinContent(ieta+1, ipt+1, num/dnm);
      
      num = tagNegDat->GetBinContent(ieta+1, ipt+1) / amcNegDat->GetBinContent(ieta+1, ipt+1);
      dnm = tagNegMC ->GetBinContent(ieta+1, ipt+1) / amcNegMC ->GetBinContent(ieta+1, ipt+1);
			hTagPos->SetBinContent(ieta+1, ipt+1, num/dnm);
		}
	}
  char outfilename[500];
  sprintf(outfilename,"%s/SysUnc_%s.root",outDir.Data(),eType.Data());
	TFile *f = new TFile(outfilename, "RECREATE");
	hFSRNeg->Write();
	hFSRPos->Write();
	hMCNeg->Write();
	hMCPos->Write();
	hBkgNeg->Write();
	hBkgPos->Write();
  hTagNeg->Write();
	hTagPos->Write();
	f->Close();
}
