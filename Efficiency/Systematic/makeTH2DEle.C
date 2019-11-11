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
// const int NBeta = 8;
// const float ptrange[NBpt+1]   = {25., 40., 8000.};
// const float etarange[NBeta+1] = {-2.4, -2.1, -1.2, -0.9, 0., 0.9, 1.2, 2.1, 2.4};

const int NBeta = 12;
const float etarange[NBeta+1] = {-2.4,-2.0,-1.566,-1.4442,-1.0,-0.5,0,0.5,1.0,1.4442,1.566,2.0,2.4};
const int NBpt  = 3;
const float ptrange[NBpt+1]   = {25., 35.0, 50., 10000.};

// const vector<TString> charges{"Negative","Positive"};
const vector<TString> charges{"Combined","Combined"};

// TString masterDir="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/testReweights_v2_2/results/";
TString masterDir="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results/TOYS/";
TString DirFSR="_POWxPythia_POWxPhotos/"; // change this once i finish all my shit
TString DirMC="_aMCxPythia_minloxPythia/";
TString bkgDir="_aMCxPythia_POWBKG/"; // should be exp vs Powerlaw


// TString tagDir   ="_aMCxPythia_aMCxPythia_tagPt/"; // should be exp vs Powerlaw
TString effDirD = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results/Zee/Data/";
TString effDirM = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results/Zee/MC/";
TString amc="_aMCxPythia/";
TString tagpt="_aMCxPythia_tagPt/"; // should be exp vs Powerlaw

// TString DirFSR="_aMCxPythia_v2_aMCxPythia_v2/"; // change this once i finish all my shit
// TString DirMC="_aMCxPythia_v2_aMCxPythia_v2/";
// TString bkgDir="_aMCxPythia_v2_aMCxPythia_v2/"; // should be exp vs Powerlaw

TString outDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics";
TString subf="";

void makeTH2DEle(TString eType = "EleGSFSelEff"){
    // alternate input option is HLTEff
    vector<double> vMCNeg;
    vector<double> vMCPos;
    vector<double> vFSRNeg;
    vector<double> vFSRPos;
    vector<double> vBkgNeg;
    vector<double> vBkgPos;
     std::cout << "hello" << std::endl;
     
     
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
     
     
  for(int i = 0; i < NBeta*NBpt+1; ++i){
    // double value=0;
    // char infilename[250];
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),eType.Data(),DirFSR.Data(),charges[0].Data(),i);
    // std::cout << infilename << std::endl;
    
    // std::cout << "what" << std::endl;
    ifstream infile1(infilename); assert (infile1);  infile1>>value; vFSRNeg.push_back(value); value=0;
    std::cout << infilename << vFSRNeg.back() << std::endl;
    // std::cout << "blah" << std::endl;
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),eType.Data(),DirFSR.Data(),charges[1].Data(),i);
    ifstream infile2(infilename); infile2>>value; vFSRPos.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),eType.Data(),DirMC.Data(),charges[0].Data(),i);
    ifstream infile3(infilename); infile3>>value; vMCNeg.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),eType.Data(),DirMC.Data(),charges[1].Data(),i);
    ifstream infile4(infilename); infile4>>value; vMCPos.push_back(value); value=0;    
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),eType.Data(),bkgDir.Data(),charges[0].Data(),i);
    ifstream infile5(infilename); infile5>>value; vBkgNeg.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),eType.Data(),bkgDir.Data(),charges[1].Data(),i);
    ifstream infile6(infilename); infile6>>value; vBkgPos.push_back(value); value=0; 
  }
  std::cout << "done loop" << std::endl;
    
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
			hMCNeg->SetBinContent(ieta+1, ipt+1, 1.+vMCNeg[ipt*NBeta+ieta]);
			hMCPos->SetBinContent(ieta+1, ipt+1, 1.+vMCPos[ipt*NBeta+ieta]);
			hBkgNeg->SetBinContent(ieta+1, ipt+1, 1.+vBkgNeg[ipt*NBeta+ieta]);
			hBkgPos->SetBinContent(ieta+1, ipt+1, 1.+vBkgPos[ipt*NBeta+ieta]);
      
      double num = tagPosDat->GetBinContent(ieta+1, ipt+1) / amcPosDat->GetBinContent(ieta+1, ipt+1);
      double dnm = tagPosMC ->GetBinContent(ieta+1, ipt+1) / amcPosMC ->GetBinContent(ieta+1, ipt+1);
      hTagNeg->SetBinContent(ieta+1, ipt+1, num/dnm);
      
      num = tagNegDat->GetBinContent(ieta+1, ipt+1) / amcNegDat->GetBinContent(ieta+1, ipt+1);
      dnm = tagNegMC ->GetBinContent(ieta+1, ipt+1) / amcNegMC ->GetBinContent(ieta+1, ipt+1);
			hTagPos->SetBinContent(ieta+1, ipt+1, num/dnm);
      // double num = 1;
      // double dnm = 1;
      // hTagNeg->SetBinContent(ieta+1, ipt+1, num/dnm);
      
      // num = 1;
      // dnm = 1;
			// hTagPos->SetBinContent(ieta+1, ipt+1, num/dnm);
      
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
