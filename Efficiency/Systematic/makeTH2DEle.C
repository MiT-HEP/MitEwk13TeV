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
const int NBpt = 8;
const float ptrange[NBpt+1] = {25,30,35,40,45,50,60,80,8000};

TString masterDir="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/testReweights_v2_2/results/";
// TString masterDir="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/LowPU2017ID_13TeV_v0//results/";
TString sigDirMC="_POWxPythia_v1_POWxPhotos_v1/"; // change this once i finish all my shit
TString sigDirFSR="_aMCxPythia_v1_minloxPythia_v1/";
TString bkgDir="_aMCxPythia_v1_POWBKG_v1/"; // should be exp vs Powerlaw

void makeTH2DEle(TString effTypeSig = "GSFSelEff"){
    // alternate input option is HLTEff
    vector<double> vSigMCNeg;
    vector<double> vSigMCPos;
    vector<double> vSigFSRNeg;
    vector<double> vSigFSRPos;
    vector<double> vBkgNeg;
    vector<double> vBkgPos;
    
  for(int i = 0; i < NBeta*NBpt+1; ++i){
    double value=0;
    char infilename[200];
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),sigDirFSR.Data(),"Negative/_v2",i);
    // std::cout << infilename << std::endl;
    ifstream infile1(infilename); infile1>>value; vSigFSRNeg.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),sigDirFSR.Data(),"Positive/_v2",i);
    ifstream infile2(infilename); infile2>>value; vSigFSRPos.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),sigDirMC.Data(),"Negative/_v2",i);
    ifstream infile3(infilename); infile3>>value; vSigMCNeg.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),sigDirMC.Data(),"Positive/_v2",i);
    ifstream infile4(infilename); infile4>>value; vSigMCPos.push_back(value); value=0;    
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),bkgDir.Data(),"Negative/_v2",i);
    ifstream infile5(infilename); infile5>>value; vBkgNeg.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),bkgDir.Data(),"Positive/_v2",i);
    ifstream infile6(infilename); infile6>>value; vBkgPos.push_back(value); value=0; 
  }
    
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
			hSigFSRNeg->SetBinContent(ieta+1, ipt+1, 1.+vSigFSRNeg[ipt*NBeta+ieta]);
			hSigFSRPos->SetBinContent(ieta+1, ipt+1, 1.+vSigFSRPos[ipt*NBeta+ieta]);
			hSigMCNeg->SetBinContent(ieta+1, ipt+1, 1.+vSigMCNeg[ipt*NBeta+ieta]);
			hSigMCPos->SetBinContent(ieta+1, ipt+1, 1.+vSigMCPos[ipt*NBeta+ieta]);
			hBkgNeg->SetBinContent(ieta+1, ipt+1, 1.+vBkgNeg[ipt*NBeta+ieta]);
			hBkgPos->SetBinContent(ieta+1, ipt+1, 1.+vBkgPos[ipt*NBeta+ieta]);
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
