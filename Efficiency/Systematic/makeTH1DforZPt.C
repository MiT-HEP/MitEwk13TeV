#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMinuit.h"

const int muEtaNB = 28;
const int elEtaNB = 29;
const float muEtaRange[muEtaNB+1] = {-2.4,-2.3,-2.2,-2.1,-2.0,-1.7,-1.6,-1.5,-1.4,-1.2,-0.8,-0.5,-0.3,-0.2,0.0,0.2,0.3,0.5,0.8,1.2,1.4,1.5,1.6,1.7,2.0,2.1,2.2,2.3,2.4};
const float elEtaRange[elEtaNB+1] = {-2.4,-2.3,-2.2,-2.1,-2.0,-1.7,-1.6,-1.566,-1.4442,-1.4,-1.2,-0.8,-0.5,-0.3,-0.2,0.0,0.2,0.3,0.5,0.8,1.2,1.4442,1.566,1.6,1.7,2.0,2.1,2.2,2.3,2.4};

void makeTH1DforZPt(
){

	TH1D *h_muEffSigSys = new TH1D("h_muEffSigSys","h_muEffSigSys", muEtaNB, muEtaRange);
	TH1D *h_muEffBkgSys = new TH1D("h_muEffBkgSys","h_muEffBkgSys", muEtaNB, muEtaRange);
	TH1D *h_elEffSigSys = new TH1D("h_elEffSigSys","h_elEffSigSys", elEtaNB, elEtaRange);
	TH1D *h_elEffBkgSys = new TH1D("h_elEffBkgSys","h_elEffBkgSys", elEtaNB, elEtaRange);

	for(int ieta=0; ieta<muEtaNB; ieta++){
		char inputfilesig[100];
		sprintf(inputfilesig,"/afs/cern.ch/user/x/xniu/public/Systematics/Results/Sig_Mu_Trk_etapt_%d.txt", ieta);
		std::ifstream infilesig(inputfilesig);
		float syssig = 0;
		while (infilesig >> syssig)
		h_muEffSigSys->SetBinContent(ieta+1, syssig);

                char inputfilebkg[100];
                sprintf(inputfilebkg,"/afs/cern.ch/user/x/xniu/public/Systematics/Results/Bkg_Mu_Trk_etapt_%d.txt", ieta);
                std::ifstream infilebkg(inputfilebkg);
                float sysbkg = 0;
                while (infilebkg >> sysbkg)
                h_muEffBkgSys->SetBinContent(ieta+1, sysbkg);
	}

        for(int ieta=0; ieta<elEtaNB; ieta++){
                char inputfilesig[100];
                sprintf(inputfilesig,"/afs/cern.ch/user/x/xniu/public/Systematics/Results/Sig_Ele_Gsf_etapt_%d.txt", ieta);
                std::ifstream infilesig(inputfilesig);
                float syssig = 0;
                while (infilesig >> syssig)
                h_elEffSigSys->SetBinContent(ieta+1, syssig);

                char inputfilebkg[100];
                sprintf(inputfilebkg,"/afs/cern.ch/user/x/xniu/public/Systematics/Results/Bkg_Ele_Gsf_etapt_%d.txt", ieta);
                std::ifstream infilebkg(inputfilebkg);
                float sysbkg = 0;
                while (infilebkg >> sysbkg)
                h_elEffBkgSys->SetBinContent(ieta+1, sysbkg);
        }

	TCanvas *c = new TCanvas("c","c");
        h_muEffSigSys->SetTitle("");
        h_muEffSigSys->SetStats(0);
        h_muEffSigSys->GetXaxis()->SetTitle("#eta");
	h_muEffSigSys->GetYaxis()->SetTitle("systematics");
        h_muEffSigSys->GetXaxis()->SetTitleSize(0.05);
        h_muEffSigSys->GetYaxis()->SetTitleSize(0.05);
        h_muEffSigSys->SetMarkerStyle(20);
        h_muEffSigSys->SetLineColor(1);
        h_muEffSigSys->SetLineWidth(2);
        h_muEffSigSys->Draw();
	c->SaveAs("muEffSigSys.png");

        h_muEffBkgSys->SetTitle("");
        h_muEffBkgSys->SetStats(0);
        h_muEffBkgSys->GetXaxis()->SetTitle("#eta");
        h_muEffBkgSys->GetYaxis()->SetTitle("systematics");
        h_muEffBkgSys->GetXaxis()->SetTitleSize(0.05);
        h_muEffBkgSys->GetYaxis()->SetTitleSize(0.05);
        h_muEffBkgSys->SetMarkerStyle(20);
        h_muEffBkgSys->SetLineColor(1);
        h_muEffBkgSys->SetLineWidth(2);
        h_muEffBkgSys->Draw();
        c->SaveAs("muEffBkgSys.png");

        h_elEffSigSys->SetTitle("");
        h_elEffSigSys->SetStats(0);
        h_elEffSigSys->GetXaxis()->SetTitle("#eta");
        h_elEffSigSys->GetYaxis()->SetTitle("systematics");
        h_elEffSigSys->GetXaxis()->SetTitleSize(0.05);
        h_elEffSigSys->GetYaxis()->SetTitleSize(0.05);
        h_elEffSigSys->SetMarkerStyle(20);
        h_elEffSigSys->SetLineColor(1);
        h_elEffSigSys->SetLineWidth(2);
        h_elEffSigSys->Draw();
        c->SaveAs("elEffSigSys.png");

        h_elEffBkgSys->SetTitle("");
        h_elEffBkgSys->SetStats(0);
        h_elEffBkgSys->GetXaxis()->SetTitle("#eta");
        h_elEffBkgSys->GetYaxis()->SetTitle("systematics");
        h_elEffBkgSys->GetXaxis()->SetTitleSize(0.05);
        h_elEffBkgSys->GetYaxis()->SetTitleSize(0.05);
        h_elEffBkgSys->SetMarkerStyle(20);
        h_elEffBkgSys->SetLineColor(1);
        h_elEffBkgSys->SetLineWidth(2);
        h_elEffBkgSys->Draw();
        c->SaveAs("elEffBkgSys.png");

	TFile *f = new TFile("LeptonRecoEffUnc.root", "RECREATE");
        h_muEffSigSys->Write();
        h_muEffBkgSys->Write();
        h_elEffSigSys->Write();
        h_elEffBkgSys->Write();
	f->Close();
}
