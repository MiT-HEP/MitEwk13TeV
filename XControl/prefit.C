#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "RooMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

Double_t myfunction(Double_t *x, Double_t *par){

	Double_t erf = RooMath::erfc((par[0] - x[0])*par[1]);
	Double_t u = (x[0] - 91.1876)*par[2];

	if(u < -70) u = 1e20;
	else if( u>70 ) u = 0;
	else u = exp(-u);   //exponential decay
	return erf*u;
} 

void prefit(){

	TFile *f = TFile::Open("/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/MitEwk13TeV/XControl/f_bkgfail.root"); 
	TFile *fbkg = new TFile("fbkg.root","RECREATE");
	TH1D *h[24];
	TH1D *hbkg[24];

//	TF1 *fq = new TF1("fq","[0]+[1]*x+[2]*x*x",60.,120.);
//	TF1 *fl = new TF1("fl","[0]+[1]*x",60.,120.);
	
//	TF1 *fq = new TF1("fq",myfunction,60.,120.,3);

	TF1 *fq = new TF1("fq","x**(-1.*[0])",60.,120.);

	for(int i = 0; i < 24; i++){
		char name[20];
		char namebkg[20];
		sprintf(name, "histFail_%d", i);
		sprintf(namebkg, "histbkgFail_%d", i);

		h[i] = (TH1D*) f->Get(name);
		hbkg[i] = new TH1D(namebkg, "", 60, 60., 120.);

		for(int ib = 0; ib < 20; ib++){
		        hbkg[i]->SetBinContent(ib+1, h[i]->GetBinContent(ib+1));
		        hbkg[i]->SetBinError(ib+1, h[i]->GetBinError(ib+1));
		}

                for(int ib = 40; ib < 60; ib++){
                        hbkg[i]->SetBinContent(ib+1, h[i]->GetBinContent(ib+1));
                        hbkg[i]->SetBinError(ib+1, h[i]->GetBinError(ib+1));
                }
		//hbkg[i]->Write();
		hbkg[i]->Fit("fq");
		std::cout<<"ibin = "<<i<<": "<<fq->GetParameter(0)<<","<<fq->GetParError(0)<<" ,"<<fq->GetParameter(1)<<","<<fq->GetParError(1)<<","<<fq->GetParameter(2)<<","<<fq->GetParError(2)<<std::endl;
//		if(i==10){
//		hbkg[i]->Fit("fl");
//		std::cout<<"ibin = "<<i<<": "<<fl->GetParameter(0)<<","<<fl->GetParError(0)<<","<<fl->GetParameter(1)<<","<<fl->GetParError(1)<<std::endl;
//		}
	}

	fbkg->Write();
	fbkg->Close();
	f->Close();
}
