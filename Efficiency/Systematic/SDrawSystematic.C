#include <iostream>
#include <TLegend.h>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include "TH1D.h"
#include "TH2D.h"
#include <THStack.h>
#include "TProfile.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFractionFitter.h"
#include <string>
#include <vector>
#include <math.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPave.h>
#include <TPaveStats.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <fstream>
#include <TString.h>

using namespace std;

string int2string(int i){
  stringstream ss;
  string ret;
  ss<<i;
  ss>>ret;
  return ret;
}

void SDrawSystematic(const TString inputFile,
		     const TString outputDir,
		     const TString outputName)
{
	char inputfile[100];
	sprintf(inputfile,"%s",inputFile.Data());

	TTree *t = new TTree("t","t");
	t->ReadFile(inputfile,"Pull");
	t->SetEstimate(t->GetEntries());

//	TH1D *h = new TH1D("h", "", 100, t->GetMinimum("Pull"), t->GetMaximum("Pull"));
	TH1D *h = new TH1D("h", "", 100, -15., 15.);
//	TH1D *h = new TH1D("h", "", 100, -50., 50.);
//	TH1D *h = new TH1D("h", "", 100, -0.1, 0.1);
//	TH1D *h = new TH1D("h", "", 100, 0., 1.);
	t->Project("h", "Pull");
//	h->Rebin(10);

	h->Fit("gaus");
	float meanhist = abs(h->GetMean());
	float meangaus = abs(h->GetFunction("gaus")->GetParameter(1));
//	float mean = meangaus < meanhist ? meangaus : meanhist; 
	//float mean = meantemp > 0? meantemp : -1.*meantemp;
	float mean = meangaus;
	float sigma = h->GetFunction("gaus")->GetParameter(2);

	char outputfile[100];
	sprintf(outputfile,"%s/%s.txt",outputDir.Data(), outputName.Data());
	ofstream meanfile;
	meanfile.open(outputfile);
	meanfile<<mean<<" "<<sigma<<endl;
	meanfile.close();

	char outputpng[100];
	sprintf(outputpng, "%s/%s.png",outputDir.Data(), outputName.Data());
        TCanvas *c = new TCanvas("c","c");
        gStyle->SetOptFit();
	h->Draw();
	c->SaveAs(outputpng);
	c->Clear();

	delete c;
	delete h;
	delete t;
}
