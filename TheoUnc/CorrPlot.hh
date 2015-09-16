#ifndef CORRPLOT_HH
#define CORRPLOT_HH

#include <TROOT.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TEllipse.h>
#include <TPaveText.h>
#include <vector>

using namespace std;

class CorrPlotItem {
public:
  CorrPlotItem():graph(0),ellipse(0){}
  ~CorrPlotItem(){}

  TGraph* graph;
  TEllipse* ellipse;

};

class CorrPlot {
public:
  CorrPlot();
  CorrPlot(TString name, TString title, TString xtitle, TString ytitle, Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax);
  ~CorrPlot(){}
  
  void AddCorrPlot(TGraph *gr, TEllipse *el, TString label, int color=kBlack, int marksty=kFullDotLarge, int linesty=1, int fillsty=0);
  void Draw(TCanvas *c, TString fname, Int_t lumi);
  
protected:
  TString fName;
  TString fTitle;
  TString fXtitle;
  TString fYtitle;
  Float_t fXmin;
  Float_t fXmax;
  Float_t fYmin;
  Float_t fYmax;
  TLegend *fLeg=0;
  
  vector<CorrPlotItem> fItems;
};

#endif
