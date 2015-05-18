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
  CorrPlot(TString name, TString title, TString xtitle, TString ytitle);
  ~CorrPlot(){}
  
  void AddCorrPlot(TGraph *gr, TEllipse *el, TString label, int color=kBlack, int marksty=kFullDotLarge, int linesty=1);
  void Draw(TCanvas *c, TString fname);
  
protected:
  TString fName;
  TString fTitle;
  TString fXtitle;
  TString fYtitle;
  TLegend *fLeg; 	     
  
  vector<CorrPlotItem> fItems;
};

#endif
