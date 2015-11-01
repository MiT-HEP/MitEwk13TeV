#include <CorrPlot.hh>
#include <iostream>

CorrPlot::CorrPlot()
{
  CorrPlot("wtf", "is", "going", "on",0,10,0,10);
}

CorrPlot::CorrPlot(TString name, TString title, TString xtitle, TString ytitle, Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax):
  fName(name),
  fTitle(title),
  fXtitle(xtitle),
  fYtitle(ytitle),
  fXmin(xmin),
  fXmax(xmax),
  fYmin(ymin),
  fYmax(ymax)
{}

void CorrPlot::AddCorrPlot(TGraph *gr, TEllipse *el, TString label, int color, int marksty, int linesty, int fillsty) {
  if (!gr || !el) return;
  if (fLeg==0) {
    fLeg = new TLegend(0.6, 0.2, 0.93, 0.4);
    fLeg->SetFillColor(0);
    fLeg->SetLineColor(0);
    fLeg->SetShadowColor(0);
    fLeg->AddEntry(gr,label,"F");
  }
  else
    fLeg->AddEntry(gr,label,"PL");

  gr->SetMarkerColor(color);
  gr->SetLineColor  (color);
  gr->SetFillColor  (color);
  gr->SetLineStyle  (linesty);
  gr->SetLineWidth  (2);
  gr->SetMarkerStyle(marksty);
  gr->SetMarkerSize(1.2);

  el->SetLineColor  (color);
  el->SetFillColor  (color);
  el->SetFillStyle  (fillsty);
  el->SetLineStyle  (linesty);
  el->SetLineWidth  (2);

  CorrPlotItem item;
  item.graph = gr;
  item.ellipse = el;
  fItems.push_back(item);
}

void CorrPlot::Draw(TCanvas *c, TString fname, Int_t nlumi) {

  fItems[0].graph->SetTitle(fTitle);
  fItems[0].graph->GetXaxis()->SetTitle(fXtitle);
  fItems[0].graph->GetYaxis()->SetTitle(fYtitle);

  fItems[0].graph->GetXaxis()->SetLimits(fXmin,fXmax);
  fItems[0].graph->GetYaxis()->SetRangeUser(fYmin,fYmax);

  fItems[0].graph->Draw("ap");
  fItems[0].ellipse->Draw("same s");

  char lumitxt[150];
  sprintf(lumitxt,"#scale[0.75]{%i{#bf{pb^{-1} (13 TeV)}}}", nlumi);
  TPaveText *lumi = new TPaveText(0.65,0.93,0.97,0.99,"NDC");
  lumi->SetFillStyle(0); lumi->SetShadowColor(0); lumi->SetLineColor(0);
  lumi->SetTextFont(62);
  lumi->AddText(lumitxt);
  TPaveText *prelim = new TPaveText(0.18,0.93,0.44,0.99,"NDC");
  prelim->SetFillStyle(0); prelim->SetShadowColor(0); prelim->SetLineColor(0);
  //prelim->SetTextFont(62);
  prelim->AddText("#bf{CMS} #scale[0.75]{#it{Preliminary}}");
  TPaveText *fewz = new TPaveText(0.20,0.80,0.58,0.89,"NDC");
  fewz->SetFillStyle(0); fewz->SetShadowColor(0); fewz->SetLineColor(0);
  fewz->SetTextFont(62);
  fewz->AddText("Acc. #times FEWZ NNLO Prediction");
  lumi->Draw();
  prelim->Draw();
  fewz->Draw();

  for (UInt_t i=1; i<fItems.size(); i++) {
    fItems[i].ellipse->Draw("same s");
    fItems[i].graph->Draw("same p");
  }
  
  fLeg->Draw();
  
  c->SaveAs(fname);
  
}
