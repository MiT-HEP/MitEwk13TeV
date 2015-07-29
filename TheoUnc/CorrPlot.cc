#include <CorrPlot.hh>
#include <iostream>

CorrPlot::CorrPlot()
{
  CorrPlot("wtf", "is", "going", "on");
}

CorrPlot::CorrPlot(TString name, TString title, TString xtitle, TString ytitle):
  fName(name),
  fTitle(title),
  fXtitle(xtitle),
  fYtitle(ytitle)
{}

void CorrPlot::AddCorrPlot(TGraph *gr, TEllipse *el, TString label, int color, int marksty, int linesty) {
  if (!gr || !el) return;
  if (!fLeg) {
    fLeg = new TLegend(0.6, 0.2, 0.93, 0.4);
    fLeg->SetFillColor(0);
    fLeg->SetLineColor(0);
    fLeg->SetShadowColor(0);
  }
  fLeg->AddEntry(gr,label,"P");


  gr->SetMarkerColor(color);
  gr->SetLineColor  (color);
  gr->SetFillColor  (color);
  gr->SetLineStyle  (linesty);
  gr->SetLineWidth  (2);
  gr->SetMarkerStyle(marksty);
  gr->SetMarkerSize(1.2);

  el->SetLineColor  (color);
  el->SetFillColor  (0);
  el->SetFillStyle  (4000);
  el->SetLineStyle  (linesty);
  el->SetLineWidth  (2);

  CorrPlotItem item;
  item.graph = gr;
  item.ellipse = el;
  fItems.push_back(item);
}

void CorrPlot::Draw(TCanvas *c, TString fname) {

  fItems[0].graph->SetTitle(fTitle);
  fItems[0].graph->GetXaxis()->SetTitle(fXtitle);
  fItems[0].graph->GetYaxis()->SetTitle(fYtitle);

  fItems[0].graph->GetXaxis()->SetLimits(8.2,8.6);
  fItems[0].graph->GetYaxis()->SetRangeUser(11.0,11.6);

  fItems[0].graph->Draw("ap");
  fItems[0].ellipse->Draw("same s");
  
  for (UInt_t i=1; i<fItems.size(); i++) {
    fItems[i].graph->Draw("same p");
    fItems[i].ellipse->Draw("same s");
  }
  
  fLeg->Draw();
  
  c->SaveAs(fname);
  
}
