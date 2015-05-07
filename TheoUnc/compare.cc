#include "compare.hh"

void compare() { return; }

void drawThreePt(Channel chan, TString var, Config fB, Config f1, Config f2, Int_t nbinsPt, Int_t xminPt, Int_t xmaxPt) {

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);

  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.5,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.05);
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetRightMargin(0.07);
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);

  c1->cd(2)->SetPad(0,0,1.0,0.5);
  c1->cd(2)->SetTopMargin(0.05);
  c1->cd(2)->SetBottomMargin(0.25);
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(1);

  c1->cd(1);

  TH1D *hB = returnPlot(chan, fB, nbinsPt, xminPt, xmaxPt, var);
  TH1D *h1 = returnPlot(chan, f1, nbinsPt, xminPt, xmaxPt, var);
  TH1D *h2 = returnPlot(chan, f2, nbinsPt, xminPt, xmaxPt, var);

  hB->SetLineColor(kBlack);
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);

  hB->SetTitle("");

  TString xtitle="";
  if (var=="genVf_pt") {
    if (chan<=wpe) {
      xtitle+="W";
      if (chan==wmm || chan==wme) xtitle+="-";
      else xtitle+="+";
    }
    else xtitle+="Z";
  }
  else if (var=="genL1f_pt") {
    if (chan==wpm || chan==zmm) xtitle+="#mu+";
    else if (chan==wpe || chan==zee) xtitle+="e+";
  }
  else if (var=="genL2f_pt") {
    if (chan==wmm || chan==zmm) xtitle+="#mu-";
    else if (chan==wme || chan==zee) xtitle+="e-";
  }
  xtitle +=" p_{T} (GeV)";

  hB->GetXaxis()->SetTitle(xtitle);

  hB->GetYaxis()->SetRangeUser(0, 1.3*max( max(hB->GetMaximum(), h1->GetMaximum()), h2->GetMaximum()));

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");
  h2->DrawCopy("histsame");

  TLegend *leg = new TLegend(0.5, 0.65, 0.9, 0.85);
  leg->SetShadowColor(0); leg->SetFillColor(0); leg->SetLineColor(0);
  leg->AddEntry(hB, fB.getLabel(), "l");
  leg->AddEntry(h1, f1.getLabel(), "l");
  leg->AddEntry(h2, f2.getLabel(), "l");
  leg->Draw();

  c1->cd(2);

  h2->Divide(hB);
  h1->Divide(hB);
  hB->Divide(hB);

  hB->GetYaxis()->SetRangeUser(0.5,1.5);

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");
  h2->DrawCopy("histsame");

  c1->SaveAs(chan_name[chan]+"/"+fB.getSampleName()+"_"+f1.getSampleName()+"_"+f2.getSampleName()+"_"+var+".png");

  delete h1; delete h2; delete hB; delete c1; delete leg;
  h1=0; h2=0; hB=0; c1=0; leg=0;

}
void drawTwoPt(Channel chan, TString var, Config fB, Config f1, Int_t nbinsPt, Int_t xminPt, Int_t xmaxPt) {

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);

  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.5,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.05);
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetRightMargin(0.07);
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);

  c1->cd(2)->SetPad(0,0,1.0,0.5);
  c1->cd(2)->SetTopMargin(0.05);
  c1->cd(2)->SetBottomMargin(0.25);
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(1);

  c1->cd(1);

  TH1D *hB = returnPlot(chan, fB, nbinsPt, xminPt, xmaxPt, var);
  TH1D *h1 = returnPlot(chan, f1, nbinsPt, xminPt, xmaxPt, var);

  hB->SetLineColor(kBlack);
  h1->SetLineColor(kRed);

  hB->SetTitle("");

  TString xtitle="";
  if (var=="genVf_pt") {
    if (chan<=wpe) {
      xtitle+="W";
      if (chan==wmm || chan==wme) xtitle+="-";
      else xtitle+="+";
    }
    else xtitle+="Z";
  }
  else if (var=="genL1f_pt") {
    if (chan==wpm || chan==zmm) xtitle+="#mu+";
    else if (chan==wpe || chan==zee) xtitle+="e+";
  }
  else if (var=="genL2f_pt") {
    if (chan==wmm || chan==zmm) xtitle+="#mu-";
    else if (chan==wme || chan==zee) xtitle+="e-";
  }
  xtitle +=" p_{T} (GeV)";

  hB->GetXaxis()->SetTitle(xtitle);

  hB->GetYaxis()->SetRangeUser(0, 1.3*max(hB->GetMaximum(), h1->GetMaximum()));

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");

  TLegend *leg = new TLegend(0.5, 0.65, 0.9, 0.85);
  leg->SetShadowColor(0); leg->SetFillColor(0); leg->SetLineColor(0);
  leg->AddEntry(hB, fB.getLabel(), "l");
  leg->AddEntry(h1, f1.getLabel(), "l");
  leg->Draw();

  c1->cd(2);

  h1->Divide(hB);
  hB->Divide(hB);

  hB->GetYaxis()->SetRangeUser(0.5,1.5);

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");

  c1->SaveAs(chan_name[chan]+"/"+fB.getSampleName()+"_"+f1.getSampleName()+"_"+var+".png");

  delete h1; delete hB; delete c1; delete leg;
  h1=0; hB=0; c1=0; leg=0;

}

void drawThreeEta(Channel chan, TString var, Config fB, Config f1, Config f2, Int_t nbinsEta, Int_t xminEta, Int_t xmaxEta) {

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);

  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.5,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.05);
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetRightMargin(0.07);
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);

  c1->cd(2)->SetPad(0,0,1.0,0.5);
  c1->cd(2)->SetTopMargin(0.05);
  c1->cd(2)->SetBottomMargin(0.25);
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(1);

  c1->cd(1);

  TH1D *hB = returnPlot(chan, fB, nbinsEta, xminEta, xmaxEta, var);
  TH1D *h1 = returnPlot(chan, f1, nbinsEta, xminEta, xmaxEta, var);
  TH1D *h2 = returnPlot(chan, f2, nbinsEta, xminEta, xmaxEta, var);

  hB->SetLineColor(kBlack);
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);

  hB->SetTitle("");

  TString xtitle="";
  if (var=="genVf_eta") {
    if (chan<=wpe) {
      xtitle+="W";
      if (chan==wmm || chan==wme) xtitle+="-";
      else xtitle+="+";
    }
    else xtitle+="Z";
  }
  else if (var=="genL1f_eta") {
    if (chan==wpm || chan==zmm) xtitle+="#mu+";
    else if (chan==wpe || chan==zee) xtitle+="e+";
  }
  else if (var=="genL2f_eta") {
    if (chan==wmm || chan==zmm) xtitle+="#mu-";
    else if (chan==wme || chan==zee) xtitle+="e-";
  }
  xtitle +=" #eta";

  hB->GetXaxis()->SetTitle(xtitle);

  hB->GetYaxis()->SetRangeUser(0, 1.3*max( max(hB->GetMaximum(), h1->GetMaximum()), h2->GetMaximum()));

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");
  h2->DrawCopy("histsame");

  c1->cd(2);

  h2->Divide(hB);
  h1->Divide(hB);
  hB->Divide(hB);

  hB->GetYaxis()->SetRangeUser(0.5,1.5);

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");
  h2->DrawCopy("histsame");

  TLegend *leg = new TLegend(0.35, 0.3, 0.75, 0.5);
  leg->SetShadowColor(0); leg->SetFillColor(0); leg->SetLineColor(0);
  leg->AddEntry(hB, fB.getLabel(), "l");
  leg->AddEntry(h1, f1.getLabel(), "l");
  leg->AddEntry(h2, f2.getLabel(), "l");
  leg->Draw();

  c1->SaveAs(chan_name[chan]+"/"+fB.getSampleName()+"_"+f1.getSampleName()+"_"+f2.getSampleName()+"_"+var+".png");

  delete h1; delete h2; delete hB; delete c1; delete leg;
  h1=0; h2=0; hB=0; c1=0; leg=0;

}

void drawTwoEta(Channel chan, TString var, Config fB, Config f1, Int_t nbinsEta, Int_t xminEta, Int_t xmaxEta) {

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);

  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.5,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.05);
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetRightMargin(0.07);
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);

  c1->cd(2)->SetPad(0,0,1.0,0.5);
  c1->cd(2)->SetTopMargin(0.05);
  c1->cd(2)->SetBottomMargin(0.25);
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(1);

  c1->cd(1);

  TH1D *hB = returnPlot(chan, fB, nbinsEta, xminEta, xmaxEta, var);
  TH1D *h1 = returnPlot(chan, f1, nbinsEta, xminEta, xmaxEta, var);

  hB->SetLineColor(kBlack);
  h1->SetLineColor(kRed);

  hB->SetTitle("");

  TString xtitle="";
  if (var=="genVf_eta") {
    if (chan<=wpe) {
      xtitle+="W";
      if (chan==wmm || chan==wme) xtitle+="-";
      else xtitle+="+";
    }
    else xtitle+="Z";
  }
  else if (var=="genL1f_eta") {
    if (chan==wpm || chan==zmm) xtitle+="#mu+";
    else if (chan==wpe || chan==zee) xtitle+="e+";
  }
  else if (var=="genL2f_eta") {
    if (chan==wmm || chan==zmm) xtitle+="#mu-";
    else if (chan==wme || chan==zee) xtitle+="e-";
  }
  xtitle +=" #eta";

  hB->GetXaxis()->SetTitle(xtitle);

  hB->GetYaxis()->SetRangeUser(0, 1.3*max(hB->GetMaximum(), h1->GetMaximum()));

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");

  c1->cd(2);

  h1->Divide(hB);
  hB->Divide(hB);

  hB->GetYaxis()->SetRangeUser(0.5,1.5);

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");

  TLegend *leg = new TLegend(0.35, 0.3, 0.75, 0.5);
  leg->SetShadowColor(0); leg->SetFillColor(0); leg->SetLineColor(0);
  leg->AddEntry(hB, fB.getLabel(), "l");
  leg->AddEntry(h1, f1.getLabel(), "l");
  leg->Draw();

  c1->SaveAs(chan_name[chan]+"/"+fB.getSampleName()+"_"+f1.getSampleName()+"_"+var+".png");

  delete h1; delete hB; delete c1; delete leg;
  h1=0; hB=0; c1=0; leg=0;

}

void calcAcc(Channel chan, Config f) {
  
  TTree *t1 = (TTree*) f.file->Get("Events");

  TH1D *all = new TH1D("all", "", 1, -100, 100);
  if (chan==wpm)      t1->Draw("genVf_eta>>all", "(genV_id==24 && genL1_id==-13)*weight");
  else if (chan==wpe) t1->Draw("genVf_eta>>all", "(genV_id==24 && genL1_id==-11)*weight");
  else if (chan==wmm) t1->Draw("genVf_eta>>all", "(genV_id==-24 && genL2_id==13)*weight");
  else if (chan==wme) t1->Draw("genVf_eta>>all", "(genV_id==-24 && genL2_id==11)*weight");
  else if (chan==zmm) t1->Draw("genVf_eta>>all", "(genV_m>60 && genV_m<120)*(genL1_id==-13 && genL2_id==13)*weight");
  else if (chan==zee) t1->Draw("genVf_eta>>all", "(genV_m>60 && genV_m<120)*(genL1_id==-11 && genL2_id==11)*weight");

  TH1D *sel = new TH1D("sel", "", 1, -100, 100);
  if (chan==wpm)      t1->Draw("genVf_eta>>sel", "(genV_id==24 && genL1_id==-13)*weight*(abs(genL1f_eta)<2.1 && genL1f_pt>25)");
  else if (chan==wpe) t1->Draw("genVf_eta>>sel", "(genV_id==24 && genL1_id==-11)*weight*(abs(genL1f_eta)<2.5 && (abs(genL1f_eta)<1.4442 || abs(genL1_eta)>1.566) && genL1_pt>25)");
  else if (chan==wmm) t1->Draw("genVf_eta>>sel", "(genV_id==-24 && genL2_id==13)*weight*(abs(genL2f_eta)<2.1 && genL2f_pt>25)");
  else if (chan==wme) t1->Draw("genVf_eta>>sel", "(genV_id==-24 && genL2_id==11)*weight*(abs(genL2f_eta)<2.5 && (abs(genL2f_eta)<1.4442 || abs(genL2f_eta)>1.566) && genL2f_pt>25)");
  else if (chan==zmm) t1->Draw("genVf_eta>>sel", "(genV_m>60 && genV_m<120)*(genL1_id==-13 && genL2_id==13)*weight*(abs(genL1f_eta)<2.1 && genL1f_pt>25)*(abs(genL2f_eta)<2.1 && genL2f_pt>25)");
  else if (chan==zee) t1->Draw("genVf_eta>>sel", "(genV_m>60 && genV_m<120)*(genL1_id==-11 && genL2_id==11)*weight*(abs(genL1f_eta)<2.5 && genL1f_pt>25)*(abs(genL2f_eta)<2.5 && genL2f_pt>25)*(abs(genL1f_eta)<1.4442 || abs(genL1f_eta)>1.566)*(abs(genL2f_eta)<1.4442 || abs(genL2f_eta)>1.566)");

  Double_t a_t = all->Integral();
  Double_t a_s = sel->Integral();

  //cout << setprecision(4) << a_s/a_t << " \pm " << setprecision(2) << sqrt((a_s/a_t)*(1-a_s/a_t)/all->GetEntries());
  std::cout << chan_name[chan] << std::endl;
  std::cout << "|  " << f.getSampleName() << "  |  " << std::setprecision(4) << a_s/a_t << " +/- " << std::setprecision(2) << sqrt((a_s/a_t)*(1-a_s/a_t)/all->GetEntries()) << "  |" << std::endl;

  delete sel; delete all; delete t1;

}


TH1D* returnPlot(Channel chan, Config f, Int_t nbins, Double_t xmin, Double_t xmax, TString var) {
  
  TTree *t1 = (TTree*) f.file->Get("Events");
  TH1D *h1 = new TH1D(f.getSampleName(), f.getSampleName(), nbins, xmin, xmax);

  if      (chan==wmm) t1->Draw(var+">>"+f.getSampleName(), "weight*(genV_id==-24 && genL2_id== 13)");
  else if (chan==wme) t1->Draw(var+">>"+f.getSampleName(), "weight*(genV_id==-24 && genL2_id== 11)");
  else if (chan==wpm) t1->Draw(var+">>"+f.getSampleName(), "weight*(genV_id== 24 && genL1_id==-13)");
  else if (chan==wpe) t1->Draw(var+">>"+f.getSampleName(), "weight*(genV_id== 24 && genL1_id==-11)");
  else if (chan==zmm) t1->Draw(var+">>"+f.getSampleName(), "weight*(genV_id== 23 && genL1_id==-13 && genL2_id==13)");
  else if (chan==zee) t1->Draw(var+">>"+f.getSampleName(), "weight*(genV_id== 23 && genL1_id==-11 && genL2_id==11)");

  h1->Scale(1.0/h1->Integral());
  
  return h1;

}
