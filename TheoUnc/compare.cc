#include "compare.hh"

void compare() { return; }

void drawFourPt(Channel chan, TString var, Config fB, Config f1, Config f2, Config f3, Int_t nbinsPt, Int_t xminPt, Int_t xmaxPt) {

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);

  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.3,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.01); //0.01
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetRightMargin(0.07);
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);

  c1->cd(2)->SetPad(0,0,1.0,0.3);
  c1->cd(2)->SetTopMargin(0.01);
  c1->cd(2)->SetBottomMargin(0.45);//0.25
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.100,"Y");

  c1->cd(1);

  TH1D *hB = returnPlot(chan, fB, nbinsPt, xminPt, xmaxPt, var);
  TH1D *h1 = returnPlot(chan, f1, nbinsPt, xminPt, xmaxPt, var);
  TH1D *h2 = returnPlot(chan, f2, nbinsPt, xminPt, xmaxPt, var);
  TH1D *h3 = returnPlot(chan, f3, nbinsPt, xminPt, xmaxPt, var);

  hB->SetLineColor(kBlack);
  h1->SetLineColor(kBlue);h1->SetLineStyle(9);
  h2->SetLineColor(kRed);h2->SetLineStyle(7);
  h3->SetLineColor(kGreen+2);h3->SetLineStyle(3);

  hB->SetLineWidth(3);h1->SetLineWidth(3);h2->SetLineWidth(3);h3->SetLineWidth(3);

  hB->SetTitle("");
  hB->GetYaxis()->SetTitle("Unit Norm.");

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
  xtitle +=" p_{T}  [GeV]";

  hB->GetYaxis()->SetRangeUser(0.0001, 1.3*max(max( max(hB->GetMaximum(), h1->GetMaximum()), h2->GetMaximum()), h3->GetMaximum()));
  hB->GetYaxis()->SetNdivisions(310,kTRUE);

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");
  h2->DrawCopy("histsame");
  h3->DrawCopy("histsame");

  TLegend *leg = new TLegend(0.5, 0.6, 0.9, 0.85);
  leg->SetShadowColor(0); leg->SetFillColor(0); leg->SetLineColor(0);
  leg->AddEntry(hB, fB.getLabel(), "l");
  leg->AddEntry(h1, f1.getLabel(), "l");
  leg->AddEntry(h2, f2.getLabel(), "l");
  leg->AddEntry(h3, f3.getLabel(), "l");
  leg->Draw();

  c1->cd(2);

  TH1D *r3 = returnRelDiff(h3,hB,"h3B");
  TH1D *r2 = returnRelDiff(h2,hB,"h2B");
  TH1D *r1 = returnRelDiff(h1,hB,"h1B");
  TH1D *rB = returnRelDiff(hB,hB,"hBB");

  rB->GetYaxis()->SetRangeUser(-0.5,0.5);
  rB->GetXaxis()->SetTitle(xtitle);
  rB->GetYaxis()->SetTitle("Rel. Diff");

  rB->DrawCopy("hist");
  r1->DrawCopy("histsame");
  r2->DrawCopy("histsame");
  r3->DrawCopy("histsame");

  c1->SaveAs(chan_name[chan]+"/"+fB.getSampleName()+"_compare_"+var+".png");

  delete rB; delete r1; delete r2; delete r3;
  delete h3; delete h1; delete h2; delete hB; 
  delete c1; delete leg;
  rB=0; r1=0; r2=0; r3=0;
  h1=0; h2=0; hB=0; h3=0;
  c1=0; leg=0;

}

void drawThreePt(Channel chan, TString var, Config fB, Config f1, Config f2, Int_t nbinsPt, Int_t xminPt, Int_t xmaxPt) {

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);

  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.3,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.01); //0.01
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetRightMargin(0.07);
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);

  c1->cd(2)->SetPad(0,0,1.0,0.3);
  c1->cd(2)->SetTopMargin(0.01);
  c1->cd(2)->SetBottomMargin(0.45);//0.25
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.100,"Y");

  c1->cd(1);

  TH1D *hB = returnPlot(chan, fB, nbinsPt, xminPt, xmaxPt, var);
  TH1D *h1 = returnPlot(chan, f1, nbinsPt, xminPt, xmaxPt, var);
  TH1D *h2 = returnPlot(chan, f2, nbinsPt, xminPt, xmaxPt, var);

  hB->SetLineColor(kBlack);
  h1->SetLineColor(kBlue);h1->SetLineStyle(9);
  h2->SetLineColor(kRed);h2->SetLineStyle(7);

  hB->SetLineWidth(3);h1->SetLineWidth(3);h2->SetLineWidth(3);

  hB->SetTitle("");
  hB->GetYaxis()->SetTitle("Unit Norm.");

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
  xtitle +=" p_{T}  [GeV]";

  hB->GetYaxis()->SetRangeUser(0.0001, 1.3*max( max(hB->GetMaximum(), h1->GetMaximum()), h2->GetMaximum()));
  hB->GetYaxis()->SetNdivisions(310,kTRUE);

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");
  h2->DrawCopy("histsame");

  TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.85);
  leg->SetShadowColor(0); leg->SetFillColor(0); leg->SetLineColor(0);
  leg->AddEntry(hB, fB.getLabel(), "l");
  leg->AddEntry(h1, f1.getLabel(), "l");
  leg->AddEntry(h2, f2.getLabel(), "l");
  leg->Draw();

  c1->cd(2);

  TH1D *r2 = returnRelDiff(h2,hB,"h2B");
  TH1D *r1 = returnRelDiff(h1,hB,"h1B");
  TH1D *rB = returnRelDiff(hB,hB,"hBB");

  rB->GetYaxis()->SetRangeUser(-0.5,0.5);
  rB->GetXaxis()->SetTitle(xtitle);
  rB->GetYaxis()->SetTitle("Rel. Diff");

  rB->DrawCopy("hist");
  r1->DrawCopy("histsame");
  r2->DrawCopy("histsame");

  c1->SaveAs(chan_name[chan]+"/"+fB.getSampleName()+"_"+f1.getSampleName()+"_"+f2.getSampleName()+"_"+var+".png");

  delete rB; delete r1; delete r2;
  delete h1; delete h2; delete hB; 
  delete c1; delete leg;
  rB=0; r1=0; r2=0;
  h1=0; h2=0; hB=0; 
  c1=0; leg=0;

}
void drawTwoPt(Channel chan, TString var, Config fB, Config f1, Int_t nbinsPt, Int_t xminPt, Int_t xmaxPt) {

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);

  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.3,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.01); //0.01
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetRightMargin(0.07);
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);

  c1->cd(2)->SetPad(0,0,1.0,0.3);
  c1->cd(2)->SetTopMargin(0.01);
  c1->cd(2)->SetBottomMargin(0.45);//0.25
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.100,"Y");

  c1->cd(1);

  TH1D *hB = returnPlot(chan, fB, nbinsPt, xminPt, xmaxPt, var);
  TH1D *h1 = returnPlot(chan, f1, nbinsPt, xminPt, xmaxPt, var);

  hB->SetLineColor(kBlack);
  h1->SetLineColor(kBlue);h1->SetLineStyle(9);

  hB->SetLineWidth(3);h1->SetLineWidth(3);

  hB->SetTitle("");
  hB->GetYaxis()->SetTitle("Unit Norm.");

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

  hB->GetYaxis()->SetRangeUser(0.0001, 1.3*max(hB->GetMaximum(), h1->GetMaximum()));
  hB->GetYaxis()->SetNdivisions(310,kTRUE);

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");

  TLegend *leg = new TLegend(0.5, 0.65, 0.9, 0.85);
  leg->SetShadowColor(0); leg->SetFillColor(0); leg->SetLineColor(0);
  leg->AddEntry(hB, fB.getLabel(), "l");
  leg->AddEntry(h1, f1.getLabel(), "l");
  leg->Draw();

  c1->cd(2);

  TH1D *r1 = returnRelDiff(h1,hB,"h1B");
  TH1D *rB = returnRelDiff(hB,hB,"hBB");

  rB->GetYaxis()->SetRangeUser(-0.5,0.5);
  rB->GetXaxis()->SetTitle(xtitle);
  rB->GetYaxis()->SetTitle("Rel. Diff");

  rB->DrawCopy("hist");
  r1->DrawCopy("histsame");

  c1->SaveAs(chan_name[chan]+"/"+fB.getSampleName()+"_"+f1.getSampleName()+"_"+var+".png");

  delete h1; delete hB; delete c1; delete leg;
  delete r1; delete rB;
  h1=0; hB=0; c1=0; leg=0;
  r1=0; rB=0;

}

void drawThreeEta(Channel chan, TString var, Config fB, Config f1, Config f2, Int_t nbinsEta, Int_t xminEta, Int_t xmaxEta) {

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);

  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.3,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.01); //0.01
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetRightMargin(0.07);
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);

  c1->cd(2)->SetPad(0,0,1.0,0.3);
  c1->cd(2)->SetTopMargin(0.01);
  c1->cd(2)->SetBottomMargin(0.45);//0.25
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.100,"Y");

  c1->cd(1);

  TH1D *hB = returnPlot(chan, fB, nbinsEta, xminEta, xmaxEta, var);
  TH1D *h1 = returnPlot(chan, f1, nbinsEta, xminEta, xmaxEta, var);
  TH1D *h2 = returnPlot(chan, f2, nbinsEta, xminEta, xmaxEta, var);

  hB->SetLineColor(kBlack);
  h1->SetLineColor(kBlue);h1->SetLineStyle(9);
  h2->SetLineColor(kRed);h2->SetLineStyle(7);

  hB->SetLineWidth(3);h1->SetLineWidth(3);h2->SetLineWidth(3);

  hB->SetTitle("");
  hB->GetYaxis()->SetTitle("Unit Norm.");

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

  hB->GetYaxis()->SetRangeUser(0.0001, 1.3*max( max(hB->GetMaximum(), h1->GetMaximum()), h2->GetMaximum()));
  hB->GetYaxis()->SetNdivisions(310,kTRUE);

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");
  h2->DrawCopy("histsame");

  c1->cd(2);

  TH1D *r2 = returnRelDiff(h2,hB,"h2B");
  TH1D *r1 = returnRelDiff(h1,hB,"h1B");
  TH1D *rB = returnRelDiff(hB,hB,"hBB");

  rB->GetYaxis()->SetRangeUser(-0.5,0.5);
  rB->GetXaxis()->SetTitle(xtitle);
  rB->GetYaxis()->SetTitle("Rel. Diff");

  rB->DrawCopy("hist");
  r1->DrawCopy("histsame");
  r2->DrawCopy("histsame");

  TLegend *leg = new TLegend(0.35, 0.3, 0.75, 0.5);
  leg->SetShadowColor(0); leg->SetFillColor(0); leg->SetLineColor(0);
  leg->AddEntry(rB, fB.getLabel(), "l");
  leg->AddEntry(r1, f1.getLabel(), "l");
  leg->AddEntry(r2, f2.getLabel(), "l");
  leg->Draw();

  c1->SaveAs(chan_name[chan]+"/"+fB.getSampleName()+"_"+f1.getSampleName()+"_"+f2.getSampleName()+"_"+var+".png");

  delete h1; delete h2; delete hB; delete c1; delete leg;
  delete r1; delete r2; delete rB;
  h1=0; h2=0; hB=0; c1=0; leg=0;
  r1=0; r2=0; rB=0;

}

void drawTwoEta(Channel chan, TString var, Config fB, Config f1, Int_t nbinsEta, Int_t xminEta, Int_t xmaxEta) {

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);

  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.3,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.01); //0.01
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetRightMargin(0.07);
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);

  c1->cd(2)->SetPad(0,0,1.0,0.3);
  c1->cd(2)->SetTopMargin(0.01);
  c1->cd(2)->SetBottomMargin(0.45);//0.25
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.100,"Y");

  c1->cd(1);

  TH1D *hB = returnPlot(chan, fB, nbinsEta, xminEta, xmaxEta, var);
  TH1D *h1 = returnPlot(chan, f1, nbinsEta, xminEta, xmaxEta, var);

  hB->SetLineColor(kBlack);
  h1->SetLineColor(kBlue);h1->SetLineStyle(9);

  hB->SetLineWidth(3);h1->SetLineWidth(3);

  hB->SetTitle("");
  hB->GetYaxis()->SetTitle("Unit Norm.");

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

  hB->GetYaxis()->SetRangeUser(0.0001, 1.3*max(hB->GetMaximum(), h1->GetMaximum()));
  hB->GetYaxis()->SetNdivisions(310,kTRUE);

  hB->DrawCopy("hist");
  h1->DrawCopy("histsame");

  c1->cd(2);

  TH1D *r1 = returnRelDiff(h1,hB,"h1B");
  TH1D *rB = returnRelDiff(hB,hB,"hBB");

  rB->GetYaxis()->SetRangeUser(-0.5,0.5);
  rB->GetXaxis()->SetTitle(xtitle);
  rB->GetYaxis()->SetTitle("Rel. Diff");

  rB->DrawCopy("hist");
  r1->DrawCopy("histsame");

  TLegend *leg = new TLegend(0.35, 0.3, 0.75, 0.5);
  leg->SetShadowColor(0); leg->SetFillColor(0); leg->SetLineColor(0);
  leg->AddEntry(rB, fB.getLabel(), "l");
  leg->AddEntry(r1, f1.getLabel(), "l");
  leg->Draw();

  c1->SaveAs(chan_name[chan]+"/"+fB.getSampleName()+"_"+f1.getSampleName()+"_"+var+".png");

  delete h1; delete hB; delete c1; delete leg;
  delete r1; delete rB;
  h1=0; hB=0; c1=0; leg=0;
  r1=0; rB=0;

}

Double_t calcAcc(Channel chan, Config f) {
  
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
  //if (chan==wpm)      t1->Draw("genVf_eta>>sel", "(genV_id==24 && genL1_id==-13)*weight*(abs(genL1_eta)<2.1 && genL1_pt>25)");
  else if (chan==wpe) t1->Draw("genVf_eta>>sel", "(genV_id==24 && genL1_id==-11)*weight*(abs(genL1f_eta)<2.5 && (abs(genL1f_eta)<1.4442 || abs(genL1f_eta)>1.566) && genL1f_pt>25)");
  //else if (chan==wpe) t1->Draw("genVf_eta>>sel", "(genV_id==24 && genL1_id==-11)*weight*(abs(genL1_eta)<2.5 && (abs(genL1_eta)<1.4442 || abs(genL1_eta)>1.566) && genL1_pt>25)");
  else if (chan==wmm) t1->Draw("genVf_eta>>sel", "(genV_id==-24 && genL2_id==13)*weight*(abs(genL2f_eta)<2.1 && genL2f_pt>25)");
  //else if (chan==wmm) t1->Draw("genVf_eta>>sel", "(genV_id==-24 && genL2_id==13)*weight*(abs(genL2_eta)<2.1 && genL2_pt>25)");
  else if (chan==wme) t1->Draw("genVf_eta>>sel", "(genV_id==-24 && genL2_id==11)*weight*(abs(genL2f_eta)<2.5 && (abs(genL2f_eta)<1.4442 || abs(genL2f_eta)>1.566) && genL2f_pt>25)");
  //else if (chan==wme) t1->Draw("genVf_eta>>sel", "(genV_id==-24 && genL2_id==11)*weight*(abs(genL2_eta)<2.5 && (abs(genL2_eta)<1.4442 || abs(genL2_eta)>1.566) && genL2_pt>25)");
  else if (chan==zmm) t1->Draw("genVf_eta>>sel", "(genV_m>60 && genV_m<120)*(genL1_id==-13 && genL2_id==13)*weight*(abs(genL1f_eta)<2.1 && genL1f_pt>25)*(abs(genL2f_eta)<2.1 && genL2f_pt>25)");
  //else if (chan==zmm) t1->Draw("genVf_eta>>sel", "(genV_m>60 && genV_m<120)*(genL1_id==-13 && genL2_id==13)*weight*(abs(genL1_eta)<2.1 && genL1_pt>25)*(abs(genL2_eta)<2.1 && genL2_pt>25)");
  else if (chan==zee) t1->Draw("genVf_eta>>sel", "(genV_m>60 && genV_m<120)*(genL1_id==-11 && genL2_id==11)*weight*(abs(genL1f_eta)<2.5 && genL1f_pt>25)*(abs(genL2f_eta)<2.5 && genL2f_pt>25)*(abs(genL1f_eta)<1.4442 || abs(genL1f_eta)>1.566)*(abs(genL2f_eta)<1.4442 || abs(genL2f_eta)>1.566)");
  //else if (chan==zee) t1->Draw("genVf_eta>>sel", "(genV_m>60 && genV_m<120)*(genL1_id==-11 && genL2_id==11)*weight*(abs(genL1_eta)<2.5 && genL1_pt>25)*(abs(genL2_eta)<2.5 && genL2_pt>25)*(abs(genL1_eta)<1.4442 || abs(genL1_eta)>1.566)*(abs(genL2_eta)<1.4442 || abs(genL2_eta)>1.566)");

  Double_t a_t = all->Integral();
  Double_t a_s = sel->Integral();

  std::cout << chan_name[chan] << ", " << f.getSampleName() << ",  " << std::setprecision(4) << a_s/a_t << ", " << std::setprecision(2) << sqrt((a_s/a_t)*(1-a_s/a_t)/all->GetEntries())  << std::endl;
  delete sel; delete all; delete t1;

  return a_s/a_t;

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

TH1D* returnRelDiff(TH1D* h, TH1D* b, TString name) {
  TH1D* hRelDiff = new TH1D(name, "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  hRelDiff->SetLineColor(h->GetLineColor());
  hRelDiff->SetLineStyle(h->GetLineStyle());
  hRelDiff->SetLineWidth(h->GetLineWidth());

  hRelDiff->GetYaxis()->SetTitleOffset(0.42);
  hRelDiff->GetYaxis()->SetTitleSize(0.13);
  hRelDiff->GetYaxis()->SetLabelSize(0.10);
  hRelDiff->GetXaxis()->SetTitleOffset(1.2);
  hRelDiff->GetXaxis()->SetTitleSize(0.13);
  hRelDiff->GetXaxis()->SetLabelSize(0.12);
  //hRelDiff->GetXaxis()->CenterTitle();
  hRelDiff->GetYaxis()->CenterTitle();
  hRelDiff->GetYaxis()->SetNdivisions(303,kTRUE);

  for (Int_t i=1; i<h->GetNbinsX()+1; i++) {
    Double_t val = h->GetBinContent(i) - b->GetBinContent(i);
    if (b->GetBinContent(i)>0) hRelDiff->SetBinContent(i, val/b->GetBinContent(i));
    else hRelDiff->SetBinContent(i, 0);
  }
  return hRelDiff;
}
