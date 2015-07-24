{

  TFile *f = new TFile("test.root", "read");
  
  TTree *RawWm = (TTree*) f->Get("RawWm");
  TH1D *hRaw = new TH1D("hRaw", "hRaw", 40, 0, 80);
  RawWm->Draw("out_met>>hRaw", "weight");

  TTree *CorrWm = (TTree*) f->Get("CorrWm");
  TH1D *hCorr = new TH1D("hCorr", "hCorr", 40, 0, 80);
  CorrWm->Draw("out_met>>hCorr", "weight");

  TTree *CorrUpWm = (TTree*) f->Get("CorrUpWm");
  TH1D *hCorrUp = new TH1D("hCorrUp", "hCorrUp", 40, 0, 80);
  CorrUpWm->Draw("out_met>>hCorrUp", "weight");
  TTree *CorrDownWm = (TTree*) f->Get("CorrDownWm");
  TH1D *hCorrDown = new TH1D("hCorrDown", "hCorrDown", 40, 0, 80);
  CorrDownWm->Draw("out_met>>hCorrDown", "weight");

  TTree *LepScaleUpWm = (TTree*) f->Get("LepScaleUpWm");
  TH1D *hLepScaleUp = new TH1D("hLepScaleUp", "hLepScaleUp", 40, 0, 80);
  LepScaleUpWm->Draw("out_met>>hLepScaleUp", "weight");
  TTree *LepScaleDownWm = (TTree*) f->Get("LepScaleDownWm");
  TH1D *hLepScaleDown = new TH1D("hLepScaleDown", "hLepScaleDown", 40, 0, 80);
  LepScaleDownWm->Draw("out_met>>hLepScaleDown", "weight");

  TTree *LepResUpWm = (TTree*) f->Get("LepResUpWm");
  TH1D *hLepResUp = new TH1D("hLepResUp", "hLepResUp", 40, 0, 80);
  LepResUpWm->Draw("out_met>>hLepResUp", "weight");
  TTree *LepResDownWm = (TTree*) f->Get("LepResDownWm");
  TH1D *hLepResDown = new TH1D("hLepResDown", "hLepResDown", 40, 0, 80);
  LepResDownWm->Draw("out_met>>hLepResDown", "weight");

  TCanvas *c1 = MakeCanvas("c1", "", 800, 600);
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetShadowColor(0); leg->SetLineColor(0);
  hRaw->SetLineWidth(2); hCorr->SetLineWidth(2);
  hCorr->SetLineColor(kRed);
  leg->AddEntry(hRaw, "Uncorrected", "l");
  leg->AddEntry(hCorr, "Corrected", "l");
  
  hRaw->SetTitle("");
  hRaw->GetYaxis()->SetTitle("A.U.");
  hRaw->GetYaxis()->SetRangeUser(0, 1.2*TMath::Max(hRaw->GetMaximum(), hCorr->GetMaximum()));
  hRaw->GetXaxis()->SetTitle("MET");
  hRaw->Draw("hist");
  hCorr->Draw("histsame");
  leg->Draw();

  c1->SaveAs("rawVsCorr.png");

  hCorr->SetLineColor(kBlack); hCorrUp->SetLineColor(kRed); hCorrDown->SetLineColor(kBlue);
  hCorrUp->SetLineStyle(2); hCorrUp->SetLineWidth(2);
  hCorrDown->SetLineStyle(2); hCorrDown->SetLineWidth(2);
  leg->Clear();
  leg->AddEntry(hCorrUp, "+1#sigma", "l");
  leg->AddEntry(hCorr, "Nominal", "l");
  leg->AddEntry(hCorrDown, "-1#sigma", "l");

  hCorr->SetTitle("");
  hCorr->GetYaxis()->SetTitle("A.U.");
  hCorr->GetYaxis()->SetRangeUser(0, 1.2*hCorr->GetMaximum()); //1.2*TMath::Max(hCorr->GetMaximum(), hCorrUp->GetMaximum()));
  hCorr->GetXaxis()->SetTitle("MET");
  hCorr->Draw("hist");
  hCorrUp->Draw("histsame");
  hCorrDown->Draw("histsame");
  leg->Draw();

  c1->SaveAs("recoilVariations.png");

  hLepScaleUp->SetLineColor(kRed); hLepScaleDown->SetLineColor(kBlue);
  hLepScaleUp->SetLineStyle(2); hLepScaleUp->SetLineWidth(2);
  hLepScaleDown->SetLineStyle(2); hLepScaleDown->SetLineWidth(2);
  leg->Clear();
  leg->AddEntry(hLepScaleUp, "+1#sigma", "l");
  leg->AddEntry(hCorr, "Nominal", "l");
  leg->AddEntry(hLepScaleDown, "-1#sigma", "l");

  hCorr->SetTitle("");
  hCorr->GetYaxis()->SetTitle("A.U.");
  hCorr->GetYaxis()->SetRangeUser(0, hCorr->GetMaximum()); //1.2*TMath::Max(hCorr->GetMaximum(), hLepScaleUp->GetMaximum()));
  hCorr->GetXaxis()->SetTitle("MET");
  hCorr->Draw("hist");
  hLepScaleUp->Draw("histsame");
  hLepScaleDown->Draw("histsame");
  leg->Draw();

  c1->SaveAs("lepScaleVariations.png");

  hLepResUp->SetLineColor(kRed); hLepResDown->SetLineColor(kBlue);
  hLepResUp->SetLineStyle(2); hLepResUp->SetLineWidth(2);
  hLepResDown->SetLineStyle(2); hLepResDown->SetLineWidth(2);
  leg->Clear();
  leg->AddEntry(hLepResUp, "+1#sigma", "l");
  leg->AddEntry(hCorr, "Nominal", "l");
  leg->AddEntry(hLepResDown, "-1#sigma", "l");

  hCorr->SetTitle("");
  hCorr->GetYaxis()->SetTitle("A.U.");
  hCorr->GetYaxis()->SetRangeUser(0, hCorr->GetMaximum()); //1.2*TMath::Max(hCorr->GetMaximum(), hLepResUp->GetMaximum()));
  hCorr->GetXaxis()->SetTitle("MET");
  hCorr->Draw("hist");
  hLepResUp->Draw("histsame");
  hLepResDown->Draw("histsame");
  leg->Draw();

  c1->SaveAs("lepResVariations.png");

}
