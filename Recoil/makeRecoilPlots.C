{

  TFile *f = new TFile("testWe.root", "read");
  
  TTree *RawWe = (TTree*) f->Get("RawWe");
  TH1D *hRaw = new TH1D("hRaw", "hRaw", 40, 0, 80);
  RawWe->Draw("out_met>>hRaw", "weight");

  TTree *CorrWe = (TTree*) f->Get("CorrWe");
  TH1D *hCorr = new TH1D("hCorr", "hCorr", 40, 0, 80);
  CorrWe->Draw("out_met>>hCorr", "weight");

  TTree *CorrUpWe = (TTree*) f->Get("CorrUpWe");
  TH1D *hCorrUp = new TH1D("hCorrUp", "hCorrUp", 40, 0, 80);
  CorrUpWe->Draw("out_met>>hCorrUp", "weight");
  TTree *CorrDownWe = (TTree*) f->Get("CorrDownWe");
  TH1D *hCorrDown = new TH1D("hCorrDown", "hCorrDown", 40, 0, 80);
  CorrDownWe->Draw("out_met>>hCorrDown", "weight");

  TTree *LepScaleUpWe = (TTree*) f->Get("LepScaleUpWe");
  TH1D *hLepScaleUp = new TH1D("hLepScaleUp", "hLepScaleUp", 40, 0, 80);
  LepScaleUpWe->Draw("out_met>>hLepScaleUp", "weight");
  TTree *LepScaleDownWe = (TTree*) f->Get("LepScaleDownWe");
  TH1D *hLepScaleDown = new TH1D("hLepScaleDown", "hLepScaleDown", 40, 0, 80);
  LepScaleDownWe->Draw("out_met>>hLepScaleDown", "weight");

  TTree *LepResUpWe = (TTree*) f->Get("LepResUpWe");
  TH1D *hLepResUp = new TH1D("hLepResUp", "hLepResUp", 40, 0, 80);
  LepResUpWe->Draw("out_met>>hLepResUp", "weight");
  TTree *LepResDownWe = (TTree*) f->Get("LepResDownWe");
  TH1D *hLepResDown = new TH1D("hLepResDown", "hLepResDown", 40, 0, 80);
  LepResDownWe->Draw("out_met>>hLepResDown", "weight");

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

  c1->SaveAs("rawVsCorrWe.png");

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

  c1->SaveAs("recoilVariationsWe.png");

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

  c1->SaveAs("lepScaleVariationsWe.png");

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

  c1->SaveAs("lepResVariationsWe.png");

}
