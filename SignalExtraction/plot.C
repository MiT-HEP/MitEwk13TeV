void plot()
{
  TFile *ff = TFile::Open("Zmumu_pdfTemplates.root", "read");
  RooWorkspace* combine_workspace = (RooWorkspace*) ff->Get("combine_workspace");
  RooAbsData *data = combine_workspace->data("dataMetp");
  RooAbsPdf *metp = combine_workspace->pdf("wmp");
  RooAbsPdf *metpup = combine_workspace->pdf("wmp_RecoilUp");
  RooAbsPdf *metpdown = combine_workspace->pdf("wmp_RecoilDown");
  RooAbsPdf *metp_no = combine_workspace->pdf("wm");
  RooPlot *plot3 = combine_workspace->var("pfmet")->frame();
  data->plotOn(plot3,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),RooFit::Name("a"));
  metp->plotOn(plot3,LineColor(kBlue),RooFit::Name("b"));
  metpup->plotOn(plot3,LineColor(kGreen),RooFit::Name("c"));
//   metpdown->plotOn(plot3,LineColor(kRed),RooFit::Name("d"));
  data->plotOn(plot3,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),RooFit::Name("a"));
  //metp_no->plotOn(plot3,LineColor(kBlack));
  plot3->Draw();

  TLegend* leg = new TLegend(0.65, 0.65, 0.95, 0.90);
  leg->AddEntry(plot3->findObject("a")  , "data", "LP" );
  leg->AddEntry(plot3->findObject("b")  , "Z#rightarrow#mu#mu corrected", "L" );
  leg->AddEntry(plot3->findObject("c")  , "Z#rightarrow#mu#mu no correction", "L" );
//   leg->AddEntry(plot3->findObject("d")  , "Z#rightarrow#mu#mu down", "L" );
  //leg->AddEntry(plot3->findObject("res_sig")  , "HH->bb#gamma#gamma", "L" );
  plot3->SetTitle("");
  plot3->GetXaxis()->SetTitle("PF MET [GeV]");
  plot3->GetYaxis()->SetTitle("Events / 2.0 GeV");
  leg->Draw();
}
