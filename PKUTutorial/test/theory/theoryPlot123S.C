void theoryPlot123S(){

  double norm = 0.25;

  //gROOT->SetStyle("Plain");
  gROOT->LoadMacro("../xsection/setTDRStyle_modified.C");
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.14);
  TString ytitle = TString::Format("d^{2}#sigma/dp_{T}dy #times #Beta(#mu#mu)  (nb/(GeV/c))");

  enum enuth {_data=0,_pythia,_cem,_cascade};
  int grcolor[4] = {1,4, 2,8};
  for(int i=0;i<4;i++) grcolor[i]=1; //remove color

  TCanvas *c = new TCanvas("c","title",1500,600);
  c->Divide(3,1,0,0);
  
  c->cd(1);
  TH1F *frame = gPad->DrawFrame(0,0.0002,31,1.1);
  frame->SetStats(0);
  frame->GetXaxis()->SetTitle("p_{T}^{Y} (GeV/c)");
  frame->GetYaxis()->SetTitle(ytitle);
  frame->GetYaxis()->SetTitleOffset(1.1);
  gPad->SetLogy();
  frame->Draw();
  
  TFile* fp = new TFile("../GenCrossSection/genXSupsilon1S.root"); // open the file
  TH1F* genPt = (TH1F*)gDirectory->Get("genPt"); // get the hist
  TH1F* genPtLargeBins = (TH1F*)gDirectory->Get("genPtLargeBins"); // get the hist
  TGraphAsymmErrors* genPtLargeBinsGraph = (TGraphAsymmErrors*)gDirectory->Get("genPtLargeBinsGraph");
  genPt->SetLineColor(grcolor[_pythia]);
  genPt->SetLineWidth(1);
  genPt->SetLineStyle(1);
  genPt->Scale(7.37/15.17*norm);
  genPt->Draw("same HIST c");
  
  TFile cmsFile("xsection_1s_rap_0_2.root");
  TGraphAsymmErrors* cmsdata = ((TMultiGraph*)cmsFile.Get("mg1"))->GetListOfGraphs()->At(1);
  for(int i=0; i<cmsdata->GetN(); i++){
    cmsdata->GetX()[i+1] = genPtLargeBinsGraph->GetX()[i];
    cmsdata->GetEXhigh()[i+1] = 0;
    cmsdata->GetEXlow()[i+1] = 0;
  }
  cmsdata->SetLineWidth(1);
  cmsdata->SetLineColor(grcolor[_data]);
  cmsdata->SetMarkerSize(.5);
  normalize(cmsdata,norm);
  cmsdata->Draw("P");
  leg = new TLegend(0.4,0.75,0.95,0.9);
  leg->AddEntry(cmsdata, "CMS data", "PML");
  leg->AddEntry(genPt,"PYTHIA (normalized)", "L");
  leg->Draw();
  TLatex latex;
  latex.DrawLatex(21,0.05,"|y^{#Upsilon}| < 2");
  latex.DrawLatex(2,0.0015,"#Upsilon(1S)");
  latex.DrawLatex(2,0.0004 , "#sqrt{s} = 7 TeV, L = 3 pb^{ -1}");
  
  
  c->cd(2);
  gPad->SetFillStyle(4000);
  TH1F *frame2 = gPad->DrawFrame(0.001,0.0002,31,1.1);
  frame2->SetStats(0);
  frame2->GetXaxis()->SetTitle("p_{T}^{Y} (GeV/c)");
//  frame2->GetXaxis()->CenterTitle();
  frame2->GetYaxis()->SetTitle("");
  frame2->GetYaxis()->SetTitleOffset(0);
  frame2->GetYaxis()->SetLabelSize(0);
  frame2->GetYaxis()->SetLabelOffset(0);
  gPad->SetLogy();
  frame2->Draw();

  TFile* fp2 = new TFile("../GenCrossSection/genXSupsilon2S.root"); // open the file
  TH1F* genPt2 = (TH1F*)gDirectory->Get("genPt"); // get the hist
  TH1F* genPtLargeBins2 = (TH1F*)gDirectory->Get("genPtLargeBins"); // get the hist
  TGraphAsymmErrors* genPtLargeBinsGraph2 = (TGraphAsymmErrors*)gDirectory->Get("genPtLargeBinsGraph");
  genPt2->SetLineColor(grcolor[_pythia]);
  genPt2->SetLineWidth(1);
  genPt2->SetLineStyle(1);
  genPt2->Scale(2.16/6.47*norm);
  genPt2->Draw("same HIST c");
  
  TFile cmsFile2("xsection_2s_rap_0_2.root");
  TGraphAsymmErrors* cmsdata2 = ((TMultiGraph*)cmsFile2.Get("mg1"))->GetListOfGraphs()->At(1);
  for(int i=0; i<cmsdata2->GetN(); i++){
    cmsdata2->GetX()[i+1] = genPtLargeBinsGraph2->GetX()[i];
    cmsdata2->GetEXhigh()[i+1] = 0;
    cmsdata2->GetEXlow()[i+1] = 0;
  }
  cmsdata2->SetLineWidth(1);
  cmsdata2->SetLineColor(grcolor[_data]);
  cmsdata2->SetMarkerSize(.5);
  normalize(cmsdata2,norm);
  cmsdata2->Draw("P");
  leg2 = new TLegend(0.4,0.75,0.95,0.9);
  leg2->AddEntry(cmsdata2, "CMS data", "PML");
  leg2->AddEntry(genPt2,"PYTHIA (normalized)", "L");
  leg2->Draw();
  TLatex latex;
  latex.DrawLatex(21,0.05,"|y^{#Upsilon}| < 2");
  latex.DrawLatex(2,0.0015,"#Upsilon(2S)");
  latex.DrawLatex(2,0.0004 , "#sqrt{s} = 7 TeV, L = 3 pb^{ -1}");
  
  
  c->cd(3);
  TH1F *frame3 = gPad->DrawFrame(0.001,0.0002,31,1.1);
  frame3->SetStats(0);
  frame3->GetXaxis()->SetTitle("p_{T}^{Y} (GeV/c)");
  frame3->GetYaxis()->SetTitle("");
  frame3->GetYaxis()->SetTitleOffset(0);
  frame3->GetYaxis()->SetLabelSize(0);
  frame3->GetYaxis()->SetLabelOffset(0);
  gPad->SetLogy();
  frame3->Draw();
  
  TFile* fp3 = new TFile("../GenCrossSection/genXSupsilon3S.root"); // open the file
  TH1F* genPt3 = (TH1F*)gDirectory->Get("genPt"); // get the hist
  TH1F* genPtLargeBins3 = (TH1F*)gDirectory->Get("genPtLargeBins"); // get the hist
  TGraphAsymmErrors* genPtLargeBinsGraph3 = (TGraphAsymmErrors*)gDirectory->Get("genPtLargeBinsGraph");
  genPt3->SetLineColor(grcolor[_pythia]);
  genPt3->SetLineWidth(1);
  genPt3->SetLineStyle(1);
  genPt3->Scale(1.02/1.7*norm);
  genPt3->Draw("same HIST c");
  
  TFile cmsFile3("xsection_3s_rap_0_2.root");
  TGraphAsymmErrors* cmsdata3 = ((TMultiGraph*)cmsFile3.Get("mg1"))->GetListOfGraphs()->At(1);
  for(int i=0; i<cmsdata3->GetN(); i++){
    cmsdata3->GetX()[i+1] = genPtLargeBinsGraph3->GetX()[i];
    cmsdata3->GetEXhigh()[i+1] = 0;
    cmsdata3->GetEXlow()[i+1] = 0;
  }
  cmsdata3->SetLineWidth(1);
  cmsdata3->SetLineColor(grcolor[_data]);
  cmsdata3->SetMarkerSize(.5);
  normalize(cmsdata3,norm);
  cmsdata3->Draw("P");
  leg3 = new TLegend(0.4,0.75,0.95,0.9);
  leg3->AddEntry(cmsdata3, "CMS data", "PML");
  leg3->AddEntry(genPt3,"PYTHIA (normalized)", "L");
  leg3->Draw();
  latex.DrawLatex(21,0.05,"|y^{#Upsilon}| < 2");
  latex.DrawLatex(2,0.0015,"#Upsilon(3S)");
  latex.DrawLatex(2,0.0004 , "#sqrt{s} = 7 TeV, L = 3 pb^{ -1}");

  
  c->Print("theoryPlot123Sxpos.gif");
  c->Print("theoryPlot123Sxpos.eps");
  c->Print("theoryPlot123Sxpos.pdf");
}

void normalize(TGraphAsymmErrors *gr, double norm) {
  double nbin    = gr->GetN();
  double* xval   = gr->GetX();
  double* yval   = gr->GetY();
  double* yerr_h = gr->GetEYhigh();
  double* yerr_l = gr->GetEYlow ();
  for(int i=0; i<nbin; i++) {
    gr->SetPoint(i,xval[i],yval[i]*norm);
    gr->SetPointEYhigh(i,yerr_h[i]*norm);
    gr->SetPointEYlow (i,yerr_l[i]*norm);
  }
}
