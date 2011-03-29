{
  TString ytitle = TString::Format("d^{2}#sigma/dp_{T}d|y| #times BR(#mu#mu)  (nb/(GeV/c))");

  TFile cmsFile("xsection_1s_rap_1_2.root");
  TGraphAsymmErrors* cmsdata = ((TMultiGraph*)cmsFile.Get("mg1"))->GetListOfGraphs()->At(0);

  TCanvas *c2 = new TCanvas();//"c1","title",800,600);
  TH1F *frame = gPad->DrawFrame(0,0.000005,35,2);
  //TH2F * frame = new TH2F("frame","", 1000,0.,30., 0, 0.001,10);
  frame->SetStats(0);
  frame->GetXaxis()->SetTitle("p_{T}^{Y} (GeV/c)");
  frame->GetYaxis()->SetTitle("d#sigma/dp_{T} * BR (nb/(GeV/c))");
  frame->GetXaxis()->SetTitle("p_{T}^{Y} (GeV/c)");
  frame->GetYaxis()->SetTitle(ytitle);
  frame->GetYaxis()->SetTitleOffset(1.1);
  frame->Draw();

  float scale=1.0/0.5; //the feed down fraction factor

  float Nlo_pt[8]={3,5,10,15,20,25,30,35};
  float Nlo_pt_err[8]={0,0,0,0,0,0,0,0};
  float Nlo_high[8]={0.16,0.095,0.011,0.0015,0.0003,0.00006,0.00002,0.000006};
  float Nlo_low[8] ={0.065,0.04,0.005,0.00075,0.00013,0.00003,0.00001,0.0000035};
  float Nlo_mean[8];
  float Nlo_error[8];
  for (int i=0; i<8; i++){
    Nlo_mean[i]=(Nlo_high[i]+Nlo_low[i])/2.0;
    Nlo_error[i]=(Nlo_high[i]-Nlo_low[i])/2.0;
    Nlo_mean[i]=Nlo_mean[i]*scale;
    Nlo_error[i]=Nlo_error[i]*scale;
  }
  TGraphErrors* NLO = new TGraphErrors(8,Nlo_pt,Nlo_mean,Nlo_pt_err,Nlo_error);

  float NNlo_pt[7]={5,10,15,20,25,30,35};
  float NNlo_pt_err[7]={0,0,0,0,0,0,0};
  float NNlo_high[7]={0.24,0.048,0.0089,0.0022,0.00058,0.00018,0.000062};
  float NNlo_low[7] ={0.042,0.006,0.0012,0.00026,0.00007,0.000026,0.0000097};
  float NNlo_mean[7];
  float NNlo_error[7];
  for (int i=0; i<7; i++){
    NNlo_mean[i]=(NNlo_high[i]+NNlo_low[i])/2.0;
    NNlo_error[i]=(NNlo_high[i]-NNlo_low[i])/2.0;
    NNlo_mean[i]=NNlo_mean[i]*scale;
    NNlo_error[i]=NNlo_error[i]*scale;
  }
  TGraphErrors* NNLO = new TGraphErrors(7,NNlo_pt,NNlo_mean,NNlo_pt_err,NNlo_error);

  c2->cd();
  c2->SetLeftMargin(0.15);
  c2->SetLogy();
  frame->Draw();
  
  NLO->SetFillColor(7);
  NLO->SetLineColor(7);
  NLO->Draw("3");

  NNLO->SetFillColor(6);
  NNLO->SetFillStyle(3244);
  NNLO->SetLineColor(6);
  NNLO->Draw("same 3");

  cmsdata->SetLineWidth(3);
  //cmsdata->SetMarkerColor(1);
  //cmsdata->SetMarkerStyle(20);
  //cmsdata->SetMarkerSize(0.5);

  cmsdata->SetLineColor(1);
  cmsdata->Draw("same P");

  TLatex latex;
  latex.DrawLatex(2,0.00001,"1<|y^{#Upsilon}|<2");

  TLatex latex2;
  latex2.DrawLatex(2,0.0002,"#Upsilon(1S)");

  TLatex latex3;
  latex3.DrawLatex(2,0.00004 , "#sqrt{s}=7 TeV, L= 3 pb^{-1}");

  TLatex latex4;
  latex4.DrawLatex(8,0.002 , "Theory curves have been scaled up");
  TLatex latex5;
  latex5.DrawLatex(8,0.001 , "with the feed down fractions at CDF");

  TLegend *leg = new TLegend(0.65,0.65,0.89,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);

  leg->AddEntry(cmsdata, "CMS data", "PML");
  leg->AddEntry(NLO, "NLO prediction", "FL");
  leg->AddEntry(NNLO,"NNLO* prediction", "FL");
  leg->Draw();

  //c2->Print("U1S_CSM_y_1-2_noscale.gif");
  //c2->Print("U1S_CSM_y_1-2_noscale.eps");
  //c2->Print("U1S_CSM_y_1-2_noscale.pdf");
  c2->Print("U1S_CSM_y_1-2_scaled.gif");
  c2->Print("U1S_CSM_y_1-2_scaled.eps");
  c2->Print("U1S_CSM_y_1-2_scaled.pdf");

}
