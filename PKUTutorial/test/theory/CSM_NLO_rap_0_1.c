{
  TString ytitle = TString::Format("d^{2}#sigma/dp_{T}d|y| #times BR(#mu#mu)  (nb/(GeV/c))");

  TFile cmsFile("xsection_1s_rap_0_1.root");
  TGraphAsymmErrors* cmsdata = ((TMultiGraph*)cmsFile.Get("mg1"))->GetListOfGraphs()->At(0);
  float cmsY[6];
  float cmsX[6];
  for (int i=1; i<7; i++){
    cmsY[i-1]=cmsdata->GetY()[i]/2;
    cmsX[i-1]=cmsdata->GetX()[i];
    cout<<cmsY[i-1]<<" "<<cmsX[i-1]<<endl;
  }
  TGraph *cmsdata_direct = new TGraph(6,cmsX,cmsY); 

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
  float Nlo_high[8]={0.18,0.11, 0.013,0.0018,0.00036,0.000075,0.000025,0.000008};
  float Nlo_low[8] ={0.07,0.043,0.006,0.001, 0.0002, 0.000045,0.000015,0.000005};
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
  float NNlo_high[7]={0.27, 0.056, 0.01,  0.0026, 0.0007,  0.00021, 0.000076};
  float NNlo_low[7] ={0.045,0.0072,0.0015,0.00036,0.000093,0.000035,0.000012};
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
  //cmsdata_direct->Draw("same *");//directly produecd Y(1S)


  TLatex latex;
  latex.DrawLatex(2,0.00001,"|y^{#Upsilon}|<1");

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

  //c2->Print("U1S_CSM_y_0-1_noscale.gif");
  //c2->Print("U1S_CSM_y_0-1_noscale.eps");
  //c2->Print("U1S_CSM_y_0-1_noscale.pdf");
  c2->Print("U1S_CSM_y_0-1_scaled.gif");
  c2->Print("U1S_CSM_y_0-1_scaled.eps");
  c2->Print("U1S_CSM_y_0-1_scaled.pdf");

}
