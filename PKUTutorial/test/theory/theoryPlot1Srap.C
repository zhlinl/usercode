theoryPlot1Srap(){
  gROOT->LoadMacro("../xsection/setTDRStyle_modified.C");
  setTDRStyle();
  enum enuth {_data=0,_pythia,_cem,_cascade};
  int grcolor[4] = {1,4, 2,8};
  for(int i=0;i<4;i++) grcolor[i]=1; //remove color
  TFile cmsFile("xsection_1s_y0.root");
  TGraphAsymmErrors* cmsdata = ((TMultiGraph*)cmsFile.Get("mg1"))->GetListOfGraphs()->At(0);
  TCanvas *c2 = new TCanvas();//"c1","title",800,600);
  TH1F *frame = gPad->DrawFrame(0,0.1,2,50);
  //TH2F * frame = new TH2F("frame","", 1000,0.,30., 0, 0.001,10);
  frame->SetStats(0);

  frame->GetXaxis()->SetTitle("|y^{#Upsilon}|");
  TString ytitle = TString::Format("d#sigma/d|y|* BR(#Upsilon(1S)#rightarrow#mu#mu)  (nb/(GeV/c))");
  frame->GetYaxis()->SetTitle(ytitle);

  frame->GetYaxis()->SetTitleOffset(0.95);
  c2->cd();
  c2->SetLeftMargin(0.15);
  c2->SetLogy();
  frame->Draw();

  cmsdata->SetLineWidth(3);
  cmsdata->SetLineColor(grcolor[_data]);

  TFile* fp = new TFile("../GenCrossSection/genXSupsilon1Srap.root"); // open the file
  TH1F* hist1 = (TH1F*)gDirectory->Get("genY"); // get the hist
  hist1  ->SetLineColor(grcolor[_pythia]);
  hist1->SetLineWidth(2);
  hist1->SetLineStyle(1);

///////////////////////////////////////////////
  ////// Color octet Pierre
  ///////////////////////////////////////////////
  /*
  ifstream in;
  if(opt==1){
  in.open("DirectJpsi_CMSfullrap_clean.dat");
  } 
  if(opt==2){
    in.open("DirectJpsi_CMSrap_0_1.4_clean.dat");
  }
  if(opt==3){
    in.open("DirectJpsi_CMSrap_1.4_2.4_clean.dat");
  }
  Float_t pt[200], xscen[200], xsmin[200], xsmax[200], xserrlow[200], xserrhigh[200],pterr[200], xscen_corr[200], xserrlow_corr[200], xserrhigh_corr[200];
  for(int i=0; i<39; i++){
    in >> pt[i] >> xscen[i] >> xserrlow[i] >> xserrhigh[i]; 
    cout << "pt[i]= " << pt[i] << " xscen[i]= "<< xscen[i] <<"+" << xserrlow[i]<<  "-" << xserrhigh[i]<<  endl;
    xserrlow[i]=fabs(xserrhigh[i]);
    xserrhigh[i]=xserrlow[i];
    pterr[i]=0;
    Float_t factor=1000*0.0598;
    xscen_corr[i]=xscen[i]*factor;
    xserrlow_corr[i]=xserrlow[i]*factor;
    xserrhigh_corr[i]=xserrhigh[i]*factor;
  }
  
  TGraphAsymmErrors* gr;
  gr= new TGraphAsymmErrors(39,pt,xscen_corr,pterr,pterr, xserrlow_corr,xserrhigh_corr);
  gr->SetFillColor(7);
  gr->SetLineColor(7);
  TGraph * grline;
  grline= new TGraph(39,pt,xscen_corr);
  grline->SetLineWidth(5);
  grline->SetLineColor(7);
  */

  ///////////////////////////////////////////
  ///// CASCADE
  ///////////////////////////////////////////
/*
  TFile* f1 = new TFile("cascade-pp-7000-upsilon.root"); // open the file
  f1->cd("BPH10_003");
  TH1F* hh=(TH1F*)gDirectory->Get("h101");
  TFile* f2 = new TFile("cascade-pp-7000-chi_b.root"); // open the file
  f2->cd("BPH10_003");
  TH1F* gg=(TH1F*)gDirectory->Get("h101");
  gg->Add(hh);
  gg->SetLineWidth(2);
  gg->SetLineColor(8);
*/
  ////////////////////////
  /// CEM
  ////////////////////////
/*
  TH1F* histcem = new TH1F("htot1","dsigma/dptjpsi*BR(jpsi->mumu) (pb/GeV) vs ptjpsi ",98,1.,50.);
  ifstream in;
  in.open("upsilon_cem.dat");

  Float_t pt[200], pterr[200], xscen[200];
  for(int i=1; i<97; i++){
    Float_t cem;
    in >> pt[i] >> cem;
    // from Ramona Vogt
    xscen[i]=cem*2*3.14*pt[i]*0.0248;
//    cout << "CEM: pt[i]= " << pt[i] << " xscen[i]= "<< xscen[i]<<  endl;
    pterr[i]=0; // not known
  }
  for(int i=1; i<97; i++){
    histcem->SetBinContent(i, xscen[i]);
  }
  histcem->SetLineWidth(2);
  histcem->SetLineColor(2);

*/
  ////////////////////////////
  /// NOW ALL TOGETHER!
  ////////////////////////////
  ///// pierre
  //gr->Draw("E4");
  ///// cascade 
  //gg->Draw("sameHIST");
  ///// cem
  //histcem->Draw("sameHIST");
  //// /pyt
  hist1->Draw("sameHIST");
  cmsdata->Draw("P");
  //cmsdata1->Draw("Psame");

  TLatex latex3;
  latex3.DrawLatex(0.2,0.2 , "#sqrt{s}=7 TeV, L= 3 pb^{-1}");

  TLegend *leg = new TLegend(0.55,0.75,0.89,0.90);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);

  leg->AddEntry(cmsdata, "CMS preliminary", "PML");
  
  //leg->AddEntry(gr, "Direct #Upsilon, LO NRQCD, Artoisenet", "FL");
  leg->AddEntry(hist1,"PYTHIA", "L");
  //leg->AddEntry(histcem, "CEM", "L");
  //leg->AddEntry(gg, "CASCADE", "L");

  leg->Draw();
  c2->Print("theoryPlot1Srap.gif");
  c2->Print("theoryPlot1Srap.eps");
  c2->Print("theoryPlot1Srap.pdf");

  hist1->Scale(7.52/15.17);
  hist1->Draw("sameHIST");
  c2->Print("theoryPlot1SrapScaled.gif");
  c2->Print("theoryPlot1SrapScaled.eps");
  c2->Print("theoryPlot1SrapScaled.pdf");
}
