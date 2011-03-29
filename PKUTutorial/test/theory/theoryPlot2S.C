void theoryPlot2S(){

  double norm = 0.25;

  gROOT->LoadMacro("../xsection/setTDRStyle_modified.C");
  setTDRStyle();

   TString ytitle = TString::Format("d^{2}#sigma/dp_{T}dy #times #Beta(#mu#mu)  (nb/(GeV/c))");

  enum enuth {_data=0,_pythia,_cem,_cascade};
  int grcolor[4] = {1,4, 2,8};
  for(int i=0;i<4;i++) grcolor[i]=1; //remove color
 
  TFile cmsFile("xsection_2s_rap_0_2.root");

  TGraphAsymmErrors* cmsdata = ((TMultiGraph*)cmsFile.Get("mg1"))->GetListOfGraphs()->At(1);
   //TGraphAsymmErrors* cmsdata1 = ((TMultiGraph*)cmsFile.Get("mg1"))->GetListOfGraphs()->At(1);
 
   TCanvas *c2 = new TCanvas();//"c1","title",800,600);
   TH1F *frame = gPad->DrawFrame(0,0.0002,30,1);
   //TH2F * frame = new TH2F("frame","", 1000,0.,30., 0, 0.001,10);

  frame->SetStats(0);
  frame->GetXaxis()->SetTitle("p_{T}^{Y} (GeV/c)");
  frame->GetYaxis()->SetTitle("d#sigma/dp_{T}(pp->#Upsilon(2S)X) * B(#Upsilon(2S)->#mu#mu) (nb/(GeV/c))");
  frame->GetYaxis()->SetTitle(ytitle);
  frame->GetYaxis()->SetTitleOffset(1.1);
  frame->Draw();
  //////////////////////////////////////////////
  ///// PYTHIA
  //////////////////////////////////////////////

  TFile* fp = new TFile("../GenCrossSection/genXSupsilon2S.root"); // open the file
  TH1F* hist1 = (TH1F*)gDirectory->Get("genPt"); // get the hist
  TH1F* genPtLargeBins = (TH1F*)gDirectory->Get("genPtLargeBins"); // get the hist
  TGraphAsymmErrors* genPtLargeBinsGraph = (TGraphAsymmErrors*)gDirectory->Get("genPtLargeBinsGraph");
//  hist1->SetLineWidth(2);
//  hist1->SetLineColor(4);
  //  hist1->Draw();

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
  TH1F* histcascade_dir = new TH1F("hcas1","dsigma/dptjpsi*BR(jpsi->mumu) (pb/GeV) vs ptjpsi ",100,0.,100.);
  TH1F* histcascade_chi2p = new TH1F("hcas2","dsigma/dptjpsi*BR(jpsi->mumu) (pb/GeV) vs ptjpsi ",100,0.,100.);

  ifstream in,in1;
  in.open("Upsilon-txt/upsilon2s-pt-1eta2-dir2s.txt");
  in1.open("Upsilon-txt/upsilon2s-pt-1eta2-chi2p.txt");

  Float_t pt_cas_dir[200], pterr_cas_dir[200], xsc_cas_dir[200];
  Float_t pt_cas_chi2p[200], pterr_cas_chi2p[200], xsc_cas_chi2p[200];
  for(int i=1; i<101; i++){
    in >> xsc_cas_dir[i] >> pt_cas_dir[i] ;//cascade_dir;
    in1 >> xsc_cas_chi2p[i] >> pt_cas_chi2p[i];//cascade_chi2p

    // from Ramona Vogt
//    xscen[i]=cem*2*3.14*pt[i]*0.0248;
    cout << "CASCADE_dir: pt[i]= " << pt_cas_dir[i] << " xsc_cas_dir[i]= "<< xsc_cas_dir[i]<<  endl;
    cout << "CASCADE_chi2: pt[i]= " << pt_cas_chi2p[i] << " xsc_cas_chi2p[i]"<< xsc_cas_chi2p[i]<<  endl;
//    pterr[i]=0; // not known
  }
  for(int i=1; i<101; i++){
    histcascade_dir->SetBinContent(i, xsc_cas_dir[i]);
    histcascade_chi2p->SetBinContent(i, xsc_cas_chi2p[i]);
  }
   histcascade_dir->Add(histcascade_chi2p);
  //histcem->SetLineWidth(2);
  //histcem->SetLineColor(2);
  cout<<"CASCADE intgrated cross section = "<<histcascade_dir->Integral("width")<<endl;

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

  TH1F* histcem = new TH1F("htot1","dsigma/dptjpsi*BR(jpsi->mumu) (pb/GeV) vs ptjpsi ",98,1.,50.);
  ifstream in;
  in.open("upsilon_cem.dat");

  Float_t pt[200], pterr[200], xscen[200];
  for(int i=1; i<97; i++){
    Float_t cem;
    in >> pt[i] >> cem;
    // from Ramona Vogt
    xscen[i]=cem*2*3.14*pt[i]*0.0248 *0.33;
//    cout << "CEM: pt[i]= " << pt[i] << " xscen[i]= "<< xscen[i]<<  endl;
    pterr[i]=0; // not known
  }
  for(int i=1; i<97; i++){
    histcem->SetBinContent(i, xscen[i]);
  }
//  histcem->SetLineWidth(2);
//  histcem->SetLineColor(2);

 
  ////////////////////////////
  /// NOW ALL TOGETHER!
  ////////////////////////////
   c2->cd();
   c2->SetLeftMargin(0.15);
   c2->SetLogy();
   frame->Draw();
   cmsdata->SetLineWidth(3);
   //cmsdata->SetMarkerColor(1);
   //cmsdata->SetMarkerStyle(20);
   //cmsdata->SetMarkerSize(0.5);
 
   cmsdata->SetLineColor(grcolor[_data]);
   hist1  ->SetLineColor(grcolor[_pythia]);
   histcem->SetLineColor(grcolor[_cem]);
   histcascade_dir   ->SetLineColor(grcolor[_cascade]);
   
   histcascade_dir->SetLineWidth(3);
   histcascade_dir->SetLineStyle(3);
   
   histcem->SetLineWidth(2);
   histcem->SetLineStyle(6);
 
   hist1->SetLineWidth(2);
   hist1->SetLineStyle(1);
 
  ///// cascade 
  histcascade_dir->Draw("sameHIST");
  ///// cem
  histcem->Draw("sameL");
  //// /pyt
  hist1->Draw("same HIST 9c");
  cmsdata->Draw("P");
//  cmsdata1->Draw("Psame");

  //cout<<"wrong2"<<endl;
  TLatex latex;
  latex.DrawLatex(21,0.03,"|y^{#Upsilon}|<2");

  TLatex latex2;
  latex2.DrawLatex(2,0.003,"#Upsilon(2S)");

  TLatex latex3;
  latex3.DrawLatex(2,0.001 , "#sqrt{s}=7 TeV, L= 3 pb^{-1}");

  TLegend *leg = new TLegend(0.65,0.65,0.9,0.9);
//  leg->SetBorderSize(0);
//  leg->SetFillColor(10);

  leg->AddEntry(cmsdata, "CMS data", "PML");
  
  //leg->AddEntry(gr, "Direct #Upsilon, LO NRQCD, Artoisenet", "FL");
  leg->AddEntry(hist1,"PYTHIA", "L");
  leg->AddEntry(histcem, "CEM", "L");
  leg->AddEntry(histcascade_dir, "CASCADE", "L");
 

//normalize
  histcascade_dir     ->Scale(norm);
  hist1  ->Scale(norm);
  histcem->Scale(norm);
  normalize(cmsdata,norm);


  leg->Draw();
  //cout<<"wrong 4"<<endl; 
  c2->Print("theoryPlot2S.gif");
  c2->Print("theoryPlot2S.eps");
  c2->Print("theoryPlot2S.pdf");
  //cout<<"wrong44"<<endl;
  hist1->Scale(2.16/6.47);
  hist1->Draw("sameHIST");
  //cout<<"wrong5"<<endl;
  c2->Print("theoryPlot2Sscaled.gif");
  c2->Print("theoryPlot2Sscaled.eps");
  c2->Print("theoryPlot2Sscaled.pdf");

  for(int i=0; i<cmsdata->GetN(); i++){
    cmsdata->GetEXhigh()[i+1] = 0;//genPtLargeBins->GetEXhigh()[i];
    cmsdata->GetEXlow()[i+1] = 0;// genPtLargeBins->GetEXlow()[i];
  }
  for(int i=0; i<genPtLargeBinsGraph->GetN(); i++){
    genPtLargeBinsGraph->GetY()[i] = genPtLargeBinsGraph->GetY()[i]*2.16/6.47*norm;
  }
  frame->Draw();
  genPtLargeBinsGraph->SetMarkerStyle(1);
  genPtLargeBinsGraph->SetMarkerSize(0.0);
  genPtLargeBinsGraph->Draw("same PEz");
  cmsdata->Draw("same P");
  leg = new TLegend(0.4,0.75,0.9,0.9);
  leg->AddEntry(cmsdata, "CMS data", "PMLE");
  leg->AddEntry(genPtLargeBinsGraph,"PYTHIA (normalized)", "L");
  leg->Draw();
  latex.DrawLatex(21,0.08,"|y^{#Upsilon}| < 2");
  latex2.DrawLatex(2,0.0015,"#Upsilon(2S)");
  latex3.DrawLatex(2,0.0004 , "#sqrt{s} = 7 TeV, L = 3 pb^{ -1}");
  c2->Print("theoryPlot2Sbinned.gif");
  c2->Print("theoryPlot2Sbinned.eps");
  c2->Print("theoryPlot2Sbinned.pdf");
  hist1->Draw("same HIST c");
  c2->Print("theoryPlot2Scombined.gif");
  c2->Print("theoryPlot2Scombined.eps");
  c2->Print("theoryPlot2Scombined.pdf");

  for(int i=0; i<cmsdata->GetN(); i++){
    cmsdata->GetX()[i+1] = genPtLargeBinsGraph->GetX()[i];
    cmsdata->GetEXhigh()[i+1] = 0;//genPtLargeBins->GetEXhigh()[i];
    cmsdata->GetEXlow()[i+1] = 0;//genPtLargeBins->GetEXlow()[i];
  }
  frame->Draw();
  hist1->Draw("same HIST c");
  cmsdata->Draw("same P");
  leg = new TLegend(0.4,0.75,0.9,0.9);
  leg->AddEntry(cmsdata, "CMS data", "PML");
  leg->AddEntry(hist1,"PYTHIA (normalized)", "L");
  leg->Draw();
  latex.DrawLatex(21,0.08,"|y^{#Upsilon}| < 2");
  latex2.DrawLatex(2,0.0015,"#Upsilon(2S)");
  latex3.DrawLatex(2,0.0004 , "#sqrt{s} = 7 TeV, L = 3 pb^{ -1}");
  c2->Print("theoryPlot2Sxpos.gif");
  c2->Print("theoryPlot2Sxpos.eps");
  c2->Print("theoryPlot2Sxpos.pdf");

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
  //cout << "aa " << nbin << endl;
}
