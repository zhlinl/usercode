bool absrap = false;

void normalize(TGraphAsymmErrors *gr, double norm);

void  peakOverlay(bool allpk=0, int ipeak=3) {

  double normalization = allpk?0.5:1.0;  // rapidity normalization -- note normalization values are set by hand at the moment 

  gROOT->LoadMacro("setTDRStyle_modified.C");
  setTDRStyle();

  if(allpk)
    assert(!ipeak);
  else
    assert(ipeak);
    
  const int ngr = allpk?3:2;
  
  TString path("mode");
  if(allpk) path+="0";
  else    path+="1";
  path+="/";
  

  TString title[ngr];
  if(allpk) {
    title[0] = "#Upsilon(1S)"; title[1] = "#Upsilon(2S)"; title[2] = "#Upsilon(3S)";
  } else {
    title[0] = "    |y|<1"; title[1]="1<|y|<2";
  }
  
  TString nfile[ngr];
  for(int i=0; i<ngr; i++)
    if(allpk)
      nfile[i] = TString::Format("%sxsection_%ds_y0.root",path.Data(),i+1);
    else 
      nfile[i] = TString::Format("%sxsection_%ds_y%d.root",path.Data(),ipeak,i);

  TFile* file[ngr];
  for(int i=0; i<ngr; i++)
    file[i] = TFile::Open(nfile[i],"read");
  
  TGraphAsymmErrors* gr[ngr], *grsta[ngr], *grsys[ngr];
  enum unc {all,sys,sta,cen};

  TFile* file[ngr];
  for(int i=0; i<ngr; i++) {
    file [i] = TFile::Open(nfile[i],"read");  
    gr   [i] = (TGraphAsymmErrors*)((TMultiGraph*)file[i]->Get("mg1"))->GetListOfGraphs()->At(all);//all
    grsys[i] = (TGraphAsymmErrors*)((TMultiGraph*)file[i]->Get("mg1"))->GetListOfGraphs()->At(sys);//all but lumi
    grsta[i] = (TGraphAsymmErrors*)((TMultiGraph*)file[i]->Get("mg1"))->GetListOfGraphs()->At(sta);//sta
  }
  
  for(int i=0; i<ngr; i++) {
    normalize(gr   [i],normalization);
    normalize(grsta[i],normalization);
    normalize(grsys[i],normalization);
  }

  TPaveText *tp = new TPaveText(0.5,0.8,0.8,0.9,"brNDC");
  tp->SetBorderSize(0);
  tp->AddText("CMS,  #sqrt{s} = 7 TeV");
  if(allpk)
    tp->AddText("L = 3 pb^{-1},  |y|<2");
  else {
    tp->AddText("L = 3 pb^{-1}");
  }

  TPaveText *tpp = new TPaveText(0.75,0.5,0.86,0.6,"brNDC");
  tpp->SetBorderSize(0);
  tpp->AddText(TString::Format("#Upsilon(%dS)",ipeak)); 

  
  TLegend *leg;
  if(allpk)
    leg = new TLegend(0.7,0.6,0.9,0.75);	
  else
    leg = new TLegend(0.7,0.65,0.9,0.75);
	
  int marker[]={21,24,22};
  for(int i=0; i<ngr; i++) {
    //these are needed, unless we re-run xsecresults with newstyle
    gr[i]->SetMarkerSize (1.4);
    if(i==2)
      gr[i]->SetMarkerSize (1.7);
    gr[i]->SetMarkerColor(13);
    //gr[i]->SetMarkerColor(2+2*i);
    //gr[i]->SetLineColor  (2+2*i);
    gr[i]->SetMarkerStyle(marker[i]);
    //gr[i]->SetMarkerStyle(27-3*i);
    leg->AddEntry(gr[i],title[i],"lpe");
  }

  if(!allpk) {
    //gr[0]->SetMarkerSize (1.4);
    //  gr[1]->SetMarkerSize (1.4);
  //gr[1]->SetFillColor(0.);
  }
  TCanvas a;
  TH1F *hh = gPad->DrawFrame(0.001,0,30,1);

  for(int i=0; i<ngr; i++) {
    grsta[i]->SetLineColor  (38);
    grsta[i]->SetMarkerColor(38);
    grsta[i]->SetFillColor  (38);
    grsta[i]->SetMarkerStyle(21);
    grsta[i]->SetMarkerSize(0.8);

    grsys[i]->SetLineColor  (kRed);
    grsys[i]->SetMarkerColor(kRed);
    grsys[i]->SetMarkerStyle(25);	
    grsys[i]->SetMarkerSize(0.8);	
    grsys[i]->SetFillStyle(0);    	
    grsys[i]->SetLineWidth(2);
  }

  for(int i=0; i<ngr; i++) {
    gr   [i]->Draw("e0");
    grsys[i]->Draw("[]");
    gr   [i]->Draw("pxy");
    //grsta[i]->Draw("e2"); //disabled, to simplify plot
  }

  double dummyx[1] = {0.}, dummyy[1] = {0.};
  dummy = new TGraph(1, dummyx, dummyy);

  TLegend *uleg = new TLegend(0.15,0.2,0.45,0.3);	
  //uleg->AddEntry(gr[0],"total uncertainty","p");
  //uleg->AddEntry(grsta[0],"stat. unc.","p");
  uleg->AddEntry(grsys[0],"total unc. except","p");
  //uleg->AddEntry(dummy,"global lumi unc. (11%)","e0");
  uleg->AddEntry(dummy,"luminosity (11%)","e0");
  uleg->Draw();

  leg->Draw();
  tp->Draw();
  if(!allpk)
    tpp->Draw();

  hh->GetYaxis()->SetTitleOffset(1.08); //needed
  hh->GetYaxis()->SetLabelSize(0.0475); //0.05
  //hh->GetXaxis()->SetLabelSize(0.045); //for consistency with y
  hh->GetYaxis()->SetTitleSize(0.055); //needed to fit...

  //TString ytitle = 
  TString ytitle;
  if(allpk)
    if(absrap)
      ytitle = TString::Format("d^{2}#sigma/dp_{T}d|y| #times #Beta(#mu#mu) (nb/(GeV/c))");
    else
      ytitle = TString::Format("d^{2}#sigma/dp_{T}dy #times #Beta(#mu#mu) (nb/(GeV/c))");
  else
    if(absrap)
      ytitle = TString::Format("d^{2}#sigma/dp_{T}d|y| #times #Beta(#mu#mu) (nb/(GeV/c))");
    else
      ytitle = TString::Format("d^{2}#sigma/dp_{T}dy #times #Beta(#mu#mu) (nb/(GeV/c))");

  hh->GetYaxis()->SetTitle(ytitle);
  hh->GetXaxis()->SetTitle("p_{T}^{#Upsilon} (GeV/c)");
  if(absrap)
    hh->SetMinimum(0.0005);
  else
    hh->SetMinimum(0.0002);
  hh->SetMaximum(1.0001);
  hh->GetXaxis()->SetRangeUser(0,30);
  
  gPad->SetLogy();
  
  TString pname;
  if(allpk)
    pname = TString::Format("xsec_overlay");
  else
    pname = TString::Format("xsec_%ds_2ybin",ipeak);

  a.SaveAs(path+pname + ".pdf");
  a.SaveAs(path+pname + ".gif");
  
  for(int i=0; i<ngr; i++)
    file[i]->Close();

  return;
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
