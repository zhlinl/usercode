void plotGrPtDiff();
void plotGrYDiff();
void plotTwoGrs();
void plotTreeGrs();
void plotTwoGrsRatio();

const int ipeak = 1;

void plotgrs() {
  //plotTwoGrsRatio();
  //plotTwoGrs();
  //plotTreeGrs();
  //plotGrPtDiff();
  plotGrYDiff();
}


void plotGrYDiff() {

  const double ptbinw  =  30;//need to un-normalize
  //const double rapbinw = 0.4;//need to normalize (nb: bins are 0.2 wide and merged +/-, thus 0.4)
  const double rapbinw = 0.4*2;//need to normalize (nb: bins are 0.2 wide and merged +/-, thus 0.4)
  const double renorm = ptbinw/rapbinw;

  TFile file("xsection.root");
  gDirectory->Cd("plots");
  //const int nybins = 10;
  //double rapval[nybins] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9};//note: bin centers
  const int nybins = 5;
  double rapval[nybins] = {0.2, 0.6, 1.0, 1.4, 1.8};//note: bin centers
  TGraphAsymmErrors* grin[nybins];
  double xsecy [nybins] = {0.};
  double xsecye[nybins] = {0.};
  double raperr[nybins] = {0.};
  for (int i=0; i<nybins; i++) {
    TString gname = TString::Format("xsection_ups%dS_rap%d_pt",ipeak,i);
    grin  [i] = (TGraphAsymmErrors*) gROOT->FindObject(gname);
    xsecy [i] = (grin[i]->GetY())[0];
    xsecye[i] = (grin[i]->GetEYhigh())[0];
    raperr[i] = 0.2;
    xsecy [i] *= renorm;
    xsecye[i] *= renorm;
    printf("%s  %d->%7.5f +/- %7.5f\n",gname.Data(), i,xsecy[i],xsecye[i]);
  }
    
  TGraphAsymmErrors* thegr = new TGraphAsymmErrors(nybins, rapval, xsecy,
						   raperr,raperr,xsecye, xsecye);
						   
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0); 
  gStyle->SetFillColor(0);
  
  TString ytitle = TString::Format("#frac{d#sigma}{dy} . BR(#Upsilon(%dS)#rightarrow#mu#mu) (nb)",ipeak);
  thegr->SetTitle( "" );
  
  TAxis* ax = thegr->GetXaxis(); 
  ax->SetRangeUser(0,2);
  ax->SetTitle("#Upsilon rapidity");
  thegr->GetYaxis()->SetTitle(ytitle);
 
  thegr->SetMarkerColor(kBlue);
  thegr->SetLineColor(kBlue);
  thegr->SetLineWidth(2);

  TPaveText *tp = new TPaveText(0.5,0.3,0.9,0.5,"brNDC");
  tp->SetFillColor(0);
  tp->SetBorderSize(0);
  //tp->SetTextColor(kGray);
  tp->SetTextAlign(22);
  tp->SetTextSize(0.04);
  tp->AddText(TString::Format("#Upsilon(%dS) cross section",ipeak));


  thegr->SetMinimum(0);
  TCanvas a;
  thegr->Draw("ap");
  //tp->Draw("same");
  a.SaveAs(TString::Format("xsec%ds_ydiff.gif",ipeak));
}

void plotGrPtDiff() {

  bool single = true;
  bool posneg = false;

  TFile file("xsection.root");
  gDirectory->Cd("plots");
 
  TGraphAsymmErrors *gr0, *gr1;
  gr0 = (TGraphAsymmErrors*) gROOT->FindObject(TString::Format("xsection_ups%dS_rap0_pt",ipeak));
  if(!single)
    gr1 = (TGraphAsymmErrors*) gROOT->FindObject(TString::Format("xsection_ups%dS_rap1_pt",ipeak));
  else
    gr1 =  (TGraphAsymmErrors*)gr0->Clone();

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr0);
  mg->Add(gr1);
  
  TString ytitle = TString::Format("#frac{d#sigma}{dp_{T}} . BR(#Upsilon(%dS)#rightarrow#mu#mu) [nb/GeV/c]",ipeak);
  mg->SetTitle( "" );

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0); 
  gStyle->SetFillColor(0);

  TAxis* ax = mg->GetXaxis(); 
  //???ax->SetTitle("p_{T} (#mu#mu) [GeV/c]");
  //???mg->GetYaxis()->SetTitle(ytitle);
 
  gr0->SetMarkerColor(kPink);
  gr0->SetLineColor(kPink);
  gr0->SetLineWidth(2);
  //totGr->SetMarkerStyle(21);

  gr1->SetMarkerColor(kBlue);
  gr1->SetLineColor(kBlue);
  gr1->SetLineWidth(2);

   TLegend *leg = new TLegend(0.7,0.6,0.8,0.7);	
   if(single) {
     leg->AddEntry(gr1,"|y|<2","l");
   } else {
     if(posneg) {
       leg->AddEntry(gr0,"y<0","l");
       leg->AddEntry(gr1,"y>0","l");
     } else {
       leg->AddEntry(gr0,"0<|y|<1","l");
       leg->AddEntry(gr1,"1<|y|<2","l");
     }
   }

   double lumi = 3;
   TPaveText *tp = new TPaveText(0.5,0.75,0.9,0.85,"brNDC");
   tp->SetFillColor(0);
   tp->SetBorderSize(0);
   //tp->SetTextColor(kGray);
   tp->SetTextAlign(22);
   tp->SetTextSize(0.04);
   tp->AddText(TString::Format("#Upsilon(%dS) cross section",ipeak));
   tp->AddText(TString::Format("#sqrt{s}=7TeV preliminary, L= %2.1fpb^{-1}",lumi));


  TCanvas a;
  mg->Draw("ap");
  leg->Draw();
  tp->Draw();
  gPad->SetLogy();
  //a.SaveAs(TString::Format("gr_%ds.gif",ipeak));
  TString pname = TString::Format("xsec%ds",ipeak);
  if(posneg)
    pname += "_negative_positive_y";
  else
    pname += "_barrel_forward";
  pname += ".gif";
  a.SaveAs(pname);
  return;
}



void  plotTreeGrs() {

  const int ipeak =2;
  const int ngr = 2;
  TString title[] = {"0<|y^{#Upsilon}|<1","1<|y^{#Upsilon}|<2"};
  //TString title[] = {"#Upsilon(1S)", "#Upsilon(2S)", "#Upsilon(3S)"};
  

  TString nfile[ngr];

  nfile[0] = TString::Format("figs/xsection_%ds_y0.root",ipeak);
  nfile[1] = TString::Format("figs/xsection_%ds_y1.root",ipeak);

  //nfile[0] = "figs/xsection_1s_y0.root";
  //nfile[1] = "figs/xsection_2s_y0.root";
  //nfile[2] = "figs/xsection_3s_y0.root";

  TGraphAsymmErrors* gr[ngr];
  TMultiGraph *mg = new TMultiGraph();

  TFile* file[ngr];
  for(int i=0; i<ngr; i++) {
    file[i] = TFile::Open(nfile[i],"read");  
    gr[i] = (TGraphAsymmErrors*)((TMultiGraph*)file[i]->Get("mg1"))->GetListOfGraphs()->At(0);
    mg->Add(gr[i]);
  }


  //setTDRStyle();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0); 
  gStyle->SetFillColor(0);

  //TString ytitle = TString::Format("d#sigma/dp_{T} . BR(#Upsilon(%dS)#rightarrow#mu#mu)",ipeak);
  TString ytitle = TString::Format("d#sigma/dp_{T} . BR(#Upsilon#rightarrow#mu#mu)");

  TLegend *leg = new TLegend(0.7,0.5,0.85,0.7);	

  TPaveText *tp = new TPaveText(0.5,0.73,0.9,0.86,"brNDC");
  tp->SetFillColor(0);
  tp->SetBorderSize(0);
  tp->SetTextAlign(22);
  tp->SetTextColor(12);
  tp->SetTextSize(0.06);
  tp->AddText("CMS,  #sqrt{s} = 7 TeV");
  //tp->AddText("L=3pb^{-1},  |y(#Upsilon)|<2");
  tp->AddText("L=3pb^{-1}");

  for(int i=0; i<ngr; i++) {
    gr[i]->SetMarkerSize(0.7);
    gr[i]->SetMarkerColor(2+2*i);
    gr[i]->SetLineColor  (2+2*i);
    gr[i]->SetMarkerStyle(24+i);
    gr[i]->SetLineWidth(2);
    leg->AddEntry(gr[i],title[i],"lp");
  }

  TCanvas a;
  mg->Draw("ap");
  leg->Draw();
  tp->Draw();

  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetTitleSize(0.05);

  mg->GetYaxis()->SetTitle(ytitle);
  mg->GetXaxis()->SetTitle("p_{T} (#mu#mu) [GeV/c]");
  mg->SetMinimum(0.0005);
  mg->SetMaximum(5);
  mg->GetXaxis()->SetRangeUser(0,30);

  gPad->SetLogy();
  a.SaveAs(TString::Format("figs/xsec_%ds_2ybin.pdf",ipeak));
  //a.SaveAs("figs/xsec_overlay.gif");
  
}


void plotTwoGrs() {

  //TString nfile1("fitres_1s_upsilonYieldWeighted_cowboys/xsection.root");
  //TString nfile2("fitres_1s_upsilonYieldWeighted_seaguls/xsection.root");
  //TString nfile1("figs/xsection_1s_y0.root");
  //TString nfile2("figs/xsection_2s_y0.root");

  TFile* file1 = TFile::Open(nfile1,"read");
  TFile* file2 = TFile::Open(nfile2,"read");


  TGraphAsymmErrors *gr0, *gr1;
  gDirectory->Cd(nfile1+":/plots");
  gr0 = (TGraphAsymmErrors*) gROOT->FindObject(TString::Format("xsection_ups%dS_rap0_pt",ipeak));
  gDirectory->Cd(nfile2+":/plots");
  gr1 = (TGraphAsymmErrors*) gROOT->FindObject(TString::Format("xsection_ups%dS_rap0_pt",ipeak));

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr0);
  mg->Add(gr1);
  
  TString ytitle = TString::Format("#frac{d#sigma}{dp_{T}} . BR(#Upsilon(%dS)#rightarrow#mu#mu) [nb/GeV/c]",ipeak);
  mg->SetTitle( "" );

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0); 
  gStyle->SetFillColor(0);
  TAxis* ax = mg->GetXaxis(); 
  //???ax->SetTitle("p_{T} (#mu#mu) [GeV/c]");
  //???mg->GetYaxis()->SetTitle(ytitle);
 
  gr0->SetMarkerColor(kPink);
  gr0->SetLineColor(kPink);
  gr0->SetLineWidth(2);
  //totGr->SetMarkerStyle(21);

  gr1->SetMarkerColor(kBlue);
  gr1->SetLineColor(kBlue);
  gr1->SetLineWidth(2);

   TLegend *leg = new TLegend(0.7,0.6,0.8,0.7);	
   leg->AddEntry(gr0,title[0],"l");
   leg->AddEntry(gr1,title[1],"l");

   double lumi = 3;
   TPaveText *tp = new TPaveText(0.5,0.75,0.9,0.85,"brNDC");
   tp->SetFillColor(0);
   tp->SetBorderSize(0);
   //tp->SetTextColor(kGray);
   tp->SetTextAlign(22);
   tp->SetTextSize(0.04);
   tp->AddText(TString::Format("#Upsilon(%dS) cross section",ipeak));
   tp->AddText(TString::Format("#sqrt{s}=7TeV preliminary, L= %2.1fpb^{-1}",lumi));


  TCanvas a;
  mg->Draw("ap");
  leg->Draw();
  tp->Draw();
  gPad->SetLogy();
  //a.SaveAs(TString::Format("gr_%ds.gif",ipeak));
  TString pname = TString::Format("xsec%ds",ipeak);
  pname += "_" + title[0] + "_" + title[1];
  pname += ".gif";
  a.SaveAs(pname);
  return;
}


void plotTwoGrsRatio() {
  
  int peak1 = 3; 
  int peak2 = 0;
  const int nn = 3;

  TString nfile[nn];
  for(int i=0; i<nn; i++)
    nfile[i] = TString::Format("figs/xsection_%ds_y0.root",i+1);
  //TString nfile1("fitres_1s_upsilonYieldWeighted_nominal/xsection.root");
  //TString nfile2("fitres_2s_upsilonYieldWeighted_2s/xsection.root");
  //TString nfile3("fitres_3s_upsilonYieldWeighted_3s/xsection.root");

  TFile* file[nn];
  for(int i=0; i<nn; i++)
    file[i] = TFile::Open(nfile[i],"read");


  TGraphAsymmErrors *gr[nn], *gg[nn];

  for(int i=0; i<nn; i++) {
    gr[i] = (TGraphAsymmErrors*)((TMultiGraph*)file[i]->Get("mg1"))->GetListOfGraphs()->At(0);//-1-30
    gg[i] = (TGraphAsymmErrors*)((TMultiGraph*)file[i]->Get("mg1"))->GetListOfGraphs()->At(1);//0-30
  }
  //gDirectory->Cd(nfile1+":/plots");
  //gr1 = (TGraphAsymmErrors*) gROOT->FindObject(TString::Format("xsection_ups%dS_rap0_pt",1));
  //gDirectory->Cd(nfile2+":/plots");
  //gr2 = (TGraphAsymmErrors*) gROOT->FindObject(TString::Format("xsection_ups%dS_rap0_pt",2));
  //gDirectory->Cd(nfile3+":/plots");
  //gr3 = (TGraphAsymmErrors*) gROOT->FindObject(TString::Format("xsection_ups%dS_rap0_pt",3));


  //TCanvas c; gPad->SetLogy(); gg[0]->Draw("ap"); c.SaveAs("aa.gif");return;

  TGraphAsymmErrors *gr31, *gr32, *gr21;

  const int nn = gr[0]->GetN();;
  const int nbin = nn;

  double rat31[nbin], rat31_el[nbin], rat31_eh[nbin];
  double rat32[nbin], rat32_el[nbin], rat32_eh[nbin];
  double rat21[nbin], rat21_el[nbin], rat21_eh[nbin];
  double val1, err1l,err1h, val2,err2l,err2h, val3,err3l,err3h;
  for(int i=0; i<nbin; i++) {
    val1 = (gr[0]->GetY())[i];
    val2 = (gr[1]->GetY())[i];
    val3 = (gr[2]->GetY())[i];
    err1l = (gr[0]->GetEYlow())[i];
    err2l = (gr[1]->GetEYlow())[i];
    err3l = (gr[2]->GetEYlow())[i];
    err1h = (gr[0]->GetEYhigh())[i];
    err2h = (gr[1]->GetEYhigh())[i];
    err3h = (gr[2]->GetEYhigh())[i];
    rat31  [i] = val3/val1;
    rat32  [i] = val3/val2;
    rat21  [i] = val2/val1;
    rat31_el[i] = rat31[i]*sqrt(pow(err3l/val3,2)+pow(err1l/val1,2));
    rat32_el[i] = rat32[i]*sqrt(pow(err3l/val3,2)+pow(err2l/val2,2));
    rat21_el[i] = rat21[i]*sqrt(pow(err2l/val2,2)+pow(err1l/val1,2));
    rat31_eh[i] = rat31[i]*sqrt(pow(err3h/val3,2)+pow(err1h/val1,2));
    rat32_eh[i] = rat32[i]*sqrt(pow(err3h/val3,2)+pow(err2h/val2,2));
    rat21_eh[i] = rat21[i]*sqrt(pow(err2h/val2,2)+pow(err1h/val1,2));
  }

  char ss[50];
  ofstream flatex;
  flatex.open("xsec_ratio.tex");
  flatex << "\\begin{tabular}{cccc} \\hline \n";
  flatex << "\\pt (\\GeVc) & \\upsiii/\\upsi & \\upsiii/\\upsii & \\upsii/\\upsi \\\\ \\hline\n";
  bool xbars = 0;
  double* ptl = gr[0]->GetEXlow();
  double* pth = gr[0]->GetEXhigh();
  for(int i=0; i<nbin; i++) {
    if(!xbars) continue;
    cen =  gr[0]->GetX()[i];
    if(i) {
      ptl[i] = (gg[0]->GetEXlow()) [i-1];
      pth[i] = (gg[0]->GetEXhigh())[i-1];
    } else {
      cen = 0 ;
      ptl[i] = 0;      
      pth[i] = gr[0]->GetX()[nbin-1] + (gg[0]->GetEXhigh())[nbin-2];
    }
    //printf("cen:%5.2f ptl:%5.2f pth:%5.2f [%5.2f,%5.2f] \n", cen,ptl[i],pth[i],cen-ptl[i],cen+pth[i]);
    sprintf(ss,"%3.0f - %3.0f & $%4.2f^{+%3.2f}_{-%3.2f}$ & $%4.2f^{+%3.2f}_{-%3.2f}$ & $%4.2f^{+%3.2f}_{-%3.2f}$ \\\\ \\hline\n", 
	    cen-ptl[i],cen+pth[i],  
	    rat31[i], rat31_eh[i], rat31_el[i],
	    rat32[i], rat32_eh[i], rat32_el[i], 
	    rat21[i], rat21_eh[i], rat21_el[i]); 
    flatex << ss;
  }
  sprintf(ss,"\\end{tabular}\n");
  flatex << ss; 
  flatex.close();

  gr31 = new TGraphAsymmErrors(nbin, gr[0]->GetX(), rat31, ptl,pth, rat31_el, rat31_eh);
  gr32 = new TGraphAsymmErrors(nbin, gr[0]->GetX(), rat32, ptl,pth, rat32_el, rat32_eh);
  gr21 = new TGraphAsymmErrors(nbin, gr[0]->GetX(), rat21, ptl,pth, rat21_el, rat21_eh);

  TMultiGraph *mg = new TMultiGraph();
  if      (peak1 == 3 && peak2 == 1) mg->Add(gr31);
  else if (peak1 == 3 && peak2 == 2) mg->Add(gr32);
  else if (peak1 == 2 && peak2 == 1) mg->Add(gr21);
  else {mg->Add(gr31); mg->Add(gr32); mg->Add(gr21);}

  mg->SetMinimum(0);
  mg->SetMaximum(2);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0); 
  gStyle->SetFillColor(0);

  TString ytitle = TString::Format("#Upsilon(nS) cross section ratios");
  //TString ytitle = TString::Format("(d#sigma/dp_{T}.BR(#Upsilon(%dS)) / (d#sigma/dp_{T}.BR(#Upsilon(%dS))",peak1,peak2);
  mg->SetTitle( "" );

  gr31->SetMarkerColor(6);
  gr31->SetLineColor(6);
  gr31->SetLineWidth(2);
  gr32->SetMarkerColor(2);
  gr32->SetLineColor(2);
  gr32->SetLineWidth(2);
  gr21->SetMarkerColor(4);
  gr21->SetLineColor(4);
  gr21->SetLineWidth(2);
  gr31->SetMarkerStyle(24);//20);
  gr32->SetMarkerStyle(25);//25);
  gr21->SetMarkerStyle(26);//22);
  gr31->SetMarkerSize(0.7);
  gr32->SetMarkerSize(0.7);
  gr21->SetMarkerSize(0.7);

   TLegend *leg = new TLegend(0.2,0.45,0.5,0.65);	
   //TLegend *leg = new TLegend(0.2,0.55,0.5,0.75);	//tdr
   leg->AddEntry(gr31,"#Upsilon(3S)/#Upsilon(1S)","lp");
   leg->AddEntry(gr32,"#Upsilon(3S)/#Upsilon(2S)","lp");
   leg->AddEntry(gr21,"#Upsilon(2S)/#Upsilon(1S)","lp");

   TPaveText *tp = new TPaveText(0.2,0.7,0.45,0.85,"brNDC");
   //TPaveText *tp = new TPaveText(0.2,0.8,0.5,0.9,"brNDC");//tdr
   tp->SetBorderSize(0);
   tp->SetFillColor(0);
   tp->SetTextAlign(22);
   tp->SetTextColor(12);
   tp->SetTextSize(0.06);
   tp->AddText("CMS,  #sqrt{s} = 7 TeV");
   tp->AddText("L=3pb^{-1},  |y(#Upsilon)|<2");
   //tp->AddText("L=3pb^{-1}");

   TCanvas a;
   //mg->SetMaximum(1);
   mg->Draw("ap");

   TAxis* ax = mg->GetXaxis(); 
   ax->SetTitle("p_{T} (#mu#mu) (GeV/c)");
   mg->GetYaxis()->SetTitle(ytitle);
   mg->GetYaxis()->SetTitleSize(0.05);
   mg->GetXaxis()->SetTitleSize(0.05);
   mg->GetXaxis()->SetRangeUser(0,30);

   tp->Draw();
   TString pname = TString::Format("xsecratio");
   if(peak1 * peak2==0) 
     leg->Draw();
   else
     pname.Append(TString::Format("_%ds_%ds",peak1,peak2));
   pname += ".pdf";
   a.SaveAs(pname);
   return;
}
