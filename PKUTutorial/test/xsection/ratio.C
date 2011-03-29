//double lumirelerr = 0.11; 

double* getptbinedges(int &nbin);

void ratio() {
  
  gROOT->LoadMacro("setTDRStyle_modified.C");
  setTDRStyle();

  TString path = "mode3/";
  TString pname = TString::Format("xsec_ratio");
  TString ytitle("#sigma #times #Beta (#mu#mu)  ratio");
  //TString ytitle = TString::Format("#sigma #times BR(#mu#mu) ratio");
  //TString ytitle = TString::Format("cross section ratio");
  TString xtitle("p_{T}^{#Upsilon} (GeV/c)");

  const int npeak = 3; //3 peaks

  //get pt bins
  const  int nbin;
  double* pte = getptbinedges(nbin);

  TString nfile[npeak];
  for(int i=0; i<npeak; i++)
    nfile[i] = TString::Format("mode3/xsection_%ds_y0.root",i+1);

  TFile* file[npeak];
  for(int i=0; i<npeak; i++)
    file[i] = TFile::Open(nfile[i],"read");
  
  //TGraphAsymmErrors *grall[npeak],*grsys[npeak],*grsta[npeak];
  TGraphAsymmErrors *gr[3][npeak];

  enum err {lo=0,hi};
  enum unc {all=0,sys,sta,cen};
  enum peak {ups1s=0,ups2s,ups3s};

  for(int i=0; i<npeak; i++) {
    gDirectory->Cd( TString::Format("%s:",nfile[i].Data()));
    //gDirectory->pwd();
    //gDirectory->ls();
    gr[all][i] = (TGraphAsymmErrors*)((TMultiGraph*)file[i]->Get("mg1"))->GetListOfGraphs()->At(all);
    gr[sys][i] = (TGraphAsymmErrors*)((TMultiGraph*)file[i]->Get("mg1"))->GetListOfGraphs()->At(sys);
    gr[sta][i] = (TGraphAsymmErrors*)((TMultiGraph*)file[i]->Get("mg1"))->GetListOfGraphs()->At(sta);
  }
    
  const int nb = gr[0][0]->GetN();
  assert(nb==nbin);

  const int npeak=3;//1s/2s/3s
  const int nunc =3+1; //central, total, syst, sta 
  double xsec[nbin][npeak][nunc][2];
  
  for(int i=0; i<nbin; i++) {
    for(int p=0; p<npeak; p++) {
      xsec[i][p][cen][0 ] =      (gr[all][p]->GetY())[i];      //central
      xsec[i][p][all][hi] = fabs((gr[all][p]->GetEYhigh())[i]);//total unc hi
      xsec[i][p][all][lo] = fabs((gr[all][p]->GetEYlow ())[i]);//total unc lo
      xsec[i][p][sys][hi] = fabs((gr[sys][p]->GetEYhigh())[i]);//systm unc hi
      xsec[i][p][sys][lo] = fabs((gr[sys][p]->GetEYlow ())[i]);//systm unc lo
      xsec[i][p][sta][hi] = fabs((gr[sta][p]->GetEYhigh())[i]);//stats unc hi
      xsec[i][p][sta][lo] = fabs((gr[sta][p]->GetEYlow ())[i]);//stats unc lo
      //remove stat from syst+stat
      xsec[i][p][sys][hi] = sqrt( pow((gr[sys][p]->GetEYhigh())[i],2) - pow((gr[sta][p]->GetEYhigh())[i],2) );//systm unc hi
      xsec[i][p][sys][lo] = sqrt( pow((gr[sys][p]->GetEYlow ())[i],2) - pow((gr[sta][p]->GetEYlow ())[i],2) );//systm unc lo
    }
  }
  
  double a_rat31[nbin], a_rat31EtotHi[nbin], a_rat31EtotLo[nbin], a_rat31EsysHi[nbin], a_rat31EsysLo[nbin], a_rat31EstaHi[nbin], a_rat31EstaLo[nbin];
  double a_rat21[nbin], a_rat21EtotHi[nbin], a_rat21EtotLo[nbin], a_rat21EsysHi[nbin], a_rat21EsysLo[nbin], a_rat21EstaHi[nbin], a_rat21EstaLo[nbin];
  double a_rat32[nbin], a_rat32EtotHi[nbin], a_rat32EtotLo[nbin], a_rat32EsysHi[nbin], a_rat32EsysLo[nbin], a_rat32EstaHi[nbin], a_rat32EstaLo[nbin];

  for(int i=0; i<nbin; i++) {

    //ratio ups3s / ups1s
    a_rat31       [i] = xsec[i][ups3s][cen][0] / xsec[i][ups1s][cen][0];
    a_rat31EtotHi [i] = -1*a_rat31[i] + (xsec[i][ups3s][cen][0]+xsec[i][ups3s][all][hi]) / (xsec[i][ups1s][cen][0]-xsec[i][ups1s][all][lo]);
    a_rat31EtotLo [i] =    a_rat31[i] - (xsec[i][ups3s][cen][0]-xsec[i][ups3s][all][lo]) / (xsec[i][ups1s][cen][0]+xsec[i][ups1s][all][hi]);
    a_rat31EsysHi [i] = -1*a_rat31[i] + (xsec[i][ups3s][cen][0]+xsec[i][ups3s][sys][hi]) / (xsec[i][ups1s][cen][0]-xsec[i][ups1s][sys][lo]);
    a_rat31EsysLo [i] =    a_rat31[i] - (xsec[i][ups3s][cen][0]-xsec[i][ups3s][sys][lo]) / (xsec[i][ups1s][cen][0]+xsec[i][ups1s][sys][hi]);
    a_rat31EstaHi [i] = -1*a_rat31[i] + (xsec[i][ups3s][cen][0]+xsec[i][ups3s][sta][hi]) / (xsec[i][ups1s][cen][0]-xsec[i][ups1s][sta][lo]);
    a_rat31EstaLo [i] =    a_rat31[i] - (xsec[i][ups3s][cen][0]-xsec[i][ups3s][sta][lo]) / (xsec[i][ups1s][cen][0]+xsec[i][ups1s][sta][hi]);
    
    //ratio ups2s / ups1s
    a_rat21       [i] = xsec[i][ups2s][cen][0] / xsec[i][ups1s][cen][0];
    a_rat21EtotHi [i] = -1*a_rat21[i] + (xsec[i][ups2s][cen][0]+xsec[i][ups2s][all][hi]) / (xsec[i][ups1s][cen][0]-xsec[i][ups1s][all][lo]);
    a_rat21EtotLo [i] =    a_rat21[i] - (xsec[i][ups2s][cen][0]-xsec[i][ups2s][all][lo]) / (xsec[i][ups1s][cen][0]+xsec[i][ups1s][all][hi]);
    a_rat21EsysHi [i] = -1*a_rat21[i] + (xsec[i][ups2s][cen][0]+xsec[i][ups2s][sys][hi]) / (xsec[i][ups1s][cen][0]-xsec[i][ups1s][sys][lo]);
    a_rat21EsysLo [i] =    a_rat21[i] - (xsec[i][ups2s][cen][0]-xsec[i][ups2s][sys][lo]) / (xsec[i][ups1s][cen][0]+xsec[i][ups1s][sys][hi]);
    a_rat21EstaHi [i] = -1*a_rat21[i] + (xsec[i][ups2s][cen][0]+xsec[i][ups2s][sta][hi]) / (xsec[i][ups1s][cen][0]-xsec[i][ups1s][sta][lo]);
    a_rat21EstaLo [i] =    a_rat21[i] - (xsec[i][ups2s][cen][0]-xsec[i][ups2s][sta][lo]) / (xsec[i][ups1s][cen][0]+xsec[i][ups1s][sta][hi]);

    //ratio ups3s / ups2s
    a_rat32       [i] = xsec[i][ups3s][cen][0] / xsec[i][ups2s][cen][0];
    a_rat32EtotHi [i] = -2*a_rat32[i] + (xsec[i][ups3s][cen][0]+xsec[i][ups3s][all][hi]) / (xsec[i][ups2s][cen][0]-xsec[i][ups2s][all][lo]);
    a_rat32EtotLo [i] =    a_rat32[i] - (xsec[i][ups3s][cen][0]-xsec[i][ups3s][all][lo]) / (xsec[i][ups2s][cen][0]+xsec[i][ups2s][all][hi]);
    a_rat32EsysHi [i] = -2*a_rat32[i] + (xsec[i][ups3s][cen][0]+xsec[i][ups3s][sys][hi]) / (xsec[i][ups2s][cen][0]-xsec[i][ups2s][sys][lo]);
    a_rat32EsysLo [i] =    a_rat32[i] - (xsec[i][ups3s][cen][0]-xsec[i][ups3s][sys][lo]) / (xsec[i][ups2s][cen][0]+xsec[i][ups2s][sys][hi]);
    a_rat32EstaHi [i] = -2*a_rat32[i] + (xsec[i][ups3s][cen][0]+xsec[i][ups3s][sta][hi]) / (xsec[i][ups2s][cen][0]-xsec[i][ups2s][sta][lo]);
    a_rat32EstaLo [i] =    a_rat32[i] - (xsec[i][ups3s][cen][0]-xsec[i][ups3s][sta][lo]) / (xsec[i][ups2s][cen][0]+xsec[i][ups2s][sta][hi]);
  }
  
  double *ptl[3], *pth[3];
  for(int k=0;k<3;k++) {
    ptl[k] = gr[0][k]->GetEXlow();
    pth[k] = gr[0][k]->GetEXhigh();
  }

  double ptn[nbin],pty[nbin];
  for(int i=0; i<nbin; i++) {
    ptn[i] = 0.2;
    pty[i] = 0.;
  }

  grall31 = new TGraphAsymmErrors(nbin, gr[0][ups3s]->GetX(), a_rat31, ptl[ups3s],pth[ups3s], a_rat31EtotLo, a_rat31EtotHi);
  grsys31 = new TGraphAsymmErrors(nbin, gr[0][ups3s]->GetX(), a_rat31, pty,pty, a_rat31EsysLo, a_rat31EsysHi);
  grsta31 = new TGraphAsymmErrors(nbin, gr[0][ups3s]->GetX(), a_rat31, ptn,ptn, a_rat31EstaLo, a_rat31EstaHi);

  grall21 = new TGraphAsymmErrors(nbin, gr[0][ups2s]->GetX(), a_rat21, ptl[ups2s],pth[ups2s], a_rat21EtotLo, a_rat21EtotHi);
  grsys21 = new TGraphAsymmErrors(nbin, gr[0][ups2s]->GetX(), a_rat21, pty,pty, a_rat21EsysLo, a_rat21EsysHi);
  grsta21 = new TGraphAsymmErrors(nbin, gr[0][ups2s]->GetX(), a_rat21, ptn,ptn, a_rat21EstaLo, a_rat21EstaHi);

  grall32 = new TGraphAsymmErrors(nbin, gr[0][ups3s]->GetX(), a_rat32, ptl[ups3s],pth[ups3s], a_rat32EtotLo, a_rat32EtotHi);
  grsys32 = new TGraphAsymmErrors(nbin, gr[0][ups3s]->GetX(), a_rat32, pty,pty, a_rat32EsysLo, a_rat32EsysHi);
  grsta32 = new TGraphAsymmErrors(nbin, gr[0][ups3s]->GetX(), a_rat32, ptn,ptn, a_rat32EstaLo, a_rat32EstaHi);

  TCanvas x;
  TH1F *hh = gPad->DrawFrame(0,0,30,0.8);

  hh->GetYaxis()->SetTitleOffset(1.1); //needed
  hh->SetXTitle(xtitle);
  hh->SetYTitle(ytitle);
  
  grall31->SetMarkerStyle(20);   grall21->SetMarkerStyle(25);   grall32->SetMarkerStyle(26); 
  grall31->SetMarkerSize(1.5);   grall21->SetMarkerSize(1.7);   grall32->SetMarkerSize(2); 
  grall31->SetMarkerColor(13);   grall21->SetMarkerColor(13);   grall32->SetMarkerColor(13); 


  grsta31->SetLineColor  (38);   grsta21->SetLineColor  (38);   grsta32->SetLineColor  (38); 
  grsta31->SetMarkerColor(38);   grsta21->SetMarkerColor(38);   grsta32->SetMarkerColor(38); 
  grsta31->SetFillColor  (38);   grsta21->SetFillColor  (40);   grsta32->SetFillColor  (38); 
  grsta31->SetMarkerStyle(21);   grsta21->SetMarkerStyle(21);   grsta32->SetMarkerStyle(21); 
  grsta31->SetMarkerSize(0.8);   grsta21->SetMarkerSize(0.8);   grsta32->SetMarkerSize(0.8); 
  //grsta31->SetFillStyle(3004);    	 grsta21->SetFillStyle(0);    	grsta32->SetFillStyle(3005);    

  grsys31->SetLineColor  (kRed); grsys21->SetLineColor  (kRed); grsys32->SetLineColor  (kRed);
  grsys31->SetMarkerColor(kRed); grsys21->SetMarkerColor(kRed);	grsys32->SetMarkerColor(kRed);
  grsys31->SetMarkerStyle(25);	 grsys21->SetMarkerStyle(25); 	grsys32->SetMarkerStyle(25); 
  grsys31->SetMarkerSize(0.8);	 grsys21->SetMarkerSize(0.8); 	grsys32->SetMarkerSize(0.8); 
  grsys31->SetFillStyle(0);    	 grsys21->SetFillStyle(0);    	grsys32->SetFillStyle(0);    
  grsys31->SetLineWidth  (3); grsys21->SetLineWidth  (3); grsys32->SetLineWidth  (3);

  grall31->Draw("e0");
  grsta31->Draw("e2");
  grsys31->Draw("[]");
  grall31->Draw("pxy");

  grall21->Draw("e0");
  grsta21->Draw("e2");
  grsys21->Draw("[]");
  grall21->Draw("pxy");
 
  grsys31->Draw("[]");
  grsys21->Draw("[]");

  //z (no end bars); x (no error bars), || and [] (only ending of error bars)

  //grsta32->Draw("e2");
  //grsys32->Draw("e2");
  //grall32->Draw("e0");

  TPaveText *ttp = new TPaveText(0.2,0.8,0.5,0.9,"brNDC");//tdr
  ttp->SetBorderSize(0);
  ttp->AddText("CMS,  #sqrt{s} = 7 TeV");
  ttp->AddText("L = 3 pb^{-1},  |y|<2");

  TLegend *uleg = new TLegend(0.7,0.2,0.9,0.3);	
  //uleg->AddEntry(grall31,"total unc.","P");
  uleg->AddEntry(grsta31,"stat. unc.","p");
  uleg->AddEntry(grsys31,"syst. unc.","p");
  uleg->Draw();

  TLegend *lleg = new TLegend(0.2,0.6,0.5,0.75);
  lleg->AddEntry(grall31,"#Upsilon(3S)/#Upsilon(1S)","lpe");
  lleg->AddEntry(grall21,"#Upsilon(2S)/#Upsilon(1S)","lpe");
  lleg->Draw();
  ttp->Draw();

  //mgu->Draw("ap");
  //mgu->GetXaxis()->SetRangeUser(0,30);
  //mgu->GetYaxis()->SetRangeUser(0,0.8);

  //x.SaveAs("aa.gif");
  x.SaveAs(path+pname + "_detailed.pdf");
  x.SaveAs(path+pname + "_detailed.gif");


  TGraphAsymmErrors *gr31, *gr32, *gr21;


  /*
    Note:
    - luminosity uncertainty set to zero on the input (produced by xsecresults.C)
    Tbd:
    - display information (esp. in table) re stat and syst errors separately
    -  need to save relevant info in advance to input graphs (ie in xsecresults.C)
   */

  double rat31[nbin], rat31_el[nbin], rat31_eh[nbin];
  double rat32[nbin], rat32_el[nbin], rat32_eh[nbin];
  double rat21[nbin], rat21_el[nbin], rat21_eh[nbin];
  double val1, err1l,err1h, val2,err2l,err2h, val3,err3l,err3h;
  for(int i=0; i<nbin; i++) {
    val1  = (gr[0][0]->GetY())[i];
    val2  = (gr[0][1]->GetY())[i];
    val3  = (gr[0][2]->GetY())[i];
    err1l = fabs((gr[0][0]->GetEYlow())[i] );
    err2l = fabs((gr[0][1]->GetEYlow())[i] );
    err3l = fabs((gr[0][2]->GetEYlow())[i] );
    err1h = fabs((gr[0][0]->GetEYhigh())[i]);
    err2h = fabs((gr[0][1]->GetEYhigh())[i]);
    err3h = fabs((gr[0][2]->GetEYhigh())[i]);
    rat31  [i] = val3/val1;
    rat32  [i] = val3/val2;
    rat21  [i] = val2/val1;
    rat31_el[i] =    rat31[i] - (val3-err3l)/(val1+err1h);
    rat32_el[i] =    rat32[i] - (val3-err3l)/(val2+err2h);
    rat21_el[i] =    rat21[i] - (val2-err2l)/(val1+err1h);
    rat31_eh[i] = -1*rat31[i] + (val3+err3h)/(val1-err1l);
    rat32_eh[i] = -1*rat32[i] + (val3+err3h)/(val2-err2l);
    rat21_eh[i] = -1*rat21[i] + (val2+err2h)/(val1-err1l);
    double xx_rat31_el = rat31[i]*sqrt(pow(err3l/val3,2)+pow(err1l/val1,2));
    double xx_rat32_el = rat32[i]*sqrt(pow(err3l/val3,2)+pow(err2l/val2,2));
    double xx_rat21_el = rat21[i]*sqrt(pow(err2l/val2,2)+pow(err1l/val1,2));
    double xx_rat31_eh = rat31[i]*sqrt(pow(err3h/val3,2)+pow(err1h/val1,2));
    double xx_rat32_eh = rat32[i]*sqrt(pow(err3h/val3,2)+pow(err2h/val2,2));
    double xx_rat21_eh = rat21[i]*sqrt(pow(err2h/val2,2)+pow(err1h/val1,2));
    printf("ratio-31:  %4.2f err(+%4.2f-%4.2f) err'(+%4.2f-%4.2f)\n",
	   rat31[i], rat31_eh[i],rat31_el[i], xx_rat31_eh, xx_rat31_el);
  }

  char ss[50];
  ofstream flatex, fdetai;
  flatex.open(path+pname+".tex");
  fdetai.open(path+pname+"_detailed.tex");
  sprintf(ss,"\\begin{tabular}{ccc} \\hline \n");
  flatex << ss;
  fdetai << ss;
  sprintf(ss,"\\pt (\\GeVc) & \\upsiii/\\upsi & \\upsii/\\upsi \\\\ \\hline\n");
  flatex << ss;
  fdetai << ss;
  double* aptl = gr[0][0]->GetEXlow();
  double* apth = gr[0][0]->GetEXhigh();
  //double cen=-1;
  double pt1,pt2;
  for(int i=0; i<nbin; i++) {
    //cen =  gr[0]->GetX()[i];
    if(i==0) {
      pt1 = pte[0];
      pt2 = pte[nbin-1];
    } else {
      pt1 = pte[i-1];
      pt2 = pte[i];
    }   
    sprintf(ss,"$%3.0f - %3.0f$ & $%4.2f \\pm %3.2f\\,(%3.2f)$ & $%4.2f \\pm %3.2f\\,(%3.2f)$ \\\\ \n", 
	    pt1, pt2, 
	    rat31[i], rat31_eh[i], rat31_el[i],
	    //rat32[i], rat32_eh[i], rat32_el[i], 
	    rat21[i], rat21_eh[i], rat21_el[i]); 
    flatex << ss;
    //    sprintf(ss,"$%3.0f - %3.0f$ & $%4.2f \\pm %3.2f\\,(%3.2f) \\pm %3.2f\\,(%3.2f)$ & $%4.2f \\pm %3.2f\\,(%3.2f) \\pm %3.2f\\,(%3.2f)$ \\\\ \n", 
    // sprintf(ss,"$%3.0f - %3.0f$ & $%4.2f {}^{+%3.2f}_{-%3.2f} {}^{%3.2f}_{%3.2f}$ & $%4.2f {}^{+%3.2f}_{-%3.2f} {}^{+%3.2f}_{-%3.2f}$ \\\\ \n", 
    sprintf(ss,"$%3.0f - %3.0f$ & $%4.2f \\pm^{%3.2f}_{%3.2f} \\pm^{%3.2f}_{%3.2f}$ & $%4.2f \\pm^{%3.2f}_{%3.2f} \\pm^{%3.2f}_{%3.2f}$ \\\\ \n", 
	    pt1, pt2, 
	    a_rat31[i], a_rat31EstaHi[i], a_rat31EstaLo[i],a_rat31EsysHi[i], a_rat31EsysLo[i],
	    a_rat21[i], a_rat21EstaHi[i], a_rat21EstaLo[i],a_rat21EsysHi[i], a_rat21EsysLo[i]);
    fdetai << ss;
    
    //cout << "pt:" << cen << " " << ss << endl;
  }
  sprintf(ss,"\\hline \n \\end{tabular}\n");
  flatex << ss; 
  fdetai << ss; 
  flatex.close();
  fdetai.close();

  gr31 = new TGraphAsymmErrors(nbin, gr[0][0]->GetX(), rat31, aptl,apth, rat31_el, rat31_eh);
  gr32 = new TGraphAsymmErrors(nbin, gr[0][0]->GetX(), rat32, aptl,apth, rat32_el, rat32_eh);
  gr21 = new TGraphAsymmErrors(nbin, gr[0][0]->GetX(), rat21, aptl,apth, rat21_el, rat21_eh);

  TMultiGraph *mg = new TMultiGraph();
  //if      (peak1 == 3 && peak2 == 1) mg->Add(gr31);
  //else if (peak1 == 3 && peak2 == 2) mg->Add(gr32);
  //else if (peak1 == 2 && peak2 == 1) mg->Add(gr21);
  //else 
  {mg->Add(gr31); mg->Add(gr21);} //mg->Add(gr32); 


  //mg->SetMinimum(0);
  //mg->SetMaximum(0.8);


  //TString ytitle = TString::Format("(d#sigma/dp_{T}.BR(#Upsilon(%dS)) / (d#sigma/dp_{T}.BR(#Upsilon(%dS))",peak1,peak2);
  mg->SetTitle( "" );

  gr31->SetMarkerColor(13);
//  gr31->SetMarkerColor(2);
//  gr31->SetLineColor  (2);
//  //gr31->SetLineWidth(2);
//  gr32->SetMarkerColor(6);
//  gr32->SetLineColor  (6);
//  //gr32->SetLineWidth(2);
//  gr21->SetMarkerColor(4);
//  gr21->SetLineColor  (4);
//  //gr21->SetLineWidth(2);
  gr31->SetMarkerStyle(20);//24
  gr32->SetMarkerStyle(26);
  gr21->SetMarkerStyle(25);
  gr31->SetMarkerSize(1.4);
  gr32->SetMarkerSize(1.4);
  gr21->SetMarkerSize(1.4);

  //TLegend *leg = new TLegend(0.2,0.55,0.5,0.75);a
   TLegend *leg = new TLegend(0.2,0.6,0.5,0.75);
   leg->AddEntry(gr31,"#Upsilon(3S)/#Upsilon(1S)","lp");
   //leg->AddEntry(gr32,"#Upsilon(3S)/#Upsilon(2S)","lp");
   leg->AddEntry(gr21,"#Upsilon(2S)/#Upsilon(1S)","lp");

   //TPaveText *tp = new TPaveText(0.2,0.7,0.45,0.85,"brNDC");
   TPaveText *tp = new TPaveText(0.2,0.8,0.5,0.9,"brNDC");//tdr
   tp->SetBorderSize(0);
   tp->AddText("CMS,  #sqrt{s} = 7 TeV");
   tp->AddText("L = 3 pb^{-1},  |y|<2");
   //tp->AddText("L=3pb^{-1}");

   TCanvas a;
   //mg->SetMaximum(1);
   mg->Draw("ap");

   TAxis* ax = mg->GetXaxis(); 
   ax->SetTitle(xtitle);
   mg->GetYaxis()->SetTitle(ytitle);
   //mg->GetYaxis()->SetTitleSize(0.05);
   //mg->GetXaxis()->SetTitleSize(0.05);
   mg->GetXaxis()->SetRangeUser(0,30);
   mg->GetYaxis()->SetRangeUser(0,0.8);

   tp->Draw();
   //if(peak1 * peak2==0) 
   leg->Draw();
   //else
     //pname.Append(TString::Format("_%ds_%ds",peak1,peak2));
   a.SaveAs(path+pname + ".pdf");
   a.SaveAs(path+pname + ".gif");
   
   for(int i=0; i<npeak; i++)
     file[i]->Close();

   return;
}


double* getptbinedges(int &nbin) {
  TString fnvar = TString::Format("mode3/fitres_1s_upsilonYieldWeighted_nominal/xsection.root");
  TFile fvar; fvar.Open(fnvar,"READ");
  gDirectory->Cd(fnvar+":/plots");
  difvar = (TH2F*)gROOT->FindObject(TString::Format("pt_y_bins"));
  double* diff = difvar->GetXaxis()->GetXbins()->GetArray();
  nbin = difvar->GetXaxis()->GetNbins() + 1;
  fvar.Close();
  for(int i=0;i<nbin;i++)  printf("pt-edges[%d]:%3.1f\t",i,diff[i]); printf("\nnbin:%d\n",nbin);
  return diff;
}
