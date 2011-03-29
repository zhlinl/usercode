#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TMath.h>
#include <TF1.h>
#include <iostream>
#endif
using namespace std;

void combine_accp_binned() {

  gROOT->Reset();

  TFile* f = new TFile("acceptance_filter.root");

//  f->cd("IdEfficiencyTnPrelMC");
  f->cd("");
  TH2F* ha;
  gDirectory->GetObject("acceptance2D_binned",ha);
  TH1F* hc;
  gDirectory->GetObject("acceptance1D_binned",hc);

  TFile* f2 = new TFile("acceptance_tracking.root");//acceptance_p1_mixed.root");
  f2->cd("");
  TH2F* hb;
  gDirectory->GetObject("acceptance2D_binned",hb);
  TH1F* hd;
  gDirectory->GetObject("acceptance1D_binned",hd);

  TH2F* accp_comb_binned = (TH2F*)hb->Clone("accp_comb_binned");
  TH1F* accp1d_comb = (TH1F*)hd->Clone("accp1d_comb");

  cout<<"binx"<<ha->GetNbinsX()<<endl;
  cout<<"biny"<<ha->GetNbinsY()<<endl; 
  for(int j=1;j<=17;j++){//pt
     accp1d_comb->SetBinContent(j,(hc->GetBinContent(j))*(hd->GetBinContent(j)));
     cout<<"j"<<j<<"hc"<<hc->GetBinContent(j)<<"hd"<<hd->GetBinContent(j)<<"content"<<hc->GetBinContent(j)*hd->GetBinContent(j)<<endl;
     accp1d_comb->SetBinError(j,sqrt( pow((hc->GetBinError(j))* (hd->GetBinContent(j)),2) + pow((hc->GetBinError(j))* (hd->GetBinContent(j)),2)));
  }

     for(int i=1;i<=17;i++){//pt
       for(int j=1;j<=9;j++){//eta
        accp_comb_binned->SetBinContent(i,j,(ha->GetBinContent(i,j))*(hb->GetBinContent(i,j)));
        cout<<"i"<<i<<"j"<<j<<"ha"<<ha->GetBinContent(i,j)<<"hb"<<hb->GetBinContent(i,j)<<"tol"<<(ha->GetBinContent(i,j))*(hb->GetBinContent(i,j))<<endl;
        accp_comb_binned->SetBinError(i,j,sqrt( pow((ha->GetBinError(i,j))* (hb->GetBinContent(i,j)),2) + pow((hb->GetBinError(i,j))* (ha->GetBinContent(i,j)),2)));
     }
  } 

  TFile *output =  new TFile("accp_comb_mixed.root","RECREATE");
  output->cd();
  accp_comb_binned->Write();
  accp1d_comb->Write();

  output->Write();
  output->Close();
} 
      
