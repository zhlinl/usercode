/*
determine the upsilon efficiencies in two ways:
- from MC truth (reco|gen, trig|reco)
- from multiplication of single muon efficiencies
- estimates the difference as a systematics source
 */

#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif // ! M_PI
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1D.h"
#include "TLorentzVector.h"

bool passAcc(double muPlusEta, double muPlusPt, double muMinusEta, double muMinusPt) {
//  cout<<"muPlusEta"<<muPlusEta<<endl;
  bool pass = true;
  pass &= ((fabs(muPlusEta)<1.6&&muPlusPt>3.5) || (fabs(muPlusEta)<2.4&&fabs(muPlusEta)>1.6&&muPlusPt>2.5));
  pass &= ((fabs(muMinusEta)<1.6&&muMinusPt>3.5) || (fabs(muMinusEta)<2.4&&fabs(muMinusEta)>1.6&&muMinusPt>2.5));
  return pass;
}

Float_t get(TH2F* h, Float_t pt, Float_t eta){
  int i = h->FindBin(pt,fabs(eta));//eta,pt);
//  cout<<"i"<<i<<endl;
  return h->GetBinContent(i);
}

void dimuoneff(TString inFile="testmceff_total2.root",TString outFile="mceff.root"){

  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.13);

  gStyle->SetPaintTextFormat("6.3f");


  TString title("(#epsilon_{#mu}#epsilon_{#mu} - #epsilon_{#mu#mu} ) / #epsilon_{#mu#mu}");
  TString xtitle("|y^{#Upsilon}|");
  TString ytitle("p_{T}^{#Upsilon}  (GeV/c)");

  TFile *f1 = TFile::Open("effMaps_Map_Upsilon_MC_MCTRUTH_MuonID.root");
  f1->cd();
  TH2F* reco_mu_tnp_eta_pt = (TH2F*)gROOT->FindObject("TH2F_tagAndProbe_TM_pt_abseta_mcTrue_cnt_eff");
  if(!reco_mu_tnp_eta_pt){ 
    cout<<"Could not access tnp map!"<<endl;
    return;
  }

  TFile *f2 = TFile::Open("effMaps_Map_Upsilon_MC_MCTRUTH_Trigger.root");
  f2->cd(); 
  TH2F* trigger_mu_tnp_eta_pt = (TH2F*)gROOT->FindObject("TH2F_tagAndProbe_Trigger_TM_pt_abseta_mcTrue_cnt_eff");
  if(!trigger_mu_tnp_eta_pt){
    cout<<"Could not access tnp map!"<<endl;
    return;
  }
  
  TFile infile(inFile);
  gDirectory->Cd("");
  TTree* tree_trackid = (TTree*)gROOT->FindObject("inTree");
  if(!tree_trackid){
    cout<<"Could not access yield tree!"<<endl;
    return;
  }
  

  Float_t genPosMuPt,  genPosMuEta,  genNegMuPt,  genNegMuEta, genUpsPt,  genUpsRap,  
    recoPosMuPt,  recoPosMuEta,  recoNegMuPt,  recoNegMuEta,  recoUpsPt,  recoUpsRap;
  Int_t HLT2MuonOpen; 
  inTree->SetBranchAddress("genUpsPt"    ,&genUpsPt);
  inTree->SetBranchAddress("genUpsRap"   ,&genUpsRap);
  inTree->SetBranchAddress("recoUpsPt"   ,&recoUpsPt);
  inTree->SetBranchAddress("recoUpsRap"  ,&recoUpsRap);
  inTree->SetBranchAddress("genPosMuPt"  ,&genPosMuPt);
  inTree->SetBranchAddress("genPosMuEta" ,&genPosMuEta);
  inTree->SetBranchAddress("genNegMuPt"  ,&genNegMuPt);
  inTree->SetBranchAddress("genNegMuEta" ,&genNegMuEta);
  inTree->SetBranchAddress("recoPosMuPt" ,&recoPosMuPt);
  inTree->SetBranchAddress("recoPosMuEta",&recoPosMuEta);
  inTree->SetBranchAddress("recoNegMuPt" ,&recoNegMuPt);
  inTree->SetBranchAddress("recoNegMuEta",&recoNegMuEta);
  inTree->SetBranchAddress("HLT2MuonOpen",&HLT2MuonOpen);


  //const int nyups   = 1;  double yBinsUps  [nyups  +1] = {0, 2};
  //const int nptups  = 1;  double ptBinsUps [nptups +1] = {0, 30};

  //const int nyups   = 5;  double yBinsUps  [nyups  +1] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0};
  //const int nptups  = 1;  double ptBinsUps [nptups +1] = {0, 30};

  //const int nyups   = 2;  double yBinsUps  [nyups  +1] = {0,1,2};
  //const int nptups  = 6;  double ptBinsUps [nptups +1] = {0, 2, 5, 8, 11, 15, 30};

  const int nyups   = 1;  double yBinsUps  [nyups  +1] = {0,2};
  const int nptups = 15;  double ptBinsUps [nptups +1] = {0,1,2,3,4,5,6,7,8,9,10,12,14,17,20,30};//1s
  //const int nptups =  8;  double ptBinsUps [nptups +1] = {0, 2, 4, 6, 9,12,16,20,30};//2S
  //const int nptups =  6;  double ptBinsUps [nptups +1] = {0, 3, 6, 9, 14, 20, 30};//3s

  const int netamuo = 6;  double etaBinsMuo[netamuo+1] = {0,0.4,0.8,1.2,1.6,2.0,2.4};
  const int nptmuo  = 8;  double ptBinsMuo  [nptmuo+1] = {2.5,3.0,3.5,4.0,4.5,5.0,6.0,8.0,50};

  double yBinsUps_[nyups+1];
  for(int i=0; i<nyups+1; i++) yBinsUps_[i] = yBinsUps[i]*1.2; //+ 0.1*(yBinsUps[i+1]-yBinsUps[i]);

  h_DiMuReco       = new TH2F("h_DiMuReco","h_DiMuReco",            nyups,yBinsUps,nptups,ptBinsUps);
  h_DiMuTrigger    = new TH2F("h_DiMuTrigger","h_DiMuTrigger",      nyups,yBinsUps,nptups,ptBinsUps);
  h_DiMuTrigReco   = new TH2F("h_DiMuTrigReco","h_DiMuTrigReco",    nyups,yBinsUps,nptups,ptBinsUps);
  recoEffRatio_tnp = new TH2F("recoEffRatio_tnp","recoEffRatio_tnp",nyups,yBinsUps,nptups,ptBinsUps);
  trigEffRatio_tnp = new TH2F("trigEffRatio_tnp","trigEffRatio_tnp",nyups,yBinsUps,nptups,ptBinsUps);
  trcoEffRatio_tnp = new TH2F("trcoEffRatio_tnp","trcoEffRatio_tnp",nyups,yBinsUps,nptups,ptBinsUps);
  //recoEffRatio_mu  = new TH2F("recoEffRatio_mu","recoEffRatio_mu",netamuo,etaBinsMuo,nptmuo,ptBinsMuo);
  recoEffRatio_mu  = new TH2F("recoEffRatio_mu","recoEffRatio_mu",  nptmuo,ptBinsMuo,netamuo,etaBinsMuo);

  recoEffRatio_err = new TH2F("recoEffRatio_err","recoEffRatio_err",nyups,yBinsUps_,nptups,ptBinsUps);
  trigEffRatio_err = new TH2F("trigEffRatio_err","trigEffRatio_err",nyups,yBinsUps_,nptups,ptBinsUps);
  trcoEffRatio_err = new TH2F("trcoEffRatio_err","trcoEffRatio_err",nyups,yBinsUps_,nptups,ptBinsUps);


  h_GenUpsPtRap    = new TH2F("GenUpsPtRap","GenUpsPtRap",          nyups,yBinsUps,nptups,ptBinsUps); //gen distr (for normalization of dimu eff)
  h_RecUpsPtRap    = new TH2F("RecUpsPtRap","RecUpsPtRap",          nyups,yBinsUps,nptups,ptBinsUps); //reco distr, for reco eff (given gen)
  h_TriUpsPtRap    = new TH2F("TriUpsPtRap","TriUpsPtRap",          nyups,yBinsUps,nptups,ptBinsUps); //reco+trig distr, for trig eff given reco
  h_TrcoUpsPtRap   = new TH2F("TrcoUpsPtRap","TrcoUpsPtRap",        nyups,yBinsUps,nptups,ptBinsUps);//reco+trig distr, for reco+trig (given gen)
  h_GenUpsPt       = new TH1F("GenUpsPt","GenUpsPt",                nptups,ptBinsUps);
  h_RecUpsPt       = new TH1F("RecUpsPt","RecUpsPt",                nptups,ptBinsUps);
  h_TriUpsPt       = new TH1F("TriUpsPt","TriUpsPt",                nptups,ptBinsUps);
  
  h_GenPosMuPt     = new TH1F("GenPosMuPt","GenPosMuPt",            nptmuo, ptBinsMuo);
//  h_GenPosMuPtEta  = new TH2F("GenPosMuPtEta","GenPosMuPtEta",      netamuo,etaBinsMuo,nptmuo,ptBinsMuo); 
  h_GenPosMuPtEta  = new TH2F("GenPosMuPtEta","GenPosMuPtEta",      nptmuo,ptBinsMuo,netamuo,etaBinsMuo);

  h_GenNegMuPt     = new TH1F("GenNegMuPt","GenNegMuPt",            nptmuo, ptBinsMuo);
//  h_GenNegMuPtEta  = new TH2F("GenNegMuPtEta","GenNegMuPtEta",      netamuo,etaBinsMuo,nptmuo,ptBinsMuo);
  h_GenNegMuPtEta  = new TH2F("GenNegMuPtEta","GenNegMuPtEta",      nptmuo,ptBinsMuo,netamuo,etaBinsMuo);

  h_PosMuPt        = new TH1F("PosMuPt","PosMuPt",                  nptmuo, ptBinsMuo);
//  h_PosMuPtEta     = new TH2F("PosMuPtEta","PosMuPtEta",            netamuo,etaBinsMuo,nptmuo,ptBinsMuo);
  h_PosMuPtEta     = new TH2F("PosMuPtEta","PosMuPtEta",            nptmuo,ptBinsMuo,netamuo,etaBinsMuo);

  
  h_NegMuPt        = new TH1F("NegMuPt","NegMuPt",                  nptmuo, ptBinsMuo);
//  h_NegMuPtEta     = new TH2F("NegMuPtEta","NegMuPtEta",            netamuo,etaBinsMuo,nptmuo,ptBinsMuo);
  h_NegMuPtEta     = new TH2F("NegMuPtEta","NegMuPtEta",            nptmuo,ptBinsMuo,netamuo,etaBinsMuo);
  
  /*
    if(hlt_L1DoubleMuOpen == 1){
    hltMuOpen_Y_pt->Fill(genUpsPt);
    hltMuOpen_Y_rap_pt->Fill(fabs(genUpsRap),genUpsPt);

    hltMuOpen_Y_rap_pt->Divide(gen_Y_rap_pt);
    effDiMuMuOpen = getCJ(hltMuOpen_mu_tnp_eta_pt,muPlusPt,muPlusEta)*getCJ(hltMuOpen_mu_tnp_eta_pt,muMinusPt,muMinusEta);

    h_DiMuMuOpen_tnp  ->Fill(upsRapidity,upsPt,effDiMuMuOpen);
    h_DiMuMuOpen_tnp  ->Add(hltMuOpen_Y_rap_pt,-1);
    MuOpeneffRatio_tnp->Divide(h_DiMuMuOpen_tnp,hltMuOpen_Y_rap_pt,1,1,"B");

  */

  h_DiMuReco->Sumw2();
  h_DiMuTrigger->Sumw2();
  h_DiMuTrigReco->Sumw2();
  recoEffRatio_tnp->Sumw2();
  trigEffRatio_tnp->Sumw2(); 
  trcoEffRatio_tnp->Sumw2(); 
  recoEffRatio_mu->Sumw2();
  h_GenPosMuPt->Sumw2();
  h_GenPosMuPtEta->Sumw2();
  h_GenNegMuPt->Sumw2();
  h_GenNegMuPtEta->Sumw2();
  h_GenUpsPt->Sumw2();
  h_GenUpsPtRap->Sumw2();
  h_PosMuPt->Sumw2();
  h_PosMuPtEta->Sumw2();
  h_NegMuPt->Sumw2();
  h_NegMuPtEta->Sumw2();
  h_RecUpsPt->Sumw2();
  h_RecUpsPtRap->Sumw2();
  h_TriUpsPt->Sumw2();
  h_TriUpsPtRap->Sumw2();
  h_TrcoUpsPtRap->Sumw2();
 
  //loop through the candidates
  for (Int_t i=0;i<inTree->GetEntries(); i++) {
    if(i%10000==0) std::cout<<i<<std::endl;
    inTree->GetEntry(i);

    //if(i==100 ) break;

    if(!passAcc(genPosMuEta,genPosMuPt,genNegMuEta,genNegMuPt)) continue;


    h_GenPosMuPt->Fill(genPosMuPt);
    h_GenPosMuPtEta->Fill(genPosMuPt,genPosMuEta);
    if(recoPosMuPt != 9999){
      h_PosMuPt->Fill(recoPosMuPt);
      h_PosMuPtEta->Fill(recoPosMuPt,recoPosMuEta);
    }
    h_GenNegMuPt->Fill(genNegMuPt);
    h_GenNegMuPtEta->Fill(genNegMuPt,genNegMuEta);
    if (recoNegMuPt != 9999){
      h_NegMuPt->Fill(recoNegMuPt);
      h_NegMuPtEta->Fill(recoNegMuPt,recoNegMuEta);
    }
    h_GenUpsPt->Fill(genUpsPt);
    h_GenUpsPtRap->Fill(fabs(genUpsRap),genUpsPt);
    if(recoNegMuPt != 9999 && recoPosMuPt != 9999){
      h_RecUpsPt->Fill(recoUpsPt);
      h_RecUpsPtRap->Fill(fabs(recoUpsRap),recoUpsPt);
      if(HLT2MuonOpen==1){
        h_TrcoUpsPtRap->Fill(fabs(recoUpsRap),recoUpsPt);
        h_TriUpsPtRap->Fill(fabs(recoUpsRap),recoUpsPt);
        h_TriUpsPt->Fill(recoUpsPt);
      }


      double effDiMuTnP = get(reco_mu_tnp_eta_pt,recoPosMuPt,recoPosMuEta) * get(reco_mu_tnp_eta_pt,recoNegMuPt,recoNegMuEta);
      double eff2MuOpen = get(trigger_mu_tnp_eta_pt,recoPosMuPt,recoPosMuEta) * get(trigger_mu_tnp_eta_pt,recoNegMuPt,recoNegMuEta);
      h_DiMuReco    ->Fill(fabs(recoUpsRap),recoUpsPt,effDiMuTnP);
      h_DiMuTrigger ->Fill(fabs(recoUpsRap),recoUpsPt,eff2MuOpen);
      h_DiMuTrigReco->Fill(fabs(recoUpsRap),recoUpsPt,eff2MuOpen*effDiMuTnP);
      
      if(i%10000==0)
	printf("  reco:%6.3f  trig:%6.3f   reco*trig:%6.3f \n", effDiMuTnP, eff2MuOpen,eff2MuOpen*effDiMuTnP);
      
    }


  }

  //close everything
  TFile outfile(outFile,"recreate");


  h_GenPosMuPt->Write();
  h_GenPosMuPtEta->Write();
  h_PosMuPt->Write();
  h_PosMuPtEta->Write();
  h_GenNegMuPt->Write();
  h_GenNegMuPtEta->Write();
  h_NegMuPt->Write();
  h_NegMuPtEta->Write();
  h_GenUpsPt->Write();
  h_GenUpsPtRap->Write();
  h_RecUpsPt->Write();
  h_RecUpsPtRap->Write();

  h_TriUpsPtRap->Divide(h_RecUpsPtRap);
  h_TriUpsPt->SetName("ups_2MuOpen_eff_pt");
  h_TriUpsPtRap->SetName("ups_2MuOpen_eff_pt_rap");
  h_TriUpsPtRap->Write();
  h_TriUpsPt->Write();

  // the dimuon reco efficiency from tnp
//  h_DiMuReco->Write();
//  h_DiMuTrigger->Write();
  h_DiMuReco    ->Divide(h_RecUpsPtRap);
  h_DiMuTrigger ->Divide(h_RecUpsPtRap);
  h_DiMuTrigReco->Divide(h_RecUpsPtRap);
  h_DiMuReco     ->Write();
  h_DiMuTrigger ->Write();
  h_DiMuTrigReco->Write();

  // reco eff vs pt
  h_RecUpsPt   ->Divide(h_GenUpsPt);

  // reco eff vs pt and rap
  h_RecUpsPtRap->Divide(h_GenUpsPtRap);

  // reco+trig eff vs pt eta
  h_TrcoUpsPtRap->Divide(h_GenUpsPtRap);

  h_RecUpsPt   ->SetName("ups_reco_eff_pt");
  h_RecUpsPtRap->SetName("ups_reco_eff_pt_rap");
  h_RecUpsPt   ->Write();
  h_RecUpsPtRap->Write();

  h_RecUpsPtRap ->SetTitle("#epsilon_{#mu#mu}   Muon Id | Gen");
  h_TriUpsPtRap ->SetTitle("#epsilon_{#mu#mu}   Trigger | Reco");  
  h_TrcoUpsPtRap->SetTitle("#epsilon_{#mu#mu}   Muon ID & Trigger | Gen");  
  h_RecUpsPtRap ->SetXTitle(xtitle);
  h_TriUpsPtRap ->SetXTitle(xtitle);
  h_TrcoUpsPtRap->SetXTitle(xtitle);
  h_RecUpsPtRap ->SetYTitle(ytitle);
  h_TriUpsPtRap ->SetYTitle(ytitle);
  h_TrcoUpsPtRap->SetYTitle(ytitle);

  h_DiMuReco    ->SetTitle("#epsilon_{#mu}#epsilon_{#mu}   Muon Id | Gen");
  h_DiMuTrigger ->SetTitle("#epsilon_{#mu}#epsilon_{#mu}   Trigger | Reco");  
  h_DiMuTrigReco->SetTitle("#epsilon_{#mu}#epsilon_{#mu}   Muon ID & Trigger | Gen");  
  h_DiMuReco    ->SetXTitle(xtitle);
  h_DiMuTrigger ->SetXTitle(xtitle);
  h_DiMuTrigReco->SetXTitle(xtitle);
  h_DiMuReco    ->SetYTitle(ytitle);
  h_DiMuTrigger ->SetYTitle(ytitle);
  h_DiMuTrigReco->SetYTitle(ytitle);

  //tmp, jus to display erros
  for(int i=1; i<recoEffRatio_tnp->GetNbinsX()+1; i++) {
    for(int j=1; j<recoEffRatio_tnp->GetNbinsY()+1; j++) { 
      recoEffRatio_err->SetBinContent(i,j,h_RecUpsPtRap ->GetBinError(i,j));
      trigEffRatio_err->SetBinContent(i,j,h_TriUpsPtRap ->GetBinError(i,j));
      trcoEffRatio_err->SetBinContent(i,j,h_TrcoUpsPtRap->GetBinError(i,j));
    }
  }
  h_RecUpsPtRap   ->SetMarkerSize(1.5);
  h_TriUpsPtRap   ->SetMarkerSize(1.5);
  h_TrcoUpsPtRap  ->SetMarkerSize(1.5);
  h_DiMuReco      ->SetMarkerSize(1.5);
  h_DiMuTrigger   ->SetMarkerSize(1.5);
  h_DiMuTrigReco  ->SetMarkerSize(1.5);
  recoEffRatio_err->SetMarkerSize(1.1);
  trigEffRatio_err->SetMarkerSize(1.1);
  trcoEffRatio_err->SetMarkerSize(1.1);
  recoEffRatio_err->SetMarkerColor(15);
  trigEffRatio_err->SetMarkerColor(15);
  trcoEffRatio_err->SetMarkerColor(15);

  TCanvas c1; h_RecUpsPtRap ->Draw("colz text"); recoEffRatio_err->Draw("text same"); c1.SaveAs("effDiMuReco.gif");
  TCanvas c2; h_TriUpsPtRap ->Draw("colz text"); trigEffRatio_err->Draw("text same"); c2.SaveAs("effDiMuTrig.gif");
  TCanvas c3; h_TrcoUpsPtRap->Draw("colz text"); trcoEffRatio_err->Draw("text same"); c3.SaveAs("effDiMuTotal.gif");

  //tmp, jus to display erros
  for(int i=1; i<recoEffRatio_tnp->GetNbinsX()+1; i++) {
    for(int j=1; j<recoEffRatio_tnp->GetNbinsY()+1; j++) { 
      recoEffRatio_err->SetBinContent(i,j,h_DiMuReco    ->GetBinError(i,j));
      trigEffRatio_err->SetBinContent(i,j,h_DiMuTrigger ->GetBinError(i,j));
      trcoEffRatio_err->SetBinContent(i,j,h_DiMuTrigReco->GetBinError(i,j));
    }
  }

  TCanvas c10; h_DiMuReco    ->Draw("colz text"); recoEffRatio_err->Draw("text same"); c10.SaveAs("effDiMuRecoTnp.gif");
  TCanvas c20; h_DiMuTrigger ->Draw("colz text"); trigEffRatio_err->Draw("text same"); c20.SaveAs("effDiMuTrigTnp.gif");
  TCanvas c30; h_DiMuTrigReco->Draw("colz text"); trcoEffRatio_err->Draw("text same"); c30.SaveAs("effDiMuTotalTnp.gif");

  // difference reco of dimuon effs, tnp - true
  h_DiMuReco->Add(h_RecUpsPtRap,-1);

  // difference divided by ratio, for dimuon eff tnp vs true
  recoEffRatio_tnp ->Divide(h_DiMuReco, h_RecUpsPtRap,   1,1,"B");
  recoEffRatio_tnp ->Write();  
  
  h_DiMuTrigger->Add(h_TriUpsPtRap,-1);
  //trigEffRatio_tnp ->Divide(h_DiMuTrigger,h_RecUpsPtRap,1,1,"B");//WRONG!!!
  trigEffRatio_tnp ->Divide(h_DiMuTrigger,h_TriUpsPtRap,1,1,"B");
  trigEffRatio_tnp->Write();

  h_DiMuTrigReco->Add(h_TrcoUpsPtRap,-1);
  trcoEffRatio_tnp ->Divide(h_DiMuTrigReco,h_TrcoUpsPtRap,1,1,"B");
  trcoEffRatio_tnp->Write();

  h_GenPosMuPt->Add(h_GenNegMuPt);
  h_PosMuPt->Add(h_NegMuPt);
  h_PosMuPt->Divide(h_GenPosMuPt);
  h_PosMuPt->SetName("mu_reco_eff_pt");

  h_GenPosMuPtEta->Add(h_GenNegMuPtEta);
  h_PosMuPtEta   ->Add(h_NegMuPtEta);
  h_PosMuPtEta->Divide(h_GenPosMuPtEta);
  h_PosMuPtEta->SetName("mu_reco_eff_pt_eta");
  h_PosMuPt->Write();
  h_PosMuPtEta->Write();
  
  cout<<"tnp bins"<<reco_mu_tnp_eta_pt->GetXaxis()->GetNbins()<<endl;
  cout<<"mc bins"<<h_PosMuPtEta->GetXaxis()->GetNbins()<<endl;
  reco_mu_tnp_eta_pt->Add(h_PosMuPtEta,-1);
  recoEffRatio_mu->Divide(reco_mu_tnp_eta_pt,h_PosMuPtEta,1,1,"B");
  recoEffRatio_mu->Write();

  //plot
  TString title("(#epsilon_{#mu}#epsilon_{#mu} - #epsilon_{#mu#mu} ) / #epsilon_{#mu#mu}");
  TString xtitle("|y^{#Upsilon}|");
  TString ytitle("p_{T}^{#Upsilon}  (GeV/c)");

  recoEffRatio_tnp->SetTitle(title + "   Muon Id | Gen");
  recoEffRatio_tnp->SetXTitle(xtitle);
  recoEffRatio_tnp->SetYTitle(ytitle);

  trigEffRatio_tnp->SetTitle(title + "   Trigger | Reco");  
  trigEffRatio_tnp->SetXTitle(xtitle);
  trigEffRatio_tnp->SetYTitle(ytitle);

  trcoEffRatio_tnp->SetTitle(title + "   Muon ID & Trigger | Gen");  
  trcoEffRatio_tnp->SetXTitle(xtitle);
  trcoEffRatio_tnp->SetYTitle(ytitle);

  recoEffRatio_tnp->SetMarkerSize(1.5);
  trigEffRatio_tnp->SetMarkerSize(1.5);
  trcoEffRatio_tnp->SetMarkerSize(1.5);
  recoEffRatio_err->SetMarkerSize(1.1);
  trigEffRatio_err->SetMarkerSize(1.1);
  trcoEffRatio_err->SetMarkerSize(1.1);
  recoEffRatio_err->SetMarkerColor(15);
  trigEffRatio_err->SetMarkerColor(15);
  trcoEffRatio_err->SetMarkerColor(15);

  for(int i=1; i<recoEffRatio_tnp->GetNbinsX()+1; i++) {
    for(int j=1; j<recoEffRatio_tnp->GetNbinsY()+1; j++) { 
      recoEffRatio_err->SetBinContent(i,j,recoEffRatio_tnp->GetBinError(i,j));
      trigEffRatio_err->SetBinContent(i,j,trigEffRatio_tnp->GetBinError(i,j));
      trcoEffRatio_err->SetBinContent(i,j,trcoEffRatio_tnp->GetBinError(i,j));
    }
  }

  TString nn("");//("_global");

  TCanvas r1;  
  gPad->SetLogy(); 
  //recoEffRatio_tnp->DrawCopy("colz texte");  
  recoEffRatio_tnp->Draw("text colz");  
  recoEffRatio_err->Draw("text same");  
  r1.SaveAs("recoEffRatio" + nn + ".gif"); r1.SaveAs("recoEffRatio" + nn + ".pdf");

  TCanvas r2;  
  gPad->SetLogy(); 
  //trigEffRatio_tnp->Draw("colz texte");  
  trigEffRatio_tnp->Draw("colz text");  
  trigEffRatio_err->Draw("text same");  
  r2.SaveAs("trigEffRatio" + nn + ".gif"); r2.SaveAs("trigEffRatio" + nn + ".pdf");

  TCanvas r3;  
  gPad->SetLogy(); 
  //trcoEffRatio_tnp->Draw("colz texte");  
  trcoEffRatio_tnp->Draw("colz text");  
  trcoEffRatio_err->Draw("text same");  
  r3.SaveAs("trcoEffRatio" + nn + ".gif"); r3.SaveAs("trcoEffRatio" + nn + ".pdf");


  /// dump values
  for(int irap=0; irap<nyups; irap++ ) {
    printf("%2.0f<|y|<%2.0f\n",yBinsUps[irap],yBinsUps[irap+1]);
      for(int ipt=0; ipt<nptups; ipt++ ) {
	printf("\t%2.0f<|pt|<%2.0f\t",
	       ptBinsUps[ipt],ptBinsUps[ipt+1]); 
	printf("muid:%6.3f+/-%5.3f \t",
	       recoEffRatio_tnp->GetBinContent(irap+1,ipt+1),
	       recoEffRatio_err->GetBinContent(irap+1,ipt+1));
	printf("trig:%6.3f+/-%5.3f\t",
	       trigEffRatio_tnp->GetBinContent(irap+1,ipt+1),
	       trigEffRatio_err->GetBinContent(irap+1,ipt+1));
	printf("reco+trig:%6.3f+/-%5.3f\n",
	       trcoEffRatio_tnp->GetBinContent(irap+1,ipt+1),
	       trcoEffRatio_err->GetBinContent(irap+1,ipt+1));
    }
  }



  outfile.Close();
  infile.Close();
}
