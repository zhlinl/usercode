#ifndef __CINT__
#include <RooGlobalFunc.h>
#endif                    

#include <TLatex.h>       
#include <RooRealVar.h>   
#include <RooDataSet.h>   
#include <RooGaussian.h>  
#include <RooFitResult.h> 
#include <RooLandau.h>    
#include <RooChebychev.h> 
#include <RooAddPdf.h>    
#include <RooPlot.h>      
#include <RooDataHist.h>  
#include <RooVoigtian.h>  
#include <RooCBShape.h>   
#include <TCanvas.h>      
#include <TROOT.h>        
#include <TStyle.h>        
#include <TLorentzVector.h>        

#include <TAxis.h>        
#include <TH1.h>          
#include <TH1F.h>          
#include <TH2F.h>          
#include <TTree.h>        
#include <TFile.h>        
#include <TH1D.h>         
#include <TH1I.h>         
#include <TCanvas.h>      
#include <TLine.h>        
#include <TMath.h>        
#include <TVector3.h>     
#include <TString.h>      
#include <TLegend.h> 
#include <TPad.h> 
#include "myStyle.h"

void effTrackMu()
{

	gROOT->Reset();
	Long_t nEntries;  
	int passedEvt=0;

	Double_t massMin,massMax;  
	massMin=2.6;               
	massMax=3.5;               

	char fileName[500];                                     
	sprintf(fileName,                                       
			"/home/zhlinl/data/muonTree_Run2010B-Nov4ReReco_v1-Onia2MuMu-v_new/muonTree_nocut_1.root"
			//"/home/zhlinl/cutAnalysis/test/macro/reducedTree.root"
			);

	Double_t JpsiMass,JpsiPt,JpsiRap;
	int HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2;
	//TClonesArray* muPosP=new TClonesArray("TLorentzVector");
	//TClonesArray* muNegP=new TClonesArray("TLorentzVector");
	//TClonesArray* JpsiP=new TClonesArray("TLorentzVector");
	TLorentzVector *muPosP, *muNegP, *JpsiP;
	int muPos_arbitrated, muPos_oneStationTight, muPos_lastStationAngTight,
			muPos_lastStationTight, muPos_lastStationOptimizedLowPtTight,
			muPos_lastStationOptimizedBarrelLowPtTight,muPos_oneStationAngTight;                                    
	//int muPos_found, muPos_pixeLayers, muPos_nValidMuHits;
	int muNeg_arbitrated, muNeg_oneStationTight, muNeg_lastStationAngTight,
			muNeg_lastStationTight, muNeg_lastStationOptimizedLowPtTight, 
			muNeg_lastStationOptimizedBarrelLowPtTight, muNeg_oneStationAngTight;                                   
	//int muNeg_found, muNeg_pixeLayers, muNeg_nValidMuHits;                                                      

	TFile *infile=new TFile(fileName,"R");        
	TTree *tree=(TTree*)infile->Get("data");   
	tree->SetBranchAddress("muPosP",&muPosP);
	tree->SetBranchAddress("muNegP",&muNegP);
	tree->SetBranchAddress("JpsiP",&JpsiP);
	tree->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2",&HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2);

	tree->SetBranchAddress("muPos_arbitrated",&muPos_arbitrated);
	tree->SetBranchAddress("muPos_oneStationTight", &muPos_oneStationTight);                                    
	tree->SetBranchAddress("muPos_lastStationAngTight",&muPos_lastStationAngTight); 
	tree->SetBranchAddress("muPos_lastStationTight",&muPos_lastStationTight);  
	tree->SetBranchAddress("muPos_lastStationOptimizedLowPtTight",&muPos_lastStationOptimizedLowPtTight);   
	tree->SetBranchAddress("muPos_lastStationOptimizedBarrelLowPtTight",&muPos_lastStationOptimizedBarrelLowPtTight);
	tree->SetBranchAddress("muPos_oneStationAngTight",&muPos_oneStationAngTight);   

	tree->SetBranchAddress("muNeg_arbitrated",&muNeg_arbitrated);
	tree->SetBranchAddress("muNeg_oneStationTight", &muNeg_oneStationTight);                                    
	tree->SetBranchAddress("muNeg_lastStationAngTight",&muNeg_lastStationAngTight);                            
	tree->SetBranchAddress("muNeg_lastStationTight",&muNeg_lastStationTight);                                   
	tree->SetBranchAddress("muNeg_lastStationOptimizedLowPtTight",&muNeg_lastStationOptimizedLowPtTight);    
	tree->SetBranchAddress("muNeg_lastStationOptimizedBarrelLowPtTight",&muNeg_lastStationOptimizedBarrelLowPtTight);
	tree->SetBranchAddress("muNeg_oneStationAngTight",&muNeg_oneStationAngTight);                               

	const int NBinPt=10;
	const int NBinEta=26;
	double maxPt=10;
	TH2F *histMu_pt_eta_raw=new TH2F("histMu_pt_eta_Raw","Mu_pt_eta_Raw",NBinEta,-2.6,2.6,NBinPt,0,maxPt);
	histMu_pt_eta_raw->SetXTitle("#eta");
	histMu_pt_eta_raw->SetYTitle("P_{T} (GeV/c)");
	TH2F *histMu_pt_eta_arbitr=new TH2F("histMu_pt_eta_arbitr","Mu_pt_eta_arbitr",NBinEta,-2.6,2.6,NBinPt,0,maxPt);
	histMu_pt_eta_arbitr->SetXTitle("#eta");
	histMu_pt_eta_arbitr->SetYTitle("P_{T} (GeV/c)");
	TH2F *histMu_pt_eta_oneSta=new TH2F("histMu_pt_eta_oneSta","Mu_pt_eta_oneSta",NBinEta,-2.6,2.6,NBinPt,0,maxPt);
	histMu_pt_eta_oneSta->SetXTitle("#eta");
	histMu_pt_eta_oneSta->SetYTitle("P_{T} (GeV/c)");
	TH2F *histMu_pt_eta_lastSta=new TH2F("histMu_pt_eta_lastSta","Mu_pt_eta_lastSta",NBinEta,-2.6,2.6,NBinPt,0,maxPt);
	histMu_pt_eta_lastSta->SetXTitle("#eta");
	histMu_pt_eta_lastSta->SetYTitle("P_{T} (GeV/c)");

	TH2F *effMu_pt_eta_arbitr=new TH2F("effMu_pt_eta_arbitr","effMu_pt_eta_TMArbitrated",NBinEta,-2.6,2.6,NBinPt,0,maxPt);
	effMu_pt_eta_arbitr->SetXTitle("#eta");
	effMu_pt_eta_arbitr->SetYTitle("P_{T} (GeV/c)");
	TH2F *effMu_pt_eta_oneSta=new TH2F("effMu_pt_eta_oneSta","effMu_pt_eta_OneStationTight",NBinEta,-2.6,2.6,NBinPt,0,maxPt);
	effMu_pt_eta_oneSta->SetXTitle("#eta");
	effMu_pt_eta_oneSta->SetYTitle("P_{T} (GeV/c)");
	TH2F *effMu_pt_eta_lastSta=new TH2F("effMu_pt_eta_lastSta","effMu_pt_eta_LastStationAngTight",NBinEta,-2.6,2.6,NBinPt,0,maxPt);
	effMu_pt_eta_lastSta->SetXTitle("#eta");
	effMu_pt_eta_lastSta->SetYTitle("P_{T} (GeV/c)");

	TH1F *histMu_pt_raw=new TH1F("histMu_pt_raw","Mu_pt_raw",NBinPt,0,maxPt);
	histMu_pt_raw->SetXTitle("P_{T} (GeV/c)");
	histMu_pt_raw->SetYTitle("Events");
	TH1F *histMu_pt_arbitr=new TH1F("histMu_pt_arbitr","Mu_pt_arbitr",NBinPt,0,maxPt);
	histMu_pt_arbitr->SetXTitle("P_{T} (GeV/c)");
	histMu_pt_arbitr->SetYTitle("Events");
	TH1F *histMu_pt_oneSta=new TH1F("histMu_pt_oneSta","Mu_pt_oneSta",NBinPt,0,maxPt);
	histMu_pt_oneSta->SetXTitle("P_{T} (GeV/c)");
	histMu_pt_oneSta->SetYTitle("Events");
	TH1F *histMu_pt_lastSta=new TH1F("histMu_pt_lastSta","Mu_pt_lastSta",NBinPt,0,maxPt);
	histMu_pt_lastSta->SetXTitle("P_{T} (GeV/c)");
	histMu_pt_lastSta->SetYTitle("Events");

	TH1F *effMu_pt_arbitr=new TH1F("effMu_pt_arbitr","effMu_pt_TMArbitrated",NBinPt,0,maxPt);
	effMu_pt_arbitr->SetXTitle("P_{T} (GeV/c)");
	effMu_pt_arbitr->SetYTitle("efficiency");
	TH1F *effMu_pt_oneSta=new TH1F("effMu_pt_oneSta","effMu_pt_OneStationTight",NBinPt,0,maxPt);
	effMu_pt_oneSta->SetXTitle("P_{T} (GeV/c)");
	effMu_pt_oneSta->SetYTitle("efficiency");
	TH1F *effMu_pt_lastSta=new TH1F("effMu_pt_lastSta","effMu_pt_LastStationAngTight",NBinPt,0,maxPt);
	effMu_pt_lastSta->SetXTitle("P_{T} (GeV/c)");
	effMu_pt_lastSta->SetYTitle("efficiency");

	TH1F *histMu_eta_raw=new TH1F("histMu_eta_raw","Mu_eta_raw",NBinEta,-2.6,2.6);
	histMu_eta_raw->SetXTitle("#eta");//Rapidity");
	histMu_eta_raw->SetYTitle("Events");
	TH1F *histMu_eta_arbitr=new TH1F("histMu_eta_arbitr","Mu_eta_arbitr",NBinEta,-2.6,2.6);
	histMu_eta_arbitr->SetXTitle("#eta");//Rapidity");
	histMu_eta_arbitr->SetYTitle("Events");
	TH1F *histMu_eta_oneSta=new TH1F("histMu_eta_oneSta","Mu_eta_oneSta",NBinEta,-2.6,2.6);
	histMu_eta_oneSta->SetXTitle("#eta");//Rapidity");
	histMu_eta_oneSta->SetYTitle("Events");
	TH1F *histMu_eta_lastSta=new TH1F("histMu_eta_lastSta","Mu_eta_lastSta",NBinEta,-2.6,2.6);
	histMu_eta_lastSta->SetXTitle("#eta");//Rapidity");
	histMu_eta_lastSta->SetYTitle("Events");
	
	TH1F *effMu_eta_arbitr=new TH1F("effMu_eta_arbitr","effMu_eta_TMArbitrated",NBinEta,-2.6,2.6);
	effMu_eta_arbitr->SetXTitle("#eta");
	effMu_eta_arbitr->SetYTitle("efficiency");
	TH1F *effMu_eta_oneSta=new TH1F("effMu_eta_oneSta","effMu_eta_OneStationTight",NBinEta,-2.6,2.6);
	effMu_eta_oneSta->SetXTitle("#eta");
	effMu_eta_oneSta->SetYTitle("efficiency");
	TH1F *effMu_eta_lastSta=new TH1F("effMu_eta_lastSta","effMu_eta_LastStationAngTight",NBinEta,-2.6,2.6);
	effMu_eta_lastSta->SetXTitle("#eta");
	effMu_eta_lastSta->SetYTitle("efficiency");

	nEntries=tree->GetEntries();

	cout<<"===============total Entries in Tree: "<<nEntries<<"==============="<<endl;
	for(int i=0; i<nEntries/1; i++)
	{
		tree->GetEntry(i);
		if((i<=100 && i%20==0)
				||(i>100 && i<=1000 && i%100==0)
				||(i>1000 && i%1000==0) )
		{   
			cout<<"entryID "<<i<<" in Total Entries "<<nEntries<<endl;
		}
		JpsiMass=JpsiP->M();
		if( 1
				//&& JpsiMass>massMin && JpsiMass<massMax
				&& HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2 ==1
				//&& muPos_arbitrated ==1 && muNeg_arbitrated ==1
				//&& cut 
			)
		{
			histMu_pt_eta_raw->Fill(muPosP->Eta(),muPosP->Pt());
			histMu_pt_eta_raw->Fill(muNegP->Eta(),muNegP->Pt());
			histMu_pt_raw->Fill(muPosP->Pt());
			histMu_pt_raw->Fill(muNegP->Pt());
			histMu_eta_raw->Fill(muPosP->Eta());
			histMu_eta_raw->Fill(muNegP->Eta());
			if(muPos_arbitrated ==1 && muNeg_arbitrated ==1)
			{
				histMu_pt_eta_arbitr->Fill(muPosP->Eta(),muPosP->Pt());
				histMu_pt_eta_arbitr->Fill(muNegP->Eta(),muNegP->Pt());
				histMu_pt_arbitr->Fill(muPosP->Pt());
				histMu_pt_arbitr->Fill(muNegP->Pt());
				histMu_eta_arbitr->Fill(muPosP->Eta());
				histMu_eta_arbitr->Fill(muNegP->Eta());
				passedEvt++;
			}
			if(muPos_arbitrated ==1 && muNeg_arbitrated ==1
					&&muPos_oneStationTight ==1 && muNeg_oneStationTight ==1 )
			{
			histMu_pt_eta_oneSta->Fill(muPosP->Eta(),muPosP->Pt());
			histMu_pt_eta_oneSta->Fill(muNegP->Eta(),muNegP->Pt());
				histMu_pt_oneSta->Fill(muPosP->Pt());
				histMu_pt_oneSta->Fill(muNegP->Pt());
				histMu_eta_oneSta->Fill(muPosP->Eta());
				histMu_eta_oneSta->Fill(muNegP->Eta());
				//passedEvt++;
			}
			if(muPos_arbitrated ==1 && muNeg_arbitrated ==1
					&& muPos_lastStationAngTight ==1 && muNeg_lastStationAngTight ==1)
			{
			histMu_pt_eta_lastSta->Fill(muPosP->Eta(),muPosP->Pt());
			histMu_pt_eta_lastSta->Fill(muNegP->Eta(),muNegP->Pt());
				histMu_pt_lastSta->Fill(muPosP->Pt());
				histMu_pt_lastSta->Fill(muNegP->Pt());
				histMu_eta_lastSta->Fill(muPosP->Eta());
				histMu_eta_lastSta->Fill(muNegP->Eta());
			}

		}
	}
	cout<<"========totally passed Events: "<<passedEvt<<"============"<<endl;

	effMu_pt_eta_arbitr->Divide(histMu_pt_eta_arbitr,histMu_pt_eta_raw,1,1,"B");
	effMu_pt_arbitr->Divide(histMu_pt_arbitr,histMu_pt_raw,1,1,"B");
	effMu_eta_arbitr->Divide(histMu_eta_arbitr,histMu_eta_raw,1,1,"B");
	bool select=1;
	if(select==1)
	{
		effMu_pt_eta_oneSta->Divide(histMu_pt_eta_oneSta,histMu_pt_eta_arbitr,1,1,"B");
		effMu_pt_eta_lastSta->Divide(histMu_pt_eta_lastSta,histMu_pt_eta_arbitr,1,1,"B");
		effMu_pt_oneSta->Divide(histMu_pt_oneSta,histMu_pt_arbitr,1,1,"B");
		effMu_pt_lastSta->Divide(histMu_pt_lastSta,histMu_pt_arbitr,1,1,"B");
		effMu_eta_oneSta->Divide(histMu_eta_oneSta,histMu_eta_arbitr,1,1,"B");
		effMu_eta_lastSta->Divide(histMu_eta_lastSta,histMu_eta_arbitr,1,1,"B");
	}

	if(select==0)
	{
		effMu_pt_eta_oneSta->Divide(histMu_pt_eta_oneSta,histMu_pt_eta_raw,1,1,"B");
		effMu_pt_eta_lastSta->Divide(histMu_pt_eta_lastSta,histMu_pt_eta_raw,1,1,"B");
		effMu_pt_oneSta->Divide(histMu_pt_oneSta,histMu_pt_raw,1,1,"B");
		effMu_pt_lastSta->Divide(histMu_pt_lastSta,histMu_pt_raw,1,1,"B");
		effMu_eta_oneSta->Divide(histMu_eta_oneSta,histMu_eta_raw,1,1,"B");
		effMu_eta_lastSta->Divide(histMu_eta_lastSta,histMu_eta_raw,1,1,"B");
	}

	gStyle->SetPalette(1);
	TCanvas *cv00=new TCanvas("cv00","cv00",100,100,900,500);
	SetMyStyle();
	gStyle->SetOptStat(0);	
	//effMu_pt_eta_arbitr->SetMarkerStyle(22);
	effMu_pt_eta_arbitr->GetZaxis()->SetRangeUser(0.,1.);
	effMu_pt_eta_arbitr->Draw("colz");
	cv00->Modified();
	TCanvas *cv01=new TCanvas("cv01","cv01",100,100,900,500);
	SetMyStyle();
	gStyle->SetOptStat(0);	
	//effMu_pt_eta_oneSta->SetMarkerStyle(22);
	effMu_pt_eta_oneSta->GetZaxis()->SetRangeUser(0.,1.);
	effMu_pt_eta_oneSta->Draw("colz");
	cv01->Modified();
	TCanvas *cv02=new TCanvas("cv02","cv02",100,100,900,500);
	SetMyStyle();
	gStyle->SetOptStat(0);	
	//effMu_pt_eta_lastSta->SetMarkerStyle(22);
	effMu_pt_eta_lastSta->GetZaxis()->SetRangeUser(0.,1.);
	effMu_pt_eta_lastSta->Draw("colz");
	cv02->Modified();
	cv00->SaveAs("effMu_pt_eta_arbitr.eps");
	cv01->SaveAs("effMu_pt_eta_oneSta.eps");
	cv02->SaveAs("effMu_pt_eta_lastSta.eps");
	cv00->SaveAs("effMu_pt_eta_arbitr.gif");
	cv01->SaveAs("effMu_pt_eta_oneSta.gif");
	cv02->SaveAs("effMu_pt_eta_lastSta.gif");

	TCanvas *cv10=new TCanvas("cv10","cv10",100,100,900,500);
	SetMyStyle();
	gStyle->SetOptStat(0);	
	effMu_pt_arbitr->SetMarkerStyle(22);
	effMu_pt_arbitr->SetMarkerColor(kBlue);
	effMu_pt_arbitr->SetLineColor(kBlue);
	effMu_pt_arbitr->GetYaxis()->SetRangeUser(0.,1.);
	effMu_pt_arbitr->Draw("LP");
	cv10->Modified();
	TCanvas *cv11=new TCanvas("cv11","cv11",100,100,900,500);
	SetMyStyle();
	gStyle->SetOptStat(0);	
	effMu_pt_oneSta->SetMarkerStyle(22);
	effMu_pt_oneSta->SetMarkerColor(kBlue);
	effMu_pt_oneSta->SetLineColor(kBlue);
	effMu_pt_oneSta->GetYaxis()->SetRangeUser(0.,1.);
	effMu_pt_oneSta->Draw("LP");
	cv11->Modified();
	TCanvas *cv12=new TCanvas("cv12","cv12",100,100,900,500);
	SetMyStyle();
	gStyle->SetOptStat(0);	
	effMu_pt_lastSta->SetMarkerStyle(22);
	effMu_pt_lastSta->SetMarkerColor(kBlue);
	effMu_pt_lastSta->SetLineColor(kBlue);
	effMu_pt_lastSta->GetYaxis()->SetRangeUser(0.,1.);
	effMu_pt_lastSta->Draw("LP");
	cv12->Modified();
	cv10->SaveAs("effMu_pt_arbitr.eps");
	cv11->SaveAs("effMu_pt_oneSta.eps");
	cv12->SaveAs("effMu_pt_lastSta.eps");
	cv10->SaveAs("effMu_pt_arbitr.gif");
	cv11->SaveAs("effMu_pt_oneSta.gif");
	cv12->SaveAs("effMu_pt_lastSta.gif");
	TCanvas *cv20=new TCanvas("cv20","cv20",100,100,900,500);
	SetMyStyle();
	gStyle->SetOptStat(0);	
	effMu_eta_arbitr->SetMarkerStyle(22);
	effMu_eta_arbitr->SetMarkerColor(kBlue);
	effMu_eta_arbitr->SetLineColor(kBlue);
	effMu_eta_arbitr->GetYaxis()->SetRangeUser(0.,1.);
	effMu_eta_arbitr->Draw("LP");
	cv20->Modified();
	TCanvas *cv21=new TCanvas("cv21","cv21",100,100,900,500);
	SetMyStyle();
	gStyle->SetOptStat(0);	
	effMu_eta_oneSta->SetMarkerStyle(22);
	effMu_eta_oneSta->SetMarkerColor(kBlue);
	effMu_eta_oneSta->SetLineColor(kBlue);
	effMu_eta_oneSta->GetYaxis()->SetRangeUser(0.,1.);
	effMu_eta_oneSta->Draw("LP");
	cv21->Modified();
	TCanvas *cv22=new TCanvas("cv22","cv22",100,100,900,500);
	SetMyStyle();
	gStyle->SetOptStat(0);	
	effMu_eta_lastSta->SetMarkerStyle(22);
	effMu_eta_lastSta->SetMarkerColor(kBlue);
	effMu_eta_lastSta->SetLineColor(kBlue);
	effMu_eta_lastSta->GetYaxis()->SetRangeUser(0.,1.);
	effMu_eta_lastSta->Draw("LP");
	cv22->Modified();
	cv20->SaveAs("effMu_eta_arbitr.eps");
	cv21->SaveAs("effMu_eta_oneSta.eps");
	cv22->SaveAs("effMu_eta_lastSta.eps");
	cv20->SaveAs("effMu_eta_arbitr.gif");
	cv21->SaveAs("effMu_eta_oneSta.gif");
	cv22->SaveAs("effMu_eta_lastSta.gif");
}


