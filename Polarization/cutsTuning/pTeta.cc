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
//#include <vector.h>       
#include <TString.h>      
//#include <iostream.h>     
#include <TLegend.h> 
#include <TPad.h> 
#include "/home/zhlinl/work/codes/myStyle.h"
void pTeta()//int which=0)
{

	gROOT->Reset();
	//gSystem->Load("$ROOTSYS/lib/libPhysics.so");
	//gStyle->SetOptStat(1);
	Long_t nEntries;  
	int passedEvt=0;

	int which=0;  //0,1,2,3
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

	TH2F *histMu_pt_eta=new TH2F("histMu_pt_eta","Mu_pt_eta",100,-2.6,2.6,100,0,25);
	histMu_pt_eta->SetXTitle("#eta");
	histMu_pt_eta->SetYTitle("P_{T} (GeV/c)");
	TH1F *histMu_pt=new TH1F("histMu_pt","Mu_pt",100,0,25);
	histMu_pt->SetXTitle("P_{T} (GeV/c)");
	histMu_pt->SetYTitle("Events");
	TH1F *histMu_eta=new TH1F("histMu_eta","Mu_eta",100,-2.6,2.6);
	histMu_eta->SetXTitle("#eta");//Rapidity");
	histMu_eta->SetYTitle("Events");

	TH1F *accMu_pt=new TH1F("accMu_pt","accMu_pt",100,0,25);
	TH1F *accMu_eta=new TH1F("accMu_eta","accMu_eta",100,-2.6,2.6);
	

	TH2F *histMuPos_pt_eta=new TH2F("histMuPos_pt_eta","MuPos_pt_eta",100,-2.6,2.6,100,0,25);
	histMuPos_pt_eta->SetXTitle("#eta");
	histMuPos_pt_eta->SetYTitle("P_{T} (GeV/c)");
	TH1F *histMuPos_pt=new TH1F("histMuPos_pt","MuPos_pt",100,0,25);
	histMuPos_pt->SetXTitle("P_{T} (GeV/c)");
	histMuPos_pt->SetYTitle("Events");
	TH1F *histMuPos_eta=new TH1F("histMuPos_eta","MuPos_eta",100,-2.6,2.6);
	histMuPos_eta->SetXTitle("#eta");//Rapidity");
	histMuPos_eta->SetYTitle("Events");
	
	TH2F *histMuNeg_pt_eta=new TH2F("histMuNeg_pt_eta","MuNeg_pt_eta",100,-2.6,2.6,100,0,25);
	histMuNeg_pt_eta->SetXTitle("#eta");
	histMuNeg_pt_eta->SetYTitle("P_{T} (GeV/c)");
	TH1F *histMuNeg_pt=new TH1F("histMuNeg_pt","MuNeg_pt",100,0,25);
	histMuNeg_pt->SetXTitle("P_{T} (GeV/c)");
	histMuNeg_pt->SetYTitle("Events");
	TH1F *histMuNeg_eta=new TH1F("histMuNeg_eta","MuNeg_eta",100,-2.6,2.6);
	histMuNeg_eta->SetXTitle("#eta");//Rapidity");
	histMuNeg_eta->SetYTitle("Events");

	TH2F *histJpsi_pt_eta=new TH2F("histJpsi_pt_eta","Jpsi_pt_eta",100,-2.9,2.9,100,0,25);
	histJpsi_pt_eta->SetXTitle("#eta");
	histJpsi_pt_eta->SetYTitle("P_{T} (GeV/c)");
	TH1F *histJpsi_pt=new TH1F("histJpsi_pt","Jpsi_pt",100,0,25);
	histJpsi_pt->SetXTitle("P_{T} (GeV/c)");
	histJpsi_pt->SetYTitle("Events");
	TH1F *histJpsi_eta=new TH1F("histJpsi_eta","Jpsi_eta",100,-2.9,2.9);
	histJpsi_eta->SetXTitle("#eta");//Rapidity");
	histJpsi_eta->SetYTitle("Events");
	

	nEntries=tree->GetEntries();

	cout<<"===============total Entries in Tree: "<<nEntries<<"==============="<<endl;
	//for(int which=0; which<4; which++)
	//{
		bool cut=kFALSE;
		for(int i=0; i<nEntries/1; i++)
		{
			tree->GetEntry(i);
			//bool cut;
			if(which==0) cut=kTRUE ;
			if(which==1) cut= muPos_arbitrated ==1 && muNeg_arbitrated ==1 ;
			if(which==2) cut= muPos_arbitrated ==1 && muNeg_arbitrated ==1
				&& muPos_oneStationTight ==1 && muNeg_oneStationTight ==1 ;
			if(which==3) cut= muPos_arbitrated ==1 && muNeg_arbitrated ==1
				&& muPos_lastStationAngTight ==1 && muNeg_lastStationAngTight ==1 ;
			if((i<=100 && i%20==0)
					||(i>100 && i<=1000 && i%100==0)
					||(i>1000 && i%1000==0) )
			{   
				cout<<"entryID "<<i<<" in Total Entries "<<nEntries<<endl;
			}

			JpsiMass=JpsiP->M();
			JpsiPt=JpsiP->Pt();
			//JpsiRap=JpsiP->Rapidity();
			JpsiRap=JpsiP->Eta();
			//histMu_pt_eta->Fill(muPos.Rap(),muPos.Pt());
			if( 1
					//&& JpsiMass>massMin && JpsiMass<massMax
					&& HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2 ==1
					//&& muPos_arbitrated ==1 && muNeg_arbitrated ==1
					//&& cut 
				)
			{
				//histMu_pt_eta->Fill(muPosP->Rapidity(),muPosP->Pt());
				//histMu_pt_eta->Fill(muNegP->Rapidity(),muNegP->Pt());
				//histMu_eta->Fill(muPosP->Rapidity());
				//histMu_eta->Fill(muNegP->Rapidity());
				histMu_pt_eta->Fill(muPosP->Eta(),muPosP->Pt());
				histMu_pt_eta->Fill(muNegP->Eta(),muNegP->Pt());
				histMu_pt->Fill(muPosP->Pt());
				histMu_pt->Fill(muNegP->Pt());
				histMu_eta->Fill(muPosP->Eta());
				histMu_eta->Fill(muNegP->Eta());
				
				histMuPos_pt_eta->Fill(muPosP->Eta(),muPosP->Pt());
				histMuNeg_pt_eta->Fill(muNegP->Eta(),muNegP->Pt());
				histMuPos_pt->Fill(muPosP->Pt());
				histMuNeg_pt->Fill(muNegP->Pt());
				histMuPos_eta->Fill(muPosP->Eta());
				histMuNeg_eta->Fill(muNegP->Eta());
				histJpsi_pt_eta->Fill(JpsiRap,JpsiPt);
				histJpsi_pt->Fill(JpsiPt);
				histJpsi_eta->Fill(JpsiRap);
				
				passedEvt++;
			}
		}
	
		cout<<"========totally passed Events: "<<passedEvt<<"============"<<endl;
		TCanvas *cv10=new TCanvas("cv10","cv10",100,100,900,500);
		//SetMyStyle();
		gStyle->SetPalette(1);
		gStyle->SetOptStat(0);
		histMu_pt_eta->Draw("colz");
		cv10->Modified();
		TCanvas *cv20=new TCanvas("cv20","cv20",100,100,900,500);
		cv20->SetLogy(1);
		//SetMyStyle();
		gStyle->SetOptStat(0);	
		//gStyle->SetTitleYOffset(1.2);
		histMu_pt->GetYaxis()->SetTitleOffset(1.2);
		histMu_pt->Draw();
		cv20->Modified();
		TCanvas *cv30=new TCanvas("cv30","cv30",100,100,900,500);
		//SetMyStyle();
		gStyle->SetOptStat(0);	
		//gStyle->SetTitleYOffset(1.2);
		histMu_eta->GetYaxis()->SetTitleOffset(1.2);
		histMu_eta->Draw();
		cv30->Modified();


		TCanvas *cv11=new TCanvas("cv11","cv11",100,100,900,500);
		//SetMyStyle();
		gStyle->SetPalette(1);
		gStyle->SetOptStat(0);
		histMuPos_pt_eta->Draw("colz");
		cv11->Modified();
		TCanvas *cv21=new TCanvas("cv21","cv21",100,100,900,500);
		//SetMyStyle();
		gStyle->SetOptStat(0);	
		//gStyle->SetTitleYOffset(1.2);
		histMuPos_pt->GetYaxis()->SetTitleOffset(1.2);
		histMuPos_pt->Draw();
		cv21->Modified();
		TCanvas *cv31=new TCanvas("cv31","cv31",100,100,900,500);
		//SetMyStyle();
		gStyle->SetOptStat(0);	
		//gStyle->SetTitleYOffset(1.2);
		histMuPos_eta->GetYaxis()->SetTitleOffset(1.2);
		histMuPos_eta->Draw();
		cv31->Modified();

		TCanvas *cv12=new TCanvas("cv12","cv12",100,100,900,500);
		//SetMyStyle();
		gStyle->SetPalette(1);
		gStyle->SetOptStat(0);
		histMuNeg_pt_eta->Draw("colz");
		cv12->Modified();
		TCanvas *cv22=new TCanvas("cv22","cv22",100,100,900,500);
		//SetMyStyle();
		gStyle->SetOptStat(0);	
		//gStyle->SetTitleYOffset(1.2);
		histMuNeg_pt->GetYaxis()->SetTitleOffset(1.2);
		histMuNeg_pt->Draw();
		cv22->Modified();
		TCanvas *cv32=new TCanvas("cv32","cv32",100,100,900,500);
		//SetMyStyle();
		gStyle->SetOptStat(0);	
		//gStyle->SetTitleYOffset(1.2);
		histMuNeg_eta->GetYaxis()->SetTitleOffset(1.2);
		histMuNeg_eta->Draw();
		cv32->Modified();

		TCanvas *cv4=new TCanvas("cv4","cv4",100,100,900,500);
		//SetMyStyle();
		gStyle->SetOptStat(0); 
		histJpsi_pt_eta->Draw("colz");
		cv4->Modified();
		TCanvas *cv5=new TCanvas("cv5","cv5",100,100,900,500);
		//SetMyStyle();
		gStyle->SetOptStat(0); 
		histJpsi_pt->Draw();
		cv5->Modified();
		TCanvas *cv6=new TCanvas("cv6","cv6",100,100,900,500);
		//SetMyStyle();
		gStyle->SetOptStat(0); 
		//gStyle->SetStatX(0.6);
		//gStyle->SetStatY(0.9);
		histJpsi_eta->Draw();
		cv6->Modified();
		if(which==0)
		{
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta.eps");
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta.gif");
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta.C");
			cv20->SaveAs("./pic/arbitrated/mu_pt.eps");
			cv20->SaveAs("./pic/arbitrated/mu_pt.gif");
			cv20->SaveAs("./pic/arbitrated/mu_pt.C");
			cv30->SaveAs("./pic/arbitrated/mu_eta.eps");
			cv30->SaveAs("./pic/arbitrated/mu_eta.gif");
			cv30->SaveAs("./pic/arbitrated/mu_eta.C");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta.eps");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta.gif");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta.C");
			cv21->SaveAs("./pic/arbitrated/muPos_pt.eps");
			cv21->SaveAs("./pic/arbitrated/muPos_pt.gif");
			cv21->SaveAs("./pic/arbitrated/muPos_pt.C");
			cv31->SaveAs("./pic/arbitrated/muPos_eta.eps");
			cv31->SaveAs("./pic/arbitrated/muPos_eta.gif");
			cv31->SaveAs("./pic/arbitrated/muPos_eta.C");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta.eps");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta.gif");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta.C");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt.eps");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt.gif");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt.C");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta.eps");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta.gif");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta.C");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta.eps");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta.gif");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta.C");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt.eps");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt.gif");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt.C");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta.eps");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta.gif");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta.C");
		}
		if(which==1)
		{
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta_arbitr.eps");
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta_arbitr.gif");
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta_arbitr.C");
			cv20->SaveAs("./pic/arbitrated/mu_pt_arbitr.eps");
			cv20->SaveAs("./pic/arbitrated/mu_pt_arbitr.gif");
			cv20->SaveAs("./pic/arbitrated/mu_pt_arbitr.C");
			cv30->SaveAs("./pic/arbitrated/mu_eta_arbitr.eps");
			cv30->SaveAs("./pic/arbitrated/mu_eta_arbitr.gif");
			cv30->SaveAs("./pic/arbitrated/mu_eta_arbitr.C");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta_arbitr.eps");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta_arbitr.gif");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta_arbitr.C");
			cv21->SaveAs("./pic/arbitrated/muPos_pt_arbitr.eps");
			cv21->SaveAs("./pic/arbitrated/muPos_pt_arbitr.gif");
			cv21->SaveAs("./pic/arbitrated/muPos_pt_arbitr.C");
			cv31->SaveAs("./pic/arbitrated/muPos_eta_arbitr.eps");
			cv31->SaveAs("./pic/arbitrated/muPos_eta_arbitr.gif");
			cv31->SaveAs("./pic/arbitrated/muPos_eta_arbitr.C");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta_arbitr.eps");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta_arbitr.gif");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta_arbitr.C");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt_arbitr.eps");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt_arbitr.gif");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt_arbitr.C");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta_arbitr.eps");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta_arbitr.gif");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta_arbitr.C");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta_arbitr.eps");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta_arbitr.gif");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta_arbitr.C");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt_arbitr.eps");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt_arbitr.gif");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt_arbitr.C");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta_arbitr.eps");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta_arbitr.gif");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta_arbitr.C");
		}
		if(which==2)
		{
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta_oneSta.eps");
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta_oneSta.gif");
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta_oneSta.C");
			cv20->SaveAs("./pic/arbitrated/mu_pt_oneSta.eps");
			cv20->SaveAs("./pic/arbitrated/mu_pt_oneSta.gif");
			cv20->SaveAs("./pic/arbitrated/mu_pt_oneSta.C");
			cv30->SaveAs("./pic/arbitrated/mu_eta_oneSta.eps");
			cv30->SaveAs("./pic/arbitrated/mu_eta_oneSta.gif");
			cv30->SaveAs("./pic/arbitrated/mu_eta_oneSta.C");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta_oneSta.eps");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta_oneSta.gif");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta_oneSta.C");
			cv21->SaveAs("./pic/arbitrated/muPos_pt_oneSta.eps");
			cv21->SaveAs("./pic/arbitrated/muPos_pt_oneSta.gif");
			cv21->SaveAs("./pic/arbitrated/muPos_pt_oneSta.C");
			cv31->SaveAs("./pic/arbitrated/muPos_eta_oneSta.eps");
			cv31->SaveAs("./pic/arbitrated/muPos_eta_oneSta.gif");
			cv31->SaveAs("./pic/arbitrated/muPos_eta_oneSta.C");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta_oneSta.eps");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta_oneSta.gif");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta_oneSta.C");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt_oneSta.eps");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt_oneSta.gif");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt_oneSta.C");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta_oneSta.eps");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta_oneSta.gif");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta_oneSta.C");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta_oneSta.eps");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta_oneSta.gif");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta_oneSta.C");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt_oneSta.eps");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt_oneSta.gif");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt_oneSta.C");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta_oneSta.eps");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta_oneSta.gif");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta_oneSta.C");
		}
		if(which==3)
		{
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta_lastSta.eps");
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta_lastSta.gif");
			cv10->SaveAs("./pic/arbitrated/mu_pt_eta_lastSta.C");
			cv20->SaveAs("./pic/arbitrated/mu_pt_lastSta.eps");
			cv20->SaveAs("./pic/arbitrated/mu_pt_lastSta.gif");
			cv20->SaveAs("./pic/arbitrated/mu_pt_lastSta.C");
			cv30->SaveAs("./pic/arbitrated/mu_eta_lastSta.eps");
			cv30->SaveAs("./pic/arbitrated/mu_eta_lastSta.gif");
			cv30->SaveAs("./pic/arbitrated/mu_eta_lastSta.C");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta_lastSta.eps");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta_lastSta.gif");
			cv11->SaveAs("./pic/arbitrated/muPos_pt_eta_lastSta.C");
			cv21->SaveAs("./pic/arbitrated/muPos_pt_lastSta.eps");
			cv21->SaveAs("./pic/arbitrated/muPos_pt_lastSta.gif");
			cv21->SaveAs("./pic/arbitrated/muPos_pt_lastSta.C");
			cv31->SaveAs("./pic/arbitrated/muPos_eta_lastSta.eps");
			cv31->SaveAs("./pic/arbitrated/muPos_eta_lastSta.gif");
			cv31->SaveAs("./pic/arbitrated/muPos_eta_lastSta.C");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta_lastSta.eps");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta_lastSta.gif");
			cv12->SaveAs("./pic/arbitrated/muNeg_pt_eta_lastSta.C");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt_lastSta.eps");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt_lastSta.gif");
			cv22->SaveAs("./pic/arbitrated/muNeg_pt_lastSta.C");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta_lastSta.eps");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta_lastSta.gif");
			cv32->SaveAs("./pic/arbitrated/muNeg_eta_lastSta.C");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta_lastSta.eps");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta_lastSta.gif");
			cv4->SaveAs("./pic/arbitrated/jpsi_pt_eta_lastSta.C");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt_lastSta.eps");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt_lastSta.gif");
			cv5->SaveAs("./pic/arbitrated/jpsi_pt_lastSta.C");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta_lastSta.eps");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta_lastSta.gif");
			cv6->SaveAs("./pic/arbitrated/jpsi_eta_lastSta.C");
		}
	//}
}
