#ifndef __CINT__
#include <RooGlobalFunc.h>
#endif

#include "TLatex.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooLandau.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "TAxis.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMath.h"
#include "TVector3.h"
#include "vector.h"
#include "TString.h"
#include "iostream.h"
#include "fstream.h"
#include "TLegend.h"
#include "TStyle.h"

#include "fitPeak.cc"
#include "myStyle.h"

using namespace RooFit;
void cutSBG()
{
	gROOT->Reset();
	const int nCut=9;
	const int nBin=3;
	double SBG[nBin][nCut];
	double nSig[nBin][nCut];
	double cutVar[nCut]={0.,1.,2.,3.,4.,5.,6.,7.,8.};
	char *label[nCut]={"noCut","innerTrackerHits","pixelLayers","innerTrackerChi2","dxy","dz","globalChi2","validMuHits","vertexProb"};
	ofstream printTxt("SBG.txt");

	Long_t nEntries, nData=0;
	char fileName[500];
	sprintf(fileName,
			"/home/zhlinl/cutAnalysis/test/macro/reducedTree.root"
			);
	Double_t massMin,massMax;
	massMin=2.6;
	massMax=3.5;
	RooRealVar mass("mass","dimuon mass(GeV/c^{2})",massMin,massMax);
	RooDataSet *dataSet[nBin][nCut];

	for(int i=0; i<nBin; i++)
	{
		for(int j=0; j<nCut; j++)
		{
			dataSet[i][j]=new RooDataSet(Form("dataSet%d%d",i,j),"signal+background",mass);
		}
	}
	//Jpsi Variables
	Double_t JpsiMass,JpsiPt,JpsiRap;
	Double_t JpsiVprob;
	//(1).Positive Muon                                     
	double muPos_nchi2In, muPos_dxy, muPos_dz, muPos_nchi2Gl;
	int muPos_arbitrated, muPos_oneStationTight;
	int muPos_found, muPos_pixeLayers, muPos_nValidMuHits;
	//(2).Negative Muon                                     
	double muNeg_nchi2In, muNeg_dxy, muNeg_dz, muNeg_nchi2Gl;
	int muNeg_arbitrated, muNeg_oneStationTight;
	int muNeg_found, muNeg_pixeLayers, muNeg_nValidMuHits;
	//Trigger Information                                   
	int HLT_L1DoubleMuOpen, HLT_Mu0_TkMu0_Jpsi, HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2;

	TFile *infile=new TFile(fileName,"R");
	TTree *tree=(TTree*)infile->Get("data");
	//Jpsi Variables
	tree->SetBranchAddress("JpsiMass",&JpsiMass);
	tree->SetBranchAddress("JpsiPt",&JpsiPt);
	tree->SetBranchAddress("JpsiRap",&JpsiRap);
	tree->SetBranchAddress("JpsiVprob",&JpsiVprob);
	//Trigger information
	tree->SetBranchAddress("HLT_L1DoubleMuOpen",&HLT_L1DoubleMuOpen);
	tree->SetBranchAddress("HLT_Mu0_TkMu0_Jpsi",&HLT_Mu0_TkMu0_Jpsi);
	tree->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2",&HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2);
	//1) Positive Muon                                      
	tree->SetBranchAddress("muPos_nchi2In", &muPos_nchi2In);
	tree->SetBranchAddress("muPos_dxy", &muPos_dxy);
	tree->SetBranchAddress("muPos_dz", &muPos_dz);
	tree->SetBranchAddress("muPos_nchi2Gl", &muPos_nchi2Gl);
	tree->SetBranchAddress("muPos_arbitrated", &muPos_arbitrated);
	tree->SetBranchAddress("muPos_oneStationTight", &muPos_oneStationTight);
	tree->SetBranchAddress("muPos_found", &muPos_found);
	tree->SetBranchAddress("muPos_pixeLayers", &muPos_pixeLayers);
	tree->SetBranchAddress("muPos_nValidMuHits", &muPos_nValidMuHits);
	//2) Negative Muon                                      
	tree->SetBranchAddress("muNeg_nchi2In", &muNeg_nchi2In);
	tree->SetBranchAddress("muNeg_dxy", &muNeg_dxy);
	tree->SetBranchAddress("muNeg_dz", &muNeg_dz);
	tree->SetBranchAddress("muNeg_nchi2Gl", &muNeg_nchi2Gl);
	tree->SetBranchAddress("muNeg_arbitrated", &muNeg_arbitrated);
	tree->SetBranchAddress("muNeg_oneStationTight", &muNeg_oneStationTight);
	tree->SetBranchAddress("muNeg_found", &muNeg_found);
	tree->SetBranchAddress("muNeg_pixeLayers", &muNeg_pixeLayers);
	tree->SetBranchAddress("muNeg_nValidMuHits", &muNeg_nValidMuHits);

	nEntries=tree->GetEntries();
	cout<<"===============total Entries in Tree: "<<nEntries<<"==============="<<endl;

	Int_t nPeak=1;
	vector<Double_t> fMean,dfMean,fSigma,dfSigma;
	fMean.push_back(3.09);
	dfMean.push_back(0.01);
	fSigma.push_back(0.04);
	dfSigma.push_back(0.04);

	vector<double> sigBkg[nBin][nCut];

	bool cutPass[nCut];
	for(int i=0; i<nCut; i++)
	{

		for(int k=0; k<nEntries; k++)
		{
			tree->GetEntry(k);

			if(i==0) cutPass[i]= 1;
			if(i==1) cutPass[i]= muPos_found > 11 && muNeg_found > 11;
			if(i==2) cutPass[i]=  muPos_found > 11 && muNeg_found > 11 
				&& muPos_pixeLayers > 1 && muNeg_pixeLayers > 1;
			if(i==3) cutPass[i]= muPos_found > 11 && muNeg_found > 11 
				&& muPos_pixeLayers > 1 && muNeg_pixeLayers > 1 
					&& muPos_nchi2In < 1.8 && muNeg_nchi2In < 1.8;
			if(i==4) cutPass[i]=muPos_found > 11 && muNeg_found > 11 
				&& muPos_pixeLayers > 1 && muNeg_pixeLayers > 1 
					&& muPos_nchi2In < 1.8 && muNeg_nchi2In < 1.8
					&& fabs(muPos_dxy) < 0.1 && fabs(muNeg_dxy) < 0.1;
			if(i==5) cutPass[i]= muPos_found > 11 && muNeg_found > 11 
				&& muPos_pixeLayers > 1 && muNeg_pixeLayers > 1 
					&& muPos_nchi2In < 1.8 && muNeg_nchi2In < 1.8
					&& fabs(muPos_dxy) < 0.1 && fabs(muNeg_dxy) < 0.1
					&& fabs(muPos_dz) < 15. && fabs(muNeg_dz) < 15.;
			if(i==6) cutPass[i]= muPos_found > 11 && muNeg_found > 11 
				&& muPos_pixeLayers > 1 && muNeg_pixeLayers > 1 
					&& muPos_nchi2In < 1.8 && muNeg_nchi2In < 1.8
					&& fabs(muPos_dxy) < 0.1 && fabs(muNeg_dxy) < 0.1
					&& fabs(muPos_dz) < 15. && fabs(muNeg_dz) < 15.
					&& (muPos_nchi2Gl < -100. || (muPos_nchi2Gl > -100. && muPos_nchi2Gl < 20. ))
					&& (muNeg_nchi2Gl < -100. || (muNeg_nchi2Gl > -100. && muNeg_nchi2Gl < 20. ));
			if(i==7) cutPass[i]=muPos_found > 11 && muNeg_found > 11 
				&& muPos_pixeLayers > 1 && muNeg_pixeLayers > 1 
					&& muPos_nchi2In < 1.8 && muNeg_nchi2In < 1.8
					&& fabs(muPos_dxy) < 0.1 && fabs(muNeg_dxy) < 0.1
					&& fabs(muPos_dz) < 15. && fabs(muNeg_dz) < 15.
					&& (muPos_nchi2Gl < -100. || (muPos_nchi2Gl > -100. && muPos_nchi2Gl < 20. ))
					&& (muNeg_nchi2Gl < -100. || (muNeg_nchi2Gl > -100. && muNeg_nchi2Gl < 20. ))
					&& (muPos_nValidMuHits < -100 || (muPos_nValidMuHits > -100 && muPos_nValidMuHits > 0))
					&& (muNeg_nValidMuHits < -100 || (muNeg_nValidMuHits > -100 && muNeg_nValidMuHits > 0)); 
			if(i==8) cutPass[i]= muPos_found > 11 && muNeg_found > 11 
				&& muPos_pixeLayers > 1 && muNeg_pixeLayers > 1 
					&& muPos_nchi2In < 1.8 && muNeg_nchi2In < 1.8
					&& fabs(muPos_dxy) < 0.1 && fabs(muNeg_dxy) < 0.1
					&& fabs(muPos_dz) < 15. && fabs(muNeg_dz) < 15.
					&& (muPos_nchi2Gl < -100. || (muPos_nchi2Gl > -100. && muPos_nchi2Gl < 20. ))
					&& (muNeg_nchi2Gl < -100. || (muNeg_nchi2Gl > -100. && muNeg_nchi2Gl < 20. ))
					&& (muPos_nValidMuHits < -100 || (muPos_nValidMuHits > -100 && muPos_nValidMuHits > 0))
					&& (muNeg_nValidMuHits < -100 || (muNeg_nValidMuHits > -100 && muNeg_nValidMuHits > 0))
					&& JpsiVprob > 0.01;


			if( HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2 == 1 
					&& muPos_arbitrated ==1 && muNeg_arbitrated ==1 
					&& muPos_oneStationTight ==1 && muNeg_oneStationTight ==1  
				)
			{	
				mass.setVal(JpsiMass);
				if(fabs(JpsiRap)<=0.9 && cutPass[i])
				{
					dataSet[0][i]->add(mass);
					nData++;
				}  
				if(fabs(JpsiRap)>0.9 && fabs(JpsiRap)<1.5 && cutPass[i])
				{
					dataSet[1][i]->add(mass);
					nData++;
				}
				if(fabs(JpsiRap)>1.5 && fabs(JpsiRap)<2.1 && cutPass[i])
				{
					dataSet[2][i]->add(mass);
					nData++;
				}
			}
		}
		cout<<"======total matched Entries in Pt<6GeV && 0<|y|<2.1 : "<<nData<<" ======"<<endl;
		nData=0;

		TCanvas *c1[nBin];
		for(int j=0; j<nBin; j++)
		{
			c1[j]=new TCanvas(Form("c1_%d",j),"");
			c1[j]->SetFillColor(10);
			sigBkg[j][i].clear();
			sigBkg[j][i]=fitPeak(mass,0,nPeak,0,dataSet[j][i],fMean,dfMean,fSigma,dfSigma,0.2);
			c1[j]->Modified();
			c1[j]->Update();
			SBG[j][i]=sigBkg[j][i][2];
			nSig[j][i]=sigBkg[j][i][6];
		}
	}

	for(int i=0; i<nBin; i++)
	{
		for(int j=0; j<nCut; j++)
		{
			cout<<"SBG["<<i<<"]["<<j<<"]:"<<SBG[i][j]<<" "<<endl;
			printTxt<<"SBG["<<i<<"]["<<j<<"]:"<<SBG[i][j]<<" "<<endl;
		}
		cout<<endl;
		printTxt<<endl;
		for(int j=0; j<nCut; j++)
		{
			cout<<"nSig["<<i<<"]["<<j<<"]:"<<nSig[i][j]<<" "<<endl;
			printTxt<<"nSig["<<i<<"]["<<j<<"]:"<<nSig[i][j]<<" "<<endl;
		}
		cout<<endl;
		printTxt<<endl;
	}

	printTxt.close();

	TCanvas *c5[nBin];
	for(int i=0; i<nBin; i++)
	{
		c5[i]=new TCanvas(Form("c5_%d",i),"");
		c5[i]->SetFillColor(10);

		SetMyStyle();
		gStyle->SetOptStat(0);
		gStyle->SetTitleXOffset(1.2); 
		gStyle->SetTitleYOffset(1.3); 
		gPad->SetLeftMargin(0.12);  
		gPad->SetBottomMargin(0.12);
		gPad->SetRightMargin(0.12); 
		gPad->SetTopMargin(0.12);

		TH1F *histSB=new TH1F("histSB","histSB",nCut,0,nCut);
		histSB->SetMarkerColor(2);
		histSB->SetMarkerStyle(23);
		histSB->SetMarkerSize(1.5);
		histSB->SetLineColor(2);
		histSB->GetXaxis()->CenterTitle();
		histSB->GetYaxis()->SetTitle("S/BG");
		histSB->GetYaxis()->CenterTitle();
		if(i==0) histSB->SetTitle(Form("|y|<%.1f, P_{T}>%.1f",0.9,6.0));
		if(i==1) histSB->SetTitle(Form("%.1f<|y|<%.1f, P_{T}>%.1f",0.9,1.5,6.0));
		if(i==2) histSB->SetTitle(Form("%.1f<|y|<%.1f, P_{T}>%.1f",1.5,2.1,6.0));

		for(int k=0; k<nCut; k++)
		{
			histSB->SetBinContent(k+1,SBG[i][k]);
			histSB->GetXaxis()->SetBinLabel(k+1,label[k]);
		}
		histSB->Draw("LP");
		c5[i]->SaveAs(Form("./pic/SBG/SBG_EtaBin%d.eps",i));
		c5[i]->SaveAs(Form("./pic/SBG/SBG_EtaBin%d.gif",i));
		c5[i]->SaveAs(Form("./pic/SBG/SBG_EtaBin%d.C",i));
	}
	return;

}
