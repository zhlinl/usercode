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

void reduceTree()
{
	gROOT->Reset();
	Long_t nEntries;
	
	char fileName[500];
	
	Double_t massMin,massMax;
	massMin=2.6;
	massMax=3.5;
	
	sprintf(fileName,
	"/scratch/zhlinl/muonTree_Run2010B-Nov4ReReco_v1-Onia2MuMu-v_new/muonTree_nocut_1.root"
	);
	
	//Jpsi Variables
	Double_t JpsiMass,JpsiPt,JpsiRap;
	Double_t JpsiVprob;
	TLorentzVector *muPosP,*muNegP,*JpsiP;
	//(1).Positive Muon                                     
	double muPos_nchi2In, muPos_dxy, muPos_dz, muPos_nchi2Gl;
	int muPos_arbitrated, muPos_oneStationTight, muPos_lastStationAngTight,
			muPos_lastStationTight, muPos_lastStationOptimizedLowPtTight,
			muPos_lastStationOptimizedBarrelLowPtTight,muPos_oneStationAngTight;
	int muPos_found, muPos_pixeLayers, muPos_nValidMuHits;
	//(2).Negative Muon                                     
	double muNeg_nchi2In, muNeg_dxy, muNeg_dz, muNeg_nchi2Gl;
	int muNeg_arbitrated, muNeg_oneStationTight, muNeg_lastStationAngTight,
			muNeg_lastStationTight, muNeg_lastStationOptimizedLowPtTight,
			muNeg_lastStationOptimizedBarrelLowPtTight, muNeg_oneStationAngTight;
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
	tree->SetBranchAddress("JpsiP",&JpsiP);
	tree->SetBranchAddress("muPosP",&muPosP);
	tree->SetBranchAddress("muNegP",&muNegP);
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
	tree->SetBranchAddress("muPos_lastStationAngTight",&muPos_lastStationAngTight);
	tree->SetBranchAddress("muPos_lastStationTight",&muPos_lastStationTight);
	tree->SetBranchAddress("muPos_lastStationOptimizedLowPtTight",&muPos_lastStationOptimizedLowPtTight);
	tree->SetBranchAddress("muPos_lastStationOptimizedBarrelLowPtTight",&muPos_lastStationOptimizedBarrelLowPtTight);
	tree->SetBranchAddress("muPos_oneStationAngTight",&muPos_oneStationAngTight);
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
	tree->SetBranchAddress("muNeg_lastStationAngTight",&muNeg_lastStationAngTight);
	tree->SetBranchAddress("muNeg_lastStationTight",&muNeg_lastStationTight);
	tree->SetBranchAddress("muNeg_lastStationOptimizedLowPtTight",&muNeg_lastStationOptimizedLowPtTight);
	tree->SetBranchAddress("muNeg_lastStationOptimizedBarrelLowPtTight",&muNeg_lastStationOptimizedBarrelLowPtTight);
	tree->SetBranchAddress("muNeg_oneStationAngTight",&muNeg_oneStationAngTight);

	TFile *outfile=new TFile("reducedTree.root","recreate");
	TTree *data=new TTree("data","reduced data Tree");
	
	data->Branch("HLT_L1DoubleMuOpen",&HLT_L1DoubleMuOpen,"HLT_L1DoubleMuOpen/I");
	data->Branch("HLT_Mu0_TkMu0_Jpsi",&HLT_Mu0_TkMu0_Jpsi,"HLT_Mu0_TkMu0_Jpsi/I");
	data->Branch("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2",&HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2,"HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2/I");
	//Jpsi Variables	
	data->Branch("JpsiMass",   &JpsiMass,  "JpsiMass/D");
	data->Branch("JpsiPt",     &JpsiPt,    "JpsiPt/D");
	data->Branch("JpsiRap",    &JpsiRap,   "JpsiRap/D");
	data->Branch("JpsiVprob",  &JpsiVprob, "JpsiVprob/D");
	data->Branch("JpsiP","TLorentzVector",&JpsiP);
	data->Branch("muPosP","TLorentzVector",&muPosP);
	data->Branch("muNegP","TLorentzVector",&muNegP);
	//1). Positive Muon
	data->Branch("muPos_nchi2In", &muPos_nchi2In, "muPos_nchi2In/D");
	data->Branch("muPos_dxy", &muPos_dxy, "muPos_dxy/D");
	data->Branch("muPos_dz", &muPos_dz, "muPos_dz/D");
	data->Branch("muPos_nchi2Gl", &muPos_nchi2Gl, "muPos_nchi2Gl/D");
	data->Branch("muPos_arbitrated", &muPos_arbitrated, "muPos_arbitrated/I");
	data->Branch("muPos_oneStationTight", &muPos_oneStationTight, "muPos_oneStationTight/I");
	data->Branch("muPos_found", &muPos_found, "muPos_found/I");
	data->Branch("muPos_pixeLayers", &muPos_pixeLayers, "muPos_pixeLayers/I");
	data->Branch("muPos_nValidMuHits", &muPos_nValidMuHits, "muPos_nValidMuHits/I");
	data->Branch("muPos_lastStationAngTight",&muPos_lastStationAngTight, "muPos_lastStationAngTight/I");
	data->Branch("muPos_lastStationTight",&muPos_lastStationTight, "muPos_lastStationTight/I");
	data->Branch("muPos_lastStationOptimizedLowPtTight",&muPos_lastStationOptimizedLowPtTight, "muPos_lastStatio nOptimizedLowPtTight/I");
	data->Branch("muPos_lastStationOptimizedBarrelLowPtTight",&muPos_lastStationOptimizedBarrelLowPtTight, "muPo s_lastStationOptimizedBarrelLowPtTight/I");
	data->Branch("muPos_oneStationAngTight",&muPos_oneStationAngTight,"muPos_oneStationAngTight/I");
																      //2). Negative Muon
	//2). Negative Muon
	data->Branch("muNeg_nchi2In", &muNeg_nchi2In, "muNeg_nchi2In/D");
	data->Branch("muNeg_dxy", &muNeg_dxy, "muNeg_dxy/D");
	data->Branch("muNeg_dz", &muNeg_dz, "muNeg_dz/D");
	data->Branch("muNeg_nchi2Gl", &muNeg_nchi2Gl, "muNeg_nchi2Gl/D");
	data->Branch("muNeg_arbitrated", &muNeg_arbitrated, "muNeg_arbitrated/I");
	data->Branch("muNeg_oneStationTight", &muNeg_oneStationTight, "muNeg_oneStationTight/I");
	data->Branch("muNeg_found", &muNeg_found, "muNeg_found/I");
	data->Branch("muNeg_pixeLayers", &muNeg_pixeLayers, "muNeg_pixeLayers/I");
	data->Branch("muNeg_nValidMuHits", &muNeg_nValidMuHits, "muNeg_nValidMuHits/I");
	data->Branch("muNeg_lastStationAngTight",&muNeg_lastStationAngTight, "muNeg_lastStationAngTight/I");
	data->Branch("muNeg_lastStationTight",&muNeg_lastStationTight, "muNeg_lastStationTight/I");
	data->Branch("muNeg_lastStationOptimizedLowPtTight",&muNeg_lastStationOptimizedLowPtTight, "muNeg_lastStatio nOptimizedLowPtTight/I");
	data->Branch("muNeg_lastStationOptimizedBarrelLowPtTight",&muNeg_lastStationOptimizedBarrelLowPtTight, "muNe g_lastStationOptimizedBarrelLowPtTight/I");
	data->Branch("muNeg_oneStationAngTight",&muNeg_oneStationAngTight,"muNeg_oneStationAngTight/I");
	
	nEntries=tree->GetEntries();
	for(int i=0; i<nEntries; i++)
	{
		tree->GetEntry(i);
		if( JpsiMass>massMin && JpsiMass<massMax && JpsiPt>6.) data->Fill();
	}
	cout<<"Entries of Old Tree: "<<tree->GetEntries()<<endl;
	cout<<"Entries of New Tree: "<<data->GetEntries()<<endl;
	outfile->cd();
	data->Write();
	outfile->Close();
	infile->Close();

	
}
