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

#include <TAxis.h>        
#include <TH1.h>          
#include <TTree.h>        
#include <TFile.h>        
#include <TH1D.h>         
#include <TH1I.h>         
#include <TCanvas.h>      
#include <TLine.h>        
#include <TMath.h>        
#include <TVector3.h>     
#include <vector.h>       
#include <TString.h>      
#include <iostream.h>     
#include <TLegend.h>      

//using namespace RooFit;

void drawDistri()
{

	Long_t nEntries, nData=0;  

	Double_t massMin,massMax;  
	massMin=2.6;               
	massMax=3.5;               

	char fileName[500];                                     
	sprintf(fileName,                                       
			//"/home/zhlinl/work/CMSSW_3_8_4/src/HeavyFlavorAnalysis/Onia2MuMu/test/rootfiles/tree/tree_nocut.root"         
			//"/home/zhlinl/data/tree_nocut_runs140042-141914_MuOnia1.root"
			//"/home/zhlinl/data/tree_nocut_runs140042-141914_MuOnia.root"
				"/home/zhlinl/data/muonTree_Run2010B-Nov4ReReco_v1-Onia2MuMu-v_new/muonTree_nocut_1.root"
			);                                                  

	//Jpsi Variables                                        
	Double_t JpsiMass,JpsiPt,JpsiRap;                       
	Double_t JpsiVprob;                                     
	//(1).Positive Muon                                     
	double muPos_nchi2In, muPos_dxy, muPos_dz, muPos_nchi2Gl;                                                         
	int muPos_arbitrated, muPos_oneStaTight;               
	int muPos_found, muPos_pixeLayers, muPos_nValidMuHits;  
	//(2).Negative Muon                                     
	double muNeg_nchi2In, muNeg_dxy, muNeg_dz, muNeg_nchi2Gl;                                                         
	int muNeg_arbitrated, muNeg_oneStaTight;               
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
	//Triger information
	tree->SetBranchAddress("HLT_L1DoubleMuOpen",&HLT_L1DoubleMuOpen);                                                 
	tree->SetBranchAddress("HLT_Mu0_TkMu0_Jpsi",&HLT_Mu0_TkMu0_Jpsi);                                     
	tree->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2",&HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2);
	//1) Positive Muon                                      
	tree->SetBranchAddress("muPos_nchi2In", &muPos_nchi2In);
	tree->SetBranchAddress("muPos_dxy", &muPos_dxy);        
	tree->SetBranchAddress("muPos_dz", &muPos_dz);          
	tree->SetBranchAddress("muPos_nchi2Gl", &muPos_nchi2Gl);
	tree->SetBranchAddress("muPos_arbitrated", &muPos_arbitrated);                                                    
	tree->SetBranchAddress("muPos_oneStaTight", &muPos_oneStaTight);                                                  
	tree->SetBranchAddress("muPos_found", &muPos_found);    
	tree->SetBranchAddress("muPos_pixeLayers", &muPos_pixeLayers);                                                    
	tree->SetBranchAddress("muPos_nValidMuHits", &muPos_nValidMuHits);                                                
	//2) Negative Muon                                      
	tree->SetBranchAddress("muNeg_nchi2In", &muNeg_nchi2In);
	tree->SetBranchAddress("muNeg_dxy", &muNeg_dxy);        
	tree->SetBranchAddress("muNeg_dz", &muNeg_dz);          
	tree->SetBranchAddress("muNeg_nchi2Gl", &muNeg_nchi2Gl);
	tree->SetBranchAddress("muNeg_arbitrated", &muNeg_arbitrated);                                                    
	tree->SetBranchAddress("muNeg_oneStaTight", &muNeg_oneStaTight);                                                  
	tree->SetBranchAddress("muNeg_found", &muNeg_found);    
	tree->SetBranchAddress("muNeg_pixeLayers", &muNeg_pixeLayers);                                                    
	tree->SetBranchAddress("muNeg_nValidMuHits", &muNeg_nValidMuHits);                                                

	TH1D *vProb=new TH1D("vProb","dimuon vertex probability",100,0.,0.02);


	TH1D *muonPos_nchi2In=new TH1D("muonPos_nchi2In","mu inner tracker normalized chi2",200,0.,10.);
	TH1D *muonPos_dxy=new TH1D("muonPos_dxy","mu dxy",100,-4.0,4.0);
	TH1D *muonPos_dz=new TH1D("muonPos_dz","mu dz",100,-20.,20.);
	TH1D *muonPos_nchi2Gl=new TH1D("muonPos_nchi2Gl","mu global tracker normalized chi2",100,0.,30.);
	TH1I *muonPos_arbitrated=new TH1I("muonPos_arbitrated","TrackerMuonArbitrated",3,0,3);
	TH1I *muonPos_oneStaTight=new TH1I("muonPos_oneStaTight","TMOneStationTight",3,0,3);
	TH1I *muonPos_found=new TH1I("muonPos_found","mu innerTrack Hits ",100,0,35);
	TH1I *muonPos_pixeLayers=new TH1I("muonPos_pixeLayers","mu pixeLayers",100,0,6);
	TH1I *muonPos_nValidMuHits=new TH1I("muonPos_nValidMuHits","mu nValidMuHits",100,0,55);
	//2) Negative Muon                                    
	TH1D *muonNeg_nchi2In=new TH1D("muonNeg_nchi2In","mu inner tracker normalized chi2",200,0.,10.);
	TH1D *muonNeg_dxy=new TH1D("muonNeg_dxy","mu dxy",100,-4.0,4.0);                                  
	TH1D *muonNeg_dz=new TH1D("muonNeg_dz","mu dz",100,-20.,20.);
	TH1D *muonNeg_nchi2Gl=new TH1D("muonNeg_nchi2Gl","mu global tracker normalized chi2",100,0.,30.);
	TH1I *muonNeg_arbitrated=new TH1I("muonNeg_arbitrated","TrackerMuonArbitrated",3,0,3);
	TH1I *muonNeg_oneStaTight=new TH1I("muonNeg_oneStaTight","TMOneStationTight",3,0,3);              
	TH1I *muonNeg_found=new TH1I("muonNeg_found","mu innerTrack Hits",100,0,35);
	TH1I *muonNeg_pixeLayers=new TH1I("muonNeg_pixeLayers","mu pixeLayers",100,0,6);
	TH1I *muonNeg_nValidMuHits=new TH1I("muonNeg_nValidMuHits","mu nValidMuHits",100,0,55);
	cout<<"===============end definition of Histo==============="<<endl;

	nEntries=tree->GetEntries();

	cout<<"===============total Entries in Tree: "<<nEntries<<"==============="<<endl;
	for(int i=0; i<nEntries; i++)
	{
		tree->GetEntry(i);
		if( //JpsiMass>massMin && JpsiMass<massMax
				//&& HLT_L1DoubleMuOpen == 1 
			 	//&& HLT_Mu0_TkMu0_Jpsi == 2
				HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2 ==1
				&& JpsiPt > 6
				//&& muPos_arbitrated ==1 && muNeg_arbitrated ==1
				//&& muPos_oneStaTight ==1 && muNeg_oneStaTight ==1
				)

		{
			if(muPos_nchi2In > -5000.)
			{
				vProb->Fill(JpsiVprob);
				
				muonPos_nchi2In->Fill(muPos_nchi2In);
				muonPos_dxy->Fill(muPos_dxy);
				muonPos_dz->Fill(muPos_dz);
				if(muPos_nchi2Gl>-5000.)
				{
					muonPos_nchi2Gl->Fill(muPos_nchi2Gl);
					muonPos_nValidMuHits->Fill(muPos_nValidMuHits);
				}
				muonPos_arbitrated->Fill(muPos_arbitrated);
				muonPos_oneStaTight->Fill(muPos_oneStaTight);
				muonPos_found->Fill(muPos_found);
				muonPos_pixeLayers->Fill(muPos_pixeLayers);

				muonNeg_nchi2In->Fill(muNeg_nchi2In);
				muonNeg_dxy->Fill(muNeg_dxy);
				muonNeg_dz->Fill(muNeg_dz);
				if(muNeg_nchi2Gl>-5000.)
				{
					muonNeg_nchi2Gl->Fill(muNeg_nchi2Gl);
					muonNeg_nValidMuHits->Fill(muNeg_nValidMuHits);
				}
				muonNeg_arbitrated->Fill(muNeg_arbitrated);
				muonNeg_oneStaTight->Fill(muNeg_oneStaTight);
				muonNeg_found->Fill(muNeg_found);
				muonNeg_pixeLayers->Fill(muNeg_pixeLayers);
			}
			nData++;
		}
	}
	cout<<"======total matched Entries: "<<nData<<" ======"<<endl;

	TCanvas *c0=new TCanvas("c0","",100,100,900,500);
	c0->cd();
	c0->SetFillColor(10);
	c0->SetLogy();
	vProb->SetXTitle("dimuon Vertex Prob");
	vProb->Draw();
	c0->Modified();
	c0->Update();
	c0->SaveAs("pic/dist/vProb.eps");
	c0->SaveAs("pic/dist/vProb.gif");

	TCanvas *c1=new TCanvas("c1","",100,100,900,500);
	c1->cd();
	c1->SetFillColor(10);
	c1->SetLogy();
	muonPos_nchi2In->SetXTitle("innerTrack chi2/ndof");
	muonPos_nchi2In->Draw();
	muonNeg_nchi2In->SetLineColor(kRed);
	muonNeg_nchi2In->Draw("same");
	TLegend *leg1=new TLegend(0.7,0.6,0.9,0.8);
	leg1->AddEntry(muonPos_nchi2In,"Pos muon");
	leg1->AddEntry(muonNeg_nchi2In,"Neg muon");
	leg1->Draw();
	c1->Modified();
	c1->Update();
	c1->SaveAs("pic/dist/mu_nchi2In.eps");
	c1->SaveAs("pic/dist/mu_nchi2In.gif");

	TCanvas *c2=new TCanvas("c2","",100,100,900,500);
	c2->cd();
	c2->SetFillColor(10);
	c2->SetLogy();
	muonPos_dxy->SetXTitle("mu_dxy");
	muonPos_dxy->Draw();
	muonNeg_dxy->SetLineColor(kRed);
	muonNeg_dxy->Draw("same");
	TLegend *leg2=new TLegend(0.7,0.6,0.9,0.8);
	leg2->AddEntry(muonPos_dxy,"Pos muon");
	leg2->AddEntry(muonNeg_dxy,"Neg muon");
	leg2->Draw();
	c2->Modified();
	c2->Update();
	c2->SaveAs("pic/dist/mu_dxy.eps");
	c2->SaveAs("pic/dist/mu_dxy.gif");

	TCanvas *c3=new TCanvas("c3","",100,100,900,500);
	c3->cd();
	c3->SetFillColor(10);
	c3->SetLogy();
	muonPos_dz->SetXTitle("mu_dz");
	muonPos_dz->Draw();
	muonNeg_dz->SetLineColor(kRed);
	muonNeg_dz->Draw("same");
	TLegend *leg3=new TLegend(0.7,0.6,0.9,0.8);
	leg3->AddEntry(muonPos_dz,"Pos muon");
	leg3->AddEntry(muonNeg_dz,"Neg muon");
	leg3->Draw();
	c3->Modified();
	c3->Update();
	c3->SaveAs("pic/dist/mu_dz.eps");
	c3->SaveAs("pic/dist/mu_dz.gif");

	TCanvas *c4=new TCanvas("c4","",100,100,900,500);
	c4->cd();
	c4->SetFillColor(10);
	c4->SetLogy();
	muonPos_nchi2Gl->SetXTitle("globalTrack chi2/ndof");
	muonPos_nchi2Gl->Draw();
	muonNeg_nchi2Gl->SetLineColor(kRed);
	muonNeg_nchi2Gl->Draw("same");
	TLegend *leg4=new TLegend(0.7,0.6,0.9,0.8);
	leg4->AddEntry(muonPos_nchi2Gl,"Pos muon");
	leg4->AddEntry(muonNeg_nchi2Gl,"Neg muon");
	leg4->Draw();
	c4->Modified();
	c4->Update();
	c4->SaveAs("pic/dist/mu_nchi2Gl.eps");
	c4->SaveAs("pic/dist/mu_nchi2Gl.gif");

	TCanvas *c5=new TCanvas("c5","",100,100,900,500);
	c5->cd();
	c5->SetFillColor(10);
	muonPos_arbitrated->SetXTitle("ID TrackerMuonArbitrated");
	muonPos_arbitrated->Draw();
	muonNeg_arbitrated->SetLineColor(kRed);
	muonNeg_arbitrated->Draw("same");
	TLegend *leg5=new TLegend(0.7,0.6,0.9,0.8);
	leg5->AddEntry(muonPos_arbitrated,"Pos muon");
	leg5->AddEntry(muonNeg_arbitrated,"Neg muon");
	leg5->Draw();
	c5->Modified();
	c5->Update();
	c5->SaveAs("pic/dist/mu_arbitrated.eps");
	c5->SaveAs("pic/dist/mu_arbitrated.gif");

	TCanvas *c6=new TCanvas("c6","",100,100,900,500);
	c6->cd();
	c6->SetFillColor(10);
	muonPos_oneStaTight->SetXTitle("ID OneStationTight");
	muonPos_oneStaTight->Draw();
	muonNeg_oneStaTight->SetLineColor(kRed);
	muonNeg_oneStaTight->Draw("same");
	TLegend *leg6=new TLegend(0.7,0.6,0.9,0.8);
	leg6->AddEntry(muonPos_oneStaTight,"Pos muon");
	leg6->AddEntry(muonNeg_oneStaTight,"Neg muon");
	leg6->Draw();
	c6->Modified();
	c6->Update();
	c6->SaveAs("pic/dist/mu_oneStaTight.eps");
	c6->SaveAs("pic/dist/mu_oneStaTight.gif");

	TCanvas *c7=new TCanvas("c7","",100,100,900,500);
	c7->cd();
	c7->SetFillColor(10);
	muonPos_found->SetXTitle("innerTrack Hits");
	muonPos_found->Draw();
	muonNeg_found->SetLineColor(kRed);
	muonNeg_found->Draw("same");
	TLegend *leg7=new TLegend(0.7,0.6,0.9,0.8);
	leg7->AddEntry(muonPos_found,"Pos muon");
	leg7->AddEntry(muonNeg_found,"Neg muon");
	leg7->Draw();
	c7->Modified();
	c7->Update();
	c7->SaveAs("pic/dist/mu_found.eps");
	c7->SaveAs("pic/dist/mu_found.gif");

	TCanvas *c8=new TCanvas("c8","",100,100,900,500);
	c8->cd();
	c8->SetFillColor(10);
	muonPos_pixeLayers->SetXTitle("mu_PixelLayers");
	muonPos_pixeLayers->Draw();
	muonNeg_pixeLayers->SetLineColor(kRed);
	muonNeg_pixeLayers->Draw("same");
	TLegend *leg8=new TLegend(0.7,0.6,0.9,0.8);
	leg8->AddEntry(muonPos_pixeLayers,"Pos muon");
	leg8->AddEntry(muonNeg_pixeLayers,"Neg muon");
	leg8->Draw();
	c8->Modified();
	c8->Update();
	c8->SaveAs("pic/dist/mu_pixeLayers.eps");
	c8->SaveAs("pic/dist/mu_pixeLayers.gif");

	TCanvas *c9=new TCanvas("c9","",100,100,900,500);
	c9->cd();
	c9->SetFillColor(10);
	c9->SetLogy();
	muonPos_nValidMuHits->SetXTitle("nValidMuHits");
	muonPos_nValidMuHits->Draw();
	muonNeg_nValidMuHits->SetLineColor(kRed);
	muonNeg_nValidMuHits->Draw("same");
	//muonNeg_nValidMuHits->Draw();
	//muonPos_nValidMuHits->Draw("same");
	TLegend *leg9=new TLegend(0.7,0.6,0.9,0.8);
	leg9->AddEntry(muonPos_nValidMuHits,"Pos muon");
	leg9->AddEntry(muonNeg_nValidMuHits,"Neg muon");
	leg9->Draw();
	c9->Modified();
	c9->Update();
	c9->SaveAs("pic/dist/mu_nValidMuHits.eps");
	c9->SaveAs("pic/dist/mu_nValidMuHits.gif");

}
