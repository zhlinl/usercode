#include <iostream>
//#include <string>

//#include <sstream>
//J/Psi common vars
#include "commonVar.h"

//Fitting routine
#include "CompositeModelBuilder.h"

// RooFit Includes
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooExponential.h"
#include "RooAbsPdf.h"
#include "RooGaussModel.h"
#include "RooDecay.h"



//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TLegend.h"



int main(int argc, char**argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;


  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
  RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,1.5);//0=GG,1=GT,2=TT
  RooRealVar MCweight("MCweight","MCweight",0,1.1);
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);

  TChain *dataTreesNP = new TChain("data");
  dataTreesNP->Add("/scratch/knuenz/Polarization/RootInput/TTree_NP.root");

  TChain *dataTreesPR = new TChain("data");
  dataTreesPR->Add("/scratch/knuenz/Polarization/RootInput/TTree_red_PR.root");

//  TChain *dataTreesRealData = new TChain("data");
//  dataTreesRealData->Add("/scratch/knuenz/Polarization/RootInput/TTree_pol_Mu0TkMu0Jpsi_data_11Sep2010.root");

   Char_t *fileNameInput = "jPsiFit.root";
   TFile* fInput = new TFile(fileNameInput);



  RooDataSet *dataNP = new RooDataSet("dataNP","Supplied Data NonPrompt",dataTreesNP,varlist);
  RooDataSet *dataPR = new RooDataSet("dataPR","Supplied Data Prompt",dataTreesPR,varlist);
//  RooDataSet *dataBK = new RooDataSet("dataBK","Supplied Data Background",dataTreesBK,varlist);

  RooDataSet *dataNPbin, *dataPRbin, *dataBKbin, *dataREALbin;

  CompositeModelBuilder* cs;
  CompositeModelBuilder* hx;

  RooPlot* CS_cos_Frame;
  RooPlot* CS_phi_Frame;
  RooPlot* HX_cos_Frame;
  RooPlot* HX_phi_Frame;

	char outputfilename[200];
	sprintf(outputfilename,"Results/PolarizationParameterResults.txt");
	printf("output filename is: %s\n", outputfilename);
	FILE *outputFile = fopen(outputfilename,"w");



  for(int ptBin = 3; ptBin < jpsi::kNbPTBins-2; ++ptBin) {
    for(int yBin = 2; yBin < jpsi::kNbRapForPTBins; ++yBin) {




      cout<<"pT"<<ptBin+1<<"rapidity"<<yBin+1<<endl;

	  char reduceNP[200];
	  sprintf(reduceNP,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
	  char reducePR[200];
	  sprintf(reducePR,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
	  char reduceBK[200];
	  sprintf(reduceBK,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f && JpsiMass < %f | JpsiMass > %f",
		  jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],
		  jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgLow,
		  jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgHigh);

	  dataNPbin = (RooDataSet*)dataNP->reduce(reduceNP);
	  dataPRbin = (RooDataSet*)dataPR->reduce(reducePR);
//	  dataBKbin = (RooDataSet*)dataBK->reduce(reduceBK);
//	  dataREALbin = (RooDataSet*)dataBK->reduce(reducePR);

	  TTree *TreedataPRbin = (TTree*)dataPRbin->tree();

	  costh_CS.setBins(jpsi::kNbBinsCosT);
	  phi_CS.setBins(jpsi::kNbBinsPhiPol);

	  costh_HX.setBins(jpsi::kNbBinsCosT);
	  phi_HX.setBins(jpsi::kNbBinsPhiPol);

	  TH2F * hist_CS_PR = dataPRbin->createHistogram(costh_CS,phi_CS);//"hist_CS_PR",costh_CS,Binning(jpsi::kNbBinsCosT),YVar(phi_CS,Binning(jpsi::kNbBinsPhiPol))) ;
	  TH2F * hist_HX_PR = dataPRbin->createHistogram(costh_HX,phi_HX);//"hist_CS_PR",costh_CS,Binning(jpsi::kNbBinsCosT),YVar(phi_CS,Binning(jpsi::kNbBinsPhiPol))) ;

	  char DirectoryPathPol[200];
	  sprintf(DirectoryPathPol,"/pt%d_rapidity%d/PolarizationModel",ptBin+1,yBin+1);

	  TDirectory *InputDirectoryPol = (TDirectory*)fInput->Get(DirectoryPathPol);

	  cs = new CompositeModelBuilder("CS");
	  cs->setUsePrompt(true);
	  cs->setUseNonPrompt(false);
	  cs->setUseMass(false);
	  cs->setUseLifetime(false);
	  cs->setUsePol(true);
//	  cs->getPromptPolarizationModel()->loadParameters(*InputDirectoryPol);
//	  cs->getNonPromptPolarizationModel()->loadParameters(*InputDirectoryPol);
	  cs->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);

	  TH1 * hist_CS_PR_model = cs->model()->createHistogram("hist_CS_PR_model",costh_CS,Binning(jpsi::kNbBinsCosT),YVar(phi_CS,Binning(jpsi::kNbBinsPhiPol))) ;

      hx = new CompositeModelBuilder("HX");
 	  hx->setUsePrompt(true);
 	  hx->setUseNonPrompt(false);
 	  hx->setUseMass(false);
 	  hx->setUseLifetime(false);
 	  hx->setUsePol(true);
//	  hx->getPromptPolarizationModel()->loadParameters(*InputDirectoryPol);
//	  hx->getNonPromptPolarizationModel()->loadParameters(*InputDirectoryPol);
 	  hx->initModel(JpsiMass,Jpsict,costh_HX,phi_HX);

 	  TH1 * hist_HX_PR_model = hx->model()->createHistogram("hist_HX_PR_model",costh_HX,Binning(jpsi::kNbBinsCosT),YVar(phi_HX,Binning(jpsi::kNbBinsPhiPol))) ;


  CS_cos_Frame = new RooPlot;
  CS_cos_Frame = costh_CS.frame(Bins(50)) ;
  dataPRbin->plotOn(CS_cos_Frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
//  cs->model()->plotOn(CS_cos_Frame, Normalization(1.0), LineWidth(3),LineColor(kBlue));
  double chi2_CS_cos_Frame=CS_cos_Frame->chiSquare();
  CS_cos_Frame->SetMinimum(1);
  char CS_cos_Title[200];
  sprintf(CS_cos_Title,"Cos#theta_{CS} MC Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_CS_cos_Frame);
  CS_cos_Frame->SetTitle(CS_cos_Title);

  CS_phi_Frame = new RooPlot;
  CS_phi_Frame = phi_CS.frame(Bins(50)) ;
  dataPRbin->plotOn(CS_phi_Frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
//  cs->model()->plotOn(CS_phi_Frame, Normalization(1.0), LineWidth(3),LineColor(kBlue));
  double chi2_CS_phi_Frame=CS_phi_Frame->chiSquare();
  CS_phi_Frame->SetMinimum(1);
  char CS_phi_Title[200];
  sprintf(CS_phi_Title,"#phi_{CS} MC Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_CS_phi_Frame);
  CS_phi_Frame->SetTitle(CS_phi_Title);


  char CS2DTitle[200];
  sprintf(CS2DTitle,"MC %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
  char CS2DFitTitle[200];
  sprintf(CS2DFitTitle,"Fit %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);


  TCanvas* CS_cos_Canvas = new TCanvas("CS_cos_Canvas","CS_cos_Canvas",1400, 600);
  CS_cos_Canvas->Divide(2);  CS_cos_Canvas->SetFillColor(kWhite);
  CS_cos_Canvas->cd(1)->SetLogy(1);
  CS_cos_Canvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); CS_cos_Frame->GetYaxis()->SetTitleOffset(2) ; CS_cos_Frame->Draw();
  CS_cos_Canvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); CS_cos_Frame->GetYaxis()->SetTitleOffset(2) ; CS_cos_Frame->Draw(); //PromptLegend->Draw();

  TCanvas* CS_phi_Canvas = new TCanvas("CS_phi_Canvas","CS_phi_Canvas",1400, 600);
  CS_phi_Canvas->Divide(2);  CS_phi_Canvas->SetFillColor(kWhite);
  CS_phi_Canvas->cd(1)->SetLogy(1);
  CS_phi_Canvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); CS_phi_Frame->GetYaxis()->SetTitleOffset(2) ; CS_phi_Frame->Draw();
  CS_phi_Canvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); CS_phi_Frame->GetYaxis()->SetTitleOffset(2) ; CS_phi_Frame->Draw(); //PromptLegend->Draw();


  TCanvas* CSPol2DCanvas = new TCanvas("CSPol2DCanvas","CSPol2DCanvas",1400, 700);
  CSPol2DCanvas->Divide(2);  CSPol2DCanvas->SetFillColor(kWhite);
  CSPol2DCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); hist_CS_PR->SetStats(0); hist_CS_PR->GetXaxis()->SetTitle("cos #theta_{CS}"); hist_CS_PR->GetYaxis()->SetTitle("#phi_{CS} [deg]"); hist_CS_PR->GetYaxis()->SetTitleOffset(2); hist_CS_PR->SetTitle(CS2DTitle); hist_CS_PR->Draw("colz");
  CSPol2DCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); hist_CS_PR_model->SetStats(0); hist_CS_PR_model->SetTitle(CS2DFitTitle); hist_CS_PR_model->GetYaxis()->SetTitleOffset(2); hist_CS_PR_model->Draw("colz");

  char FilenameCScos[200];
  sprintf(FilenameCScos,"Plots/Polarization/CS_costh_pt%d_rapidity%d.png",ptBin+1,yBin+1);
  char FilenameCSphi[200];
  sprintf(FilenameCSphi,"Plots/Polarization/CS_phi_pt%d_rapidity%d.png",ptBin+1,yBin+1);
  char FilenameCS2D[200];
  sprintf(FilenameCS2D,"Plots/Polarization/CS_2D_pt%d_rapidity%d.png",ptBin+1,yBin+1);

  CS_phi_Canvas->SaveAs(FilenameCSphi);
  CS_cos_Canvas->SaveAs(FilenameCScos);
  CSPol2DCanvas->SaveAs(FilenameCS2D);

  CS_phi_Canvas->Close();
  CS_cos_Canvas->Close();
  CSPol2DCanvas->Close();



  HX_cos_Frame = new RooPlot;
   HX_cos_Frame = costh_HX.frame(Bins(50)) ;
   dataPRbin->plotOn(HX_cos_Frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
 //  HX->model()->plotOn(HX_cos_Frame, Normalization(1.0), LineWidth(3),LineColor(kBlue));
   double chi2_HX_cos_Frame=HX_cos_Frame->chiSquare();
   HX_cos_Frame->SetMinimum(1);
   char HX_cos_Title[200];
   sprintf(HX_cos_Title,"Cos#theta_{HX} MC Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_HX_cos_Frame);
   HX_cos_Frame->SetTitle(HX_cos_Title);

   HX_phi_Frame = new RooPlot;
   HX_phi_Frame = phi_HX.frame(Bins(50)) ;
   dataPRbin->plotOn(HX_phi_Frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
 //  HX->model()->plotOn(HX_phi_Frame, Normalization(1.0), LineWidth(3),LineColor(kBlue));
   double chi2_HX_phi_Frame=HX_phi_Frame->chiSquare();
   HX_phi_Frame->SetMinimum(1);
   char HX_phi_Title[200];
   sprintf(HX_phi_Title,"#phi_{HX} MC Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_HX_phi_Frame);
   HX_phi_Frame->SetTitle(HX_phi_Title);


   char HX2DTitle[200];
   sprintf(HX2DTitle,"MC %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
   char HX2DFitTitle[200];
   sprintf(HX2DFitTitle,"Fit %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);


   TCanvas* HX_cos_Canvas = new TCanvas("HX_cos_Canvas","HX_cos_Canvas",1400, 600);
   HX_cos_Canvas->Divide(2);  HX_cos_Canvas->SetFillColor(kWhite);
   HX_cos_Canvas->cd(1)->SetLogy(1);
   HX_cos_Canvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); HX_cos_Frame->GetYaxis()->SetTitleOffset(2) ; HX_cos_Frame->Draw();
   HX_cos_Canvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); HX_cos_Frame->GetYaxis()->SetTitleOffset(2) ; HX_cos_Frame->Draw(); //PromptLegend->Draw();

   TCanvas* HX_phi_Canvas = new TCanvas("HX_phi_Canvas","HX_phi_Canvas",1400, 600);
   HX_phi_Canvas->Divide(2);  HX_phi_Canvas->SetFillColor(kWhite);
   HX_phi_Canvas->cd(1)->SetLogy(1);
   HX_phi_Canvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); HX_phi_Frame->GetYaxis()->SetTitleOffset(2) ; HX_phi_Frame->Draw();
   HX_phi_Canvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); HX_phi_Frame->GetYaxis()->SetTitleOffset(2) ; HX_phi_Frame->Draw(); //PromptLegend->Draw();


   TCanvas* HXPol2DCanvas = new TCanvas("HXPol2DCanvas","HXPol2DCanvas",1400, 700);
   HXPol2DCanvas->Divide(2);  HXPol2DCanvas->SetFillColor(kWhite);
   HXPol2DCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); hist_HX_PR->SetStats(0); hist_HX_PR->GetXaxis()->SetTitle("cos #theta_{HX}"); hist_HX_PR->GetYaxis()->SetTitle("#phi_{HX} [deg]"); hist_HX_PR->GetYaxis()->SetTitleOffset(2); hist_HX_PR->SetTitle(HX2DTitle); hist_HX_PR->Draw("colz");
   HXPol2DCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); hist_HX_PR_model->SetStats(0); hist_HX_PR_model->SetTitle(HX2DFitTitle); hist_HX_PR_model->GetYaxis()->SetTitleOffset(2); hist_HX_PR_model->Draw("colz");

   char FilenameHXcos[200];
   sprintf(FilenameHXcos,"Plots/Polarization/HX_costh_pt%d_rapidity%d.png",ptBin+1,yBin+1);
   char FilenameHXphi[200];
   sprintf(FilenameHXphi,"Plots/Polarization/HX_phi_pt%d_rapidity%d.png",ptBin+1,yBin+1);
   char FilenameHX2D[200];
   sprintf(FilenameHX2D,"Plots/Polarization/HX_2D_pt%d_rapidity%d.png",ptBin+1,yBin+1);

   HX_phi_Canvas->SaveAs(FilenameHXphi);
   HX_cos_Canvas->SaveAs(FilenameHXcos);
   HXPol2DCanvas->SaveAs(FilenameHX2D);

   HX_phi_Canvas->Close();
   HX_cos_Canvas->Close();
   HXPol2DCanvas->Close();


	  }
	}




	fclose(outputFile);



  delete dataTreesNP;
  delete dataNP;
  delete dataTreesPR;
  delete dataPR;
//  delete dataTreesBK;
//  delete dataBK;
  delete fInput;



  return 0;
}
