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
#include "RooCBShape.h"
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

  TChain *dataTreesPR = new TChain("data");
  dataTreesPR->Add("/Users/lindseygray/jpsi_workdir/Spring10/promptJpsi/TTree_pol_Mu0Track0Jpsi_MCprompt.root");

  TChain *dataTreesRealData = new TChain("data");
  dataTreesRealData->Add("/Users/lindseygray/jpsi_workdir/data/TTree_pol_Mu0TkMu0Jpsi_dataR_30Aug2010.root");


  Char_t *fileNameInput = "jPsiFit.root";
  TFile* fInput = new TFile(fileNameInput);


  RooDataSet *dataPR = new RooDataSet("dataPR","Supplied Data Prompt",dataTreesPR,varlist);
  RooDataSet *dataReal = new RooDataSet("dataPR","Supplied Data Prompt",dataTreesRealData,varlist);

  RooDataSet *dataPRbin;
  RooDataSet *dataRealbin;

  RooDataHist *MassHist;

  MassModel* massinp_;
  CompositeModelBuilder* theModel;

  RooPlot* REALMassFrame;
  RooPlot* FSRMassFrame;

	char outputfilename[200];
	sprintf(outputfilename,"Results/MassParameterResults.txt");
	printf("output filename is: %s\n", outputfilename);
	FILE *outputFile = fopen(outputfilename,"w");


  RooRealVar *CBm, *CBs, *CBa, *CBn, *FracGauss, *gSig;
  RooRealVar *errCBm, *errCBs, *errCBa, *errCBn, *errFracGauss, *errgSig;

  double CBm_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], CBs_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], CBa_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], CBn_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], FracGauss_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], gSig_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
  double errCBm_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errCBs_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errCBa_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errCBn_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errFracGauss_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errgSig_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];


  for(int ptBin = 0; ptBin < jpsi::kNbPTBins; ++ptBin) {
    for(int yBin = 0; yBin < jpsi::kNbRapForPTBins; ++yBin) {


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

	  dataPRbin = (RooDataSet*)dataPR->reduce(reducePR);
	  dataRealbin = (RooDataSet*)dataReal->reduce(reducePR);

	  JpsiMass.setBins(50);
	  MassHist = new RooDataHist("MassHist","Binned Mass",RooArgSet(JpsiMass),*dataRealbin);

	  char DirectoryPathMass[200];
	  sprintf(DirectoryPathMass,"/pt%d_rapidity%d/MassModel",ptBin+1,yBin+1);

	  TDirectory *InputDirectoryMass = (TDirectory*)fInput->Get(DirectoryPathMass);

	  theModel = new CompositeModelBuilder("CS");
	  theModel->setUseLifetime(false);
	  theModel->setUsePol(false);

	  theModel->getMassModel()->loadParameters(*InputDirectoryMass);
	  theModel->getMassModel()->fix("CBa");
	  theModel->getMassModel()->fix("CBn");
	  theModel->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);

	  //  theModel->getMassModel()->fix("CBm");
	  //  theModel->getMassModel()->fix("CBs");

	  theModel->model()->fitTo(*MassHist,RooFit::NumCPU(2),RooFit::Timer(true),
				   RooFit::Extended(true),RooFit::SumW2Error(true));

	  massinp_ = new MassModel();
	  massinp_->loadParameters(*InputDirectoryMass);
	  massinp_->initModels(JpsiMass);

  FSRMassFrame = new RooPlot;
  FSRMassFrame = JpsiMass.frame(Bins(50)) ;
  dataPRbin->plotOn(FSRMassFrame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
  massinp_->prompt()->plotOn(FSRMassFrame,Normalization(1.0), LineWidth(2));
  double chi2_FSRMassFrame=FSRMassFrame->chiSquare();
  FSRMassFrame->SetMinimum(1);
  char FSRMassTitle[200];
  sprintf(FSRMassTitle,"MC MassFit %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_FSRMassFrame);
  FSRMassFrame->SetTitle(FSRMassTitle);

  REALMassFrame = new RooPlot;
  REALMassFrame = JpsiMass.frame(Bins(50)) ;
  dataRealbin->plotOn(REALMassFrame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
  theModel->model()->plotOn(REALMassFrame, Normalization(1.0), LineWidth(2));
  double chi2_REALMassFrame=REALMassFrame->chiSquare();
  REALMassFrame->SetMinimum(1);
  char REALMassTitle[200];
  sprintf(REALMassTitle,"MC MassFit with real data %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_REALMassFrame);
  REALMassFrame->SetTitle(REALMassTitle);

  TH1* legendBlue = dataPR->createHistogram("legendBlue",JpsiMass,Binning(100)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(3) ;
  TH1* legendRed = dataPR->createHistogram("legendRed",JpsiMass,Binning(100)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(kSolid) ; legendRed->SetLineWidth(2) ;
  TH1* legendBlack = dataPR->createHistogram("legendBlack",JpsiMass,Binning(100)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
  TH1* legendGreen = dataPR->createHistogram("legendGreen",JpsiMass,Binning(100)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;

  TLegend* PromptLegend=new TLegend(0.6,0.6,0.9,0.9);
  PromptLegend->SetFillColor(kWhite);
  PromptLegend->SetTextFont(72);
  PromptLegend->SetTextSize(0.025);
  PromptLegend->AddEntry(legendBlue,"MC Prompt fit","l");
  PromptLegend->AddEntry(legendBlack,"Gauss1","l");
  PromptLegend->AddEntry(legendRed,"Gauss2","l");
  PromptLegend->AddEntry(legendGreen,"Gauss3","l");

  TLegend* NonPromptLegend=new TLegend(0.6,0.6,0.9,0.9);
  NonPromptLegend->SetFillColor(kWhite);
  NonPromptLegend->SetTextFont(72);
  NonPromptLegend->SetTextSize(0.025);
  NonPromptLegend->AddEntry(legendBlue,"MC Non Prompt fit","l");
  NonPromptLegend->AddEntry(legendBlack,"SingleSided Decay","l");
  NonPromptLegend->AddEntry(legendRed,"Flipped Decay","l");

  TLegend* BackgroundLegend=new TLegend(0.6,0.6,0.9,0.9);
  BackgroundLegend->SetFillColor(kWhite);
  BackgroundLegend->SetTextFont(72);
  BackgroundLegend->SetTextSize(0.025);
  BackgroundLegend->AddEntry(legendBlue,"Data background fit","l");
  BackgroundLegend->AddEntry(legendBlack,"SingleSided Decay","l");
  BackgroundLegend->AddEntry(legendRed,"Flipped Decay","l");
  BackgroundLegend->AddEntry(legendGreen,"DoubleSided Decay","l");

  TCanvas* FSRMassCanvas = new TCanvas("FSRMassCanvas","FSRMassCanvas",1400, 600);
  FSRMassCanvas->Divide(2);  FSRMassCanvas->SetFillColor(kWhite);
  FSRMassCanvas->cd(1)->SetLogy(1);
  FSRMassCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); FSRMassFrame->GetYaxis()->SetTitleOffset(2) ; FSRMassFrame->Draw();
  FSRMassCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); FSRMassFrame->GetYaxis()->SetTitleOffset(2) ; FSRMassFrame->Draw();

  TCanvas* REALMassCanvas = new TCanvas("REALMassCanvas","REALMassCanvas",1400, 600);
  REALMassCanvas->Divide(2);  REALMassCanvas->SetFillColor(kWhite);
  REALMassCanvas->cd(1)->SetLogy(1);
  REALMassCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); REALMassFrame->GetYaxis()->SetTitleOffset(2) ; REALMassFrame->Draw();
  REALMassCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); REALMassFrame->GetYaxis()->SetTitleOffset(2) ; REALMassFrame->Draw();

  char FilenamePRMass[200];
  sprintf(FilenamePRMass,"Plots/Mass/massPrompt_pt%d_rapidity%d.png",ptBin+1,yBin+1);
  char FilenameREALMass[200];
  sprintf(FilenameREALMass,"Plots/Mass/REALmassPrompt_pt%d_rapidity%d.png",ptBin+1,yBin+1);

  FSRMassCanvas->SaveAs(FilenamePRMass);
  REALMassCanvas->SaveAs(FilenameREALMass);

  FSRMassCanvas->Close();
  REALMassCanvas->Close();

///////////////////////////////////// EXTRACT AND PRINT PARAMETERS //////////////////////////////////////////////////////

  cout<<"here"<<endl;

	  CBm = (RooRealVar*)InputDirectoryMass->Get("CBm");
	  CBs = (RooRealVar*)InputDirectoryMass->Get("CBs");
	  CBa = (RooRealVar*)InputDirectoryMass->Get("CBa");
	  CBn = (RooRealVar*)InputDirectoryMass->Get("CBn");
	  //FracGauss = (RooRealVar*)InputDirectoryMass->Get("FracGauss");
	  //gSig = (RooRealVar*)InputDirectoryMass->Get("gSig");

	  CBm_[ptBin+1][yBin+1]=CBm->getVal();
	  CBs_[ptBin+1][yBin+1]=CBs->getVal();
	  CBa_[ptBin+1][yBin+1]=CBa->getVal();
	  CBn_[ptBin+1][yBin+1]=CBn->getVal();
	  //FracGauss_[ptBin+1][yBin+1]=FracGauss->getVal();
	  //gSig_[ptBin+1][yBin+1]=gSig->getVal();

	  errCBm_[ptBin+1][yBin+1]=CBm->getError();
	  errCBs_[ptBin+1][yBin+1]=CBs->getError();
	  errCBa_[ptBin+1][yBin+1]=CBa->getError();
	  errCBn_[ptBin+1][yBin+1]=CBn->getError();
	  //errFracGauss_[ptBin+1][yBin+1]=FracGauss->getError();
	  //errgSig_[ptBin+1][yBin+1]=gSig->getError();


	  cout<<"here"<<endl;
	  cout<<"there"<<endl;

  delete REALMassFrame;
  delete FSRMassFrame;
  delete massinp_;
  delete theModel;
  delete MassHist;

	  }
	}

  cout<<"there"<<endl;

	fprintf(outputFile, "\nMC FSR Mass Fit:\n");

	  fprintf(outputFile, "\n  CBm        |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],CBm_[1][1],errCBm_[1][1],CBm_[2][1],errCBm_[2][1],CBm_[3][1],errCBm_[3][1],CBm_[4][1],errCBm_[4][1],CBm_[5][1],errCBm_[5][1],CBm_[6][1],errCBm_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],CBm_[1][2],errCBm_[1][2],CBm_[2][2],errCBm_[2][2],CBm_[3][2],errCBm_[3][2],CBm_[4][2],errCBm_[4][2],CBm_[5][2],errCBm_[5][2],CBm_[6][2],errCBm_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],CBm_[1][3],errCBm_[1][3],CBm_[2][3],errCBm_[2][3],CBm_[3][3],errCBm_[3][3],CBm_[4][3],errCBm_[4][3],CBm_[5][3],errCBm_[5][3],CBm_[6][3],errCBm_[6][3]);

	  fprintf(outputFile, "\n  CBs         |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],CBs_[1][1],errCBs_[1][1],CBs_[2][1],errCBs_[2][1],CBs_[3][1],errCBs_[3][1],CBs_[4][1],errCBs_[4][1],CBs_[5][1],errCBs_[5][1],CBs_[6][1],errCBs_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],CBs_[1][2],errCBs_[1][2],CBs_[2][2],errCBs_[2][2],CBs_[3][2],errCBs_[3][2],CBs_[4][2],errCBs_[4][2],CBs_[5][2],errCBs_[5][2],CBs_[6][2],errCBs_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],CBs_[1][3],errCBs_[1][3],CBs_[2][3],errCBs_[2][3],CBs_[3][3],errCBs_[3][3],CBs_[4][3],errCBs_[4][3],CBs_[5][3],errCBs_[5][3],CBs_[6][3],errCBs_[6][3]);

	  fprintf(outputFile, "\n  CBa        |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],CBa_[1][1],errCBa_[1][1],CBa_[2][1],errCBa_[2][1],CBa_[3][1],errCBa_[3][1],CBa_[4][1],errCBa_[4][1],CBa_[5][1],errCBa_[5][1],CBa_[6][1],errCBa_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],CBa_[1][2],errCBa_[1][2],CBa_[2][2],errCBa_[2][2],CBa_[3][2],errCBa_[3][2],CBa_[4][2],errCBa_[4][2],CBa_[5][2],errCBa_[5][2],CBa_[6][2],errCBa_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],CBa_[1][3],errCBa_[1][3],CBa_[2][3],errCBa_[2][3],CBa_[3][3],errCBa_[3][3],CBa_[4][3],errCBa_[4][3],CBa_[5][3],errCBa_[5][3],CBa_[6][3],errCBa_[6][3]);

	  fprintf(outputFile, "\n  CBn       |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],CBn_[1][1],errCBn_[1][1],CBn_[2][1],errCBn_[2][1],CBn_[3][1],errCBn_[3][1],CBn_[4][1],errCBn_[4][1],CBn_[5][1],errCBn_[5][1],CBn_[6][1],errCBn_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],CBn_[1][2],errCBn_[1][2],CBn_[2][2],errCBn_[2][2],CBn_[3][2],errCBn_[3][2],CBn_[4][2],errCBn_[4][2],CBn_[5][2],errCBn_[5][2],CBn_[6][2],errCBn_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],CBn_[1][3],errCBn_[1][3],CBn_[2][3],errCBn_[2][3],CBn_[3][3],errCBn_[3][3],CBn_[4][3],errCBn_[4][3],CBn_[5][3],errCBn_[5][3],CBn_[6][3],errCBn_[6][3]);
/*
	  fprintf(outputFile, "\n  FracGauss        |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],FracGauss_[1][1],errFracGauss_[1][1],FracGauss_[2][1],errFracGauss_[2][1],FracGauss_[3][1],errFracGauss_[3][1],FracGauss_[4][1],errFracGauss_[4][1],FracGauss_[5][1],errFracGauss_[5][1],FracGauss_[6][1],errFracGauss_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],FracGauss_[1][2],errFracGauss_[1][2],FracGauss_[2][2],errFracGauss_[2][2],FracGauss_[3][2],errFracGauss_[3][2],FracGauss_[4][2],errFracGauss_[4][2],FracGauss_[5][2],errFracGauss_[5][2],FracGauss_[6][2],errFracGauss_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],FracGauss_[1][3],errFracGauss_[1][3],FracGauss_[2][3],errFracGauss_[2][3],FracGauss_[3][3],errFracGauss_[3][3],FracGauss_[4][3],errFracGauss_[4][3],FracGauss_[5][3],errFracGauss_[5][3],FracGauss_[6][3],errFracGauss_[6][3]);

	  fprintf(outputFile, "\n  gSig        |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],gSig_[1][1],errgSig_[1][1],gSig_[2][1],errgSig_[2][1],gSig_[3][1],errgSig_[3][1],gSig_[4][1],errgSig_[4][1],gSig_[5][1],errgSig_[5][1],gSig_[6][1],errgSig_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],gSig_[1][2],errgSig_[1][2],gSig_[2][2],errgSig_[2][2],gSig_[3][2],errgSig_[3][2],gSig_[4][2],errgSig_[4][2],gSig_[5][2],errgSig_[5][2],gSig_[6][2],errgSig_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],gSig_[1][3],errgSig_[1][3],gSig_[2][3],errgSig_[2][3],gSig_[3][3],errgSig_[3][3],gSig_[4][3],errgSig_[4][3],gSig_[5][3],errgSig_[5][3],gSig_[6][3],errgSig_[6][3]);
*/


	fclose(outputFile);

	  cout<<"there"<<endl;

  delete dataTreesRealData;
  delete dataTreesPR;
  delete dataPR;
  delete fInput;



  return 0;
}
