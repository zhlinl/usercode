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

  bool pereverr(false);

  for(int i=0;i<argc; ++i) {
    if(std::string(argv[i]).find("--perEventErrors") != std::string::npos) pereverr = true;
    //if(std::string(argv[i]).find("--newTTree") != std::string::npos) newdata = true;
  }

  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
  RooRealVar JpsictErr("JpsictErr","Error on l_{J/#psi} [mm]",1e-6,(1-1e-6));
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
  RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,2.5);//0=GG,1=GT,2=TT
  RooRealVar MCweight("MCweight","MCweight",0,1.1);
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  if(pereverr)
    varlist.add(JpsictErr);


  TChain *realData = new TChain("data");
  realData->Add("/Users/lindseygray/jpsi_workdir/data/TTree_pol_Mu0TkMu0Jpsi_dataR_30Aug2010.root");

  RooDataSet *dataBK =  new RooDataSet("dataBK","Supplied Data Background",realData,varlist);


  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);
  varlist.add(MCweight);

  TChain *dataTreesNP = new TChain("data");
  dataTreesNP->Add("/Users/lindseygray/jpsi_workdir/Spring10/nonPromptJpsi/TTree_pol_Mu0Track0Jpsi_B0ToJPsiMumu.root");
  dataTreesNP->Add("/Users/lindseygray/jpsi_workdir/Spring10/nonPromptJpsi/TTree_pol_Mu0Track0Jpsi_BpToJPsiMuMu.root");
  dataTreesNP->Add("/Users/lindseygray/jpsi_workdir/Spring10/nonPromptJpsi/TTree_pol_Mu0Track0Jpsi_BsToJPsiMuMu.root");
  dataTreesNP->Add("/Users/lindseygray/jpsi_workdir/Spring10/nonPromptJpsi/TTree_pol_Mu0Track0Jpsi_LambdaBToJPsiMuMu.root");

  TChain *dataTreesPR = new TChain("data");
  dataTreesPR->Add("/Users/lindseygray/jpsi_workdir/Spring10/promptJpsi/TTree_pol_Mu0Track0Jpsi_MCprompt.root");

   Char_t *fileNameInput = "jPsiFit.root";
   TFile* fInput = new TFile(fileNameInput);



   RooDataSet *dataNP = new RooDataSet("dataNP","Supplied Data NonPrompt",dataTreesNP,varlist,0,"MCweight");
   RooDataSet *dataPR = new RooDataSet("dataPR","Supplied Data Prompt",dataTreesPR,varlist,0,"MCweight");
 
  RooDataSet *dataNPbin, *dataPRbin, *dataBKbin, *dataREALbin;

  LifetimeModel* lifetinp_;

  RooPlot* BackgroundFrame;
  RooPlot* NonPromptFrame;
  RooPlot* PromptFrame;

	char outputfilename[200];
	sprintf(outputfilename,"Results/LifetimeParameterResults.txt");
	printf("output filename is: %s\n", outputfilename);
	FILE *outputFile = fopen(outputfilename,"w");


  RooRealVar *ctCoefPrompt1, *ctCoefPrompt2, *ctCoefPrompt3, *ctCoefPromptTot, *ctMeanPrompt1, *ctMeanPrompt2, 
    *ctMeanPrompt3, *ctSigmaPrompt1, *ctSigmaPrompt2, *ctSigmaPrompt3;
  RooRealVar *ssTauNonPrompt, *fTauNonPrompt, *dsTauNonPrompt, *ssCoefNonPrompt, *fCoefNonPrompt;
  RooRealVar *ssTauBkg, *fTauBkg, *dsTauBkg, *ssCoefBkg, *fCoefBkg, *dsCoefBkg;

  double ctCoefPrompt1_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], ctCoefPrompt2_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    ctCoefPrompt3_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], ctCoefPromptTot_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    ctMeanPrompt1_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], ctMeanPrompt2_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    ctMeanPrompt3_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], ctSigmaPrompt1_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    ctSigmaPrompt2_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], ctSigmaPrompt3_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
  double errctCoefPrompt1_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errctCoefPrompt2_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    errctCoefPrompt3_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errctCoefPromptTot_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    errctMeanPrompt1_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errctMeanPrompt2_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    errctMeanPrompt3_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errctSigmaPrompt1_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    errctSigmaPrompt2_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errctSigmaPrompt3_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
  double ssTauNonPrompt_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], fTauNonPrompt_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    dsTauNonPrompt_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], ssCoefNonPrompt_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    fCoefNonPrompt_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
  double errssTauNonPrompt_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errfTauNonPrompt_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    errdsTauNonPrompt_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errssCoefNonPrompt_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    errfCoefNonPrompt_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
  double ssTauBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], fTauBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    dsTauBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], ssCoefBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    fCoefBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], dsCoefBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];
  double errssTauBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errfTauBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    errdsTauBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errssCoefBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], 
    errfCoefBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1], errdsCoefBkg_[jpsi::kNbPTBins+1][jpsi::kNbRapForPTBins+1];


  for(int ptBin = 0; ptBin < jpsi::kNbPTBins; ++ptBin) {
    for(int yBin = 0; yBin < jpsi::kNbRapForPTBins; ++yBin) {

      //if(ptBin == 0 && yBin == 0) continue;

      cout<<"pT"<<ptBin+1<<"rapidity"<<yBin+1<<endl;

	  char reduceNP[200];
	  sprintf(reduceNP,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[ptBin],
		  jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
	  char reducePR[200];
	  sprintf(reducePR,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[ptBin],
		  jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
	  char reduceBK[200];
	  sprintf(reduceBK,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f && JpsiMass < %f | JpsiMass > %f",
		  jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],
		  jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgLow,
		  jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgHigh);

	  dataNPbin = (RooDataSet*)dataNP->reduce(reduceNP);
	  dataPRbin = (RooDataSet*)dataPR->reduce(reducePR);
	  dataBKbin = (RooDataSet*)dataBK->reduce(reduceBK);
//	  dataREALbin = (RooDataSet*)dataBK->reduce(reducePR);

	  char DirectoryPath[200];
	  sprintf(DirectoryPath,"/pt%d_rapidity%d/LifetimeModel",ptBin+1,yBin+1);

	  TDirectory *InputDirectory = (TDirectory*)fInput->Get(DirectoryPath);

	  lifetinp_ = new LifetimeModel();
	  lifetinp_->loadParameters(*InputDirectory);
	  if(pereverr) 
	    lifetinp_->initModels(Jpsict,JpsictErr);
	  else
	    lifetinp_->initModels(Jpsict);

	  PromptFrame = new RooPlot;
	  PromptFrame = Jpsict.frame(Bins(50)) ;
	  dataPRbin->plotOn(PromptFrame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	  if(pereverr) {	    
	    lifetinp_->prompt()->plotOn(PromptFrame, Normalization(1.0), 
					RooFit::ProjWData(RooArgSet(JpsictErr),*dataPRbin),
					LineWidth(3),LineColor(kBlue));
	  } else {
	    lifetinp_->prompt()->plotOn(PromptFrame, Normalization(1.0), LineWidth(3),LineColor(kBlue));
	  }
	  double chi2_PromptFrame=PromptFrame->chiSquare();
	  if(pereverr) {
	    /*
	      lifetinp_->prompt()->plotOn(PromptFrame, Components(*lifetinp_->promptGauss1()),  
	      RooFit::ProjWData(RooArgSet(JpsictErr),*dataPRbin),
	      RooFit::LineStyle(kSolid),RooFit::LineColor(kBlack),LineWidth(2));
	      lifetinp_->prompt()->plotOn(PromptFrame, Components(*lifetinp_->promptGauss2()),  
	      RooFit::ProjWData(RooArgSet(JpsictErr),*dataPRbin),
	      RooFit::LineStyle(kSolid),RooFit::LineColor(kRed),LineWidth(2));
	      lifetinp_->prompt()->plotOn(PromptFrame, Components(*lifetinp_->promptGauss3()),  
	      RooFit::ProjWData(RooArgSet(JpsictErr),*dataPRbin),
	      RooFit::LineStyle(kSolid),RooFit::LineColor(kGreen),LineWidth(2));
	    */
	  } else {
	    lifetinp_->prompt()->plotOn(PromptFrame, Components(*lifetinp_->promptGauss1()), 
					RooFit::LineStyle(kSolid),RooFit::LineColor(kBlack),LineWidth(2));
	    lifetinp_->prompt()->plotOn(PromptFrame, Components(*lifetinp_->promptGauss2()), 
					RooFit::LineStyle(kSolid),RooFit::LineColor(kRed),LineWidth(2));
	    lifetinp_->prompt()->plotOn(PromptFrame, Components(*lifetinp_->promptGauss3()), 
					RooFit::LineStyle(kSolid),RooFit::LineColor(kGreen),LineWidth(2));
	  }
	  PromptFrame->SetMinimum(1e-3);
	  char PromptTitle[200];
	  sprintf(PromptTitle,"MC Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_PromptFrame);
	  PromptFrame->SetTitle(PromptTitle);

	  NonPromptFrame = new RooPlot;
	  NonPromptFrame = Jpsict.frame(Bins(50)) ;
	  dataNPbin->plotOn(NonPromptFrame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	  if(pereverr) {
	    lifetinp_->nonPrompt()->plotOn(NonPromptFrame,Normalization(1.0),  					    
					   RooFit::ProjWData(RooArgSet(JpsictErr),*dataNPbin),
					   LineWidth(3));
	  } else {
	    lifetinp_->nonPrompt()->plotOn(NonPromptFrame,Normalization(1.0), LineWidth(3));
	  }

	  double chi2_NonPromptFrame=NonPromptFrame->chiSquare();
//	  lifetinp_->nonPrompt()->plotOn(NonPromptFrame, Components(*lifetinp_->nonPromptDecayss()), RooFit::LineStyle(kSolid),RooFit::LineColor(kBlack),LineWidth(2));
//	  lifetinp_->nonPrompt()->plotOn(NonPromptFrame, Components(*lifetinp_->nonPromptDecayf()), RooFit::LineStyle(kSolid),RooFit::LineColor(kRed),LineWidth(2));
	  NonPromptFrame->SetMinimum(1e-3);
	  char NonPromptTitle[200];
	  sprintf(NonPromptTitle,"MC Non Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_NonPromptFrame);
	  NonPromptFrame->SetTitle(NonPromptTitle);

	  BackgroundFrame = new RooPlot;
	  BackgroundFrame = Jpsict.frame(Bins(50)) ;
	  dataBKbin->plotOn(BackgroundFrame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	  if(pereverr) {
	    lifetinp_->background()->plotOn(BackgroundFrame,Normalization(1.0),  					    
					    RooFit::ProjWData(RooArgSet(JpsictErr),*dataBKbin),
					    LineWidth(3));
	  } else {
	    lifetinp_->background()->plotOn(BackgroundFrame,Normalization(1.0), 
					    LineWidth(3));
	  }
	  double chi2_BackgroundFrame=BackgroundFrame->chiSquare();
	  
	  lifetinp_->background()->plotOn(BackgroundFrame, Components(*lifetinp_->backgroundDecayss()), 					    
					  RooFit::ProjWData(RooArgSet(JpsictErr),*dataBKbin),
					  RooFit::LineStyle(kSolid),RooFit::LineColor(kBlack),LineWidth(2));
	  lifetinp_->background()->plotOn(BackgroundFrame, Components(*lifetinp_->backgroundDecayf()),  					    
					  RooFit::ProjWData(RooArgSet(JpsictErr),*dataBKbin),
					  RooFit::LineStyle(kSolid),RooFit::LineColor(kRed),LineWidth(2));
	  lifetinp_->background()->plotOn(BackgroundFrame, Components(*lifetinp_->backgroundDecayds()), 					    
					  RooFit::ProjWData(RooArgSet(JpsictErr),*dataBKbin),
					  RooFit::LineStyle(kSolid),RooFit::LineColor(kGreen),LineWidth(2));
	  BackgroundFrame->SetMinimum(1e-3);
	  char BackgroundTitle[200];
	  sprintf(BackgroundTitle,"Data Background Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_BackgroundFrame);
	  BackgroundFrame->SetTitle(BackgroundTitle);


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


  TCanvas* NonPromptCanvas = new TCanvas("NonPromptCanvas","NonPromptCanvas",1400, 600);
  NonPromptCanvas->Divide(2);  NonPromptCanvas->SetFillColor(kWhite);
  NonPromptCanvas->cd(1)->SetLogy(1);
  NonPromptCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); NonPromptFrame->GetYaxis()->SetTitleOffset(2) ; NonPromptFrame->Draw();
  NonPromptCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); NonPromptFrame->GetYaxis()->SetTitleOffset(2) ; NonPromptFrame->Draw(); NonPromptLegend->Draw();

  TCanvas* PromptCanvas = new TCanvas("PromptCanvas","PromptCanvas",1400, 600);
  PromptCanvas->Divide(2);  PromptCanvas->SetFillColor(kWhite);
  PromptCanvas->cd(1)->SetLogy(1);
  PromptCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); PromptFrame->GetYaxis()->SetTitleOffset(2) ; PromptFrame->Draw();
  PromptCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); PromptFrame->GetYaxis()->SetTitleOffset(2) ; PromptFrame->Draw(); PromptLegend->Draw();

  TCanvas* BackgroundCanvas = new TCanvas("BackgroundCanvas","BackgroundCanvas",1400, 600);
  BackgroundCanvas->Divide(2);  BackgroundCanvas->SetFillColor(kWhite);
  BackgroundCanvas->cd(1)->SetLogy(1);
  BackgroundCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); BackgroundFrame->GetYaxis()->SetTitleOffset(2) ; BackgroundFrame->Draw();
  BackgroundCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); BackgroundFrame->GetYaxis()->SetTitleOffset(2) ; BackgroundFrame->Draw(); BackgroundLegend->Draw();


  char FilenameNP[200];
  sprintf(FilenameNP,"Plots/Lifetime/ctNonPrompt_pt%d_rapidity%d.png",ptBin+1,yBin+1);
  char FilenamePR[200];
  sprintf(FilenamePR,"Plots/Lifetime/ctPrompt_pt%d_rapidity%d.png",ptBin+1,yBin+1);
  char FilenameBK[200];
  sprintf(FilenameBK,"Plots/Lifetime/data_ctBackground_pt%d_rapidity%d.png",ptBin+1,yBin+1);
  NonPromptCanvas->SaveAs(FilenameNP);
  PromptCanvas->SaveAs(FilenamePR);
  BackgroundCanvas->SaveAs(FilenameBK);

  NonPromptCanvas->Close();
  PromptCanvas->Close();
  BackgroundCanvas->Close();

////////////////////////////////////// EXTRACT AND PRINT PARAMETERS //////////////////////////////////////////////////////

  cout<<"here"<<endl;
  /*
	  ctCoefPrompt1 = (RooRealVar*)InputDirectory->Get("ctCoefPrompt1");
	  ctCoefPrompt2 = (RooRealVar*)InputDirectory->Get("ctCoefPrompt2");
	  ctCoefPrompt3 = (RooRealVar*)InputDirectory->Get("ctCoefPrompt3");
	  ctCoefPromptTot = (RooRealVar*)InputDirectory->Get("ctCoefPromptTot");
	  ctMeanPrompt1 = (RooRealVar*)InputDirectory->Get("ctMeanPrompt1");
	  ctMeanPrompt2 = (RooRealVar*)InputDirectory->Get("ctMeanPrompt2");
	  ctMeanPrompt3 = (RooRealVar*)InputDirectory->Get("ctMeanPrompt3");
	  ctSigmaPrompt1 = (RooRealVar*)InputDirectory->Get("ctSigmaPrompt1");
	  ctSigmaPrompt2 = (RooRealVar*)InputDirectory->Get("ctSigmaPrompt2");
	  ctSigmaPrompt3 = (RooRealVar*)InputDirectory->Get("ctSigmaPrompt3");

	  ssTauNonPrompt = (RooRealVar*)InputDirectory->Get("ssTauNonPrompt");
	  fTauNonPrompt = (RooRealVar*)InputDirectory->Get("fTauNonPrompt");
	  dsTauNonPrompt = (RooRealVar*)InputDirectory->Get("dsTauNonPrompt");
	  ssCoefNonPrompt = (RooRealVar*)InputDirectory->Get("ssCoefNonPrompt");
	  fCoefNonPrompt = (RooRealVar*)InputDirectory->Get("fCoefNonPrompt");

	  ssTauBkg = (RooRealVar*)InputDirectory->Get("ssTauBkg");
	  fTauBkg = (RooRealVar*)InputDirectory->Get("fTauBkg");
	  dsTauBkg = (RooRealVar*)InputDirectory->Get("dsTauBkg");
	  ssCoefBkg = (RooRealVar*)InputDirectory->Get("ssCoefBkg");
	  fCoefBkg = (RooRealVar*)InputDirectory->Get("fCoefBkg");
	  dsCoefBkg = (RooRealVar*)InputDirectory->Get("dsCoefBkg");

	  ctCoefPrompt1_[ptBin+1][yBin+1]=ctCoefPrompt1->getVal();
	  ctCoefPrompt2_[ptBin+1][yBin+1]=ctCoefPrompt2->getVal();
	  ctCoefPrompt3_[ptBin+1][yBin+1]=ctCoefPrompt3->getVal();
	  ctCoefPromptTot_[ptBin+1][yBin+1]=ctCoefPromptTot->getVal();
	  ctMeanPrompt1_[ptBin+1][yBin+1]=ctMeanPrompt1->getVal();
	  ctMeanPrompt3_[ptBin+1][yBin+1]=ctMeanPrompt3->getVal();
	  ctSigmaPrompt1_[ptBin+1][yBin+1]=ctSigmaPrompt1->getVal();
	  ctSigmaPrompt2_[ptBin+1][yBin+1]=ctSigmaPrompt2->getVal();
	  ctSigmaPrompt3_[ptBin+1][yBin+1]=ctSigmaPrompt3->getVal();

	  ssTauNonPrompt_[ptBin+1][yBin+1]=ssTauNonPrompt->getVal();
	  fTauNonPrompt_[ptBin+1][yBin+1]=fTauNonPrompt->getVal();
	  dsTauNonPrompt_[ptBin+1][yBin+1]=dsTauNonPrompt->getVal();
	  ssCoefNonPrompt_[ptBin+1][yBin+1]=ssCoefNonPrompt->getVal();
	  fCoefNonPrompt_[ptBin+1][yBin+1]=fCoefNonPrompt->getVal();

	  ssTauBkg_[ptBin+1][yBin+1]=ssTauBkg->getVal();
	  fTauBkg_[ptBin+1][yBin+1]=fTauBkg->getVal();
	  dsTauBkg_[ptBin+1][yBin+1]=dsTauBkg->getVal();
	  ssCoefBkg_[ptBin+1][yBin+1]=ssCoefBkg->getVal();
	  fCoefBkg_[ptBin+1][yBin+1]=fCoefBkg->getVal();
	  dsCoefBkg_[ptBin+1][yBin+1]=dsCoefBkg->getVal();

	  cout<<"here"<<endl;

	  errctCoefPrompt1_[ptBin+1][yBin+1]=ctCoefPrompt1->getError();
	  errctCoefPrompt2_[ptBin+1][yBin+1]=ctCoefPrompt2->getError();
	  errctCoefPrompt3_[ptBin+1][yBin+1]=ctCoefPrompt3->getError();
	  errctCoefPromptTot_[ptBin+1][yBin+1]=ctCoefPromptTot->getError();
	  errctMeanPrompt1_[ptBin+1][yBin+1]=ctMeanPrompt1->getError();
	  errctMeanPrompt2_[ptBin+1][yBin+1]=ctMeanPrompt1->getError();
	  errctMeanPrompt3_[ptBin+1][yBin+1]=ctMeanPrompt3->getError();
	  errctSigmaPrompt1_[ptBin+1][yBin+1]=ctSigmaPrompt1->getError();
	  errctSigmaPrompt2_[ptBin+1][yBin+1]=ctSigmaPrompt2->getError();
	  errctSigmaPrompt3_[ptBin+1][yBin+1]=ctSigmaPrompt3->getError();

	  errssTauNonPrompt_[ptBin+1][yBin+1]=ssTauNonPrompt->getError();
	  errfTauNonPrompt_[ptBin+1][yBin+1]=fTauNonPrompt->getError();
	  errdsTauNonPrompt_[ptBin+1][yBin+1]=dsTauNonPrompt->getError();
	  errssCoefNonPrompt_[ptBin+1][yBin+1]=ssCoefNonPrompt->getError();
	  errfCoefNonPrompt_[ptBin+1][yBin+1]=fCoefNonPrompt->getError();

	  errssTauBkg_[ptBin+1][yBin+1]=ssTauBkg->getError();
	  errfTauBkg_[ptBin+1][yBin+1]=fTauBkg->getError();
	  errdsTauBkg_[ptBin+1][yBin+1]=dsTauBkg->getError();
	  errssCoefBkg_[ptBin+1][yBin+1]=ssCoefBkg->getError();
	  errfCoefBkg_[ptBin+1][yBin+1]=fCoefBkg->getError();
	  errdsCoefBkg_[ptBin+1][yBin+1]=dsCoefBkg->getError();
  */

	  cout<<"here"<<endl;
	  cout<<"there"<<endl;

  delete BackgroundFrame;
  delete NonPromptFrame;
  delete PromptFrame;

  delete lifetinp_;
  delete PromptLegend;

	  }
	}

  cout<<"there"<<endl; // everywhere?

	  fprintf(outputFile, "Parameter Results of Lifetime Component fits\n");
	  fprintf(outputFile, "\nPrompt Component (MC):\n");

	  fprintf(outputFile, "\n  ctCoefPrompt1   |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ctCoefPrompt1_[1][1],errctCoefPrompt1_[1][1],ctCoefPrompt1_[2][1],errctCoefPrompt1_[2][1],ctCoefPrompt1_[3][1],errctCoefPrompt1_[3][1],ctCoefPrompt1_[4][1],errctCoefPrompt1_[4][1],ctCoefPrompt1_[5][1],errctCoefPrompt1_[5][1],ctCoefPrompt1_[6][1],errctCoefPrompt1_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ctCoefPrompt1_[1][2],errctCoefPrompt1_[1][2],ctCoefPrompt1_[2][2],errctCoefPrompt1_[2][2],ctCoefPrompt1_[3][2],errctCoefPrompt1_[3][2],ctCoefPrompt1_[4][2],errctCoefPrompt1_[4][2],ctCoefPrompt1_[5][2],errctCoefPrompt1_[5][2],ctCoefPrompt1_[6][2],errctCoefPrompt1_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ctCoefPrompt1_[1][3],errctCoefPrompt1_[1][3],ctCoefPrompt1_[2][3],errctCoefPrompt1_[2][3],ctCoefPrompt1_[3][3],errctCoefPrompt1_[3][3],ctCoefPrompt1_[4][3],errctCoefPrompt1_[4][3],ctCoefPrompt1_[5][3],errctCoefPrompt1_[5][3],ctCoefPrompt1_[6][3],errctCoefPrompt1_[6][3]);

	  fprintf(outputFile, "\n  ctCoefPrompt2   |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ctCoefPrompt2_[1][1],errctCoefPrompt2_[1][1],ctCoefPrompt2_[2][1],errctCoefPrompt2_[2][1],ctCoefPrompt2_[3][1],errctCoefPrompt2_[3][1],ctCoefPrompt2_[4][1],errctCoefPrompt2_[4][1],ctCoefPrompt2_[5][1],errctCoefPrompt2_[5][1],ctCoefPrompt2_[6][1],errctCoefPrompt2_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ctCoefPrompt2_[1][2],errctCoefPrompt2_[1][2],ctCoefPrompt2_[2][2],errctCoefPrompt2_[2][2],ctCoefPrompt2_[3][2],errctCoefPrompt2_[3][2],ctCoefPrompt2_[4][2],errctCoefPrompt2_[4][2],ctCoefPrompt2_[5][2],errctCoefPrompt2_[5][2],ctCoefPrompt2_[6][2],errctCoefPrompt2_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ctCoefPrompt2_[1][3],errctCoefPrompt2_[1][3],ctCoefPrompt2_[2][3],errctCoefPrompt2_[2][3],ctCoefPrompt2_[3][3],errctCoefPrompt2_[3][3],ctCoefPrompt2_[4][3],errctCoefPrompt2_[4][3],ctCoefPrompt2_[5][3],errctCoefPrompt2_[5][3],ctCoefPrompt2_[6][3],errctCoefPrompt2_[6][3]);

	  fprintf(outputFile, "\n  ctCoefPrompt3   |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ctCoefPrompt3_[1][1],errctCoefPrompt3_[1][1],ctCoefPrompt3_[2][1],errctCoefPrompt3_[2][1],ctCoefPrompt3_[3][1],errctCoefPrompt3_[3][1],ctCoefPrompt3_[4][1],errctCoefPrompt3_[4][1],ctCoefPrompt3_[5][1],errctCoefPrompt3_[5][1],ctCoefPrompt3_[6][1],errctCoefPrompt3_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ctCoefPrompt3_[1][2],errctCoefPrompt3_[1][2],ctCoefPrompt3_[2][2],errctCoefPrompt3_[2][2],ctCoefPrompt3_[3][2],errctCoefPrompt3_[3][2],ctCoefPrompt3_[4][2],errctCoefPrompt3_[4][2],ctCoefPrompt3_[5][2],errctCoefPrompt3_[5][2],ctCoefPrompt3_[6][2],errctCoefPrompt3_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ctCoefPrompt3_[1][3],errctCoefPrompt3_[1][3],ctCoefPrompt3_[2][3],errctCoefPrompt3_[2][3],ctCoefPrompt3_[3][3],errctCoefPrompt3_[3][3],ctCoefPrompt3_[4][3],errctCoefPrompt3_[4][3],ctCoefPrompt3_[5][3],errctCoefPrompt3_[5][3],ctCoefPrompt3_[6][3],errctCoefPrompt3_[6][3]);
/*
	  fprintf(outputFile, "\n ctCoefPromptTot  |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ctCoefPromptTot_[1][1],errctCoefPromptTot_[1][1],ctCoefPromptTot_[2][1],errctCoefPromptTot_[2][1],ctCoefPromptTot_[3][1],errctCoefPromptTot_[3][1],ctCoefPromptTot_[4][1],errctCoefPromptTot_[4][1],ctCoefPromptTot_[5][1],errctCoefPromptTot_[5][1],ctCoefPromptTot_[6][1],errctCoefPromptTot_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ctCoefPromptTot_[1][2],errctCoefPromptTot_[1][2],ctCoefPromptTot_[2][2],errctCoefPromptTot_[2][2],ctCoefPromptTot_[3][2],errctCoefPromptTot_[3][2],ctCoefPromptTot_[4][2],errctCoefPromptTot_[4][2],ctCoefPromptTot_[5][2],errctCoefPromptTot_[5][2],ctCoefPromptTot_[6][2],errctCoefPromptTot_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ctCoefPromptTot_[1][3],errctCoefPromptTot_[1][3],ctCoefPromptTot_[2][3],errctCoefPromptTot_[2][3],ctCoefPromptTot_[3][3],errctCoefPromptTot_[3][3],ctCoefPromptTot_[4][3],errctCoefPromptTot_[4][3],ctCoefPromptTot_[5][3],errctCoefPromptTot_[5][3],ctCoefPromptTot_[6][3],errctCoefPromptTot_[6][3]);
*/
	  fprintf(outputFile, "\n  ctMeanPrompt1   |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ctMeanPrompt1_[1][1],errctMeanPrompt1_[1][1],ctMeanPrompt1_[2][1],errctMeanPrompt1_[2][1],ctMeanPrompt1_[3][1],errctMeanPrompt1_[3][1],ctMeanPrompt1_[4][1],errctMeanPrompt1_[4][1],ctMeanPrompt1_[5][1],errctMeanPrompt1_[5][1],ctMeanPrompt1_[6][1],errctMeanPrompt1_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ctMeanPrompt1_[1][2],errctMeanPrompt1_[1][2],ctMeanPrompt1_[2][2],errctMeanPrompt1_[2][2],ctMeanPrompt1_[3][2],errctMeanPrompt1_[3][2],ctMeanPrompt1_[4][2],errctMeanPrompt1_[4][2],ctMeanPrompt1_[5][2],errctMeanPrompt1_[5][2],ctMeanPrompt1_[6][2],errctMeanPrompt1_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ctMeanPrompt1_[1][3],errctMeanPrompt1_[1][3],ctMeanPrompt1_[2][3],errctMeanPrompt1_[2][3],ctMeanPrompt1_[3][3],errctMeanPrompt1_[3][3],ctMeanPrompt1_[4][3],errctMeanPrompt1_[4][3],ctMeanPrompt1_[5][3],errctMeanPrompt1_[5][3],ctMeanPrompt1_[6][3],errctMeanPrompt1_[6][3]);

	  fprintf(outputFile, "\n  ctMeanPrompt2   |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ctMeanPrompt2_[1][1],errctMeanPrompt2_[1][1],ctMeanPrompt2_[2][1],errctMeanPrompt2_[2][1],ctMeanPrompt2_[3][1],errctMeanPrompt2_[3][1],ctMeanPrompt2_[4][1],errctMeanPrompt2_[4][1],ctMeanPrompt2_[5][1],errctMeanPrompt2_[5][1],ctMeanPrompt2_[6][1],errctMeanPrompt2_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ctMeanPrompt2_[1][2],errctMeanPrompt2_[1][2],ctMeanPrompt2_[2][2],errctMeanPrompt2_[2][2],ctMeanPrompt2_[3][2],errctMeanPrompt2_[3][2],ctMeanPrompt2_[4][2],errctMeanPrompt2_[4][2],ctMeanPrompt2_[5][2],errctMeanPrompt2_[5][2],ctMeanPrompt2_[6][2],errctMeanPrompt2_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ctMeanPrompt2_[1][3],errctMeanPrompt2_[1][3],ctMeanPrompt2_[2][3],errctMeanPrompt2_[2][3],ctMeanPrompt2_[3][3],errctMeanPrompt2_[3][3],ctMeanPrompt2_[4][3],errctMeanPrompt2_[4][3],ctMeanPrompt2_[5][3],errctMeanPrompt2_[5][3],ctMeanPrompt2_[6][3],errctMeanPrompt2_[6][3]);


	  fprintf(outputFile, "\n  ctMeanPrompt3   |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ctMeanPrompt3_[1][1],errctMeanPrompt3_[1][1],ctMeanPrompt3_[2][1],errctMeanPrompt3_[2][1],ctMeanPrompt3_[3][1],errctMeanPrompt3_[3][1],ctMeanPrompt3_[4][1],errctMeanPrompt3_[4][1],ctMeanPrompt3_[5][1],errctMeanPrompt3_[5][1],ctMeanPrompt3_[6][1],errctMeanPrompt3_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ctMeanPrompt3_[1][2],errctMeanPrompt3_[1][2],ctMeanPrompt3_[2][2],errctMeanPrompt3_[2][2],ctMeanPrompt3_[3][2],errctMeanPrompt3_[3][2],ctMeanPrompt3_[4][2],errctMeanPrompt3_[4][2],ctMeanPrompt3_[5][2],errctMeanPrompt3_[5][2],ctMeanPrompt3_[6][2],errctMeanPrompt3_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ctMeanPrompt3_[1][3],errctMeanPrompt3_[1][3],ctMeanPrompt3_[2][3],errctMeanPrompt3_[2][3],ctMeanPrompt3_[3][3],errctMeanPrompt3_[3][3],ctMeanPrompt3_[4][3],errctMeanPrompt3_[4][3],ctMeanPrompt3_[5][3],errctMeanPrompt3_[5][3],ctMeanPrompt3_[6][3],errctMeanPrompt3_[6][3]);

	  fprintf(outputFile, "\n  ctSigmaPrompt1  |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ctSigmaPrompt1_[1][1],errctSigmaPrompt1_[1][1],ctSigmaPrompt1_[2][1],errctSigmaPrompt1_[2][1],ctSigmaPrompt1_[3][1],errctSigmaPrompt1_[3][1],ctSigmaPrompt1_[4][1],errctSigmaPrompt1_[4][1],ctSigmaPrompt1_[5][1],errctSigmaPrompt1_[5][1],ctSigmaPrompt1_[6][1],errctSigmaPrompt1_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ctSigmaPrompt1_[1][2],errctSigmaPrompt1_[1][2],ctSigmaPrompt1_[2][2],errctSigmaPrompt1_[2][2],ctSigmaPrompt1_[3][2],errctSigmaPrompt1_[3][2],ctSigmaPrompt1_[4][2],errctSigmaPrompt1_[4][2],ctSigmaPrompt1_[5][2],errctSigmaPrompt1_[5][2],ctSigmaPrompt1_[6][2],errctSigmaPrompt1_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ctSigmaPrompt1_[1][3],errctSigmaPrompt1_[1][3],ctSigmaPrompt1_[2][3],errctSigmaPrompt1_[2][3],ctSigmaPrompt1_[3][3],errctSigmaPrompt1_[3][3],ctSigmaPrompt1_[4][3],errctSigmaPrompt1_[4][3],ctSigmaPrompt1_[5][3],errctSigmaPrompt1_[5][3],ctSigmaPrompt1_[6][3],errctSigmaPrompt1_[6][3]);

	  fprintf(outputFile, "\n  ctSigmaPrompt2  |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ctSigmaPrompt2_[1][1],errctSigmaPrompt2_[1][1],ctSigmaPrompt2_[2][1],errctSigmaPrompt2_[2][1],ctSigmaPrompt2_[3][1],errctSigmaPrompt2_[3][1],ctSigmaPrompt2_[4][1],errctSigmaPrompt2_[4][1],ctSigmaPrompt2_[5][1],errctSigmaPrompt2_[5][1],ctSigmaPrompt2_[6][1],errctSigmaPrompt2_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ctSigmaPrompt2_[1][2],errctSigmaPrompt2_[1][2],ctSigmaPrompt2_[2][2],errctSigmaPrompt2_[2][2],ctSigmaPrompt2_[3][2],errctSigmaPrompt2_[3][2],ctSigmaPrompt2_[4][2],errctSigmaPrompt2_[4][2],ctSigmaPrompt2_[5][2],errctSigmaPrompt2_[5][2],ctSigmaPrompt2_[6][2],errctSigmaPrompt2_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ctSigmaPrompt2_[1][3],errctSigmaPrompt2_[1][3],ctSigmaPrompt2_[2][3],errctSigmaPrompt2_[2][3],ctSigmaPrompt2_[3][3],errctSigmaPrompt2_[3][3],ctSigmaPrompt2_[4][3],errctSigmaPrompt2_[4][3],ctSigmaPrompt2_[5][3],errctSigmaPrompt2_[5][3],ctSigmaPrompt2_[6][3],errctSigmaPrompt2_[6][3]);

	  fprintf(outputFile, "\n  ctSigmaPrompt3  |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ctSigmaPrompt3_[1][1],errctSigmaPrompt3_[1][1],ctSigmaPrompt3_[2][1],errctSigmaPrompt3_[2][1],ctSigmaPrompt3_[3][1],errctSigmaPrompt3_[3][1],ctSigmaPrompt3_[4][1],errctSigmaPrompt3_[4][1],ctSigmaPrompt3_[5][1],errctSigmaPrompt3_[5][1],ctSigmaPrompt3_[6][1],errctSigmaPrompt3_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ctSigmaPrompt3_[1][2],errctSigmaPrompt3_[1][2],ctSigmaPrompt3_[2][2],errctSigmaPrompt3_[2][2],ctSigmaPrompt3_[3][2],errctSigmaPrompt3_[3][2],ctSigmaPrompt3_[4][2],errctSigmaPrompt3_[4][2],ctSigmaPrompt3_[5][2],errctSigmaPrompt3_[5][2],ctSigmaPrompt3_[6][2],errctSigmaPrompt3_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ctSigmaPrompt3_[1][3],errctSigmaPrompt3_[1][3],ctSigmaPrompt3_[2][3],errctSigmaPrompt3_[2][3],ctSigmaPrompt3_[3][3],errctSigmaPrompt3_[3][3],ctSigmaPrompt3_[4][3],errctSigmaPrompt3_[4][3],ctSigmaPrompt3_[5][3],errctSigmaPrompt3_[5][3],ctSigmaPrompt3_[6][3],errctSigmaPrompt3_[6][3]);

	  fprintf(outputFile, "\nNon Prompt Component (MC):\n");

	  fprintf(outputFile, "\n  ssTauNonPrompt  |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ssTauNonPrompt_[1][1],errssTauNonPrompt_[1][1],ssTauNonPrompt_[2][1],errssTauNonPrompt_[2][1],ssTauNonPrompt_[3][1],errssTauNonPrompt_[3][1],ssTauNonPrompt_[4][1],errssTauNonPrompt_[4][1],ssTauNonPrompt_[5][1],errssTauNonPrompt_[5][1],ssTauNonPrompt_[6][1],errssTauNonPrompt_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ssTauNonPrompt_[1][2],errssTauNonPrompt_[1][2],ssTauNonPrompt_[2][2],errssTauNonPrompt_[2][2],ssTauNonPrompt_[3][2],errssTauNonPrompt_[3][2],ssTauNonPrompt_[4][2],errssTauNonPrompt_[4][2],ssTauNonPrompt_[5][2],errssTauNonPrompt_[5][2],ssTauNonPrompt_[6][2],errssTauNonPrompt_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ssTauNonPrompt_[1][3],errssTauNonPrompt_[1][3],ssTauNonPrompt_[2][3],errssTauNonPrompt_[2][3],ssTauNonPrompt_[3][3],errssTauNonPrompt_[3][3],ssTauNonPrompt_[4][3],errssTauNonPrompt_[4][3],ssTauNonPrompt_[5][3],errssTauNonPrompt_[5][3],ssTauNonPrompt_[6][3],errssTauNonPrompt_[6][3]);

	  fprintf(outputFile, "\n  fTauNonPrompt   |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],fTauNonPrompt_[1][1],errfTauNonPrompt_[1][1],fTauNonPrompt_[2][1],errfTauNonPrompt_[2][1],fTauNonPrompt_[3][1],errfTauNonPrompt_[3][1],fTauNonPrompt_[4][1],errfTauNonPrompt_[4][1],fTauNonPrompt_[5][1],errfTauNonPrompt_[5][1],fTauNonPrompt_[6][1],errfTauNonPrompt_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],fTauNonPrompt_[1][2],errfTauNonPrompt_[1][2],fTauNonPrompt_[2][2],errfTauNonPrompt_[2][2],fTauNonPrompt_[3][2],errfTauNonPrompt_[3][2],fTauNonPrompt_[4][2],errfTauNonPrompt_[4][2],fTauNonPrompt_[5][2],errfTauNonPrompt_[5][2],fTauNonPrompt_[6][2],errfTauNonPrompt_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],fTauNonPrompt_[1][3],errfTauNonPrompt_[1][3],fTauNonPrompt_[2][3],errfTauNonPrompt_[2][3],fTauNonPrompt_[3][3],errfTauNonPrompt_[3][3],fTauNonPrompt_[4][3],errfTauNonPrompt_[4][3],fTauNonPrompt_[5][3],errfTauNonPrompt_[5][3],fTauNonPrompt_[6][3],errfTauNonPrompt_[6][3]);
/*
	  fprintf(outputFile, "\n  dsTauNonPrompt  |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],dsTauNonPrompt_[1][1],errdsTauNonPrompt_[1][1],dsTauNonPrompt_[2][1],errdsTauNonPrompt_[2][1],dsTauNonPrompt_[3][1],errdsTauNonPrompt_[3][1],dsTauNonPrompt_[4][1],errdsTauNonPrompt_[4][1],dsTauNonPrompt_[5][1],errdsTauNonPrompt_[5][1],dsTauNonPrompt_[6][1],errdsTauNonPrompt_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],dsTauNonPrompt_[1][2],errdsTauNonPrompt_[1][2],dsTauNonPrompt_[2][2],errdsTauNonPrompt_[2][2],dsTauNonPrompt_[3][2],errdsTauNonPrompt_[3][2],dsTauNonPrompt_[4][2],errdsTauNonPrompt_[4][2],dsTauNonPrompt_[5][2],errdsTauNonPrompt_[5][2],dsTauNonPrompt_[6][2],errdsTauNonPrompt_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],dsTauNonPrompt_[1][3],errdsTauNonPrompt_[1][3],dsTauNonPrompt_[2][3],errdsTauNonPrompt_[2][3],dsTauNonPrompt_[3][3],errdsTauNonPrompt_[3][3],dsTauNonPrompt_[4][3],errdsTauNonPrompt_[4][3],dsTauNonPrompt_[5][3],errdsTauNonPrompt_[5][3],dsTauNonPrompt_[6][3],errdsTauNonPrompt_[6][3]);
*/
	  fprintf(outputFile, "\n ssCoefNonPrompt  |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ssCoefNonPrompt_[1][1],errssCoefNonPrompt_[1][1],ssCoefNonPrompt_[2][1],errssCoefNonPrompt_[2][1],ssCoefNonPrompt_[3][1],errssCoefNonPrompt_[3][1],ssCoefNonPrompt_[4][1],errssCoefNonPrompt_[4][1],ssCoefNonPrompt_[5][1],errssCoefNonPrompt_[5][1],ssCoefNonPrompt_[6][1],errssCoefNonPrompt_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ssCoefNonPrompt_[1][2],errssCoefNonPrompt_[1][2],ssCoefNonPrompt_[2][2],errssCoefNonPrompt_[2][2],ssCoefNonPrompt_[3][2],errssCoefNonPrompt_[3][2],ssCoefNonPrompt_[4][2],errssCoefNonPrompt_[4][2],ssCoefNonPrompt_[5][2],errssCoefNonPrompt_[5][2],ssCoefNonPrompt_[6][2],errssCoefNonPrompt_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ssCoefNonPrompt_[1][3],errssCoefNonPrompt_[1][3],ssCoefNonPrompt_[2][3],errssCoefNonPrompt_[2][3],ssCoefNonPrompt_[3][3],errssCoefNonPrompt_[3][3],ssCoefNonPrompt_[4][3],errssCoefNonPrompt_[4][3],ssCoefNonPrompt_[5][3],errssCoefNonPrompt_[5][3],ssCoefNonPrompt_[6][3],errssCoefNonPrompt_[6][3]);

	  fprintf(outputFile, "\n fCoefNonPrompt  |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],fCoefNonPrompt_[1][1],errfCoefNonPrompt_[1][1],fCoefNonPrompt_[2][1],errfCoefNonPrompt_[2][1],fCoefNonPrompt_[3][1],errfCoefNonPrompt_[3][1],fCoefNonPrompt_[4][1],errfCoefNonPrompt_[4][1],fCoefNonPrompt_[5][1],errfCoefNonPrompt_[5][1],fCoefNonPrompt_[6][1],errfCoefNonPrompt_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],fCoefNonPrompt_[1][2],errfCoefNonPrompt_[1][2],fCoefNonPrompt_[2][2],errfCoefNonPrompt_[2][2],fCoefNonPrompt_[3][2],errfCoefNonPrompt_[3][2],fCoefNonPrompt_[4][2],errfCoefNonPrompt_[4][2],fCoefNonPrompt_[5][2],errfCoefNonPrompt_[5][2],fCoefNonPrompt_[6][2],errfCoefNonPrompt_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],fCoefNonPrompt_[1][3],errfCoefNonPrompt_[1][3],fCoefNonPrompt_[2][3],errfCoefNonPrompt_[2][3],fCoefNonPrompt_[3][3],errfCoefNonPrompt_[3][3],fCoefNonPrompt_[4][3],errfCoefNonPrompt_[4][3],fCoefNonPrompt_[5][3],errfCoefNonPrompt_[5][3],fCoefNonPrompt_[6][3],errfCoefNonPrompt_[6][3]);

	  fprintf(outputFile, "\nBackground Component (Data Sidebands):\n");

	  fprintf(outputFile, "\n  ssTauBkg        |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ssTauBkg_[1][1],errssTauBkg_[1][1],ssTauBkg_[2][1],errssTauBkg_[2][1],ssTauBkg_[3][1],errssTauBkg_[3][1],ssTauBkg_[4][1],errssTauBkg_[4][1],ssTauBkg_[5][1],errssTauBkg_[5][1],ssTauBkg_[6][1],errssTauBkg_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ssTauBkg_[1][2],errssTauBkg_[1][2],ssTauBkg_[2][2],errssTauBkg_[2][2],ssTauBkg_[3][2],errssTauBkg_[3][2],ssTauBkg_[4][2],errssTauBkg_[4][2],ssTauBkg_[5][2],errssTauBkg_[5][2],ssTauBkg_[6][2],errssTauBkg_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ssTauBkg_[1][3],errssTauBkg_[1][3],ssTauBkg_[2][3],errssTauBkg_[2][3],ssTauBkg_[3][3],errssTauBkg_[3][3],ssTauBkg_[4][3],errssTauBkg_[4][3],ssTauBkg_[5][3],errssTauBkg_[5][3],ssTauBkg_[6][3],errssTauBkg_[6][3]);

	  fprintf(outputFile, "\n  fTauBkg         |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],fTauBkg_[1][1],errfTauBkg_[1][1],fTauBkg_[2][1],errfTauBkg_[2][1],fTauBkg_[3][1],errfTauBkg_[3][1],fTauBkg_[4][1],errfTauBkg_[4][1],fTauBkg_[5][1],errfTauBkg_[5][1],fTauBkg_[6][1],errfTauBkg_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],fTauBkg_[1][2],errfTauBkg_[1][2],fTauBkg_[2][2],errfTauBkg_[2][2],fTauBkg_[3][2],errfTauBkg_[3][2],fTauBkg_[4][2],errfTauBkg_[4][2],fTauBkg_[5][2],errfTauBkg_[5][2],fTauBkg_[6][2],errfTauBkg_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],fTauBkg_[1][3],errfTauBkg_[1][3],fTauBkg_[2][3],errfTauBkg_[2][3],fTauBkg_[3][3],errfTauBkg_[3][3],fTauBkg_[4][3],errfTauBkg_[4][3],fTauBkg_[5][3],errfTauBkg_[5][3],fTauBkg_[6][3],errfTauBkg_[6][3]);

	  fprintf(outputFile, "\n  dsTauBkg        |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],dsTauBkg_[1][1],errdsTauBkg_[1][1],dsTauBkg_[2][1],errdsTauBkg_[2][1],dsTauBkg_[3][1],errdsTauBkg_[3][1],dsTauBkg_[4][1],errdsTauBkg_[4][1],dsTauBkg_[5][1],errdsTauBkg_[5][1],dsTauBkg_[6][1],errdsTauBkg_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],dsTauBkg_[1][2],errdsTauBkg_[1][2],dsTauBkg_[2][2],errdsTauBkg_[2][2],dsTauBkg_[3][2],errdsTauBkg_[3][2],dsTauBkg_[4][2],errdsTauBkg_[4][2],dsTauBkg_[5][2],errdsTauBkg_[5][2],dsTauBkg_[6][2],errdsTauBkg_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],dsTauBkg_[1][3],errdsTauBkg_[1][3],dsTauBkg_[2][3],errdsTauBkg_[2][3],dsTauBkg_[3][3],errdsTauBkg_[3][3],dsTauBkg_[4][3],errdsTauBkg_[4][3],dsTauBkg_[5][3],errdsTauBkg_[5][3],dsTauBkg_[6][3],errdsTauBkg_[6][3]);

	  fprintf(outputFile, "\n  ssCoefBkg       |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],ssCoefBkg_[1][1],errssCoefBkg_[1][1],ssCoefBkg_[2][1],errssCoefBkg_[2][1],ssCoefBkg_[3][1],errssCoefBkg_[3][1],ssCoefBkg_[4][1],errssCoefBkg_[4][1],ssCoefBkg_[5][1],errssCoefBkg_[5][1],ssCoefBkg_[6][1],errssCoefBkg_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],ssCoefBkg_[1][2],errssCoefBkg_[1][2],ssCoefBkg_[2][2],errssCoefBkg_[2][2],ssCoefBkg_[3][2],errssCoefBkg_[3][2],ssCoefBkg_[4][2],errssCoefBkg_[4][2],ssCoefBkg_[5][2],errssCoefBkg_[5][2],ssCoefBkg_[6][2],errssCoefBkg_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],ssCoefBkg_[1][3],errssCoefBkg_[1][3],ssCoefBkg_[2][3],errssCoefBkg_[2][3],ssCoefBkg_[3][3],errssCoefBkg_[3][3],ssCoefBkg_[4][3],errssCoefBkg_[4][3],ssCoefBkg_[5][3],errssCoefBkg_[5][3],ssCoefBkg_[6][3],errssCoefBkg_[6][3]);

	  fprintf(outputFile, "\n  fCoefBkg        |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],fCoefBkg_[1][1],errfCoefBkg_[1][1],fCoefBkg_[2][1],errfCoefBkg_[2][1],fCoefBkg_[3][1],errfCoefBkg_[3][1],fCoefBkg_[4][1],errfCoefBkg_[4][1],fCoefBkg_[5][1],errfCoefBkg_[5][1],fCoefBkg_[6][1],errfCoefBkg_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],fCoefBkg_[1][2],errfCoefBkg_[1][2],fCoefBkg_[2][2],errfCoefBkg_[2][2],fCoefBkg_[3][2],errfCoefBkg_[3][2],fCoefBkg_[4][2],errfCoefBkg_[4][2],fCoefBkg_[5][2],errfCoefBkg_[5][2],fCoefBkg_[6][2],errfCoefBkg_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],fCoefBkg_[1][3],errfCoefBkg_[1][3],fCoefBkg_[2][3],errfCoefBkg_[2][3],fCoefBkg_[3][3],errfCoefBkg_[3][3],fCoefBkg_[4][3],errfCoefBkg_[4][3],fCoefBkg_[5][3],errfCoefBkg_[5][3],fCoefBkg_[6][3],errfCoefBkg_[6][3]);

	  fprintf(outputFile, "\n  dsCoefBkg        |   %1.1f < pT < %1.1f   |    %1.1f < pT < %1.1f  |    %1.1f < pT < %1.1f  |   %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |  %1.1f < pT < %1.1f  |\n",jpsi::pTRange[0],jpsi::pTRange[1],jpsi::pTRange[1],jpsi::pTRange[2],jpsi::pTRange[2],jpsi::pTRange[3],jpsi::pTRange[3],jpsi::pTRange[4],jpsi::pTRange[4],jpsi::pTRange[5],jpsi::pTRange[5],jpsi::pTRange[6]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[0],jpsi::rapForPTRange[1],dsCoefBkg_[1][1],errdsCoefBkg_[1][1],dsCoefBkg_[2][1],errdsCoefBkg_[2][1],dsCoefBkg_[3][1],errdsCoefBkg_[3][1],dsCoefBkg_[4][1],errdsCoefBkg_[4][1],dsCoefBkg_[5][1],errdsCoefBkg_[5][1],dsCoefBkg_[6][1],errdsCoefBkg_[6][1]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[1],jpsi::rapForPTRange[2],dsCoefBkg_[1][2],errdsCoefBkg_[1][2],dsCoefBkg_[2][2],errdsCoefBkg_[2][2],dsCoefBkg_[3][2],errdsCoefBkg_[3][2],dsCoefBkg_[4][2],errdsCoefBkg_[4][2],dsCoefBkg_[5][2],errdsCoefBkg_[5][2],dsCoefBkg_[6][2],errdsCoefBkg_[6][2]);
	  fprintf(outputFile, "%1.1f < |rap| < %1.1f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f | %1.5f +- %1.5f |\n", jpsi::rapForPTRange[2],jpsi::rapForPTRange[3],dsCoefBkg_[1][3],errdsCoefBkg_[1][3],dsCoefBkg_[2][3],errdsCoefBkg_[2][3],dsCoefBkg_[3][3],errdsCoefBkg_[3][3],dsCoefBkg_[4][3],errdsCoefBkg_[4][3],dsCoefBkg_[5][3],errdsCoefBkg_[5][3],dsCoefBkg_[6][3],errdsCoefBkg_[6][3]);


	fclose(outputFile);

	  cout<<"there"<<endl;


  delete dataTreesNP;
  delete dataNP;
  delete dataTreesPR;
  delete dataPR;
//  delete dataTreesBK;
//  delete dataBK;
  delete fInput;



  return 0;
}
