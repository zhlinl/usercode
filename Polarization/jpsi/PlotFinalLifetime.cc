#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
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
  
  bool pereverr(false),newdata(false);

  for(int i=0;i<argc; ++i) {
    if(std::string(argv[i]).find("--perEventErrors") != std::string::npos) pereverr = true;
    if(std::string(argv[i]).find("--newTTree") != std::string::npos) newdata = true;
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
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,1.5);//0=GG,1=GT,2=TT
  RooRealVar MCweight("MCweight","MCweight",0,1.1);
  RooRealVar HLT_Mu0_TkMu0_Jpsi("HLT_Mu0_TkMu0_Jpsi","Passes HLT_Mu0_TkMu0_Jpsi Trigger",.5,1.5);
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  if(pereverr)
    varlist.add(JpsictErr);
  varlist.add(MCType_idx);
  if(newdata)
    varlist.add(HLT_Mu0_TkMu0_Jpsi);
  varlist.add(JpsiType_idx);

  TChain *dataTrees = new TChain("data");
  dataTrees->Add("/Users/lindseygray/jpsi_workdir/data/TTree_pol_Mu0TkMu0Jpsi_dataR_30Aug2010.root");
  
  Char_t *fileNameInput = "jPsiFit.root";
  TFile* fInput = new TFile(fileNameInput);
  
  Char_t *fileNameInputFinal = "jPsiFitFinal_cs.root";
  TFile* fInputFinal = new TFile(fileNameInputFinal);
  
  RooDataSet *data = new RooDataSet("data","Supplied Data",dataTrees,varlist);

  RooDataSet *databin;

  LifetimeModel* lifetinp_;

  RooPlot* LifetimeFrame;


 for(int ptBin = 0; ptBin < jpsi::kNbPTBins; ++ptBin) {
    for(int yBin = 0; yBin < jpsi::kNbRapForPTBins; ++yBin) {

      JpsiMass.setRange("lowBand",2.7,
			jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgLow);
      JpsiMass.setRange("signalRegion",
			jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass,
			jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass);
      JpsiMass.setRange("highBand",
			jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgHigh,3.5);
      
      cout<<"pT"<<ptBin+1<<"rapidity"<<yBin+1<<endl;

	  char reduce[200];
	  sprintf(reduce,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f && JpsiMass < %f | JpsiMass > %f",
		  jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],
		  jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgLow,
		  jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgHigh);
	  

	  databin = (RooDataSet*)data->reduce(reduce);

	  char DirectoryPath[200];
	  sprintf(DirectoryPath,"/pt%d_rapidity%d",ptBin+1,yBin+1);

	  TDirectory *InputDirectory = (TDirectory*)fInputFinal->GetDirectory(DirectoryPath);

	  CompositeModelBuilder *LifetimeModel = new CompositeModelBuilder("","");
	  LifetimeModel->setUseLifetime(true);
	  LifetimeModel->setUseMass(false);
	  LifetimeModel->setUsePol(false);

	  LifetimeModel->loadParameters(*InputDirectory);

	  if(pereverr)
	    LifetimeModel->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS);
	  else
	    LifetimeModel->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);

//	  lifetinp_ = new LifetimeModel();
//	  lifetinp_->loadParameters(*InputDirectory);
//	  lifetinp_->initModels(Jpsict);

	  double N_PR = 2000;
	  double N_NP = 500;
	  double N_BK = 100;
	  double N_T = N_PR+N_NP+N_BK;

	  LifetimeFrame = new RooPlot;
	  LifetimeFrame = Jpsict.frame(Bins(50)) ;
	  LifetimeFrame->setPadFactor(0.5);
	  databin->plotOn(LifetimeFrame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	  if(pereverr)
	    LifetimeModel->model()->plotOn(LifetimeFrame, RooFit::LineWidth(3),
					   RooFit::LineColor(kBlue),RooFit::Normalization(1.0),RooFit::Range("signalRegion"),
					   RooFit::ProjWData(RooArgSet(JpsictErr),*databin));
	  else
	     LifetimeModel->model()->plotOn(LifetimeFrame, RooFit::LineWidth(3),
					    RooFit::LineColor(kBlue),RooFit::Normalization(1.0),RooFit::Range("signalRegion"));

	  std::cout << "Total expected events in Signal Region: " << LifetimeModel->model()->expectedEvents(RooArgSet(Jpsict));
	  double chi2_LifetimeFrame=LifetimeFrame->chiSquare();
	  if(pereverr)
	    LifetimeModel->model()->plotOn(LifetimeFrame, Components(*LifetimeModel->getPromptModel()), RooFit::LineStyle(kSolid),
					   RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::Normalization(1.0),RooFit::Range("signalRegion"),
					   RooFit::ProjWData(RooArgSet(JpsictErr),*databin));
	  else
	    LifetimeModel->model()->plotOn(LifetimeFrame, Components(*LifetimeModel->getPromptModel()), RooFit::LineStyle(kSolid),
					   RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::Normalization(1.0),RooFit::Range("signalRegion"));
	  std::cout << "Expected Prompt Events: " << LifetimeModel->getPromptModel()->expectedEvents(RooArgSet(Jpsict));
	  if(pereverr)
	    LifetimeModel->model()->plotOn(LifetimeFrame, Components(*LifetimeModel->getNonPromptModel()), RooFit::LineStyle(kSolid),RooFit::LineColor(kRed),
					   RooFit::LineWidth(2),RooFit::Normalization(1.0),RooFit::Range("signalRegion"),
					   RooFit::ProjWData(RooArgSet(JpsictErr),*databin));
	  else
	    LifetimeModel->model()->plotOn(LifetimeFrame, Components(*LifetimeModel->getNonPromptModel()), RooFit::LineStyle(kSolid),RooFit::LineColor(kRed),
					   RooFit::LineWidth(2),RooFit::Normalization(1.0),RooFit::Range("signalRegion"));

	  std::cout << "Expected Non-Prompt Events: " << LifetimeModel->getNonPromptModel()->expectedEvents(RooArgSet(Jpsict));
	  if(pereverr)
	    LifetimeModel->model()->plotOn(LifetimeFrame, Components(*LifetimeModel->getBkgModel()), RooFit::LineStyle(kSolid),RooFit::LineColor(kGreen),
					   RooFit::LineWidth(2),RooFit::Normalization(1.0),RooFit::Range("signalRegion"),
					   RooFit::ProjWData(RooArgSet(JpsictErr),*databin));
	  else
	    LifetimeModel->model()->plotOn(LifetimeFrame, Components(*LifetimeModel->getNonPromptModel()), RooFit::LineStyle(kSolid),RooFit::LineColor(kRed),
					   RooFit::LineWidth(2),RooFit::Normalization(1.0),RooFit::Range("signalRegion"));

	  std::cout << "Expected Background Events: " << LifetimeModel->getBkgModel()->expectedEvents(RooArgSet(Jpsict));
	  LifetimeFrame->SetMinimum(.5);
	  char LifetimeTitle[200];
	  sprintf(LifetimeTitle,"Data Lifetime Fit %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_LifetimeFrame);
	  LifetimeFrame->SetTitle(LifetimeTitle);

  TH1* legendBlue = data->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(3) ;
  TH1* legendRed = data->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(kSolid) ; legendRed->SetLineWidth(2) ;
  TH1* legendBlack = data->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
  TH1* legendGreen = data->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;

  TLegend* LifetimeLegend=new TLegend(0.6,0.6,0.9,0.9);
  LifetimeLegend->SetFillColor(kWhite);
  LifetimeLegend->SetTextFont(72);
  LifetimeLegend->SetTextSize(0.025);
  LifetimeLegend->AddEntry(legendBlue,"Lifetime Fit","l");
  LifetimeLegend->AddEntry(legendBlack,"Prompt Component","l");
  LifetimeLegend->AddEntry(legendRed,"Non Prompt Component","l");
  LifetimeLegend->AddEntry(legendGreen,"Background Component","l");

  TCanvas* LifetimeCanvas = new TCanvas("LifetimeCanvas","LifetimeCanvas",700, 600);
  LifetimeCanvas->Divide(1);  LifetimeCanvas->SetFillColor(kWhite);
  LifetimeCanvas->cd(1)->SetLogy(1);
  LifetimeCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); LifetimeFrame->GetYaxis()->SetTitleOffset(2) ; LifetimeFrame->Draw(); LifetimeLegend->Draw();

  char FilenameLifetime[200];
  sprintf(FilenameLifetime,"Plots/Lifetime/FinalLifetime_pt%d_rapidity%d.png",ptBin+1,yBin+1);

  LifetimeCanvas->SaveAs(FilenameLifetime);

  LifetimeCanvas->Close();


  delete LifetimeFrame;

  delete LifetimeModel;
  delete LifetimeLegend;

	  }
	}


  delete dataTrees;
  delete data;

  delete fInput;
  delete fInputFinal;



  return 0;
}
