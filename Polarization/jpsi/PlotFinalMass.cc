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
#include "RooChebychev.h"
#include "RooAbsPdf.h"
#include "RooGaussModel.h"
#include "RooCBShape.h"
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

  bool newdata(false);
  
  for(int i=0;i<argc; ++i) {
    if(std::string(argv[i]).find("--newTTree") != std::string::npos) newdata = true;
  }

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
  RooRealVar HLT_Mu0_TkMu0_Jpsi("HLT_Mu0_TkMu0_Jpsi","Passes HLT_Mu0_TkMu0_Jpsi Trigger",.5,1.5);
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);
  if(newdata)
    varlist.add(HLT_Mu0_TkMu0_Jpsi);

  TChain *dataTrees = new TChain("data");
  dataTrees->Add("/Users/lindseygray/jpsi_workdir/data/TTree_pol_Mu0TkMu0Jpsi_dataR_30Aug2010.root");

   Char_t *fileNameInput = "jPsiFit.root";
   TFile* fInput = new TFile(fileNameInput);

   Char_t *fileNameInputFinal = "jPsiFitFinal_cs.root";
      TFile* fInputFinal = new TFile(fileNameInputFinal);

  RooDataSet *data = new RooDataSet("data","Supplied Data",dataTrees,varlist);

  RooDataSet *databin;


  RooPlot* MassFrame;


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
      sprintf(reduce,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",
	      jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
      
      databin = (RooDataSet*)data->reduce(reduce);
      
      char DirectoryPath[200];
      sprintf(DirectoryPath,"/pt%d_rapidity%d",ptBin+1,yBin+1);
      
      TDirectory *InputDirectory = (TDirectory*)fInputFinal->GetDirectory(DirectoryPath);
      
      CompositeModelBuilder *MassModel = new CompositeModelBuilder("","");
      MassModel->setUseLifetime(false);
      MassModel->setUseMass(true);
      MassModel->setUsePol(false);
      
      MassModel->loadParameters(*InputDirectory);
      
      MassModel->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);
      
      
      MassFrame = new RooPlot;
      MassFrame = JpsiMass.frame(Bins(50)) ;
      MassFrame->setPadFactor(0.5);
      databin->plotOn(MassFrame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
      MassModel->model()->plotOn(MassFrame, RooFit::LineWidth(3),RooFit::LineColor(kBlue),
				 RooFit::Normalization(1.0),RooFit::Range("lowBand,signalRegion,highBand"));
      double chi2_MassFrame=MassFrame->chiSquare();
      MassModel->model()->plotOn(MassFrame, Components(*MassModel->getPromptModel()), RooFit::LineStyle(kSolid),
				 RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::Normalization(1.0),RooFit::Range("lowBand,signalRegion,highBand"));
      MassModel->model()->plotOn(MassFrame, Components(*MassModel->getNonPromptModel()), RooFit::LineStyle(kSolid),
				 RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::Normalization(1.0),RooFit::Range("lowBand,signalRegion,highBand"));
      MassModel->model()->plotOn(MassFrame, Components(*MassModel->getBkgModel()), RooFit::LineStyle(kSolid),
				 RooFit::LineColor(kGreen),RooFit::LineWidth(2),RooFit::Normalization(1.0),RooFit::Range("lowBand,signalRegion,highBand"));
      MassFrame->SetMinimum(0);
      char MassTitle[200];
      sprintf(MassTitle,"Data Mass Fit %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",
	      jpsi::pTRange[ptBin],jpsi::pTRange[ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],
	      chi2_MassFrame);
      MassFrame->SetTitle(MassTitle);
      
      TH1* legendBlue = data->createHistogram("legendBlue",JpsiMass,Binning(100)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(3) ;
      TH1* legendRed = data->createHistogram("legendRed",JpsiMass,Binning(100)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(kSolid) ; legendRed->SetLineWidth(2) ;
      TH1* legendBlack = data->createHistogram("legendBlack",JpsiMass,Binning(100)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
      TH1* legendGreen = data->createHistogram("legendGreen",JpsiMass,Binning(100)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
      
      TLegend* MassLegend=new TLegend(0.6,0.6,0.9,0.9);
      MassLegend->SetFillColor(kWhite);
      MassLegend->SetTextFont(72);
      MassLegend->SetTextSize(0.025);
      MassLegend->AddEntry(legendBlue,"Mass Fit","l");
      MassLegend->AddEntry(legendBlack,"Mass Fit (prompt)","l");
      MassLegend->AddEntry(legendRed,"Mass Fit (non-prompt)","l");
      MassLegend->AddEntry(legendGreen,"Background Component","l");
      
      TCanvas* MassCanvas = new TCanvas("LifetimeCanvas","LifetimeCanvas",700, 600);
      MassCanvas->Divide(1);  MassCanvas->SetFillColor(kWhite);  
      MassCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); MassFrame->GetYaxis()->SetTitleOffset(2) ; MassFrame->Draw(); MassLegend->Draw();
      
      char FilenameMass[200];
      sprintf(FilenameMass,"Plots/Mass/FinalMass_pt%d_rapidity%d.png",ptBin+1,yBin+1);
      
      MassCanvas->SaveAs(FilenameMass);
      
      MassCanvas->Close();
      
      
      delete MassFrame;
      
      delete MassModel;
      delete MassLegend;
      
    }
 }
 
 
 delete dataTrees;
 delete data;
 
 delete fInput;
 delete fInputFinal;
 
 
 
 return 0;
}
