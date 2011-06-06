#include <iostream>
//#include <stdio.h>
//#include <direct.h>
//#include <sys/stat.h>
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
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TStyle.h"



int main(int argc, char**argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;


   gStyle->SetPalette(1,0);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetPadRightMargin(0.05);
//  gStyle->SetPadTopMargin(0.15);

   gStyle->SetTickLength(-0.02, "xyz");
   gStyle->SetLabelOffset(0.02, "x");
   gStyle->SetLabelOffset(0.02, "y");
   gStyle->SetTitleOffset(1.3, "x");
   gStyle->SetTitleOffset(1.4, "y");

  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
  RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,2.5);//0=GG,1=GT,2=TT
  RooRealVar MCweight("MCweight","MCweight",0,500);
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);

//  const char Direc[200];
//  sprintf(Direc,"Directory");
//  int stat = mkdir("Direc");

   Char_t *fileNameInputCS = "jPsiFitFinal_cs.root";
   TFile* fInputCS = new TFile(fileNameInputCS);

   Char_t *fileNameInputHX = "jPsiFitFinal_hx.root";
   TFile* fInputHX = new TFile(fileNameInputHX);


	RooRealVar *nNonPromptCS, *nPromptCS, *nBackground, *nNonPromptHX, *nPromptHX;
	RooRealVar *lambda_theta_CS_PR, *lambda_phi_CS_PR, *lambda_theta_HX_PR, *lambda_phi_HX_PR;

	double nNonPromptCS_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errnNonPromptCS_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double nPromptCS_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errnPromptCS_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double nBackground_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errnBackground_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double BfractionCS_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errBfractionCS_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double nNonPromptHX_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errnNonPromptHX_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double nPromptHX_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errnPromptHX_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double BfractionHX_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errBfractionHX_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double sigOvbkg_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errsigOvbkg_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double lambda_theta_CS_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errlambda_theta_CS_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double lambda_phi_CS_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errlambda_phi_CS_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double lambda_theta_HX_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errlambda_theta_HX_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double lambda_phi_HX_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errlambda_phi_HX_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double invariant_lambda_CS_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errinvariant_lambda_CS_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double invariant_lambda_HX_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errinvariant_lambda_HX_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double invariant_F_CS_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errinvariant_F_CS_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double invariant_F_HX_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];
	double errinvariant_F_HX_PR_[jpsi::kNbPTBins][jpsi::kNbRapForPTBins];


	double ptMean[jpsi::kNbPTBins];
	double errptMean[jpsi::kNbPTBins];
	double BfractionCS_rap[jpsi::kNbPTBins];
	double errBfractionCS_rap[jpsi::kNbPTBins];
	double BfractionHX_rap[jpsi::kNbPTBins];
	double errBfractionHX_rap[jpsi::kNbPTBins];
	double nBackground_rap[jpsi::kNbPTBins];
	double errnBackground_rap[jpsi::kNbPTBins];
	double sigOvbkg_rap[jpsi::kNbPTBins];
	double errsigOvbkg_rap[jpsi::kNbPTBins];
	double lambda_theta_CS_PR_rap[jpsi::kNbPTBins];
	double errlambda_theta_CS_PR_rap[jpsi::kNbPTBins];
	double lambda_phi_CS_PR_rap[jpsi::kNbPTBins];
	double errlambda_phi_CS_PR_rap[jpsi::kNbPTBins];
	double lambda_theta_HX_PR_rap[jpsi::kNbPTBins];
	double errlambda_theta_HX_PR_rap[jpsi::kNbPTBins];
	double lambda_phi_HX_PR_rap[jpsi::kNbPTBins];
	double errlambda_phi_HX_PR_rap[jpsi::kNbPTBins];
	double invariant_lambda_CS_PR_rap[jpsi::kNbPTBins];
	double errinvariant_lambda_CS_PR_rap[jpsi::kNbPTBins];
	double invariant_lambda_HX_PR_rap[jpsi::kNbPTBins];
	double errinvariant_lambda_HX_PR_rap[jpsi::kNbPTBins];
	double invariant_F_CS_PR_rap[jpsi::kNbPTBins];
	double errinvariant_F_CS_PR_rap[jpsi::kNbPTBins];
	double invariant_F_HX_PR_rap[jpsi::kNbPTBins];
	double errinvariant_F_HX_PR_rap[jpsi::kNbPTBins];


	Double_t XS_Bfrac_mid[3] = {0.162, 0.257, 0.369};
	Double_t errXS_Bfrac_mid[3] = {0.071, 0.033, 0.041};
	Double_t XS_pTMean_mid[3] = {5.59, 7.88, 13.53};
	Double_t XS_Bfrac_fwd[5] = {0.098, 0.112, 0.165, 0.203, 0.331};
	Double_t errXS_Bfrac_fwd[5] = {0.058, 0.024, 0.029, 0.029, 0.057};
	Double_t XS_pTMean_fwd[5] = {1.23, 2.86, 4.89, 7.59, 13.14};



	for(int ptBin = 0; ptBin < jpsi::kNbPTBins; ++ptBin) {
	  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
      
      
      
      
      cout<<"pT"<<ptBin+1<<"rapidity"<<yBin+1<<endl;

	  char DirectoryPath[200];
	  sprintf(DirectoryPath,"pt%d_rapidity%d",ptBin+1,yBin+1);
	  char DirectoryPathPol[200];
	  sprintf(DirectoryPathPol,"PromptPolarizationModel",ptBin+1,yBin+1);



	  
	  TDirectory *InputDirectoryCS = (TDirectory*)fInputCS->GetDirectory(DirectoryPath);
	  TDirectory *InputDirectoryHX = (TDirectory*)fInputHX->GetDirectory(DirectoryPath);

	  if(!InputDirectoryCS || !InputDirectoryHX) continue;

	  TDirectory *InputDirectoryPolCS = (TDirectory*)InputDirectoryCS->GetDirectory(DirectoryPathPol);
	  TDirectory *InputDirectoryPolHX = (TDirectory*)InputDirectoryHX->GetDirectory(DirectoryPathPol);

	  if(!InputDirectoryCS || !InputDirectoryHX) continue;

	  nPromptCS = (RooRealVar*)InputDirectoryCS->Get("nPrompt");
	  std::cout << "nPrompt CS: " <<  nPromptCS << std::endl;
	  nPromptCS_[ptBin][yBin]=nPromptCS->getVal();
	  errnPromptCS_[ptBin][yBin]=nPromptCS->getError();

	  nNonPromptCS = (RooRealVar*)InputDirectoryCS->Get("nNonPrompt");
	  std::cout << "nNonPrompt CS: " << nNonPromptCS << std::endl;
	  nNonPromptCS_[ptBin][yBin]=nNonPromptCS->getVal();
	  errnNonPromptCS_[ptBin][yBin]=nNonPromptCS->getError();

	  nPromptHX = (RooRealVar*)InputDirectoryHX->Get("nPrompt");
	  std::cout << "nPrompt HX: " << nPromptHX << std::endl;
	  nPromptHX_[ptBin][yBin]=nPromptHX->getVal();
	  errnPromptHX_[ptBin][yBin]=nPromptHX->getError();

	  nNonPromptHX = (RooRealVar*)InputDirectoryHX->Get("nNonPrompt");
	  std::cout << "nNonPrompt HX: " << nNonPromptHX << std::endl;
	  nNonPromptHX_[ptBin][yBin]=nNonPromptHX->getVal();
	  errnNonPromptHX_[ptBin][yBin]=nNonPromptHX->getError();

	  nBackground = (RooRealVar*)InputDirectoryCS->Get("nBackground");
	  std::cout <<"nBackground: " <<  nBackground << std::endl;
	  nBackground_[ptBin][yBin]=nBackground->getVal();
	  errnBackground_[ptBin][yBin]=nBackground->getError();

	  BfractionCS_[ptBin][yBin]=nNonPromptCS_[ptBin][yBin]/(nNonPromptCS_[ptBin][yBin]+nPromptCS_[ptBin][yBin]);
	  errBfractionCS_[ptBin][yBin]=(1.0/(nPromptCS_[ptBin][yBin]+nNonPromptCS_[ptBin][yBin])+nNonPromptCS_[ptBin][yBin]/pow(nPromptCS_[ptBin][yBin]+nNonPromptCS_[ptBin][yBin],2))*errnNonPromptCS_[ptBin][yBin]+nNonPromptCS_[ptBin][yBin]/pow(nPromptCS_[ptBin][yBin]+nNonPromptCS_[ptBin][yBin],2)*errnPromptCS_[ptBin][yBin];

	  std::cout << "B Fraction CS: " << BfractionCS_[ptBin][yBin] << std::endl;
	  std::cout << "B Fraction CS Error: " << errBfractionCS_[ptBin][yBin] << std::endl;

	  BfractionHX_[ptBin][yBin]=nNonPromptHX_[ptBin][yBin]/(nNonPromptHX_[ptBin][yBin]+nPromptHX_[ptBin][yBin]);
	  errBfractionHX_[ptBin][yBin]=(1.0/(nPromptHX_[ptBin][yBin]+nNonPromptHX_[ptBin][yBin])+nNonPromptHX_[ptBin][yBin]/pow(nPromptHX_[ptBin][yBin]+nNonPromptHX_[ptBin][yBin],2))*errnNonPromptHX_[ptBin][yBin]+nNonPromptHX_[ptBin][yBin]/pow(nPromptHX_[ptBin][yBin]+nNonPromptHX_[ptBin][yBin],2)*errnPromptHX_[ptBin][yBin];

	  std::cout << "B Fraction HX: " << BfractionHX_[ptBin][yBin] << std::endl;
	  std::cout << "B Fraction HX Error: " << errBfractionHX_[ptBin][yBin] << std::endl;

	  sigOvbkg_[ptBin][yBin]=(nNonPromptCS_[ptBin][yBin]+nPromptCS_[ptBin][yBin])/nBackground_[ptBin][yBin];
	  errsigOvbkg_[ptBin][yBin] = sigOvbkg_[ptBin][yBin] * sqrt(pow((errnNonPromptCS_[ptBin][yBin]+errnPromptCS_[ptBin][yBin])/(nNonPromptCS_[ptBin][yBin]+nPromptCS_[ptBin][yBin]),2) + pow(errnBackground_[ptBin][yBin]/nBackground_[ptBin][yBin],2));

	  std::cout << "Signal/Background: " << sigOvbkg_[ptBin][yBin] << std::endl;
	  std::cout << "S/B Error: " << errsigOvbkg_[ptBin][yBin] << std::endl;
	  if(InputDirectoryPolCS) {
	    lambda_theta_CS_PR = (RooRealVar*)InputDirectoryPolCS->Get("promptlambda_theta_CS");
	    std::cout << "Prompt Lambda Theta CS: " << lambda_theta_CS_PR << std::endl;
	    lambda_theta_CS_PR_[ptBin][yBin]=lambda_theta_CS_PR->getVal();
	    errlambda_theta_CS_PR_[ptBin][yBin]=lambda_theta_CS_PR->getError();
	    
	    lambda_phi_CS_PR = (RooRealVar*)InputDirectoryPolCS->Get("promptlambda_phi_CS");
	    std::cout <<  "Prompt Lambda Phi CS: " << lambda_phi_CS_PR << std::endl;
	    lambda_phi_CS_PR_[ptBin][yBin]=lambda_phi_CS_PR->getVal();
	    errlambda_phi_CS_PR_[ptBin][yBin]=lambda_phi_CS_PR->getError();
	  }

	  if(InputDirectoryPolHX) {
	    lambda_theta_HX_PR = (RooRealVar*)InputDirectoryPolHX->Get("promptlambda_theta_HX");
	    std::cout <<  "Prompt Lambda Theta HX: " <<lambda_theta_HX_PR << std::endl;
	    lambda_theta_HX_PR_[ptBin][yBin]=lambda_theta_HX_PR->getVal();
	    errlambda_theta_HX_PR_[ptBin][yBin]=lambda_theta_HX_PR->getError();
	    
	    lambda_phi_HX_PR = (RooRealVar*)InputDirectoryPolHX->Get("promptlambda_phi_HX");
	    std::cout <<  "Prompt Lambda Phi HX: " <<lambda_phi_HX_PR << std::endl;
	    lambda_phi_HX_PR_[ptBin][yBin]=lambda_phi_HX_PR->getVal();
	    errlambda_phi_HX_PR_[ptBin][yBin]=lambda_phi_HX_PR->getError();
	  }
	  invariant_lambda_CS_PR_[ptBin][yBin]=(lambda_theta_CS_PR_[ptBin][yBin]+3*lambda_phi_CS_PR_[ptBin][yBin])/(1-lambda_phi_CS_PR_[ptBin][yBin]);
	  invariant_lambda_HX_PR_[ptBin][yBin]=(lambda_theta_HX_PR_[ptBin][yBin]+3*lambda_phi_HX_PR_[ptBin][yBin])/(1-lambda_phi_HX_PR_[ptBin][yBin]);
	  invariant_F_CS_PR_[ptBin][yBin]=(1+lambda_theta_CS_PR_[ptBin][yBin]+2*lambda_phi_CS_PR_[ptBin][yBin])/(3-lambda_theta_CS_PR_[ptBin][yBin]);
	  invariant_F_HX_PR_[ptBin][yBin]=(1+lambda_theta_HX_PR_[ptBin][yBin]+2*lambda_phi_HX_PR_[ptBin][yBin])/(3-lambda_theta_HX_PR_[ptBin][yBin]);
	  
	  
	  }	  
	}

////////////////////////////////// 	B - FRACTION CS ///////////////////////////////////////////////////////////////////////////////////

	TLegend* bfracLegendCS=new TLegend(0.15,0.7,0.65,0.89);
	bfracLegendCS->SetFillColor(kWhite);
	bfracLegendCS->SetTextFont(72);
	bfracLegendCS->SetTextSize(0.02);
	bfracLegendCS->SetBorderSize(0);

	char legendrap[200];

	TCanvas *BfractionCanvasCS = new TCanvas("fraction of J/#psi from B hadrons","fraction of J/#psi from B hadrons",1000,700);
	BfractionCanvasCS->SetFillColor(kWhite);
//	BfractionCanvasCS->SetGrid();
	BfractionCanvasCS->GetFrame()->SetFillColor(kWhite);
	BfractionCanvasCS->GetFrame()->SetBorderSize(0);
//	BfractionCanvasCS->SetLeftMargin(0.15) ;

	TH1F *BbfracHistoCS = BfractionCanvasCS->DrawFrame(0,0,20,1);
	BbfracHistoCS->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	BbfracHistoCS->SetYTitle("fraction of J/#psi from B hadrons");
	BbfracHistoCS->GetYaxis()->SetTitleOffset(1.5);

	for(int yBin=0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
	  for(int ptBin=0; ptBin < jpsi::kNbPTBins; ++ptBin) {
	    ptMean[ptBin]=jpsi::pTWCentre_rap[yBin][ptBin];
	    errptMean[ptBin]=0;
	    BfractionCS_rap[ptBin]=BfractionCS_[ptBin][yBin];
	    if (ptBin==0 && yBin==0) BfractionCS_rap[ptBin] = -1;
	    if (ptBin==1 && yBin==0) BfractionCS_rap[ptBin] = -1;
	    if (ptBin==2 && yBin==1) BfractionCS_rap[ptBin] = -1;
	    if (ptBin==3 && yBin==1) BfractionCS_rap[ptBin] = -1;
	    
	    errBfractionCS_rap[ptBin]=errBfractionCS_[ptBin][yBin];
	  }

	  TGraphErrors *bfracGraphCS = new TGraphErrors(jpsi::kNbPTBins,ptMean,BfractionCS_rap,errptMean,errBfractionCS_rap);
	  bfracGraphCS->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	  //	bfracGraphCS->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	  bfracGraphCS->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	  sprintf(legendrap,"Rapidity %d",yBin+1);
	  bfracLegendCS->AddEntry(bfracGraphCS,legendrap,"ple");
	  bfracGraphCS->Draw("P");
	  
	}
	
	TGraphErrors *bfracGraphCS = new TGraphErrors(3,XS_pTMean_mid,XS_Bfrac_mid,errptMean,errXS_Bfrac_mid);
	bfracGraphCS->SetMarkerColor(kBlack);
	//	bfracGraphCS->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	bfracGraphCS->SetMarkerStyle(13);
	sprintf(legendrap,"XS cross check, |y|<1.4 (100nb^{-1})");
	bfracLegendCS->AddEntry(bfracGraphCS,legendrap,"ple");
	bfracGraphCS->Draw("P");
	
	TGraphErrors *bfracGraphCS2 = new TGraphErrors(5,XS_pTMean_fwd,XS_Bfrac_fwd,errptMean,errXS_Bfrac_fwd);
	bfracGraphCS2->SetMarkerColor(kBlue);
	//	bfracGraphCS2->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	bfracGraphCS2->SetMarkerStyle(24);
	sprintf(legendrap,"XS cross check, 1.4<|y|<2.4 (100nb^{-1})");
	bfracLegendCS->AddEntry(bfracGraphCS2,legendrap,"ple");
	bfracGraphCS2->Draw("P");
	
	bfracLegendCS->Draw();
	
	char Filename[200];
	sprintf(Filename,"Plots/pT_BfractionCS.png");
	
	BfractionCanvasCS->SaveAs(Filename);
	
	BfractionCanvasCS->Close();

////////////////////////////////// 	B - FRACTION HX ///////////////////////////////////////////////////////////////////////////////////


	TLegend* bfracLegendHX=new TLegend(0.15,0.7,0.65,0.89);
	bfracLegendHX->SetFillColor(kWhite);
	bfracLegendHX->SetTextFont(72);
	bfracLegendHX->SetTextSize(0.02);
	bfracLegendHX->SetBorderSize(0);


	TCanvas *BfractionCanvasHX = new TCanvas("fraction of J/#psi from B hadrons","fraction of J/#psi from B hadrons",1000,700);
	BfractionCanvasHX->SetFillColor(kWhite);
	//	BfractionCanvasHX->SetGrid();
	BfractionCanvasHX->GetFrame()->SetFillColor(kWhite);
	BfractionCanvasHX->GetFrame()->SetBorderSize(0);
	//	BfractionCanvasHX->SetLeftMargin(0.15) ;
	
	TH1F *BbfracHistoHX = BfractionCanvasHX->DrawFrame(0,0,20,1);
	BbfracHistoHX->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	BbfracHistoHX->SetYTitle("fraction of J/#psi from B hadrons");
	BbfracHistoHX->GetYaxis()->SetTitleOffset(1.5);
	
	for(int yBin=0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
	  for(int ptBin=0; ptBin < jpsi::kNbPTBins; ++ptBin) {
	    ptMean[ptBin]=jpsi::pTWCentre_rap[yBin][ptBin];
	    errptMean[ptBin]=0;
	    BfractionHX_rap[ptBin]=BfractionHX_[ptBin][yBin];
	    if (ptBin==0 && yBin==0) BfractionHX_rap[ptBin] = -1;
	    if (ptBin==1 && yBin==0) BfractionHX_rap[ptBin] = -1;
	    if (ptBin==2 && yBin==1) BfractionHX_rap[ptBin] = -1;
	    if (ptBin==3 && yBin==1) BfractionHX_rap[ptBin] = -1;
	    
	    errBfractionHX_rap[ptBin]=errBfractionHX_[ptBin][yBin];
	  }
	  
	  TGraphErrors *bfracGraphHX = new TGraphErrors(jpsi::kNbPTBins,ptMean,BfractionHX_rap,errptMean,errBfractionHX_rap);
	  bfracGraphHX->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	  //	bfracGraphHX->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	  bfracGraphHX->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	  sprintf(legendrap,"Rapidity %d",yBin+1);
	  bfracLegendHX->AddEntry(bfracGraphHX,legendrap,"ple");
	  bfracGraphHX->Draw("P");
	  
	}
	


	TGraphErrors *bfracGraphHX3 = new TGraphErrors(3,XS_pTMean_mid,XS_Bfrac_mid,errptMean,errXS_Bfrac_mid);
	bfracGraphHX3->SetMarkerColor(kBlack);
	//	bfracGraphHX3->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	bfracGraphHX3->SetMarkerStyle(13);
	sprintf(legendrap,"XS cross check, |y|<1.4 (100nb^{-1})");
	bfracLegendHX->AddEntry(bfracGraphHX3,legendrap,"ple");
	bfracGraphHX3->Draw("P");
	
	TGraphErrors *bfracGraphHX4 = new TGraphErrors(5,XS_pTMean_fwd,XS_Bfrac_fwd,errptMean,errXS_Bfrac_fwd);
	bfracGraphHX4->SetMarkerColor(kBlue);
	//	bfracGraphHX4->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	bfracGraphHX4->SetMarkerStyle(24);
	sprintf(legendrap,"XS cross check, 1.4<|y|<2.4 (100nb^{-1})");
	bfracLegendHX->AddEntry(bfracGraphHX4,legendrap,"ple");
	bfracGraphHX4->Draw("P");
	
	
	
	bfracLegendHX->Draw();
	
	sprintf(Filename,"Plots/pT_BfractionHX.png");
	
	BfractionCanvasHX->SaveAs(Filename);
	
	BfractionCanvasHX->Close();
	
////////////////////////////////// 	Number of Background Events ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* BackgroundLegend=new TLegend(0.8,0.6,0.9,0.89);
	 BackgroundLegend->SetFillColor(kWhite);
	 BackgroundLegend->SetTextFont(72);
	 BackgroundLegend->SetTextSize(0.02);
	 BackgroundLegend->SetBorderSize(0);

	 TCanvas *BackgroundCanvas = new TCanvas("number of Background events","",1000,700);
	 BackgroundCanvas->SetFillColor(kWhite);
//   BackgroundCanvas->SetGrid();
	 BackgroundCanvas->GetFrame()->SetFillColor(kWhite);
	 BackgroundCanvas->GetFrame()->SetBorderSize(0);
	 BackgroundCanvas->SetLeftMargin(0.15) ;

	 TH1F *BackgroundhistoHisto = BackgroundCanvas->DrawFrame(0,0,30,65000);
	 BackgroundhistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 BackgroundhistoHisto->SetYTitle("number of Background events");
	 BackgroundhistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 	  		for(int yBin=0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
	 	  		for(int ptBin=0; ptBin < jpsi::kNbPTBins; ++ptBin) {
	 	  					   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin][ptBin];
		  					   errptMean[ptBin]=0;
	 	  					nBackground_rap[ptBin]=nBackground_[ptBin][yBin];
	 	  					errnBackground_rap[ptBin]=errnBackground_[ptBin][yBin];
	 	  					  }

	 TGraphErrors *BackgroundGraph = new TGraphErrors(jpsi::kNbPTBins,ptMean,nBackground_rap,errptMean,errnBackground_rap);
	 BackgroundGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
//	 BackgroundGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 BackgroundGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 BackgroundLegend->AddEntry(BackgroundGraph,legendrap,"ple");
	 BackgroundGraph->Draw("P");

	 	  		}

	 	  		BackgroundLegend->Draw();

	 	  		  sprintf(Filename,"Plots/pT_nBackground.png");

	 	  		BackgroundCanvas->SaveAs(Filename);

	 	  		BackgroundCanvas->Close();

////////////////////////////////// 	SIG/BKG RATIO ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* sigOvbkgLegend=new TLegend(0.8,0.6,0.9,0.89);
	 sigOvbkgLegend->SetFillColor(kWhite);
	 sigOvbkgLegend->SetTextFont(72);
	 sigOvbkgLegend->SetTextSize(0.02);
	 sigOvbkgLegend->SetBorderSize(0);

	 TCanvas *sigOvbkgCanvas = new TCanvas("Sig/Bkg ratio","",1000,700);
	 sigOvbkgCanvas->SetFillColor(kWhite);
//   sigOvbkgCanvas->SetGrid();
	 sigOvbkgCanvas->GetFrame()->SetFillColor(kWhite);
	 sigOvbkgCanvas->GetFrame()->SetBorderSize(0);
//	 sigOvbkgCanvas->SetLeftMargin(0.15) ;

	 TH1F *sigOvbkghistoHisto = sigOvbkgCanvas->DrawFrame(0,0,30,20);
	 sigOvbkghistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 sigOvbkghistoHisto->SetYTitle("Sig/Bkg ratio");
	 sigOvbkghistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 	  		for(int yBin=0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
	 	  		for(int ptBin=0; ptBin < jpsi::kNbPTBins; ++ptBin) {
	 	  					   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin][ptBin];
		  					   errptMean[ptBin]=0;
	 	  					sigOvbkg_rap[ptBin]=sigOvbkg_[ptBin][yBin];
	 	  					errsigOvbkg_rap[ptBin]=errsigOvbkg_[ptBin][yBin];
	 	  					  }

	 TGraphErrors *sigOvbkgGraph = new TGraphErrors(jpsi::kNbPTBins,ptMean,sigOvbkg_rap,errptMean,errsigOvbkg_rap);
	 sigOvbkgGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
//	 sigOvbkgGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 sigOvbkgGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 sigOvbkgLegend->AddEntry(sigOvbkgGraph,legendrap,"ple");
	 sigOvbkgGraph->Draw("P");

	 	  		}

	 	  		sigOvbkgLegend->Draw();

	 	  		  sprintf(Filename,"Plots/pT_sigOvbkg.png");

	 	  		sigOvbkgCanvas->SaveAs(Filename);

	 	  		sigOvbkgCanvas->Close();

////////////////////////////////// LAMBDA THETA PROMPT CS ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* LambdaThetaPRCSLegend=new TLegend(0.8,0.6,0.9,0.89);
	 LambdaThetaPRCSLegend->SetFillColor(kWhite);
	 LambdaThetaPRCSLegend->SetTextFont(72);
	 LambdaThetaPRCSLegend->SetTextSize(0.02);
	 LambdaThetaPRCSLegend->SetBorderSize(0);

	 TCanvas *LambdaThetaPRCSCanvas = new TCanvas("Prompt #lambda_{#theta_{CS}}","",1000,700);
	 LambdaThetaPRCSCanvas->SetFillColor(kWhite);
     LambdaThetaPRCSCanvas->SetGrid(0,1);
	 LambdaThetaPRCSCanvas->GetFrame()->SetFillColor(kWhite);
	 LambdaThetaPRCSCanvas->GetFrame()->SetBorderSize(0);
//	 LambdaThetaPRCSCanvas->SetLeftMargin(0.15) ;

	 TH1F *LambdaThetaPRCShistoHisto = LambdaThetaPRCSCanvas->DrawFrame(0,-1.3,30,1.3);
	 LambdaThetaPRCShistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 LambdaThetaPRCShistoHisto->SetYTitle("Prompt #lambda_{#theta_{CS}}");
	 LambdaThetaPRCShistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 	  		for(int yBin=0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
	 	  		for(int ptBin=0; ptBin < jpsi::kNbPTBins; ++ptBin) {
	 	  					   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin][ptBin];
		  					   errptMean[ptBin]=0;
		  					 lambda_theta_CS_PR_rap[ptBin]=lambda_theta_CS_PR_[ptBin][yBin];
	 	  					errlambda_theta_CS_PR_rap[ptBin]=errlambda_theta_CS_PR_[ptBin][yBin];
	 	  					  }

	 TGraphErrors *LambdaThetaPRCSGraph = new TGraphErrors(jpsi::kNbPTBins,ptMean,lambda_theta_CS_PR_rap,errptMean,errlambda_theta_CS_PR_rap);
	 LambdaThetaPRCSGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPRCSGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPRCSGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 LambdaThetaPRCSLegend->AddEntry(LambdaThetaPRCSGraph,legendrap,"ple");
	 LambdaThetaPRCSGraph->Draw("P");

	 	  		}

	 	  		LambdaThetaPRCSLegend->Draw();

	 	  		  sprintf(Filename,"Plots/pT_LambdaThetaPRCS.png");

	 	  		LambdaThetaPRCSCanvas->SaveAs(Filename);

	 	  		LambdaThetaPRCSCanvas->Close();

////////////////////////////////// LAMBDA PHI PROMPT CS ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* LambdaPhiPRCSLegend=new TLegend(0.8,0.6,0.9,0.89);
	 LambdaPhiPRCSLegend->SetFillColor(kWhite);
	 LambdaPhiPRCSLegend->SetTextFont(72);
	 LambdaPhiPRCSLegend->SetTextSize(0.02);
	 LambdaPhiPRCSLegend->SetBorderSize(0);

	 TCanvas *LambdaPhiPRCSCanvas = new TCanvas("Prompt #lambda_{#phi_{CS}}","",1000,700);
	 LambdaPhiPRCSCanvas->SetFillColor(kWhite);
     LambdaPhiPRCSCanvas->SetGrid(0,1);
	 LambdaPhiPRCSCanvas->GetFrame()->SetFillColor(kWhite);
	 LambdaPhiPRCSCanvas->GetFrame()->SetBorderSize(0);
//	 LambdaPhiPRCSCanvas->SetLeftMargin(0.15) ;

	 TH1F *LambdaPhiPRCShistoHisto = LambdaPhiPRCSCanvas->DrawFrame(0,-1.3,30,1.3);
	 LambdaPhiPRCShistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 LambdaPhiPRCShistoHisto->SetYTitle("Prompt #lambda_{#phi_{CS}}");
	 LambdaPhiPRCShistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 	  		for(int yBin=0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
	 	  		for(int ptBin=0; ptBin < jpsi::kNbPTBins; ++ptBin) {
	 	  					   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin][ptBin];
		  					   errptMean[ptBin]=0;
		  					 lambda_phi_CS_PR_rap[ptBin]=lambda_phi_CS_PR_[ptBin][yBin];
	 	  					errlambda_phi_CS_PR_rap[ptBin]=errlambda_phi_CS_PR_[ptBin][yBin];
	 	  					  }

	 TGraphErrors *LambdaPhiPRCSGraph = new TGraphErrors(jpsi::kNbPTBins,ptMean,lambda_phi_CS_PR_rap,errptMean,errlambda_phi_CS_PR_rap);
	 LambdaPhiPRCSGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaPhiPRCSGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaPhiPRCSGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 LambdaPhiPRCSLegend->AddEntry(LambdaPhiPRCSGraph,legendrap,"ple");
	 LambdaPhiPRCSGraph->Draw("P");

	 	  		}

	 	  		LambdaPhiPRCSLegend->Draw();

	 	  		  sprintf(Filename,"Plots/pT_LambdaPhiPRCS.png");

	 	  		LambdaPhiPRCSCanvas->SaveAs(Filename);

	 	  		LambdaPhiPRCSCanvas->Close();

////////////////////////////////// LAMBDA THETA PROMPT HX ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* LambdaThetaPRHXLegend=new TLegend(0.8,0.6,0.9,0.89);
	 LambdaThetaPRHXLegend->SetFillColor(kWhite);
	 LambdaThetaPRHXLegend->SetTextFont(72);
	 LambdaThetaPRHXLegend->SetTextSize(0.02);
	 LambdaThetaPRHXLegend->SetBorderSize(0);

	 TCanvas *LambdaThetaPRHXCanvas = new TCanvas("Prompt #lambda_{#theta_{HX}}","",1000,700);
	 LambdaThetaPRHXCanvas->SetFillColor(kWhite);
     LambdaThetaPRHXCanvas->SetGrid(0,1);
	 LambdaThetaPRHXCanvas->GetFrame()->SetFillColor(kWhite);
	 LambdaThetaPRHXCanvas->GetFrame()->SetBorderSize(0);
//	 LambdaThetaPRHXCanvas->SetLeftMargin(0.15) ;

	 TH1F *LambdaThetaPRHXhistoHisto = LambdaThetaPRHXCanvas->DrawFrame(0,-1.3,30,1.3);
	 LambdaThetaPRHXhistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 LambdaThetaPRHXhistoHisto->SetYTitle("Prompt #lambda_{#theta_{HX}}");
	 LambdaThetaPRHXhistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 	  		for(int yBin=0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
	 	  		for(int ptBin=0; ptBin < jpsi::kNbPTBins; ++ptBin) {
	 	  					   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin][ptBin];
		  					   errptMean[ptBin]=0;
		  					 lambda_theta_HX_PR_rap[ptBin]=lambda_theta_HX_PR_[ptBin][yBin];
	 	  					errlambda_theta_HX_PR_rap[ptBin]=errlambda_theta_HX_PR_[ptBin][yBin];
	 	  					  }

	 TGraphErrors *LambdaThetaPRHXGraph = new TGraphErrors(jpsi::kNbPTBins,ptMean,lambda_theta_HX_PR_rap,errptMean,errlambda_theta_HX_PR_rap);
	 LambdaThetaPRHXGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPRHXGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPRHXGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 LambdaThetaPRHXLegend->AddEntry(LambdaThetaPRHXGraph,legendrap,"ple");
	 LambdaThetaPRHXGraph->Draw("P");

	 	  		}

	 	  		LambdaThetaPRHXLegend->Draw();

	 	  		  sprintf(Filename,"Plots/pT_LambdaThetaPRHX.png");

	 	  		LambdaThetaPRHXCanvas->SaveAs(Filename);

	 	  		LambdaThetaPRHXCanvas->Close();

////////////////////////////////// LAMBDA PHI PROMPT HX ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* LambdaPhiPRHXLegend=new TLegend(0.8,0.6,0.9,0.89);
	 LambdaPhiPRHXLegend->SetFillColor(kWhite);
	 LambdaPhiPRHXLegend->SetTextFont(72);
	 LambdaPhiPRHXLegend->SetTextSize(0.02);
	 LambdaPhiPRHXLegend->SetBorderSize(0);

	 TCanvas *LambdaPhiPRHXCanvas = new TCanvas("Prompt #lambda_{#phi_{HX}}","",1000,700);
	 LambdaPhiPRHXCanvas->SetFillColor(kWhite);
     LambdaPhiPRHXCanvas->SetGrid(0,1);
	 LambdaPhiPRHXCanvas->GetFrame()->SetFillColor(kWhite);
	 LambdaPhiPRHXCanvas->GetFrame()->SetBorderSize(0);
//	 LambdaPhiPRHXCanvas->SetLeftMargin(0.15) ;

	 TH1F *LambdaPhiPRHXhistoHisto = LambdaPhiPRHXCanvas->DrawFrame(0,-1.3,30,1.3);
	 LambdaPhiPRHXhistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 LambdaPhiPRHXhistoHisto->SetYTitle("Prompt #lambda_{#phi_{HX}}");
	 LambdaPhiPRHXhistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 	  		for(int yBin=0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
	 	  		for(int ptBin=0; ptBin < jpsi::kNbPTBins; ++ptBin) {
	 	  					   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin][ptBin];
		  					   errptMean[ptBin]=0;
		  					 lambda_phi_HX_PR_rap[ptBin]=lambda_phi_HX_PR_[ptBin][yBin];
	 	  					errlambda_phi_HX_PR_rap[ptBin]=errlambda_phi_HX_PR_[ptBin][yBin];
	 	  					  }

	 TGraphErrors *LambdaPhiPRHXGraph = new TGraphErrors(jpsi::kNbPTBins,ptMean,lambda_phi_HX_PR_rap,errptMean,errlambda_phi_HX_PR_rap);
	 LambdaPhiPRHXGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaPhiPRHXGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaPhiPRHXGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 LambdaPhiPRHXLegend->AddEntry(LambdaPhiPRHXGraph,legendrap,"ple");
	 LambdaPhiPRHXGraph->Draw("P");

	 	  		}

	 	  		LambdaPhiPRHXLegend->Draw();

	 	  		  sprintf(Filename,"Plots/pT_LambdaPhiPRHX.png");

	 	  		LambdaPhiPRHXCanvas->SaveAs(Filename);

	 	  		LambdaPhiPRHXCanvas->Close();

////////////////////////////////// INVARIANT LAMBDA PROMPT ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* PromptInvariantLambdaLegend=new TLegend(0.7,0.6,0.9,0.89);
	 PromptInvariantLambdaLegend->SetFillColor(kWhite);
	 PromptInvariantLambdaLegend->SetTextFont(72);
	 PromptInvariantLambdaLegend->SetTextSize(0.02);
	 PromptInvariantLambdaLegend->SetBorderSize(0);

	 TCanvas *PromptInvariantLambdaCanvas = new TCanvas("Invariant Prompt #lambda_{CS}","",1000,700);
	 PromptInvariantLambdaCanvas->SetFillColor(kWhite);
     PromptInvariantLambdaCanvas->SetGrid(0,1);
	 PromptInvariantLambdaCanvas->GetFrame()->SetFillColor(kWhite);
	 PromptInvariantLambdaCanvas->GetFrame()->SetBorderSize(0);
//	 PromptInvariantLambdaCanvas->SetLeftMargin(0.15) ;

	 TH1F *PromptInvariantLambdahistoHisto = PromptInvariantLambdaCanvas->DrawFrame(0,-1.3,30,1.3);
	 PromptInvariantLambdahistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 PromptInvariantLambdahistoHisto->SetYTitle("Invariant Prompt #lambda_{CS}");
	 PromptInvariantLambdahistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 	  		for(int yBin=0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
	 	  		for(int ptBin=0; ptBin < jpsi::kNbPTBins; ++ptBin) {
	 	  					   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin][ptBin];
		  					   errptMean[ptBin]=0;
		  					invariant_lambda_CS_PR_rap[ptBin]=invariant_lambda_CS_PR_[ptBin][yBin];
	 	  					errinvariant_lambda_CS_PR_rap[ptBin]=0;
	 	  					invariant_lambda_HX_PR_rap[ptBin]=invariant_lambda_HX_PR_[ptBin][yBin];
	 	  					errinvariant_lambda_HX_PR_rap[ptBin]=0;
	 	  					  }

	 TGraphErrors *PromptInvariantLambdaCSGraph = new TGraphErrors(jpsi::kNbPTBins,ptMean,invariant_lambda_CS_PR_rap,errptMean,errinvariant_lambda_CS_PR_rap);
	 PromptInvariantLambdaCSGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantLambdaCSGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantLambdaCSGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d, CS",yBin+1);
	 PromptInvariantLambdaLegend->AddEntry(PromptInvariantLambdaCSGraph,legendrap,"ple");
	 PromptInvariantLambdaCSGraph->Draw("P");

	 TGraphErrors *PromptInvariantLambdaHXGraph = new TGraphErrors(jpsi::kNbPTBins,ptMean,invariant_lambda_HX_PR_rap,errptMean,errinvariant_lambda_HX_PR_rap);
	 PromptInvariantLambdaHXGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantLambdaHXGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantLambdaHXGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]+3);
	 sprintf(legendrap,"Rapidity %d, HX",yBin+1);
	 PromptInvariantLambdaLegend->AddEntry(PromptInvariantLambdaHXGraph,legendrap,"ple");
	 PromptInvariantLambdaHXGraph->Draw("P");

	 	  		}

	 	  		PromptInvariantLambdaLegend->Draw();

	 	  		  sprintf(Filename,"Plots/pT_PromptInvariantLambda.png");

	 	  		PromptInvariantLambdaCanvas->SaveAs(Filename);

	 	  		PromptInvariantLambdaCanvas->Close();


////////////////////////////////// INVARIANT F PROMPT ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* PromptInvariantFLegend=new TLegend(0.7,0.6,0.9,0.89);
	 PromptInvariantFLegend->SetFillColor(kWhite);
	 PromptInvariantFLegend->SetTextFont(72);
	 PromptInvariantFLegend->SetTextSize(0.02);
	 PromptInvariantFLegend->SetBorderSize(0);

	 TCanvas *PromptInvariantFCanvas = new TCanvas("Invariant Prompt F_{CS}","",1000,700);
	 PromptInvariantFCanvas->SetFillColor(kWhite);
     PromptInvariantFCanvas->SetGrid(0,1);
	 PromptInvariantFCanvas->GetFrame()->SetFillColor(kWhite);
	 PromptInvariantFCanvas->GetFrame()->SetBorderSize(0);
//	 PromptInvariantFCanvas->SetLeftMargin(0.15) ;

	 TH1F *PromptInvariantFhistoHisto = PromptInvariantFCanvas->DrawFrame(0,-1.3,30,1.3);
	 PromptInvariantFhistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 PromptInvariantFhistoHisto->SetYTitle("Invariant Prompt F_{CS}");
	 PromptInvariantFhistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 	  		for(int yBin=0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
	 	  		for(int ptBin=0; ptBin < jpsi::kNbPTBins; ++ptBin) {
	 	  					   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin][ptBin];
		  					   errptMean[ptBin]=0;
		  					 invariant_F_CS_PR_rap[ptBin]=invariant_F_CS_PR_[ptBin][yBin];
	 	  					errinvariant_F_CS_PR_rap[ptBin]=0;
	 	  					invariant_F_HX_PR_rap[ptBin]=invariant_F_HX_PR_[ptBin][yBin];
	 	  					errinvariant_F_HX_PR_rap[ptBin]=0;
	 	  					  }

	 TGraphErrors *PromptInvariantFCSGraph = new TGraphErrors(jpsi::kNbPTBins,ptMean,invariant_F_CS_PR_rap,errptMean,errinvariant_F_CS_PR_rap);
	 PromptInvariantFCSGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantFCSGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantFCSGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d, CS",yBin+1);
	 PromptInvariantFLegend->AddEntry(PromptInvariantFCSGraph,legendrap,"ple");
	 PromptInvariantFCSGraph->Draw("P");

	 TGraphErrors *PromptInvariantFHXGraph = new TGraphErrors(jpsi::kNbPTBins,ptMean,invariant_F_HX_PR_rap,errptMean,errinvariant_F_HX_PR_rap);
	 PromptInvariantFHXGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantFHXGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantFHXGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]+3);
	 sprintf(legendrap,"Rapidity %d, HX",yBin+1);
	 PromptInvariantFLegend->AddEntry(PromptInvariantFHXGraph,legendrap,"ple");
	 PromptInvariantFHXGraph->Draw("P");

	 	  		}

	 	  		PromptInvariantFLegend->Draw();

	 	  		  sprintf(Filename,"Plots/pT_PromptInvariantF.png");

	 	  		PromptInvariantFCanvas->SaveAs(Filename);

	 	  		PromptInvariantFCanvas->Close();





  delete fInputCS;
  delete fInputHX;



  return 0;
}
