/*
 * BinExamplePlots.C
 *
 *  Created on: Apr 22, 2012
 *      Author: valentinknuenz
 */

#include <iostream>
#include <sstream>
#include <iomanip>


//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TFormula.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"


void plotSummary( TTree* Results, char parameter[200], char parameterTitle[200], int iState , int irapBin, int ipTbin, char DataID[200] ,double nSigma) {

	int nBins=100;
	TH1D* h_param   = new TH1D( "h_param", "h_param", nBins, -1,1);
	char DrawChar[200];
	sprintf(DrawChar,"%s>>h_param",parameter);
	Results->Draw(DrawChar);
	TH1D* h_param_   = new TH1D( "h_param_", "h_param_", nBins, h_param->GetMean()-nSigma*h_param->GetRMS(), h_param->GetMean()+nSigma*h_param->GetRMS());
	sprintf(DrawChar,"%s>>h_param_",parameter);
	Results->Draw(DrawChar);

	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	h_param_->SetXTitle(parameterTitle);
	gPad->SetFillColor(kWhite);
	h_param_->SetTitle(0);
	h_param_->Draw();

	char savename[200];
	sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/%s_rap%d_pT%d.pdf",DataID,iState,parameter,irapBin,ipTbin);
	c2->SaveAs(savename);

}


void BinExamplePlots(){

	gROOT->SetBatch();

	  char DataID[200];
	  sprintf(DataID,"Data_TowardsPRL_Aug11_FinalResults_1Sigma_AlteredPPD_BKGlinPLUSRestSquaredGauss_5nRand_BiasCorrection_1S2S3SAug12_1Sig");

	  bool PlotAcceptance=true;

	  char directory[200];
	  char Fig_directory[200];
	  sprintf(Fig_directory,"FigBuffer/BinExamplePlots");
	  gSystem->mkdir(Fig_directory);
	  sprintf(Fig_directory,"FigBuffer/BinExamplePlots/%s",DataID);
	  gSystem->mkdir(Fig_directory);

	  int pTmin=6;
	  int pTmax=10;
	  int rapmin=1;
	  int rapmax=2;
	  int nStatemin=1;
	  int nStatemax=3;


	  for(iState=nStatemin;iState<nStatemax+1;iState++){
		  sprintf(directory,"%s/Ups%dS",Fig_directory,iState);
		  gSystem->mkdir(directory);

		  for(irapBin=rapmin;irapBin<rapmax+1;irapBin++){
			  for(ipTbin=pTmin;ipTbin<pTmax+1;ipTbin++){


	  char Filename[200];
	  sprintf(Filename,"/scratch/knuenz/Polarization/Upsilon/Data/%s/results_%dSUps_rap%d_pT%d.root",DataID,iState,irapBin,ipTbin);
      TFile *resultFile = new TFile(Filename,"READ");

      TTree* Results=(TTree*)resultFile->Get("Results");
      TH1D* SubtractedBG_test=(TH1D*)resultFile->Get("SubtractedBG_test");

      double nSigma=4;

      for(iLam=1;iLam<13;iLam++){

      char parameter[200];
      char parameterTitle[200];

      if(iLam==1) {sprintf(parameter,"lthCS"); sprintf(parameterTitle,"#lambda_{#theta}^{CS}");}
      if(iLam==2) {sprintf(parameter,"lphCS"); sprintf(parameterTitle,"#lambda_{#phi}^{CS}");}
      if(iLam==3) {sprintf(parameter,"ltpCS"); sprintf(parameterTitle,"#lambda_{#theta#phi}^{CS}");}
      if(iLam==4) {sprintf(parameter,"ltildeCS"); sprintf(parameterTitle,"#tilde{#lambda}^{CS}");}

      if(iLam==5) {sprintf(parameter,"lthHX"); sprintf(parameterTitle,"#lambda_{#theta}^{HX}");}
      if(iLam==6) {sprintf(parameter,"lphHX"); sprintf(parameterTitle,"#lambda_{#phi}^{HX}");}
      if(iLam==7) {sprintf(parameter,"ltpHX"); sprintf(parameterTitle,"#lambda_{#theta#phi}^{HX}");}
      if(iLam==8) {sprintf(parameter,"ltildeHX"); sprintf(parameterTitle,"#tilde{#lambda}^{HX}");}

      if(iLam==9) {sprintf(parameter,"lthPX"); sprintf(parameterTitle,"#lambda_{#theta}^{PX}");}
      if(iLam==10) {sprintf(parameter,"lphPX"); sprintf(parameterTitle,"#lambda_{#phi}^{PX}");}
      if(iLam==11) {sprintf(parameter,"ltpPX"); sprintf(parameterTitle,"#lambda_{#theta#phi}^{PX}");}
      if(iLam==12) {sprintf(parameter,"ltildePX"); sprintf(parameterTitle,"#tilde{#lambda}^{PX}");}

      plotSummary( Results, parameter, parameterTitle, iState , irapBin, ipTbin, DataID , nSigma);

      }



      double lth_min=-1.1;
      double lth_max=1.1;
      double lph_min=-1.1;
      double lph_max=1.1;

      TCanvas *c1 = new TCanvas("c1", "c1", 10, 28, 580,571);
      c1->Range(-237.541,-66.47556,187.377,434.8609);
      c1->SetFillColor(0);
      c1->SetBorderMode(0);
      c1->SetBorderSize(0);
      c1->SetLeftMargin(0.1354167);
      c1->SetRightMargin(0.01736111);
      c1->SetTopMargin(0.01841621);
      c1->SetBottomMargin(0.1325967);
      c1->SetFrameBorderMode(0);

      TH2D* h_lph_vs_lth_CS   = new TH2D( "h_lph_vs_lth_CS", "h_lph_vs_lth_CS", 100, lth_min,lth_max,100,lph_min,lph_max);
      Results->Draw("lphCS:lthCS>>h_lph_vs_lth_CS");
      TH2D* h_lph_vs_lth_HX   = new TH2D( "h_lph_vs_lth_HX", "h_lph_vs_lth_HX", 100,lth_min,lth_max,100,lph_min,lph_max);
      Results->Draw("lphHX:lthHX>> h_lph_vs_lth_HX");
      TH2D* h_lph_vs_lth_PX   = new TH2D( "h_lph_vs_lth_PX", "h_lph_vs_lth_PX", 100, lth_min,lth_max,100,lph_min,lph_max);
      Results->Draw("lphPX:lthPX>> h_lph_vs_lth_PX");

      h_lph_vs_lth_CS->SetXTitle("#lambda_{#theta}");
      h_lph_vs_lth_CS->SetYTitle("#lambda_{#phi}");


      h_lph_vs_lth_CS->SetTitle(0);
      h_lph_vs_lth_HX->SetTitle(0);
      h_lph_vs_lth_PX->SetTitle(0);
      h_lph_vs_lth_CS->SetStats(0);
      h_lph_vs_lth_HX->SetStats(0);
      h_lph_vs_lth_PX->SetStats(0);

      h_lph_vs_lth_CS->SetMarkerStyle(20);
      h_lph_vs_lth_HX->SetMarkerStyle(20);
      h_lph_vs_lth_PX->SetMarkerStyle(20);
      h_lph_vs_lth_CS->SetMarkerSize(0.5);
      h_lph_vs_lth_HX->SetMarkerSize(0.5);
      h_lph_vs_lth_PX->SetMarkerSize(0.5);
      h_lph_vs_lth_CS->SetMarkerColor(kBlue);
      h_lph_vs_lth_HX->SetMarkerColor(kRed);
      h_lph_vs_lth_PX->SetMarkerColor(kGreen);

      TH2D* h_lph_vs_lth = (TH2D*)h_lph_vs_lth_CS->Clone();
      h_lph_vs_lth->SetName("h_lph_vs_lth");
      h_lph_vs_lth->Reset();

      h_lph_vs_lth->GetXaxis()->SetTitle("#lambda_{#vartheta}");
      h_lph_vs_lth->GetXaxis()->SetLabelOffset(0.028);
      h_lph_vs_lth->GetXaxis()->SetTitleSize(0.05);
      h_lph_vs_lth->GetXaxis()->SetTickLength(-0.03);
      h_lph_vs_lth->GetXaxis()->SetTitleOffset(1.15);
      h_lph_vs_lth->GetYaxis()->SetTitle("#lambda_{#varphi}");
      h_lph_vs_lth->GetYaxis()->SetLabelOffset(0.032);
      h_lph_vs_lth->GetYaxis()->SetTitleSize(0.05);
      h_lph_vs_lth->GetYaxis()->SetTickLength(-0.03);
      h_lph_vs_lth->GetYaxis()->SetTitleOffset(1.3);
      h_lph_vs_lth->Draw("");

      double x_background_ph_vs_th[4] = { lth_min+0.01, lth_min+0.01, lth_max-0.01, lth_max-0.01 };
      double y_background_ph_vs_th[4] = { lph_min+0.01, lph_max-0.01, lph_max-0.01, lph_min+0.01 };
      TPolyLine *background_ph_vs_th = new TPolyLine( 4, x_background_ph_vs_th, y_background_ph_vs_th );
      background_ph_vs_th->SetFillColor(kGray);
      background_ph_vs_th->SetLineStyle(0);
      background_ph_vs_th->Draw("f same");

      double x_triangle_ph_vs_th[3] = {-1.,  1., 1.};
      double y_triangle_ph_vs_th[3] = { 0., -1., 1.};
      TPolyLine *triangle_ph_vs_th = new TPolyLine( 3, x_triangle_ph_vs_th, y_triangle_ph_vs_th );
      triangle_ph_vs_th->SetFillColor(kWhite);
      triangle_ph_vs_th->SetLineStyle(0);
      triangle_ph_vs_th->Draw("f same");

      TLine* lh_p1 = new TLine( lth_min, 1., lth_max, 1. );
      lh_p1->SetLineWidth( 1 );
      lh_p1->SetLineStyle( 2 );
      lh_p1->SetLineColor( kGray+2 );
      lh_p1->Draw( "same" );

      TLine* lh_m1 = new TLine( lth_min, -1., lth_max, -1. );
      lh_m1->SetLineWidth( 1 );
      lh_m1->SetLineStyle( 2 );
      lh_m1->SetLineColor( kGray+2 );
      lh_m1->Draw( "same" );

      TLine* lh_0 = new TLine( lth_min, 0., lth_max, 0. );
      lh_0->SetLineWidth( 1 );
      lh_0->SetLineStyle( 2 );
      lh_0->SetLineColor( kGray+2 );
      lh_0->Draw( "same" );

      TLine* lv_p1 = new TLine( 1., lph_min, 1., lph_max );
      lv_p1->SetLineWidth( 1 );
      lv_p1->SetLineStyle( 2 );
      lv_p1->SetLineColor( kGray+2 );
      lv_p1->Draw( "same" );

      TLine* lv_m1 = new TLine( -1., lph_min, -1., lph_max );
      lv_m1->SetLineWidth( 1 );
      lv_m1->SetLineStyle( 2 );
      lv_m1->SetLineColor( kGray+2 );
      lv_m1->Draw( "same" );

      TLine* lv_0 = new TLine( 0., lph_min, 0., lph_max );
      lv_0->SetLineWidth( 1 );
      lv_0->SetLineStyle( 2 );
      lv_0->SetLineColor( kGray+2 );
      lv_0->Draw( "same" );



      // CS frame

/*      setContourHistogram ( h_lph_vs_lth_CS );
      h_lph_vs_lth_CS->SetLineColor( kBlue );
      h_lph_vs_lth_CS->SetLineWidth( 2 );
      h_lph_vs_lth_CS->SetLineStyle( 1 );
      h_lph_vs_lth_CS->Draw( "same" );

      // HX frame

      setContourHistogram ( h_lph_vs_lth_HX );
      h_lph_vs_lth_HX->SetLineColor( kRed );
      h_lph_vs_lth_HX->SetLineWidth( 2 );
      h_lph_vs_lth_HX->SetLineStyle( 1 );
      h_lph_vs_lth_HX->Draw( "cont2, same" );

      // PX frame

      setContourHistogram ( h_lph_vs_lth_PX );
      h_lph_vs_lth_PX->SetLineColor( kGreen+2 );
      h_lph_vs_lth_PX->SetLineWidth( 2 );
      h_lph_vs_lth_PX->SetLineStyle( 1 );
      h_lph_vs_lth_PX->Draw( "cont2, same" );
*/

      h_lph_vs_lth_CS->Draw("same");
      h_lph_vs_lth_HX->Draw("same");
      h_lph_vs_lth_PX->Draw("same");

 	 TLegend* plotLegend=new TLegend(0.15,0.8,0.25,0.95);
    	 plotLegend->SetFillColor(kWhite);
//    	 plotLegend->SetTextFont(72);
    	 plotLegend->SetTextSize(0.035);
    	 plotLegend->SetBorderSize(1);
    	 char legendentry[200];
    	 sprintf(legendentry,"CS");
    	 plotLegend->AddEntry(h_lph_vs_lth_CS,legendentry,"p");
    	 sprintf(legendentry,"HX");
    	 plotLegend->AddEntry(h_lph_vs_lth_HX,legendentry,"p");
    	 sprintf(legendentry,"PX");
    	 plotLegend->AddEntry(h_lph_vs_lth_PX,legendentry,"p");
    	 plotLegend->Draw(); plotLegend->Draw();

         char savename[200];
         sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/lth_vs_lph_rap%d_pT%d.pdf",DataID,iState,irapBin,ipTbin);
         c1->SaveAs(savename);



/*      TCanvas *c2 = new TCanvas("c2","c2",800,800);
      gPad->SetFillColor(kWhite);
      h_lthlphCS->Draw();
      h_lthlphHX->Draw("same");
      h_lthlphPX->Draw("same");

*/

c2 = new TCanvas("c2","c2",1200,800);
SubtractedBG_test->SetXTitle("N_{BG,subtracted} / N_{BG,actual}");
gPad->SetFillColor(kWhite);
SubtractedBG_test->SetTitle(0);
SubtractedBG_test->Draw();

sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/SubtractedBackground_rap%d_pT%d.pdf",DataID,iState,irapBin,ipTbin);
c2->SaveAs(savename);

if(PlotAcceptance){
c2 = new TCanvas("c2","c2",1200,800);
SubtractedBG_test->SetXTitle("Metr.Hast. acceptance, CS");
gPad->SetFillColor(kWhite);
MetropolisHastingsAcceptanceCS->SetTitle(0);
MetropolisHastingsAcceptanceCS->Draw();

sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/MetropolisHastingsAcceptanceCS_rap%d_pT%d.pdf",DataID,iState,irapBin,ipTbin);
c2->SaveAs(savename);


c2 = new TCanvas("c2","c2",1200,800);
SubtractedBG_test->SetXTitle("Metr.Hast. acceptance, HX");
gPad->SetFillColor(kWhite);
MetropolisHastingsAcceptanceHX->SetTitle(0);
MetropolisHastingsAcceptanceHX->Draw();

sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/MetropolisHastingsAcceptanceHX_rap%d_pT%d.pdf",DataID,iState,irapBin,ipTbin);
c2->SaveAs(savename);


c2 = new TCanvas("c2","c2",1200,800);
SubtractedBG_test->SetXTitle("Metr.Hast. acceptance, PX");
gPad->SetFillColor(kWhite);
MetropolisHastingsAcceptancePX->SetTitle(0);
MetropolisHastingsAcceptancePX->Draw();

sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/MetropolisHastingsAcceptancePX_rap%d_pT%d.pdf",DataID,iState,irapBin,ipTbin);
c2->SaveAs(savename);
}
resultFile->Close();

			  }
		  }
	  }
		  return;
}
