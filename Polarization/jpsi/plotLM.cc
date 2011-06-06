#ifndef __CINT__
#include <RooGlobalFunc.h>
#endif
//#include "/home/zhlinl/cluster1/work/codes/myStyle.h"
#include "/afs/ihep.ac.cn/users/z/zhangll/workspace/work/codes/myStyle.h"
#include "commonVar.h"

//#include <TLatex.h>
#include <RooRealVar.h> 
#include <RooDataSet.h>
#include <RooGaussian.h>
//#include <RooFitResult.h>
//#include <RooLandau.h>
//#include <RooChebychev.h>
//#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooDataHist.h>
#include <RooVoigtian.h>
#include <RooBreitWigner.h>
#include <RooKeysPdf.h>
#include <RooFFTConvPdf.h>
//#include <RooCBShape.h>
//#include <TCanvas.h>
#include <TROOT.h>
#include <RooRealVar.h> 
#include <RooCategory.h> 
#include <RooSimultaneous.h>
#include <RooWorkspace.h>
#include <RooArgSet.h>

#include <TBenchmark.h>
//#include <TAxis.h>
//#include <TH1.h>
//#include <TLine.h>
#include <TMath.h>
#include <TFile.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TH1.h>
#include <TLatex.h>
#include <TIterator.h>
//#include <TVector3.h>
#include <vector.h>
//#include <TString.h>
#include <iostream.h>
using namespace RooFit;

void plotLM(int rapBin=2, int ptBin=8) {

	char path[100];
	char scen[4]="CS";
	char dataType[10]="Run2010AB";
	int nbins=50;
	sprintf(path,Form("pic/%s",dataType));

	gSystem->mkdir(path,kTRUE);

	char inName[200];
	sprintf(inName, Form("/afs/ihep.ac.cn/users/z/zhangll/workspace/work/polarization/Polarization/project/%s/%s/pola%s_rap%d_pt%d.root",dataType,scen,dataType,rapBin,ptBin)); cout<<"inName: "<<inName<<endl;
	TFile *inFile=new TFile(inName,"R");
	if(!inFile) return;
	RooWorkspace *ws=(RooWorkspace *)inFile->Get(Form("pola%s_rap%d_pt%d",dataType,rapBin,ptBin));
	if(!ws) 
	{
		cout<<">>=======Error: no workspace in root file=========="<<endl;
		return;
	}
	//RooRealVar *JpsiMass=(RooRealVar*)ws->var("JpsiMass");
	RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);

	RooPlot *ctauFrameSig=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins));
	ctauFrameSig->SetName(Form("ctausig_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrameSig->SetTitle("");
	RooPlot *ctauFrameBkg=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins));
	ctauFrameBkg->SetName(Form("ctaubkg_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrameBkg->SetTitle("");
	RooPlot *massFrame=((RooRealVar*)ws->var("JpsiMass"))->frame(Bins(nbins));
	massFrame->SetName(Form("mass_plot_rap%d_pt%d",rapBin,ptBin));
	massFrame->SetTitle("");
	//massFrame->SetYTitle("");
	massFrame->GetYaxis()->SetTitle("");

	RooDataSet *data=(RooDataSet *)ws->data(Form("data_rap%d_pt%d",rapBin,ptBin));
	
	RooArgSet *dataVars=(RooArgSet *)data->get();
	dataVars->Print();
	TIterator* di=(TIterator*)dataVars->createIterator();

	RooRealVar *nNP_=(RooRealVar *)ws->var("nNonPromptSignal");
	RooRealVar *nP_=(RooRealVar *)ws->var("nPromptSignal");
	double nNP=nNP_->getVal();
	double nP=nP_->getVal();
	cout<<">>==nNP: "<<nNP<<endl;
	cout<<">>==nP: "<<nP<<endl;

	//return;

	RooAbsData *dsig = data->reduce("massRegion == massRegion::signal");
	RooAbsData *dpr = data->reduce("mlRegion == mlRegion::promptSignal");
	RooAbsData *dnp = data->reduce("mlRegion == mlRegion::nonPromptSignal");
	RooAbsData *dbkg = data->reduce("massRegion == massRegion::rightMassSideBand || massRegion == massRegion::leftMassSideBand");
	RooAbsData *dbkgl = data->reduce("massRegion == massRegion::leftMassSideBand");
	RooAbsData *dbkgr = data->reduce("massRegion == massRegion::rightMassSideBand");

	RooSimultaneous *mpdf=(RooSimultaneous *)ws->pdf("MPdf");
	RooSimultaneous *lpdf=(RooSimultaneous *)ws->pdf("LPdf");

	dsig->plotOn(ctauFrameSig,MarkerSize(0.8),MarkerStyle(20));
	lpdf->plotOn(ctauFrameSig,
			LineWidth(2),
			ProjWData(*dsig));

	lpdf->plotOn(ctauFrameSig,
			Components("promptExtCTau"),
			LineStyle(5),
			LineColor(kBlue),
			LineWidth(2),
			ProjWData(*dsig));

	lpdf->plotOn(ctauFrameSig,
			Components("nonpromptExtCTau"),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2),
			ProjWData(*dsig));

	lpdf->plotOn(ctauFrameSig,
			Components("backgroundExtCTau"),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2),
			ProjWData(*dsig));               

	dbkg->plotOn(ctauFrameBkg,MarkerSize(0.8));
	lpdf->plotOn(ctauFrameBkg,
			LineWidth(2),
			ProjWData(*dbkg));


	data->plotOn(massFrame,MarkerSize(0.8));
	mpdf->plotOn(massFrame,
			LineWidth(2),
			ProjWData(*data));

	mpdf->plotOn(massFrame,
			Components("backgroundExtMass"),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2),
			ProjWData(*data));
	mpdf->plotOn(massFrame,
			Components("nonpromptExtMass"),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2),
			ProjWData(*data));

	mpdf->plotOn(massFrame,
			Components("promptExtMass"),
			LineStyle(5),
			LineColor(kBlue),
			LineWidth(2),
			ProjWData(*data));

	TH1* legendBlue = data->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = data->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = data->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = data->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = data->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
	TH1* legendPink = data->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	//TLegend* LifetimeLegend=new TLegend(0.6,0.6,0.9,0.9);
	TLegend* LifetimeLegend=new TLegend(0.65,0.7,0.88,0.88);
	LifetimeLegend->SetFillColor(kWhite);
	//LifetimeLegend->SetTextFont(72);  
	LifetimeLegend->SetTextFont(42);  
	LifetimeLegend->SetTextSize(0.025);
	LifetimeLegend->SetBorderSize(0.);
	LifetimeLegend->AddEntry(legendBlue,"sum","l"); 
 	LifetimeLegend->AddEntry(legendBlueDash,"prompt J/#psi","l");
	LifetimeLegend->AddEntry(legendRed,"non-Prompt J/#psi","l");
	LifetimeLegend->AddEntry(legendPink,"background","l");

	TLegend* MassLegend=new TLegend(0.65,0.7,0.88,0.88);
	MassLegend->SetFillColor(kWhite);
	MassLegend->SetTextFont(42);  
	MassLegend->SetTextSize(0.025);
	//MassLegend->SetLineWidth(0.);
	MassLegend->SetBorderSize(0.);
	MassLegend->AddEntry(legendBlue,"sum","l"); 
 	MassLegend->AddEntry(legendBlueDash,"prompt J/#psi","l");
	MassLegend->AddEntry(legendRed,"non-Prompt J/#psi","l");
	MassLegend->AddEntry(legendPink,"background","l");

	double left=0.65, top=0.65, textSize=0.032;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	Double_t step=textSize*1.3;


	//SetMyStyle();
	TCanvas *c1=new TCanvas("c1","");
	c1->SetLogy();
	ctauFrameSig->Draw();
	LifetimeLegend->Draw();
	if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",jpsi::rapForPTRange[rapBin]));
	else 
		latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",jpsi::rapForPTRange[rapBin-1],jpsi::rapForPTRange[rapBin]));
	top-=step;
	latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",jpsi::pTRange[rapBin][ptBin-1],jpsi::pTRange[rapBin][ptBin]));

	TCanvas *c2=new TCanvas("c2","");
	c2->SetLogy();
	ctauFrameBkg->Draw();
	top=0.75;
	if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",jpsi::rapForPTRange[rapBin]));
	else
	 	latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",jpsi::rapForPTRange[rapBin-1],jpsi::rapForPTRange[rapBin]));
	top-=step;
	latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",jpsi::pTRange[rapBin][ptBin-1],jpsi::pTRange[rapBin][ptBin]));

	TCanvas *c3=new TCanvas("c3","");
	massFrame->Draw();
	MassLegend->Draw();
	top=0.65;
	if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",jpsi::rapForPTRange[rapBin]));
	else 
		latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",jpsi::rapForPTRange[rapBin-1],jpsi::rapForPTRange[rapBin]));
	top-=step;
	latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",jpsi::pTRange[rapBin][ptBin-1],jpsi::pTRange[rapBin][ptBin]));
	//c1->SetLogy(0);
	c1->SaveAs(Form("%s/lifeTimeSig_rap%d_pT%d_%s.png",path,rapBin,ptBin,scen));
	c2->SaveAs(Form("%s/lifeTimeBkg_rap%d_pT%d_%s.png",path,rapBin,ptBin,scen));
	c3->SaveAs(Form("%s/mass_rap%d_pT%d_%s.png",path,rapBin,ptBin,scen));
}
