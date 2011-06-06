#ifndef __CINT__
#include <RooGlobalFunc.h>
#endif
//#include "/home/zhlinl/cluster1/work/codes/myStyle.h"
#include "/afs/ihep.ac.cn/users/z/zhangll/workspace/work/codes/myStyle.h"
#include "commonVar.h"

//#include <TCanvas.h>
#include <TROOT.h>

#include <RooRealVar.h>
#include <RooDataSet.h>

#include <TBenchmark.h>
//#include <TAxis.h>
#include <TMath.h>
#include <TFile.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLine.h>
#include <TLatex.h>
#include <TIterator.h>
//#include <TVector3.h>
#include <vector.h>
//#include <TString.h>
#include <iostream.h>
#include <TGraphAsymmErrors.h>
using namespace std;

//void plotFinal(char *dataType="Run2010AB",char *scen="CS",char *sub="lambda_tilde"){
void plotFinal(char *dataType="Run2010AB",char *scen="CS",char sub[20]="lambda_phi",int temp=0){

	const int NRaps=2;
	char path[100];
	//char *test="lambda_tilde";
	char test[20];
	sprintf(test,"lambda_tilde");
	sprintf(path,Form("picLM/%s/%s/total",dataType,sub));
	gSystem->mkdir(path);
	//char dataType[10];
	//char scen[5];
	TGraphAsymmErrors *bFrac_pt[NRaps];
	TGraphAsymmErrors *lth_pt_p[NRaps];
	TGraphAsymmErrors *lphi_pt_p[NRaps];
	TGraphAsymmErrors *lthphi_pt_p[NRaps];
	TGraphAsymmErrors *lthtilde_pt_p[2][NRaps];
	TGraphAsymmErrors *f_pt_p[2][NRaps];
	TGraphAsymmErrors *lth_pt_np[NRaps];
	TGraphAsymmErrors *lphi_pt_np[NRaps];
	TGraphAsymmErrors *lthphi_pt_np[NRaps];
	TGraphAsymmErrors *lthtilde_pt_np[2][NRaps];
	TGraphAsymmErrors *f_pt_np[2][NRaps];
	//TCanvas *cv[NRaps][11];
	//TH1F *hg[NRaps][11];
	TCanvas *cv[11];
	TH1F *hg[11];
	char inName[200];

	sprintf(inName, Form("/afs/ihep.ac.cn/users/z/zhangll/workspace/work/polarization/Polarization/PolPt_%s_%s_%s.root",dataType,scen,sub)); cout<<"inName: "<<inName<<endl;
	TFile *infile=new TFile(inName,"R");
	if(!infile){
		cout<<">>===Error: no input File"<<endl;
		return;
	}
	cout<<">>===Open File successfully!"<<endl;

	if(test=="lambda_tilde") cout<<">>==debug002xx"<<endl;
	if(strcmp(sub,"lambda_tilde")) cout<<">>==debug002"<<endl;
	for(int rapBin=0;rapBin<NRaps;rapBin++){
		bFrac_pt[rapBin]=(TGraphAsymmErrors *)infile->Get(Form("bFrac_rap%d",rapBin+1));

		lth_pt_p[rapBin]=(TGraphAsymmErrors *)infile->Get(Form("lth_p_rap%d",rapBin+1));
		lphi_pt_p[rapBin]=(TGraphAsymmErrors *)infile->Get(Form("lphi_p_rap%d",rapBin+1));
		lthphi_pt_p[rapBin]=(TGraphAsymmErrors *)infile->Get(Form("lthphi_p_rap%d",rapBin+1));

		lth_pt_np[rapBin]=(TGraphAsymmErrors *)infile->Get(Form("lth_np_rap%d",rapBin+1));
		lphi_pt_np[rapBin]=(TGraphAsymmErrors *)infile->Get(Form("lphi_np_rap%d",rapBin+1));
		lthphi_pt_np[rapBin]=(TGraphAsymmErrors *)infile->Get(Form("lthphi_np_rap%d",rapBin+1));

		//if(sub=="lambda_tilde" && scen=="CS"){
		if(temp==1){
			cout<<">>==debug001"<<endl;
			lthtilde_pt_p[0][rapBin]=(TGraphAsymmErrors *)infile->Get(Form("lthtilde_p_rap%d",rapBin+1));
			lthtilde_pt_np[0][rapBin]=(TGraphAsymmErrors *)infile->Get(Form("lthtilde_np_rap%d",rapBin+1));
		}
		//if(sub=="F_invar" && scen=="CS"){
		if(temp==2){
			f_pt_p[0][rapBin]=(TGraphAsymmErrors *)infile->Get(Form("f_p_rap%d",rapBin+1));
			f_pt_np[0][rapBin]=(TGraphAsymmErrors *)infile->Get(Form("f_np_rap%d",rapBin+1));
		}

	}
	TLine *line0=new TLine(0.,0.,30.,0.);
	line0->SetLineStyle(9);
	TLine *line1=new TLine(0.,1./3.,30.,1./3.);
	line1->SetLineStyle(9);

	RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
	RooDataSet *data=new RooDataSet("data","data",JpsiMass);
	//TCanvas c1; TH1* legendRed =c1.DrawFrame(0.,0.,30.,0.8); legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(kSolid) ; legendRed->SetLineWidth(2) ; legendRed->SetMarkerStyle(20); legendRed->SetMarkerColor(kRed);
	TH1* legendRed = data->createHistogram("legendRed",JpsiMass) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(kSolid) ; legendRed->SetLineWidth(2) ; legendRed->SetMarkerStyle(20); legendRed->SetMarkerColor(kRed);
	TH1* legendBlack = data->createHistogram("legendBlack",JpsiMass) ; legendBlack->SetLineColor(kBlack); legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ; legendBlack->SetMarkerStyle(21); legendBlack->SetMarkerColor(kBlack);

	TLegend* legend=new TLegend(0.7,0.75,0.80,0.84);
	legend->SetFillColor(kWhite);
	legend->SetTextFont(42);  
	legend->SetTextSize(0.035);
	legend->SetBorderSize(0);
	legend->AddEntry(legendRed,"|y| < 0.9","p");
	legend->AddEntry(legendBlack,"0.9 < |y| < 1.2","p");

	double left=0.7, top=0.85, textSize=0.032;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	Double_t step=textSize*1.3;

	cv[0]=new TCanvas("cv_0","");
	hg[0]=cv[0]->DrawFrame(0.,0.,30.,0.8);
	hg[0]->GetYaxis()->SetTitle("B-fraction");
	hg[0]->GetXaxis()->SetTitle("p_{T} [GeV]");
	bFrac_pt[0]->Draw("P");
	bFrac_pt[1]->Draw("P");
	legend->Draw();

	cv[1]=new TCanvas("cv_1","");
	hg[1]=cv[1]->DrawFrame(0.,-1.,30.,1.);
	hg[1]->GetYaxis()->SetTitle(Form("#lambda^{%s}_{#theta}",scen));
	hg[1]->GetXaxis()->SetTitle("p_{T} [GeV]");
	lth_pt_p[0]->Draw("P");
	lth_pt_p[1]->Draw("P");
	line0->Draw();
	legend->Draw();
	latex->DrawLatex(left,top,"Prompt");

	cv[2]=new TCanvas("cv_2","");
	hg[2]=cv[2]->DrawFrame(0.,-1.,30.,1.);
	hg[2]->GetYaxis()->SetTitle(Form("#lambda^{%s}_{#phi}",scen));
	hg[2]->GetXaxis()->SetTitle("p_{T} [GeV]");
	lphi_pt_p[0]->Draw("P");
	lphi_pt_p[1]->Draw("P");
	line0->Draw();
	legend->Draw();
	latex->DrawLatex(left,top,"Prompt");

	cv[3]=new TCanvas("cv_3","");
	hg[3]=cv[3]->DrawFrame(0.,-1.,30.,1.);
	hg[3]->GetYaxis()->SetTitle(Form("#lambda^{%s}_{#theta#phi}",scen));
	hg[3]->GetXaxis()->SetTitle("p_{T} [GeV]");
	lthphi_pt_p[0]->Draw("P");
	lthphi_pt_p[1]->Draw("P");
	line0->Draw();
	legend->Draw();
	latex->DrawLatex(left,top,"Prompt");

	cv[4]=new TCanvas("cv_4","");
	hg[4]=cv[4]->DrawFrame(0.,-1.,30.,1.);
	hg[4]->GetYaxis()->SetTitle(Form("#lambda^{%s}_{#theta}",scen));
	hg[4]->GetXaxis()->SetTitle("p_{T} [GeV]");
	lth_pt_np[0]->Draw("P");
	lth_pt_np[1]->Draw("P");
	line0->Draw();
	legend->Draw();
	latex->DrawLatex(left,top,"non-Prompt");

	cv[5]=new TCanvas("cv_5","");
	hg[5]=cv[5]->DrawFrame(0.,-1.,30.,1.);
	hg[5]->GetYaxis()->SetTitle(Form("#lambda^{%s}_{#phi}",scen));
	hg[5]->GetXaxis()->SetTitle("p_{T} [GeV]");
	lphi_pt_np[0]->Draw("P");
	lphi_pt_np[1]->Draw("P");
	line0->Draw();
	legend->Draw();
	latex->DrawLatex(left,top,"non-Prompt");

	cv[6]=new TCanvas("cv_6","");
	hg[6]=cv[6]->DrawFrame(0.,-1.,30.,1.);
	hg[6]->GetYaxis()->SetTitle(Form("#lambda^{%s}_{#theta#phi}",scen));
	hg[6]->GetXaxis()->SetTitle("p_{T} [GeV]");
	lthphi_pt_np[0]->Draw("P");
	lthphi_pt_np[1]->Draw("P");
	line0->Draw();
	legend->Draw();
	latex->DrawLatex(left,top,"non-Prompt");

	cv[0]->SaveAs(Form("%s/bFrac_pt-%s.png",path,scen));
	cv[1]->SaveAs(Form("%s/lth_pt_p-%s.png",path,scen));
	cv[2]->SaveAs(Form("%s/lphi_pt_p-%s.png",path,scen));
	cv[3]->SaveAs(Form("%s/lthphi_pt_p-%s.png",path,scen));
	cv[4]->SaveAs(Form("%s/lth_pt_np-%s.png",path,scen));
	cv[5]->SaveAs(Form("%s/lphi_pt_np-%s.png",path,scen));
	cv[6]->SaveAs(Form("%s/lthphi_pt_np-%s.png",path,scen));
	cv[0]->SaveAs(Form("%s/bFrac_pt-%s.eps",path,scen));
	cv[1]->SaveAs(Form("%s/lth_pt_p-%s.eps",path,scen));
	cv[2]->SaveAs(Form("%s/lphi_pt_p-%s.eps",path,scen));
	cv[3]->SaveAs(Form("%s/lthphi_pt_p-%s.eps",path,scen));
	cv[4]->SaveAs(Form("%s/lth_pt_np-%s.eps",path,scen));
	cv[5]->SaveAs(Form("%s/lphi_pt_np-%s.eps",path,scen));
	cv[6]->SaveAs(Form("%s/lthphi_pt_np-%s.eps",path,scen));


	left=0.13;
	//if(sub=="lambda_tilde" && scen=="CS"){
	if(temp==1){
		sprintf(inName, Form("/afs/ihep.ac.cn/users/z/zhangll/workspace/work/polarization/Polarization/PolPt_%s_HX_%s.root",dataType,sub)); cout<<"inName: "<<inName<<endl;
		TFile *inFile=new TFile(inName,"R");
		if(!inFile){
			cout<<">>===Error: no input File"<<endl;
			return;
		}
		cout<<"<<===Open File successfully!"<<endl;
		for(int rapBin=0;rapBin<NRaps;rapBin++){
			lthtilde_pt_p[1][rapBin]=(TGraphAsymmErrors *)inFile->Get(Form("lthtilde_p_rap%d",rapBin+1));
			lthtilde_pt_np[1][rapBin]=(TGraphAsymmErrors *)inFile->Get(Form("lthtilde_np_rap%d",rapBin+1));
		}
		TLegend* legend0=new TLegend(left,0.75,0.80,0.84);
		legend0->SetFillColor(kWhite);
		legend0->SetTextFont(42);  
		legend0->SetTextSize(0.035);
		legend0->SetBorderSize(0);
		legend0->AddEntry(legendRed,"CS","p");
		legend0->AddEntry(legendBlack,"HX","p");

		cv[7]=new TCanvas("cv_7","");
		hg[7]=cv[7]->DrawFrame(0.,-1.,30.,1.);
		hg[7]->GetYaxis()->SetTitle("#lambdatilde");
		hg[7]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lthtilde_pt_p[0][0]->SetMarkerStyle(20);
		lthtilde_pt_p[0][0]->SetMarkerColor(kRed);
		lthtilde_pt_p[1][0]->SetMarkerStyle(21);
		lthtilde_pt_p[1][0]->SetMarkerColor(kBlack);
		lthtilde_pt_p[0][0]->Draw("P");
		lthtilde_pt_p[1][0]->Draw("P");
		line0->Draw();
		legend0->Draw();
		latex->DrawLatex(left,top,"Prompt J/#psi, |y| < 0.9");
		cv[7]->SaveAs(Form("%s/lthtilde_pt_p-rap1.png",path));
		cv[7]->SaveAs(Form("%s/lthtilde_pt_p-rap1.eps",path));

		cv[8]=new TCanvas("cv_8","");
		hg[8]=cv[8]->DrawFrame(0.,-1.,30.,1.);
		hg[8]->GetYaxis()->SetTitle("#lambdatilde");
		hg[8]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lthtilde_pt_p[0][1]->SetMarkerStyle(20);
		lthtilde_pt_p[0][1]->SetMarkerColor(kRed);
		lthtilde_pt_p[1][1]->SetMarkerStyle(21);
		lthtilde_pt_p[1][1]->SetMarkerColor(kBlack);
		lthtilde_pt_p[0][1]->Draw("P");
		lthtilde_pt_p[1][1]->Draw("P");
		line0->Draw();
		legend0->Draw();
		latex->DrawLatex(left,top,"Prompt J/#psi, 0.9 < |y| < 1.2");
		cv[8]->SaveAs(Form("%s/lthtilde_pt_p-rap2.png",path));
		cv[8]->SaveAs(Form("%s/lthtilde_pt_p-rap2.eps",path));

		cv[9]=new TCanvas("cv_9","");
		hg[9]=cv[9]->DrawFrame(0.,-1.,30.,1.);
		hg[9]->GetYaxis()->SetTitle("#lambdatilde");
		hg[9]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lthtilde_pt_np[0][0]->SetMarkerStyle(20);
		lthtilde_pt_np[0][0]->SetMarkerColor(kRed);
		lthtilde_pt_np[1][0]->SetMarkerStyle(21);
		lthtilde_pt_np[1][0]->SetMarkerColor(kBlack);
		lthtilde_pt_np[0][0]->Draw("P");
		lthtilde_pt_np[1][0]->Draw("P");
		line0->Draw();
		legend0->Draw();
		latex->DrawLatex(left,top,"non-Prompt J/#psi, |y| < 0.9");
		cv[9]->SaveAs(Form("%s/lthtilde_pt_np-rap1.png",path));
		cv[9]->SaveAs(Form("%s/lthtilde_pt_np-rap1.eps",path));

		cv[10]=new TCanvas("cv_10","");
		hg[10]=cv[10]->DrawFrame(0.,-1.,30.,1.);
		hg[10]->GetYaxis()->SetTitle("#lambdatilde");
		hg[10]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lthtilde_pt_np[0][1]->SetMarkerStyle(20);
		lthtilde_pt_np[0][1]->SetMarkerColor(kRed);
		lthtilde_pt_np[1][1]->SetMarkerStyle(21);
		lthtilde_pt_np[1][1]->SetMarkerColor(kBlack);
		lthtilde_pt_np[0][1]->Draw("P");
		lthtilde_pt_np[1][1]->Draw("P");
		line0->Draw();
		legend0->Draw();
		latex->DrawLatex(left,top,"non-Prompt J/#psi, 0.9 < |y| < 1.2");
		cv[10]->SaveAs(Form("%s/lthtilde_pt_np-rap2.png",path));
		cv[10]->SaveAs(Form("%s/lthtilde_pt_np-rap2.eps",path));
	}


	//if(sub=="F_invar" && scen=="CS"){
	if(temp==2){
		sprintf(inName, Form("/afs/ihep.ac.cn/users/z/zhangll/workspace/work/polarization/Polarization/PolPt_%s_HX_%s.root",dataType,sub)); cout<<"inName: "<<inName<<endl;
		TFile *inFile=new TFile(inName,"R");
		if(!inFile){
			cout<<">>===Error: no input File"<<endl;
			return;
		}
		cout<<"<<===Open File successfully!"<<endl;
		for(int rapBin=0;rapBin<NRaps;rapBin++){
			f_pt_p[1][rapBin]=(TGraphAsymmErrors *)inFile->Get(Form("f_p_rap%d",rapBin+1));
			f_pt_np[1][rapBin]=(TGraphAsymmErrors *)inFile->Get(Form("f_p_rap%d",rapBin+1));
		}

		TLegend* legend1=new TLegend(left,0.75,0.80,0.84);
		legend1->SetFillColor(kWhite);
		legend1->SetTextFont(42);  
		legend1->SetTextSize(0.035);
		legend1->SetBorderSize(0);
		legend1->AddEntry(legendRed,"CS","p");
		legend1->AddEntry(legendBlack,"HX","p");

		cv[7]=new TCanvas("cv_7","");
		hg[7]=cv[7]->DrawFrame(0.,0.,30.,1.);
		hg[7]->GetYaxis()->SetTitle("F");
		hg[7]->GetXaxis()->SetTitle("p_{T} [GeV]");
		f_pt_p[0][0]->SetMarkerStyle(20);
		f_pt_p[0][0]->SetMarkerColor(kRed);
		f_pt_p[1][0]->SetMarkerStyle(21);
		f_pt_p[1][0]->SetMarkerColor(kBlack);
		f_pt_p[0][0]->Draw("P");
		f_pt_p[1][0]->Draw("P");
		line1->Draw();
		legend1->Draw();
		latex->DrawLatex(left,top,"Prompt J/#psi, |y| < 0.9");
		cv[7]->SaveAs(Form("%s/f_pt_p-rap1.png",path));
		cv[7]->SaveAs(Form("%s/f_pt_p-rap1.eps",path));

		cv[8]=new TCanvas("cv_8","");
		hg[8]=cv[8]->DrawFrame(0.,0.,30.,1.);
		hg[8]->GetYaxis()->SetTitle("F");
		hg[8]->GetXaxis()->SetTitle("p_{T} [GeV]");
		f_pt_p[0][1]->SetMarkerStyle(20);
		f_pt_p[0][1]->SetMarkerColor(kRed);
		f_pt_p[1][1]->SetMarkerStyle(21);
		f_pt_p[1][1]->SetMarkerColor(kBlack);
		f_pt_p[0][1]->Draw("P");
		f_pt_p[1][1]->Draw("P");
		line1->Draw();
		legend1->Draw();
		latex->DrawLatex(left,top,"Prompt J/#psi, 0.9 < |y| < 1.2");
		cv[8]->SaveAs(Form("%s/f_pt_p-rap2.png",path));
		cv[8]->SaveAs(Form("%s/f_pt_p-rap2.eps",path));

		cv[9]=new TCanvas("cv_9","");
		hg[9]=cv[9]->DrawFrame(0.,0.,30.,1.);
		hg[9]->GetYaxis()->SetTitle("F");
		hg[9]->GetXaxis()->SetTitle("p_{T} [GeV]");
		f_pt_np[0][0]->SetMarkerStyle(20);
		f_pt_np[0][0]->SetMarkerColor(kRed);
		f_pt_np[1][0]->SetMarkerStyle(21);
		f_pt_np[1][0]->SetMarkerColor(kBlack);
		f_pt_np[0][0]->Draw("P");
		f_pt_np[1][0]->Draw("P");
		line1->Draw();
		legend1->Draw();
		latex->DrawLatex(left,top,"non-Prompt J/#psi, |y| < 0.9");
		cv[9]->SaveAs(Form("%s/f_pt_np-rap1.png",path));
		cv[9]->SaveAs(Form("%s/f_pt_np-rap1.eps",path));

		cv[10]=new TCanvas("cv_10","");
		hg[10]=cv[10]->DrawFrame(0.,0.,30.,1.);
		hg[10]->GetYaxis()->SetTitle("F");
		hg[10]->GetXaxis()->SetTitle("p_{T} [GeV]");
		f_pt_np[0][1]->SetMarkerStyle(20);
		f_pt_np[0][1]->SetMarkerColor(kRed);
		f_pt_np[1][1]->SetMarkerStyle(21);
		f_pt_np[1][1]->SetMarkerColor(kBlack);
		f_pt_np[0][1]->Draw("P");
		f_pt_np[1][1]->Draw("P");
		line1->Draw();
		legend1->Draw();
		latex->DrawLatex(left,top,"non-Prompt J/#psi, 0.9 < |y| < 1.2");
		cv[10]->SaveAs(Form("%s/f_pt_np-rap2.png",path));
		cv[10]->SaveAs(Form("%s/f_pt_np-rap2.eps",path));
	}

}

