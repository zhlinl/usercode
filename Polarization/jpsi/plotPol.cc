#ifndef __CINT__
#include <RooGlobalFunc.h>
#endif
//#include "/home/zhlinl/cluster1/work/codes/myStyle.h"
#include "/afs/ihep.ac.cn/users/z/zhangll/workspace/work/codes/myStyle.h"
#include "commonVar.h"

//#include <TLatex.h>
//#include <TCanvas.h>
#include <TROOT.h>

#include <TBenchmark.h>
//#include <TAxis.h>
//#include <TLine.h>
#include <TMath.h>
#include <TFile.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TIterator.h>
//#include <TVector3.h>
#include <vector.h>
//#include <TString.h>
#include <iostream.h>
#include <TGraphAsymmErrors.h>

void plotPol(char *dataType="Run2010AB",char *scen="CS",char *sub="lambda_phi"){

	//double nNP[2][jpsi::kNbPTMaxBins];
	//double nP[2][jpsi::kNbPTMaxBins];
	const int NRaps=2;
	int total_nonzero_points=0;
	char path[100];
	//char dataType[10];
	//char scen[5];
	int nbins=50;
	int NbinsX=0;
	double x=0,y=0;
	double exh=0,exl=0,eyh=0,eyl=0;
	double lamdaTheta=0,lamdaPhi=0,lThtilde=0,f=0;
	//sprintf(scen,"HX");
	//sprintf(dataType,"Run2010A");
	sprintf(path,Form("pic/%s/%s",dataType,sub));
	gSystem->mkdir(path,kTRUE);
	TH1F *hist_bFrac_pt[NRaps];
	TGraphAsymmErrors *graphThis;
	TCanvas *cv[NRaps][11];
	TH1F *hg[NRaps][11];
	//TFile *outFile[NRaps];
	TFile *outFile=new TFile(Form("PolPt_%s_%s_%s.root",dataType,scen,sub),"recreate");

	TGraphAsymmErrors *bFrac_pt[NRaps];
	TGraphAsymmErrors *lth_pt_p[NRaps];
	TGraphAsymmErrors *lphi_pt_p[NRaps];
	TGraphAsymmErrors *lthphi_pt_p[NRaps];
	TGraphAsymmErrors *lthtilde_pt_p[NRaps];
	TGraphAsymmErrors *f_pt_p[NRaps];
	TGraphAsymmErrors *lth_pt_np[NRaps];
	TGraphAsymmErrors *lphi_pt_np[NRaps];
	TGraphAsymmErrors *lthphi_pt_np[NRaps];
	TGraphAsymmErrors *lthtilde_pt_np[NRaps];
	TGraphAsymmErrors *f_pt_np[NRaps];

	char inName[200];
	for(int rapBin=0;rapBin<NRaps;rapBin++){

		total_nonzero_points=0;
		bFrac_pt[rapBin]=new TGraphAsymmErrors();
		lth_pt_p[rapBin]=new TGraphAsymmErrors();
		lphi_pt_p[rapBin]=new TGraphAsymmErrors();
		lthphi_pt_p[rapBin]=new TGraphAsymmErrors();
		lthtilde_pt_p[rapBin]=new TGraphAsymmErrors();
		f_pt_p[rapBin]=new TGraphAsymmErrors();
		lth_pt_np[rapBin]=new TGraphAsymmErrors();
		lphi_pt_np[rapBin]=new TGraphAsymmErrors();
		lthphi_pt_np[rapBin]=new TGraphAsymmErrors();
		lthtilde_pt_np[rapBin]=new TGraphAsymmErrors();
		f_pt_np[rapBin]=new TGraphAsymmErrors();
		//hist_bFrac_pt[rapBin]= new TH1F(Form("hist_bFrac_pt%d",rapBin),"",100,0,35);

		for(int ptBin=5;ptBin<jpsi::kNbPTBins[rapBin+1];ptBin++){ 
			sprintf(inName, Form("/afs/ihep.ac.cn/users/z/zhangll/workspace/work/polarization/Polarization/project/%s/%s/%s/pola%s_rap%d_pt%d-%s_test.root",dataType,scen,sub,dataType,rapBin+1,ptBin+1,scen)); cout<<"inName: "<<inName<<endl;
			TFile *inFile=new TFile(inName,"R");
			if(!inFile){
				cout<<">>===Error: no input File"<<endl;
				continue;
			}
			//gDirectory->ls();
			graphThis=(TGraphAsymmErrors *)inFile->Get(Form("bFrac_rap%d",rapBin+1));
			if(!graphThis){
				cout<<">>===Error: no object found!"<<endl;
				continue;
			}
			lamdaTheta=0;lamdaPhi=0;lThtilde=0,f=0;

			NbinsX=graphThis->GetN();
			//TMath::LocMax to locate the maximum element of an array 
			int point=TMath::LocMax(graphThis->GetN(),graphThis->GetY());
			// //the corresponding x,y value can be obtained with
			x=graphThis->GetX()[point];
			y=graphThis->GetY()[point];
			if(x!=0 && y!=0){
				cout<<">>x: "<<x<<" y: "<<y<<endl;

				bFrac_pt[rapBin]->SetName(Form("bFrac_rap%d",rapBin+1));
				x=graphThis->GetX()[point];
				y=graphThis->GetY()[point];
				bFrac_pt[rapBin]->SetPoint(total_nonzero_points,x,y);
				//exh=graphThis->GetErrorXhigh(point);
				//exl=graphThis->GetErrorXlow(point);
				exh=0.;exl=0.;
				eyh=graphThis->GetErrorYhigh(point);
				eyl=graphThis->GetErrorYlow(point);
				bFrac_pt[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);

				//prompt
				graphThis=(TGraphAsymmErrors *)inFile->Get(Form("lth_p_rap%d",rapBin+1));
				x=graphThis->GetX()[point];
				y=graphThis->GetY()[point];
				lth_pt_p[rapBin]->SetName(Form("lth_p_rap%d",rapBin+1));
				lth_pt_p[rapBin]->SetPoint(total_nonzero_points,x,y);
				exh=graphThis->GetErrorXhigh(point);
				exl=graphThis->GetErrorXlow(point);
				eyh=graphThis->GetErrorYhigh(point);
				eyl=graphThis->GetErrorYlow(point);
				lth_pt_p[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);
				lamdaTheta=y;

				graphThis=(TGraphAsymmErrors *)inFile->Get(Form("lthphi_p_rap%d",rapBin+1));
				x=graphThis->GetX()[point];
				y=graphThis->GetY()[point];
				lthphi_pt_p[rapBin]->SetName(Form("lthphi_p_rap%d",rapBin+1));
				lthphi_pt_p[rapBin]->SetPoint(total_nonzero_points,x,y);
				exh=graphThis->GetErrorXhigh(point);
				exl=graphThis->GetErrorXlow(point);
				eyh=graphThis->GetErrorYhigh(point);
				eyl=graphThis->GetErrorYlow(point);
				lthphi_pt_p[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);


				graphThis=(TGraphAsymmErrors *)inFile->Get(Form("lthtilde_p_rap%d",rapBin+1));
				x=graphThis->GetX()[point];
				y=graphThis->GetY()[point];
				lthtilde_pt_p[rapBin]->SetName(Form("lthtilde_p_rap%d",rapBin+1));
				lthtilde_pt_p[rapBin]->SetPoint(total_nonzero_points,x,y);
				exh=graphThis->GetErrorXhigh(point);
				exl=graphThis->GetErrorXlow(point);
				eyh=graphThis->GetErrorYhigh(point);
				eyl=graphThis->GetErrorYlow(point);
				lthtilde_pt_p[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);
				lThtilde=y;

				if(lThtilde!=-3)	lamdaPhi=(lThtilde-lamdaTheta)/(3.+lThtilde);
				graphThis=(TGraphAsymmErrors *)inFile->Get(Form("f_p_rap%d",rapBin+1));
				graphThis=(TGraphAsymmErrors *)inFile->Get(Form("lphi_p_rap%d",rapBin+1));
				x=graphThis->GetX()[point];
				y=graphThis->GetY()[point];
				lphi_pt_p[rapBin]->SetName(Form("lphi_p_rap%d",rapBin+1));
				lphi_pt_p[rapBin]->SetPoint(total_nonzero_points,x,y);
				//lphi_pt_p[rapBin]->SetPoint(total_nonzero_points,x,lamdaPhi);
				exh=graphThis->GetErrorXhigh(point);
				exl=graphThis->GetErrorXlow(point);
				eyh=graphThis->GetErrorYhigh(point);
				eyl=graphThis->GetErrorYlow(point);
				lphi_pt_p[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);

				if(lamdaTheta!=-3) f=(1.+lamdaTheta+2.*lamdaPhi)/(3.+lamdaTheta);
				graphThis=(TGraphAsymmErrors *)inFile->Get(Form("f_p_rap%d",rapBin+1));
				x=graphThis->GetX()[point];
				y=graphThis->GetY()[point];
				f_pt_p[rapBin]->SetName(Form("f_p_rap%d",rapBin+1));
				f_pt_p[rapBin]->SetPoint(total_nonzero_points,x,y);
				//f_pt_p[rapBin]->SetPoint(total_nonzero_points,x,f);
				exh=graphThis->GetErrorXhigh(point);
				exl=graphThis->GetErrorXlow(point);
				eyh=graphThis->GetErrorYhigh(point);
				eyl=graphThis->GetErrorYlow(point);
				f_pt_p[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);

				//non-prompt
				lamdaTheta=0;lamdaPhi=0;lThtilde=0;f=0;
				graphThis=(TGraphAsymmErrors *)inFile->Get(Form("lth_np_rap%d",rapBin+1));
				x=graphThis->GetX()[point];
				y=graphThis->GetY()[point];
				lth_pt_np[rapBin]->SetName(Form("lth_np_rap%d",rapBin+1));
				lth_pt_np[rapBin]->SetPoint(total_nonzero_points,x,y);
				exh=graphThis->GetErrorXhigh(point);
				exl=graphThis->GetErrorXlow(point);
				eyh=graphThis->GetErrorYhigh(point);
				eyl=graphThis->GetErrorYlow(point);
				lth_pt_np[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);
				lamdaTheta=y;

				graphThis=(TGraphAsymmErrors *)inFile->Get(Form("lthphi_np_rap%d",rapBin+1));
				x=graphThis->GetX()[point];
				y=graphThis->GetY()[point];
				lthphi_pt_np[rapBin]->SetName(Form("lthphi_np_rap%d",rapBin+1));
				lthphi_pt_np[rapBin]->SetPoint(total_nonzero_points,x,y);
				exh=graphThis->GetErrorXhigh(point);
				exl=graphThis->GetErrorXlow(point);
				eyh=graphThis->GetErrorYhigh(point);
				eyl=graphThis->GetErrorYlow(point);
				lthphi_pt_np[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);

				graphThis=(TGraphAsymmErrors *)inFile->Get(Form("lthtilde_np_rap%d",rapBin+1));
				x=graphThis->GetX()[point];
				y=graphThis->GetY()[point];
				lthtilde_pt_np[rapBin]->SetName(Form("lthtilde_np_rap%d",rapBin+1));
				lthtilde_pt_np[rapBin]->SetPoint(total_nonzero_points,x,y);
				exh=graphThis->GetErrorXhigh(point);
				exl=graphThis->GetErrorXlow(point);
				eyh=graphThis->GetErrorYhigh(point);
				eyl=graphThis->GetErrorYlow(point);
				lthtilde_pt_np[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);
				lThtilde=y;

				if(lThtilde!=-3)	lamdaPhi=(lThtilde-lamdaTheta)/(3.+lThtilde);
				graphThis=(TGraphAsymmErrors *)inFile->Get(Form("lphi_np_rap%d",rapBin+1));
				x=graphThis->GetX()[point];
				y=graphThis->GetY()[point];
				lphi_pt_np[rapBin]->SetName(Form("lphi_np_rap%d",rapBin+1));
				lphi_pt_np[rapBin]->SetPoint(total_nonzero_points,x,y);
				//lphi_pt_np[rapBin]->SetPoint(total_nonzero_points,x,lamdaPhi);
				exh=graphThis->GetErrorXhigh(point);
				exl=graphThis->GetErrorXlow(point);
				eyh=graphThis->GetErrorYhigh(point);
				eyl=graphThis->GetErrorYlow(point);
				lphi_pt_np[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);

				if(lamdaTheta!=-3) f=(1.+lamdaTheta+2.*lamdaPhi)/(3.+lamdaTheta);
				graphThis=(TGraphAsymmErrors *)inFile->Get(Form("f_np_rap%d",rapBin+1));
				x=graphThis->GetX()[point];
				y=graphThis->GetY()[point];
				f_pt_np[rapBin]->SetName(Form("f_np_rap%d",rapBin+1));
				f_pt_np[rapBin]->SetPoint(total_nonzero_points,x,y);
				//f_pt_np[rapBin]->SetPoint(total_nonzero_points,x,f);
				exh=graphThis->GetErrorXhigh(point);
				exl=graphThis->GetErrorXlow(point);
				eyh=graphThis->GetErrorYhigh(point);
				eyl=graphThis->GetErrorYlow(point);
				f_pt_np[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);

				total_nonzero_points++;
			}
	//		for(int point=0;point<NbinsX;point++){
	//			graphThis->GetPoint(point,x,y);
	//			if(x!=0 && y!=0){

	//				cout<<">>x: "<<x<<" y: "<<y<<endl;
	//				bFrac_pt[rapBin]->SetPoint(total_nonzero_points,x,y);
	//				exh=graphThis->GetErrorXhigh(point);
	//				exl=graphThis->GetErrorXlow(point);
	//				eyh=graphThis->GetErrorYhigh(point);
	//				eyl=graphThis->GetErrorYlow(point);
	//				bFrac_pt[rapBin]->SetPointError(total_nonzero_points,exl,exh,eyl,eyh);
	//				total_nonzero_points++;
	//			}
	//		}
		}
		outFile->cd();
		bFrac_pt[rapBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
		bFrac_pt[rapBin]->SetMarkerStyle(20+rapBin);
		bFrac_pt[rapBin]->SetMarkerColor(kRed+4*rapBin);
		bFrac_pt[rapBin]->SetLineColor(kRed+4*rapBin);
		bFrac_pt[rapBin]->SetLineWidth(2.5);
		//bFrac_pt[rapBin]->GetYaxis()->SetRangeUser(0.,0.8);
		//bFrac_pt[rapBin]->GetXaxis()->SetRangeUser(0.,30.);

		//prompt
		lth_pt_p[rapBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lth_pt_p[rapBin]->SetMarkerStyle(20+rapBin);
		lth_pt_p[rapBin]->SetMarkerColor(kRed+4*rapBin);
		lth_pt_p[rapBin]->SetLineColor(kRed+4*rapBin);
		//lth_pt_p[rapBin]->GetYaxis()->SetRangeUser(-1.,1.);
		//lth_pt_p[rapBin]->GetXaxis()->SetRangeUser(0.,30.);
		//lphi_pt_p[rapBin]->GetYaxis()->SetTitle(Form("#lamda^{%s}_{#theta}",scen));
		lphi_pt_p[rapBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lphi_pt_p[rapBin]->SetMarkerStyle(20+rapBin);
		lphi_pt_p[rapBin]->SetMarkerColor(kRed+4*rapBin);
		lphi_pt_p[rapBin]->SetLineColor(kRed+4*rapBin);
		//lphi_pt_p[rapBin]->GetYaxis()->SetRangeUser(-1.,1.);
		//lphi_pt_p[rapBin]->GetXaxis()->SetRangeUser(0.,30.);
		lthphi_pt_p[rapBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lthphi_pt_p[rapBin]->SetMarkerStyle(20+rapBin);
		lthphi_pt_p[rapBin]->SetMarkerColor(kRed+4*rapBin);
		lthphi_pt_p[rapBin]->SetLineColor(kRed+4*rapBin);
		//lthphi_pt_p[rapBin]->GetYaxis()->SetRangeUser(-1.,1.);
		//lthphi_pt_p[rapBin]->GetXaxis()->SetRangeUser(0.,30.);
		lthtilde_pt_p[rapBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lthtilde_pt_p[rapBin]->SetMarkerStyle(20+rapBin);
		lthtilde_pt_p[rapBin]->SetMarkerColor(kRed+4*rapBin);
		lthtilde_pt_p[rapBin]->SetLineColor(kRed+4*rapBin);
		//lthtilde_pt_p[rapBin]->GetYaxis()->SetRangeUser(-1.,1.);
		//lthtilde_pt_p[rapBin]->GetXaxis()->SetRangeUser(0.,30.);
		f_pt_p[rapBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
		f_pt_p[rapBin]->SetMarkerStyle(20+rapBin);
		f_pt_p[rapBin]->SetMarkerColor(kRed+4*rapBin);
		f_pt_p[rapBin]->SetLineColor(kRed+4*rapBin);
		//f_pt_p[rapBin]->GetYaxis()->SetRangeUser(0.,1.);
		//f_pt_p[rapBin]->GetXaxis()->SetRangeUser(0.,30.);

		//non-prompt
		lth_pt_np[rapBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lth_pt_np[rapBin]->SetMarkerStyle(20+rapBin);
		lth_pt_np[rapBin]->SetMarkerColor(kRed+4*rapBin);
		lth_pt_np[rapBin]->SetLineColor(kRed+4*rapBin);
		//lth_pt_np[rapBin]->GetYaxis()->SetRangeUser(-1.,1.);
		//lth_pt_np[rapBin]->GetXaxis()->SetRangeUser(0.,30.);
		lphi_pt_np[rapBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lphi_pt_np[rapBin]->SetMarkerStyle(20+rapBin);
		lphi_pt_np[rapBin]->SetMarkerColor(kRed+4*rapBin);
		lphi_pt_np[rapBin]->SetLineColor(kRed+4*rapBin);
		//lphi_pt_np[rapBin]->GetYaxis()->SetRangeUser(-1.,1.);
		//lphi_pt_np[rapBin]->GetXaxis()->SetRangeUser(0.,30.);
		lthphi_pt_np[rapBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lthphi_pt_np[rapBin]->SetMarkerStyle(20+rapBin);
		lthphi_pt_np[rapBin]->SetMarkerColor(kRed+4*rapBin);
		lthphi_pt_np[rapBin]->SetLineColor(kRed+4*rapBin);
		//lthphi_pt_np[rapBin]->GetYaxis()->SetRangeUser(-1.,1.);
		//lthphi_pt_np[rapBin]->GetXaxis()->SetRangeUser(0.,30.);
		lthtilde_pt_np[rapBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
		lthtilde_pt_np[rapBin]->SetMarkerStyle(20+rapBin);
		lthtilde_pt_np[rapBin]->SetMarkerColor(kRed+4*rapBin);
		lthtilde_pt_np[rapBin]->SetLineColor(kRed+4*rapBin);
		//lthtilde_pt_np[rapBin]->GetYaxis()->SetRangeUser(-1.,1.);
		//lthtilde_pt_np[rapBin]->GetXaxis()->SetRangeUser(0.,30.);
		f_pt_np[rapBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
		f_pt_np[rapBin]->SetMarkerStyle(20+rapBin);
		f_pt_np[rapBin]->SetMarkerColor(kRed+4*rapBin);
		f_pt_np[rapBin]->SetLineColor(kRed+4*rapBin);
		//f_pt_np[rapBin]->GetYaxis()->SetRangeUser(0.,1.);
		//f_pt_np[rapBin]->GetXaxis()->SetRangeUser(0.,30.);


		bFrac_pt[rapBin]->Write();
		lth_pt_p[rapBin]->Write();
		lphi_pt_p[rapBin]->Write();
		lthphi_pt_p[rapBin]->Write();
		lthtilde_pt_p[rapBin]->Write();
		f_pt_p[rapBin]->Write();
		lth_pt_np[rapBin]->Write();
		lphi_pt_np[rapBin]->Write();
		lthphi_pt_np[rapBin]->Write();
		lthtilde_pt_np[rapBin]->Write();
		f_pt_np[rapBin]->Write();

		if(total_nonzero_points==0) continue;

		cv[rapBin][0]=new TCanvas(Form("cv%d_0",rapBin),"");
		hg[rapBin][0]=cv[rapBin][0]->DrawFrame(0.,0.,30.,0.8);
		bFrac_pt[rapBin]->Draw("P");
		cv[rapBin][1]=new TCanvas(Form("cv%d_1",rapBin),"");
		hg[rapBin][1]=cv[rapBin][1]->DrawFrame(0.,-1.,30.,1.);
		lth_pt_p[rapBin]->Draw("P");
		cv[rapBin][2]=new TCanvas(Form("cv%d_2",rapBin),"");
		hg[rapBin][2]=cv[rapBin][2]->DrawFrame(0.,-1.,30.,1.);
		lphi_pt_p[rapBin]->Draw("P");
		cv[rapBin][3]=new TCanvas(Form("cv%d_3",rapBin),"");
		hg[rapBin][3]=cv[rapBin][3]->DrawFrame(0.,-1.,30.,1.);
		lthphi_pt_p[rapBin]->Draw("P");
		cv[rapBin][4]=new TCanvas(Form("cv%d_4",rapBin),"");
		hg[rapBin][4]=cv[rapBin][4]->DrawFrame(0.,-1.,30.,1.);
		lthtilde_pt_p[rapBin]->Draw("P");
		cv[rapBin][5]=new TCanvas(Form("cv%d_5",rapBin),"");
		hg[rapBin][5]=cv[rapBin][5]->DrawFrame(0.,0.,30.,1.);
		f_pt_p[rapBin]->Draw("P");
		cv[rapBin][6]=new TCanvas(Form("cv%d_6",rapBin),"");
		hg[rapBin][6]=cv[rapBin][6]->DrawFrame(0.,-1.,30.,1.);
		lth_pt_np[rapBin]->Draw("P");
		cv[rapBin][7]=new TCanvas(Form("cv%d_7",rapBin),"");
		hg[rapBin][7]=cv[rapBin][7]->DrawFrame(0.,-1.,30.,1.);
		lphi_pt_np[rapBin]->Draw("P");
		cv[rapBin][8]=new TCanvas(Form("cv%d_8",rapBin),"");
		hg[rapBin][8]=cv[rapBin][8]->DrawFrame(0.,-1.,30.,1.);
		lthphi_pt_np[rapBin]->Draw("P");
		cv[rapBin][9]=new TCanvas(Form("cv%d_9",rapBin),"");
		hg[rapBin][9]=cv[rapBin][9]->DrawFrame(0.,-1.,30.,1.);
		lthtilde_pt_np[rapBin]->Draw("P");
		cv[rapBin][10]=new TCanvas(Form("cv%d_10",rapBin),"");
		hg[rapBin][10]=cv[rapBin][10]->DrawFrame(0.,0.,30.,1.);
		f_pt_np[rapBin]->Draw("P");

		cv[rapBin][0]->SaveAs(Form("%s/bFrac_pt_rap%d-%s.png",path,rapBin,scen));
		cv[rapBin][1]->SaveAs(Form("%s/lth_pt_p_rap%d-%s.png",path,rapBin,scen));
		cv[rapBin][2]->SaveAs(Form("%s/lphi_pt_p_rap%d-%s.png",path,rapBin,scen));
		cv[rapBin][3]->SaveAs(Form("%s/lthphi_pt_p_rap%d-%s.png",path,rapBin,scen));
		cv[rapBin][4]->SaveAs(Form("%s/lthtilde_pt_p_rap%d-%s.png",path,rapBin,scen));
		cv[rapBin][5]->SaveAs(Form("%s/f_pt_p_rap%d-%s.png",path,rapBin,scen));
		cv[rapBin][6]->SaveAs(Form("%s/lth_pt_np_rap%d-%s.png",path,rapBin,scen));
		cv[rapBin][7]->SaveAs(Form("%s/lphi_pt_np_rap%d-%s.png",path,rapBin,scen));
		cv[rapBin][8]->SaveAs(Form("%s/lthphi_pt_np_rap%d-%s.png",path,rapBin,scen));
		cv[rapBin][9]->SaveAs(Form("%s/lthtilde_pt_np_rap%d-%s.png",path,rapBin,scen));
		cv[rapBin][10]->SaveAs(Form("%s/f_pt_np_rap%d-%s.png",path,rapBin,scen));
		
		//bFrac_pt[rapBin]->SetHistogram(hist_bFrac_pt[rapBin]);
		//hist_bFrac_pt[rapBin]=(TH1F *)bFrac_pt[rapBin]->GetHistogram();
	}
	outFile->Close();

}


