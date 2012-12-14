#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"

#include "TH2F.h"

enum{L,R};
const char *bgLabel[2] = {"L", "R"};
TH2F *hCosThetaPhi[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::kNbFrames][2];
TH2F *hCosThetaPhiHighct[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::kNbFrames][2];
TH2F *hCosThetaPhiNPBG[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::kNbFrames];
//TH1F *events_SR[onia::kNbRapForPTBins][onia::kNbPTMaxBins];
void LoadHistos(Int_t iRapBin, Int_t iPTBin, Int_t nState);
void PlotHistos(Int_t iRapBin, Int_t iPTBin, Int_t iFrame, Int_t iWindow);
//===========================

//========================================================
// code to read input arguments
	template<typename T>
void fromSplit(const std::string& key, const std::string &arg, T& out)
{
	const char delim = '=';
	// Skip if key or delimiter not there
	if ((arg.find(key) == std::string::npos) ||
			(arg.find(delim) == std::string::npos))
		return;

	std::string skey, sval;
	std::stringstream sstr(arg);
	std::getline(sstr, skey, delim); // Dummy read to skip key
	std::getline(sstr, sval, delim); // Get value
	T tout;
	if (!(std::istringstream(sval) >> std::boolalpha >> tout))
		return;
	out = tout;
	std::cout << std::boolalpha << skey << ": "  << out << std::endl;
}
// Special version for string without the conversion 
	template<>
void fromSplit(const std::string& key, const std::string &arg, std::string &out)
{
	const char delim = '=';
	// Skip if key or delimiter not there
	if ((arg.find(key) == std::string::npos) ||
			(arg.find(delim) == std::string::npos))
		return;
	std::string skey, sval;
	std::stringstream sstr(arg);
	std::getline(sstr, skey, delim); // Dummy read to skip key
	std::getline(sstr, sval, delim); // Get value
	out = sval;
	std::cout << skey << ": "  << out << std::endl;
}


int main(int argc, char* argv[]){

	// set default values
	int nState = 999;

	// Loop over argument list                                                          
	for (int i=1; i < argc; i++){
		std::string arg = argv[i];
		fromSplit("nState", arg, nState);
	}

	gStyle->SetPalette(1);

	for(int iRap = 1; iRap < onia::kNbRapForPTBins+1; iRap++){
		for(int iPT = 1; iPT < onia::kNbPTBins[iRap]+1; iPT++){
			LoadHistos(iRap, iPT, nState);
			for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
				PlotHistos(iRap, iPT, iFrame, L);
				PlotHistos(iRap, iPT, iFrame, R);
			}
		}
	}
	return 1;
}

//===========================
void PlotHistos(Int_t iRapBin, Int_t iPTBin, Int_t iFrame, Int_t iWindow){

	double lvalue = 0.28, tvalue = 0.92;
	double left=lvalue, top=tvalue, textSize=0.035;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	gStyle->SetPadRightMargin(0.12);
	gStyle->SetOptStat(0);
	gStyle->SetFrameBorderMode(0);

	Char_t name[200];
	sprintf(name, "c1_%s_rap%d_pT%d_%s", onia::frameLabel[iFrame], iRapBin, iPTBin, bgLabel[iWindow]);
	TCanvas *c1 = new TCanvas(name, "", 500, 500);
	c1->SetFillColor(kWhite);

	double yOffset=1.4;
	hCosThetaPhi[iRapBin][iPTBin][iFrame][iWindow]->GetYaxis()->SetTitleOffset(yOffset);
	gPad->SetLeftMargin(0.125);

	hCosThetaPhi[iRapBin][iPTBin][iFrame][iWindow]->Draw("colz");

	if(iRapBin==1) 
		latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
					onia::rapForPTRange[iRapBin],
					onia::pTRange[iRapBin][iPTBin-1],onia::pTRange[iRapBin][iPTBin]));
	else
		latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
					onia::rapForPTRange[iRapBin-1],onia::rapForPTRange[iRapBin],
					onia::pTRange[iRapBin][iPTBin-1],onia::pTRange[iRapBin][iPTBin]));

	sprintf(name, "Figures/cosThetaPhi_%s_rap%d_pT%d_%s.pdf", onia::frameLabel[iFrame], iRapBin, iPTBin, bgLabel[iWindow]);
	c1->Print(name);

	sprintf(name, "c2_%s_rap%d_pT%d_%s", onia::frameLabel[iFrame], iRapBin, iPTBin, bgLabel[iWindow]);
	TCanvas *c2 = new TCanvas(name, "", 500, 500);
	c2->SetFillColor(kWhite);
	hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][iWindow]->GetYaxis()->SetTitleOffset(yOffset);
	hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][iWindow]->Draw("colz");

	if(iRapBin==1) 
		latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
					onia::rapForPTRange[iRapBin],
					onia::pTRange[iRapBin][iPTBin-1],onia::pTRange[iRapBin][iPTBin]));
	else
		latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
					onia::rapForPTRange[iRapBin-1],onia::rapForPTRange[iRapBin],
					onia::pTRange[iRapBin][iPTBin-1],onia::pTRange[iRapBin][iPTBin]));

	sprintf(name, "Figures/cosThetaPhi_%s_rap%d_pT%d_highct_%s.pdf", onia::frameLabel[iFrame], iRapBin, iPTBin, bgLabel[iWindow]);
	c2->Print(name);

	if(iWindow==0){
		sprintf(name, "c3_%s_rap%d_pT%d", onia::frameLabel[iFrame], iRapBin, iPTBin);
		TCanvas *c3 = new TCanvas(name, "", 500, 500);
		c3->SetFillColor(kWhite);
		hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame]->GetYaxis()->SetTitleOffset(yOffset);
		hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame]->Draw("colz");

		if(iRapBin==1) 
			latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
						onia::rapForPTRange[iRapBin],
						onia::pTRange[iRapBin][iPTBin-1],onia::pTRange[iRapBin][iPTBin]));
		else
			latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
						onia::rapForPTRange[iRapBin-1],onia::rapForPTRange[iRapBin],
						onia::pTRange[iRapBin][iPTBin-1],onia::pTRange[iRapBin][iPTBin]));

		sprintf(name, "Figures/cosThetaPhiNPBG_%s_rap%d_pT%d.pdf", onia::frameLabel[iFrame], iRapBin, iPTBin);
		c3->Print(name);

		delete c3;

	}

	delete c1;
	delete c2;
}

//===========================
void LoadHistos(Int_t iRapBin, Int_t iPTBin, Int_t nState){
	cout<<"rap "<<iRapBin<<" pt "<<iPTBin<<endl;

	Char_t name[200];
	sprintf(name, "tmpFiles/data_Psi%dS_rap%d_pT%d.root", nState-3, iRapBin, iPTBin);

	TFile *fIn = new TFile(name);

	//sprintf(name, "events_SR");
	//events_SR[iRapBin][iPTBin] = (TH1F*) fIn->Get(name);
	//cout<<"rap"<<iRapBin<<" pt"<<iPTBin<<" Prompt signal in prompt region: "<<events_SR[iRapBin][iPTBin]->GetBinContent(1)<<endl;

	for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
		sprintf(name, "hBG_cosThetaPhi_%s_L", onia::frameLabel[iFrame]);
		cout<<"name: "<<name<<endl;
		hCosThetaPhi[iRapBin][iPTBin][iFrame][L] = (TH2F *) fIn->Get(name);
		//hCosThetaPhi[iRapBin][iPTBin][iFrame][L]->Print();
		sprintf(name, "hCosThetaPhi_%s_pT%d_rap%d_L", onia::frameLabel[iFrame], iRapBin, iPTBin);
		hCosThetaPhi[iRapBin][iPTBin][iFrame][L]->SetName(name);

		sprintf(name, "hBG_cosThetaPhi_%s_R", onia::frameLabel[iFrame]);
		hCosThetaPhi[iRapBin][iPTBin][iFrame][R] = (TH2F *) fIn->Get(name);
		sprintf(name, "hCosThetaPhi_%s_pT%d_rap%d_R", onia::frameLabel[iFrame], iRapBin, iPTBin);
		hCosThetaPhi[iRapBin][iPTBin][iFrame][R]->SetName(name);

		sprintf(name, "hBGinNP_cosThetaPhi_%s_L", onia::frameLabel[iFrame]);
		hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][L] = (TH2F *) fIn->Get(name);
		//hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][L]->Print();
		sprintf(name, "hCosThetaPhi_%s_pT%d_rap%d_highct_L", onia::frameLabel[iFrame], iRapBin, iPTBin);
		hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][L]->SetName(name);

		sprintf(name, "hBGinNP_cosThetaPhi_%s_R", onia::frameLabel[iFrame]);
		hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][R] = (TH2F *) fIn->Get(name);
		sprintf(name, "hCosThetaPhi_%s_pT%d_rap%d_highct_R", onia::frameLabel[iFrame], iRapBin, iPTBin);
		hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][R]->SetName(name);

		sprintf(name, "hNPBG_cosThetaPhi_%s", onia::frameLabel[iFrame]);
		hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame] = (TH2F *) fIn->Get(name);
		//hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame]->Print();
		sprintf(name, "hCosThetaPhiNPBG_%s_pT%d_rap%d", onia::frameLabel[iFrame], iRapBin, iPTBin);
		hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame]->SetName(name);

		if(iFrame==2){
			std::cout<<"rap "<<iRapBin<<" pt "<<iPTBin<<" nBinsCosthBG: "<<
				hCosThetaPhi[iRapBin][iPTBin][iFrame][L]->GetNbinsX()<<std::endl;
			std::cout<<"rap "<<iRapBin<<" pt "<<iPTBin<<" nBinsPhiBG:   "<<
				hCosThetaPhi[iRapBin][iPTBin][iFrame][L]->GetNbinsY()<<std::endl;
			std::cout<<"rap "<<iRapBin<<" pt "<<iPTBin<<" nBinsCosthNPBG: "<<
				hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame]->GetNbinsX()<<std::endl;
			std::cout<<"rap "<<iRapBin<<" pt "<<iPTBin<<" nBinsPhiNPBG:   "<<
				hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame]->GetNbinsY()<<std::endl;
		}

	}

}
