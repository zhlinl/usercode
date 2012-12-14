/*
 * EvaluateSyst.cc
 *
 *  Created on: Dec 3, 2011
 *      Author: valentinknuenz
 */

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "../../interface/rootIncludes.inc"
#include "../../interface/commonVar.h"
#include "ToyMC.h"

#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"

#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"



int main(int argc, char** argv) {

	Char_t *storagedir = "Default"; //Storage Directory
	Char_t *basedir = "Default"; //Storage Directory
  	Char_t *JobID1 = "Default";
  	Char_t *JobID2 = "Default";
  	Char_t *SystDir = "Default";

	int ptBinMin=1;
	int ptBinMax=1;
	int nState=1;

	bool statErrConsideration=false;

	  for( int i=0;i < argc; ++i ) {

		    if(std::string(argv[i]).find("JobID1") != std::string::npos) {char* JobID1char = argv[i]; char* JobID1char2 = strtok (JobID1char, "="); JobID1 = JobID1char2; cout<<"JobID1 = "<<JobID1<<endl;}
		    if(std::string(argv[i]).find("JobID2") != std::string::npos) {char* JobID2char = argv[i]; char* JobID2char2 = strtok (JobID2char, "="); JobID2 = JobID2char2; cout<<"JobID2 = "<<JobID2<<endl;}
		    if(std::string(argv[i]).find("SystDir") != std::string::npos) {char* SystDirchar = argv[i]; char* SystDirchar2 = strtok (SystDirchar, "="); SystDir = SystDirchar2; cout<<"SystDir = "<<SystDir<<endl;}
		    if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
		    if(std::string(argv[i]).find("basedir") != std::string::npos) {char* basedirchar = argv[i]; char* basedirchar2 = strtok (basedirchar, "="); basedir = basedirchar2; cout<<"basedir = "<<basedir<<endl;}

		    if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
		    if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
		    if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}

		    if(std::string(argv[i]).find("statErrConsideration=true") != std::string::npos) {statErrConsideration=true; cout<<"statErrConsideration"<<endl;}

	    }

	    char tmpfilename[200];
		sprintf(tmpfilename,"%s/macros/polFit/%s/TGraphResults_%dSUps.root",basedir,SystDir,nState);
		gSystem->Unlink(tmpfilename);

		char filename[200];
		sprintf(filename,"%s/%s/TGraphResults_%dSUps.root",storagedir,JobID1,nState);
//		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/Systematics/TotalSyst/TotalSquaredSystApr27/TGraphResults_%dSUps.root",nState);// ifCentralsWithTotalSyst
		TFile *infile1 = new TFile(filename,"READ");
		sprintf(filename,"%s/%s/TGraphResults_%dSUps.root",storagedir,JobID2,nState);
//		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/Systematics/TotalSyst/Apr25CentralWithStatSystSquared/TGraphResults_2SUps.root");
		TFile *infile2 = new TFile(filename,"READ");

//		sprintf(filename,"%s/MCclosure_CombinedStates_MCtruthFineEta_3Sig/TGraphResults_%dSUps.root",storagedir,nState);//ifAddInforFromThirdTGraph
//		TFile *infile3 = new TFile(filename,"READ");//ifAddInforFromThirdTGraph

		char GraphName[200];

		for(int iLam = 1; iLam<19; iLam++){

		for(int rapBin = 1; rapBin < 3; rapBin++){

		sprintf(filename,"%s/macros/polFit/%s/TGraphResults_%dSUps.root",basedir,SystDir,nState);
		TFile *outfile = new TFile(filename,"UPDATE");


		if(iLam==1)  sprintf(GraphName,"lth_CS_rap%d",rapBin);
		if(iLam==2)  sprintf(GraphName,"lph_CS_rap%d",rapBin);
		if(iLam==3)  sprintf(GraphName,"ltp_CS_rap%d",rapBin);
		if(iLam==4)  sprintf(GraphName,"lthstar_CS_rap%d",rapBin);
		if(iLam==5)  sprintf(GraphName,"lphstar_CS_rap%d",rapBin);
		if(iLam==6)  sprintf(GraphName,"ltilde_CS_rap%d",rapBin);

		if(iLam==7)  sprintf(GraphName,"lth_HX_rap%d",rapBin);
		if(iLam==8)  sprintf(GraphName,"lph_HX_rap%d",rapBin);
		if(iLam==9)  sprintf(GraphName,"ltp_HX_rap%d",rapBin);
		if(iLam==10) sprintf(GraphName,"lthstar_HX_rap%d",rapBin);
		if(iLam==11) sprintf(GraphName,"lphstar_HX_rap%d",rapBin);
		if(iLam==12) sprintf(GraphName,"ltilde_HX_rap%d",rapBin);

		if(iLam==13) sprintf(GraphName,"lth_PX_rap%d",rapBin);
		if(iLam==14) sprintf(GraphName,"lph_PX_rap%d",rapBin);
		if(iLam==15) sprintf(GraphName,"ltp_PX_rap%d",rapBin);
		if(iLam==16) sprintf(GraphName,"lthstar_PX_rap%d",rapBin);
		if(iLam==17) sprintf(GraphName,"lphstar_PX_rap%d",rapBin);
		if(iLam==18) sprintf(GraphName,"ltilde_PX_rap%d",rapBin);

		int nFrame=0;

		if(iLam>0&&iLam<7) nFrame=1;
		if(iLam>6&&iLam<13) nFrame=2;
		if(iLam>12&&iLam<19) nFrame=3;

		TGraphAsymmErrors* graph1 = (TGraphAsymmErrors*) infile1->Get(GraphName);
		TGraphAsymmErrors* graph2 = (TGraphAsymmErrors*) infile2->Get(GraphName);

//		TGraphAsymmErrors* graph3 = (TGraphAsymmErrors*) infile3->Get(GraphName);//ifAddInforFromThirdTGraph


		int nBinspT=ptBinMax-ptBinMin+1;
		double ptCentre_[nBinspT];
		double ptCentreErr_low[nBinspT];
		double ptCentreErr_high[nBinspT];
		double lmean[nBinspT];

		double ptCentre1_[nBinspT];
		double ptCentreErr1_low[nBinspT];
		double ptCentreErr1_high[nBinspT];
		double ptCentre2_[nBinspT];
		double ptCentreErr2_low[nBinspT];
		double ptCentreErr2_high[nBinspT];

		double ptCentre3_[nBinspT];
		double lmean3[nBinspT];

		double lmean1[nBinspT];
		double lmean2[nBinspT];

		double lmeanErr1_low[nBinspT];
		double lmeanErr1_high[nBinspT];
		double lmeanErr2_low[nBinspT];
		double lmeanErr2_high[nBinspT];

		double lmeanErr_change_low[nBinspT];
		double lmeanErr_change_high[nBinspT];
		double lmeanErr_change[nBinspT];



		int pt=0;
		for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

			graph1->GetPoint(pt,ptCentre1_[pt],lmean1[pt]);
			ptCentreErr1_high[pt]=graph1->GetErrorXhigh(pt);
			ptCentreErr1_low[pt]=graph1->GetErrorXlow(pt);
			graph2->GetPoint(pt,ptCentre2_[pt],lmean2[pt]);
			ptCentreErr2_high[pt]=graph2->GetErrorXhigh(pt);
			ptCentreErr2_low[pt]=graph2->GetErrorXlow(pt);

			if(lmean1[pt]>998) lmean1[pt]=9999;

			ptCentre_[pt]=(ptCentre1_[pt]+ptCentre1_[pt])/2.;
			ptCentreErr_low[pt]=(ptCentreErr1_low[pt]+ptCentreErr2_low[pt])/2.;
			ptCentreErr_high[pt]=(ptCentreErr1_high[pt]+ptCentreErr2_high[pt])/2.;

			lmeanErr1_high[pt]=graph1->GetErrorYhigh(pt);
			lmeanErr1_low[pt]=graph1->GetErrorYlow(pt);
			lmeanErr2_high[pt]=graph2->GetErrorYhigh(pt);
			lmeanErr2_low[pt]=graph2->GetErrorYlow(pt);

			lmean[pt]=lmean1[pt]-lmean2[pt];

			lmeanErr_change_high[pt]=TMath::Sqrt(lmeanErr2_high[pt]*lmeanErr2_high[pt]+lmean1[pt]*lmean1[pt]);
			lmeanErr_change_low[pt]=TMath::Sqrt(lmeanErr2_low[pt]*lmeanErr2_low[pt]+lmean1[pt]*lmean1[pt]);
			lmeanErr_change[pt]=(lmeanErr_change_low[pt]+lmeanErr_change_high[pt])/2.;



//			lmean[pt]=lmean1[pt]+lmean2[pt];//ifCorrectCentralResultsForBias

//BiasCorrectionError

			lmeanErr_change[pt]=999;
//			lmeanErr_change[pt]=TMath::Sqrt(TMath::Power((lmeanErr1_high[pt]+lmeanErr1_low[pt])/2.,2.)+TMath::Power((lmeanErr2_high[pt]+lmeanErr2_low[pt])/2.,2.));
			if((lmeanErr1_high[pt]+lmeanErr1_low[pt])/2.>(lmeanErr2_high[pt]+lmeanErr2_low[pt])/2.) lmeanErr_change[pt]=TMath::Sqrt(TMath::Power((lmeanErr1_high[pt]+lmeanErr1_low[pt])/2.,2.)-TMath::Power((lmeanErr2_high[pt]+lmeanErr2_low[pt])/2.,2.));
			if((lmeanErr1_high[pt]+lmeanErr1_low[pt])/2.<(lmeanErr2_high[pt]+lmeanErr2_low[pt])/2.) lmeanErr_change[pt]=TMath::Sqrt(TMath::Power((lmeanErr2_high[pt]+lmeanErr2_low[pt])/2.,2.)-TMath::Power((lmeanErr1_high[pt]+lmeanErr1_low[pt])/2.,2.));

			double fracErrHigh2=lmeanErr2_high[pt]/(lmeanErr2_high[pt]+lmeanErr2_low[pt]);


//	unc		lmeanErr1_high[pt]=fracErrHigh2*(lmeanErr1_high[pt]+lmeanErr1_low[pt]);
//	unc		lmeanErr1_low[pt]=(1-fracErrHigh2)*(lmeanErr1_high[pt]+lmeanErr1_low[pt]);


//			graph3->GetPoint(pt,ptCentre3_[pt],lmean3[pt]);//ifAddInforFromThirdTGraph
//			lmean[pt]=lmean1[pt]-lmean2[pt]+lmean3[pt];

/*			if(ptBin<6){
				lmean[pt]=lmean1[pt]-lmean3[pt];//ifAddInforFromThirdTGraph
			}
*/
//			lmean[pt]=lmean[pt]/TMath::Sqrt(12.);//ifSQRT12

			double SigComp=1;

			double lmeanBuff=lmean[pt];
			if(statErrConsideration){
				cout<<"StatErrCheck"<<endl;
				bool getSqrt=true;
				if(lmeanBuff*lmeanBuff-SigComp*SigComp*lmeanErr2_low[pt]*lmeanErr2_low[pt]-SigComp*SigComp*lmeanErr1_high[pt]*lmeanErr1_high[pt]<0 && lmeanBuff < 0) {lmean[pt]=0;getSqrt=false;}
				if(lmeanBuff*lmeanBuff-SigComp*SigComp*lmeanErr1_low[pt]*lmeanErr1_low[pt]-SigComp*SigComp*lmeanErr2_high[pt]*lmeanErr2_high[pt]<0 && lmeanBuff > 0) {lmean[pt]=0;getSqrt=false;}
				if(lmeanBuff==0) {lmean[pt]=0;getSqrt=false;}
				if(lmeanBuff < 0 && getSqrt) lmean[pt]=-TMath::Sqrt(lmeanBuff*lmeanBuff-SigComp*SigComp*lmeanErr2_low[pt]*lmeanErr2_low[pt]-SigComp*SigComp*lmeanErr1_high[pt]*lmeanErr1_high[pt]);
				if(lmeanBuff > 0 && getSqrt) lmean[pt]=TMath::Sqrt(lmeanBuff*lmeanBuff-SigComp*SigComp*lmeanErr1_low[pt]*lmeanErr1_low[pt]-SigComp*SigComp*lmeanErr2_high[pt]*lmeanErr2_high[pt]);
			}

//			if(lmeanBuff < 0) lmean[pt]=lmeanBuff/TMath::Sqrt(lmeanErr2_low[pt]*lmeanErr2_low[pt]+lmeanErr1_high[pt]*lmeanErr1_high[pt]); //ifPull
//			else lmean[pt]=lmeanBuff/TMath::Sqrt(lmeanErr1_low[pt]*lmeanErr1_low[pt]+lmeanErr2_high[pt]*lmeanErr2_high[pt]); //ifPull

		pt++;
		}




//		TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,0,0);
//		TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre2_,lmean2,ptCentreErr2_low,ptCentreErr2_high,lmean1,lmean1);// ifCentralsWithTotalSyst
//		TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,lmeanErr_change,lmeanErr_change);// ifCalculateBiasCorrection
//		TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,lmeanErr1_low,lmeanErr1_high);// ifCorrectCentralResultsForBias
		TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean2,ptCentreErr_low,ptCentreErr_high,lmeanErr1_low,lmeanErr1_high);// if 'take central value from JobID2, take error from JobID1'
//		TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,lmeanErr1_low,lmeanErr1_high);// ifAddInforFromThirdTGraph
		graphSyst->SetMarkerColor(ToyMC::MarkerColor[nFrame]);
		graphSyst->SetLineColor(ToyMC::MarkerColor[nFrame]);
		graphSyst->SetMarkerStyle(ToyMC::MarkerStyle[nState][rapBin]);
		graphSyst->SetMarkerSize(ToyMC::MarkerSize[nState][rapBin]);
		graphSyst->SetName(GraphName);


		outfile->cd();
		graphSyst->Draw("P");
		graphSyst->Write();

		outfile->Write();
		outfile->Close();
		delete outfile;
		outfile = NULL;


		}


		}





	return 0;
}
