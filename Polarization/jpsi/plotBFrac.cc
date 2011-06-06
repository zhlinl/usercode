
#include "commonVar.h"

void plotBFrac(){

	//double nNP[2][jpsi::kNbPTMaxBins];
	//double nP[2][jpsi::kNbPTMaxBins];
	const int NRaps=2;
	double nNP,nP;
	double errNP,errP;
	double bFraction[NRaps][jpsi::kNbPTMaxBins];
	double bFractionErrl[NRaps][jpsi::kNbPTMaxBins];
	double bFractionErrh[NRaps][jpsi::kNbPTMaxBins];
	double ptMean[NRaps][jpsi::kNbPTMaxBins],errptMeanl[NRaps][jpsi::kNbPTMaxBins],errptMeanh[NRaps][jpsi::kNbPTMaxBins];

	char path[100];
	int nbins=50;
	sprintf(path,"picLM/Run2010A");
	//gSystem->mkdir(path);

	for(int rapBin=0;rapBin<NRaps;rapBin++){
	for(int ptBin=0;ptBin<jpsi::jpsi::kNbPTBins[rapBin+1];ptBin++){

		ptMean[rapBin][ptBin]=(jpsi::pTRange[rapBin+1][ptBin]+jpsi::pTRange[rapBin+1][ptBin+1])/2;
		errptMeanl[rapBin][ptBin]=fabs(ptMean[rapBin][ptBin]-jpsi::pTRange[rapBin+1][ptBin]);
		errptMeanh[rapBin][ptBin]=fabs(jpsi::pTRange[rapBin+1][ptBin+1]-ptMean[rapBin][ptBin]);
		bFraction[rapBin][ptBin]=0;
		bFractionErrl[rapBin][ptBin]=0;
		bFractionErrh[rapBin][ptBin]=0;
	}
	}

	char scen[4]="CS";
	char inName[200];
	for(int rapBin=0;rapBin<NRaps;rapBin++){
		for(int ptBin=0;ptBin<jpsi::kNbPTBins[rapBin+1];ptBin++){

			//char inName[200];
			sprintf(inName,Form("/afs/ihep.ac.cn/users/z/zhangll/workspace/work/polarization/Polarization/project/polaRun2010A_rap%d_pt%d-CS.root",rapBin+1,ptBin+1));
			cout<<"inName: "<<inName<<endl;
			TFile *inFile=new TFile(inName,"R");
			if(!inFile){
				cout<<">>===Error: no input File"<<endl;
				continue;
			}
			RooWorkspace *ws=(RooWorkspace *)inFile->Get(Form("polaRun2010A_rap%d_pt%d",rapBin+1,ptBin+1));
			if(!ws){
				cout<<">>===Error: no workspace in root file"<<endl;
				continue;
			}
			//RooRealVar *JpsiMass=(RooRealVar*)ws->var("JpsiMass");
			RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);


			RooDataSet *data=(RooDataSet *)ws->data(Form("data_rap%d_pt%d",rapBin+1,ptBin+1));

			RooArgSet *dataVars=(RooArgSet *)data->get();
			dataVars->Print();
			TIterator* di=(TIterator*)dataVars->createIterator();

			RooRealVar *nNP_=(RooRealVar *)ws->var("nNonPromptSignal");
			RooRealVar *nP_=(RooRealVar *)ws->var("nPromptSignal");
			nNP=nNP_->getVal();
			nP=nP_->getVal();
			errNP=nNP_->getError();
			errP=nP_->getError();
			bFraction[rapBin][ptBin]=nNP[rapBin][ptBin]/(nNP[rapBin][ptBin]+nP[rapBin][ptBin]);
			if(bFraction[rapBin][ptBin]>1.) bFraction[rapBin][ptBin]=0;
			bFractionErrl[rapBin][ptBin] = ((1.0/(nNP+nP)+nNP/((nNP+nP)*(nNP+nP)))*errNP + nNP/((nNP+nP)*(nNP+nP))*errP);
			bFractionErrh[rapBin][ptBin] = bFractionErrl[rapBin][ptBin];

			cout<<">>==nNP_rap"<<rapBin+1<<"_pt_"<<ptBin+1<<": "<<nNP<<endl;
			cout<<">>==nP_rap"<<rapBin+1<<"_pt_"<<ptBin+1<<": "<<nP<<endl;
			cout<<">>==bFraction_rap"<<rapBin+1<<"_pt_"<<ptBin+1<<": "<<bFraction[rapBin][ptBin]<<endl;

			//delete inFile;delete ws;delete data;delete dataVars;

		}
	}

	TGraphAsymmErrors *bFracGraph[2];
	bFracGraph[0]=new TGraphAsymmErrors(jpsi::kNbPTMaxBins,ptMean[0],bFraction[0],errptMeanl[0],errptMeanh[0],bFractionErrl[0],bFractionErrh[0]);
	bFracGraph[0]->SetMarkerStyle(20);
	bFracGraph[0]->SetMarkerColor(kRed);
	bFracGraph[0]->SetTitle("B-fraction rapidity bin 1");
	bFracGraph[0]->GetXaxis()->SetTitle("p_{T} [GeV]");

	bFracGraph[1]=new TGraphAsymmErrors(jpsi::kNbPTMaxBins,ptMean[1],bFraction[1],errptMeanl[1],errptMeanh[1],bFractionErrl[1],bFractionErrh[1]);
	bFracGraph[1]->SetMarkerStyle(21);
	bFracGraph[1]->SetMarkerColor(kRed);
	bFracGraph[1]->SetTitle("B-fraction rapidity bin 2");
	bFracGraph[1]->GetXaxis()->SetTitle("p_{T} [GeV]");

	TCanvas *c1=new TCanvas("c1","");
	bFracGraph[0]->Draw("AP");
	TCanvas *c2=new TCanvas("c2","");
	bFracGraph[1]->Draw("AP");
	c1->SaveAs(Form("%s/bFrac_rap1.png",path));
	c2->SaveAs(Form("%s/bFrac_rap2.png",path));



}
