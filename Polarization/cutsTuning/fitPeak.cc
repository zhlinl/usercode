#ifndef __CINT__
#include <RooGlobalFunc.h>
#endif

#include <TLatex.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooLandau.h>
#include <RooChebychev.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooDataHist.h>
#include <RooVoigtian.h>
#include <RooCBShape.h>
#include <TCanvas.h>
#include <TROOT.h>

#include <TAxis.h>
#include <TH1.h>
#include <TLine.h>
#include <TMath.h>
#include <TVector3.h>
#include <vector.h>
#include <TString.h>
#include <iostream.h>
using namespace RooFit;

Int_t getNumDigs(Double_t num){
	Int_t numDig;
	numDig = 0;
	if(num>=5) numDig = 0;
	else if(num>=0.5) numDig = 1;
	else if(num>=0.05) numDig = 2;
	else if(num>=0.005) numDig = 3;
	else if(num>=0.0005) numDig = 4;
	else if(num>=0.00005) numDig = 5;
	else if(num>=0.000005) numDig = 6;
	else if(num>=0.0000005) numDig = 7;
	else if(num>=0.00000005) numDig = 8;
	return numDig;
}

vector<double> fitPeak(RooRealVar x, Int_t whichBkg=0,const Int_t nPeak=1, TH1* hist=0, RooDataSet* dataSet, vector<Double_t> fMean, vector<Double_t> dfMean,vector<Double_t> fSigma, vector<Double_t> dfSigma,Double_t fbkg=0.5)

{
	//---------------Background Define---------------------------
	Bool_t showRlt=1;
	Bool_t getSB=1;
	Long_t nEntries;
	//ofstream log("log.txt");
	vector<Double_t> cInit;
	RooChebychev *bkg;
	RooArgList listCf("listCf");
	if(whichBkg==0)
	{
		Int_t nOrder=3;
		cout<<"===Now Using the 'RooChebychev' as Background==="<<endl;
		//cout<<"===Please input the orders using===:"<<endl;
		//cin>>nOrder;	
		cInit.clear();
		if(nOrder>=0 && nOrder<=7)
		{
			RooRealVar *cf;
			for(Int_t i=0; i<nOrder; i++)
			{		
				cInit.push_back(0.1);
				cf=new RooRealVar(Form("cf%d",i),Form("cf%d",i),cInit[i],-1,1);
				listCf.add(*cf);
			}
			bkg = new RooChebychev("bkg","Background",x,listCf);
		}
		else
		{
			cout<<"Error: orders must be 'int' type between [0-7]"<<endl;
		}
	}	
	else 
	{
		cout<<"Error:  there is no Background function you input, please input the number between(0--3)"<<endl;
	}

	//-----------------Peak Define-------------------------------
	if(nPeak<0 || nPeak>1000)
	{
		cout<<"No Peak or nPeak is too large to fit, quit"<<endl;
		return;
	}

	RooAbsData *dataHist;
	if(hist)
	{ 
		dataHist=new RooDataHist("dataHist","dataHist",x,hist);
		nEntries=(Long_t(hist->GetEntries()));
	}
	else if(dataSet)
	{ 
		//RooDataSet *dataSet=new RooDataSet("dataSet","dataSet",x,tree);
		nEntries=dataSet->numEntries();
		//nEntries=dataSet->sumEntries();
	}

	RooRealVar *mean[nPeak];
	RooRealVar *sigma[nPeak];
	RooArgList listPeak("listPeak");
	//RooVoigtian *Peak[nPeak];
	RooGaussian *peak[nPeak];
	RooArgList listArea("listArea");
	RooRealVar *areaPeak[nPeak];
	RooRealVar *areaBkg;

	Double_t val,valL,valR;
	for(Int_t i=0; i<nPeak; i++)
	{
		val=fMean[i];
		valL=fMean[i]-dfMean[i];
		valR=fMean[i]+dfMean[i];
		mean[i]=new RooRealVar(Form("mean%d",i),Form("mean%d",i),val,valL,valR);
		//log<<"mean:"<<val<<endl;
		//log<<"meanLo:"<<valL<<endl;
		//log<<"meanHi:"<<valR<<endl;

		val=fSigma[i];
		valL=fSigma[i]-dfSigma[i];
		valR=fSigma[i]+dfSigma[i];
		sigma[i]=new RooRealVar(Form("sigma%d",i),Form("sigma%d",i),val,valL,valR);
		peak[i]=new RooGaussian(Form("peak%d",i),Form("peak%d",i),x,*mean[i],*sigma[i]);
		listPeak.add(*peak[i]);

		if(fbkg<0||fbkg>1)
		{
			cout<<"======Input the bkg fraction error, Using the default one(0.5)"<<endl;
			fbkg=0.5;
		}
		val=(1-fbkg)*nEntries/nPeak;
		//log<<"number of signals:"<<val<<endl;
		areaPeak[i]=new RooRealVar(Form("areaPeak%d",i),Form("areaPeak%d",i),val,0,nEntries);
		listArea.add(*areaPeak[i]);
	}
	Double_t nBkg;
	nBkg=fbkg*nEntries;
	//log<<"number of background:"<<nBkg<<endl;
	areaBkg=new RooRealVar("areaBkg","areaBkg",nBkg,0,nEntries);
	listPeak.add(*bkg);
	listArea.add(*areaBkg);
	RooAddPdf *model=new RooAddPdf("model","fit model",listPeak,listArea);

	RooFitResult *fitRlt;
	if(hist)
	{
		fitRlt=model->fitTo(*dataHist,Extended(kTRUE),Minos(kTRUE),Save(kTRUE),NumCPU(7));
	} 
	else if(dataSet)
	{   
		fitRlt=model->fitTo(*dataSet,Extended(kTRUE),Minos(kTRUE),Save(kTRUE),NumCPU(7));
	} 

	/*
		 TCanvas *c1 = new TCanvas("c1","");
		 c1->SetFillColor(10);
		 c1->SetBorderMode(0);
		 c1->SetBorderSize(2);
		 c1->SetLeftMargin(0.1012121);
		 c1->SetRightMargin(0.02525252);
	//c1->SetTopMargin(0.02118644);
	c1->SetTopMargin(1);
	c1->SetFrameFillColor(0);
	c1->SetFrameBorderMode(0);
	c1->SetFrameBorderMode(0);
	c1->SetLineWidth(1);
	c1->SetTickx(1);
	c1->SetTicky(1);
	 */
	RooPlot* frame=x.frame();
	frame->GetXaxis()->SetTitle("M_{#mu#mu} (GeV/c^{2})");
	frame->GetXaxis()->CenterTitle();
	frame->GetYaxis()->CenterTitle();
	frame->SetTitle("Gaussian & Chebychev");
	if(hist)
	{
		dataHist->plotOn(frame);
	}
	else if(dataSet)
	{ 
		dataSet->plotOn(frame);
	} 
	//model->plotOn(frame,Components(*bkg),LineColor(kRed));
	model->plotOn(frame,Components("bkg"),LineColor(kRed),LineStyle(kDashed));
	model->plotOn(frame);
	//------get chi2------------------
	Int_t parsFit=(fitRlt->floatParsFinal()).getSize();
	Int_t nBins=frame->GetNbinsX();
	Double_t chi2Pre=frame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	Int_t ndof=nBins-parsFit;  //num of degree of freedom
	Double_t chi2=chi2Pre*ndof;

	frame->Draw();

	vector<double> sigBkg;
	if(getSB)
	{
		cout<<"==================================calculate S/B ratio=================================="<<endl;
		double Mean=0,Sigma=0,low0,high0, low1, high1;
		Mean=mean[0]->getVal();
		Sigma=sigma[0]->getVal();
		low0=Mean-2.*Sigma;
		high0=Mean+2*Sigma;
		low1=Mean-4.*Sigma;
		high1=Mean+4*Sigma;
		x.setRange("sig0",low0,high0);
		x.setRange("sig1",(double)x.getMin(),low1);
		x.setRange("sig2",high1,(double)x.getMax());

		RooAbsReal* frac0=model->createIntegral(x,NormSet(x),Range("sig0"));
		RooAbsReal* frac1=model->createIntegral(x,NormSet(x),Range("sig1"));
		RooAbsReal* frac2=model->createIntegral(x,NormSet(x),Range("sig2"));
		RooAbsReal* igPeak=peak[0]->createIntegral(x,NormSet(x),Range("sig0"));
		RooAbsReal* igBkg=bkg->createIntegral(x,NormSet(x),Range("sig0"));
		double ratioPeak=(double)areaPeak[0]->getVal()*(double)igPeak->getVal();
		double ratioBkg=(double)areaBkg->getVal()*(double)igBkg->getVal();
		double ratioSB0, ratioSB1;
		ratioSB0=ratioPeak/ratioBkg;                  //Signal-to-Backgrond S/BG
		ratioSB1=ratioPeak/sqrt(ratioBkg+ratioPeak);  //FOM(Figure Of Merit)
		double num0,num1;
		num0=(double)nEntries*(double)frac0->getVal();
		num1=(double)nEntries*((double)frac1->getVal()+(double)frac2->getVal());
		cout<<"============="<<igPeak->getVal()<<"=================="<<endl;
		cout<<"============="<<igBkg->getVal()<<"==================="<<endl;
		cout<<"========num Peak : "<<ratioPeak<<"==================="<<endl;
		cout<<"========num Bkg : "<<ratioBkg<<"==================="<<endl;
		cout<<"========num Model (+-2sigma) : "<<num0<<"=============="<<endl;
		cout<<"========num Model (out4sigma): "<<num1<<"=============="<<endl;
		cout<<"========S/B(2sigma): "<<ratioSB0<<"=============="<<endl;
		cout<<"========S/sqrt(S+B)(2sigma): "<<ratioSB1<<"=============="<<endl;

		sigBkg.push_back(Mean);
		sigBkg.push_back(Sigma);

		sigBkg.push_back(ratioSB0);
		sigBkg.push_back(ratioSB1);

		sigBkg.push_back(num0);
		sigBkg.push_back(num1);

		sigBkg.push_back((double)areaPeak[0]->getVal());
		sigBkg.push_back((double)areaBkg->getVal());

		//sigBkg.push_back(ratioPeak);
		//sigBkg.push_back(ratioBkg);

		cout<<"======================================================================================"<<endl;
	}
	//------------------------------------------------------------------------
	//----------------------------------Draw the Latex------------------------
	//------------------------------------------------------------------------
	if(showRlt)
	{
		Double_t left=0.68, top=0.80, textSize=0.032;
		Char_t formatStr[1000],str[1000];

		TLatex *latex=new TLatex();
		latex->SetTextFont(42);
		latex->SetNDC(kTRUE);
		latex->SetTextSize(textSize);
		Double_t step=textSize*1.3;
		latex->DrawLatex(left,top,Form("Entries=%ld",nEntries));
		top-=step;

		Int_t numDig=getNumDigs(chi2);
		numDig=2;
		//sprintf(formatStr,"#chi^{2}/nDOF=%%.%df",numDig);
		sprintf(formatStr,"#chi^{2}/ndof = %%.%df / %d",numDig,ndof);
		latex->DrawLatex(left,top,Form(formatStr,chi2));
		top-=step;

		/*sprintf(formatStr,"-LogL=%g",fitRlt->minNll());
			latex->DrawLatex(left,top,formatStr);
			top-=step;*/
		top=top-0.3*step;

		//-----Print the peak information-------
		for(int i=0; i<nPeak; i++)
		{
			Double_t val=mean[i]->getVal()*1000;  //GeV -> MeV
			Double_t error=mean[i]->getError()*1000;
			numDig=getNumDigs(error);
			sprintf(formatStr,"mean_{%d}=%%.%df #pm%%.%df\n",i,numDig,numDig);
			sprintf(str,formatStr,val,error);
			latex->SetTextColor(kRed);
			latex->DrawLatex(left,top,str);
			top-=step;

			val=sigma[i]->getVal()*1000;
			error=sigma[i]->getError()*1000;
			numDig=getNumDigs(error);
			sprintf(formatStr,"#sigma_{%d}=%%.%df #pm%%.%df\n",i,numDig,numDig);
			sprintf(str,formatStr,val,error);
			latex->SetTextColor(kRed);
			latex->DrawLatex(left,top,str);
			top-=step;

			val=areaPeak[i]->getVal();                             
			error=areaPeak[i]->getError();                         
			Double_t errL=areaPeak[i]->getAsymErrorLo();           
			Double_t errH=areaPeak[i]->getAsymErrorHi();           
			numDig=getNumDigs(error);                           
			sprintf(formatStr,"Nsig_{%d}=%%.%df #pm%%.%df\n",i,numDig,numDig);
			sprintf(str,formatStr,val,error);
			//sprintf(formatStr,"Nsig_{%d}=%%.%df ^{%%.%df}_{% %.%df}\n",i,numDig,numDig,numDig);         
			//sprintf(str,formatStr,val,errL,errH);         
			latex->SetTextColor(kRed);                          
			latex->DrawLatex(left,top,str);                     
			top-=step;                                          
		}

		//-----Print the background information----
		Double_t valBkg=areaBkg->getVal();                  
		Double_t error=areaBkg->getError();                         
		Double_t errL=areaBkg->getAsymErrorLo();                    
		Double_t errH=areaBkg->getAsymErrorHi();                    
		numDig=getNumDigs(error);                           
		sprintf(formatStr,"Nb_{%dChe}=%%.%df #pm%%.%df\n",nOrder,numDig,numDig); 
		sprintf(str,formatStr,valBkg,error);      
		//sprintf(formatStr,"Nb_{%iChe}=%%.%df ^{%%.%df}_{%%.%df}\n",nOrder,numDig,numDig,numDig);     
		//sprintf(str,formatStr,valBkg,errL,errH);      
		latex->SetTextColor(kBlue);                         
		latex->DrawLatex(left,top,str);                  
	}

	//c1->Modified();
	//c1->Update();
	//c1->Print("Unbinned.eps");
	RooArgSet* params = model->getVariables() ;
	cout<<endl;
	params->Print("v");
	//cout<<"Entries: "<<nEntries<<endl;
	//return fitRlt;
	return sigBkg;

}
