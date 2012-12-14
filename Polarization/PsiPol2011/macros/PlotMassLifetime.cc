#include "rootIncludes.inc"
#include "commonVar.h"
#include "RooUtils.h"

using namespace RooFit;

void plotMass(RooWorkspace *ws, int rapBin, int ptBin, int nState);
void plotMassLog(RooWorkspace *ws, int rapBin, int ptBin, int nState);
void plotLifeBg(RooWorkspace *ws, int rapBin, int ptBin, int nState);
void plotLifeSig(RooWorkspace *ws, int rapBin, int ptBin, int nState);
void plotLifeSig_linear(RooWorkspace *ws, int rapBin, int ptBin, int nState);
//==============================================

void PlotMassLifetime(const std::string &infilename, int rapBin, int ptBin, int nState, int Plotting){
	RooWorkspace* ws = getFromTFile<RooWorkspace>(infilename, "ws_masslifetime");

	switch (Plotting) {
		case 1:
			std::cout << ">>>>Plotting mass" << std::endl;
			plotMass(ws, rapBin, ptBin, nState);
			plotMassLog(ws, rapBin, ptBin, nState);
			std::cout << ">>>>Plotting lifetime sidebands" << std::endl;
			plotLifeBg(ws, rapBin, ptBin, nState);
			std::cout << ">>>>Plotting lifetime signal region" << std::endl;
			plotLifeSig(ws, rapBin, ptBin, nState);
			plotLifeSig_linear(ws, rapBin, ptBin, nState);
			break;
		case 2:
			std::cout << ">>>>Plotting mass" << std::endl;
			plotMass(ws, rapBin, ptBin, nState);
			plotMassLog(ws, rapBin, ptBin, nState);
			break;
		case 3:
			std::cout << ">>>>Plotting lifetime sidebands" << std::endl;
			plotLifeBg(ws, rapBin, ptBin, nState);
			break;
		case 4:
			std::cout << ">>>>Plotting lifetime signal region" << std::endl;
			plotLifeSig(ws, rapBin, ptBin, nState);
			//plotLifeSig_linear(ws, rapBin, ptBin, nState);
			break;
		default:
			std::cerr << "I dont know what do do with this value of Plotting" << std::endl;
	}
	delete ws;
}

//==============================================

//==============================================
void plotMass(RooWorkspace *ws, int rapBin, int ptBin, int nState){

	int  nbins=90; //0.005 bin size
	TGaxis::SetMaxDigits(3);

	if(nState == 4) nbins=90;
	if(nState == 5) nbins=120;

	RooRealVar *JpsiMass = ws->var("JpsiMass");
	assert( 0 != JpsiMass );
	//RooRealVar JpsiMass(*ws->var("JpsiMass"));

	RooPlot *massFrame = JpsiMass->frame(Bins(nbins));
	assert ( 0 != massFrame );
	massFrame->SetName(Form("mass_plot_rap%d_pt%d",rapBin,ptBin));
	massFrame->SetTitle("");
	massFrame->GetYaxis()->SetTitle("Events / 5 MeV");
	massFrame->GetYaxis()->SetTitleOffset(1.3);

	RooPlot *massFramePull = JpsiMass->frame(Bins(nbins));
	assert ( 0 != massFramePull );
	massFramePull->SetName(Form("pullmass_plot_rap%d_pt%d",rapBin,ptBin));
	massFramePull->SetTitle("");
	massFramePull->GetYaxis()->SetTitle("pull");
	massFramePull->GetXaxis()->SetTitleSize(0.08);
	massFramePull->GetYaxis()->SetTitleSize(0.08);
	massFramePull->GetXaxis()->SetLabelSize(0.08);
	massFramePull->GetYaxis()->SetLabelSize(0.08);
	//massFramePull->GetXaxis()->SetTitleOffset(1.0);
	massFramePull->GetYaxis()->SetTitleOffset(0.4);
	massFramePull->GetYaxis()->SetRangeUser(-5.5,5.5);

	RooAbsData *data= ws->data(Form("data_rap%d_pt%d",rapBin,ptBin));
	assert ( 0 != data );
	RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("m_fitresult_rap%d_pt%d",rapBin,ptBin)));
	assert ( 0 != fitRlt);

	double Mean = -1.0, MeanErr = -1.0;
	getVarFromWorkspace(ws, "CBmass", Mean, MeanErr);
	double Sigma = -1.0, SigmaErr = -1.0;
	getVarFromWorkspace(ws, "CBsigma", Sigma, SigmaErr);
	double Sigma2 = -1.0, Sigma2Err = -1.0;
	getVarFromWorkspace(ws, "CBsigma2", Sigma2, Sigma2Err);
	double Alpha = -1.0, AlphaErr = -1.0;
	getVarFromWorkspace(ws, "CBalpha", Alpha, AlphaErr);
	double cbN = -1.0, cbNErr = -1.0;
	getVarFromWorkspace(ws, "CBn", cbN, cbNErr);
	double lambda = -1.0, lambdaErr = -1.0;
	getVarFromWorkspace(ws, "bkgLambda", lambda, lambdaErr);
	double fracCB1 = -1.0, fracCB1Err = -1.0;
	getVarFromWorkspace(ws, "fracCB1", fracCB1, fracCB1Err);
	double fracBkg = -1.0, fracBkgErr = -1.0;
	getVarFromWorkspace(ws, "fracBkg", fracBkg, fracBkgErr);

	double SigmaWei = sqrt(pow(Sigma,2)*fracCB1+pow(Sigma2,2)*(1-fracCB1));
	double SigmaWeiErr =  (1./(2*SigmaWei))*
		sqrt(pow((pow(Sigma,2)-pow(Sigma2,2))*fracCB1Err,2) + pow(2*fracCB1*Sigma*SigmaErr,2) + pow(2*(1-fracCB1)*Sigma2*Sigma2Err,2));

	RooAbsPdf *massPdf = ws->pdf("massModel");
	assert ( 0 != massPdf );
	RooAbsPdf *bkgMassShape = ws->pdf("bkgMassShape");
	assert ( 0 != bkgMassShape );
	double sigMaxMass = Mean+SigmaWei*onia::nSigMass;
	double sigMinMass = Mean-SigmaWei*onia::nSigMass;
	double sbHighMass = Mean+SigmaWei*onia::nSigBkgHigh;
	double sbLowMass =  Mean-SigmaWei*onia::nSigBkgLow;
	cout<<"sigMaxMass: "<<sigMaxMass<<endl;
	cout<<"sigMinMass: "<<sigMinMass<<endl;
	cout<<"sbHighMass: "<<sbHighMass<<endl;
	cout<<"sbLowMass: "<<sbLowMass<<endl;

	int nEntries = data->numEntries();
	JpsiMass->setRange("SigRegion",sigMinMass,sigMaxMass);

	RooRealVar* fracFull3Sig=(RooRealVar*)massPdf->createIntegral(*JpsiMass,NormSet(*JpsiMass),Range("SigRegion"));
	RooRealVar* fracBkg3Sig=(RooRealVar*)bkgMassShape->createIntegral(*JpsiMass,NormSet(*JpsiMass),Range("SigRegion"));
	double evtFull3Sig = nEntries*fracFull3Sig->getVal();
	double evtBkg3Sig = nEntries*fracBkg*fracBkg3Sig->getVal();
	double BkgRatio3Sig = evtBkg3Sig/evtFull3Sig;

	cout<<"fracFull3Sig: "<<fracFull3Sig->getVal()<<endl;
	cout<<"fracBkg3Sig: "<<fracBkg3Sig->getVal()<<endl;
	cout<<"evtFull3Sig: "<<evtFull3Sig<<endl;
	cout<<"evtBkg3Sig: "<<evtBkg3Sig<<endl;
	cout<<"BkgRatio3Sig: "<<BkgRatio3Sig<<endl;

	data->plotOn(massFrame,MarkerSize(0.8));
	massPdf->plotOn(massFrame,
			LineWidth(2),
			ProjWData(*data));
	//------get chi2------------SHOULD DONE after PLOTTING------
	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_Mass=massFrame->GetNbinsX();
	double chi2Pre_Mass=massFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_Mass=nBins_Mass-parsFit;  //num of degree of freedom
	double chi2_Mass=chi2Pre_Mass*ndof_Mass;

	RooHist* hpull_mass = massFrame->pullHist() ;
	hpull_mass->SetMarkerSize(0.8);
	for(int i=0;i<hpull_mass->GetN();i++){
		hpull_mass->SetPointEYlow(i,0.);
		hpull_mass->SetPointEYhigh(i,0.);
	}
	massFramePull->addPlotable(hpull_mass,"P");

	massPdf->plotOn(massFrame,
			Components("bkgMassShape"),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2),
			ProjWData(*data));

	//if(fracCB1 < 0.5){
	if(Sigma > Sigma2){

		fracCB1 = 1-fracCB1;
		double temp = 0.;
		temp = Sigma; Sigma = Sigma2; Sigma2 = temp;
		temp = SigmaErr; SigmaErr = Sigma2Err; Sigma2Err = temp;
		massPdf->plotOn(massFrame,
				Components("massCBShape"),
				LineStyle(2),
				LineColor(kGreen),
				LineWidth(2),
				ProjWData(*data));
		massPdf->plotOn(massFrame,
				Components("massCBShape2"),
				LineStyle(2),
				LineColor(kRed),
				LineWidth(2),
				ProjWData(*data));

	}  else{
		massPdf->plotOn(massFrame,
				Components("massCBShape"),
				LineStyle(2),
				LineColor(kRed),
				LineWidth(2),
				ProjWData(*data));
		massPdf->plotOn(massFrame,
				Components("massCBShape2"),
				LineStyle(2),
				LineColor(kGreen),
				LineWidth(2),
				ProjWData(*data));
	}


	double minY = 0.;

	double maxY = 0.;
	if(nState == 4) maxY = massFrame->GetMaximum()*0.3;
	if(nState == 5) maxY = massFrame->GetMaximum()*0.4;
	double lineWidth = 2.0;
	TLine *lineSBLow = new TLine(sbLowMass, minY, sbLowMass, maxY);
	TLine *lineSBHigh = new TLine(sbHighMass, minY, sbHighMass, maxY);
	TLine *lineSigLow = new TLine(sigMinMass, minY, sigMinMass, maxY);
	TLine *lineSigHigh = new TLine(sigMaxMass, minY, sigMaxMass, maxY);
	lineSBLow->SetLineWidth(lineWidth);lineSBHigh->SetLineWidth(lineWidth);
	lineSigLow->SetLineWidth(lineWidth);lineSigHigh->SetLineWidth(lineWidth);
	lineSBLow->SetLineColor(kBlue);lineSBHigh->SetLineColor(kBlue);
	lineSigLow->SetLineColor(kRed);lineSigHigh->SetLineColor(kRed);
	lineSBLow->SetLineStyle(7);lineSBHigh->SetLineStyle(7);
	lineSigLow->SetLineStyle(5);lineSigHigh->SetLineStyle(5);

	TH1* legendBlue = data->createHistogram("legendBlue",*JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = data->createHistogram("legendBlueDash",*JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = data->createHistogram("legendRed",*JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = data->createHistogram("legendBlack",*JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = data->createHistogram("legendGreen",*JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(2) ; legendGreen->SetLineWidth(2) ;
	TH1* legendGreenDash = data->createHistogram("legendGreenDash",*JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
	TH1* legendPink = data->createHistogram("legendPink",*JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	//TLegend* MassLegend=new TLegend(0.65,0.7,0.88,0.88);
	TLegend* MassLegend=new TLegend(0.65,0.5,0.88,0.7);
	MassLegend->SetFillColor(kWhite);
	MassLegend->SetTextFont(42);
	MassLegend->SetTextSize(0.025);
	MassLegend->SetBorderSize(0.);
	MassLegend->AddEntry(legendBlue,"sum","l");
	MassLegend->AddEntry(legendRed,"signal CB_{1}","l");
	MassLegend->AddEntry(legendGreen,"signal CB_{2}","l");
	MassLegend->AddEntry(legendPink,"background","l");

	double left=0.65, top=0.65, textSize=0.025;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	TCanvas *c1=new TCanvas("c1","",800,700);
	c1->SetTickx();
	c1->SetTicky();

	massFrame->Draw(); MassLegend->Draw();
	lineSBLow->Draw("same"); lineSBHigh->Draw("same"); lineSigLow->Draw("same"); lineSigHigh->Draw("same");
	top=0.85; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5)
		latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);


	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));

	top-=step;
	top-=0.5*step;
	latex->DrawLatex(left,top,Form("n_{CB}  =  %.1f",cbN));

	left=0.15; top=0.85; textSize=0.020;
	latex->SetTextSize(textSize);
	latex->DrawLatex(left,top,Form("#chi^{2}/ndof = %.2f / %d", chi2_Mass, ndof_Mass));
	top-=step;
	latex->DrawLatex(left,top,Form("mean   =  %.3f #pm %.3f MeV",Mean*1000, MeanErr*1000));
	top-=step;
	latex->DrawLatex(left,top,Form("#sigma_{1}  =  %.3f #pm %.3f MeV",Sigma*1000, SigmaErr*1000));
	top-=step;
	latex->DrawLatex(left,top,Form("#sigma_{2}  =  %.3f #pm %.3f MeV",Sigma2*1000, Sigma2Err*1000));
	top-=step;
	latex->DrawLatex(left,top,Form("#alpha  =  %.3f #pm %.3f",Alpha, AlphaErr));
	top-=step;
	latex->DrawLatex(left,top,Form("#lambda =  %.3f #pm %.3f",lambda, lambdaErr));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{CB_{1}}  =  %.3f #pm %.3f",fracCB1, fracCB1Err));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Bg}  =  %.3f #pm %.3f",fracBkg, fracBkgErr));
	top-=step;
	top-=step;
	top-=step;
	latex->DrawLatex(left,top,Form("effective #sigma =  %.3f MeV",SigmaWei*1000));
	top-=step;
	top-=0.5*step;
	latex->SetTextSize(textSize);
	latex->DrawLatex(left,top,Form("#frac{B}{B+S} (#pm3#sigma)  =  %.3f",BkgRatio3Sig));

	std::stringstream saveMass;
	saveMass << "Fit/mass_rap" << rapBin << "_pt" << ptBin << ".pdf";
	c1->SaveAs(saveMass.str().c_str());

	delete c1;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendGreenDash;
	delete legendPink;
	return;
}

//==============================================
void plotMassLog(RooWorkspace *ws, int rapBin, int ptBin, int nState){

	int  nbins=90; //0.005 bin size
	TGaxis::SetMaxDigits(3);

	if(nState == 4) nbins=90;
	if(nState == 5) nbins=120;
	RooRealVar JpsiMass(*ws->var("JpsiMass"));

	RooPlot *massFrame=((RooRealVar*)ws->var("JpsiMass"))->frame(Bins(nbins));
	massFrame->SetName(Form("mass_plot_rap%d_pt%d",rapBin,ptBin));
	massFrame->SetTitle("");
	massFrame->GetYaxis()->SetTitle("Events / 5 MeV");
	massFrame->GetYaxis()->SetTitleOffset(1.3);

	RooPlot *massFramePull=((RooRealVar*)ws->var("JpsiMass"))->frame(Bins(nbins));
	massFramePull->SetName(Form("pullmass_plot_rap%d_pt%d",rapBin,ptBin));
	massFramePull->SetTitle("");
	massFramePull->GetYaxis()->SetTitle("pull");
	massFramePull->GetXaxis()->SetTitleSize(0.08);
	massFramePull->GetYaxis()->SetTitleSize(0.08);
	massFramePull->GetXaxis()->SetLabelSize(0.08);
	massFramePull->GetYaxis()->SetLabelSize(0.08);
	//massFramePull->GetXaxis()->SetTitleOffset(1.0);
	massFramePull->GetYaxis()->SetTitleOffset(0.4);
	massFramePull->GetYaxis()->SetRangeUser(-5.5,5.5);

	RooDataSet *data=(RooDataSet *)ws->data(Form("data_rap%d_pt%d",rapBin,ptBin));
	RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("m_fitresult_rap%d_pt%d",rapBin,ptBin));

	RooRealVar *CBmass=(RooRealVar *)ws->var("CBmass");
	RooRealVar *CBsigma=(RooRealVar *)ws->var("CBsigma");
	RooRealVar *CBsigma2=(RooRealVar *)ws->var("CBsigma2");
	RooRealVar *CBalpha=(RooRealVar *)ws->var("CBalpha");
	RooRealVar *CBn=(RooRealVar *)ws->var("CBn");
	RooRealVar *bkgLambda=(RooRealVar *)ws->var("bkgLambda");
	RooRealVar *fracCB1_=(RooRealVar *)ws->var("fracCB1");
	RooRealVar *fracBkg_=(RooRealVar *)ws->var("fracBkg");

	double Mean = CBmass->getVal();
	double MeanErr = CBmass->getError();
	double Sigma = CBsigma->getVal();
	double SigmaErr = CBsigma->getError();
	double Sigma2 = CBsigma2->getVal();
	double Sigma2Err = CBsigma2->getError();
	double Alpha = CBalpha->getVal();
	double AlphaErr = CBalpha->getError();
	double cbN = CBn->getVal();
	double cbNErr = CBn->getError();
	double lambda = bkgLambda->getVal();
	double lambdaErr = bkgLambda->getError();
	double fracCB1 = fracCB1_->getVal();
	double fracCB1Err = fracCB1_->getError();
	double fracBkg = fracBkg_->getVal();
	double fracBkgErr = fracBkg_->getError();

	double SigmaWei = sqrt(pow(Sigma,2)*fracCB1+pow(Sigma2,2)*(1-fracCB1));
	double SigmaWeiErr =  (1./(2*SigmaWei))*
		sqrt(pow((pow(Sigma,2)-pow(Sigma2,2))*fracCB1Err,2) + pow(2*fracCB1*Sigma*SigmaErr,2) + pow(2*(1-fracCB1)*Sigma2*Sigma2Err,2));
	cout<<">>=====Mean: "<<Mean<<endl;
	cout<<">>=====Sigma: "<<Sigma<<endl;
	cout<<">>=====Sigma2: "<<Sigma2<<endl;
	cout<<">>=====Alpha: "<<Alpha<<endl;
	cout<<">>=====cbN: "<<cbN<<endl;
	cout<<">>=====lambda: "<<lambda<<endl;
	cout<<">>=====fracCB1: "<<fracCB1<<endl;
	cout<<">>=====fracBkg: "<<fracBkg<<endl;
	cout<<">>=====SigmaWei: "<<SigmaWei<<endl;

	RooAddPdf *massPdf = (RooAddPdf*)ws->pdf("massModel");
	RooAddPdf *bkgMassShape = (RooAddPdf*)ws->pdf("bkgMassShape");

	double sigMaxMass = Mean+SigmaWei*onia::nSigMass;
	double sigMinMass = Mean-SigmaWei*onia::nSigMass;
	double sbHighMass = Mean+SigmaWei*onia::nSigBkgHigh;
	double sbLowMass =  Mean-SigmaWei*onia::nSigBkgLow;
	cout<<"sigMaxMass: "<<sigMaxMass<<endl;
	cout<<"sigMinMass: "<<sigMinMass<<endl;
	cout<<"sbHighMass: "<<sbHighMass<<endl;
	cout<<"sbLowMass: "<<sbLowMass<<endl;

	int nEntries = data->numEntries();
	JpsiMass.setRange("SigRegion",sigMinMass,sigMaxMass);

	RooRealVar* fracFull3Sig=(RooRealVar*)massPdf->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SigRegion"));
	RooRealVar* fracBkg3Sig=(RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SigRegion"));
	double evtFull3Sig = nEntries*fracFull3Sig->getVal();
	double evtBkg3Sig = nEntries*fracBkg*fracBkg3Sig->getVal();
	double BkgRatio3Sig = evtBkg3Sig/evtFull3Sig;
	cout<<"fracFull3Sig: "<<fracFull3Sig->getVal()<<endl;
	cout<<"fracBkg3Sig: "<<fracBkg3Sig->getVal()<<endl;
	cout<<"evtFull3Sig: "<<evtFull3Sig<<endl;
	cout<<"evtBkg3Sig: "<<evtBkg3Sig<<endl;
	cout<<"BkgRatio3Sig: "<<BkgRatio3Sig<<endl;
	data->plotOn(massFrame,MarkerSize(0.8));
	massPdf->plotOn(massFrame,
			LineWidth(2),
			ProjWData(*data));
	//------get chi2------------SHOULD DONE after PLOTTING------
	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_Mass=massFrame->GetNbinsX();
	double chi2Pre_Mass=massFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_Mass=nBins_Mass-parsFit;  //num of degree of freedom
	double chi2_Mass=chi2Pre_Mass*ndof_Mass;

	RooHist* hpull_mass = massFrame->pullHist() ;
	hpull_mass->SetMarkerSize(0.8);
	for(int i=0;i<hpull_mass->GetN();i++){
		hpull_mass->SetPointEYlow(i,0.);
		hpull_mass->SetPointEYhigh(i,0.);
	}
	massFramePull->addPlotable(hpull_mass,"P");

	massPdf->plotOn(massFrame,
			Components("bkgMassShape"),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2),
			ProjWData(*data));

	//if(fracCB1 < 0.5){
	if(Sigma > Sigma2){

		fracCB1 = 1-fracCB1;
		double temp = 0.;
		temp = Sigma; Sigma = Sigma2; Sigma2 = temp;
		temp = SigmaErr; SigmaErr = Sigma2Err; Sigma2Err = temp;
		massPdf->plotOn(massFrame,
				Components("massCBShape"),
				LineStyle(2),
				LineColor(kGreen),
				LineWidth(2),
				ProjWData(*data));
		massPdf->plotOn(massFrame,
				Components("massCBShape2"),
				LineStyle(2),
				LineColor(kRed),
				LineWidth(2),
				ProjWData(*data));

	}  else{
		massPdf->plotOn(massFrame,
				Components("massCBShape"),
				LineStyle(2),
				LineColor(kRed),
				LineWidth(2),
				ProjWData(*data));
		massPdf->plotOn(massFrame,
				Components("massCBShape2"),
				LineStyle(2),
				LineColor(kGreen),
				LineWidth(2),
				ProjWData(*data));
	}
	double minY = 0.;
	double maxY = 0.;
	if(nState == 4) maxY = massFrame->GetMaximum()*0.1;
	if(nState == 5) maxY = massFrame->GetMaximum()*0.3;
	double lineWidth = 2.0;
	TLine *lineSBLow = new TLine(sbLowMass, minY, sbLowMass, maxY);
	TLine *lineSBHigh = new TLine(sbHighMass, minY, sbHighMass, maxY);
	TLine *lineSigLow = new TLine(sigMinMass, minY, sigMinMass, maxY);
	TLine *lineSigHigh = new TLine(sigMaxMass, minY, sigMaxMass, maxY);
	lineSBLow->SetLineWidth(lineWidth);lineSBHigh->SetLineWidth(lineWidth);
	lineSigLow->SetLineWidth(lineWidth);lineSigHigh->SetLineWidth(lineWidth);
	lineSBLow->SetLineColor(kBlue);lineSBHigh->SetLineColor(kBlue);
	lineSigLow->SetLineColor(kRed);lineSigHigh->SetLineColor(kRed);
	lineSBLow->SetLineStyle(7);lineSBHigh->SetLineStyle(7);
	lineSigLow->SetLineStyle(5);lineSigHigh->SetLineStyle(5);

	double Ymax = massFrame->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;

	if(nState == 4) massFrame->SetMaximum(90000.);
	if(nState == 5) massFrame->SetMaximum(7000.);
	massFrame->SetMinimum(2.);
	TH1* legendBlue = data->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = data->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = data->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = data->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = data->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(2) ; legendGreen->SetLineWidth(2) ;
	TH1* legendPink = data->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	TLegend* MassLegend=new TLegend(0.65,0.75,0.88,0.88);
	MassLegend->SetFillColor(kWhite);
	MassLegend->SetTextFont(42);
	MassLegend->SetTextSize(0.025);
	MassLegend->SetBorderSize(0.);
	MassLegend->AddEntry(legendBlue,"sum","l");
	MassLegend->AddEntry(legendRed,"signal CB_{1}","l");
	MassLegend->AddEntry(legendGreen,"signal CB_{2}","l");
	MassLegend->AddEntry(legendPink,"background","l");

	double left=0.15, top=0.85, textSize=0.025;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	TCanvas *c1=new TCanvas("c1","",800,700);
	c1->SetTickx();
	c1->SetTicky();

	c1->cd();

	c1->SetLogy(1);
	massFrame->Draw(); MassLegend->Draw();
	lineSBLow->Draw("same"); lineSBHigh->Draw("same"); lineSigLow->Draw("same"); lineSigHigh->Draw("same");
	textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5)
		latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));

	top-=step;
	top-=step;

	std::stringstream saveMasslog;
	saveMasslog << "Fit/mass_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
	c1->SaveAs(saveMasslog.str().c_str());

	delete c1;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendPink;
	return;
}

//==============================================
void plotLifeBg(RooWorkspace *ws, int rapBin, int ptBin, int nState){

	int nbins=140; //bin size 0.025 mm

	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	RooRealVar *JpsictErr = ws->var("JpsictErr");

	RooPlot *ctauFrameBkgSBL=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(-1., 2.5));
	ctauFrameBkgSBL->SetName(Form("ctaubkg_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrameBkgSBL->SetYTitle("Events per 25 micron");
	ctauFrameBkgSBL->SetTitle("");

	RooPlot *ctauFrameBkgSBLPull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(-1., 2.5));
	ctauFrameBkgSBLPull->SetName(Form("pullctaubkg_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrameBkgSBLPull->SetTitle("");
	if(nState == 4) ctauFrameBkgSBLPull->GetXaxis()->SetTitle("c_{#tau}^{J/#psi} [mm]");
	if(nState == 5) ctauFrameBkgSBLPull->GetXaxis()->SetTitle("c_{#tau}^{#psi(2S)} [mm]");
	ctauFrameBkgSBLPull->GetYaxis()->SetTitle("pull");
	ctauFrameBkgSBLPull->GetXaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBLPull->GetYaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBLPull->GetXaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBLPull->GetYaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBLPull->GetYaxis()->SetTitleOffset(0.4);
	ctauFrameBkgSBLPull->GetYaxis()->SetRangeUser(-5.5,5.5);

	RooPlot *ctauFrameBkgSBR=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(-1., 2.5));
	ctauFrameBkgSBR->SetName(Form("ctaubkg_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrameBkgSBR->SetYTitle("Events per 25 micron");
	ctauFrameBkgSBR->SetTitle("");

	RooPlot *ctauFrameBkgSBRPull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(-1., 2.5));
	ctauFrameBkgSBRPull->SetName(Form("pullctaubkg_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrameBkgSBRPull->SetTitle("");
	if(nState == 4) ctauFrameBkgSBRPull->GetXaxis()->SetTitle("c_{#tau}^{J/#psi} [mm]");
	if(nState == 5) ctauFrameBkgSBRPull->GetXaxis()->SetTitle("c_{#tau}^{#psi(2S)} [mm]");
	ctauFrameBkgSBRPull->GetYaxis()->SetTitle("pull");
	ctauFrameBkgSBRPull->GetXaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBRPull->GetYaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBRPull->GetXaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBRPull->GetYaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBRPull->GetYaxis()->SetTitleOffset(0.4);
	ctauFrameBkgSBRPull->GetYaxis()->SetRangeUser(-5.5,5.5);

	RooDataSet *dataSBL=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_SBL",rapBin,ptBin));
	RooDataSet *dataSBR=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_SBR",rapBin,ptBin));

	RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("l_fitresult_rap%d_pt%d",rapBin,ptBin));
	if(!fitRlt) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; return; }
	fitRlt->Print();

	RooRealVar *promptMean_=(RooRealVar*)ws->var("promptMean");
	RooRealVar *ctResolution_=(RooRealVar*)ws->var("ctResolution");
	RooRealVar *ctResolution2_=(RooRealVar*)ws->var("ctResolution2");
	RooRealVar *fracGauss2_=(RooRealVar*)ws->var("fracGauss2");

	RooRealVar *bkgTauSSD_SBL_=(RooRealVar*)ws->var("bkgTauSSD_SBL");
	RooRealVar *bkgTauFD_=(RooRealVar*)ws->var("bkgTauFD");
	RooRealVar *bkgTauDSD_=(RooRealVar*)ws->var("bkgTauDSD");

	RooRealVar *bkgTauSSD_SBR_=(RooRealVar*)ws->var("bkgTauSSD_SBR");
	//RooRealVar *bkgTauFD_=(RooRealVar*)ws->var("bkgTauFD");
	//RooRealVar *bkgTauDSD_=(RooRealVar*)ws->var("bkgTauDSD");

	RooRealVar *fBkgSSDR_SBL_ = (RooRealVar*)ws->var("fBkgSSDR_SBL");
	RooRealVar *fBkgDSD_SBL_ = (RooRealVar*)ws->var("fBkgDSD_SBL");
	RooRealVar *fBkgSSDR_SBR_ = (RooRealVar*)ws->var("fBkgSSDR_SBR");
	RooRealVar *fBkgDSD_SBR_ = (RooRealVar*)ws->var("fBkgDSD_SBR");

	RooRealVar *fBkg_ = (RooRealVar*)ws->var("fBkg");

	double promptMean = promptMean_->getVal();
	double promptMeanErr = promptMean_->getError();
	double promptCtRe = ctResolution_->getVal();
	double promptCtReErr = ctResolution_->getError();
	double promptCtRe2 = ctResolution2_->getVal();
	double promptCtRe2Err = ctResolution2_->getError();
	double fracGauss2 = fracGauss2_->getVal();
	double fracGauss2Err = fracGauss2_->getError();

	double bkgTauSSD_SBL = bkgTauSSD_SBL_->getVal();
	double bkgTauSSD_SBLErr = bkgTauSSD_SBL_->getError();
	double bkgTauFD = bkgTauFD_->getVal();
	double bkgTauFDErr = bkgTauFD_->getError();
	double bkgTauDSD = bkgTauDSD_->getVal();
	double bkgTauDSDErr = bkgTauDSD_->getError();

	double bkgTauSSD_SBR = bkgTauSSD_SBR_->getVal();
	double bkgTauSSD_SBRErr = bkgTauSSD_SBR_->getError();


	double fBkg = fBkg_->getVal();
	double fBkgErr = fBkg_->getError();
	cout<<"fBkg: "<<fBkg<<endl;
	cout<<"fBkgErr: "<<fBkgErr<<endl;

	double fBkgSSDR_SBL = fBkgSSDR_SBL_->getVal();
	double fBkgSSDRErr_SBL = fBkgSSDR_SBL_->getError();
	double fBkgDSD_SBL = fBkgDSD_SBL_->getVal();
	double fBkgDSDErr_SBL = fBkgDSD_SBL_->getError();
	double fBkgSSDL_SBL = 1-fBkgSSDR_SBL-fBkgDSD_SBL;

	double fBkgSSDR_SBR = fBkgSSDR_SBR_->getVal();
	double fBkgSSDRErr_SBR = fBkgSSDR_SBR_->getError();
	double fBkgDSD_SBR = fBkgDSD_SBR_->getVal();
	double fBkgDSDErr_SBR = fBkgDSD_SBR_->getError();
	double fBkgSSDL_SBR = 1-fBkgSSDR_SBR-fBkgDSD_SBR;

	cout<<"promptMean: "<<promptMean<<endl;
	cout<<"promptCtRe: "<<promptCtRe<<endl;
	cout<<"fBkgSSDR_SBL: "<<fBkgSSDR_SBL<<endl;
	cout<<"fBkgSSDL_SBL: "<<fBkgSSDL_SBL<<endl;
	cout<<"fBkgSSDR_SBR: "<<fBkgSSDR_SBR<<endl;
	cout<<"fBkgSSDL_SBR: "<<fBkgSSDL_SBR<<endl;

	RooAddPdf *ModelLifeSBL = (RooAddPdf*)ws->pdf("backgroundlifetimeL");
	RooAddPdf *ModelLifeSBR = (RooAddPdf*)ws->pdf("backgroundlifetimeR");

	int parsFit;

	//ploting background SBL
	dataSBL->plotOn(ctauFrameBkgSBL,MarkerSize(0.8));
	ModelLifeSBL->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			LineWidth(2),
			NumCPU(1)
			);

	parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_LSBL=ctauFrameBkgSBL->GetNbinsX();
	double chi2Pre_LSBL=ctauFrameBkgSBL->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_LSBL=nBins_LSBL-parsFit;  //num of degree of freedom
	double chi2_LSBL=chi2Pre_LSBL*ndof_LSBL;

	RooHist* hpull_ctauBkgSBL = ctauFrameBkgSBL->pullHist() ;
	hpull_ctauBkgSBL->SetMarkerSize(0.8);
	for(int i=0;i<hpull_ctauBkgSBL->GetN();i++){
		hpull_ctauBkgSBL->SetPointEYlow(i,0.);
		hpull_ctauBkgSBL->SetPointEYhigh(i,0.);
	}
	cout<<"hpull_ctauBkgSBL->GetN(): "<<hpull_ctauBkgSBL->GetN()<<endl;

	ctauFrameBkgSBLPull->addPlotable(hpull_ctauBkgSBL,"P");

	ModelLifeSBL->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			Components("backgroundSSD_SBL"),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2),
			NumCPU(1)
			);

	ModelLifeSBL->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			Components("backgroundFD"),
			LineStyle(2),
			LineColor(kGreen),
			LineWidth(2),
			NumCPU(1)
			);

	ModelLifeSBL->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			Components("backgroundDSD"),
			LineStyle(1),
			LineColor(kBlack),
			LineWidth(2),
			NumCPU(1)
			);

	double Ymax = ctauFrameBkgSBL->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;

	ctauFrameBkgSBL->SetMaximum(5*Ymax);
	ctauFrameBkgSBL->SetMinimum(0.001*Ymax);

	dataSBR->plotOn(ctauFrameBkgSBR,MarkerSize(0.8));
	ModelLifeSBR->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			LineWidth(2),
			NumCPU(1)
			);

	parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_LSBR=ctauFrameBkgSBR->GetNbinsX();
	double chi2Pre_LSBR=ctauFrameBkgSBR->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_LSBR=nBins_LSBR-parsFit;  //num of degree of freedom
	double chi2_LSBR=chi2Pre_LSBR*ndof_LSBR;

	RooHist* hpull_ctauBkgSBR = ctauFrameBkgSBR->pullHist() ;
	hpull_ctauBkgSBR->SetMarkerSize(0.8);
	for(int i=0;i<hpull_ctauBkgSBR->GetN();i++){
		hpull_ctauBkgSBR->SetPointEYlow(i,0.);
		hpull_ctauBkgSBR->SetPointEYhigh(i,0.);
	}
	ctauFrameBkgSBRPull->addPlotable(hpull_ctauBkgSBR,"P");

	ModelLifeSBR->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			Components("backgroundSSD_SBR"),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2),
			NumCPU(1)
			);

	ModelLifeSBR->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			Components("backgroundFD"),
			LineStyle(2),
			LineColor(kGreen),
			LineWidth(2),
			NumCPU(1)
			);

	ModelLifeSBR->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			Components("backgroundDSD"),
			LineStyle(1),
			LineColor(kBlack),
			LineWidth(2),
			NumCPU(1)
			);

	Ymax = ctauFrameBkgSBR->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;
	ctauFrameBkgSBR->SetMaximum(5*Ymax);
	ctauFrameBkgSBR->SetMinimum(0.001*Ymax);

	TH1* legendBlue = dataSBL->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = dataSBL->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = dataSBL->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = dataSBL->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = dataSBL->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
	TH1* legendGreenDash = dataSBL->createHistogram("legendGreenDash",JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
	TH1* legendPink = dataSBL->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	//TLegend* LifetimeLegendBkgSBL=new TLegend(0.65,0.75,0.88,0.88);
	TLegend* LifetimeLegendBkgSBL=new TLegend(0.12,0.75,0.25,0.88);
	LifetimeLegendBkgSBL->SetFillColor(kWhite);
	LifetimeLegendBkgSBL->SetTextFont(42);
	LifetimeLegendBkgSBL->SetTextSize(0.025);
	LifetimeLegendBkgSBL->SetBorderSize(0.);
	LifetimeLegendBkgSBL->AddEntry(legendBlue,"sum","l");
	LifetimeLegendBkgSBL->AddEntry(legendGreenDash,"SSL","l");
	LifetimeLegendBkgSBL->AddEntry(legendRed,"SSR ","l");
	LifetimeLegendBkgSBL->AddEntry(legendBlack,"DS","l");

	TLegend* LifetimeLegendBkgSBR=new TLegend(0.12,0.75,0.25,0.88);
	LifetimeLegendBkgSBR->SetFillColor(kWhite);
	LifetimeLegendBkgSBR->SetTextFont(42);
	LifetimeLegendBkgSBR->SetTextSize(0.025);
	LifetimeLegendBkgSBR->SetBorderSize(0.);
	LifetimeLegendBkgSBR->AddEntry(legendBlue,"sum","l");
	LifetimeLegendBkgSBR->AddEntry(legendGreenDash,"SSL","l");
	LifetimeLegendBkgSBR->AddEntry(legendRed,"SSR ","l");
	LifetimeLegendBkgSBR->AddEntry(legendBlack,"DS","l");

	//TLegend* LifetimeLegendSig=new TLegend(0.65,0.7,0.88,0.88);
	TLegend* LifetimeLegendSig=new TLegend(0.12,0.75,0.25,0.88);
	LifetimeLegendSig->SetFillColor(kWhite);
	LifetimeLegendSig->SetTextFont(42);
	LifetimeLegendSig->SetTextSize(0.025);
	LifetimeLegendSig->SetBorderSize(0.);
	LifetimeLegendSig->AddEntry(legendBlue,"sum","l");
	if(nState == 4){
		LifetimeLegendSig->AddEntry(legendBlueDash,"prompt J/#psi","l");
		LifetimeLegendSig->AddEntry(legendRed,"non-Prompt J/#psi","l");
	}
	if(nState == 5){
		LifetimeLegendSig->AddEntry(legendBlueDash,"prompt #psi(2S)","l");
		LifetimeLegendSig->AddEntry(legendRed,"non-Prompt #psi(2S)","l");
	}
	LifetimeLegendSig->AddEntry(legendPink,"background","l");

	double left=0.6, top=0.85, textSize=0.025;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	TCanvas *c1=new TCanvas("c1","",800,900);
	c1->SetTickx();
	c1->SetTicky();

	c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
	pad1->SetGridy();
	pad1->SetLeftMargin(0.1);
	pad1->SetBottomMargin(0.2);
	pad1->Draw();
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
	pad2->SetLeftMargin(0.1);
	pad2->SetBottomMargin(0.1);
	pad2->Draw();

	//Bkg SBL
	pad2->cd(0);pad2->SetLogy(1);
	ctauFrameBkgSBL->Draw();
	LifetimeLegendBkgSBL->Draw();
	left=0.4; top=0.85; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));

	top-=step;
	latex->SetTextColor(kRed);
	latex->DrawLatex(left,top,"Left Mass SideBand");
	latex->SetTextColor(kBlack);
	top-=step;
	latex->DrawLatex(left,top,Form("mean  =  %.1f",promptMean));
	top-=step;
	latex->DrawLatex(left,top,Form("resolution1 =  %.1f",promptCtRe));
	top-=step;
	latex->DrawLatex(left,top,Form("resolution2  =  %.1f",promptCtRe2));
	top-=step;
	if(ptBin>7)
		latex->DrawLatex(left,top,Form("frac_{SSR}   =  %.2f",fBkgSSDR_SBL));


	left=0.6; top=0.85;
	latex->DrawLatex(left,top,Form("#chi^{2}/ndof   =  %.2f/%d",chi2_LSBL,ndof_LSBL));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{DS}   =  %.3f #pm %.3f",fBkgDSD_SBL, fBkgDSDErr_SBL));
	top-=step;
	if(ptBin<8){
		latex->DrawLatex(left,top,Form("frac_{SSR}   =  %.3f #pm %.3f",fBkgSSDR_SBL, fBkgSSDRErr_SBL));
		top-=step;
	}
	latex->DrawLatex(left,top,Form("#tau_{SSR}   =  %.3f #pm %.3f",bkgTauSSD_SBL, bkgTauSSD_SBLErr));
	top-=step;
	top-=step;

	latex->DrawLatex(left,top,Form("frac_{Gauss2}   =  %.3f #pm %.3f",fracGauss2, fracGauss2Err));
	top-=step;
	latex->DrawLatex(left,top,Form("#tau_{DS}   =  %.3f #pm %.3f",bkgTauDSD, bkgTauDSDErr));
	top-=step;
	latex->DrawLatex(left,top,Form("#tau_{SSL}  =  %.3f #pm %.3f",bkgTauFD, bkgTauFDErr));

	pad1->cd(0); pad1->SetLogy(0);
	ctauFrameBkgSBLPull->Draw();

	std::stringstream savectLSB;
	savectLSB << "Fit/ct_LSB_rap" << rapBin << "_pt" << ptBin << ".pdf";
	c1->SaveAs(savectLSB.str().c_str());

	//Bkg SBR
	pad2->cd(0);pad2->SetLogy(1);
	ctauFrameBkgSBR->Draw();
	LifetimeLegendBkgSBR->Draw();
	left=0.4; top=0.85; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));

	top-=step;
	latex->SetTextColor(kRed);
	latex->DrawLatex(left,top,"Right Mass SideBand");
	latex->SetTextColor(kBlack);
	top-=step;

	left=0.6; top=0.85;
	latex->DrawLatex(left,top,Form("#chi^{2}/ndof   =  %.2f/%d",chi2_LSBR,ndof_LSBR));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{DS}   =  %.3f #pm %.3f",fBkgDSD_SBR, fBkgDSDErr_SBR));
	top-=step;
	if(ptBin<8){
		latex->DrawLatex(left,top,Form("frac_{SSR}   =  %.3f #pm %.3f",fBkgSSDR_SBR, fBkgSSDRErr_SBR));
		top-=step;
	}
	latex->DrawLatex(left,top,Form("#tau_{SSR}   =  %.3f #pm %.3f",bkgTauSSD_SBR, bkgTauSSD_SBRErr));
	top-=step;
	top-=step;

	pad1->cd(0); pad1->SetLogy(0);
	ctauFrameBkgSBRPull->Draw();

	std::stringstream savectRSB;
	savectRSB << "Fit/ct_RSB_rap" << rapBin << "_pt" << ptBin << ".pdf";
	c1->SaveAs(savectRSB.str().c_str());

	delete c1;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendGreenDash;
	delete legendPink;
	return;

}

//==============================================
void plotLifeSig(RooWorkspace *ws, int rapBin, int ptBin, int nState){

	int nbins=230; //bin size 0.01 mm

	bool plotRegionLine=false;

	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	RooRealVar Jpsict(*ws->var("Jpsict"));
	RooRealVar *JpsictErr = ws->var("JpsictErr");

	RooPlot *ctauFrameSig=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(-.3, 2.));
	ctauFrameSig->SetName(Form("ctausig_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrameSig->SetYTitle("Events per 10 micron");
	ctauFrameSig->SetTitle("");
	RooPlot *ctauFrameSigPull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(-.3, 2.));
	ctauFrameSigPull->SetName(Form("pullctausig_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrameSigPull->SetTitle("");
	if(nState == 4) ctauFrameSigPull->SetXTitle("c_{#tau}^{J/#psi} [mm]");
	if(nState == 5) ctauFrameSigPull->SetXTitle("c_{#tau}^{#psi(2S)} [mm]");
	ctauFrameSigPull->GetYaxis()->SetTitle("pull");
	ctauFrameSigPull->GetXaxis()->SetTitleSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetTitleSize(0.08);
	ctauFrameSigPull->GetXaxis()->SetLabelSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetLabelSize(0.08);
	//ctauFrameSigPull->GetXaxis()->SetTitleOffset(1.0);
	ctauFrameSigPull->GetYaxis()->SetTitleOffset(0.4);
	ctauFrameSigPull->GetYaxis()->SetRangeUser(-5.5,5.5);


	RooDataSet *dataSR=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));

	RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("l_fitresult_rap%d_pt%d",rapBin,ptBin));
	if(!fitRlt) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; return; }
	fitRlt->Print();
	cout<<"dataSR->numEntries: "<<dataSR->numEntries()<<endl;
	//return 1;

	RooRealVar *promptMean_=(RooRealVar*)ws->var("promptMean");
	RooRealVar *ctResolution_=(RooRealVar*)ws->var("ctResolution");
	RooRealVar *ctResolution2_=(RooRealVar*)ws->var("ctResolution2");
	RooRealVar *fracGauss2_=(RooRealVar*)ws->var("fracGauss2");

	RooRealVar *nonPromptTau_=(RooRealVar*)ws->var("nonPromptTau");
	RooRealVar *fBkg_ = (RooRealVar*)ws->var("fBkg");
	RooRealVar *fPrompt_ = (RooRealVar*)ws->var("fPrompt");

	double promptMean = promptMean_->getVal();
	double promptMeanErr = promptMean_->getError();
	double promptCtRe = ctResolution_->getVal();
	double promptCtReErr = ctResolution_->getError();
	double promptCtRe2 = ctResolution2_->getVal();
	double promptCtRe2Err = ctResolution2_->getError();
	double fracGauss2 = fracGauss2_->getVal();
	double fracGauss2Err = fracGauss2_->getError();
	double nonPromptTau = nonPromptTau_->getVal();
	double nonPromptTauErr = nonPromptTau_->getError();

	double fBkg = fBkg_->getVal();
	double fBkgErr = fBkg_->getError();
	double fPrompt = fPrompt_->getVal();
	double fPromptErr = fPrompt_->getError();
	double fNonPrompt =  1.-fBkg-fPrompt;

	cout<<"promptMean: "<<promptMean<<endl;
	cout<<"promptCtRe: "<<promptCtRe<<endl;
	cout<<"nonPromptTau: "<<nonPromptTau<<endl;
	cout<<"fBkg: "<<fBkg<<endl;
	cout<<"fBkgErr: "<<fBkgErr<<endl;
	cout<<"fPrompt: "<<fPrompt<<endl;
	cout<<"fNonPrompt: "<<fNonPrompt<<endl;

	RooAddPdf *ModelLife = (RooAddPdf*)ws->pdf("fulllifetime");

	RooAddModel *Prompt = (RooAddModel*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *nonPromptSSD = (RooAbsPdf*)ws->pdf("nonPromptSSD");
	RooAbsPdf *backgroundlifetime = (RooAbsPdf*)ws->pdf("backgroundlifetime");


	int parsFit;
	//ploting signal region
	dataSR->plotOn(ctauFrameSig,MarkerSize(0.8));
	ModelLife->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			LineWidth(2), NumCPU(1));

	//------get chi2------------SHOULD DONE after PLOTTING------
	parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.

	int nBins_LSig=ctauFrameSig->GetNbinsX();
	double chi2Pre_LSig=ctauFrameSig->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_LSig=nBins_LSig-parsFit;  //num of degree of freedom
	double chi2_LSig=chi2Pre_LSig*ndof_LSig;

	RooHist* hpull_ctauSig = ctauFrameSig->pullHist() ;
	hpull_ctauSig->SetMarkerSize(0.8);
	for(int i=0;i<hpull_ctauSig->GetN();i++){
		hpull_ctauSig->SetPointEYlow(i,0.);
		hpull_ctauSig->SetPointEYhigh(i,0.);
	}
	ctauFrameSigPull->addPlotable(hpull_ctauSig,"P");

	//ModelLife->plotOn(ctauFrameSig,
	Prompt->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			//Components("promptLifetime"),
			Normalization(fPrompt),
			LineStyle(5),
			LineColor(kBlue),
			LineWidth(2), NumCPU(1));

	//ModelLife->plotOn(ctauFrameSig,
	nonPromptSSD->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			//Components("nonPromptSSD"),
			Normalization(fNonPrompt),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2), NumCPU(1));

	//ModelLife->plotOn(ctauFrameSig,
	backgroundlifetime->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			//Components("backgroundlifetime"),
			Normalization(fBkg),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2), NumCPU(1));

	double Ymax = ctauFrameSig->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;
	ctauFrameSig->SetMaximum(3*Ymax);
	ctauFrameSig->SetMinimum(0.001*Ymax);


	TH1* legendBlue = dataSR->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = dataSR->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = dataSR->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = dataSR->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = dataSR->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
	TH1* legendGreenDash = dataSR->createHistogram("legendGreenDash",JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
	TH1* legendPink = dataSR->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	//TLegend* LifetimeLegendBkgSBL=new TLegend(0.65,0.75,0.88,0.88);
	TLegend* LifetimeLegendBkgSBL=new TLegend(0.12,0.75,0.25,0.88);
	LifetimeLegendBkgSBL->SetFillColor(kWhite);
	LifetimeLegendBkgSBL->SetTextFont(42);
	LifetimeLegendBkgSBL->SetTextSize(0.025);
	LifetimeLegendBkgSBL->SetBorderSize(0.);
	LifetimeLegendBkgSBL->AddEntry(legendBlue,"sum sideband","l");
	LifetimeLegendBkgSBL->AddEntry(legendRed,"right side ","l");
	LifetimeLegendBkgSBL->AddEntry(legendGreenDash,"left side","l");
	LifetimeLegendBkgSBL->AddEntry(legendBlack,"double side","l");

	TLegend* LifetimeLegendBkgSBR=new TLegend(0.12,0.75,0.25,0.88);
	LifetimeLegendBkgSBR->SetFillColor(kWhite);
	LifetimeLegendBkgSBR->SetTextFont(42);
	LifetimeLegendBkgSBR->SetTextSize(0.025);
	LifetimeLegendBkgSBR->SetBorderSize(0.);
	LifetimeLegendBkgSBR->AddEntry(legendBlue,"sum sideband","l");
	LifetimeLegendBkgSBR->AddEntry(legendRed,"right side ","l");
	LifetimeLegendBkgSBR->AddEntry(legendGreenDash,"left side","l");
	LifetimeLegendBkgSBR->AddEntry(legendBlack,"double side","l");

	TLegend* LifetimeLegendSig=new TLegend(0.25,0.75,0.37,0.88);
	LifetimeLegendSig->SetFillColor(kWhite);
	LifetimeLegendSig->SetTextFont(42);
	LifetimeLegendSig->SetTextSize(0.025);
	LifetimeLegendSig->SetBorderSize(0.);
	LifetimeLegendSig->AddEntry(legendBlue,"sum","l");
	LifetimeLegendSig->AddEntry(legendBlueDash,"prompt","l");
	LifetimeLegendSig->AddEntry(legendRed,"non-prompt","l");
	LifetimeLegendSig->AddEntry(legendPink,"background","l");

	double left=0.4, top=0.85, textSize=0.025;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	TCanvas *c1=new TCanvas("c1","",800,900);
	c1->SetTickx();
	c1->SetTicky();

	c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
	pad1->SetGridy();
	pad1->SetLeftMargin(0.1);
	pad1->SetBottomMargin(0.2);
	pad1->Draw();
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
	pad2->SetLeftMargin(0.1);
	pad2->SetBottomMargin(0.1);
	pad2->Draw();

	double expectProbP = 0.95;
	double expectProbNP = 0.99;
	double ctBegin = 0.015;
	double ctCutP  = 0., ctCutNP = 0.;
	//if(plotRegionLine){
	//    ctBegin = 0.015;
	//    ctCutP  = evaluateCtau(expectProbP,Prompt,Jpsict,ctBegin);
	//    ctBegin = 0.018;
	//    ctCutNP = evaluateCtau(expectProbNP,Prompt,Jpsict,ctBegin);
	//}
	cout<<"ctCutP: "<<ctCutP<<endl;
	cout<<"ctCutNP: "<<ctCutNP<<endl;
	double minX=0., maxX=0.;
	double minY=0., maxY=0.;
	maxY = ctauFrameSig->GetMaximum()*0.7;
	TLine *lineP = new TLine(ctCutP, minY, ctCutP, maxY);
	lineP->SetLineWidth(1.0);
	lineP->SetLineColor(kBlue);
	lineP->SetLineStyle(7);
	TLine *lineNP = new TLine(ctCutNP, minY, ctCutNP, maxY);
	lineNP->SetLineWidth(1.0);
	lineNP->SetLineColor(kRed);
	lineNP->SetLineStyle(7);

	//Sig
	pad2->cd(0);pad2->SetLogy(1);
	ctauFrameSig->Draw();
	LifetimeLegendSig->Draw();
	if(plotRegionLine){ lineP->Draw("same");lineNP->Draw("same"); }

	left=0.4; top=0.85; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));

	top-=step;
	latex->SetTextColor(kRed);
	latex->DrawLatex(left,top,"Signal Mass Region");
	latex->SetTextColor(kBlack);
	top-=step;
	latex->DrawLatex(left,top,Form("mean  =  %.1f",promptMean));
	top-=step;
	latex->DrawLatex(left,top,Form("resolution1 =  %.1f",promptCtRe));
	top-=step;
	latex->DrawLatex(left,top,Form("resolution2  =  %.1f",promptCtRe2));

	left=0.6, top=0.85;
	latex->DrawLatex(left,top,Form("#chi^{2}/ndof   =  %.2f/%d",chi2_LSig,ndof_LSig));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Gauss2}   =  %.3f #pm %.3f",fracGauss2, fracGauss2Err));
	top-=step;
	latex->DrawLatex(left,top,Form("#tau_{NP}  =  %.3f #pm %.3f",nonPromptTau, nonPromptTauErr));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{P}   =  %.3f #pm %.3f",fPrompt, fPromptErr));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Bg}   =  %.3f #pm %.3f",fBkg, fBkgErr));
	top-=step;
	top-=step;
	latex->DrawLatex(left,top,Form("BFrac   =  %.3f ",fNonPrompt/(fNonPrompt+fPrompt)));

	pad1->cd(0); pad1->SetLogy(0);
	ctauFrameSigPull->Draw();

	std::stringstream savectlog;
	savectlog << "Fit/ct_SR_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
	c1->SaveAs(savectlog.str().c_str());

	delete c1;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendGreenDash;
	delete legendPink;
	return;
}

//==============================================
void plotLifeSig_linear(RooWorkspace *ws, int rapBin, int ptBin, int nState){

	int nbins=100; //bin size 0.005 mm
	TGaxis::SetMaxDigits(3);

	bool plotRegionLine=false;

	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	RooRealVar Jpsict(*ws->var("Jpsict"));
	RooRealVar *JpsictErr = ws->var("JpsictErr");

	RooPlot *ctauFrameSig=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(-.2, .3));
	ctauFrameSig->SetName(Form("ctausig_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrameSig->SetYTitle("Events per 5 micron");
	ctauFrameSig->SetTitle("");
	ctauFrameSig->GetYaxis()->SetTitleOffset(1.3);

	RooPlot *ctauFrameSigPull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(-.2, .3));
	ctauFrameSigPull->SetName(Form("pullctausig_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrameSigPull->SetTitle("");
	if(nState == 4) ctauFrameSigPull->SetXTitle("c_{#tau}^{J/#psi} [mm]");
	if(nState == 5) ctauFrameSigPull->SetXTitle("c_{#tau}^{#psi(2S)} [mm]");
	ctauFrameSigPull->GetYaxis()->SetTitle("pull");
	ctauFrameSigPull->GetXaxis()->SetTitleSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetTitleSize(0.08);
	ctauFrameSigPull->GetXaxis()->SetLabelSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetLabelSize(0.08);
	//ctauFrameSigPull->GetXaxis()->SetTitleOffset(1.0);
	ctauFrameSigPull->GetYaxis()->SetTitleOffset(0.4);
	ctauFrameSigPull->GetYaxis()->SetRangeUser(-5.5,5.5);


	RooDataSet *dataSR=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));

	RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("l_fitresult_rap%d_pt%d",rapBin,ptBin));
	if(!fitRlt) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; return; }
	fitRlt->Print();

	RooRealVar *promptMean_=(RooRealVar*)ws->var("promptMean");
	RooRealVar *ctResolution_=(RooRealVar*)ws->var("ctResolution");
	RooRealVar *ctResolution2_=(RooRealVar*)ws->var("ctResolution2");
	RooRealVar *fracGauss2_=(RooRealVar*)ws->var("fracGauss2");

	RooRealVar *nonPromptTau_=(RooRealVar*)ws->var("nonPromptTau");
	RooRealVar *fBkg_ = (RooRealVar*)ws->var("fBkg");
	RooRealVar *fPrompt_ = (RooRealVar*)ws->var("fPrompt");

	double promptMean = promptMean_->getVal();
	double promptMeanErr = promptMean_->getError();
	double promptCtRe = ctResolution_->getVal();
	double promptCtReErr = ctResolution_->getError();
	double promptCtRe2 = ctResolution2_->getVal();
	double promptCtRe2Err = ctResolution2_->getError();
	double fracGauss2 = fracGauss2_->getVal();
	double fracGauss2Err = fracGauss2_->getError();
	double nonPromptTau = nonPromptTau_->getVal();
	double nonPromptTauErr = nonPromptTau_->getError();
	double fBkg = fBkg_->getVal();
	double fBkgErr = fBkg_->getError();
	double fPrompt = fPrompt_->getVal();
	double fPromptErr = fPrompt_->getError();
	double fNonPrompt =  1.-fBkg-fPrompt;

	cout<<"promptMean: "<<promptMean<<endl;
	cout<<"promptCtRe: "<<promptCtRe<<endl;
	cout<<"nonPromptTau: "<<nonPromptTau<<endl;
	cout<<"fBkg: "<<fBkg<<endl;
	cout<<"fPrompt: "<<fPrompt<<endl;
	cout<<"fNonPrompt: "<<fNonPrompt<<endl;

	RooAddPdf *ModelLife = (RooAddPdf*)ws->pdf("fulllifetime");

	RooAddModel *Prompt = (RooAddModel*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *nonPromptSSD = (RooAbsPdf*)ws->pdf("nonPromptSSD");
	RooAbsPdf *backgroundlifetime = (RooAbsPdf*)ws->pdf("backgroundlifetime");

	int parsFit;
	//ploting signal region
	dataSR->plotOn(ctauFrameSig,MarkerSize(0.8));
	ModelLife->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			LineWidth(2), NumCPU(1));

	//------get chi2------------SHOULD DONE after PLOTTING------
	parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.

	int nBins_LSig=ctauFrameSig->GetNbinsX();
	double chi2Pre_LSig=ctauFrameSig->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_LSig=nBins_LSig-parsFit;  //num of degree of freedom
	double chi2_LSig=chi2Pre_LSig*ndof_LSig;

	RooHist* hpull_ctauSig = ctauFrameSig->pullHist() ;
	hpull_ctauSig->SetMarkerSize(0.8);
	for(int i=0;i<hpull_ctauSig->GetN();i++){
		hpull_ctauSig->SetPointEYlow(i,0.);
		hpull_ctauSig->SetPointEYhigh(i,0.);
	}
	ctauFrameSigPull->addPlotable(hpull_ctauSig,"P");

	//ModelLife->plotOn(ctauFrameSig,
	Prompt->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			//Components("promptLifetime"),
			Normalization(fPrompt),
			LineStyle(5),
			LineColor(kBlue),
			LineWidth(2), NumCPU(1));

	//ModelLife->plotOn(ctauFrameSig,
	nonPromptSSD->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			//Components("nonPromptSSD"),
			Normalization(fNonPrompt),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2), NumCPU(1));

	//ModelLife->plotOn(ctauFrameSig,
	backgroundlifetime->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			//Components("backgroundlifetime"),
			Normalization(fBkg),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2), NumCPU(1));
	double Ymax = ctauFrameSig->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;
	ctauFrameSig->SetMaximum(1.1*Ymax);
	//ctauFrameSig->SetMaximum(3*Ymax);
	//ctauFrameSig->SetMinimum(0.001*Ymax);


	TH1* legendBlue = dataSR->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = dataSR->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = dataSR->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = dataSR->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = dataSR->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
	TH1* legendGreenDash = dataSR->createHistogram("legendGreenDash",JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
	TH1* legendPink = dataSR->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	//TLegend* LifetimeLegendBkgSBL=new TLegend(0.65,0.75,0.88,0.88);
	TLegend* LifetimeLegendBkgSBL=new TLegend(0.12,0.75,0.25,0.88);
	LifetimeLegendBkgSBL->SetFillColor(kWhite);
	LifetimeLegendBkgSBL->SetTextFont(42);
	LifetimeLegendBkgSBL->SetTextSize(0.025);
	LifetimeLegendBkgSBL->SetBorderSize(0.);
	LifetimeLegendBkgSBL->AddEntry(legendBlue,"sum sideband","l");
	LifetimeLegendBkgSBL->AddEntry(legendRed,"right side ","l");
	LifetimeLegendBkgSBL->AddEntry(legendGreenDash,"left side","l");
	LifetimeLegendBkgSBL->AddEntry(legendBlack,"double side","l");

	TLegend* LifetimeLegendBkgSBR=new TLegend(0.12,0.75,0.25,0.88);
	LifetimeLegendBkgSBR->SetFillColor(kWhite);
	LifetimeLegendBkgSBR->SetTextFont(42);
	LifetimeLegendBkgSBR->SetTextSize(0.025);
	LifetimeLegendBkgSBR->SetBorderSize(0.);
	LifetimeLegendBkgSBR->AddEntry(legendBlue,"sum sideband","l");
	LifetimeLegendBkgSBR->AddEntry(legendRed,"right side ","l");
	LifetimeLegendBkgSBR->AddEntry(legendGreenDash,"left side","l");
	LifetimeLegendBkgSBR->AddEntry(legendBlack,"double side","l");

	////TLegend* LifetimeLegendSig=new TLegend(0.65,0.7,0.88,0.88);
	TLegend* LifetimeLegendSig=new TLegend(0.12,0.75,0.25,0.88);
	////TLegend* LifetimeLegendSig=new TLegend(0.48,0.75,0.6,0.88);
	//TLegend* LifetimeLegendSig=new TLegend(0.43,0.75,0.55,0.88);
	LifetimeLegendSig->SetFillColor(kWhite);
	LifetimeLegendSig->SetTextFont(42);
	LifetimeLegendSig->SetTextSize(0.025);
	LifetimeLegendSig->SetBorderSize(0.);
	LifetimeLegendSig->AddEntry(legendBlue,"sum","l");
	LifetimeLegendSig->AddEntry(legendBlueDash,"prompt","l");
	LifetimeLegendSig->AddEntry(legendRed,"non-prompt","l");
	LifetimeLegendSig->AddEntry(legendPink,"background","l");

	double left=0.65, top=0.85, textSize=0.032;
	//double left=0.6, top=0.85, textSize=0.025;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	TCanvas *c1=new TCanvas("c1","",800,900);
	//SetMyStyle();
	c1->SetTickx();
	c1->SetTicky();

	c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
	pad1->SetGridy();
	pad1->SetLeftMargin(0.1);
	pad1->SetBottomMargin(0.2);
	pad1->Draw();
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
	pad2->SetLeftMargin(0.1);
	pad2->SetBottomMargin(0.1);
	pad2->Draw();

	double expectProbP = 0.95;
	double expectProbNP = 0.99;
	double ctBegin = 0.015;
	double ctCutP  = 0., ctCutNP = 0.;
	//if(plotRegionLine){
	//    ctBegin = 0.015;
	//    ctCutP  = evaluateCtau(expectProbP,Prompt,Jpsict,ctBegin);
	//    ctBegin = 0.018;
	//    ctCutNP = evaluateCtau(expectProbNP,Prompt,Jpsict,ctBegin);
	//}
	cout<<"ctCutP: "<<ctCutP<<endl;
	cout<<"ctCutNP: "<<ctCutNP<<endl;
	double minX=0., maxX=0.;
	double minY=0., maxY=0.;
	maxY = ctauFrameSig->GetMaximum()*0.6;
	TLine *lineP = new TLine(ctCutP, minY, ctCutP, maxY);
	lineP->SetLineWidth(1.0);
	lineP->SetLineColor(kBlue);
	lineP->SetLineStyle(kSolid);
	maxY = ctauFrameSig->GetMaximum()*0.4;
	TLine *lineNP = new TLine(ctCutNP, minY, ctCutNP, maxY);
	lineNP->SetLineWidth(1.0);
	lineNP->SetLineColor(kRed);
	lineNP->SetLineStyle(kSolid);

	//Sig
	pad2->cd(0);pad2->SetLogy(0);
	ctauFrameSig->Draw();
	LifetimeLegendSig->Draw();
	if(plotRegionLine){ lineP->Draw("same");lineNP->Draw("same"); }

	left=0.65; top=0.85; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));

	top-=step;
	latex->SetTextColor(kRed);
	latex->DrawLatex(left,top,"Signal Mass Region");
	latex->SetTextColor(kBlack);
	top-=step;

	pad1->cd(0); pad1->SetLogy(0);
	ctauFrameSigPull->Draw();

	std::stringstream savect;
	savect << "Fit/ct_SR_rap" << rapBin << "_pt" << ptBin << ".pdf";
	c1->SaveAs(savect.str().c_str());

	delete c1;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendGreenDash;
	delete legendPink;
	return;
}
