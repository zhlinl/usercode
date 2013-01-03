#include "calculatePar.cc"
#include "calcPol.C"

#include <string>
#include <iostream>
#include <sstream>

using namespace RooFit;
using namespace std;

vector<double> calculateInte(RooWorkspace *ws, RooDataSet *dataJpsictErr, double ctCutMin, double ctCutMax);

//======================
int order(int n){
	int total=1;
	for(int i=0;i<n;i++)
		total=total*2;
	return total;
}
int findEvenNum(double number){
	int thisNum=0;
	for(int n=0;n<100;n++){
		if(number >= order(n) && number <= order(n+1)){
			thisNum=order(n+1);
			break;
		}
	}
	return thisNum;
}

//======================
TH2D* ReSetBin(TH2D* hist, int nBinX, int nBinY, const std::stringstream& name, const std::stringstream& title){
	TH2D *tempHist = (TH2D*)hist->Clone("temp_BG_cosThetaPhiL");
	delete hist;
	hist = new TH2D(name.str().c_str(), title.str().c_str(),
			nBinX, onia::cosTMin, onia::cosTMax, nBinY, onia::phiPolMin, onia::phiPolMax);
	hist->Sumw2();
	TAxis *Xold = tempHist->GetXaxis();
	TAxis *Yold = tempHist->GetYaxis();
	TAxis *Xnew = hist->GetXaxis();
	TAxis *Ynew = hist->GetYaxis();
	for(int binX = 1; binX <= Xnew->GetNbins(); binX++){
		for(int binY = 1; binY <= Ynew->GetNbins(); binY++){
			double centerX = Xnew->GetBinCenter(binX);
			double centerY = Ynew->GetBinCenter(binY);

			//find the corresponding bin and bin error
			double binCont=0.,binErr=0.;
			bool findBin=false;
			for(int BinX = 1; BinX <= Xold->GetNbins(); BinX++){
				for(int BinY = 1; BinY <= Yold->GetNbins(); BinY++){
					double lowX = Xold->GetBinLowEdge(BinX);
					double upX  = Xold->GetBinUpEdge(BinX);
					double lowY = Yold->GetBinLowEdge(BinY);
					double upY  = Yold->GetBinUpEdge(BinY);
					if(centerX > lowX && centerX < upX && centerY > lowY && centerY < upY){
						binCont = tempHist->GetBinContent(BinX,BinY);
						binErr = tempHist->GetBinError(BinX,BinY);
						findBin=true;
					}
					if(findBin) break;
				}//BinY
			}//BinX
			//done
			hist->SetBinContent(binX,binY,binCont);
			hist->SetBinError(binX,binY,binErr);
		}//binY
	}//binX

	return hist;
}

//---------------------------------------------------------------------------------------------------------------
void bkgHistos(const std::string infilename, int rapBin, int ptBin, int nState, bool MC, bool f_BG_zero, bool doCtauUncer){

	const std::string
		datafilename = "tmpFiles/selEvents_data.root",
								 treename = "selectedData",
								 wsname = "ws_masslifetime";

	// input
	TFile *datafile = TFile::Open(datafilename.c_str());
	if(!datafile){
		std::cout << "Inputfile missing" << std::endl;
		return;
	}
	TTree *intree = (TTree *)datafile->Get(treename.c_str());
	TLorentzVector *lepP = 0, *lepN = 0, *jpsi = 0;
	TFile *fitfile = TFile::Open(infilename.c_str());
	if(!fitfile){
		std::cout << "fitfile is missing" << std::endl;
		return;
	}
	RooWorkspace *ws = (RooWorkspace*)fitfile->Get(wsname.c_str());

	gStyle->SetPadRightMargin(0.2);
	gROOT->SetStyle("Plain");
	gStyle->SetTitleBorderSize(0);

	// create output
	std::stringstream outfilename;
	outfilename << "tmpFiles/data_Psi" << nState-3 << "S_rap" << rapBin << "_pT" << ptBin << ".root";
	TFile *output = TFile::Open(outfilename.str().c_str(), "RECREATE");

	TTree *outtree = intree->CloneTree(0);

	// background histos
	TH2D *hBG_cosThetaPhiL[onia::kNbFrames];
	TH2D *hBG_cosThetaPhiR[onia::kNbFrames];
	TH2D *hBG_cosThetaPhi[onia::kNbFrames];
	TH2D *hNPBG_cosThetaPhi[onia::kNbFrames];
	TH2D *hNPS_cosThetaPhi[onia::kNbFrames];
	TH2D *hBGinNP_cosThetaPhi[onia::kNbFrames];
	TH2D *hBGinNP_cosThetaPhiL[onia::kNbFrames];
	TH2D *hBGinNP_cosThetaPhiR[onia::kNbFrames];
	TH2D *hTBG_cosThetaPhi[onia::kNbFrames];

	for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
		//book the 2D (cosTheta, phi) histos for the L and R mass sideband
		std::stringstream nameL, nameR, nameNP, nameBGinNPL, nameBGinNPR, title;
		nameL << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
		nameR << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
		nameNP << "hNPBG_cosThetaPhi_" << onia::frameLabel[iFrame];
		nameBGinNPL << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
		nameBGinNPR << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
		title << ";cos#theta_{"<< onia::frameLabel[iFrame] << "};#phi_{" << onia::frameLabel[iFrame] << "} [deg]";
		hBG_cosThetaPhiL[iFrame] = new TH2D(nameL.str().c_str(), title.str().c_str(),
				onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
		hBG_cosThetaPhiL[iFrame]->Sumw2();
		hBG_cosThetaPhiR[iFrame] = new TH2D(nameR.str().c_str(), title.str().c_str(),
				onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
		hBG_cosThetaPhiR[iFrame]->Sumw2();

		hNPBG_cosThetaPhi[iFrame] = new TH2D(nameNP.str().c_str(), title.str().c_str(),
				onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
		hNPBG_cosThetaPhi[iFrame]->Sumw2();

		hBGinNP_cosThetaPhiL[iFrame] = new TH2D(nameBGinNPL.str().c_str(), title.str().c_str(),
				onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
		hBGinNP_cosThetaPhiL[iFrame]->Sumw2();
		hBGinNP_cosThetaPhiR[iFrame] = new TH2D(nameBGinNPR.str().c_str(), title.str().c_str(),
				onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
		hBGinNP_cosThetaPhiR[iFrame]->Sumw2();

	} // iFrame

	// mean pT and y histos (background-subtracted)
	int nBins = 100;
	TH1D* pT_L   = new TH1D( "pTLSB", "pTLSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	TH1D* pT_R   = new TH1D( "pTRSB", "pTRSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	TH1D* pT_highct_L   = new TH1D( "pT_highct_LSB", "pT_highct_LSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	TH1D* pT_highct_R   = new TH1D( "pT_highcta_RSB", "pT_highct_RSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	TH1D* pT_NP   = new TH1D( "pTNP", "pTNP", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	TH1D* pT_PSR   = new TH1D( "pTPSR", "pTPSR", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	TH1D* rap_L   = new TH1D( "rapLSB", "rapLSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
	TH1D* rap_R   = new TH1D( "rapRSB", "rapRSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
	TH1D* rap_highct_L   = new TH1D( "rap_highct_LSB", "rap_highct_LSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
	TH1D* rap_highct_R   = new TH1D( "rap_highct_RSB", "rap_highct_RSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
	TH1D* rap_NP   = new TH1D( "rapNP", "rapNP", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
	TH1D* rap_PSR   = new TH1D( "rapPSR", "rapPSR", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);

	//------------------------------------------------------------------------------------------------
	// store pT and y borders
	TVectorD* pTBorder = new TVectorD(1, 2, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin], "END");
	TVectorD* yBorder = new TVectorD(1, 2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin], "END");
	pTBorder->Print();
	yBorder->Print();
	output->cd();
	pTBorder->Write();
	yBorder->Write();

	///////// CPU consuming plots? ///////
	//////// Dimuon mass plots? /////////
	/////// pT distribution plots? ////////
	/////// mass vs pT, vs rap, vs mean pT, vs mean rap ///////

	// set branches
	intree->SetBranchAddress("lepP", &lepP);
	intree->SetBranchAddress("lepN", &lepN);
	intree->SetBranchAddress("JpsiP", &jpsi);
	double jpsict = 0;
	intree->SetBranchAddress("Jpsict", &jpsict);

	//---------------------------------------------------------------------------------------------
	// INPUT FROM FIT
	// MASS FIT
	// reading fitting parameters
	std::stringstream mfitresult;
	mfitresult << "m_fitresult_rap" << rapBin << "_pt" << ptBin;
	RooFitResult *mresult = (RooFitResult*)ws->genobj(mfitresult.str().c_str());
	RooArgList mvarlist = mresult->floatParsFinal();

	// parameters
	RooRealVar *CBmass=(RooRealVar*)mvarlist.find("CBmass");
	RooRealVar *CBsigma1=(RooRealVar*)mvarlist.find("CBsigma");
	RooRealVar *CBsigma2=(RooRealVar*)mvarlist.find("CBsigma2");
	RooRealVar *fracCB1=(RooRealVar*)mvarlist.find("fracCB1");
	double mass = CBmass->getVal();
	double Sigma1 = CBsigma1->getVal();
	double Sigma2 = CBsigma2->getVal();
	double fCB1 = fracCB1->getVal();
	double sigma = sqrt( pow(Sigma1,2)*fCB1 + pow(Sigma2,2)*(1-fCB1) );

	// variables
	RooRealVar *m = ws->var("JpsiMass");
	RooRealVar *bkgLambda = ws->var("bkgLambda");
	RooRealVar *CBalpha = ws->var("CBalpha");
	RooRealVar *CBn = ws->var("CBn");
	RooRealVar* ct = ws->var("Jpsict");

	// pdf
	RooAbsPdf *bkgMass = (RooAbsPdf*)ws->pdf("bkgMassShape");
	RooAbsPdf *signalMass = (RooAbsPdf*)ws->pdf("massFull");
	RooAbsPdf *PRpdf = (RooAbsPdf*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *NPpdf = (RooAbsPdf*)ws->pdf("nonPromptSSD");
	RooAbsPdf *BGpdf = (RooAbsPdf*)ws->pdf("backgroundlifetime");

	// load snapshot with all results
	std::stringstream masssnapshotname;
	masssnapshotname << "m_snapshot_rap" << rapBin << "_pt" << ptBin;
	ws->loadSnapshot(masssnapshotname.str().c_str());

	TF1* funcBG = (TF1*)bkgMass->asTF(*m, *bkgLambda, *m);
	TF1* funcSig = (TF1*)signalMass->asTF(*m, RooArgList(*CBalpha, *fracCB1, *CBn, *CBsigma1, *CBsigma2, *CBmass), *m);

	//-------------------------------------------------------------------------------------------------
	// mass edges for left and right sideband
	double massMinL = onia::massMin;
	double massMaxR = onia::massMax;
	double massMaxL = mass - onia::nSigBkgLow  * sigma;
	double massMinR = mass + onia::nSigBkgHigh * sigma;
	double massMinSR = mass - onia::nSigMass * sigma;
	double massMaxSR = mass + onia::nSigMass * sigma;
	std::cout << "-------------------------------------------------------------\n" <<
		"left  sideband: mass window " << massMinL  << " < M < " << massMaxL  << " GeV\n" <<
		"right sideband: mass window " << massMinR  << " < M < " << massMaxR  << " GeV\n" <<
		"signal  region: mass window " << massMaxL  << " < M < " << massMinR  << " GeV\n" <<
		//"signal  region: mass window " << massMinSR << " < M < " << massMaxSR << " GeV\n" <<
		"-------------------------------------------------------------\n" << std::endl;

	// LIFETIME FIT
	// get dataset and number of events in dataset
	std::stringstream dataBin, data, cutSR;
	dataBin << "data_rap" << rapBin << "_pt" << ptBin << "_SR";
	data    << "data_rap" << rapBin << "_pt" << ptBin;
	//cutSR   << "JpsiMass > "<<massMinSR<<" && JpsiMass < "<<massMaxSR;
	cutSR   << "JpsiMass > "<<massMaxL<<" && JpsiMass < "<<massMinR;

	RooAbsData* Data    = ws->data(   data.str().c_str());
	RooDataSet *binData = (RooDataSet*)Data->reduce(
			Cut(cutSR.str().c_str()),
			Name(dataBin.str().c_str()),
			Title("Data For Fitting"));
	int entries = binData->numEntries();

	// load snapshot with all results
	std::stringstream snapshotname;
	snapshotname << "l_snapshot_rap" << rapBin << "_pt" << ptBin;
	ws->loadSnapshot(snapshotname.str().c_str());

	// get parameters for calculating fraction of left sideband
	RooRealVar *mean_ws=(RooRealVar*)ws->var("MeanSR");
	double mean = mean_ws->getVal();
	RooRealVar *mean_LSB_ws=(RooRealVar*)ws->var("MeanSBL");
	double mean_LSB = mean_LSB_ws->getVal();
	RooRealVar *mean_RSB_ws=(RooRealVar*)ws->var("MeanSBR");
	double mean_RSB = mean_RSB_ws->getVal();
	double fracLSB = 1 - (mean - mean_LSB)/(mean_RSB - mean_LSB);

	RooRealVar* fBkg_ws = (RooRealVar*)ws->var("FracBkg");
	double fBkg = fBkg_ws->getVal();

	// reading fitting parameters
	std::stringstream fitresult;
	fitresult << "l_fitresult_rap" << rapBin << "_pt" << ptBin;
	RooFitResult *result = (RooFitResult*)ws->genobj(fitresult.str().c_str());
	RooArgList varlist = result->floatParsFinal();

	// get parameters
	RooRealVar* fP_ws = (RooRealVar*)varlist.find("fPrompt");
	double fP1 = fP_ws->getVal();
	double fNP1 = 1 - fBkg - fP1;

	// calculate entries in signal region
	double nP = fP1*entries;
	double nNP = fNP1*entries;
	double nBG = fBkg*entries;

	std::cout << "-------------------------------------------------------------\n" <<
		"total number of events in signal region: " << entries << "\n" <<
		"prompt events in signal region: " << nP << "\n" <<
		"non prompt events in signal region: " << nNP << "\n" <<
		"background events in signal region: " << nBG << "\n" <<
		"fraction of left sideband: " << fracLSB << "\n" <<
		"-------------------------------------------------------------\n" << std::endl;

	int NumEvt = 0., maxEvt = 100000;
	if(Data->numEntries() < maxEvt)
		NumEvt = Data->numEntries();
	else
		NumEvt = maxEvt;

	RooRealVar Jpsict(*ws->var("Jpsict"));
	RooRealVar JpsictErr(*ws->var("JpsictErr"));
	RooDataSet *dataJpsictErr = (RooDataSet*)Data->reduce(SelectVars(RooArgSet(JpsictErr)),
			EventRange(0,NumEvt),Name("dataJpsictErr"));

	double meanPt = getMeanPt(rapBin,ptBin,infilename.c_str());
	// define sigma of prompt p.d.f., got from fit the trend
	//define function y = a + b * pT
	double a = 0.073, b = 0.0027;
	//proper decay length
	double L_decay = a + b * meanPt;
	//pseudo-proper decay length
	double l_pdecay =0., scale = 0.;
	if(nState==4) l_pdecay = L_decay * mass / meanPt ;
	if(nState==5) {
		if(rapBin == 1) scale = 1.21; //1.20
		if(rapBin == 2) scale = 1.29; //1.28
		if(rapBin ==3) scale = 1.31;  //1.43
		l_pdecay = scale * L_decay * 3.092 / meanPt ;
	}
	std::cout<<"l_pdecay: "<<l_pdecay<<std::endl;

	double nSigma = 0.;
	if(nState==4) nSigma = 2.5;
	if(nState==5) nSigma = 2.0;
	cout<<"nSigma: "<<nSigma<<endl;

	double ctauCut = nSigma*l_pdecay;
	double ctCutMinPR = -ctauCut, ctCutMaxPR = ctauCut;
	double ctCutMinNP =  ctauCut, ctCutMaxNP = 6.;

	// calculate fractions in prompt signal region
	vector<double> InteRltsPR = calculateInte(ws,dataJpsictErr,ctCutMinPR,ctCutMaxPR);
	double IntePRinPR = InteRltsPR[0];
	double InteNPinPR = InteRltsPR[1];
	double InteBGinPR = InteRltsPR[2];

	double fPinP  = IntePRinPR;
	double fNPinP = InteNPinPR;
	double fBGinP = InteBGinPR;

	std::cout << "-------------------------------------------------------------\n" <<
		"ctau cut at " << ctauCut << "\n" <<
		"fraction of prompt events within ctau cut: " << fPinP << " \n" <<
		"fraction of non prompt events within ctau cut: " << fNPinP << " \n" <<
		"fraction of background events within ctau cut: " << fBGinP << " \n" <<
		"-------------------------------------------------------------\n" << std::endl;

	//------------------------------------------------------------------------------------------------
	// calculate fractions in non prompt signal region
	vector<double> InteRltsNP = calculateInte(ws,dataJpsictErr,ctCutMinNP,ctCutMaxNP);
	double IntePRinNP = InteRltsNP[0];
	double InteNPinNP = InteRltsNP[1];
	double InteBGinNP = InteRltsNP[2];

	double fPbkg = IntePRinNP;
	double fNPbkg = InteNPinNP;
	double fBGbkg = InteBGinNP;

	double fBGinNP = fBGbkg*fBkg/(fBGbkg*fBkg + fNPbkg*fNP1 + fPbkg*fP1);  // comb. background fraction in high ctau signal region
	double fP      =   fPinP*fP1/(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);  // prompt fraction in prompt signal region
	double fNPB    = fNPinP*fNP1/(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);  // non prompt background events in prompt region
	double fBGsig  = fBGinP*fBkg/(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);  // background in prompt signal region
	double fTBGsig =  fNPB + fBGsig;

	//------------------------------------------------------------------------------------------------
	// evaluate error on non-prompt fraction
	double fPerr=0., fNPerr=0., fBGerr=0.;

	if(doCtauUncer){
		int nEvents = 200;
		if(nState==4) nEvents = 100;
		double promptFrac = 0, nonpromptFrac = 0, bkgFrac = 0;

		RooRealVar* bkgTauDSD = (RooRealVar*)ws->var("bkgTauDSD");
		RooRealVar* bkgTauFD = (RooRealVar*)ws->var("bkgTauFD");
		RooRealVar* bkgTauSSD_SBL = (RooRealVar*)ws->var("bkgTauSSD_SBL");
		RooRealVar* bkgTauSSD_SBR = (RooRealVar*)ws->var("bkgTauSSD_SBR");
		RooRealVar* fBkgDSD_SBL = (RooRealVar*)ws->var("fBkgDSD_SBL");
		RooRealVar* fBkgDSD_SBR = (RooRealVar*)ws->var("fBkgDSD_SBR");
		RooRealVar* fBkgSSDR_SBL = (RooRealVar*)ws->var("fBkgSSDR_SBL");
		RooRealVar* fBkgSSDR_SBR = (RooRealVar*)ws->var("fBkgSSDR_SBR");
		RooRealVar* fPrompt = (RooRealVar*)ws->var("fPrompt");
		RooRealVar* fBKG = (RooRealVar*)ws->var("fBkg");
		RooRealVar* fracGauss2 = (RooRealVar*)ws->var("fracGauss2");
		RooRealVar* nonPromptTau = (RooRealVar*)ws->var("nonPromptTau");
		RooArgSet *paraVars = new RooArgSet(*bkgTauDSD,*bkgTauFD,*bkgTauSSD_SBL,*bkgTauSSD_SBR,
				*fBkgDSD_SBL,*fBkgDSD_SBR, *fBkgSSDR_SBL,*fBkgSSDR_SBR,
				*nonPromptTau);
		paraVars->add(RooArgSet(*fPrompt,*fBKG,*fracGauss2));

		// create Hesse pdf and generate dataset
		RooAbsPdf *multiVarPdf = (RooAbsPdf*)result->createHessePdf(*paraVars);
		RooDataSet *multiVarData = (RooDataSet*)multiVarPdf->generate(*paraVars,nEvents);

		TH1D* histPFracDist = new TH1D();
		TH1D* histNPFracDist = new TH1D();
		TH1D* histBGFracDist = new TH1D();

		double IntePR=0., InteNP=0., InteBG=0.;
		double fracP1=0., fracNP1=0., fracBkg=0.;
		for(int n = 0; n < nEvents; n++) {
			if(n%40==0) std::cout << (double)n/nEvents*100 << "%" << std::endl;
			RooArgSet* args = (RooArgSet*)multiVarData->get(n);

			double bkgTauFD_ = ((RooRealVar*)args->find("bkgTauFD"))->getVal();
			ws->var("bkgTauFD")->setVal(bkgTauFD_);
			double bkgTauSSD_SBL_ = ((RooRealVar*)args->find("bkgTauSSD_SBL"))->getVal();
			ws->var("bkgTauSSD_SBL")->setVal(bkgTauSSD_SBL_);
			double bkgTauSSD_SBR_ = ((RooRealVar*)args->find("bkgTauSSD_SBR"))->getVal();
			ws->var("bkgTauSSD_SBR")->setVal(bkgTauSSD_SBR_);
			double fBkgDSD_SBL_ = ((RooRealVar*)args->find("fBkgDSD_SBL"))->getVal();
			ws->var("fBkgDSD_SBL")->setVal(fBkgDSD_SBL_);
			double fBkgDSD_SBR_ = ((RooRealVar*)args->find("fBkgDSD_SBR"))->getVal();
			ws->var("fBkgDSD_SBR")->setVal(fBkgDSD_SBR_);
			if(!(nState==4 && ptBin > 7)){
				double fBkgSSDR_SBL_ = ((RooRealVar*)args->find("fBkgSSDR_SBL"))->getVal();
				ws->var("fBkgSSDR_SBL")->setVal(fBkgSSDR_SBL_);
				double fBkgSSDR_SBR_ = ((RooRealVar*)args->find("fBkgSSDR_SBR"))->getVal();
				ws->var("fBkgSSDR_SBR")->setVal(fBkgSSDR_SBR_);
			}
			double fPrompt_ = ((RooRealVar*)args->find("fPrompt"))->getVal();
			ws->var("fPrompt")->setVal(fPrompt_);
			double fBkg_ = ((RooRealVar*)args->find("fBkg"))->getVal();
			ws->var("fBkg")->setVal(fBkg_);
			double nonPromptTau_ = ((RooRealVar*)args->find("nonPromptTau"))->getVal();
			ws->var("nonPromptTau")->setVal(nonPromptTau_);
			double bkgTauDSD_ = ((RooRealVar*)args->find("bkgTauDSD"))->getVal();
			ws->var("bkgTauDSD")->setVal(bkgTauDSD_);
			double fracGauss2_ = ((RooRealVar*)args->find("fracGauss2"))->getVal();
			ws->var("fracGauss2")->setVal(fracGauss2_);

			fracP1 = fPrompt_; fracBkg = fBkg_; fracNP1 =  1.-fracP1-fracBkg;
			vector<double> InteRltsTemp = calculateInte(ws,dataJpsictErr,ctCutMinPR,ctCutMaxPR);
			IntePR    = InteRltsTemp[0];
			InteNP    = InteRltsTemp[1];
			InteBG    = InteRltsTemp[2];

			promptFrac    = IntePR * fracP1;
			nonpromptFrac = InteNP * fracNP1;
			bkgFrac       = InteBG * fracBkg;

			histPFracDist->Fill(promptFrac);
			histNPFracDist->Fill(nonpromptFrac);
			histBGFracDist->Fill(bkgFrac);
		}

		fPerr  = histPFracDist->GetRMS() /(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);
		fNPerr = histNPFracDist->GetRMS()/(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);
		fBGerr = histBGFracDist->GetRMS()/(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);
	}
	double fTBGerr = TMath::Sqrt(TMath::Power(fNPerr,2) + TMath::Power(fBGerr,2));

	output->cd();
	//------------------------------------------------------------------------------------------------
	// fill histogram with combinatorial background fraction
	std::stringstream bkgname;
	bkgname << ";;fraction of comb. BG in " << onia::nSigMass << "  sigma window, prompt region";
	TH1D* hFracBG = new TH1D("comb_background_fraction", bkgname.str().c_str(), 1, 0., 1.);
	if(MC == false && f_BG_zero == false){
		hFracBG->SetBinContent(1, fBGsig);
		hFracBG->SetBinError(1, fBGerr);
	}
	// for MC set background fraction to 0.001
	if(MC == true || f_BG_zero == true)  hFracBG->SetBinContent(1, 0.001);
	hFracBG->Write();

	// fill histogram with non prompt background fraction
	std::stringstream NPbkgname;
	NPbkgname << ";;fraction of non prompt BG in " << onia::nSigMass << "  sigma window, prompt region";
	TH1D* hFracNPBG = new TH1D("nonprompt_background_fraction", NPbkgname.str().c_str(), 1, 0., 1.);
	if(MC == false && f_BG_zero == false){
		hFracNPBG->SetBinContent(1, fNPB);
		hFracNPBG->SetBinError(1, fNPerr);
	}
	// for MC set background fraction to 0.001
	if(MC == true || f_BG_zero == true)  hFracNPBG->SetBinContent(1, 0.001);
	hFracNPBG->Write();

	// fill histogram with total background fraction
	std::stringstream tbkgname;
	tbkgname << ";;fraction of total BG in " << onia::nSigMass << "  sigma window, prompt region";
	TH1D* hFracTBG = new TH1D("background_fraction", tbkgname.str().c_str(), 1, 0., 1.);
	if(MC == false && f_BG_zero == false){
		hFracTBG->SetBinContent(1, fTBGsig);
		hFracTBG->SetBinError(1, fTBGerr);
	}
	if(MC == true || f_BG_zero == true)  hFracTBG->SetBinContent(1, 0.001);
	hFracTBG->Write();

	// fill histogram with prompt fraction
	std::stringstream Pname;
	Pname << ";;fraction of prompt events in " << onia::nSigMass << "  sigma window, prompt region (1-fNP-fBkg)";
	TH1D* hFracP = new TH1D("prompt_fraction", Pname.str().c_str(), 1, 0., 1.);
	hFracP->SetBinContent(1, fP);
	hFracP->SetBinError(1, fPerr);
	hFracP->Write();

	// fill histogram with fraction of LSB
	TH1D* hFracLSB = new TH1D("fraction_LSB", ";;f_{LSB}", 1, 0., 1.);
	hFracLSB->SetBinContent(1, fracLSB);
	hFracLSB->Write();

	// fill histogram with events in signal region
	std::stringstream evtSRname;
	evtSRname << ";;events in " << onia::nSigMass << "  sigma window";
	TH1D* hEvtSR = new TH1D("events_SR", evtSRname.str().c_str(), 3, 0., 3.);
	hEvtSR->SetBinContent(1, nP);
	hEvtSR->SetBinContent(2, nNP);
	hEvtSR->SetBinContent(3, nBG);
	hEvtSR->Write();

	// fill histogram with events in prompt signal region
	std::stringstream evtPSRname;
	evtPSRname << ";;events in " << onia::nSigMass << "  sigma window and low ctau region";
	TH1D* hEvtPSR = new TH1D("events_promptSR", evtPSRname.str().c_str(), 3, 0., 3.);
	hEvtPSR->SetBinContent(1, fPinP*nP);
	hEvtPSR->SetBinContent(2, fNPinP*nNP);
	hEvtPSR->SetBinContent(3, fBGinP*nBG);
	hEvtPSR->Write();

	// fill histogram with weighted sigma
	TH1D* hsigma = new TH1D("weighted_sigma", ";;weighted #sigma", 1, 0, 1);
	hsigma->SetBinContent(1, sigma);
	hsigma->Write();

	//------------------------------------------------------------------------------------------------
	// build the 3D (pT, |y|, M) histos for the L and R mass sideband
	TH3D* hBG_pTRapMass_L = new TH3D("hBG_pTRapMass_L", ";p_{T} [GeV/c]; |y|; M [GeV]",
			7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
			2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
			7, massMinL, massMaxR); // signal mass window!
	hBG_pTRapMass_L->Sumw2();

	TH3D* hBG_pTRapMass_R = new TH3D("hBG_pTRapMass_R", ";p_{T} [GeV/c]; |y|; M [GeV]",
			7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
			2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
			7, massMinL, massMaxR); // signal mass window!
	hBG_pTRapMass_R->Sumw2();

	TH3D* hBG_pTRapMass_highct_L = new TH3D("hBG_pTRapMass_highct_L", ";p_{T} [GeV/c]; |y|; M [GeV]",
			7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
			2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
			7, massMinL, massMaxR);
	hBG_pTRapMass_highct_L->Sumw2();

	TH3D* hBG_pTRapMass_highct_R = new TH3D("hBG_pTRapMass_highct_R", ";p_{T} [GeV/c]; |y|; M [GeV]",
			7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
			2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
			7, massMinL, massMaxR);
	hBG_pTRapMass_highct_R->Sumw2();

	TH3D* hNP_pTRapMass = new TH3D("hNP_pTRapMass_NP", ";p_{T} [GeV/c]; |y|; M [GeV]",
			7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
			2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
			7, massMinL, massMaxR);
	hNP_pTRapMass->Sumw2();

	//------------------------------------------------------------------------------------------------
	// loop through tree, fill background histos and save data in pt and y bins
	int index = -1;
	int n = intree->GetEntries();

	//// for 1S, rap1,pt 1 2 3 4; rap2,pt 1 2 3 4, not use the full statistics
	//if(nState==4 && ((rapBin==1 && ptBin<5)||(rapBin==2 && ptBin<5))) 
	//	n = n / (6-ptBin);

	for(int i = 0; i < n; i++){

		Long64_t iEntry = intree->LoadTree(i);
		intree->GetEntry(iEntry);
		if(i % 100000 == 0) {std::cout << "entry " << i << " out of " << n << std::endl;}

		if(jpsi->Pt() >= onia::pTRange[rapBin-1][ptBin-1] && 
				jpsi->Pt() < onia::pTRange[rapBin-1][ptBin] && 
				TMath::Abs(jpsi->Rapidity()) >= onia::rapForPTRange[rapBin-1] && 
				TMath::Abs(jpsi->Rapidity()) < onia::rapForPTRange[rapBin]){

			//store TLorenzVectors of the two muons in the given pT and rap cell
			//store only events from signal region
			if(jpsi->M() >= massMaxL && jpsi->M() < massMinR && TMath::Abs(jpsict) < ctauCut){
				outtree->Fill();
				pT_PSR->Fill(jpsi->Pt());
				rap_PSR->Fill(TMath::Abs(jpsi->Rapidity()));
			}

			// store cosTheta and phi distributions of the background
			// events gets index 0 if it is in the left sideband and 1 if it is in the right one
			// events with index 2 are from the high ctau non prompt region
			// events with index 3 and 4 are from the high ctau region from left and right sideband
			if(jpsi->M() >= massMinL && jpsi->M() < massMaxL && TMath::Abs(jpsict) < ctauCut){
				index = 0;
				hBG_pTRapMass_L->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcBG->GetRandom(massMaxL, massMinR));
				pT_L->Fill(jpsi->Pt());
				rap_L->Fill(TMath::Abs(jpsi->Rapidity()));
			}
			else if(jpsi->M() >= massMinR && jpsi->M() < massMaxR && TMath::Abs(jpsict) < ctauCut){
				index = 1;
				hBG_pTRapMass_R->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcBG->GetRandom(massMaxL, massMinR));
				pT_R->Fill(jpsi->Pt());
				rap_R->Fill(TMath::Abs(jpsi->Rapidity()));
			}
			else if(jpsi->M() >= massMaxL && jpsi->M() < massMinR && jpsict >= ctCutMinNP){
				index = 2;
				hNP_pTRapMass->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcSig->GetRandom(massMaxL, massMinR));
				pT_NP->Fill(jpsi->Pt());
				rap_NP->Fill(TMath::Abs(jpsi->Rapidity()));
			}
			else if(jpsi->M() >= massMinL && jpsi->M() < massMaxL && jpsict >= ctCutMinNP){
				index = 3;
				hBG_pTRapMass_highct_L->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcSig->GetRandom(massMaxL, massMinR));
				pT_highct_L->Fill(jpsi->Pt());
				rap_highct_L->Fill(TMath::Abs(jpsi->Rapidity()));
			}
			else if(jpsi->M() >= massMinR && jpsi->M() < massMaxR && jpsict >= ctCutMinNP){
				index = 4;
				hBG_pTRapMass_highct_R->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcSig->GetRandom(massMaxL, massMinR));
				pT_highct_R->Fill(jpsi->Pt());
				rap_highct_R->Fill(TMath::Abs(jpsi->Rapidity()));
			}
			else continue;

			calcPol(*lepP, *lepN);

			for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
				if(MC == false){
					if(index == 0) hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 1) hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 2) hNPBG_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 3) hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 4) hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
				}
				else{
					hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
				}
			} // iFrame

			if(MC || f_BG_zero){
				hBG_pTRapMass_L->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), gRandom->Uniform(massMinL, massMaxR));
				hBG_pTRapMass_R->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), gRandom->Uniform(massMinL, massMaxR));
				hBG_pTRapMass_highct_L->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), gRandom->Gaus(massMinL, massMaxR));
				hBG_pTRapMass_highct_R->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), gRandom->Gaus(massMinL, massMaxR));
			}

		} // if(onia...)
	} // i

	//---------------- binning algorithm
	int nBinsPhi = 16, nBinsCosth = 160;
	int totalBins=0, filledBins=0;
	for(int binCosth = 0; binCosth<hBGinNP_cosThetaPhiR[2]->GetNbinsX(); binCosth++){
		for(int binPhi = 0; binPhi<hBGinNP_cosThetaPhiR[2]->GetNbinsY(); binPhi++){
			totalBins++;
			int binContent = hNPBG_cosThetaPhi[2]->GetBinContent(binCosth+1,binPhi+1);//Val: Use here the NP histo (same physical coverage, but better filled, no holes -> better estimate of coverage)
			if(binContent>0) filledBins++;
		}
	}
	double coverage = 2*(double)filledBins/(double)totalBins;
	nBinsCosth = 16*2/coverage;
	cout<<"coverage: "<<coverage<<endl;
	cout<<"nBinsPhi: "<<nBinsPhi<<endl;
	cout<<"nBinsCosth: "<<nBinsCosth<<endl;

	int IntBG = hBGinNP_cosThetaPhiR[2]->Integral(); 
	if(IntBG > hBGinNP_cosThetaPhiL[2]->Integral())
		IntBG = hBGinNP_cosThetaPhiL[2]->Integral();
	if(IntBG > hBG_cosThetaPhiL[2]->Integral())
		IntBG = hBG_cosThetaPhiL[2]->Integral();
	if(IntBG > hBG_cosThetaPhiR[2]->Integral())
		IntBG = hBG_cosThetaPhiR[2]->Integral();
	//Val: calculate here the integral of the lowstatBG histo (calculate all integrals of the 4 BG regions, and use the one with the smallest integral). At the moment a few 1S RSB histos have N < 10. This will be solved by this change.

	int IntNPBG = hNPBG_cosThetaPhi[2]->Integral();

	int nBinsPhiBG = nBinsPhi, nBinsCosthBG = nBinsCosth,
			nBinsPhiNPBG = nBinsPhi, nBinsCosthNPBG = nBinsCosth;

	bool binningDone=false;
	//average events per-bin cell
	double Naverge=0;
	Naverge = (double)IntBG/((double)nBinsPhi*nBinsCosth*coverage/2.);
	if(Naverge>10){
		binningDone=true;
	}
	cout<<"binningDone: "<<binningDone<<endl;

	if(Naverge<10){
		cout<<"-----------------------------------"<<endl;
		cout<<"old nBinsCosth: "<<nBinsCosth<<endl;
		cout<<"old nBinsPhi: "<<nBinsPhi<<endl;
		nBinsCosth = findEvenNum((double)nBinsCosth);
		nBinsPhi = findEvenNum(nBinsCosth*coverage/2.);
		cout<<"new nBinsCosth: "<<nBinsCosth<<endl;
		cout<<"new nBinsPhi: "<<nBinsPhi<<endl;
		cout<<"-----------------------------------"<<endl;
		//Val: set here nBinsPhi to the lowest 2^n, such that nBinsPhi > nBinsCosth*coverage/2 (nBinsCosth*coverage/2 is what I call 'effective' bins in costh). Making nPhiBins here bigger makes ensures that there is maximally a factor of 2 in between nPhiBins and nCosthBins, as nBinsPhi is cut in half first.

		int nBinsPhiMin=4;
		int nBinsCosthMin=4;

		//BG
		nBinsCosthBG=nBinsCosth;
		nBinsPhiBG=nBinsPhi; //Val: added this line
		double NavergeBG=0.;
		for(int i=0; i<500; i++){
			if(i%50==0) cout<<"i: "<<i<<endl;
			if(nBinsPhiBG/2<nBinsPhiMin && nBinsCosthBG/2<nBinsCosthMin) break; // If the mimimum number of bins for both phi and costh are reached, stop the loop

			//Change the binning, first in phi, then in costh:
			if(nBinsPhiBG/2>=nBinsPhiMin) nBinsPhiBG=nBinsPhiBG/2;  //This ensures a mimimum number of bins in phi, e.g. 4
			NavergeBG = (double)IntBG/((double)nBinsPhiBG*nBinsCosthBG*coverage/2.);
			if(NavergeBG>10) break;

			if(nBinsCosthBG/2>=nBinsCosthMin) nBinsCosthBG=nBinsCosthBG/2; //This ensures a mimimum number of bins in costh, e.g. 4
			NavergeBG = (double)IntBG/((double)nBinsPhiBG*nBinsCosthBG*coverage/2.);
			if(NavergeBG>10) break;
		}
		cout<<"NavergeBG: "<<NavergeBG<<endl;

		//NPBG
		nBinsCosthNPBG=nBinsCosth;
		nBinsPhiNPBG=nBinsPhi; //Val: added this line
		double NavergeNPBG=0.;
		for(int i=0; i<500; i++){
			if(i%50==0) cout<<"i: "<<i<<endl;
			if(nBinsPhiNPBG/2<nBinsPhiMin && nBinsCosthNPBG/2<nBinsCosthMin) break; // If the mimimum number of bins for both phi and costh are reached, stop the loop

			//Change the binning, first in phi, then in costh:
			if(nBinsPhiNPBG/2>=nBinsPhiMin) nBinsPhiNPBG=nBinsPhiNPBG/2;  //This ensures a mimimum number of bins in phi, e.g. 4
			NavergeNPBG = (double)IntNPBG/((double)nBinsPhiNPBG*nBinsCosthNPBG*coverage/2.);
			if(NavergeNPBG>10) break; //Val

			if(nBinsCosthNPBG/2>=nBinsCosthMin) nBinsCosthNPBG=nBinsCosthNPBG/2; //This ensures a mimimum number of bins in costh, e.g. 4
			NavergeNPBG = (double)IntNPBG/((double)nBinsPhiNPBG*nBinsCosthNPBG*coverage/2.);
			if(NavergeNPBG>10) break;
		}
		cout<<"NavergeNPBG: "<<NavergeNPBG<<endl;

	}

	cout<<"final binning for NPBG....."<<endl;
	cout<<"nBinsPhiNPBG: "<<nBinsPhiNPBG<<endl;
	cout<<"nBinsCosthNPBG: "<<nBinsCosthNPBG<<endl;

	cout<<"final binning for BG....."<<endl;
	cout<<"nBinsPhiBG: "<<nBinsPhiBG<<endl;
	cout<<"nBinsCosthBG: "<<nBinsCosthBG<<endl;

	//loop again with new binning
	for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
		//book the 2D (cosTheta, phi) histos for the L and R mass sideband
		std::stringstream nameL, nameR, nameNP, nameBGinNPL, nameBGinNPR, title;
		nameL << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
		nameR << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
		nameNP << "hNPBG_cosThetaPhi_" << onia::frameLabel[iFrame];
		nameBGinNPL << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
		nameBGinNPR << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
		title << ";cos#theta_{"<< onia::frameLabel[iFrame] << "};#phi_{" << onia::frameLabel[iFrame] << "} [deg]";

		delete hBG_cosThetaPhiL[iFrame];
		hBG_cosThetaPhiL[iFrame] = new TH2D(nameL.str().c_str(), title.str().c_str(),
				nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
		hBG_cosThetaPhiL[iFrame]->Sumw2();
		delete hBG_cosThetaPhiR[iFrame];
		hBG_cosThetaPhiR[iFrame] = new TH2D(nameR.str().c_str(), title.str().c_str(),
				nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
		hBG_cosThetaPhiR[iFrame]->Sumw2();

		delete hNPBG_cosThetaPhi[iFrame];
		hNPBG_cosThetaPhi[iFrame] = new TH2D(nameNP.str().c_str(), title.str().c_str(),
				nBinsCosthNPBG, onia::cosTMin, onia::cosTMax, nBinsPhiNPBG, onia::phiPolMin, onia::phiPolMax);
		hNPBG_cosThetaPhi[iFrame]->Sumw2();

		delete hBGinNP_cosThetaPhiL[iFrame];
		hBGinNP_cosThetaPhiL[iFrame] = new TH2D(nameBGinNPL.str().c_str(), title.str().c_str(),
				nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
		hBGinNP_cosThetaPhiL[iFrame]->Sumw2();
		delete hBGinNP_cosThetaPhiR[iFrame];
		hBGinNP_cosThetaPhiR[iFrame] = new TH2D(nameBGinNPR.str().c_str(), title.str().c_str(),
				nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
		hBGinNP_cosThetaPhiR[iFrame]->Sumw2();

	} // iFrame

	for(int i = 0; i < n; i++){

		Long64_t iEntry = intree->LoadTree(i);
		intree->GetEntry(iEntry);
		if(i % 100000 == 0) {std::cout << "entry " << i << " out of " << n << std::endl;}

		if(jpsi->Pt() >= onia::pTRange[rapBin-1][ptBin-1] 
				&& jpsi->Pt() < onia::pTRange[rapBin-1][ptBin] &&
				TMath::Abs(jpsi->Rapidity()) >= onia::rapForPTRange[rapBin-1] 
				&& TMath::Abs(jpsi->Rapidity()) < onia::rapForPTRange[rapBin]){
			if(jpsi->M() >= massMinL && jpsi->M() < massMaxL && TMath::Abs(jpsict) < ctauCut)
				index = 0;
			else if(jpsi->M() >= massMinR && jpsi->M() < massMaxR && TMath::Abs(jpsict) < ctauCut)
				index = 1;
			else if(jpsi->M() >= massMaxL && jpsi->M() < massMinR && jpsict >= ctCutMinNP)
				index = 2;
			else if(jpsi->M() >= massMinL && jpsi->M() < massMaxL && jpsict >= ctCutMinNP)
				index = 3;
			else if(jpsi->M() >= massMinR && jpsi->M() < massMaxR && jpsict >= ctCutMinNP)
				index = 4;
			else continue;

			calcPol(*lepP, *lepN);
			for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
				if(MC == false){
					if(index == 0) hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 1) hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 2) hNPBG_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 3) hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 4) hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
				}
				else{
					hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
				}
			} // iFrame
		} // if(onia...)
	} // i
	//------loop finished

	//---------------- end- binning algorithm

	//----------------------------------------------------------------------------------------------------
	// write background histos to file
	// 3D (pT, |y|, M) histos
	hBG_pTRapMass_L->Write();
	hBG_pTRapMass_R->Write();
	hBG_pTRapMass_L->Scale(fracLSB/hBG_pTRapMass_L->Integral());
	hBG_pTRapMass_R->Scale((1.-fracLSB)/hBG_pTRapMass_R->Integral());
	std::string namepTrapMasslowct = "comb_background_pTrapMass";
	TH3D* hBG_pTRapMass_lowct = (TH3D*) hBG_pTRapMass_L->Clone(namepTrapMasslowct.c_str());
	hBG_pTRapMass_lowct->Add(hBG_pTRapMass_R);
	hBG_pTRapMass_lowct->Write();

	hBG_pTRapMass_highct_L->Write();
	hBG_pTRapMass_highct_R->Write();
	hBG_pTRapMass_highct_L->Scale(fracLSB/hBG_pTRapMass_highct_L->Integral());
	hBG_pTRapMass_highct_R->Scale((1.-fracLSB)/hBG_pTRapMass_highct_R->Integral());
	std::string namepTrapMasshighct = "comb_background_highct_pTrapMass";
	TH3D* hBG_pTRapMass_highct = (TH3D*) hBG_pTRapMass_highct_L->Clone(namepTrapMasshighct.c_str());
	hBG_pTRapMass_highct->Add(hBG_pTRapMass_highct_R);
	hBG_pTRapMass_highct->Write();

	hNP_pTRapMass->Scale(1./hNP_pTRapMass->Integral());
	hBG_pTRapMass_highct->Scale(fBGinNP/hBG_pTRapMass_highct->Integral());
	std::string namepTrapMassNPS = "NPS_highct_pTrapMass";
	TH3D* hNPS_pTRapMass = (TH3D*) hNP_pTRapMass->Clone(namepTrapMassNPS.c_str());
	hNPS_pTRapMass->Add(hBG_pTRapMass_highct, -1);
	hNPS_pTRapMass->Write();

	hBG_pTRapMass_lowct->Scale(fBGsig/hBG_pTRapMass_lowct->Integral());
	hNPS_pTRapMass->Scale(fNPB/hNPS_pTRapMass->Integral());
	std::string namepTrapMass = "background_pTrapMass";
	TH3D* hBG_pTRapMass = (TH3D*) hNPS_pTRapMass->Clone(namepTrapMass.c_str());
	hBG_pTRapMass->Add(hBG_pTRapMass_lowct);
	hBG_pTRapMass->Write();

	// mean pT histos
	pT_L->Scale(fracLSB/pT_L->Integral());
	pT_R->Scale((1.-fracLSB)/pT_R->Integral());
	pT_L->Add(pT_R);

	pT_highct_L->Scale(fracLSB/pT_highct_L->Integral());
	pT_highct_R->Scale((1.-fracLSB)/pT_highct_R->Integral());
	pT_highct_L->Add(pT_highct_R);

	pT_NP->Scale(1./pT_NP->Integral());
	pT_highct_L->Scale(fBGinNP/pT_highct_L->Integral());
	pT_NP->Add(pT_highct_L, -1);

	pT_L->Scale(fBGsig/pT_L->Integral());
	pT_NP->Scale(fNPB/pT_NP->Integral());
	pT_L->Add(pT_NP);

	pT_PSR->Scale(1./pT_PSR->Integral());
	pT_PSR->Add(pT_L, -1);
	double meanPT = pT_PSR->GetMean();

	std::stringstream meanPTname;
	meanPTname << ";;mean p_{T}";
	TH1D* h_meanPT = new TH1D("mean_pT", meanPTname.str().c_str(), 1, 0., 1.);
	h_meanPT->SetBinContent(1, meanPT);
	h_meanPT->Write();

	// mean y histos
	rap_L->Scale(fracLSB/rap_L->Integral());
	rap_R->Scale((1.-fracLSB)/rap_R->Integral());
	rap_L->Add(rap_R);

	rap_highct_L->Scale(fracLSB/rap_highct_L->Integral());
	rap_highct_R->Scale((1.-fracLSB)/rap_highct_R->Integral());
	rap_highct_L->Add(rap_highct_R);

	rap_NP->Scale(1./rap_NP->Integral());
	rap_highct_L->Scale(fBGinNP/rap_highct_L->Integral());
	rap_NP->Add(rap_highct_L, -1);

	rap_L->Scale(fBGsig/rap_L->Integral());
	rap_NP->Scale(fNPB/rap_NP->Integral());
	rap_L->Add(rap_NP);

	rap_PSR->Scale(1./rap_PSR->Integral());
	rap_PSR->Add(rap_L, -1);
	double meanY = rap_PSR->GetMean();

	std::stringstream meanYname;
	meanYname << ";;mean |y|";
	TH1D* h_meanY = new TH1D("mean_y", meanPTname.str().c_str(), 1, 0., 1.);
	h_meanY->SetBinContent(1, meanY);
	h_meanY->Write();

	//---rebin BG to be same NPBG
	bool ResetBin=true;
	if(ResetBin){
		cout<<"Resetting bins...."<<endl;
		for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
			std::stringstream nameL, nameR, nameNP, nameBGinNPL, nameBGinNPR, title;
			nameL << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
			nameR << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
			nameNP << "hNPBG_cosThetaPhi_" << onia::frameLabel[iFrame];
			nameBGinNPL << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
			nameBGinNPR << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
			title << ";cos#theta_{"<< onia::frameLabel[iFrame] << "};#phi_{" << onia::frameLabel[iFrame] << "} [deg]";

			hBG_cosThetaPhiL[iFrame] = ReSetBin(hBG_cosThetaPhiL[iFrame], nBinsCosthNPBG, nBinsPhiNPBG, nameL, title);
			hBG_cosThetaPhiR[iFrame] = ReSetBin(hBG_cosThetaPhiR[iFrame], nBinsCosthNPBG, nBinsPhiNPBG, nameR, title);
			hBGinNP_cosThetaPhiL[iFrame] = ReSetBin(hBGinNP_cosThetaPhiL[iFrame], nBinsCosthNPBG, nBinsPhiNPBG, nameBGinNPL, title);
			hBGinNP_cosThetaPhiR[iFrame] = ReSetBin(hBGinNP_cosThetaPhiR[iFrame], nBinsCosthNPBG, nBinsPhiNPBG, nameBGinNPR, title);
		}//iFrame
	}
	//---rebin finished

	//======================================================
	for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
		hBG_cosThetaPhiL[iFrame]->Write();
		hBG_cosThetaPhiR[iFrame]->Write();
		hBGinNP_cosThetaPhiL[iFrame]->Write();
		hBGinNP_cosThetaPhiR[iFrame]->Write();
		hNPBG_cosThetaPhi[iFrame]->Write();

		// combinatorial background in signal region (prompt region)
		// combination of left and right sideband
		hBG_cosThetaPhiL[iFrame]->Scale(fracLSB/hBG_cosThetaPhiL[iFrame]->Integral());
		hBG_cosThetaPhiR[iFrame]->Scale((1.-fracLSB)/hBG_cosThetaPhiR[iFrame]->Integral());
		std::stringstream name;
		name << "comb_background_costhphi" << onia::frameLabel[iFrame];
		hBG_cosThetaPhi[iFrame] = (TH2D *) hBG_cosThetaPhiL[iFrame]->Clone(name.str().c_str());
		hBG_cosThetaPhi[iFrame]->Add(hBG_cosThetaPhiR[iFrame]);
		hBG_cosThetaPhi[iFrame]->Write();

		// combinatorial background in high ctau region
		// combination of left and right sideband in high ctau region
		hBGinNP_cosThetaPhiL[iFrame]->Scale(fracLSB/hBGinNP_cosThetaPhiL[iFrame]->Integral());
		hBGinNP_cosThetaPhiR[iFrame]->Scale((1.-fracLSB)/hBGinNP_cosThetaPhiR[iFrame]->Integral());
		std::stringstream nameBGinNP;
		nameBGinNP << "background_NPR_costhphi" << onia::frameLabel[iFrame];
		hBGinNP_cosThetaPhi[iFrame] = (TH2D *) hBGinNP_cosThetaPhiL[iFrame]->Clone(nameBGinNP.str().c_str());
		hBGinNP_cosThetaPhi[iFrame]->Add(hBGinNP_cosThetaPhiR[iFrame]);
		hBGinNP_cosThetaPhi[iFrame]->Write();

		// non prompt background in high ctau region
		// (hNPBG_cosThetaPhi + hBGinNP_cosThetaPhi)_norm - fBGinNP * (hBGinNP_cosThetaPhi)_norm
		hNPBG_cosThetaPhi[iFrame]->Scale(1./hNPBG_cosThetaPhi[iFrame]->Integral());
		hBGinNP_cosThetaPhi[iFrame]->Scale(fBGinNP/hBGinNP_cosThetaPhi[iFrame]->Integral());
		std::stringstream nameNPS;
		nameNPS << "background_NPSR_costhphi" << onia::frameLabel[iFrame];
		hNPS_cosThetaPhi[iFrame] = (TH2D *) hNPBG_cosThetaPhi[iFrame]->Clone(nameNPS.str().c_str());
		hNPS_cosThetaPhi[iFrame]->Add(hBGinNP_cosThetaPhi[iFrame], -1);
		hNPS_cosThetaPhi[iFrame]->Write();

		// total background
		// fNPBG * (hNPS_cosThetaPhi)_norm + fBGsig * (hBG_cosThetaPhi)_norm
		hBG_cosThetaPhi[iFrame]->Scale(fBGsig/hBG_cosThetaPhi[iFrame]->Integral());
		hNPS_cosThetaPhi[iFrame]->Scale(fNPB/hNPS_cosThetaPhi[iFrame]->Integral());
		std::stringstream nameTBG;
		nameTBG << "background_costhphi" << onia::frameLabel[iFrame];
		hTBG_cosThetaPhi[iFrame] = (TH2D *) hBG_cosThetaPhi[iFrame]->Clone(nameTBG.str().c_str());
		hTBG_cosThetaPhi[iFrame]->Add(hNPS_cosThetaPhi[iFrame]);
		hTBG_cosThetaPhi[iFrame]->Write();


	}// iFrame

	outtree->Write();
	datafile->Close();
	fitfile->Close();

} // void


vector<double> calculateInte(RooWorkspace *ws, RooDataSet *dataJpsictErr, double ctCutMin, double ctCutMax){

	double ctMin = -2., ctMax = 6.;
	int fineBins = 8000;

	RooAbsPdf *PRpdf = (RooAbsPdf*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *NPpdf = (RooAbsPdf*)ws->pdf("nonPromptSSD");
	RooAbsPdf *BGpdf = (RooAbsPdf*)ws->pdf("backgroundlifetime");

	RooRealVar Jpsict(*ws->var("Jpsict"));
	RooRealVar JpsictErr(*ws->var("JpsictErr"));
	Jpsict.setMin(ctMin);   Jpsict.setMax(ctMax);

	RooDataSet *genDataPR = PRpdf->generate(Jpsict,ProtoData(*dataJpsictErr));
	RooDataSet *genDataNP = NPpdf->generate(Jpsict,ProtoData(*dataJpsictErr));
	RooDataSet *genDataBG = BGpdf->generate(Jpsict,ProtoData(*dataJpsictErr));

	TH2F* histPR2D = (TH2F*)genDataPR->createHistogram("histPR2D",Jpsict,Binning(fineBins),YVar(JpsictErr,Binning(fineBins/10)));
	TH1F* histPR   = (TH1F*)histPR2D->ProjectionX();
	TH2F* histNP2D = (TH2F*)genDataNP->createHistogram("histNP2D",Jpsict,Binning(fineBins),YVar(JpsictErr,Binning(fineBins/10)));
	TH1F* histNP   = (TH1F*)histNP2D->ProjectionX();
	TH2F* histBG2D = (TH2F*)genDataBG->createHistogram("histBG2D",Jpsict,Binning(fineBins),YVar(JpsictErr,Binning(fineBins/10)));
	TH1F* histBG   = (TH1F*)histBG2D->ProjectionX();

	histPR->SetLineColor(kRed);
	histPR->SetMarkerColor(kRed);
	histNP->SetLineColor(kBlue);
	histNP->SetMarkerColor(kBlue);
	histBG->SetLineColor(kBlack);
	histBG->SetMarkerColor(kBlack);

	histPR->Scale(1./histPR->Integral());
	histNP->Scale(1./histNP->Integral());
	histBG->Scale(1./histBG->Integral());

	int BinLow = 0, BinHigh = 0;
	bool getLow = false, getHigh = false;
	for(int bin = 0; bin < fineBins; bin++){
		if(ctCutMin > histPR->GetBinLowEdge(bin) && ctCutMin < histPR->GetBinLowEdge(bin+1)){
			BinLow = bin; getLow = true;}
		if(ctCutMax > histPR->GetBinLowEdge(bin) && ctCutMax < histPR->GetBinLowEdge(bin+1)){
			BinHigh = bin; getHigh = true;}
		if(getLow && getHigh) break;
	}

	double IntePR    = histPR->Integral(BinLow,BinHigh);
	double InteNP    = histNP->Integral(BinLow,BinHigh);
	double InteBG    = histBG->Integral(BinLow,BinHigh);

	vector<double> InteRlts;
	InteRlts.push_back(IntePR);
	InteRlts.push_back(InteNP);
	InteRlts.push_back(InteBG);

	delete histPR;
	delete histNP;
	delete histBG;
	delete histPR2D;
	delete histNP2D;
	delete histBG2D;

	return InteRlts;
}


