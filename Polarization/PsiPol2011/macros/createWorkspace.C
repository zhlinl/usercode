#include "rootIncludes.inc"
#include "commonVar.h"

#include <string>
#include <iostream>
#include <sstream>

using namespace RooFit;

void createWorkspace(const std::string &infilename, int nState){

	gROOT->SetStyle("Plain");
	gStyle->SetTitleBorderSize(0);

	// Set some strings
	const std::string
		workspacename = "ws_masslifetime",
		treename = "selectedData";

	// Get the tree from the data file
	TFile *f = TFile::Open(infilename.c_str());
	TTree *tree = (TTree*)f->Get(treename.c_str());

	// Set branch addresses in tree to be able to import tree to roofit
	TLorentzVector* jpsi = new TLorentzVector;
	tree->SetBranchAddress("JpsiP",&jpsi);
	double massErr = 0;
	tree->SetBranchAddress("JpsiMassErr",&massErr);
	double Vprob = 0;
	tree->SetBranchAddress("JpsiVprob",&Vprob);
	double lifetime = 0;
	tree->SetBranchAddress("Jpsict",&lifetime);
	double lifetimeErr = 0;
	tree->SetBranchAddress("JpsictErr",&lifetimeErr);

	// define variables necessary for J/Psi(Psi(2S)) mass,lifetime fit
	RooRealVar* JpsiMass =
		new RooRealVar("JpsiMass", "M [GeV]", onia::massMin, onia::massMax);
	RooRealVar* JpsiMassErr =
		new RooRealVar("JpsiMassErr", "#delta M [GeV]", 0, 5);
	RooRealVar* JpsiRap =
		new RooRealVar("JpsiRap", "y", -onia::rap, onia::rap);
	RooRealVar* JpsiPt =
		new RooRealVar("JpsiPt", "p_{T} [GeV]", 0. ,100.);
	RooRealVar* Jpsict =
		new RooRealVar("Jpsict", "c_{#tau}^{J/#psi} [mm]", -1., 2.5);
	RooRealVar* JpsictErr =
		new RooRealVar("JpsictErr", "Error on c_{#tau}^{J/#psi} [mm]", 0.0001, 1);
	RooRealVar* JpsiVprob =
		new RooRealVar("JpsiVprob", "", 0.01, 1.);

	if(nState==5){
		Jpsict->SetTitle("c_{#tau}^{#psi(2S)} [mm]");
		JpsictErr->SetTitle("Error on c_{#tau}^{#psi(2S)} [mm]");
	}
	// Set bins
	Jpsict->setBins(10000,"cache");
	Jpsict->setBins(100);
	JpsiMass->setBins(100);
	JpsictErr->setBins(100);

	// The list of data variables    
	RooArgList dataVars(*JpsiMass,*JpsiMassErr,*JpsiRap,*JpsiPt,*Jpsict,*JpsictErr,*JpsiVprob);

	// construct dataset to contain events
	RooDataSet* fullData = new RooDataSet("fullData","The Full Data From the Input ROOT Trees",dataVars);

	int entries = tree->GetEntries();

	// loop through events in tree and save them to dataset
	for (int ientries = 0; ientries < entries; ientries++) {
		if (ientries%100000==0) std::cout << "event " << ientries << " of " << entries <<  std::endl;

		tree->GetEntry(ientries);

		double M =jpsi->M();
		double y=jpsi->Rapidity();
		double pt=jpsi->Pt();

		if (M > JpsiMass->getMin() && M < JpsiMass->getMax()
				&& massErr > JpsiMassErr->getMin() && massErr < JpsiMassErr->getMax()
				&& pt > JpsiPt->getMin() && pt < JpsiPt->getMax()
				&& y > JpsiRap->getMin() && y < JpsiRap->getMax()
				&& lifetime > Jpsict->getMin() && lifetime < Jpsict->getMax()
				&& lifetimeErr > JpsictErr->getMin() && lifetimeErr < JpsictErr->getMax()
				&& Vprob > JpsiVprob->getMin() && Vprob < JpsiVprob->getMax()
				){

			JpsiPt->setVal(pt); 
			JpsiRap->setVal(y); 
			JpsiMass->setVal(M);
			JpsiMassErr->setVal(massErr);
			Jpsict->setVal(lifetime);
			JpsictErr->setVal(lifetimeErr);
			JpsiVprob->setVal(Vprob);

			fullData->add(dataVars);
		}
	}//ientries


	//------------------------------------------------------------------------------------------------------------------
	// Define workspace and import datasets

	// Get datasets binned in pT an y
	for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){ 

		Double_t yMin  = onia::rapForPTRange[iRap];
		Double_t yMax  = onia::rapForPTRange[iRap+1];

		for(int iPT = 0; iPT < onia::kNbPTBins[iRap]; iPT++){

			Double_t ptMin = onia::pTRange[iRap][iPT];
			Double_t ptMax = onia::pTRange[iRap][iPT+1]; 

			// output file name and workspace
			std::stringstream outfilename;
			outfilename << "tmpFiles/backupWorkSpace/fit_Psi" << nState-3 << "S_rap" << iRap+1 << "_pt" << iPT+1 << ".root";
			RooWorkspace* ws = new RooWorkspace(workspacename.c_str());

			// define pt and y cuts on dataset
			std::stringstream cutString;
			cutString << "(JpsiPt >= " << ptMin << " && JpsiPt < "<< ptMax << ") && "
				<< "(TMath::Abs(JpsiRap) >= " << yMin << " && TMath::Abs(JpsiRap) < " << yMax << ")";

			// get the dataset for the fit
			RooDataSet* binData = (RooDataSet*)fullData->reduce(cutString.str().c_str());
			std::stringstream name;
			name << "data_rap" << iRap+1 << "_pt" << iPT+1;
			binData->SetNameTitle(name.str().c_str(), "Data For Fitting");    

			// Import variables to workspace
			ws->import(*binData);
			ws->writeToFile(outfilename.str().c_str());
		}//iPT
	}//iRap

	f->Close();
}
