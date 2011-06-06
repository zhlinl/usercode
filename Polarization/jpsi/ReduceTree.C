#include "TTree.h"

Char_t *fileNameData = "/scratch/knuenz/Polarization/RootInput/TNtuple_red_PR.root";
Char_t *fileNamepseudoData = "/scratch/knuenz/Polarization/RootInput/TNtuple_red_PR_pseudo.root";

void ReduceTree(){


	RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
	RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
	RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
	RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
	RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
	RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
	RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
	RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
	RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2.5);//0=PR,1=NP,2=BK
	RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,2.5);//0=GG,1=GT,2=TT


	RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
	varlist.add(MCType_idx);
	varlist.add(JpsiType_idx);



	TFile* fInData = new TFile(fileNameData);
	TFile* fInpseudoData = new TFile(fileNamepseudoData);

	TTree *treeData = (TTree*)fInData->Get("data");
	TTree *treepseudoData = (TTree*)fInpseudoData->Get("data");


	treeData->Print();
	treepseudoData->Print();

	RooDataSet* dataset = new RooDataSet("data","data",treeData,varlist);//RooArgSet(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX,MCweight,MCType_idx,JpsiType_idx));
	RooDataSet* pseudodataset = new RooDataSet("data","data",treepseudoData,varlist);//RooArgSet(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX,MCweight,MCType_idx,JpsiType_idx));

//	RooDataSet* reddataset = (RooDataSet*)dataset->reduce("JpsiPt < 20 && JpsiPt > 12.5 && abs(JpsiRap) < 0.9 && abs(JpsiRap) > 0");


	dataset->Print();
	pseudodataset->Print();

	TTree* tree = dataset->tree();
	TTree* pseudotree = pseudodataset->tree();

	double nentries = tree->GetEntries();
	double pseudonentries = pseudotree->GetEntries();

	double MCfraction = nentries/(nentries+pseudonentries);

	cout<<"MC entries: "<<nentries<<endl;
	cout<<"MC pseudo entries: "<<pseudonentries<<endl;
	cout<<"MC fraction: "<<MCfraction<<endl;


	tree->SaveAs("/scratch/knuenz/Polarization/RootInput/TTree_red_PR.root");
	pseudotree->SaveAs("/scratch/knuenz/Polarization/RootInput/TTree_red_PR.root");



	return;
}
