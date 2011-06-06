#include "TTree.h"

Char_t *fileNameData = "/scratch/knuenz/Polarization/RootInput/TTree_pol_Mu0Track0Jpsi_MCprompt.root";

void CutTree(){

	int NMin=300000;

/*	char outputfilename[200];
	sprintf(outputfilename,"Results/ProblemBin.txt");
	printf("output filename is: %s\n", outputfilename);
	FILE *outputFile = fopen(outputfilename,"w");
*/
	TFile* fInData = new TFile(fileNameData);

	TTree *treeData = (TTree*)fInData->Get("data");

	treeData->Print();


	TTree *fChain;
	fChain=treeData;

	Long64_t nentries = fChain->GetEntries();
	cout<<nentries<<endl;
	Long64_t nb = 0;

	Double_t        JpsiMass;
	Double_t        JpsiPt;
	Double_t        JpsiRap;
	Double_t        Jpsict;
	Int_t        JpsiType_idx;
	Int_t        MCType_idx;
	Double_t        costh_CS;
	Double_t        phi_CS;
	Double_t        costh_HX;
	Double_t        phi_HX;

	TBranch        *b_JpsiMass;   //!
	TBranch        *b_JpsiPt;   //!
	TBranch        *b_JpsiRap;   //!
	TBranch        *b_Jpsict;   //!
	TBranch        *b_JpsiType_idx;   //!
	TBranch        *b_MCType_idx;   //!
	TBranch        *b_costh_CS;   //!
	TBranch        *b_phi_CS;   //!
	TBranch        *b_costh_HX;   //!
	TBranch        *b_phi_HX;   //!

	fChain->SetBranchAddress("JpsiMass", &JpsiMass, &b_JpsiMass);
	fChain->SetBranchAddress("JpsiPt", &JpsiPt, &b_JpsiPt);
	fChain->SetBranchAddress("JpsiRap", &JpsiRap, &b_JpsiRap);
	fChain->SetBranchAddress("Jpsict", &Jpsict, &b_Jpsict);
	fChain->SetBranchAddress("JpsiType_idx", &JpsiType_idx, &b_JpsiType_idx);
	fChain->SetBranchAddress("MCType_idx", &MCType_idx, &b_MCType_idx);
	fChain->SetBranchAddress("costh_CS", &costh_CS, &b_costh_CS);
	fChain->SetBranchAddress("phi_CS", &phi_CS, &b_phi_CS);
	fChain->SetBranchAddress("costh_HX", &costh_HX, &b_costh_HX);
	fChain->SetBranchAddress("phi_HX", &phi_HX, &b_phi_HX);


//	Long64_t jentry=NMin;

    TNtuple* pseudodata = new TNtuple("data","data","JpsiMass:JpsiPt:JpsiRap:Jpsict:costh_CS:phi_CS:costh_HX:phi_HX:JpsiType_idx:MCType_idx");
    TNtuple* data = new TNtuple("data","data","JpsiMass:JpsiPt:JpsiRap:Jpsict:costh_CS:phi_CS:costh_HX:phi_HX:JpsiType_idx:MCType_idx");

    for (Long64_t jentry=0; jentry<NMin;jentry++) {
       		    	 if(jentry % 1000 == 0) printf("event %d\n", (Int_t) jentry);

       	   nb = fChain->GetEntry(jentry);   //nbytes += nb;

       	pseudodata->Fill(JpsiMass,JpsiPt,JpsiRap,Jpsict,costh_CS,phi_CS,costh_HX,phi_HX,JpsiType_idx,MCType_idx);

        }

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
   		    	 if(jentry % 1000 == 0) printf("event %d\n", (Int_t) jentry);

   	   nb = fChain->GetEntry(jentry);   //nbytes += nb;

// 	  fprintf(outputFile, "%f\n",Jpsict);


   	data->Fill(JpsiMass,JpsiPt,JpsiRap,Jpsict,costh_CS,phi_CS,costh_HX,phi_HX,JpsiType_idx,MCType_idx);

    }

    pseudodata->Print();

    pseudodata->SaveAs("/scratch/knuenz/Polarization/RootInput/TNtuple_red_PR_pseudo.root");

	data->Print();

	data->SaveAs("/scratch/knuenz/Polarization/RootInput/TNtuple_red_PR.root");

//	RooDataSetData.Print();


//	fclose(outputFile);


//	treeData.Print();

	return;
}
