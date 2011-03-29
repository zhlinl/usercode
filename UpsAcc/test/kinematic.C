#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include <iostream>

#define NMAX 100

void kinematic(const char* file = "MCupsilonTree.root", TString outFile="kinematic.root") {

  Int_t genUpsSize;
  Float_t genUpsPt[NMAX];
  Float_t genUpsEta[NMAX];
  Float_t genUpsPhi[NMAX];
  Float_t genUpsRapidity[NMAX];
  Int_t genMuSize;
  Float_t genMuPt[NMAX];
  Float_t genMuEta[NMAX];
  Float_t genMuPhi[NMAX];
  Int_t genMuCharge[NMAX];

  TFile* f1 = new TFile(file);
  TTree* t = (TTree*)f1->Get("UpsTree");
  t->SetBranchAddress("genUpsSize",&genUpsSize);
  t->SetBranchAddress("genUpsPt",genUpsPt);
  t->SetBranchAddress("genUpsEta",genUpsEta);
  t->SetBranchAddress("genUpsPhi",genUpsPhi);
  t->SetBranchAddress("genUpsRapidity",genUpsRapidity);
  t->SetBranchAddress("genMuSize",&genMuSize);
  t->SetBranchAddress("genMuPt",genMuPt);
  t->SetBranchAddress("genMuEta",genMuEta);
  t->SetBranchAddress("genMuPhi",genMuPhi);
  t->SetBranchAddress("genMuCharge",genMuCharge);

  cout << "entries=" << t->GetEntries() << endl;


  TH2F* genMuEtaPt = new TH2F("genMuEtaPt","",100,-2.5,2.5,100,0,5);
  TH2F* genHMuUpsPt = new TH2F("genHMuUpsPt","",300,0,30,300,0,30);
  TH2F* genLMuUpsPt = new TH2F("genLMuUpsPt","",300,0,30,300,0,30);


  genMuEtaPt->SetTitle(";#mu #eta;#mu pT (GeV/c);");
  genHMuUpsPt->SetTitle(";#mu pT (GeV/c);#Upsilon(1S) pT (GeV/c);");
  genLMuUpsPt->SetTitle(";#mu pT (GeV/c);#Upsilon(1S) pT (GeV/c);");

  for(int i=0; i<t->GetEntries(); i++){
    if(i%10000 == 0)
      std::cout<<i<<std::endl;
    t->GetEntry(i);
  
    double genupspt=genUpsPt[0];
    double genupsrap=fabs(genUpsRapidity[0]);
        double mupt1,mupt2,mueta1,mueta2;

    
    if(genMuSize < 2) continue;

     for(int tr1=0; tr1<genMuSize; tr1++){

      mupt1=genMuPt[tr1];
      mueta1=genMuEta[tr1];

        for(int tr2=tr1+1; tr2<genMuSize; tr2++){

          mupt2=genMuPt[tr2];
          mueta2=genMuEta[tr2];

          if ( genMuCharge[tr1]*genMuCharge[tr2] == -1  ) {
             genMuEtaPt->Fill(genMuEta[tr1],genMuPt[tr1]);
             genMuEtaPt->Fill(genMuEta[tr2],genMuPt[tr2]);
             if(genMuPt[tr1]>genMuPt[tr2]){
               genHMuUpsPt->Fill(genUpsPt[0],genMuPt[tr1]);
               genLMuUpsPt->Fill(genUpsPt[0],genMuPt[tr2]);
             }else{
               genHMuUpsPt->Fill(genUpsPt[0],genMuPt[tr2]);
               genLMuUpsPt->Fill(genUpsPt[0],genMuPt[tr1]);
              }
          }
        }
    }
    

  } // loop over tree entries

  

  TFile outfile(outFile,"recreate");

  genMuEtaPt->Write();
  genHMuUpsPt->Write();
  genLMuUpsPt->Write();

  outfile.Close();
}

