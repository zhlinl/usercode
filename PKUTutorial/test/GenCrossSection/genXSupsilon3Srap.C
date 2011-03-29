#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include <iostream>
#include <math.h>


#define NMAX 100

using namespace std;

void genXSupsilon3Srap(const char* file = "genUpsilon3S.root", const char* outputName = "genXSupsilon3Srap.root") {
  Float_t pt;
  Float_t y;

  TFile* f1 = new TFile(file);
  gDirectory->cd("tree");
  TTree* t = (TTree*)gDirectory->Get("probe_tree");
  t->SetBranchAddress("pt",&pt);
  t->SetBranchAddress("y",&y);

  TH1F* genY = new TH1F("genY","",10,0,2);
  genY->Sumw2();
  genY->SetTitle(";Upsilon |y|;d#sigma(#Upsilon(3S))xBr(#mu#mu)/dy (nb)");

  Float_t yCut = 2.0;
  Float_t confXS = 2.850;//nb
  Float_t generated = 10000000;
  Float_t filter = t->GetEntries()/generated;
  Float_t totXSxBr = confXS*filter;

  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    genY->Fill(fabs(y));
  }

  float total = genY->Integral();
  float totalError = sqrt(total);
  cout<<"tot:"<< total*totXSxBr/t->GetEntries()<<" +- "<<totalError*totXSxBr/t->GetEntries()<<endl;
  genY->Scale(totXSxBr/t->GetEntries());
  cout<<"total: "<<genY->Integral()<<endl;

  for(int i=1; i<=genY->GetNbinsX(); i++){
    genY->SetBinContent( i, genY->GetBinContent(i)/genY->GetBinWidth(i) );
    genY->SetBinError( i, genY->GetBinError(i)/genY->GetBinWidth(i) );
  }

  TFile out(outputName,"recreate");
  genY->Write();
  out.Close();
}

int main(int argc, const char** argv){
  genXSupsilon3Srap(argv[1],argv[2]);
}
