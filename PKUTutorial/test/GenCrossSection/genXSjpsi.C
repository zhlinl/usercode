#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include <iostream>

#define NMAX 100

void genXSupsilon3S(const char* file = "genUpsilon3S.root", const char* outputName = "genXSupsilon3S.root") {
  Float_t pt;
  Float_t y;

  TFile* f1 = new TFile(file);
  gDirectory->cd("tree");
  TTree* t = (TTree*)gDirectory->Get("probe_tree");
  std::cout<<t<<std::endl;
  t->SetBranchAddress("pt",&pt);
  t->SetBranchAddress("y",&y);

  const int nPtBins = 6;
  double ptBinEdges[nPtBins+1] = {0, 2, 3, 4, 6, 9, 20};
  TH1F* genPt = new TH1F("genPt","",nPtBins,ptBinEdges);
  genPt->Sumw2();
  genPt->SetTitle(";Upsilon pT (GeV/c);d#sigma(#Upsilon(3S))xBr(#mu#mu)/dpT (nb/GeV), |y|<2.0;");

  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    if(fabs(y) < 2.0){
      genPt->Fill(pt);
    }
  }

  Float_t branchingFraction = 0.0248;
  Float_t totalXS = 4.e3;
  Float_t xsbr = 28.5;
  Int_t totalGeneratedEvents = 20000;
  Float_t integratedLuminosity = totalGeneratedEvents / xsbr; //totalXS / branchingFraction;
  genPt->Scale(1./integratedLuminosity);

  std::cout<<genPt->Integral()<<std::endl;

  for(int i=1; i<=genPt->GetNbinsX(); i++){
    genPt->SetBinContent( i, genPt->GetBinContent(i)/genPt->GetBinWidth(i) );
    genPt->SetBinError( i, genPt->GetBinError(i)/genPt->GetBinWidth(i) );
  }

  TFile out(outputName,"recreate");
  genPt->Write();
  out.Close();
}

int main(int argc, const char** argv){
  genXSupsilon3S(argv[1],argv[2]);
}
