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

double invert(TH1D* h, int ind, double y){
  double y1, y2;
  y1 = h->GetBinContent(ind);
  y2 = h->GetBinContent(ind+1);
  if((y1-y) * (y2-y) < 0){
    double x1, x2;
    x1 = h->GetBinCenter(ind);
    x2 = h->GetBinCenter(ind+1);
    return (y*(x2-x1)-y1*x2+y2*x1)/(y2-y1);
  }else{
    int new = (y1-y2)*(y1-y)>0 ? ind+1 : ind-1;
    if(new==0 || new==h->GetNbinsX())
      return 0;
    return invert(h, new, y);
  }
}

void genXSupsilon2S(const char* file = "genUpsilon2S.root", const char* outputName = "genXSupsilon2S.root") {
  Float_t pt;
  Float_t y;

  TFile* f1 = new TFile(file);
  gDirectory->cd("tree");
  TTree* t = (TTree*)gDirectory->Get("probe_tree");
  t->SetBranchAddress("pt",&pt);
  t->SetBranchAddress("y",&y);

  const int nPtBins = 8;
  double ptBinEdges[nPtBins+1] = {0, 2, 4, 6, 9, 12, 16, 20, 30};
  TH1D* genPtBinned = new TH1D("genPtBinned","",nPtBins,ptBinEdges);
  TH1D* genPt = new TH1D("genPt","",60,0.,30.);
  genPtBinned->Sumw2();
  genPt->Sumw2();
  genPt->SetTitle(";Upsilon pT (GeV/c);d#sigma(#Upsilon(2S))xBr(#mu#mu)/dpT (nb/GeV), |y|<2.0;");

  Float_t yCut = 2.0;
  Float_t confXS = 75.400;//nb
  Float_t generated = 100000000;
  Float_t filter = t->GetEntries()/generated;
  Float_t totXSxBr = confXS*filter;

  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    if(fabs(y) < yCut){
      genPtBinned->Fill(pt);
      genPt->Fill(pt);
    }
  }

genPtBinned->Print("range");
genPt->Print("range");

  float total = genPt->Integral();
  float totalError = sqrt(total);
  cout<<"tot:"<< total*totXSxBr/t->GetEntries()<<" +- "<<totalError*totXSxBr/t->GetEntries()<<endl;
  genPtBinned->Scale(totXSxBr/t->GetEntries());
  genPt->Scale(totXSxBr/t->GetEntries());
  cout<<"total: "<<genPt->Integral()<<endl;

  for(int i=1; i<=genPt->GetNbinsX(); i++){
    genPt->SetBinContent( i, genPt->GetBinContent(i)/genPt->GetBinWidth(i) );
    genPt->SetBinError( i, genPt->GetBinError(i)/genPt->GetBinWidth(i) );
  }
  for(int i=1; i<=genPtBinned->GetNbinsX(); i++){
    genPtBinned->SetBinContent( i, genPtBinned->GetBinContent(i)/genPtBinned->GetBinWidth(i) );
    genPtBinned->SetBinError( i, genPtBinned->GetBinError(i)/genPtBinned->GetBinWidth(i) );
  }
  TGraphAsymmErrors genPtGraph(genPtBinned);
  genPtGraph.SetMarkerStyle(20);
  cout<<"Starting inversions, make sure it finishes!"<<endl;
  for(int i=0; i<genPtGraph.GetN(); i++){
    double x, Y;
    genPtGraph.GetPoint(i, x, Y);
    int ind = genPt->FindBin(x);
    double newx = invert(genPt, ind, Y);
    cout<<i<<" "<<newx<<endl;
    genPtGraph.GetX()[i] = newx;
    genPtGraph.GetEXhigh()[i] -= newx-x;
    genPtGraph.GetEXlow()[i] += newx-x;
  }
  cout<<"Inversions finished!"<<endl;

  TFile out(outputName,"recreate");
  genPtGraph->Write("genPtLargeBinsGraph");
  genPtBinned->Write("genPtLargeBins");
  genPt->Write();
  out.Close();
}

int main(int argc, const char** argv){
  genXSupsilon2S(argv[1],argv[2]);
}
