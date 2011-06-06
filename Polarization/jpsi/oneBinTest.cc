#include <iostream>

//J/Psi common vars
#include "commonVar.h"

//Fitting routine
#include "CompositeModelBuilder.h"

// RooFit Includes
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"

int main(int argc, char**argv) {
  using namespace JPsiPolarization;

  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar JpsiRap("JpsiRap","#nu",-0.9,0.9);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",7.5,12.5);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
  RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,1.5);//0=GG,1=GT,2=TT
  RooRealVar MCweight("MCweight","MCweight",0,1.1);
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);
  
  TChain *dataTrees = new TChain("data");  
  dataTrees->Add("/Users/lindseygray/jpsi_workdir/data/TTree_pol_Mu0TkMu0Jpsi_dataR_30Aug2010.root");
  //dataTrees->Add("/Users/lindseygray/jpsi_workdir/Spring10/promptJpsi/TTree_pol_Mu0Track0Jpsi_MCprompt.root");

  TFile *accFile = new TFile("accHistos_HLT_Mu0Track0Jpsi_29Aug2010_2.root");
  TFile *outFile = new TFile("output.root","RECREATE");

  TH2F* promptAcceptance = (TH2F*)accFile->Get("hAcc2D_Onia_HX_pT4_rap1");

  RooDataSet *data = new RooDataSet("data","Supplied Data",dataTrees,varlist);

  CompositeModelBuilder *HXbuilder = new CompositeModelBuilder("HX");
  //CompositeModelBuilder *CSbuilder = new CompositeModelBuilder("CS",NULL,NULL,NULL);

  HXbuilder->setUseLifetime(false);
  //CSbuilder->setUseLifetime(false);

  HXbuilder->setUsePol(false);
  //CSbuilder->setUsePol(false);

  HXbuilder->initModel(JpsiMass,Jpsict,costh_HX,phi_HX);
  //CSbuilder->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);

  std::cout << "----MASS ONLY----" << std::endl;
  HXbuilder->Print();

  RooFitResult *res = HXbuilder->model()->fitTo(*data,RooFit::NumCPU(2),RooFit::Timer(true),
						RooFit::Extended(true),RooFit::SumW2Error(true),
						RooFit::PrintLevel(-1));

  //CSbuilder->model()->fitTo(*data,RooFit::Extended(true));
  
  HXbuilder->fix("mass");

  HXbuilder->saveParameters(*outFile);

  /*
  HXbuilder->fix("Nbkg");
  HXbuilder->setUseLifetime(true);

  HXbuilder->initModel(JpsiMass,Jpsict,costh_HX,phi_HX);

  std::cout << "----MASS + CTau----" << std::endl;
  HXbuilder->Print();

  HXbuilder->model()->fitTo(*data,RooFit::NumCPU(2),RooFit::Timer(true),
			    RooFit::Extended(true),RooFit::SumW2Error(true),
			    RooFit::PrintLevel(-1));

  HXbuilder->Print();
  //CSbuilder->Print();

  HXbuilder->fix("lifetime");
  HXbuilder->fix("Nprompt");
  HXbuilder->fix("Nnonprompt");
  HXbuilder->setUsePol(true);

  HXbuilder->setPromptAccHist(promptAcceptance);
  HXbuilder->setNonPromptAccHist(promptAcceptance);

  HXbuilder->initModel(JpsiMass,Jpsict,costh_HX,phi_HX);
  std::cout << "----MASS + CTau + Polarization----" << std::endl;
  HXbuilder->Print();

  HXbuilder->model()->fitTo(*data,RooFit::NumCPU(2),RooFit::Timer(true),
			    RooFit::Extended(false),RooFit::SumW2Error(true),
			    RooFit::PrintLevel(-1));
  */

  HXbuilder->Print();

  outFile->Write();
  outFile->Close();

  delete outFile;
  outFile = TFile::Open("output.root");
  
  delete HXbuilder;

  HXbuilder = new CompositeModelBuilder("HX");
  //CompositeModelBuilder *CSbuilder = new CompositeModelBuilder("CS",NULL,NULL,NULL);

  HXbuilder->setUseLifetime(false);
  //CSbuilder->setUseLifetime(false);

  HXbuilder->setUsePol(false);
  //CSbuilder->setUsePol(false);

  HXbuilder->loadParameters(*outFile);

  HXbuilder->initModel(JpsiMass,Jpsict,costh_HX,phi_HX);
  //CSbuilder->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);

  HXbuilder->Print();
  
  delete accFile;
  delete HXbuilder;
  //delete CSbuilder;
  delete data;
  delete dataTrees;
  delete outFile;

  return 0;
}
