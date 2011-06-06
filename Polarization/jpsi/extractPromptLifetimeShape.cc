
///////
// This program is meant to extract the prompt and non-prompt lifetime shapes from the monte carlo.
// It expects to be given a list of files with !!!appropriate weights!!!
// \author Lindsey Gray (UW Madison)
//////

#include <iostream>
#include <sstream>

//J/Psi common vars
#include "commonVar.h"

//Fitting routine
#include "CompositeModelBuilder.h"

// RooFit Includes
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"

int main(int argc, char** argv) {
  using namespace JPsiPolarization;
  
  bool pereverr = false;
  bool newdata = true;

  if(argc < 2) {
    std::cout << "Usage: extractMassShape /path/to/sample1.root /path/to/sample2.root ..." << std::endl;
    std::cout << "You must ensure that the MC samples are properly weighted so they can be combined." << std::endl;
    return 1;
  }  

  for( int i=0;i < argc; ++i ) { 
    if(std::string(argv[i]).find("--perEventErrors") != std::string::npos) pereverr = true;    
    if(std::string(argv[i]).find("--newTTree") != std::string::npos) newdata = true;
  }

  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
  RooRealVar JpsictErr("JpsictErr","Error on l_{J/#psi} [mm]",1e-3,(1-1e-6));
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
  RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,2.5);//0=GG,1=GT,2=TT
  RooRealVar MCweight("MCweight","MCweight",0,1.1);
  RooRealVar HLT_Mu0_TkMu0_Jpsi("HLT_Mu0_TkMu0_OST_Jpsi","Passes HLT_Mu0_TkMu0_Jpsi Trigger",0.5,1.5);
  
  //Jpsict.setBins(500,"lifetime");

  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  if(pereverr) 
    varlist.add(JpsictErr);
  if(newdata)
    varlist.add(HLT_Mu0_TkMu0_Jpsi);
  varlist.add(MCweight);
  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);

  TChain *samples = new TChain("data");
  TFile *output = new TFile("jPsiFit.root","UPDATE");
  RooDataSet *data = NULL;
  CompositeModelBuilder* theModel;
  
  for( unsigned int arg = 1; arg < argc; ++arg ) {
    if(std::string(argv[arg]).find(".root") != std::string::npos) {
      std::cout << "Adding: " << argv[arg] << std::endl;
      samples->Add(argv[arg]);
    }
  }
  
  data = new RooDataSet("data","Concatentated MC Samples",samples,varlist,0,"MCweight");

  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
    for(int ptBin = 0; ptBin < jpsi::kNbPTBins[yBin+1]; ++ptBin) {       
      std::stringstream binName,cutString;
      binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1];


      std::cout << cutString.str() << std::endl;
      
      TDirectory *current; 
      if(!(current = output->GetDirectory(binName.str().c_str()))) {
	current = output->mkdir(binName.str().c_str());
      }
      std::string promptCut("");

      RooAbsData *thisBinPrompt = data->reduce((cutString.str()+promptCut).c_str());

      RooDataHist *promptHist = new RooDataHist("promptHist","Binned Prompt Lifetime",RooArgSet(Jpsict),*thisBinPrompt);

      theModel = new CompositeModelBuilder();
      
      theModel->setUseBkg(false);
      theModel->setUseNonPrompt(false);
      theModel->setUseMass(false);
      theModel->setUsePol(false);      
      
      if(pereverr)
	theModel->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS);
      else
	theModel->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);
      
      theModel->Print();

      if(pereverr) 
	theModel->getLifetimeModel()->prompt()->fitTo(*thisBinPrompt,RooFit::NumCPU(2),RooFit::Timer(true),
						      RooFit::Extended(false),RooFit::SumW2Error(true),
						      RooFit::Strategy(2),
						      RooFit::InitialHesse(true),
						      RooFit::PrintEvalErrors(-1),
						      RooFit::ConditionalObservables(RooArgSet(JpsictErr)));
      else
	theModel->getLifetimeModel()->prompt()->fitTo(*thisBinPrompt,RooFit::NumCPU(2),RooFit::Timer(true),
						      RooFit::Extended(false),RooFit::SumW2Error(true),
						      RooFit::Strategy(2),
						      RooFit::InitialHesse(true),
						      RooFit::PrintEvalErrors(-1));
      
      theModel->getLifetimeModel()->fix("Prompt1");
      theModel->getLifetimeModel()->fix("Prompt2");
      theModel->getLifetimeModel()->fix("Prompt3");
      
      theModel->saveParameters(*current);
      
      delete theModel;
      delete promptHist;
      delete thisBinPrompt;
    }
  }

  output->Write();
  output->Close();

  delete output;
  delete data;
  delete samples;
  return 0;
}
