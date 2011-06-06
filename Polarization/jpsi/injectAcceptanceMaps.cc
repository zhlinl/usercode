
///////
// This program is meant to extract the background shapes from the data side bands. 
// Only accepts one file, weight doesn't matter.
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
#include "RooFitResult.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"

int main(int argc, char** argv) {
  using namespace JPsiPolarization;

  if(argc != 3) {
    std::cout << "Usage: ./injectAcceptanceMaps /path/to/promptAcceptance.root /path/to/nonPromptAcceptance.root" << std::endl;
    std::cout << "Two files accepted as input." << std::endl;
    return 1;
  }  

  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5); 
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
  
  TFile *output = new TFile("jPsiFit.root","UPDATE");
  TFile *promptMaps = NULL, *nonPromptMaps = NULL;
 
  std::cout << "Loading acceptance maps from: " << argv[1] << " and " << argv[2] << std::endl;
  promptMaps = TFile::Open(argv[1],"UPDATE");
  nonPromptMaps = TFile::Open(argv[2],"UPDATE");
  

  for(int yBin = 0; yBin < 1/*jpsi::kNbRapForPTBins*/; ++yBin) {
    for(int ptBin = 3; ptBin < 4 /*jpsi::kNbPTBins[yBin+1]*/; ++ptBin) {       
      std::stringstream binName,cutString;
      binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1];

      CompositeModelBuilder* modelHX = new CompositeModelBuilder("HX");
      CompositeModelBuilder* modelCS = new CompositeModelBuilder("CS");
           
      TDirectory *current; 
      if(!(current = output->GetDirectory(binName.str().c_str()))) {
	current = output->mkdir(binName.str().c_str());
      }
      
      modelHX->setUseLifetime(false);
      modelHX->setUsePol(false);
      modelHX->setUseMass(false);
      modelHX->setUseBkg(false);
      
      modelCS->setUseLifetime(false);
      modelCS->setUsePol(false);
      modelCS->setUseMass(false);
      modelCS->setUseBkg(false);
      
      modelCS->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);

      TH2F *promptMap_CS,*nonPromptMap_CS;
      TH2F *promptMap_HX,*nonPromptMap_HX;
      
      std::stringstream nameHXp,nameCSp;
      std::stringstream nameHXnp,nameCSnp;      
      
      nameHXp << "hAcc2D_Onia_HX_pT" << ptBin+1 << "_rap" << yBin+1;
      nameCSp << "hAcc2D_Onia_CS_pT" << ptBin+1 << "_rap" << yBin+1;
      //nameHXnp << "hAcc2D_Onia_HX_pT" << ptBin+1 << "_rap" << yBin+1 << "_NP";
      //nameCSnp << "hAcc2D_Onia_CS_pT" << ptBin+1 << "_rap" << yBin+1 << "_NP";
      
      promptMaps->GetObject(nameHXp.str().c_str(),promptMap_HX);
      //nonPromptMaps->GetObject(nameHXnp.str().c_str(),nonPromptMap_HX);
      
      promptMaps->GetObject(nameCSp.str().c_str(),promptMap_CS);
      //nonPromptMaps->GetObject(nameCSnp.str().c_str(),nonPromptMap_CS);
      
      modelCS->setPromptAccHist(promptMap_CS);
      //modelCS->setNonPromptAccHist(nonPromptMap_CS);
      
      modelHX->setPromptAccHist(promptMap_HX);
      //modelHX->setNonPromptAccHist(nonPromptMap_HX);

      modelCS->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);
      modelHX->initModel(JpsiMass,Jpsict,costh_HX,phi_HX);
          
      modelCS->saveParameters(*current);
      modelHX->saveParameters(*current);
      
      delete modelHX;
      delete modelCS;      
    }
  }

  //output->Write();
  output->Close();
  promptMaps->Close();
  nonPromptMaps->Close();
  
  delete promptMaps;
  delete nonPromptMaps;
  delete output;
  return 0;
}
