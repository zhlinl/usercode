
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
#include "RooNDKeysPdf.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"

int main(int argc, char** argv) {
  using namespace JPsiPolarization;

  bool pereverr(false),newdata(false);

  if(argc < 2) {
    std::cout << "Usage: ./extractBackgroundShapes /path/to/data_sample.root" << std::endl;
    std::cout << "One file accepted as input." << std::endl;
    return 1;
  }  

  for( int i=0;i < argc; ++i ) { 
    if(std::string(argv[i]).find("--perEventErrors") != std::string::npos) 
      { 
	std::cout << "using per event errors." << std::endl;
	pereverr = true;
      }
    if(std::string(argv[i]).find("--newTTree") != std::string::npos) {
      std::cout << "Running on newly formatted TTrees." << std::endl;
      newdata = true;
    }
  }

  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
  RooRealVar JpsictErr("JpsictErr","Error on l_{J/#psi} [mm]",1e-4,(1-1e-6));
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
  RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,2.5);//0=GG,1=GT,2=TT
  RooRealVar HLT_Mu0_TkMu0_Jpsi("HLT_Mu0_TkMu0_Jpsi","Passes HLT_Mu0_TkMu0_Jpsi Trigger",0.5,1.5); // restrict to passing j/psis
  RooRealVar MCweight("MCweight","MCweight",0,1.1);
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  if(pereverr) 
    varlist.add(JpsictErr);
  if(newdata)
    varlist.add(HLT_Mu0_TkMu0_Jpsi);
  //varlist.add(MCweight);
  //varlist.add(MCType_idx);  
  varlist.add(JpsiType_idx);

  TChain *samples = new TChain("data");
  TChain *mcsamples = new TChain("data");
  TFile *output = new TFile("jPsiFit.root","UPDATE");  
  RooDataSet *data = NULL, *mc=NULL;
  CompositeModelBuilder* theModel;
  
  //yeah it's ugly but it works for now
  std::cout << "Adding: " << argv[1] << std::endl;
  samples->Add(argv[1]);
    
  data = new RooDataSet("data","Data Samples",samples,varlist);


  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
    for(int ptBin = 0; ptBin < jpsi::kNbPTBins[yBin+1]; ++ptBin) {       
      std::stringstream binName,cutString;
      binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1];
      
      CompositeModelBuilder* modelHX = new CompositeModelBuilder("HX");
      CompositeModelBuilder* modelCS = new CompositeModelBuilder("CS");
      CompositeModelBuilder* bkgMassFit = new CompositeModelBuilder();
      CompositeModelBuilder* bkgLifetimeFit = new CompositeModelBuilder();

      

      JpsiMass.setRange("lowBand",2.7,
			jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgLow);
      JpsiMass.setRange("signalRegion",
			jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass,
			jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass);
      JpsiMass.setRange("highBand",
			jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgHigh,3.5);

      std::cout << cutString.str() << std::endl;
      
      TDirectory *current; 
      if(!(current = output->GetDirectory(binName.str().c_str()))) {
	current = output->mkdir(binName.str().c_str());
      }
      
      RooAbsData *thisBin = data->reduce(cutString.str().c_str());
      
      std::cout << "Number of Entries = " << thisBin->sumEntries() << std::endl;

      //Extract the shape and normalization of the background and save them

      bkgMassFit->setUseLifetime(false);
      bkgMassFit->setUsePol(false);

      bkgMassFit->setUsePrompt(false);
      bkgMassFit->setUseNonPrompt(false);     

      bkgMassFit->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);

      bkgMassFit->unfix("nBackground");
      bkgMassFit->getMassModel()->unfix("bkgLambda");

      bkgMassFit->Print();

      bkgMassFit->getBkgModel()->fitTo(*thisBin,RooFit::NumCPU(2),RooFit::Timer(true),
				   RooFit::Extended(true),RooFit::SumW2Error(true),
				   RooFit::PrintEvalErrors(-1));
      /* // it seems that roofit was doing the right thing after all
      RooAbsReal* bkgPdfVar = bkgMassFit->getBkgModel()->createIntegral(RooArgSet(JpsiMass),
								    RooFit::NormSet(JpsiMass),
								    RooFit::Range("signalRegion"));

      double bkgPdfVal = bkgPdfVar->getVal();

      //std::cout << "Total Number of events: " << bkgInt->expectedEvents() << std::endl;
	
      std::cout << "# of Background under signal: " << bkgPdfVal << std::endl;

      //Well, let's just hack it together then...

      double peakbkg = bkgPdfVal*bkgMassFit->backgroundNorm()->getVal();
      double peakbkg_err = bkgPdfVal*bkgMassFit->backgroundNorm()->getError();
      */
      bkgMassFit->fix("nBackground");
      bkgMassFit->saveParameter("nBackground",*current);
      bkgMassFit->getMassModel()->fix("bkgLambda"); // fix background shape parameter
      bkgMassFit->saveParameters(*current);

      //Now do the lifetime fit

      bkgLifetimeFit->setUseMass(false);
      bkgLifetimeFit->setUsePol(false);
      bkgLifetimeFit->setUseNonPrompt(false);
      
      bkgLifetimeFit->loadParameters(*current); // get resolution function parameters from file
      
      if(pereverr) {
	bkgLifetimeFit->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS);
      } else {
	bkgLifetimeFit->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);
      }

      //      bkgLifetimeFit->getLifetimeModel()->unfix("SigmaPrompt1");
      bkgLifetimeFit->getLifetimeModel()->unfix("Bkg");

      bkgLifetimeFit->Print();

      if(pereverr) 
	bkgLifetimeFit->getLifetimeModel()->background()->fitTo(*thisBin,RooFit::NumCPU(2),RooFit::Timer(true),
								RooFit::Extended(false),RooFit::SumW2Error(true),
								RooFit::PrintEvalErrors(-1),
								RooFit::ConditionalObservables(RooArgSet(JpsictErr)));
      else
	bkgLifetimeFit->getLifetimeModel()->background()->fitTo(*thisBin,RooFit::NumCPU(2),RooFit::Timer(true),
								RooFit::Extended(false),RooFit::SumW2Error(true),
								RooFit::PrintEvalErrors(-1));
      
      //      bkgLifetimeFit->getLifetimeModel()->fix("SigmaPrompt1");
      bkgLifetimeFit->getLifetimeModel()->fix("Bkg"); // fix background shape parameter
      
      bkgLifetimeFit->saveParameters(*current);

      //Next determine bkg shape maps and save them
      modelHX->setUseLifetime(false);
      modelHX->setUseMass(false);
      modelHX->setUsePol(false);
      modelHX->setUsePrompt(true);
      modelHX->setUseNonPrompt(true);

      modelCS->setUseLifetime(false);
      modelCS->setUseMass(false);
      modelCS->setUsePol(false);
      modelCS->setUsePrompt(true);      
      modelCS->setUseNonPrompt(true);      

      

      TH2F* bkgCS = (TH2F*)thisBin->createHistogram("bkgShape_CS",costh_CS,
						    RooFit::Binning(jpsi::kNbBinsCosT,jpsi::cosTMin,jpsi::cosTMax),
						    RooFit::YVar(phi_CS,RooFit::Binning(jpsi::kNbBinsPhiPol,jpsi::phiPolMin,jpsi::phiPolMax)));
      
      TH2F* bkgHX = (TH2F*)thisBin->createHistogram("bkgShape_HX",costh_HX,
						    RooFit::Binning(jpsi::kNbBinsCosT,jpsi::cosTMin,jpsi::cosTMax),
						    RooFit::YVar(phi_HX,RooFit::Binning(jpsi::kNbBinsPhiPol,jpsi::phiPolMin,jpsi::phiPolMax)));
      
      modelCS->setBackgroundShapeHist(bkgCS);
      modelHX->setBackgroundShapeHist(bkgHX);

      modelCS->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);
      modelHX->initModel(JpsiMass,Jpsict,costh_HX,phi_HX);
      /*
      if(thisBin->sumEntries()) {
	std::cout << "Making CS NDKeysPdf..." << std::endl;
	RooNDKeysPdf* bkgCSpdf = new RooNDKeysPdf("bkgShapePdf_CS","Estimated CS Background Shape Pdf",
						  costh_CS,phi_CS,*(RooDataSet*)thisBin);
	
	bkgCSpdf->fixShape(true);
	
	std::cout << "Making HX NDKeysPdf..." << std::endl;
	RooNDKeysPdf* bkgHXpdf = new RooNDKeysPdf("bkgShapePdf_HX","Estimated HX Background Shape Pdf",
						  costh_HX,phi_HX,*(RooDataSet*)thisBin);
	bkgHXpdf->fixShape(true);
	
	current->Add(bkgCSpdf);
	current->Add(bkgHXpdf);
      }
      */

      modelCS->saveParameters(*current);
      modelHX->saveParameters(*current);
      
      delete bkgMassFit;
      delete bkgLifetimeFit;
      //delete bkgPdfVar;
      delete modelHX;
      delete modelCS;
      delete bkgCS;
      delete bkgHX;      
      delete thisBin;
    }
  }

  output->Write();
  output->Close();

  delete output;
  delete data;
  delete samples;
  return 0;
}
