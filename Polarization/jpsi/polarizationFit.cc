///////
// This program is meant to extract the polarization parameters from the data.
// It takes the input parameters from the previously run "extract*" programs
// and then runs the constrained fits to extract the polarization of the prompt
// and non-prompt components.
// NOTE: For now you need to change the input file names in the code and recompile right now.... 
//       Will write a nice commandline parser later.
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
#include "RooNLLVar.h"
#include "RooMinuit.h"
#include "RooProfileLL.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TCanvas.h"

int main(int argc, char** argv) {
  using namespace JPsiPolarization;

  bool doprompt(true), dononprompt(true), dobkg(true),dopol(true),pereverr(false),newdata(false),domass(false),dolifetime(false);

  if(argc < 2) {
    std::cout << "Usage: extractMassShape /path/to/sample1.root /path/to/sample2.root ..." << std::endl;
    std::cout << "You must ensure that the MC samples are properly weighted so they can be combined." << std::endl;
    return 1;
  }  

  for( int i=0;i < argc; ++i ) { 
    if(std::string(argv[i]).find("--noNonPrompt") != std::string::npos) dononprompt = false;
    if(std::string(argv[i]).find("--noPrompt") != std::string::npos) doprompt = false;
    if(std::string(argv[i]).find("--noBackground") != std::string::npos) dobkg = false;
    if(std::string(argv[i]).find("--noPol") != std::string::npos) dopol = false;
    if(std::string(argv[i]).find("--perEventErrors") != std::string::npos) pereverr = true;
    if(std::string(argv[i]).find("--newTTree") != std::string::npos) newdata = true;
    if(std::string(argv[i]).find("--noMass") != std::string::npos) domass = false;
    if(std::string(argv[i]).find("--noLifetime") != std::string::npos) dolifetime = false;
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
  RooRealVar HLT_Mu0_TkMu0_Jpsi("HLT_Mu0_TkMu0_Jpsi","Passes HLT_Mu0_TkMu0_Jpsi Trigger",0.5,1.5);
  RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,2.5);//0=GG,1=GT,2=TT
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
  TFile *fitInput = new TFile("jPsiFit.root","UPDATE");
  TFile *output_cs = new TFile("jPsiFitFinal_cs.root","RECREATE");
  TFile *output_hx = new TFile("jPsiFitFinal_hx.root","RECREATE");
  RooDataSet *data = NULL;
  CompositeModelBuilder *hx, *cs;
  
  for( unsigned int arg = 1; arg < argc; ++arg ) {
    if(std::string(argv[arg]).find(".root") != std::string::npos) {
      std::cout << "Adding: " << argv[arg] << std::endl;
      samples->Add(argv[arg]);
    }
  }
  
  data = new RooDataSet("data","Concatenated Samples",samples,varlist);

  for(int yBin = 0; yBin < 1 /*jpsi::kNbRapForPTBins*/; ++yBin) {
    for(int ptBin = 3; ptBin < 4 /*jpsi::kNbPTBins[yBin+1]*/; ++ptBin) {       
      std::stringstream binName,cutString;
      binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1]
		<< " && JpsiMass > " << jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass 
		<< " && JpsiMass < " << jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass			
		<< " && Jpsict < .1";

      JpsiMass.setRange("lowBand",2.7,
			jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgLow);
      JpsiMass.setRange("signalRegion",
			jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass,
			jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass);
      JpsiMass.setRange("highBand",
			jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgHigh,3.5);

      std::cout << cutString.str() << std::endl;
      
      TDirectory *current_in,*current_out_hx,*current_out_cs; 
      if(current_in = fitInput->GetDirectory(binName.str().c_str())) {
            
	RooAbsData *thisBin = data->reduce((cutString.str()).c_str());
	
	hx = new CompositeModelBuilder("HX","signalRegion");
	cs = new CompositeModelBuilder("CS","signalRegion");
	
	hx->setUsePrompt(doprompt);
	hx->setUseNonPrompt(dononprompt);       
	hx->setUseBkg(dobkg); 
	hx->setUseMass(domass);
	hx->setUseLifetime(dolifetime);
	hx->setUsePol(dopol);
	hx->setUseAcceptanceMaps(dopol);

	cs->setUsePrompt(doprompt);
	cs->setUseNonPrompt(dononprompt);
	cs->setUseBkg(dobkg);
	cs->setUseMass(domass);
	cs->setUseLifetime(dolifetime);
	cs->setUsePol(dopol);
	cs->setUseAcceptanceMaps(dopol);

	cs->loadParameters(*current_in);
	hx->loadParameters(*current_in);      
	
	if(pereverr) {
	  cs->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS);
	  hx->initModel(JpsiMass,Jpsict,JpsictErr,costh_HX,phi_HX);	
	} else {
	  cs->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);
	  hx->initModel(JpsiMass,Jpsict,costh_HX,phi_HX);	
	}
	/*
	cs->getLifetimeModel()->unfix("MeanPrompt1");
	cs->getLifetimeModel()->unfix("CoefPrompt1");
	cs->getLifetimeModel()->unfix("CoefPrompt2");
	*/
	//cs->getLifetimeModel()->unfix("SigmaPrompt1");
	/*
	cs->getLifetimeModel()->unfix("SigmaPrompt2");

	hx->getLifetimeModel()->unfix("MeanPrompt1");
	hx->getLifetimeModel()->unfix("CoefPrompt1");
	hx->getLifetimeModel()->unfix("CoefPrompt2");
	*/
	//hx->getLifetimeModel()->unfix("SigmaPrompt1");
	/*
	hx->getLifetimeModel()->unfix("SigmaPrompt2");
	*/
	std::cout << "Collins-Soper model before fit:" << std::endl;
	cs->Print();
		
	RooAbsReal *csNLL = NULL;

	std::cout << "Fitting CS..." << std::endl;
	if(!(ptBin==0 && (yBin == 0))) {
	  csNLL = cs->model()->createNLL(*thisBin,RooFit::NumCPU(2),
					 RooFit::Extended(true),
					 //RooFit::Minos(1),
					 RooFit::InitialHesse(true),
					 RooFit::Range("signalRegion"),
					 RooFit::Strategy(2),
					 RooFit::PrintEvalErrors(-1));

	  if(pereverr) 
	    cs->model()->fitTo(*thisBin,RooFit::NumCPU(2),
			       RooFit::Extended(true),RooFit::SumW2Error(true),
			       //RooFit::Minos(1),
			       RooFit::InitialHesse(true),
			       RooFit::Range("signalRegion"),
			       RooFit::Strategy(2),
			       RooFit::PrintEvalErrors(-1));
			       //RooFit::ConditionalObservables(RooArgSet(JpsictErr)));
	  else
	    cs->model()->fitTo(*thisBin,RooFit::NumCPU(2),RooFit::Timer(true),
			       RooFit::Extended(true),RooFit::SumW2Error(true),
			       //RooFit::Minos(1),
			       RooFit::InitialHesse(true),
			       RooFit::Range("signalRegion"),
			       RooFit::Strategy(2),
			       RooFit::PrintEvalErrors(-1));
	}
	
	if(csNLL && false ) {
	  RooArgSet *vars = cs->model()->getParameters(RooArgSet());

	  RooProfileLL* prof = (RooProfileLL*)csNLL->createProfile(RooArgSet(*vars->find("promptlambda_theta_CS"),
									     *vars->find("promptlambda_phi_CS")));
	    
	  prof->getVal();

	  RooMinuit* csMinuit = prof->minuit();
	  csMinuit->setStrategy(2);
	  csMinuit->setPrintLevel(1);

	  double cl68 = TMath::ChisquareQuantile(.68,2)/2.0; // 68% CL for 2 degrees of freedom

	  csMinuit->migrad();
	  
	  TCanvas *c = new TCanvas("contours","",500,500);

	  RooPlot* thePlot = csMinuit->contour(*(RooRealVar*)vars->find("promptlambda_theta_CS"),
					       *(RooRealVar*)vars->find("promptlambda_phi_CS"),
					       TMath::Sqrt(2.0*cl68));

	  c->cd();
	  thePlot->Draw();

	  c->Print("cscontours.root");
	  delete thePlot;
	  delete c;
	}

	std::cout << "Collins-Soper model after fit:" << std::endl;
	cs->Print();

	std::cout << std::endl << "Helicity Frame model before fit:" << std::endl;
	hx->Print();

	std::cout << "Fitting HX..." << std::endl;

	RooAbsReal *hxNLL = NULL;

	if(!(ptBin==0 && (yBin == 0))) {
	  hxNLL = hx->model()->createNLL(*thisBin,RooFit::NumCPU(2),
					 RooFit::Extended(true),
					 //RooFit::Minos(1),
					 RooFit::InitialHesse(true),
					 RooFit::Range("signalRegion"),
					 RooFit::Strategy(2),
					 RooFit::PrintEvalErrors(-1));

	  if(pereverr)
	    hx->model()->fitTo(*thisBin,RooFit::NumCPU(2),
			       RooFit::Extended(true),RooFit::SumW2Error(true),
			       //RooFit::Minos(1),
			       RooFit::InitialHesse(true),
			       RooFit::Range("signalRegion"),			   
			       RooFit::Strategy(2),
			       RooFit::PrintEvalErrors(-1));
	  // RooFit::ConditionalObservables(RooArgSet(JpsictErr)));
	  else
	    hx->model()->fitTo(*thisBin,RooFit::NumCPU(2),RooFit::Timer(true),
			       RooFit::Extended(true),RooFit::SumW2Error(true),
			       //RooFit::Minos(1),
			       RooFit::InitialHesse(true),
			       RooFit::Range("signalRegion"),			   
			       RooFit::Strategy(2),
			       RooFit::PrintEvalErrors(-1));
	}
	
	if(hxNLL&&false) {
	  RooArgSet *vars = hx->model()->getParameters(RooArgSet());

	  RooProfileLL* prof = (RooProfileLL*)hxNLL->createProfile(RooArgSet(*vars->find("promptlambda_theta_HX"),
									     *vars->find("promptlambda_phi_HX")));
	    
	  prof->getVal();

	  RooMinuit* hxMinuit = prof->minuit();
	  hxMinuit->setStrategy(2);
	  hxMinuit->setPrintLevel(1);

	  double cl68 = TMath::ChisquareQuantile(.68,2)/2.0; // 68% CL for 2 degrees of freedom

	  hxMinuit->migrad();
	  
	  TCanvas *c = new TCanvas("contours","",500,500);

	  RooPlot* thePlot = hxMinuit->contour(*(RooRealVar*)vars->find("promptlambda_theta_HX"),
					       *(RooRealVar*)vars->find("promptlambda_phi_HX"),
					       TMath::Sqrt(2.0*cl68));

	  c->cd();
	  thePlot->Draw();

	  c->Print("hxcontours.root");

	  
	  delete thePlot;
	  delete c;
	}


	std::cout << std::endl << "Helicity Frame model after fit:" << std::endl;
	hx->Print();

	current_out_hx = output_hx->mkdir(binName.str().c_str());
	current_out_cs = output_cs->mkdir(binName.str().c_str());
	

	cs->saveParameter("nPrompt",*current_out_cs);
	cs->saveParameter("nNonPrompt",*current_out_cs);
	cs->saveParameter("nBackground",*current_out_cs);
	
	hx->saveParameter("nPrompt",*current_out_hx);
	hx->saveParameter("nNonPrompt",*current_out_hx);
	hx->saveParameter("nBackground",*current_out_hx);
	
	cs->saveParameters(*current_out_cs);
	hx->saveParameters(*current_out_hx);
      
	delete cs;
	delete hx;
	delete thisBin;
      }
    }
  }

  output_hx->Write();
  output_hx->Close();

  output_cs->Write();
  output_cs->Close();

  fitInput->Close();

  delete fitInput;
  delete output_hx;
  delete output_cs;
  delete data;
  delete samples;
  return 0;
}
