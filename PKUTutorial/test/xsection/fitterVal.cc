#include "oniafitter.C"
#include "RooCategory.h"
#include "RooThresholdCategory.h"
#include "RooSimPdfBuilder.h"
#include "RooSimultaneous.h"
#include "TStopwatch.h"

#include <iostream>
#include "TRandom3.h"


/*
  systematicFitter parameters
  nToys : number of times to do the fit
  genParams : file name for a file which contains the generation parameters 
              for the pdf.  This will be passed to RooArgSet::readFromFile()
              so it needs to be properly formated.
  fitParams : file containing the initial values for pdf parameters for the fit,
              Again using RooArgSet::readFromFile().
  outFilename : file name to store the results.  This will be a text file which
                can be parsed to get systematic results.
  randomizeTail : if the alpha and npow parameters in the CB tail should be
                  selected at random based on their covariance which is 
                  hardcoded in the routine.
  binning : 1,2 or 3 depending on if you want 1S binning, 2S binning, 
            or 3S binning respectively.

examples:

To validate the fitter using Toy MC
systematicFitter(1000, "genParams.txt", "fitParams.txt", 
                 "FitterValidation.asc", false, 1, "datafile.root");

to get the systematic for the fixed tail parameters
systematicFitter(1000, "", "fitParams.txt", "CBTailParamFits.asc", 
                 true, 1, "MuOnia.root");

to do a single fit to the data in the file
systematicFitter(1, "", "fitParams.txt", "/dev/null", false, 1, "MuOnia.root");

 */
void fitterVal(int nToys, TString genParams, TString fitParams,
	       TString outFilename = "validation.asc",
	       bool randomizeTail = false, int binning = 0) {
  TStopwatch timer;
  timer.Start();
  static TRandom3 rnd(0);
  RooWorkspace oniaWS("oniaWS", "Quarkonia mass");

  dirname_ = "yieldUpsilonTree";

  buildPdf(oniaWS);
  RooRealVar * invariantMass = oniaWS.var("invariantMass");
  RooAbsPdf * pdf = oniaWS.pdf("pdf");

  RooRealVar upsPt("upsPt", "upsPt", 0.0, 100.0, "GeV/c");

  RooDataSet * data = 0;

  RooThresholdCategory ptbint("ptbin", "ptbin", upsPt);
  RooCategory ptbin("ptbin", "ptbin");

  double const * ptBinArr = ptdiff_1ybin_ptBinArr1S;
  int nbins = sizeof(ptdiff_1ybin_ptBinArr1S)/sizeof(double);
  if (binning == 2) {
    ptBinArr = ptdiff_1ybin_ptBinArr2S;
    nbins = sizeof(ptdiff_1ybin_ptBinArr2S)/sizeof(double);
  } else if (binning == 3) {
    ptBinArr = ptdiff_1ybin_ptBinArr2S;
    nbins = sizeof(ptdiff_1ybin_ptBinArr3S)/sizeof(double);
  }
  char buf[100];
  std::cout << "n pt bins: " << nbins << '\n';
  for (int bin = 1; bin < nbins; ++bin) {
    sprintf(buf, "%0.1f-%0.1f", ptBinArr[bin-1], ptBinArr[bin]);
    ptbint.addThreshold(ptBinArr[bin], buf);
    ptbin.defineType(buf, bin);
  }
  //ptbint.Print("v");


  invariantMass->Print();
  invariantMass->setRange(8.0, 14.0);
  invariantMass->Print();
  RooArgList cols(*invariantMass, ptbin);
  RooSimPdfBuilder simb(*pdf);
  RooArgSet * config = simb.createProtoBuildConfig();

  config->setStringValue("physModels", pdf->GetName());

  TString splitStr(ptbin.GetName());
  splitStr += '(';
  for (int t = 1; t < ptbin.numTypes(); ++t) {
    ptbin.setIndex(t);
    splitStr += ptbin.getLabel();
    if (t < ptbin.numTypes() - 1)
      splitStr += ',';
  }
  splitStr += ')';
  std::cout << splitStr << '\n';
  config->setStringValue("splitCats", splitStr);

  splitStr = ptbin.GetName();
  splitStr += " : nsig1,nsig2,nsig3,nbkgd,bkg_a1,bkg_a2";
  config->setStringValue(pdf->GetName(), splitStr);

  RooSimultaneous * simpdf = simb.buildPdf(*config, cols);

  RooAbsPdf * fitter = simpdf;
  if (binning == 0) fitter = pdf;

  fitter->Print();
  RooArgSet * pars = fitter->getParameters(cols);
  if (genParams.Length() > 0) {
    pars->readFromFile(genParams.Data(), 0, 0, kTRUE);
  }
  RooArgSet * trueVals = dynamic_cast<RooArgSet *>(pars->snapshot());
  trueVals->Print("v");

  //return;

  // RooDataSet * weights = 0;

  // RooRealVar weight("weight", "weight", 1.0);
  // if (weightDataSet.Length() > 0) {
  //   weights = RooDataSet::read(weightDataSet.Data(), RooArgList(weight),
  // 			       "Q");
  // }

  ofstream fout(outFilename.Data(), ios_base::out|ios_base::trunc);
  for (int toy = 0; toy < nToys; ++toy) {
    std::cout << "Toy " << toy+1 << " of " << nToys << '\n';
    if (genParams.Length() > 0) {
      pars->readFromFile(genParams.Data(), 0, 0, kTRUE);
      //pars->Print("v");
      double expectedN = 0;
      TString fitterClass(fitter->ClassName());
      if (fitterClass == "RooSimultaneous") {
	for (int pt = 1; pt < ptbin.numTypes(); ++pt) {
	  ptbin.setIndex(pt);
	  expectedN += fitter->expectedEvents(RooArgSet(*invariantMass,ptbin));
	}
      } else {
	expectedN = fitter->expectedEvents(RooArgSet(*invariantMass));
      }
      std::cout << fitterClass << " expected N: " << expectedN << '\n';
      data = fitter->generate(RooArgSet(*invariantMass,ptbin),
			      Extended(kTRUE), NumEvents(int(expectedN+0.5)));
    }
    // if (weights) {
    //   RooDataSet * tmpdata = dynamic_cast<RooDataSet *>(weights->reduce(EventRange(0, data->numEntries())));
    //   data->merge(tmpdata);
    //   delete tmpdata;
    //   tmpdata = new RooDataSet("data", "data", RooArgSet(*invariantMass,weight),
    // 			       WeightVar("weight"), Import(*data));
    //   delete data;
    //   data = tmpdata;
    // } 

    pars->readFromFile(fitParams.Data(), 0, 0, kTRUE);
    std::cout << "data contains: " << data->numEntries() << " total weight: "
	      << data->sumEntries() << '\n';

    if (randomizeTail) {
      float const Q[2][2] = {{0.888,  0.4595},
			     {-0.4595, 0.888}};
      double saprime = rnd.Gaus(0, 0.022);
      double snprime = rnd.Gaus(0, 0.123);
      double alpha = pars->getRealValue("alpha", 1.88);
      double npow = pars->getRealValue("npow", 1.33);
      double aprime = Q[0][0]*alpha + Q[1][0]*npow + saprime;
      double nprime = Q[0][1]*alpha + Q[1][1]*npow + snprime;
      std::cout << "alpha: " << alpha << " npow: " << npow << '\n';
      alpha = Q[0][0]*aprime + Q[0][1]*nprime;
      npow  = Q[1][0]*aprime + Q[1][1]*nprime;
      pars->setRealValue("alpha", alpha);
      pars->setRealValue("npow", npow);
      std::cout << "saprime: " << saprime << " snprime: " << snprime << '\n';
      std::cout << "alpha: " << alpha << " npow: " << npow << '\n';
    }

    RooFitResult * fr = fitter->fitTo(*data, Extended(kTRUE),
				      //NumCPU(2),
				      //SumW2Error(kFALSE),
				      PrintLevel(((nToys==1) ? 1 : -1)),
				      Minos(kFALSE),
				      Save(kTRUE));
    int covQual = fr->covQual();
    double edm = fr->edm();

    // if (data->isWeighted()) {
    //   delete fr;
    //   //pars->readFromFile(fitParams.Data(), 0, 0, kTRUE);
    //   fr = fitter->fitTo(*data, Extended(kTRUE),
    //  		      SumW2Error(kTRUE), Minos(kFALSE),
    //  		      PrintLevel(-1),
    //  		      //NumCPU(2),
    //  		      Save(kTRUE));
		
    // }

    if ( (fr) && ( (covQual == 3) || (covQual == -1) ) && (edm <= 0.001) ) {
      double avgWgt = data->sumEntries()/data->numEntries(), trueVal;
      TIter nextPar(fr->floatParsFinal().createIterator());
      RooRealVar * par = 0;
      //RooRealVar * parInit = 0;
      fout << "nll " << fr->minNll() << ' ';
      while ((par = (RooRealVar *)nextPar())) {
	trueVal = trueVals->getRealValue(par->GetName(), -9999, kFALSE);
	if ( data->isWeighted() && (par->GetName()[0] == 'n') ) {
	    trueVal *= avgWgt;
	}
	fout << par->GetName() << ' '
	     << par->getVal() << ' '
	     << par->getError() << ' '
	     << trueVal << ' ';
      }
      fout << '\n';
    }
    
    if (nToys == 1) {
      //pars->writeToFile("systematicFitterPars.txt");
      RooPlot * mf = invariantMass->frame(Bins(120));

      //fitter->Print("v");
      data->plotOn(mf);
      fitter->plotOn(mf, ProjWData(RooArgSet(ptbin), *data));
      fitter->plotOn(mf, ProjWData(RooArgSet(ptbin), *data), 
		     Components("bkg*"), LineColor(kRed), LineStyle(kDashed));

      mf->Draw();
      fr->Print("v");
      //data->write("SingleDataset.out.asc");
    }
    delete fr;
    if (genParams.Length() > 0) {
      delete data;
      data = 0;
    }
    std::cout.flush();
  }

  fout.close();
  timer.Stop();
  timer.Print();
  // std::cout << "Real time elapsed: " << timer.RealTime() << " s\n"
  // 	    << "CPU time consumed: " << timer.CpuTime() << " s\n";
  std::cout.flush();

}
