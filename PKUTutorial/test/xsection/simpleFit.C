#include "oniafitter.C"

void simpleFit() {


  /// input and output files
  const TString finput  = "upsilonYieldWeighted_nominal.root";
  const TString foutput = "xsection.root";
  
  cout << "oniafitter processing"
       << "\n\tinput:  \t" << finput
       << "\n\toutput: \t" << foutput
       << "\n\tresults:\t" << figs_
       << endl;
  
  /// desired format of plots to be produced
  ext_ = ".gif";
  
  /// create output file/structure
  ifstream tmp(foutput); bool isofile("true"); if(!tmp)isofile=false; tmp.close(); 
  TFile file(foutput,"update");
  if (!isofile) {
    gDirectory->mkdir(res_);
    gDirectory->mkdir(plt_);
  }
  gDirectory->Cd(res_);

  /// define workspace owing the relevant data and fit objects 
  RooWorkspace* ws = new RooWorkspace("ws","upsilon mass");

  /// define the signal and background PDFs
  buildPdf(*ws);
  
  /// read in the yield data tree
  readData(*ws,finput);

  /// perform a  fit to the full dataset  
  fitYield(*ws,"total",ws->data("data"));

  /// plot the mass and fitted model
  plotMass("mass_raw",*ws);

  //un-comment for performing the weighted fit
  /*  
  /// produce a weighted dataset
  RooDataSet * data_w =  new RooDataSet("data_w", "data_w", *(ws->data("data")->get()), Import(*((RooDataSet *)ws->data("data"))), WeightVar("weight"));

  /// perform the weighted fit 
  fitYield(*ws,"total",data_w);
  plotMass("mass_wei",*ws, data_w, ws->pdf("pdf"));
  */

  //un-comment for performing the differential fits
  /*
  /// pt and rapidity binning (for differential measurements)
  const double ptBinArr  [] = {0,1,2,3,4,5,6,7,8,9,10,12,14,17,20,30};
  const double rapBinArr [] = {-2.0, 0.0, 2.0};
  const int npbin = sizeof(ptBinArr) / sizeof(double);
  const int nybin = sizeof(rapBinArr) / sizeof(double);

  /// perform the pt/eta differential fitting 
  xsectDif(*ws,1, ptBinArr,npbin, rapBinArr,nybin);
  */

  // close output root file  
  file.Close();
  
  ///oniafitterMain(finput, foutput, dfigs, ptBinArr,npbin, rapBinArr,nybin);

}
