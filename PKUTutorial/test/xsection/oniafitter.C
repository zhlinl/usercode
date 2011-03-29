#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooThresholdCategory.h"
#include "RooBinningCategory.h"
#include "RooSuperCategory.h"
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooSimPdfBuilder.h"
#include "RooCatType.h"
#include "Roo1DTable.h"
#include "RooCmdArg.h"
#include "RooChi2Var.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TGraphAsymmErrors.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TList.h"
#include "TRegexp.h"
#include "TLatex.h"
#include "TMath.h"

using namespace RooFit;

void readData(RooWorkspace&, TString);
void buildPdf(RooWorkspace&);

void fitYield(RooWorkspace&, TString, RooAbsData* d);
void plotMass(TString, RooWorkspace&, RooAbsData* d=0, RooAbsPdf* p=0, TString title="");

void xsectTot(RooWorkspace&);
void xsectDif(RooWorkspace&, bool, double*, int, double*, int);

void getInitialParValues(RooWorkspace&);
void setDefaultFitParValues(RooWorkspace&,RooAbsData* d=0);
double getWeightAv(RooWorkspace&, RooAbsData*);
double fitChisq(RooWorkspace& w, RooDataSet* data, RooAbsPdf* pdf);
//double fitProb(RooWorkspace& w, RooDataSet* data, RooAbsPdf* pdf);
double fitProb(RooWorkspace& w, RooDataSet* data, RooAbsPdf* pdf, double& chi2, double& ndf);

void   oniafitterMain  (const TString, const TString, const TString, double*,int, double*,int);
double oniafitterSimple(const TString, const TString, const TString, double*,int, double*,int, int);

RooArgSet * doSeqFit(RooWorkspace& oniaWS, int binning, int rapBinning,
		     TString fitParams, RooDataSet * data = 0, 
		     ostream& fout = std::cout, bool doPlots = false);
void plotDataPdf(RooDataSet * data, RooAbsPdf * const pdf, int resonance,
		 RooArgSet const& projVars, bool frameless = false);
void addPtBinningToData(RooDataSet * data, int binning);


static vector<double> initParVal_;
static vector<string> initParName_;

static bool linearbg_(false);
static double mmin_(8.0);
static double mmax_(14.0);
static int const nMassBins_ = 120;
static TString figs_("figs/"); /* store location for produced file output */
static TString res_("results"); /* root folder in file for storing results */
static TString plt_("plots"); /* root folder in file for storing plots */
static TString dirname_(""); /* tree location in input file */
static TString treeName("probe_tree");//("UpsilonTree");
static TString ext_(".pdf"); /* save plots format */
static double lumi_ = 3050.; /*sample lumi in nb-1*/
static bool normalize_ = true; /*normalize results with lumi and bin width*/

static const bool isMC        = false;//use thruth info
static const bool fitfraction = false; /*modified pdf, to fit for fractions */
static const bool restorepars = false;
static const bool printPars   = true;
//const bool absoluteRap = true;
const bool dofit_ = true;

//pt diff, 1 y bin
static const double ptdiff_1ybin_rapBinArr [] = {-2.0, 0.0, 2.0};
static const double ptdiff_1ybin_ptBinArr3S[] = {0, 3, 6, 9, 14, 20, 30};//3S
static const double ptdiff_1ybin_ptBinArr2S[] = {0, 2, 4, 6, 9,12,16,20,30};//2S
static const double ptdiff_1ybin_ptBinArr1S[] = {0,1,2,3,4,5,6,7,8,9,10,12,14,17,20,30};//1S

//pt diff, 2 y bins
static const double ptdiff_2ybin_rapBinArr [] = {-2.0, -1, 0.0, 1, 2.0};
static const double ptdiff_2ybin_ptBinArr3S[] = {0, 7, 12, 30};//3S
static const double ptdiff_2ybin_ptBinArr2S[] = {0, 3, 7, 11,15,30};//2S
static const double ptdiff_2ybin_ptBinArr1S[] = {0, 2, 5, 8, 11, 15, 30};//1S

//y diff, 1 pt bin
static const double ydiff_rapBinArr [] = {-2.0,-1.6,-1.2,-0.8,-0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0};
static const double ydiff_ptBinArr3S[] = {0, 30};
static const double ydiff_ptBinArr2S[] = {0, 30};
static const double ydiff_ptBinArr1S[] = {0, 30};


static TString PEAK_("1");

void oniafitter(
		const TString finput  = "upsilonYieldWeighted_nominal.root", 
		const TString foutput = "xsection.root", 
		const TString dfigs = "figs/", 
		const TString peak = "1",
		const TString linearbg = "false",
		const TString anamode = "mode0") {
  
  PEAK_ = peak;

  if     (linearbg=="true")  {linearbg_ = true;  mmin_=8.0; mmax_=12.0;}
  else if(linearbg=="false") {linearbg_ = false; mmin_=8.0; mmax_=14.0;}
  else {cout << "problem specifying background configuration\n\tExiting.\n"; return;}
  
  int npbin = 0;
  int nybin = 0;
  if        (anamode=="mode0") {
    if      (peak=="1") npbin = sizeof(ptdiff_1ybin_ptBinArr1S) / sizeof(double);
    else if (peak=="2") npbin = sizeof(ptdiff_1ybin_ptBinArr2S) / sizeof(double);
    else if (peak=="3") npbin = sizeof(ptdiff_1ybin_ptBinArr3S) / sizeof(double);
    nybin = sizeof(ptdiff_1ybin_rapBinArr) / sizeof(double);
  } else if (anamode=="mode1") {
    if      (peak=="1") npbin = sizeof(ptdiff_2ybin_ptBinArr1S) / sizeof(double);
    else if (peak=="2") npbin = sizeof(ptdiff_2ybin_ptBinArr2S) / sizeof(double);
    else if (peak=="3") npbin = sizeof(ptdiff_2ybin_ptBinArr3S) / sizeof(double);
    nybin = sizeof(ptdiff_2ybin_rapBinArr) / sizeof(double);
  } else if (anamode=="mode2") {
    if      (peak=="1") npbin = sizeof(ydiff_ptBinArr1S) / sizeof(double);
    else if (peak=="2") npbin = sizeof(ydiff_ptBinArr2S) / sizeof(double);
    else if (peak=="3") npbin = sizeof(ydiff_ptBinArr3S) / sizeof(double);
    nybin = sizeof(ydiff_rapBinArr) / sizeof(double);
  } else if (anamode=="mode3") {
    //all peaks fit with same binning, for computing ratios
    npbin = sizeof(ptdiff_1ybin_ptBinArr3S) / sizeof(double);
    nybin = sizeof(ptdiff_1ybin_rapBinArr)  / sizeof(double);
  }
  
  const int nptb = npbin;
  const int nrapb= nybin;
  
  //cout << "npt: " << nptb << " nrap: " << nybin << endl; return;

  double UpsPtBinEdges [nptb];
  for(int i=0; i<nptb; i++) {
    if        (anamode=="mode0") {
      if      (peak=="1") UpsPtBinEdges[i] = ptdiff_1ybin_ptBinArr1S[i];
      else if (peak=="2") UpsPtBinEdges[i] = ptdiff_1ybin_ptBinArr2S[i];
      else if (peak=="3") UpsPtBinEdges[i] = ptdiff_1ybin_ptBinArr3S[i];
    } else if (anamode=="mode1") {
      if      (peak=="1") UpsPtBinEdges[i] = ptdiff_2ybin_ptBinArr1S[i];
      else if (peak=="2") UpsPtBinEdges[i] = ptdiff_2ybin_ptBinArr2S[i];
      else if (peak=="3") UpsPtBinEdges[i] = ptdiff_2ybin_ptBinArr3S[i];
    } else if (anamode=="mode2") {
      if      (peak=="1") UpsPtBinEdges[i] = ydiff_ptBinArr1S[i];
      else if (peak=="2") UpsPtBinEdges[i] = ydiff_ptBinArr2S[i];
      else if (peak=="3") UpsPtBinEdges[i] = ydiff_ptBinArr3S[i];
    } else if (anamode=="mode3") {
      //all peaks fit with same binning, for computing ratios
      UpsPtBinEdges[i] = ptdiff_1ybin_ptBinArr3S[i];
    }
  }
  
  double UpsRapBinEdges[nrapb];
  for(int i=0; i<nrapb; i++) {
    if        (anamode=="mode0") {
      UpsRapBinEdges[i] = ptdiff_1ybin_rapBinArr[i];
    } else if (anamode=="mode1") {
      UpsRapBinEdges[i] = ptdiff_2ybin_rapBinArr[i];
    } else if (anamode=="mode2") {
      UpsRapBinEdges[i] = ydiff_rapBinArr[i];
    } else if (anamode=="mode3") {
      UpsRapBinEdges[i] = ptdiff_1ybin_rapBinArr[i];
    }
  }
  
  ///full analysis fitting 
  oniafitterMain(finput, foutput, dfigs, UpsPtBinEdges,nptb, UpsRapBinEdges,nrapb);

  ///partial/costumizable fitting
  //oniafitterSimple(finput, foutput, dfigs,UpsPtBinEdges,nptb,UpsRapBinEdges,nrapb,peak);
}


void oniafitterMain(
		     const TString finput,
		     const TString foutput,
		     const TString dfigs,
		     double *UpsPtBinEdges, int nptb,
		     double *UpsRapBinEdges, int nrapb
		     ) { 
  
  figs_ = dfigs;
  cout << "oniafitter processing"
       << "\n\tinput:  \t" << finput
       << "\n\toutput: \t" << foutput
       << "\n\tresults:\t" << figs_
       << endl;

  /// roofit workspace to manage fitting info
  RooWorkspace* ws = new RooWorkspace("ws","upsilon mass");

  /// create fitting model
  buildPdf(*ws);

  /// read the data
  readData(*ws,finput);
  //plotMass("data",*ws,ws->data("data"));


  /*
  RooDataSet * t_data = 
    new RooDataSet("t_data", "t_data", *(ws->data("data")->get()),
		   Import( *(dynamic_cast<RooDataSet *>(ws->data("data"))) ));
  addPtBinningToData(t_data, 1);
  RooCategory ptbin(*(dynamic_cast<RooCategory *>(t_data->get()->find("ptbin"))));
  Roo1DTable * tab = t_data->table(ptbin);
  tab->Print("v");
  delete tab;
  delete t_data;
  */
  /// cache the initial parameter
  getInitialParValues(*ws);

  /// open output root file
  //TFile file(foutput,"recreate");
  TFile file(foutput,"update");
  gDirectory->mkdir(res_);
  gDirectory->mkdir(plt_);
  gDirectory->Cd(res_);

  /* PLEASE UNCOMMENT BELOW ONLY THOSE FUNCTION CALLS
     YOU WANT TO BE EXECUTED EACH TIME THE MACRO IS RUN */

  /// fit raw yield
  if(dofit_) {
    fitYield(*ws,"total",ws->data("data"));
    plotMass("mass_raw",*ws);
  }

  /*
  t_data = 
    new RooDataSet("t_data", "t_data", *(ws->data("data")->get()),
		   Import( *(dynamic_cast<RooDataSet *>(ws->data("data"))) ));
  addPtBinningToData(t_data, 1);
  tab = t_data->table(ptbin);
  tab->Print("v");
  delete tab;
  delete t_data;
  */

  //return;
  ///...
  //((RooDataSet*)ws->data("data"))->setWeightVar(*ws->var("weight"));
  RooDataSet * data_w = 
    new RooDataSet("data_w", "data_w", *(ws->data("data")->get()),
		   Import(*((RooDataSet *)ws->data("data"))),
		   WeightVar("weight"));
  // fitYield(*ws,"total",ws->data("data"));
  //fitYield(*ws,"total",data_w);

  //RooDataSet* data_w  = ((RooDataSet*)ws->data("data"))->Clone();
  // RooDataSet* data_w((RooDataSet*)ws->data("data"));
  // data_w->setWeightVar(*ws->var("weight"));
  
  if(dofit_) {
    fitYield(*ws,"total",data_w);
    plotMass("mass_wei",*ws, data_w, ws->pdf("pdf"));
  }
  
  /// cross section vs pT
  //int nptb  = sizeof(UpsPtBinEdges)  / sizeof(double);
  //int nrapb = sizeof(UpsRapBinEdges) / sizeof(double);  

  /// total cross section
  
  if(fitfraction) {
    xsectDif(*ws,1, UpsPtBinEdges, nptb, UpsRapBinEdges, nrapb);
    return;
  }
  
  /*
  t_data = 
    new RooDataSet("t_data", "t_data", *(ws->data("data")->get()),
		   Import( *(dynamic_cast<RooDataSet *>(ws->data("data"))) ));
  addPtBinningToData(t_data, 1);
  tab = t_data->table(ptbin);
  tab->Print("v");
  delete tab;
  delete t_data;
  */
  //differential cross section and fits
  xsectDif(*ws,dofit_, UpsPtBinEdges, nptb, UpsRapBinEdges, nrapb);
  //xsectDif(*ws,0, UpsPtBinEdges, nptb, UpsRapBinEdges, nrapb);

  //total cross section
  xsectTot(*ws);
  
  //return;

  file.Close();
}


/* total cross section
   retrieve and normalize global fit results
*/
void xsectTot(RooWorkspace& w) {
  cout << "computing total cross section...\n" << std::flush;
  gDirectory->Cd("../"+res_);
  //cout <<   gDirectory->pwd() << endl;
  //cout <<   gDirectory->ls() << endl;

  ///raw yield
  RooFitResult* fitres_r  = (RooFitResult *) gROOT->FindObject("fit_result_raw_total");
  ///corrected yield
  RooFitResult* fitres_w  = (RooFitResult *) gROOT->FindObject("fit_result_corrected_total");

  //  cout << " average weight:" << ((RooDataSet*)w.data("data"))->mean(*w.var("weight")) << "\n";


  const int npeak = 3;
  double nsigVal_r[npeak], nsigErr_r[npeak];
  double nsigVal_w[npeak], nsigErr_w[npeak], nsigErrHi_w[npeak], nsigErrLo_w[npeak];
  double xsec[npeak], xsecE[npeak], xsecEhi[npeak], xsecElo[npeak];
  for(int j=0; j<npeak; j++) {
    TString yieldn = TString::Format("nsig%d",j+1);
    nsigVal_r[j]        = ((RooRealVar*)fitres_r->floatParsFinal().find(yieldn))->getVal();
    nsigErr_r[j]        = ((RooRealVar*)fitres_r->floatParsFinal().find(yieldn))->getError();
    nsigVal_w[j]        = ((RooRealVar*)fitres_w->floatParsFinal().find(yieldn))->getVal();
    nsigErr_w[j]        = ((RooRealVar*)fitres_w->floatParsFinal().find(yieldn))->getError();
    nsigErrHi_w[j]      = ((RooRealVar*)fitres_w->floatParsFinal().find(yieldn))->getErrorHi();
    nsigErrLo_w[j]      = ((RooRealVar*)fitres_w->floatParsFinal().find(yieldn))->getErrorLo();
  }

  //debug
  if(0) {
    double _yy, _yy_e, _yy_el, _yy_eh;
    TString yieldn = TString::Format("nsig%d",1);
    _yy    = ((RooRealVar*)fitres_w->floatParsFinal().find(yieldn))->getVal();
    _yy_e  = ((RooRealVar*)fitres_w->floatParsFinal().find(yieldn))->getError();
    _yy_eh = ((RooRealVar*)fitres_w->floatParsFinal().find(yieldn))->getErrorHi();
    _yy_el = ((RooRealVar*)fitres_w->floatParsFinal().find(yieldn))->getErrorLo();
    printf("check unc: yield:%5.3f error:%5.3f  err_hi:%5.3f   err_lo:%5.3f \n",
	   _yy, _yy_e, _yy_eh, _yy_el);
  }

  /// check & print
  //this number varies whether the dataset is weighted or not(!)
  double mean_weight = ((RooDataSet*)w.data("data"))->mean(*w.var("weight"));
  //cout << "xxx " <<  mean_weight << "\n";
  //return;

  //this fails for weighted dataset, as weight stops being a variable, by becoming a weight 
  //double mean2 = (w.data("data"))->createHistogram("weight",100) ->GetMean();
  //cout << "xxx: "<< getWeightAv(w, w.data("data")) << " mean2:" << mean2 << "\n";

  cout << "raw  yield " << nsigVal_r[0] << "+-" << nsigErr_r[0] << "\n";
  cout << "cross check raw * weight(" <<  mean_weight << ") = " 
       << nsigVal_r[0]*mean_weight << "+-" << nsigErr_r[0]*mean_weight << "\n";
  cout << "weighted yield " << nsigVal_w[0] << "+-" << nsigErr_w[0] << "\n";
  cout << "xsection:" << nsigVal_w[0]/lumi_ << "+-" << nsigErr_w[0]/lumi_ << " (lumi: "<< lumi_ <<")\n";
  //cout << "cross section:" << xsection << "+-" << xsection_err << " nb\n";


  for(int j=0; j<npeak; j++) {
    /// normalize
    if(normalize_) {
      xsec [j] = nsigVal_w[j] / lumi_;
      xsecE[j] = nsigErr_w[j] / lumi_;
      xsecEhi[j] = fabs(nsigErrHi_w[j]) / lumi_;
      xsecElo[j] = fabs(nsigErrLo_w[j]) / lumi_;
    }
    /// save value as root object
    TGraphAsymmErrors xst;
    xst.SetPoint(0,0, xsec[j]);
    xst.SetPointError(0,0.1,0.1, xsecElo[j], xsecEhi[j]);
    gDirectory->Cd("../"+plt_);
    xst.SetName(TString::Format("xsection_ups%dS_total",j+1));
    xst.Write();
  }
}

/* differential cross section, sigma vs pT
-split the dataset in pT regions
-fit the sub-samples, extract the yield, and normalize   
*/
void xsectDif(RooWorkspace& w, bool dofit, double *UpsPtBinEdges, int nptb, double *UpsRapBinEdges, int nrapb) {

  RooAbsPdf*  pdf   = w.pdf("pdf");
  RooDataSet* data  = 
    new RooDataSet("data", "data", *(w.data("data")->get()),
		   Import( *(dynamic_cast<RooDataSet *>(w.data("data"))) ));
//dynamic_cast<RooDataSet *>(w.data("data"));
  RooRealVar* upsPt = w.var("upsPt");
  //RooRealVar * upsPt = dynamic_cast<RooRealVar *>(data->get()->find("upsPt"));
  RooRealVar* upsRapidity = w.var("upsRapidity");

  /* 
  RooDataSet * t_data = 
     new RooDataSet("t_data", "t_data", *(data->get()),
		    Import( *data ));
  addPtBinningToData(t_data, 1);
  RooCategory ptbin(*(dynamic_cast<RooCategory *>(t_data->get()->find("ptbin"))));
  Roo1DTable * tab2 = t_data->table(ptbin);
  tab2->Print("v");
  delete tab2;
  delete t_data;
  
  addPtBinningToData(data, 1);
  tab2 = data->table(ptbin);
  tab2->Print("v");
  delete tab2;
  */
 //const int npte = sizeof(UpsPtBinEdges)/sizeof(double);
  const int nptbins = nptb - 1;
  //const int nrapbins = sizeof(UpsRapBinEdges)/sizeof(double)-1;
  const int nrapbins = nrapb - 1;

  bool absoluteRap = (nrapbins==1)?false:true;
  const int nrapregion = absoluteRap?nrapbins/2:nrapbins;

  cout << "nrapregion:" << nrapregion << " abs-y?" << absoluteRap << " nrapbins:" << nrapbins << " nrapbin:" << nrapb << "\n";


  double UpsPtBinCenter[nrapregion][nptbins], UpsPtBinEdgeLo[nrapregion][nptbins], UpsPtBinEdgeHi[nrapregion][nptbins];
  double WeightAv[nrapregion][nptbins], WeightAvErr[nrapregion][nptbins];

  RooThresholdCategory ptRegion("ptRegion","region of pt",*upsPt);
  for(int i=0; i<nptbins; i++) {
    TString reg = TString::Format("PtBin%d",i+1);    
    ptRegion.addThreshold(UpsPtBinEdges[i+1],reg) ;    
  }
  data->addColumn(ptRegion);

  Roo1DTable * tab = data->table(ptRegion);
  tab->Print("v");
  delete tab;

  RooThresholdCategory rapRegion("rapRegion","region of rap",*upsRapidity) ;

  for(int i=0; i<nrapbins; i++) {
    int jj = (i<nrapregion)?nrapregion-i-1:i-nrapregion;
    jj = absoluteRap?jj:i;
    TString reg = TString::Format("RapBin%d",jj+1);    
    printf("edge:%d  name:%s\n",i,reg.Data());
    rapRegion.addThreshold(UpsRapBinEdges[i+1],reg) ;    
  }
  data->addColumn(rapRegion);

  tab = data->table(rapRegion);
  tab->Print("v");
  delete tab;

  RooThresholdCategory mRegion("mRegion","mass region",*w.var("invariantMass"));
  mRegion.addThreshold( 8.6,  "bgl") ;
  mRegion.addThreshold( 9.26, "vvv") ;
  mRegion.addThreshold( 9.66, "sig1s") ;
  mRegion.addThreshold( 9.873,"vvv") ;
  mRegion.addThreshold(10.173,"sig2s");
  mRegion.addThreshold(10.205,"vvv") ;
  mRegion.addThreshold(10.505,"sig3s");
  mRegion.addThreshold(11.,    "vvv");
  mRegion.addThreshold(13.,    "bgh");

  data->addColumn(mRegion);

  data->Print();
  tab = data->table(RooArgList(ptRegion,rapRegion));
  tab->Print("v");
  delete tab;

  /// check
  int speak=-1;
  if     (PEAK_=="1") speak = 0;
  else if(PEAK_=="2") speak = 1; 
  else if(PEAK_=="3") speak = 2; 

  double countBinSumEntries[nrapregion][nptbins];
  double countBinNumEntries[nrapregion][nptbins];
  //double averageWeight[nrapregion][nptbins];
  double chisqred[nrapregion][nptbins];
  double fitprob[nrapregion][nptbins];

  double mychsq[nrapregion][nptbins];   double myndof[nrapregion][nptbins];


  double meanHistPt[nrapregion][nptbins], rmsHistPt[nrapregion][nptbins];
  double meanHistWei[nrapregion][nptbins], rmsHistWei[nrapregion][nptbins];
  /// process each pt subsample
  RooDataSet* dataPt, *dataPtSignalRaw, *dataPtRaw;
  for(int j=0; j<nrapregion; j++) {
    for(int i=0; i<nptbins; i++) {
      cout << "processing subsample pt" << i << " rap:" << j << "\n" << std::flush;
      // dataPt = (RooDataSet*)  data->reduce(TString::Format("rapRegion==rapRegion::RapBin%d",j+1))->reduce(TString::Format("ptRegion==ptRegion::PtBin%d",i+1));
      //dataPt = (RooDataSet*) data->reduce(TString::Format("ptRegion==ptRegion::PtBin%d && rapRegion::RapBin%d",i+1,j+1));
      //averageWeight[j][i]=(dataPt->createHistogram("weight",100))->GetMean();

      //dataPt->setWeightVar(*w.var("weight"));
      TString ptcut(TString::Format("(rapRegion==rapRegion::RapBin%d)", j+1));
      ptcut += TString::Format("&&(ptRegion==ptRegion::PtBin%d)", i+1);

      dataPt = new RooDataSet("dataPt", "dataPt", *(data->get()),
			      Import(*data), WeightVar("weight"),
			      Cut(ptcut));

      std::cout << ptcut << std::endl;
      // dataPt->Print();
      if (dataPt->numEntries() < 1) return;

      //sideband subtracted weighted yield
      // signal: 9.46-6*0.09+3*0.09; background: 8.8.54, 11-11.54; fraction: S/B = (6+3)*0.09 / 2*0.5 = 0.81 
      // countBinSumEntries[j][i] = dataPt->reduce(Cut("invariantMass>8.92 && invariantMass<9.73"))->sumEntries() 
      //	- 0.81 * ( dataPt->reduce(Cut("invariantMass>8.0 && invariantMass<8.5"))->sumEntries() +
      //	   dataPt->reduce(Cut("invariantMass>11.0 && invariantMass<11.5"))->sumEntries() ) ;
      
      dataPtRaw = new RooDataSet("dataPt", "dataPt", *(data->get()),Import(*data), Cut(ptcut));

      TString sigcut =  TString::Format("mRegion==mRegion::sig%ss",  PEAK_.Data());
      dataPtSignalRaw = new RooDataSet("dataPt", "dataPt", *(data->get()),Import(*data), Cut(ptcut + "&&" + sigcut));
      
      RooRealVar* meanpt  = dataPtSignalRaw->meanVar(*upsPt);
      RooRealVar* meanwei = dataPtSignalRaw->meanVar(*w.var("weight"));

      countBinSumEntries[j][i] = dataPtSignalRaw->sumEntries();
      countBinNumEntries[j][i] = dataPtSignalRaw->numEntries();

      gROOT->SetStyle("Plain");
      gStyle->SetTitleBorderSize(0.);
      TH1F* hpt = (TH1F*) dataPtSignalRaw->createHistogram("upsPt",1000);
      hpt->SetTitle(TString::Format("#Upsilon(%dS) p_{T} distribution",speak+1));
      hpt->GetXaxis()->SetRangeUser(UpsPtBinEdges[i]-1,UpsPtBinEdges[i+1]+1);
      meanHistPt[j][i] = hpt->GetMean();
      rmsHistPt [j][i] = hpt->GetRMS();
      TCanvas cp;hpt->Draw(); cp.SaveAs(TString::Format("%spt_%ds_y%d_pt%d%s",figs_.Data(),speak+1,j,i,ext_.Data()));
      delete hpt;
      TH1F* hwei;
      hwei = (TH1F*) dataPtSignalRaw->createHistogram("weight",500);
      hwei->SetTitle(TString::Format("#Upsilon(%dS) weight distribution",speak+1));
      hwei->GetXaxis()->SetRangeUser(0,10);
      meanHistWei[j][i] = hwei->GetMean();
      rmsHistWei [j][i] = hwei->GetRMS();
      hwei->Draw(); cp.SaveAs(TString::Format("%swei_%ds_y%d_pt%d%s",figs_.Data(),speak+1,j,i,ext_.Data()));
      delete hwei;


      if(0)
	cout << "TESTE  " 
	     << " bin mean  pt-all:" << dataPt->mean(*w.var("upsPt"))
	  //<< " sig mean:"    << dataPtSignalRaw->mean(*w.var("upsPt"))
	     << " pt_sig:" << meanpt->getVal()
	     << "+/-" << meanpt->getError() 
	     << " pt_hist:" << meanHistPt[j][i]
	     << " rms:" << rmsHistPt[j][i]
	     << "  weight-sig:" << meanwei->getVal()
	     << "+/-" << meanwei->getError() 
	     << endl;

      
      //WeightAv   [j][i]  = meanwei->getVal();
      //WeightAvErr[j][i]  = meanwei->getError();

      WeightAv   [j][i]  = meanHistWei[j][i];
      WeightAvErr[j][i]  = rmsHistWei [j][i];

      //UpsPtBinCenter[j][i] = meanpt->getVal();
      //UpsPtBinEdgeLo[j][i] = fabs(meanpt->getErrorLo());
      //UpsPtBinEdgeHi[j][i] = fabs(meanpt->getErrorHi());

      UpsPtBinCenter[j][i] = meanHistPt[j][i];
      UpsPtBinEdgeLo[j][i] = rmsHistPt [j][i];
      UpsPtBinEdgeHi[j][i] = rmsHistPt [j][i];

      
      //UpsPtBinCenter[j][i] =  dataPt->mean(*w.var("upsPt")); 
      //UpsPtBinEdgeLo[j][i] =  UpsPtBinCenter[j][i] - UpsPtBinEdges[i];
      //UpsPtBinEdgeHi[j][i] = -UpsPtBinCenter[j][i] + UpsPtBinEdges[i+1];
      //Cout << "NTOTAL sample " << i << " = " << dataPt->sumEntries() << "\n";
      //if(i!=0) continue; /* re-fit sub-sample */
      //dofit=0;
      
      TString bintitle(""); 
      if(!absoluteRap)    
	bintitle.Append(TString::Format("%3.1f<y^{#Upsilon}<%3.1f", UpsRapBinEdges[j], UpsRapBinEdges[j+1]));
      else
	bintitle.Append(TString::Format("%3.1f<|y^{#Upsilon}|<%3.1f", UpsRapBinEdges[nrapbins/2+j], UpsRapBinEdges[nrapbins/2+j+1]));
      
      bintitle.Append(TString::Format(",  %3.1f<p_{T}^{#Upsilon}<%3.1f", UpsPtBinEdges[i],UpsPtBinEdges[i+1]));
      
      // double aa,bb;

      if(dofit) {
	fitYield(w, TString::Format("rap%d_pt%d",j,i), dataPtRaw);
	plotMass(TString::Format("massfit_raw_rap%d_pt%d",j,i),w, dataPtRaw,pdf,bintitle);
	fitYield(w, TString::Format("rap%d_pt%d",j,i), dataPt);
	plotMass(TString::Format("massfit_rap%d_pt%d",j,i),w, dataPt,pdf,bintitle);  
	double chi2red = fitChisq(w,(RooDataSet*)dataPt,pdf);
	//RooPlot* frame = w.var("invariantMass")->frame(Bins(nMassBins_));
	//int nfloatpars = pdf->getParameters(dataPt)->selectByAttrib("Constant",kFALSE)->getSize(); 
	//double chi2red = frame->chiSquare(nfloatpars); 
	chisqred[j][i] = chi2red;
	fitprob [j][i] = fitProb(w,(RooDataSet*)dataPt,pdf,
				 mychsq[j][i], myndof[j][i]);
	
       }
    }
  }

 
  if(fitfraction) {
    
    cout << "extraction 2S+3S/1S fraction...\n";

    cout << "aaa 1 \n" << flush;

    //first: total ratio
    gDirectory->Cd("../"+res_);    
    RooFitResult* fitres_w  = (RooFitResult *) gROOT->FindObject("fit_result_corrected_total");
    //cout << fitres_w << "\n";

    double _f23os[1]={0.}; double _f23os_e[1]={0.}; double _f23o1[1]={0.}; double _f23o1_e[1]={0.};

    _f23os  [0] = ((RooRealVar*)fitres_w->floatParsFinal().find("f23os"))->getVal();
    _f23os_e[0] = ((RooRealVar*)fitres_w->floatParsFinal().find("f23os"))->getError();
    _f23o1  [0] = (_f23os[0]==1)? 0: _f23os[0]/(1.-_f23os[0]);
    _f23o1_e[0] = (_f23os[0]==1)? 0: _f23os_e[0]/(1.-_f23os[0])/(1.-_f23os[0]);

    cout << " total fraction " << _f23o1[0] << " +- " << _f23o1_e[0] << "\n";
    double xv [1]={5}; 
    double xvl[1]={5}; 
    double xvh[1]={15};
    //double xvh[1]={UpsPtBinEdges[sizeof(UpsPtBinEdges) / sizeof(double)]}; 
    TGraphAsymmErrors frxst(1,xv,_f23o1,xvl,xvh,_f23o1_e,_f23o1_e);
		   //frxst.SetPoint(0,5, f23o1);
    //frxst.SetPointError(0,UpsPtBinEdges[0],UpsPtBinEdges[sizeof(UpsPtBinEdges) / sizeof(double)], f23o1_e, f23o1_e);
    gDirectory->Cd("../"+plt_);
    frxst.SetName(TString::Format("fraction_total"));
    frxst.Write();

    TCanvas bla; bla.cd(); 
    //frxst.Draw("ap"); 
    frxst.SetFillColor(kGray);
    frxst.Draw("PE2");
    frxst.Draw("Psame");

    bla.SaveAs("bla.gif");   

    cout << "aaa 2 \n" << flush;

    double fraction23os  [nrapregion][nptbins];
    double fraction23os_e[nrapregion][nptbins];
    double fraction23o1  [nrapregion][nptbins];
    double fraction23o1_e[nrapregion][nptbins];
    gDirectory->Cd("../"+res_);
    RooFitResult* fitresPt;

    cout << "aaa 3 \n" << flush;

    for(int j=0; j<nrapregion; j++) {
      for(int i=0; i<nptbins; i++) {
	fitresPt = (RooFitResult *) gROOT->FindObject(TString::Format("fit_result_corrected_rap%d_pt%d",j,i));
	double f23os(0), f23os_e(0), f23o1(0), f23o1_e(0);
    	//f23os=0; f23os_e=0; f23o1=0; f23o1_e=0;
	f23os   = ((RooRealVar*)fitresPt->floatParsFinal().find("f23os"))->getVal();
	f23os_e = ((RooRealVar*)fitresPt->floatParsFinal().find("f23os"))->getError();
	f23o1   = (f23os==1)? 0: f23os/(1.-f23os);
	f23o1_e = (f23os==1)? 0: f23os_e/(1.-f23os)/(1.-f23os);
	fraction23os  [j][i] = f23os;
	fraction23os_e[j][i] = f23os_e;
	fraction23o1  [j][i] = f23o1;
	fraction23o1_e[j][i] = f23o1_e;
      }
    }
    
    cout << "aaa 4 \n" << flush;

    TGraphAsymmErrors *fracGr;
    
    for(int Y=0; Y<nrapregion; Y++) {
      double yield[nptbins], err[nptbins];
      double yieldb[nptbins], errb[nptbins];
      double UpsPtBinCenterY[nptbins], UpsPtBinEdgeLoY[nptbins], UpsPtBinEdgeHiY[nptbins]; 
      for(int i=0; i<nptbins; i++) {
	yield[i] = fraction23o1  [Y][i];
	err  [i] = fraction23o1_e[Y][i];
	yieldb[i] = fraction23os  [Y][i];
	errb  [i] = fraction23os_e[Y][i];
	UpsPtBinCenterY[i]=UpsPtBinCenter[Y][i];
	UpsPtBinEdgeLoY[i]=UpsPtBinEdgeLo[Y][i];
	UpsPtBinEdgeHiY[i]=UpsPtBinEdgeHi[Y][i]; 
      }

      cout << "aaa 5 \n" << flush;
      fracGr = new TGraphAsymmErrors(nptbins,   UpsPtBinCenterY, yield,
				     UpsPtBinEdgeLoY, UpsPtBinEdgeHiY, 
				     err, err);
  
      gDirectory->Cd("../"+plt_);
      TString bla = TString::Format("cvs_rap_%d",Y);
      TCanvas c(bla, bla); c.cd();
      double x[2] = {0, 18.5};
      double y[2] = {0.,1.};
      TGraph frame(2,x,y);
      TString ytitle = 
	TString::Format("#frac{#sigma(#Upsilon(2S))+#sigma(#Upsilon(3S))}{#sigma(#Upsilon(1S))}");
      //TString::Format("#frac{d#sigma}{dp_{T}} . BR(#Upsilon(%dS)#rightarrow#mu#mu) [nb/GeV]",j+1);
      frame.SetTitle( "" );
      frame.GetXaxis()->SetTitle("p_{T} (#mu#mu) [GeV]");
      frame.GetYaxis()->SetTitle(ytitle);
      frame.Draw("AP");
      frxst.SetFillColor(kGray);
      frxst.Draw("PE2");
      //frxst.Draw("Psame");
      //frxst.Draw("psame");
      //gPad->SetLogy();
      fracGr->Draw("P");
      fracGr->SetName(TString::Format("fraction_rap%d_pt",Y));
      fracGr->Write();
      c.SaveAs(TString::Format("%sfraction_rap%d%s",figs_.Data(),Y,ext_.Data()));
      
      cout << "aaa 6 \n" << flush;

      TGraphAsymmErrors* fracGrb = new TGraphAsymmErrors(nptbins,   UpsPtBinCenterY, yieldb,
							 UpsPtBinEdgeLoY, UpsPtBinEdgeHiY, 
							 errb, errb);
      
      TCanvas cb(bla+"bla",bla+"bla"); cb.cd();
      TGraph frameb(2,x,y);
      TString ytitle2 = 
	TString::Format("2S+3S signal fraction");
      //TString::Format("#frac{#sigma(#Upsilon(2S))+#sigma(#Upsilon(3S))}{#sigma(#Upsilon)}");
      //TString::Format("#frac{d#sigma}{dp_{T}} . BR(#Upsilon(%dS)#rightarrow#mu#mu) [nb/GeV]",j+1);
      frameb.SetTitle( "" );
      frameb.GetXaxis()->SetTitle("p_{T} (#mu#mu) [GeV]");
      frameb.GetYaxis()->SetTitle(ytitle2);
      frameb.Draw("AP");
      fracGrb->Draw("P");
      //gPad->SetLogy();
      fracGrb->SetName(TString::Format("fractionb_rap%d_pt",Y));
      fracGrb->Write();
      cb.SaveAs(TString::Format("%sfractionb_rap%d%s",figs_.Data(),Y,ext_.Data()));
      
    cout << "aaa 7 \n" << flush;

    }
  } ///end fit fraction


  /// store yield
  const int npeak = 3;
  double UpsYieldPt  [npeak][nrapregion][nptbins];
  double UpsYieldPt_e[npeak][nrapregion][nptbins];
  double UpsYieldPt_eh[npeak][nrapregion][nptbins];
  double UpsYieldPt_el[npeak][nrapregion][nptbins];

  double UpsYieldPtRaw  [nrapregion][nptbins];
  double UpsYieldPtRaw_e[nrapregion][nptbins];

  double UpsYieldPtTot(0.), UpsYieldPtTot_e(0.); 
  double UpsYieldPtRawTot(0.), UpsYieldPtRawTot_e(0.); 

  gDirectory->Cd("../"+res_);
  RooFitResult* fitresPt, *fitresPtRaw;
  TString yieldn(""); 
  for(int j=0; j<nrapregion; j++) {
    for(int i=0; i<nptbins; i++) {
      fitresPt = (RooFitResult *) gROOT->FindObject(TString::Format("fit_result_corrected_rap%d_pt%d",j,i));
      TString yieldn(""); 
      for(int k=0; k<npeak; k++) {
	yieldn = TString::Format("nsig%d",k+1);
	UpsYieldPt   [k][j][i] = ((RooRealVar*)fitresPt->floatParsFinal().find(yieldn))->getVal();
	UpsYieldPt_e [k][j][i] = ((RooRealVar*)fitresPt->floatParsFinal().find(yieldn))->getError();
	UpsYieldPt_eh[k][j][i] = ((RooRealVar*)fitresPt->floatParsFinal().find(yieldn))->getErrorHi();
	UpsYieldPt_el[k][j][i] = ((RooRealVar*)fitresPt->floatParsFinal().find(yieldn))->getErrorLo();
      }
      yieldn = TString::Format("nsig%s",PEAK_.Data());
      fitresPtRaw = (RooFitResult *) gROOT->FindObject(TString::Format("fit_result_raw_rap%d_pt%d",j,i));
      UpsYieldPtRaw  [j][i] = ((RooRealVar*)fitresPtRaw->floatParsFinal().find(yieldn))->getVal();
      UpsYieldPtRaw_e[j][i] = ((RooRealVar*)fitresPtRaw->floatParsFinal().find(yieldn))->getError();
      printf("TESTE  (rap:%d,pt:%d) 1S yields raw=%5.0f  wei=%5.0f\n", j,i, UpsYieldPtRaw[j][i], UpsYieldPt[0][j][i]);

    }
  }

  if(0) {
    for(int j=0; j<nrapregion; j++) {
      for(int i=0; i<nptbins; i++) {
	double wei = 0.;//averageWeight[j][i];
	double rr = countBinNumEntries[j][i]?countBinSumEntries[j][i]/countBinNumEntries[j][i]:0;
	rr = wei?rr/wei:0;
	double r = UpsYieldPt[0][j][i]/countBinSumEntries[j][i];
	printf(" bin eta%d pt%d count:%8.0f fit:%8.0f r:%5.3f\n",
	       j,i,countBinSumEntries[j][i], UpsYieldPt[0][j][i],r);
      }
    }
  }
  
  if(0) {//debug
    for(int j=0; j<nrapregion; j++) {
      for(int i=0; i<nptbins; i++) {
	for(int k=0; k<npeak; k++) {  
	  printf(" %dS yield=%5.0f err:%3.0f err_hi:%3.0f err_lo:%3.0f ", 
		 k+1, UpsYieldPt[k][j][i], UpsYieldPt_e[k][j][i],
		 UpsYieldPt_eh[k][j][i], UpsYieldPt_el[k][j][i]
		 );
	}
	cout << "\n";
      }
    }
  }

  for(int j=0; j<nrapregion; j++) {
    for(int i=0; i<nptbins; i++) {
      printf(" pt (%4.1f-%4.1f) mean=%6.2f", 
	     UpsPtBinEdges[i],UpsPtBinEdges[i+1], UpsPtBinCenter[j][i]); 
      for(int k=0; k<npeak; k++) {  
	printf("  %dS: yield=%5.0f+/-%3.0f  ", k+1, UpsYieldPt[k][j][i], UpsYieldPt_e[k][j][i]);
      }
      cout << "\n";
    }
  }
  
  
  double UpsRapRegEdges[nrapregion+1];
  for(int j=0; j<nrapregion+1; j++) {
    if(!absoluteRap)    
      UpsRapRegEdges[j] = UpsRapBinEdges[j];
    else
      UpsRapRegEdges[j] = UpsRapBinEdges[nrapbins/2+j];
  }

  //bool normalize=1;
  if(normalize_){ //only lumi for now, so we can add results per bin for checks
    /// normalize yield
    for(int j=0; j<nrapregion; j++) {
      for(int i=0; i<nptbins; i++) {
	// double binw =  UpsPtBinEdges[i+1]-UpsPtBinEdges[i];
	for(int k=0; k<npeak; k++) {
	  UpsYieldPt   [k][j][i] *= 1./lumi_;
	  UpsYieldPt_e [k][j][i] *= 1./lumi_;
	  UpsYieldPt_eh[k][j][i] *= 1./lumi_;
	  UpsYieldPt_el[k][j][i] *= 1./lumi_;
	}  
      }
    }
  }

  ///Latex TABLE
  char ss[1000];
  ofstream flatex;
  TString tabtexName[nrapregion];

  
  tabtexName[0] = TString::Format("%sfitvalues_ups%d.tex",figs_.Data(),speak+1);	    
  flatex.open (tabtexName[0]);
  
  flatex << "\\begin{tabular}{c ccc|ccc|ccc|ccc}\n\\hline\n";
  sprintf(ss,"\\multicolumn{4}{c}{$\\Upsilon(%dS) \\qquad \\pt ~((\\GeVc)$} & \\multicolumn{3}{|c|}{raw fit} &  \\multicolumn{3}{|c|}{weight, $w$} & \\multicolumn{3}{|c}{cross section}\\\\ \n",speak+1);
  flatex << ss;
  flatex << "\\multicolumn{1}{c|}";
  if(absoluteRap) flatex << "{$|\\rm{rapidity}|$}"; else flatex << "{$y$}";
  flatex << " & range &  mean &  rms & yield, $N$  & s/e & $\\chi^2$/ndf & mean & rms & $\\langle w\\rangle^{-1}$ & $\\sigma_{\\rm{av.}}\\sim N\\cdot\\langle w\\rangle $ & $\\sigma_{\\rm{fit}}$ & $\\Delta$ \\\\ \\hline \n";
  //flatex << " & range &  mean &  rms & yield, $N$  & s/e & $\\chi^2$/ndf & prob & mean & rms & $\\langle w\\rangle^{-1}$ & $\\sigma_{\\rm{av.}}\\sim N\\cdot\\langle w\\rangle $ & $\\sigma_{\\rm{fit}}$ & $\\Delta$ \\\\ \\hline \n";
  
  double xsec_av_tot(0.), xsec_av_tot_e(0.);

  for(int j=0; j<nrapregion; j++) {
    for(int i=0; i<nptbins; i++) {
      //if(i) continue;
      TString bintitle = absoluteRap ? "|y|" : "y";
      //      if(!absoluteRap)    
//	bintitle.Append(TString::Format("%3.1f<y<%3.1f", UpsRapBinEdges[j], UpsRapBinEdges[j+1]));
//      else
//	bintitle.Append(TString::Format("%3.1f<|y|<%3.1f", UpsRapBinEdges[nrapbins/2+j], UpsRapBinEdges[nrapbins/2+j+1]));
      printf(" Y(%dS) %s:(%3.1f-%3.1f) pt:(%4.1f-%4.1f) mean=%6.3f+/-%5.3f  weight: %6.3f+/-%5.3f", 
	     speak+1, 
	     bintitle.Data(), UpsRapRegEdges[j], UpsRapRegEdges[j+1], 
	     UpsPtBinEdges[i],UpsPtBinEdges[i+1], 
	     UpsPtBinCenter[j][i], UpsPtBinEdgeHi[j][i],
	     WeightAv[j][i], WeightAvErr[j][i]
	     ); 
           
      double xsec_av   = UpsYieldPtRaw[j][i] * WeightAv[j][i];
      double xsec_av_e = sqrt( pow(UpsYieldPtRaw_e[j][i]*WeightAv[j][i],2)  + pow(UpsYieldPtRaw_e[j][i]*WeightAvErr[j][i],2) );
      xsec_av   /= lumi_;//*(UpsPtBinEdges[i+1]-UpsPtBinEdges[i]);
      xsec_av_e /= lumi_;//*(UpsPtBinEdges[i+1]-UpsPtBinEdges[i]);
      printf("  raw_yield=%5.0f+/-%3.0f  xsec_av:%5.3f+/-%5.3f ", 
	     UpsYieldPtRaw[j][i], UpsYieldPtRaw_e[j][i],xsec_av,xsec_av_e);
      printf("  xsec_fit=%5.3f+/-%5.3f  s/e=%3.1f  chi2/ndf=%5.3f fitprob=%5.3f", 
	     UpsYieldPt[speak][j][i]/lumi_, UpsYieldPt_e[speak][j][i]/lumi_, 
	     UpsYieldPt[speak][j][i] / UpsYieldPt_e[speak][j][i],  
	     chisqred[j][i], fitprob[j][i] );
      cout << "\n";

      sprintf(ss,"%4.1f-%4.1f ",
	      UpsRapRegEdges[j], UpsRapRegEdges[j+1]);
      flatex << ss;
      sprintf(ss,"& %4.0f-%4.0f & %6.1f & %5.1f ",
	      UpsPtBinEdges[i],UpsPtBinEdges[i+1], 
	      UpsPtBinCenter[j][i], UpsPtBinEdgeHi[j][i]);
      flatex << ss;
      if(!dofit_)
	chisqred[j][i] = 0;
      //sprintf(ss,"& %6.0f $\\pm$ %6.0f & %3.0f & %4.0f/%3.0f & %4.2f ", 
      sprintf(ss,"& %6.0f $\\pm$ %6.0f & %3.0f & %4.1f", 
	      UpsYieldPtRaw[j][i],UpsYieldPtRaw_e[j][i], 
	      UpsYieldPt[speak][j][i] / UpsYieldPt_e[speak][j][i],  
	      chisqred[j][i]
	      //mychsq[j][i], myndof[j][i], fitprob[j][i]
	      ); 
      flatex << ss;
      
     sprintf(ss,"& %4.1f & %4.1f & %5.2f$\\pm$%5.2f   & %6.3f $\\pm$ %5.3f ",
	     WeightAv[j][i], WeightAvErr[j][i], 
	     1./WeightAv[j][i], WeightAvErr[j][i]/pow(WeightAv[j][i],2),
	     xsec_av,xsec_av_e);
     flatex << ss;
     sprintf(ss,"& %6.3f $\\pm$ %5.3f ",
	     UpsYieldPt[speak][j][i], UpsYieldPt_e[speak][j][i]);
     flatex << ss;
     sprintf(ss,"& %4.1f\\%% ",
	     (xsec_av-UpsYieldPt[speak][j][i])/UpsYieldPt[speak][j][i]*100);
     flatex << ss;
     sprintf(ss,"\\\\ ");
     flatex << ss;
     sprintf(ss,"%%Y(%dS) y:%d\n",speak+1,j);
     flatex << ss;

     xsec_av_tot   += xsec_av; 
     xsec_av_tot_e += pow(xsec_av_e,2);

     UpsYieldPtTot   += UpsYieldPt       [speak][j][i]; 
     UpsYieldPtTot_e += pow( UpsYieldPt_e[speak][j][i],2); 
     
     UpsYieldPtRawTot   += UpsYieldPtRaw       [j][i]; 
     UpsYieldPtRawTot_e += pow( UpsYieldPtRaw_e[j][i],2); 
     
     xsec_av_tot_e = sqrt(xsec_av_tot_e);
    }
  }
  
  UpsYieldPtTot_e    = sqrt(UpsYieldPtTot_e   );
  UpsYieldPtRawTot_e = sqrt(UpsYieldPtRawTot_e);

  //add sums to table
  //flatex << "\\begin{tabular}{cccc|ccc|cc|ccc}\n\\hline\n";
  //flatex << "{$y$} & range &  mean &  rms & yield  & s/e & $\\chi^2$/ndf & mean & rms & $\\sigma_{\\rm{av.}}$ & $\\sigma_{\\rm{fit}}$ & $\\Delta$ \\\\ \\hline \n";
  sprintf(ss,"%4.1f-%4.1f & %4.0f-%4.0f &        &       & %7.0f$\\pm$%7.0f &       &      &      &      &      & %6.3f $\\pm$ %5.3f & %6.3f $\\pm$ %5.3f & %4.1f\\%% \\\\ %% total\n",
	  UpsRapRegEdges[0],UpsRapRegEdges[nrapregion], 
	  UpsPtBinEdges[0],UpsPtBinEdges[nptbins], 
	  //UpsPtBinCenter[j][i],
	  //UpsPtBinEdgeHi[j][i], 
	  UpsYieldPtRawTot,UpsYieldPtRawTot_e, 
	  //UpsYieldPt[speak][j][i] / UpsYieldPt_e[speak][j][i],  
	  //chisqred[j][i]
	  //WeightAv[j][i], 
	  //WeightAvErr[j][i], 
	  xsec_av_tot, xsec_av_tot_e, 
	  UpsYieldPtTot, UpsYieldPtTot_e,
	  (xsec_av_tot-UpsYieldPtTot)/UpsYieldPtTot*100
	  );


  flatex << ss;

  flatex << "\\hline\n\\end{tabular}\n";
  flatex.close();

  //bool normalize=1;
  if(normalize_){ // note lumi normlaization done before
    /// normalize yield
    for(int j=0; j<nrapregion; j++) {
      for(int i=0; i<nptbins; i++) {
	double binw =  UpsPtBinEdges[i+1]-UpsPtBinEdges[i];
	for(int k=0; k<npeak; k++) {
	  UpsYieldPt   [k][j][i] *= 1./binw;
	  UpsYieldPt_e [k][j][i] *= 1./binw;
	  UpsYieldPt_eh[k][j][i] *= 1./binw;
	  UpsYieldPt_el[k][j][i] *= 1./binw;
	}  
      }
    }
  }
  
  //save binning to file
  gDirectory->Cd("../"+plt_);
    TH2F* hbin = new TH2F("pt_y_bins","pt_y_bins", nptbins, UpsPtBinEdges, nrapregion, UpsRapRegEdges);
  for(int i=0;i<nptbins;i++)
    for(int j=0;j<nrapregion;j++) {
      hbin->SetBinContent(i+1,j+1,UpsYieldPt[speak][j][i]);
    }
  hbin->Draw("col");
  hbin->Write();
  TCanvas b;
  //b.SaveAs("aa.gif");
  TString pname = TString::Format("%sbinning_ups%dS%s",figs_.Data(),speak+1,ext_.Data());
  b.SaveAs(pname);

  /// plot cross section
  TGraphAsymmErrors *fitgr[nrapregion];
  //for(int j=0; j<npeak; j++) {  
  for(int Y=0; Y<nrapregion; Y++) {
    
    double yield[nptbins], err[nptbins], err_h[nptbins], err_l[nptbins];
    double UpsPtBinCenterY[nptbins], UpsPtBinEdgeLoY[nptbins], UpsPtBinEdgeHiY[nptbins]; 
    for(int i=0; i<nptbins; i++) {
      yield[i] = UpsYieldPt   [speak][Y][i];
      err  [i] = UpsYieldPt_e [speak][Y][i];
      err_h[i] = UpsYieldPt_eh[speak][Y][i];
      err_l[i] = fabs(UpsYieldPt_el[speak][Y][i]);//low values come with negative sign
      UpsPtBinCenterY[i]=UpsPtBinCenter[Y][i];
      UpsPtBinEdgeLoY[i]=UpsPtBinEdgeLo[Y][i];
      UpsPtBinEdgeHiY[i]=UpsPtBinEdgeHi[Y][i]; 
    }
    
    fitgr[Y] = new TGraphAsymmErrors(nptbins,   UpsPtBinCenterY, yield,
				     UpsPtBinEdgeLoY, UpsPtBinEdgeHiY, 
				     err_l, err_h);
    gDirectory->Cd("../"+plt_);
    TString bla = TString::Format("cvs_rap_%d",Y);
    TCanvas c(bla,bla); c.cd();
    double x[2] = {0, 28};
    double y[2] = {0.0008,2};
    TGraph frame(2,x,y);
    TString ytitle = 
      TString::Format("#frac{d#sigma}{dp_{T}} . BR(#Upsilon(%dS)#rightarrow#mu#mu) [nb/GeV]",speak+1);
    frame.SetTitle( "" );
    frame.GetXaxis()->SetTitle("p_{T} (#mu#mu) [GeV]");
    frame.GetYaxis()->SetTitle(ytitle);
    frame.Draw("AP");
    fitgr[Y]->Draw("P");
    gPad->SetLogy();
    TString htmp = TString::Format("xsection_ups%dS_rap%d_pt",speak+1,Y);
    fitgr[Y]->SetName(htmp);
    fitgr[Y]->Write();
    c.SaveAs(TString::Format("%sxsecdiff_ups%dS_rap%d%s",figs_.Data(),speak+1,Y,ext_.Data()));
  }
}


/* plot the mass distribution and fitted model
 */
void plotMass(TString hname, RooWorkspace& w, RooAbsData *data, RooAbsPdf* pdf, TString bintitle) {

  gROOT->SetStyle("Plain");

  //float yield[4];
  bool dataonly = (data && !pdf);
  if(!data)
    data = (RooDataSet*)w.data("data");
  if(!pdf)
    pdf = (RooAbsPdf*)w.pdf("pdf");

  /*
  RooRealVar* psigma2 = ( RooRealVar*)pdf->getParameters(data)->find("sigma2");
  RooRealVar* psigmaFraction = ( RooRealVar*)pdf->getParameters(data)->find("sigmaFraction");
  //psigmaFraction->setVal(1.0); 
  //  psigma2->setConstant(kTRUE);  
  cout << "xxx " << psigmaFraction->getVal() << " " << psigma2->getVal() << "\n";
  RooArgSet* pars = (RooArgSet*)w.pdf("pdf")->getParameters(*w.data("data"));//->selectByAttrib("Constant",kFALSE);
  pars->Print("v");
  */

  //RooAbsPdf *pdf_ = (RooAbsPdf*)w.pdf("pdf");
  //pdf_->SetName("pdf");

  RooPlot* frame;
  if(printPars) 
    frame = w.var("invariantMass")->frame(Bins(nMassBins_));
  else
    frame = w.var("invariantMass")->frame(Bins(nMassBins_));
  //data->plotOn(frame, Name("aaa"));
  //pdf->plotOn(frame);
  //RooHist* hresid = frame->pullHist("aaa","pdf");
  //hresid->SetLineColor(2);

  //typedef RooAbsData::EType dataet;
  //dataet etype = data->isWeighted() ? RooAbsData::SumW2 : RooAbsData::Auto;
  if(!data->isWeighted()) 
    //data->plotOn(frame,Name("theData"), MarkerStyle(7), DrawOption("E")); // DataError(etype));
    data->plotOn(frame,Name("theData"), MarkerStyle(7), DrawOption("E")); // DataError(etype));
  else
    data->plotOn(frame,Name("theData"), MarkerStyle(7), DrawOption("E"), DataError(RooAbsData::SumW2));

  if(!dataonly) {

    pdf->plotOn(frame,Name("thePdf"));
    RooArgSet * pars = pdf->getParameters(data);
    if(printPars) {
      int nfloatpars = pars->selectByAttrib("Constant",kFALSE)->getSize(); 
      double mychsq = frame->chiSquare("thePdf","theData",nfloatpars); 
      // double chi2red = fitChisq(w,(RooDataSet*)data,pdf);
      double myndof = frame->GetNbinsX() - nfloatpars;
      // double fitprob = TMath::Prob(mychsq, myndof);
      //printf("xxx %f  \t %f \t %f \n", chi2red, chi2reda,fitprob);

      pdf->paramOn(frame,Layout(0.55,0.9,0.9),Format("NE")/*,Format("NEA",AutoPrecision(1))*/, 
		   //Label(Form("#chi^{2}/ndf = %5.3f", chi2red
		   Label(Form("#chi^{2}/ndf = %2.0f/%2.0f", mychsq, myndof
			      )));

      pdf->plotOn(frame,Name("theBgd"),Components("bkg"),
		  LineStyle(1),LineWidth(2),LineColor(16)) ; 
      
      pdf->plotOn(frame,Components(TString::Format("gauss%sS1",PEAK_.Data())),
		  FillStyle(1001), VLines(),FillColor(kRed-8), DrawOption("F")) ; 
      pdf->plotOn(frame,Components(TString::Format("gauss%sS2",PEAK_.Data())),
		  FillStyle(1001),FillColor(kRed-5), DrawOption("F")); 
      pdf->plotOn(frame,Components(TString::Format("gauss%sS1",PEAK_.Data())),
		  LineColor(kWhite),LineWidth(1));
      
    }
    delete pars;

  }
    
  pdf->plotOn(frame);
  if(!data->isWeighted())
    data->plotOn(frame,Name("theData"), MarkerStyle(7), DrawOption("E")); // DataError(etype));
  else
    data->plotOn(frame,Name("theData"), MarkerStyle(7), DrawOption("E"), DataError(RooAbsData::SumW2));

  //data->plotOn(frame,Name("theData"), MarkerStyle(7), DrawOption("E"), DataError(etype));
  //  data->plotOn(frame,Name("theData"), MarkerStyle(7), DrawOption("E"), 
  //	       DataError(RooAbsData::SumW2)); //needed in case of weights

  //  frame->SetNdivisions(1020,"X");
  //  frame->SetNdivisions(1020,"Y");

  frame->SetTitle( "" );
  frame->GetXaxis()->SetTitle("#mu^{+}#mu^{-} mass (GeV/c^{2})");
  frame->GetYaxis()->SetTitleOffset(1.65);
  frame->GetYaxis()->SetTitleSize(0.033);

  ///fit residuals
  RooPlot* rframe = w.var("invariantMass")->frame();
  RooHist* hresid = frame->pullHist();
  hresid->SetLineColor(2);
  rframe->addObject(hresid);
  int nresid=5;
  rframe->SetTitle("");
  rframe->GetXaxis()->SetTitle("");
  rframe->GetYaxis()->SetLabelSize(0.12);
  rframe->GetYaxis()->SetNdivisions(2*nresid+1);
  rframe->SetAxisRange(-1*nresid,nresid,"Y");
  rframe->SetMarkerStyle(7);

  TCanvas c; c.cd();

  TPaveLabel *t2 = new TPaveLabel(0.01,0.1,0.07,0.9, "Residual pull", "brNDC"); 
  t2->SetFillColor(0);
  t2->SetTextAngle(90);
  t2->SetTextSize(0.2);
  t2->SetBorderSize(0);

  TPad *p1 = new TPad("p1","p1",0.1,0.25,0.9,0.901);// ,0,0,0); 
  //p1->SetBottomMargin(0.);
  p1->SetBorderMode(0); 
  p1->Draw(); 
   
  TPad *p2 = new TPad("p2","p2",0.1,0.1,0.9,0.24); 
  p2->SetTopMargin(0.);    
  p2->SetBorderMode(0); 
  p2->SetTicks(1,0); 
  p2->Draw(); 

 TPaveText *cms;
  if(printPars) {
    cms = new TPaveText(0.15,0.84,0.85,0.9,"brNDC");
    cms->AddText(TString::Format("CMS preliminary    #sqrt{s}=7TeV    L=%2.0fpb^{-1} %12s %s",(lumi_/1000.),"",bintitle.Data()));
    cms->SetTextSize(0.03);
  } else {
    cms = new TPaveText(0.6,0.7,0.85,0.85,"brNDC");  
    cms->AddText("CMS");
    cms->AddText(TString::Format("#sqrt{s}=7TeV,  L=%2.0f pb^{  -1} ",(lumi_/1000.)));
    cms->SetTextSize(0.045);
  }
  cms->SetFillColor(0);
  cms->SetBorderSize(0);
  cms->SetTextAlign(12);
  cms->Draw();

  if(printPars) {
    p1->cd();
    frame->Draw();
    p2->cd();
    rframe->Draw();
    t2->Draw();
  } else {
    frame->Draw();
    cms->Draw();
  }

  //gPad->SetLogy();
  //frame->Draw();
  c.SaveAs(figs_+hname+ext_);
  gDirectory->Cd("../"+plt_);
  frame->SetName(hname);
  frame->Write(hname);
}

/* fit the mass distribution
- save fit results to file
 */
void fitYield(RooWorkspace& w, TString name, RooAbsData* data) {
  
  cout << "fitting the upsilon mass...\n" << std::flush;

  //if(!data)
  //  data = (RooDataSet*)w.data("data");

  /// reset the fit parameters, retune the yields
  setDefaultFitParValues(w,data);

  //RooFitResult* fitres = w.pdf("pdf")->fitTo(*data, Save(), Extended(kTRUE), Minos(kTRUE), SumW2Error(data->isWeighted()));//,Range(mmin_,mmax_));
  RooFitResult* fitres = w.pdf("pdf")->fitTo(*data, Save(), Extended(kTRUE), SumW2Error(data->isWeighted()));//,Range(mmin_,mmax_));
			  
  TString fres_n("fit_result_");
  fres_n += data->isWeighted()?"corrected":"raw";
  fres_n += ("_"+name);
  fitres->SetName(fres_n);
  gDirectory->Cd("../"+res_);
  fitres->Write();
  fitres->Print();
  cout << "\tsaved results in " << fres_n << "\n" << std::flush;
  //ws->pdf("pdf")->fitTo(ws->data("data")) ;//, Save(), Extended(kTRUE));
}

/* read the data from ttree
 */
void readData(RooWorkspace& w, TString fname) {
  TFile f(fname,"read");
  gDirectory->Cd(fname+":/"+dirname_);
  TTree* theTree     = (TTree*)gROOT->FindObject(treeName);
  //theTree->Print();
  RooRealVar* mass   = w.var("invariantMass");
  RooRealVar* upsPt  = new RooRealVar("upsPt","p_{T}(#Upsilon)",0,100,"GeV"); 
  RooRealVar* weight = new RooRealVar("weight",  "weight"  ,0,100);
  RooRealVar* dummy  = new RooRealVar("upsRapidity", "upsRapidity"  ,-2.5,2.5);
  RooRealVar* dummyb = new RooRealVar("maxMuEta",  "maxMuEta"  ,-2.7,2.7);
  // RooRealVar* dummyc = new RooRealVar("SampleFlag",  "SampleFlag"  ,0,4);
  RooDataSet* data0, *data;
  //if(isMC) {
  //data0 = new RooDataSet("data","data",theTree,RooArgSet(*mass,*upsPt,*weight,*dummy,*dummyb,*dummyc)); 
  //data = ( RooDataSet*) data0->reduce(Cut("SampleFlag>0"));//,EventRange(0,500000));
  //data = ( RooDataSet*) data0;//->reduce(EventRange(0,1000),Cut("invariantMass<9.8"));
  //}
  //else {
  RooArgSet cols(*mass,*upsPt,*weight,*dummy,*dummyb);
  data0 = new RooDataSet("data","data",theTree,cols); 
  //  }
  //data = ( RooDataSet*) data0->reduce(EventRange(0,100000));//,Cut("invariantMass<11"));
  TString mcut = TString::Format("invariantMass>%f && invariantMass<%f",mmin_,mmax_);
  data =  ( RooDataSet*) data0->reduce(Cut(mcut));
  w.import(*data);
  data->Print();
  f.Close();
}


/* define the fit model
+ signal
for each upsilon 1S/2S/3S: double gaussian, or gaussian + 'crystal ball' 
at each peak, the two gaussians have common center
default fitting parameters:
- mass of upsilon 1S floats
- mass shift between three peaks are fixed to pdg (or common scale factor allowed to float)
- the three peaks have a common shape 
- relative rates of three peaks float
+ background
- second order polynominal
*/
void buildPdf(RooWorkspace& w) {

  RooRealVar* mass  = new RooRealVar("invariantMass","#mu#mu mass",mmin_,mmax_,"GeV/c^{2}"); 

  const double M1S =  9.460;  //upsilon 1S pgd mass value
  const double M2S = 10.023;  //upsilon 2S pgd mass value
  const double M3S = 10.355;  //upsilon 3S pgd mass value

  RooRealVar *mean    = new RooRealVar("mass_mean","#Upsilon mean",M1S,M1S-0.3,M1S+0.3,"GeV");
  RooRealVar *shift21 = new RooRealVar("shift2","mass diff #Upsilon(1,2S)",M2S-M1S);
  RooRealVar *shift31 = new RooRealVar("shift3","mass diff #Upsilon(1,3S)",M3S-M1S);
  //RooRealVar *shift21 = new RooRealVar("shift2","mass diff #Upsilon(1,2S)",M2S-M1S,M2S-M1S-0.5,M2S-M1S+0.5);
  //RooRealVar *shift31 = new RooRealVar("shift3","mass diff #Upsilon(1,3S)",M3S-M1S,M3S-M1S-0.5,M3S-M1S+0.5);
  //shift21->setConstant(kTRUE); 
  //shift31->setConstant(kTRUE); 

  RooRealVar *mscale  = new RooRealVar("mscale","mass scale factor",1.,0.7,1.3);
  mscale->setConstant(kTRUE); /* the def. parameter value is fixed in the fit */

  RooFormulaVar *mean1S = new RooFormulaVar("mean1S","@0",
  					    RooArgList(*mean));
  RooFormulaVar *mean2S = new RooFormulaVar("mean2S","@0+@1*@2",
					    RooArgList(*mean,*mscale,*shift21));
  RooFormulaVar *mean3S = new RooFormulaVar("mean3S","@0+@1*@2",
					    RooArgList(*mean,*mscale,*shift31));

  RooRealVar *sigma1 = new RooRealVar("sigma1","Sigma_1",0.06,0.01,0.30);
  RooRealVar *sigma2 = new RooRealVar("sigma2","Sigma_2",0.10,0.01,0.30);
  //RooRealVar *sigma2 = new RooRealVar("sigma2","Sigma_2",0.05,0.0,0.30);

  ///resolution
  RooFormulaVar *reso1S = new RooFormulaVar("reso1S","@0"             ,RooArgList(*sigma1));
  RooFormulaVar *reso2S = new RooFormulaVar("reso2S","@0*10.023/9.460",RooArgList(*sigma1));
  RooFormulaVar *reso3S = new RooFormulaVar("reso3S","@0*10.355/9.460",RooArgList(*sigma1));
  
  /// to describe final state radiation tail on the left of the peaks
  RooRealVar *alpha  = new RooRealVar("alpha","tail shift",2.0,0.2,4);
  //RooRealVar *alpha  = new RooRealVar("alpha","tail shift",1.7,0.2,4);
  RooRealVar *npow   = new RooRealVar("npow","power order",1,1,3);
  //npow ->setConstant(kTRUE);
  //alpha->setConstant(kTRUE);
 
  /// relative fraction of the two peak components 
  RooRealVar *sigmaFraction = new RooRealVar("sigmaFraction","Sigma Fraction",0.3,0.,1.);

  /// Upsilon 1S
  //  RooGaussian *gauss1S1 = new RooGaussian("gauss1S1","1S Gaussian_1",*mass,*mean1S,*sigma1);
  //  RooGaussian *gauss1S2 = new RooGaussian("gauss1S2","1S Gaussian_2",*mass,*mean1S,*sigma2);
  RooCBShape  *gauss1S1 = new RooCBShape ("gauss1S1", "FSR cb 1s", 
					  *mass,*mean1S,*reso1S,*alpha,*npow); 
  RooCBShape  *gauss1S2 = new RooCBShape ("gauss1S2", "FSR cb 1s", 
					  *mass,*mean1S,*sigma2,*alpha,*npow); 
  RooAddPdf *sig1S      = new RooAddPdf  ("sig1S","1S mass pdf",
					  RooArgList(*gauss1S1,*gauss1S2),*sigmaFraction);

  /// Upsilon 2S
  //  RooGaussian *gauss2S1 = new RooGaussian("gauss2S1","2S Gaussian_1",*mass,*mean2S,*sigma1);
  //  RooGaussian *gauss2S2 = new RooGaussian("gauss2S2","2S Gaussian_2",*mass,*mean2S,*sigma2);
  RooCBShape  *gauss2S1 = new RooCBShape ("gauss2S1", "FSR cb 2s", 
					  *mass,*mean2S,*reso2S,*alpha,*npow); 
  RooCBShape  *gauss2S2 = new RooCBShape ("gauss2S2", "FSR cb 2s", 
  					  *mass,*mean2S,*sigma2,*alpha,*npow); 
  RooAddPdf *sig2S      = new RooAddPdf  ("sig2S","2S mass pdf",
					  RooArgList(*gauss2S1,*gauss2S2),*sigmaFraction);

  /// Upsilon 3S
  //  RooGaussian *gauss3S1 = new RooGaussian("gauss3S1","3S Gaussian_1",*mass,*mean3S,*sigma1);
  //  RooGaussian *gauss3S2 = new RooGaussian("gauss3S2","3S Gaussian_2",*mass,*mean3S,*sigma2);
  RooCBShape  *gauss3S1 = new RooCBShape ("gauss3S1", "FSR cb 3s", 
					  *mass,*mean3S,*reso3S,*alpha,*npow); 
  RooCBShape  *gauss3S2 = new RooCBShape ("gauss3S2", "FSR cb 3s", 
					  *mass,*mean3S,*sigma2,*alpha,*npow); 
  RooAddPdf *sig3S      = new RooAddPdf  ("sig3S","3S mass pdf",
					  RooArgList(*gauss3S1,*gauss3S2),*sigmaFraction);
  
  /// Background
  RooRealVar *bkg_a1  = new RooRealVar("bkg_a1", "background a1", 0, -1, 1);
  RooRealVar *bkg_a2  = new RooRealVar("bkg_a2", "background a2", 0, -1, 1);
  RooAbsPdf  *bkgPdf  = new RooChebychev("bkg","linear background",
					 *mass, RooArgList(*bkg_a1,*bkg_a2));
  /// Combined pdf
  int nt = 30000;

  //bool fitfraction = true;
  RooFormulaVar *nsig1f,*nsig2f,*nsig3f;

  if(!fitfraction) {
    RooRealVar *nsig1 = new RooRealVar("nsig1","signal 1s yield",nt*0.25*0.6,0,10*nt); 
    RooRealVar *nsig2 = new RooRealVar("nsig2","signal 2s yield",nt*0.25*0.3,0,10*nt); 
    RooRealVar *nsig3 = new RooRealVar("nsig3","signal 3s yield",nt*0.25*0.1,0,10*nt); 
    //RooRealVar *nbkgd = new RooRealVar("nbkgd","brackground yield",nt*0.75,0,10*nt); 

    nsig1f = new RooFormulaVar("nsig1f","@0", RooArgList(*nsig1));
    nsig2f = new RooFormulaVar("nsig2f","@0", RooArgList(*nsig2));
    nsig3f = new RooFormulaVar("nsig3f","@0", RooArgList(*nsig3));
  //  RooAbsPdf  *pdf   = new RooAddPdf ("pdf","total signal+background pdf", 
  //				    RooArgList(*sig1S,*sig2S,*sig3S,*bkgPdf),
  //				    RooArgList(*nsig1,*nsig2,*nsig3,*nbkgd));

  // pdf = N1.S1 + N2.S2 + N3.S3 + Nb.Sb
  } else { 
    // pdf = Ns.[ (1-h).S1 + h.g.S2 + h.(1-g).S3 ] + Nb.Sb
    // Ns=N1+N2+N3, h=(N2+N3)/Ns==f23os, g=N2/(N2+N3)==f2o23
    RooRealVar *nsig  = new RooRealVar("nsig","signal 1+2+3 yield",nt*0.25,0,10*nt); 
    
    RooRealVar *f23os = new RooRealVar("f23os","2+3/1+2+3 sig fration",0.4,0,1); 
    RooRealVar *f2o23 = new RooRealVar("f2o23","2/2+3 sig fration",0.75,0,1); 
    nsig1f = new RooFormulaVar("nsig1f","@0*(1.-@1)", RooArgList(*nsig,*f23os));
    nsig2f = new RooFormulaVar("nsig2f","@0*@1*@2", RooArgList(*nsig,*f23os,*f2o23));
    nsig3f = new RooFormulaVar("nsig3f","@0*@1*(1.-@2)", RooArgList(*nsig,*f23os,*f2o23));
  }

  RooRealVar *nbkgd = new RooRealVar("nbkgd","brackground yield",nt*0.75,0,10*nt); 
  
  RooAbsPdf  *pdf   = new RooAddPdf ("pdf","total signal+background pdf", 
				     RooArgList(*sig1S,*sig2S,*sig3S,*bkgPdf),
				     RooArgList(*nsig1f,*nsig2f,*nsig3f,*nbkgd));
 
  bool nobg = isMC;
  bool varycb = false;
  bool singlecb = true;
  //bool linearbg = true;

  if(linearbg_)
    bkg_a2->setConstant(kTRUE);    

  if(!varycb) {
    alpha->setConstant(kTRUE);
  }
  npow ->setConstant(kTRUE);

  if(nobg) {
    nbkgd->setVal(0); nbkgd->setConstant(kTRUE);
    bkg_a1->setConstant(kTRUE);    
    bkg_a2->setConstant(kTRUE);    
  }

  if(singlecb) {
    sigmaFraction->setVal(1.); sigmaFraction->setConstant(kTRUE);
    sigma2->setConstant(kTRUE);
  }

  w.import(*pdf);
}



/* cache the default fit parameters
 */
void getInitialParValues(RooWorkspace& w) {
  RooArgSet* pars = (RooArgSet*)w.pdf("pdf")->getParameters(*w.data("data"));//->selectByAttrib("Constant",kFALSE);
  TIterator* coefIter = pars->createIterator() ;
  RooRealVar* coef ;
  while((coef = (RooRealVar*)coefIter->Next())) {
    cout << "caching parameter " << coef->GetName() << "\t" << flush;
    coef->Print();
    initParVal_.push_back(coef->getVal());
    initParName_.push_back(coef->GetName());
  }  
}

/* reset initial parameters of the fit pdf
- to their default values
- if the dataset is specified, retune the yield to counting-based estimates
 */
void setDefaultFitParValues(RooWorkspace& w,RooAbsData* data) {
  if(!restorepars)
    return;
  RooArgSet* pars = (RooArgSet*)w.pdf("pdf")->getParameters(*w.data("data"));
  //->selectByAttrib("Constant",kFALSE);
  RooRealVar* par = 0;
  for(int i=0; i<pars->getSize();i++) {
    par = (RooRealVar*)pars->find(initParName_[i].data());
    cout << "reseting parameter " << par->GetName() << " " << par->getVal() << " ->\t" << flush;
    par->setVal(initParVal_[i]);
    par->Print();
  }
  if(data) {
    double nsig = data->sumEntries("invariantMass>9.2 & invariantMass<9.7");
    ((RooRealVar*)w.pdf("pdf")->getParameters(data)->find("nsig1"))->setMax(nsig*50);
    ((RooRealVar*)w.pdf("pdf")->getParameters(data)->find("nsig2"))->setMax(nsig*50);
    ((RooRealVar*)w.pdf("pdf")->getParameters(data)->find("nsig3"))->setMax(nsig*50);
    ((RooRealVar*)w.pdf("pdf")->getParameters(data)->find("nbkgd"))->setMax(nsig*50);
  }
}

double getWeightAv(RooWorkspace& /*w*/, RooAbsData* data) {
  RooAbsData* data0(data);
  TH1* hall = (TH1*)data0->createHistogram("weight",100);//,0,10) ;
  double wei = hall->GetMean();
  std::cout << "average signal weight " << wei << "\n";
  return wei;
}


/*
run modified configurations
 */
double oniafitterSimple(const TString finput, const TString foutput, const TString /*dfigs*/, 
			 double* UpsPtBinEdges, int nptb, 
			 double* UpsRapBinEdges, int nrapb, 
			 int ipeak) {
  //assert(nrapb==3);

  RooWorkspace* ws = new RooWorkspace("ws","upsilon mass");

  buildPdf(*ws);
  readData(*ws,finput);
  getInitialParValues(*ws);
  TFile file(foutput,"update");
  gDirectory->mkdir(res_);
  gDirectory->mkdir(plt_);
  gDirectory->Cd(res_);

  //fitYield(*ws,"total");
  //plotMass("mass_raw",*ws);
  //((RooDataSet*)ws->data("data"))->setWeightVar(*ws->var("weight"));
  //fitYield(*ws,"total",ws->data("data"));

  //RooDataSet* data_w((RooDataSet*)ws->data("data"));
  //data_w->setWeightVar(*ws->var("weight"));
  //fitYield(*ws,"total",data_w);
  //plotMass("mass_wei",*ws, data_w, ws->pdf("pdf"));
  //xsectTot(*ws);

  //xsectDif(*ws,0, UpsPtBinEdges, nptb, UpsRapBinEdges, nrapb);
  xsectDif(*ws,1, UpsPtBinEdges, nptb, UpsRapBinEdges, nrapb);

  gDirectory->Cd("../"+res_);
  //gDirectory->pwd();
  //gDirectory->ls();
  RooFitResult * fitresPt = (RooFitResult *) gROOT->FindObject(TString::Format("fit_result_corrected_rap%d_pt%d",0,1));
  TString yieldn = TString::Format("nsig%d",ipeak);
  double yield = ((RooRealVar*)fitresPt->floatParsFinal().find(yieldn))->getVal();
  double err   = ((RooRealVar*)fitresPt->floatParsFinal().find(yieldn))->getError();
  double sig = err?yield/err:0;

  printf("Y(%dS) pt:(%4.1f-%4.1f) significance:%5.2f  [nsig:%3f+-%3f]\n",ipeak,UpsPtBinEdges[0],UpsPtBinEdges[1],sig, yield, err);

  return sig;
}

double fitChisq(RooWorkspace& w, RooDataSet* data, RooAbsPdf* pdf) {
  RooPlot* frame = w.var("invariantMass")->frame(Bins(nMassBins_));
  data->plotOn(frame,Name("theData"));// MarkerStyle(7), DrawOption("E"), DataError(RooAbsData::SumW2));
  pdf->plotOn(frame,Name("thePdf"));
  int nfloatpars = pdf->getParameters(data)->selectByAttrib("Constant",kFALSE)->getSize(); 
  double chi2red = frame->chiSquare("thePdf","theData",nfloatpars); 
  return chi2red;
}

//this function is not corrected&validated
double fitProb(RooWorkspace& w, RooDataSet* data, RooAbsPdf* pdf, double& chi2, double& ndf) {
  //return -1;
  RooPlot* frame = w.var("invariantMass")->frame(Bins(nMassBins_));
  data->plotOn(frame,Name("theData"));// MarkerStyle(7), DrawOption("E"), DataError(RooAbsData::SumW2));
  pdf->plotOn(frame,Name("thePdf"));
  int nfloatpars = pdf->getParameters(data)->selectByAttrib("Constant",kFALSE)->getSize(); 
  ndf = frame->GetNbinsX();
  ndf -= nfloatpars;
  //ndf -=1;
  // double chi2red = frame->chiSquare("thePdf","theData",nfloatpars); 
  //sdouble chi2c = chi2red * ndf;

  // data->SetBins(nMassBins_);
  RooRealVar* mass = w.var("invariantMass");
  mass->setBins(nMassBins_);
  //(w.var("invariantMass"))->SetBins(nMassBins_);
  //RooDataHist* dh = data->binnedClone("myBinnedData","myBinnedData") ;
  //RooDataHist dh("dh","dh",RooArgSet(*w.var("invariantMass")),*data) ; 
  RooDataHist dh("dh","dh",RooArgSet(*mass),*data) ; 
  // dh.Print("v") ;

  RooChi2Var mychi2 ("chi2", "chi2", *pdf, dh);//, Extended(kTRUE)); 
  //Double_t 
  chi2 = mychi2.getVal(); 

  //ATTENTION: not sure why this is not exactly equal: different binning, other?
  //cout << "xxxx " << "chi2red:" << chi2red << " chi2:" << chi2 << " ndof:" << ndf << " chi2/ndof:" <<  chi2/ndf << endl; 
  //cout << "xxx fitprob  nbins:" << frame->GetNbinsX() << " npars:" << nfloatpars << " chi2:" << chi2 << " prob:" << TMath::Prob(chi2,ndf) << endl << std::flush;
  return TMath::Prob(chi2,ndf);
}


/*
  A simple fit and plot routine using the sequential fitter.
  finput : name of the root file for the data.  Will be passed to 
    readData(...)
  binning : pt binning
    0 : pt integrated
    1 : Y1S pt bins
    2 : Y2S pt bins
    3 : Y3S pt bins
    4 : Y1S pt bins for 2 rapidity bins
    5 : Y2S pt bins for 2 rapidity bins
    6 : Y3S pt bins for 2 rapidity bins
  rapBinning : rapidity binning
     0 : |y| < 2
    -1 : 5 rapidity bins using |y|
    -2 : 2 rapidity regions |y| < 1 and 1 < |y| < 2
  fitParams : file containing initial values for parameters, could be 
    "/dev/null"
  frameless : format for frameless plots, really only works with 
    binning = 1, rapBinning = 0
 */
void sequentialFitAndPlot(TString finput, int binning, int rapBinning,
			  TString fitParams, bool frameless) {

  // RooRandom::randomGenerator()->SetSeed(0);
  linearbg_ = false;
  RooWorkspace oniaWS("oniaWS", "Quarkonia mass");

  buildPdf(oniaWS);
  readData(oniaWS, finput);
  RooDataSet * data = dynamic_cast<RooDataSet *>(oniaWS.data("data"));

  RooArgSet * fr = doSeqFit(oniaWS, binning, rapBinning, fitParams, data);

  RooRealVar * invariantMass = oniaWS.var("invariantMass");
  RooCategory ptbin(*(dynamic_cast<RooCategory *>(data->get()->find("ptbin"))));
  RooCategory rapbin(*(dynamic_cast<RooCategory *>(data->get()->find("rapbin"))));

  RooArgSet cols(*invariantMass, ptbin, rapbin);
  RooAbsPdf * plotter = oniaWS.pdf("simpdf");
  if (!plotter)
    plotter = oniaWS.pdf("pdf");

  plotDataPdf(data, plotter, binning, cols, frameless);

  delete fr;
}


void buildSimPdf(RooWorkspace& oniaWS, RooArgSet const& cats) {
  RooAbsPdf * pdf = oniaWS.pdf("pdf");

  if (!pdf) {
    buildPdf(oniaWS);
    pdf = oniaWS.pdf("pdf");
  }

  //RooRealVar * invariantMass = oniaWS.var("invariantMass");
  RooCategory ptbin(*((RooCategory *)cats.find("ptbin")));
  RooCategory rapbin(*((RooCategory *)cats.find("rapbin")));

  if ((ptbin.numTypes() < 2) && (rapbin.numTypes() < 2)) return;

  //ptbin.Print("v");
  RooSimPdfBuilder simb(*pdf);
  RooArgSet * config = simb.createProtoBuildConfig();

  //RooArgSet cols(*invariantMass);
  //cols.add(cats);

  config->setStringValue("physModels", pdf->GetName());
  TString splitStr("");
  if (ptbin.numTypes() > 1) {
    splitStr += ptbin.GetName();
    splitStr += '(';
    for (int t = 1; t <= ptbin.numTypes(); ++t) {
      ptbin.setIndex(t);
      splitStr += ptbin.getLabel();
      if (t < ptbin.numTypes())
	splitStr += ',';
    }
    splitStr += ") ";
  }
  if (rapbin.numTypes() > 1) {
    splitStr += rapbin.GetName();
    splitStr += '(';
    for (int t = 1; t <= rapbin.numTypes(); ++t) {
      rapbin.setIndex(t);
      splitStr += rapbin.getLabel();
      if (t < rapbin.numTypes())
	splitStr += ',';
    }
    splitStr += ") ";
  }

  std::cout << splitStr << '\n';
  config->setStringValue("splitCats", splitStr);
  splitStr = "";
  if (ptbin.numTypes() > 1) {
    splitStr += ptbin.GetName();
    if (rapbin.numTypes() > 1) splitStr += ',';
  }
  if (rapbin.numTypes() > 1)
    splitStr += rapbin.GetName();

  splitStr += " : nsig1,nsig2,nsig3,nbkgd,bkg_a1,bkg_a2"
    ",sigma1,mass_mean";

  std::cout << splitStr << '\n';

  config->setStringValue(pdf->GetName(), splitStr.Data());
  RooSimultaneous * simpdf = simb.buildPdf(*config, 
					   cats);
  simpdf->SetName("simpdf");
  oniaWS.import(*simpdf, RecycleConflictNodes(), Silence());
}

// static double rapPtBins1S[] = {0., 2., 5., 8., 11., 15., 30.};
// static double rapPtBins2S[] = {0., 3., 7., 11., 15., 30.};
// static double rapPtBins3S[] = {0., 7., 12., 30.};
// static double highLowRapBins[] = {0., 1., 2.};
// static double rapidityBins[] = {0., 0.4, 0.8, 1.2, 1.6, 2.0};

void addPtBinningToData(RooDataSet * data, int binning) {
  RooRealVar * upsPt = dynamic_cast<RooRealVar *>(data->get()->find("upsPt"));
  static double fullpt[2] = {0., 30.};
  double const * ptBinArr = fullpt;
  int nbins = 2;
  if (binning == 1) {
    ptBinArr = ptdiff_1ybin_ptBinArr1S;
    nbins = sizeof(ptdiff_1ybin_ptBinArr1S)/sizeof(double);
  } else if (binning == 2) {
    ptBinArr = ptdiff_1ybin_ptBinArr2S;
    nbins = sizeof(ptdiff_1ybin_ptBinArr2S)/sizeof(double);
  } else if (binning == 3) {
    ptBinArr = ptdiff_1ybin_ptBinArr3S;
    nbins = sizeof(ptdiff_1ybin_ptBinArr3S)/sizeof(double);
  }else if (binning == 4) {
    ptBinArr = ptdiff_2ybin_ptBinArr1S;
    nbins = sizeof(ptdiff_2ybin_ptBinArr1S)/sizeof(double);
  } else if (binning == 5) {
    ptBinArr = ptdiff_2ybin_ptBinArr2S;
    nbins = sizeof(ptdiff_2ybin_ptBinArr2S)/sizeof(double);
  } else if (binning == 6) {
    ptBinArr = ptdiff_2ybin_ptBinArr3S;
    nbins = sizeof(ptdiff_2ybin_ptBinArr3S)/sizeof(double);
  }
  char buf[100];
  char const * binPat = "%0.1f__%0.1f";
  //std::cout << "n pt bins: " << nbins << '\n';
  RooThresholdCategory ptbint("ptbin", "ptbin", *upsPt, 
			      TString::Format(binPat, ptBinArr[0], ptBinArr[1]),
			      1);
  for (int bin = 1; bin < nbins; ++bin) {
    sprintf(buf, binPat, ptBinArr[bin-1], ptBinArr[bin]);
    ptbint.addThreshold(ptBinArr[bin], buf, bin);
  }
  //ptbint.Print("v");

  data->addColumn(ptbint);
}

void addRapBinningToData(RooDataSet * data, int rapBinning) {
  RooRealVar * upsRapidity = 
    dynamic_cast<RooRealVar *>(data->get()->find("upsRapidity"));
  RooFormulaVar absRap("absRap", "TMath::Abs(@0)", RooArgList(*upsRapidity));

  RooAbsReal * rap = upsRapidity;
  if (rapBinning < 1)
    rap = &absRap;
  static double fullRap[2] = {0.0, 2.0};
  double const * rapBinArr = fullRap;
  int nbins = 2;
  int offset = 0;
  if (abs(rapBinning) == 1) {
    rapBinArr = ydiff_rapBinArr;
    nbins = sizeof(ydiff_rapBinArr)/sizeof(double);
  } else if (abs(rapBinning) == 2) {
    rapBinArr = ptdiff_2ybin_rapBinArr;
    nbins = sizeof(ptdiff_2ybin_rapBinArr)/sizeof(double);
  }
  if (rapBinning < 0)
    offset = nbins/2;
  
  RooThresholdCategory rapbint("rapbin", "rapbin", *rap, 
			       TString::Format("rapbin_%i", 1),
			       1);
  for (int bin = 1+offset; bin < nbins; ++bin) {
    rapbint.addThreshold(rapBinArr[bin], TString::Format("rapbin_%i", bin), 
			 bin);
  }
  // rapbint.Print("v");
  data->addColumn(rapbint);
}

/*
  doSeqFit : do a sequential fit to a dataset.
  
  oniaWS : RooWorkspace with the pdf and data in it
  binning : which resonance to bin for, 0 integrated, 1 Y1S, 2 Y2S, 3 Y3S
  rapBinning : specifiy rapidity binning
      Right now this isn't implemented corretly.  Use 0 for [-2, 2].
  fitParams : file name with the inital state parameters for the fitter.
  data : the dataset.  If data == 0 then it will get "data" from the workspace.
  fout : where to print out the results of the fit in the format
     par_name par_value par_error par_initialValue
  doPlots : if true then it will do the plotting of the fit as well.
 */
RooArgSet * doSeqFit(RooWorkspace& oniaWS, int binning, int rapBinning,
		     TString fitParams, RooDataSet * data, ostream& fout,
		     bool doPlots) {
  static TRandom3 rnd(0);

  RooRealVar * invariantMass = oniaWS.var("invariantMass");
  RooCategory * sig = oniaWS.cat("sig");
  if (!data) {
    data = dynamic_cast<RooDataSet *>(oniaWS.data("data"));
  }

  addPtBinningToData(data, binning);
  RooCategory ptbin(*(dynamic_cast<RooCategory *>(data->get()->find("ptbin"))));

  addRapBinningToData(data, rapBinning);
  RooCategory rapbin(*(dynamic_cast<RooCategory *>(data->get()->find("rapbin"))));

  RooAbsPdf * pdf = oniaWS.pdf("pdf");
  RooArgSet cols(*invariantMass, ptbin, rapbin);

  RooAbsPdf * fitter = oniaWS.pdf("simpdf");
  if (!fitter) {
    buildSimPdf(oniaWS, cols);
    fitter = oniaWS.pdf("simpdf");
  }

  Roo1DTable * trueYields = 0;
  if (!fitter) {
    fitter = pdf;
    if (sig)
      trueYields = data->table(*sig);
  } else if (sig)
    trueYields = data->table(RooArgSet(*sig,ptbin,rapbin));

  RooArgSet * pars = fitter->getParameters(cols);
  int covQual=0;
  double edm=10.;

  pars->readFromFile(fitParams.Data(), 0, 0, kTRUE);
  RooArgSet truePars;
  truePars.addClone(*pars);

  std::cout << "data contains: " << data->numEntries() << " total weight: "
	    << data->sumEntries() << '\n';
  if (trueYields)
    trueYields->Print("v");

  // return pars;

  // fitter->Print("v");
  // return pars;

  TString className(fitter->IsA()->GetName());
  std::cout << className << '\n';
  RooFitResult * fr = 0;
  //data->Print("v");
  TList * dss = 0;
  RooSuperCategory * simCat = 0;
  if (className == "RooSimultaneous") {
      RooSimultaneous * simfitter = dynamic_cast<RooSimultaneous *>(fitter);
    dss = data->split(simfitter->indexCat());
    simCat = new RooSuperCategory(dynamic_cast<RooSuperCategory const &>(simfitter->indexCat()));
  } else {
    dss = new TList();
    dss->Add(data);
  }

  // dss->Print();
  TIter nextds(dss);
  RooDataSet * fitdata;
  RooAbsPdf * tmpFitter;
  //pars->readFromFile(fitParams.Data(), 0, 0, kTRUE);
  while ((fitdata = dynamic_cast<RooDataSet *>(nextds()))) {
    //fitdata->Print();
    tmpFitter = fitter;
    ptbin.setIndex(fitdata->get()->getCatIndex("ptbin",1));
    rapbin.setIndex(fitdata->get()->getCatIndex("rapbin",1));
    if (className == "RooSimultaneous") {
      RooSimultaneous * simfitter = dynamic_cast<RooSimultaneous *>(fitter);
      tmpFitter = simfitter->getPdf(fitdata->GetName());
      simCat->setLabel(fitdata->GetName());
      simCat->Print();
      rapbin.setIndex(simCat->inputCatList().getCatIndex("rapbin", 1));
      ptbin.setIndex(simCat->inputCatList().getCatIndex("ptbin", 1));
    	//((dynamic_cast<RooSimultaneous *>(fitter))->indexCat().getLabel());
    }


    rapbin.Print();
    ptbin.Print();
    tmpFitter->Print();
    fr = tmpFitter->fitTo(*fitdata, Extended(kTRUE),
			  //NumCPU(2),
			  SumW2Error(kFALSE),
			  Minos(kFALSE),
			  Save(kTRUE));
    covQual = fr->covQual();
    edm = fr->edm();

    std::cout << "fit completed covQual: " << covQual << " edm: " << edm
	      << '\n';
    if ( (fr) && (covQual==3) && (edm <= 0.01) ) {
      if (fitdata->isWeighted()) {
	delete fr;
	fr = tmpFitter->fitTo(*fitdata, Extended(kTRUE),
			      SumW2Error(kTRUE),
			      Minos(kFALSE),
			      Save(kTRUE));
	std::cout << "weighted fit complete: " << fr->covQual() 
		  << " edm: " << fr->edm()
		  << '\n';
      }
      double /*avgWgt = data->sumEntries()/data->numEntries(),*/ trueVal = 0.;
      TIter nextPar(fr->floatParsFinal().createIterator());
      RooRealVar * par = 0;
      fout << "nll " << fr->minNll() << ' ';
      while ((par = (RooRealVar *)nextPar())) {
	trueVal = truePars.getRealValue(par->GetName(), 0.);
	if ( (trueYields) && (par->GetName()[0] == 'n') ) {
	  TString parName(par->GetName());
	  // int sep = parName.First('_');
	  TString binName = TString::Format("%s;%s", ptbin.getLabel(),
					    rapbin.getLabel());
	  TString resName = "Bkg";
	  if (parName(1) == 's') {
	    resName = parName(TRegexp("[123]"));
	    resName = "Y" + resName + "S";
	  }
	  //std::cout << parName << ' ' << binName << ' ' << resName << '\n';
	  TString rowname = resName;
	  if (className == "RooSimultaneous")
	    rowname = "{" + resName + ";" + binName + "}";
	  trueVal = trueYields->get(rowname);
	  std::cout << "True yield name: " << rowname 
		    << " yield: " << trueVal
		    << '\n';
	}
	fout << par->GetName() << ' '
	     << par->getVal() << ' '
	     << par->getError() << ' '
	     << trueVal << ' ';
      }
    }
    if (className == "RooSimultaneous") {
      delete fitdata;
      delete fr;
      fr = 0;
    }
  }
  fout << '\n';
  if (fr) delete fr;
  if (trueYields) delete trueYields;

  if (doPlots)
    plotDataPdf(data, fitter, binning, cols);
  return pars;
}

/*
  plotDataPdf : overlay data and pdf projections can do simultaneous pdf

  data : data to use
  pdf : pdf to use
  resonance : resonance that was used, 0 none, 1 Y1S, 2 Y2S, 3 Y3S

 */
void plotDataPdf(RooDataSet * data, RooAbsPdf * const pdf, int resonance,
		 RooArgSet const& projVars, bool frameless) {
  RooCategory ptbin(*(dynamic_cast<RooCategory *>(data->get()->find("ptbin"))));
  RooCategory rapbin(*(dynamic_cast<RooCategory *>(data->get()->find("rapbin"))));
  RooRealVar * invariantMass = 
    dynamic_cast<RooRealVar*>(data->get()->find("invariantMass"));

  TString className(pdf->IsA()->GetName());
  std::cout << className << '\n';
  gStyle->SetOptTitle(0);
  std::cout << "starting to plot" << '\n';
  RooArgSet sliceVars(ptbin, rapbin);
  RooPlot * mf = 0;
  TLatex l;
  l.SetNDC(true);
  l.SetTextAlign(31);
  l.SetTextFont(42);
  char buf[255];
  TCanvas * c;
  if (className == "RooSimultaneous") {
    TList * dss = 0;
    dss = data->split((dynamic_cast<RooSimultaneous *>(pdf))->indexCat());
    TIter nextds(dss);
    RooDataSet * sdata = 0;
    RooAbsPdf * plotpdf = 0;
    while ((sdata = dynamic_cast<RooDataSet *>(nextds()))) {
      plotpdf = 
	(dynamic_cast<RooSimultaneous *>(pdf))->getPdf(sdata->GetName());

      // TIter nextType((dynamic_cast<RooSimultaneous *>(pdf))->indexCat().typeIterator());
      // RooCatType * typ = 0;
      // Nexp = 0.;
      // while ( (typ = (RooCatType *)nextType()) ) {
      //   RooAbsPdf * plotpdf = 
      // 	(dynamic_cast<RooSimultaneous *>(pdf))->getPdf(typ->GetName());
      TString typName(sdata->GetName());
      typName.ReplaceAll("{", "");
      typName.ReplaceAll("}", "");
      typName.ReplaceAll(";", "_");
      std::cout << typName << '\n';
      ptbin.setLabel(typName.Data());
      //   if (ptbin.getIndex() < 1) continue;
      //   TString indexStr(ptbin.GetName());
      //   indexStr += "==";
      //   indexStr += ptbin.getIndex();
      //   //std::cout << "ptbin: " << ptbin.getLabel() << '\n';
      // RooDataSet * sdata = dynamic_cast<RooDataSet *>(data->reduce(indexStr));
      mf = invariantMass->frame(Bins(60), Range(8., 14.));
      //sdata->Print();
      sdata->plotOn(mf, //MarkerStyle(kFullDotMedium), 
		    Invisible(),
		    MarkerStyle(kFullCircle),
		    DrawOption("pz"),
		    XErrorSize(0));
      // sliceVars.Print("v");
      plotpdf->plotOn(mf, Components("bkg*"),
		      ProjWData(projVars, *sdata),
		      LineStyle(kDashed),
		      LineWidth(4.0),
		      LineColor(kGray+2));
      plotpdf->plotOn(mf, 
		      ProjWData(projVars, *sdata),
		      LineWidth(4.0),
		      LineColor(kBlue+2));
      sdata->plotOn(mf, //MarkerStyle(kFullDotLarge), 
		    MarkerStyle(kFullCircle),
		    MarkerSize(1.0),
		    DrawOption("pz"),
		    XErrorSize(0));
      c = new TCanvas("c_" + typName, typName, 650, 600);
      mf->Draw();
      if (frameless) {
	float ptlow, pthigh;
	std::cout << sdata->GetName() << '\n';
	sscanf(sdata->GetName(), "{%4f__%4f}", &ptlow, &pthigh);
	// uncomment this to make the plots for the frameless plot
	double maxy = 700.0;
	switch (ptbin.getIndex()) {
	case 1:
	case 2:
	case 3:
	  maxy = 700.0;
	  break;
	case 4:
	case 5:
	case 6:
	  maxy = 550.0;
	  break;
	case 7:
	case 8:
	case 9:
	  maxy = 325.0;
	  break;
	case 10:
	case 11:
	case 12:
	  maxy = 240.0;
	  break;
	case 13:
	case 14:
	case 15:
	  maxy = 120.0;
	  break;
	}
	mf->SetMaximum(maxy);
	sprintf(buf, "(%0.0f < p_{T}^{#mu#mu} < %0.0f) GeV/c", ptlow, pthigh);
	l.DrawLatex(1.0 - c->GetRightMargin() - 0.04, 
		    1.0 - c->GetTopMargin() - 0.09, buf);
      }
      mf->SetMinimum(1e-6);
      gPad->Modified();
      gPad->Update();
      gPad->Print(TString::Format("plots/massProj_%s_Y%iS%s.eps",
				  typName.Data(), resonance,
				  ((sdata->isWeighted()) ? "_wgt": "")));
      delete sdata;
    }

    //delete dss;
  }

  RooDataHist projData("projData", "projData", projVars, *data);
  std::cout << "plot sig+bg" << '\n';
  mf = invariantMass->frame(Bins(80), Range(8., 12.));
  mf->SetMinimum(1e-6);
  data->plotOn(mf, //MarkerStyle(kFullDotMedium), 
	       Invisible(),
	       MarkerStyle(kFullCircle),
	       DrawOption("pz"),
	       XErrorSize(0));
  pdf->plotOn(mf, Components("bkg*"),
	      //Project(RooArgSet(ptbin)),
	      //Normalization(Nexp, RooAbsReal::NumEvent),
	      ProjWData(projVars, projData), 
	      LineStyle(kDashed),
	      LineColor(kGray+2),
	      LineWidth(4.));
  pdf->plotOn(mf,
	      //Project(RooArgSet(ptbin)),
	      //Normalization(Nexp, RooAbsReal::NumEvent),
	      ProjWData(projVars, projData),
	      LineWidth(4.),
	      LineColor(kBlue+2));
  std::cout << "plot bg" << '\n';
  data->plotOn(mf, //MarkerStyle(kFullDotLarge), 
	       MarkerStyle(kFullCircle),
	       MarkerSize(0.9),
	       DrawOption("pz"), 
	       XErrorSize(0));
  c = new TCanvas("c_integrated", "integrated", 600, 600);
  gPad->SetLeftMargin(0.2);
  mf->Draw();
  mf->SetMinimum(1e-6);
  l.DrawLatex(1.0 - c->GetRightMargin() - 0.04,
	      1.0 - c->GetTopMargin() - 0.09,
	      "CMS, #sqrt{s} = 7 TeV");
  l.DrawLatex(1.0 - c->GetRightMargin() - 0.04,
	      1.0 - c->GetTopMargin() - 0.09 - l.GetTextSize()*1.2,
	      TString::Format("L = %3.0f pb^{-1}", (lumi_/1000.)));
  gPad->Update();
  mf->GetYaxis()->SetTitleOffset(1.5);
  gPad->Modified();
  gPad->Print(TString::Format("plots/massProj_integrated_Y%iS%s.eps", 
			      resonance,((data->isWeighted()) ? "_wgt": "")));
}
