#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include <stdio.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TAxis.h"
#include "TLegend.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <vector>

using namespace std;

bool verbose = 0;
bool logscale =1;
double lumirelerr = 0.11; 
const double lumiforlabel = 3.05; //pb-1

bool diffonly = true;
const double ptmax_=30;
const double ymax_=2;
const double ymin_=0.001;

TString getFileName(TString name);
double* getptbins(TString filevar, int &nbin, int center_lo_hi);
double* getUncertainty(TString filevar, int &nbin, bool hi);
double* getResults(TString filevar, int &nbin);
double* getVariation(TString filevar, TString fileref);

double* getResultsRap(TString filevar, const int nbin);
double* getUncertaintyRap(TString filevar, const int nbin, bool hi);
double* getptbinsRap(TString filevar, const int nbin, int center);

double* getptbinedges(TString filevar, int &nbin, bool ispt = false);
double* getrapbinedges(TString filevar, int &nbin);

void printline(double* diff, int nbin, TString n, double* norm=NULL);
void printptline(TString fvar);
void drawrule(int nbin);
void plotxsecgr   (TGraphAsymmErrors* xsecgr, TGraphAsymmErrors* xsecgrsta=0, TGraphAsymmErrors* xsecgrsys=0);
void plotxsecgrRap(TGraphAsymmErrors* xsecgr, TGraphAsymmErrors* xsecgrsta=0, TGraphAsymmErrors* xsecgrsys=0);
void plotpolgr(vector <TString>  polariz,  vector <TString>  polname);

int irap_;
int ipeak_;

const bool pol = true;
double npol_;

bool absrap = false;

//const double ptwidth_ = 30;  
/* un-do pt normalization when in the rapidity mode
   tbd: get remaining calls to this number dynamically from code
*/

const int nsystotable = 9;

int anamode_;
bool rapmode_;

void xsecresults(int ipeak=1, int irap=0, int anamode = 0) {

  ipeak_ = ipeak;
  irap_ = irap;

  anamode_ = anamode;
  rapmode_ = (anamode==2);
  
  if(anamode==3)
    lumirelerr = 0;
    

  gROOT->LoadMacro("setTDRStyle_modified.C");
  setTDRStyle();

  vector <TString>  sources, sysname,sysnameT;
  sources.push_back("Acc");      sysname.push_back("acceptance");     sysnameT.push_back("$A$");      
  //sources.push_back("Emuid");    sysname.push_back("muon id eff");    sysnameT.push_back("$\\varepsilon_{\\rm{muid}}$");
  //sources.push_back("Etrig");    sysname.push_back("trigger eff");    sysnameT.push_back("$\\varepsilon_{\\rm{trig}}$");    
  sources.push_back("Etreco");   sysname.push_back("trig+muid eff");  sysnameT.push_back("$\\varepsilon_{\\rm{trig,id}}$");    
  sources.push_back("ptscale");  sysname.push_back("pt scale");	      sysnameT.push_back("$S_{p}$");  
  sources.push_back("ptspec");   sysname.push_back("pt spectrum");    sysnameT.push_back("$A_{p_{T}}$");
  sources.push_back("vtxpos");   sysname.push_back("luminous region");sysnameT.push_back("$A_{\\rm vtx}$");
  sources.push_back("nofsr");    sysname.push_back("fsr");	      sysnameT.push_back("$A_{\\rm fsr}$");    
  sources.push_back("mctrue");   sysname.push_back("tnp bias");	      sysnameT.push_back("t\\&p");
  sources.push_back("tnpmcUps"); sysname.push_back("tnp jpsi/ups");   sysnameT.push_back("$\\varepsilon_{J\/\\psi,\\Upsilon}$"); 
  sources.push_back("linear812");sysname.push_back("linear bg 8-12"); sysnameT.push_back("{\\sc bg}");
  sources.push_back("ptreso");   sysname.push_back("pt resolution");  sysnameT.push_back("$\\sigma_{p}$");   
  //sources.push_back("floatmscale");   sysname.push_back("mscale");  sysnameT.push_back("$M_{\\rm scale}$");    
  //sources.push_back("floatMDiff");   sysname.push_back("mdiff");  sysnameT.push_back("$pdf_{\\Delta m}$");    
  sources.push_back("Etrk");     sysname.push_back("tracking eff");   sysnameT.push_back("$\\varepsilon_{\\rm{trk}}$");

  TString hilo   ("Acc Etrk Emuid Etrig Etreco ptscale ptreso");

  double *diff;

  TString nominal("nominal");
  TString peak;
  if(ipeak_==1) 
    peak = TString::Format("nominal");
  else
    peak = TString::Format("%ds",ipeak_);
  
  vector <TString>  polariz, polname, polnameT;
  polariz.push_back(peak);       polname.push_back("unpolarized");  polnameT.push_back("unpol");
  if(pol) {
    polariz.push_back("helT");     polname.push_back("hel. trans");   polnameT.push_back("HX T"); 
    polariz.push_back("helL");     polname.push_back("hel. long");    polnameT.push_back("HX L");  
    polariz.push_back("csT");      polname.push_back("c.s. trans");   polnameT.push_back("CS T"); 
    polariz.push_back("csL");      polname.push_back("c.s. long");    polnameT.push_back("CS L");  
  }

  npol_ = polariz.size();
  const int npol = npol_;

  for(int i=1;i<npol_;i++) {
    if(ipeak_!=1) 
      polariz[i] += TString::Format("_%ds",ipeak_);
  }  

  int nbin_;

  //cross section central values (lumi and pt normalized, as it coes from the fitter)
  double* xsec = getResults(peak,nbin_);
  const int nbin = nbin_; 

  //printline(xsec,nbin,"the cross section: ");  return;
  //cout << "number of bins:" << nbin << endl; return;

  //statistical uncertainty
  double* xsecEh = getUncertainty(peak, nbin_, 1);
  double* xsecEl = getUncertainty(peak, nbin_, 0);
  //printline(xsecEh,nbin,peak,xsec);


  //save to array
  const int nsys = sources.size(); cout << "nsys:" << nsys << endl;
  double results[nbin][nsys+2+3][2];
  for(int ibin=0; ibin<nbin; ibin++) 
    for(int isys=0;isys<nsys+5;isys++) 
      for(int i=0;i<2;i++)
	results[ibin][isys][i]=0.;

  enum err {lo=0,hi};
  for(int ibin=0; ibin<nbin; ibin++) {
    results[ibin][0][0]  = (xsec)  [ibin]; //xsection central values 
    results[ibin][1][hi] = (xsecEh)[ibin]; //stat error high
    results[ibin][1][lo] = (xsecEl)[ibin]; //stat error low
  } 

  printf("\n\nCross section results for Y(%dS)\n",ipeak_);
  drawrule(nbin);
  printptline(nominal);
  drawrule(nbin);
  //  return;

  //getptbinsA(peak, nbin_);
  //double *ptb = getptbins(peak,nbin_,0);
  //double *pt  = getptbins(peak,nbin_,1);
  printline(xsec,nbin,"cross section");
  printline(xsecEh,nbin,"statistical",xsec);
  printline(xsecEl,nbin," unc.",xsec);
  //return;  

  double syshi[nbin], syslo[nbin];
  for(int j=0; j<nbin; j++) { syshi[j]=0.; syslo[j]=0.;}

  int ishilo = lo;
  int isys=0;

  vector<TString>::iterator it,in; 
  it = sources.begin();
  in = sysname.begin();
  while(it!=sources.end()) {

    TString var, ref, name;
    ref  = nominal;
    var  = *it;
    name = *in;

    //in case of hi/lo variations, first computed lo
    if(hilo.Contains(*it)) {
      var  = *it + "Lo";
    }  

    //sources with non nominal reference
    if (*it=="mctrue") {
      ref = "tnpmc";
    } else if (*it=="tnpmcUps") {
      ref = "mctrue";
    }

    //compute the systematic from current source
    diff = getVariation(var, ref);

    // special uncertainty scaling
    for(int ibin=0; ibin<nbin; ibin++) {
      if(*it=="nofsr")
	diff[ibin]*=0.2; //take 20% of the maximal difference 
    }

    //dump results to screen
    printline(diff,nbin,name,xsec);
    //cout << "\t";
    for(int ibin=0; ibin<nbin; ibin++) {
      results[ibin][isys+2][ishilo] = (diff)[ibin]; //syst-j error low
      //printf("\t%-3.1f%%",results[ibin][isys+2][ishilo]/xsec[ibin]*100);
    }
    //cout << endl << endl;

    //combine systematic sources 
    for(int j=0; j<nbin; j++) {
      double val = diff[j]; 
      if(val>0)
	syshi[j] += pow(val,2);
      else
	syslo[j] += pow(val,2);
    }

    //in case of hi/lo variations, recompute in next iteration for hi
    if(hilo.Contains(*it)) {
      *it += "Hi";
      //*in += "Hi"; 
      *in = "";
      ishilo=hi;
    } else {
      it++;
      in++;
      isys++;
      ishilo=lo;
    }
  }

  ///extra sources (not automatically retrieved)
  double othersys[nbin];
  double lumi[nbin], external[nbin];

  //fit bias, relative errors in % !!
  double fitbias[4][3][2][16];
  //fitbias[mode][ipeak-1][irap][iptbin]
  for(int i=0;i<4;i++) 
    for(int j=0;j<3;j++) 
      for(int k=0; k<2; k++) 
	for(int l=0; l<16; l++) 
	  fitbias[i][j][k][l]=0.;

//|y|<2
//1s
fitbias[0][0][0][ 0]=0.017; // 0--30   
fitbias[0][0][0][ 1]=0.74 ; // 0--1    
fitbias[0][0][0][ 2]=0.18 ; // 1--2    
fitbias[0][0][0][ 3]=0.21 ; // 2--3    
fitbias[0][0][0][ 4]=0.27 ; // 3--4    
fitbias[0][0][0][ 5]=0.24 ; // 4--5    
fitbias[0][0][0][ 6]=0.38 ; // 5--6    
fitbias[0][0][0][ 7]=0.64 ; // 6--7    
fitbias[0][0][0][ 8]=0.93 ; // 7--8    
fitbias[0][0][0][ 9]=0.35 ; // 8--9    
fitbias[0][0][0][10]=1.0  ; // 9--10   
fitbias[0][0][0][11]=0.47 ; //10--12  
fitbias[0][0][0][12]=0.95 ; //12--14  
fitbias[0][0][0][13]=1.0  ; //14--17  
fitbias[0][0][0][14]=2.0  ; //17--20  
fitbias[0][0][0][15]=1.9  ; //20--30  
//2s
fitbias[0][1][0][ 0]=0.019; // 0--30 
fitbias[0][1][0][ 1]=0.72 ; // 0--2     
fitbias[0][1][0][ 2]=0.73 ; // 2--4  
fitbias[0][1][0][ 3]=0.72 ; // 4--6  
fitbias[0][1][0][ 4]=1.3  ; // 6--9  
fitbias[0][1][0][ 5]=1.5  ; // 9--12 
fitbias[0][1][0][ 6]=2.4  ; //12--16
fitbias[0][1][0][ 7]=5.7  ; //16--20
fitbias[0][1][0][ 8]=17   ; //20--30
//3s
fitbias[0][2][0][ 0]=0.26; // 0--30 
fitbias[0][2][0][ 1]=1.9 ; // 0--3  
fitbias[0][2][0][ 2]=4.9 ; // 3--5  
fitbias[0][2][0][ 3]=1.6 ; // 5--8  
fitbias[0][2][0][ 4]=2.9 ; // 8--12 
fitbias[0][2][0][ 5]=2.3 ; // 12--20
fitbias[0][2][0][ 6]=6.3 ; // 20--30
//|y|<1	  	  
//1s	  	  
fitbias[1][0][0][ 0]=0.12; // 0--30
fitbias[1][0][0][ 1]=0.03; // 0--2  
fitbias[1][0][0][ 2]=0.04; // 2--5  
fitbias[1][0][0][ 3]=0.08; // 5--8  
fitbias[1][0][0][ 4]=0.34; // 8--11
fitbias[1][0][0][ 5]=0.64; //11--15
fitbias[1][0][0][ 6]=0.52; //15--30
//2s	  	  
fitbias[1][1][0][ 0]=0.07; // 0--30  
fitbias[1][1][0][ 1]=0.96; // 0--3    
fitbias[1][1][0][ 2]=0.91; // 3--7    
fitbias[1][1][0][ 3]=1.2 ; // 7--11  
fitbias[1][1][0][ 4]=1.6 ; //11--15  
fitbias[1][1][0][ 5]=2.7 ; //15--30  
//3s	  	  
fitbias[1][2][0][ 0]=0.07; // 0--30  
fitbias[1][2][0][ 1]=0.33; // 0--7   
fitbias[1][2][0][ 2]=1.9 ; // 7--12  
fitbias[1][2][0][ 3]=2.7 ; //12--30  
//1<|y|<2 
//1s	  
fitbias[1][0][1][ 0]=0.08; // 0--30   
fitbias[1][0][1][ 1]=0.21; // 0--2    
fitbias[1][0][1][ 2]=0.08; // 2--5    
fitbias[1][0][1][ 3]=0.27; // 5--8    
fitbias[1][0][1][ 4]=0.76; // 8--11   
fitbias[1][0][1][ 5]=0.52; //11--15  
fitbias[1][0][1][ 6]=0.75; //15--30  
//2s	  	  
fitbias[1][1][1][ 0]=0.44; // 0--30
fitbias[1][1][1][ 1]=2.1 ; // 0--3 
fitbias[1][1][1][ 2]=1.1 ; // 3--7 
fitbias[1][1][1][ 3]=1.4 ; // 7--11
fitbias[1][1][1][ 4]=3.3 ; //11--15
fitbias[1][1][1][ 5]=6.4 ; //15--30
//3s	  	  
fitbias[1][2][1][ 0]=1.1; //   0--30
fitbias[1][2][1][ 1]=3.7; //   0--7 
fitbias[1][2][1][ 2]=3.9; //   7--12
fitbias[1][2][1][ 3]=2.4; //  12--30
//ds/dy 
//1s
fitbias[2][0][0][ 0]=0.0085; // 0.0--2.0  
fitbias[2][0][0][ 1]=0.16  ; //	0.0--0.4
fitbias[2][0][0][ 2]=0.27  ; //	0.4--0.8
fitbias[2][0][0][ 3]=0.07  ; //	0.8--1.2
fitbias[2][0][0][ 4]=0.14  ; //	1.2--1.6
fitbias[2][0][0][ 5]=0.07  ; //	1.6--2.0

 double tt;
 for(int j=0;j<3;j++) //peak
   for(int k=0; k<2; k++) //irap
     for(int l=0; l<16; l++)  //iptbin
       fitbias[3][j][k][l]=fitbias[0][j][k][l]; //syst mode3 == mode0
 
 for(int i=0;i<4;i++) //mode
   for(int j=0;j<3;j++) //peak
     for(int k=0; k<2; k++) //irap
       for(int l=0; l<16; l++) //iptbin
	 { 
	   tt = fitbias[i][j][k][l];
	   if(i!=1 && k==1)  assert(!tt);
	   if(j==1 && l>8)   assert(!tt);
	   if(j==2 && l>6)   assert(!tt);
	 }
  
 ///fsr from toy studies, relative errors
 double fsrmcsys[] = {0.008, 0.009, 0.017};
 // multiply by 3 (equiv to taking 3 sigma variation from mc)
 for(int i=0; i<3; i++)  fsrmcsys[i] *= 3.;

 double vala(0.),valb(0.),valc(0.),vald(0.);
  for(int j=0; j<nbin; j++) { 
    vala = fsrmcsys[ipeak_-1]   * xsec[j];
    valb = fitbias[anamode_][ipeak-1][irap][j]*0.01*xsec[j]; /*values given in percent*/
    //fitbias[mode][ipeak-1][irap][iptbin]
    //valb = fitbias [ipeak_-1][j]* xsec[j];
    valc = 0.016 * xsec[j]; //rho/dimuon eff
    vald = 0.006 * xsec[j]; //fit/mass constraint
    external[j] = sqrt(pow(vala,2) + pow(valb,2) + pow(valc,2) + pow(vald,2));
    lumi    [j] = lumirelerr*xsec[j];
    othersys[j] = sqrt(pow(external[j],2) + pow(lumi[j],2));
  }
  printline(external,nbin,"external",xsec);  
  printline(lumi,nbin,"luminosity",xsec);  


  //total uncertainty
  double totunchi[nbin], totunclo[nbin];
  //double totunc_butLumi_hi[nbin], totunc_butLumi_lo[nbin];
  for(int j=0; j<nbin; j++) { 
    totunchi[j]=0.; totunclo[j]=0.;
    //totunc_butLumi_hi[j]=0.; totunc_butLumi_lo[j]=0.;
  }

  drawrule(nbin);
  double extra;
  for(int j=0; j<nbin; j++) {
    extra = pow(external[j],2) + pow(lumi[j],2);
    syshi[j] +=extra;
    syslo[j] +=extra;
    totunchi[j] = sqrt(pow(syshi[j],1)+pow(xsecEh[j],2));
    totunclo[j] = sqrt(pow(syslo[j],1)+pow(xsecEl[j],2));
    syshi[j] = sqrt(syshi[j]);
    syslo[j] = sqrt(syslo[j]);
  }
  printline(syshi,nbin,"total syst",xsec);
  printline(syslo,nbin," unc.",xsec);
  drawrule(nbin);  drawrule(nbin);

  //print cross section, stat, syst, total (differential values, divided by *pt* binwidth)
  printline(xsec,nbin,"cross section");
  printline(xsecEh,nbin,"statistical",xsec);
  //printline(xsecEl,nbin," unc.",xsec);
  printline(totunchi,nbin,"total unc.",xsec);
  printline(totunclo,nbin,"",xsec);
  drawrule(nbin);  

  //polarization
  double polarization[nbin][npol][3];//ptbin,pol-type,cen/lo/hi
  double maxpol[nbin];
  for(int j=0; j<nbin; j++) { maxpol[j]=0.;}
  drawrule(nbin);
  //printline(xsec,nbin,"unpolarized");
  int ipol=0;
  it = polariz.begin();
  in = polname.begin();
  while(it!=polariz.end()) {
    double* tmpp;
    cout << *in << endl << flush;
    tmpp = getResults(*it,nbin);
    for(int ibin=0; ibin<nbin; ibin++) polarization[ibin][ipol][0] = (tmpp)[ibin]; //central
    tmpp = getUncertainty(*it, nbin_, 0);
    for(int ibin=0; ibin<nbin; ibin++) polarization[ibin][ipol][1] = (tmpp)[ibin]; //stat lo
    tmpp = getUncertainty(*it, nbin_, 1);
    for(int ibin=0; ibin<nbin; ibin++) polarization[ibin][ipol][2] = (tmpp)[ibin]; //stat hi
    tmpp = getVariation(*it, peak);
    if(*it!=peak)
      printline(tmpp,nbin,*in,xsec);
    for(int j=0; j<nbin; j++) { 
      maxpol[j]= (fabs(tmpp[j])>maxpol[j])?fabs(tmpp[j]):maxpol[j];
    }
    delete tmpp;
    it++;
    in++;
    ipol++;
  }
  printline(maxpol,nbin,"maximum",xsec);
  drawrule(nbin);  

  double *pt  = getptbins    (peak,nbin_,0);//pt bin center
  double *ptb = getptbins    (peak,nbin_,1);//pt error (lo==hi)
  double *pte = getptbinedges(peak,nbin_);  //pt/rap bin edges

  ////double *pth = getptbins(peak,nbin_,2);//hi
  //printline(pt ,nbin,"debug pt center");
  //printline(ptb,nbin,"debug pt error");
  //printline(pte,nbin,"debug pt edges");
  //printptline(nominal);
  //printline(pth,nbin,"debug pt hi");
  //cout << "number of bins:" << nbin_ << endl;

  //PRINT TABLE
  //before unc processing
  if(0)
  for (int ibin=0; ibin<nbin; ibin++) {//row
    if(ibin==0) {
      printf("pt-bin \t   pt-av.  central %4s %4s", "stat-hi","stat-lo");
      for(int isys=0;isys<nsys;isys++) {    
	printf("%15s", sysnameT[isys].Data());
      }
      printf("\n");
    } 
    if(ibin==0) {
      printf("integrated\t ");
    } else {
      printf("[%4.1f-%4.1f] %5.2f",pte[ibin-1],pte[ibin],pt[ibin]);
    }
    printf("%8.2f %8.2f%% %8.2f%%", results[ibin][0][0], results[ibin][1][0]/xsec[ibin]*100, results[ibin][1][1]/xsec[ibin]*100);

    for(int isys=0;isys<nsys;isys++) {    
      printf(" %+8.1f%%",results[ibin][isys+2][lo]/xsec[ibin]*100);
      if(hilo.Contains(sysname[isys])) 
      	printf(" %+4.1f%%",results[ibin][isys+2][hi]/xsec[ibin]*100);
    }
    cout << endl;
  }


  ///process values: normalization, error sign
  double binw[nbin], binwpt[nbin];
  for(int ibin=0; ibin<nbin; ibin++) {
    /*given the pt or rapidity bin widh depending on the analysis mode*/
    if(ibin==0)  binw[ibin]=1.;
    else         binw[ibin]=pte[ibin]-pte[ibin-1];
    binwpt[ibin]=1.;//binw[ibin]; //no need to un-normalize pt now (this is done directly when inititally reading the results)
    //if(ibin==0)  binwpt[ibin]=1.;
    //else         binwpt[ibin]=pte[ibin]-pte[ibin-1];
    if(rapmode_) binwpt[ibin]=1.; 
    /*note: pt normalization in this mode is removed in getresults() and getuncertainty()
      this quantity is used specifically to remove the normalization done in the fitter, which is pt only 
    */    
    ///printf("binwpt:%d %4.1f = %4.1f-%4.1f\n",ibin,binwpt[ibin], ptb[ibin],ptb[ibin-1]);
  }

  //printline(binw,nbin,"debug binw");
  //return;

  double hiv, lov, tmp,norm;
  for(int isys=0;isys<nsys;isys++) {//column    
    for (int ibin=0; ibin<nbin; ibin++) {//row
      norm = xsec[ibin]/100.;
      hiv = results[ibin][isys+2][hi];
      lov = results[ibin][isys+2][lo];
      //if(ibin==3) {printf("%10s a",sysname[isys].Data());printf(" hi:%+4.1f%% lo:%+4.1f%% ",hiv/norm,lov/norm);} 
      if(hilo.Contains(sysname[isys]) && hiv*lov<0) {//hi/lo with opp side
	if(hiv<lov) {tmp=hiv;hiv=lov; lov=tmp;}
      } else { //hi/lo same sign or single variation type source: assign the largest for +/-
	hiv = fabs(hiv); lov=fabs(lov);
	tmp = hiv>lov?hiv:lov; hiv=tmp; lov=tmp;
      }
      results[ibin][isys+2][lo]=lov;
      results[ibin][isys+2][hi]=hiv;
      //if(ibin==3) printf("-> hi:%+6.2f%% lo:%+6.2f%% \n",hiv/norm,lov/norm);
    }
  }

  //un-normalize differential values, for tables
  for (int ibin=0; ibin<nbin; ibin++) {
    //undo pt normalization performed in fitter
    //results[ibin][0][0] *= binwpt[ibin]; //central val
    /// Note:
    ///    re-normalization is done only for xsection central values
    ///    all uncertainties are kept unchanged: they are given in relative values, wrt to also-normalized xsec[]
    ///    hence: all errors stored in result* are normalized as they come from oniafitter
    ///           they are good to give as relative errors (once divided also by normalized xsec)
    ///           but not good to give absolute errors (this would be off by a normalization factor) 
  }

  /// compute cross section sum, and overwrite global fit results
  double xsum(0.), xsumE[2];
  double xglb(0.), xglbE[2];
  for (int k=0; k<2; k++) {xsumE[k]=0;xglbE[k]=0;}
  for (int ibin=0; ibin<nbin; ibin++) {
    if(!ibin) continue; 
    xsum      += results[ibin][0][0];
    xsumE[lo] += pow(results[ibin][1][lo],2);
    xsumE[hi] += pow(results[ibin][1][hi],2);
    //printf("xsum (bin:%d) sum:%6.4f - %6.4f + %6.4f\n",
    //   ibin,results[ibin][0][0],results[ibin][1][lo],results[ibin][1][hi]);
  }
  //printf("xsum:%6.4f - %6.4f + %6.4f\n",xsum, xsumE[hi], xsumE[lo]); 
  xsumE[lo] = sqrt(xsumE[lo]);
  xsumE[hi] = sqrt(xsumE[hi]);
  xglb      = results[0][0][0];
  xglbE[hi] = results[0][1][hi];
  xglbE[lo] = results[0][1][lo];
  printf("xsection total = %5.3f + %4.3f - %4.3f \t sum = %5.3f + %4.3f - %4.3f \t diff=%4.3f (%6.4f%%) (+%4.2f sigma, -%4.2f sigma)\n",
	 results[0][0][0],  results[0][1][hi],  results[0][1][lo],
	 xsum, xsumE[hi], xsumE[lo], 
	 (results[0][0][0]-xsum),(results[0][0][0]-xsum)/xsum*100, 
	 fabs(results[0][0][0]-xsum)/xsumE[lo], fabs(results[0][0][0]-xsum)/xsumE[hi], 
	 );

  //OVERWRITE RESULTS of GLOBAL FIT to full (pt/y-integrated) sample by SUM of partial results
  results[0][0][0]  = xsum;
  results[0][1][hi] = xsumE[hi];  
  results[0][1][lo] = xsumE[lo];

  ///add other sources and compute sum
  for (int ibin=0; ibin<nbin; ibin++) {//row
    //other sources
    for(int k=0;k<2;k++)
      results[ibin][nsys+2][k]=othersys[ibin];//lumi+external
    //total syst
    for(int k=0;k<2;k++)
      for(int isys=0;isys<nsys+1;isys++)
	results[ibin][nsys+3][k] += pow(results[ibin][isys+2][k],2);
    for(int k=0;k<2;k++)
      results[ibin][nsys+3][k] = sqrt(results[ibin][nsys+3][k]);
    //total stat+syst
    for(int k=0;k<2;k++)
      results[ibin][nsys+4][k] = sqrt(pow(results[ibin][1][k],2) + pow(results[ibin][nsys+3][k],2));
  }
    
  //check (in the works)
  if(0)
  for(int i=0; i<nsys; i++) {
    double thres = 0.1;
    bool check = true;
    double vall, valh;
    double min = 100;
    for (int ibin=0; ibin<nbin; ibin++) {
      vall = fabs( results[ibin][i+2][lo]) / xsec[ibin]*100;
      valh = fabs( results[ibin][i+2][hi]) / xsec[ibin]*100;
      if(vall<min) min=val;
      if(valh<min) min=val;
      check &= (vall > thres || valh > thres);
    }
    vall = results[0][i+2][lo] / xsec[0]*100;
    valh = results[0][i+2][hi] / xsec[0]*100;
    if(verbose)
      if(!check) cout << "exclude from table " << i 
		      << "\t" << min 
		      << "\t" << sysnameT[i] 
		      << endl; 
  }


  char ss[100];
  const int ntab = 6;
  TString tabtexName[ntab];
  tabtexName[0] = TString::Format("results_ups%d_y%d.tex",ipeak,irap);	    
  tabtexName[1] = TString::Format("syst_table_ups%d_y%d.tex",ipeak,irap);	    
  tabtexName[2] = TString::Format("summary_ups%d_y%d.tex",ipeak,irap);	    
  tabtexName[3] = TString::Format("polariz_ups%d_y%d.tex",ipeak,irap);	    
  tabtexName[4] = TString::Format("summary_polariz_ups%d_y%d.tex",ipeak,irap);
  tabtexName[5] = TString::Format("global_vs_sum_ups%d_y%d.tex",ipeak,irap);
  //if(rapmode_) 
  // for(int i=0; i<ntab; i++)
  //   tabtexName[i].ReplaceAll(".tex","_rapdiff.tex"); 
  for(int i=0; i<ntab; i++)
    tabtexName[i] = TString::Format("mode%d/",anamode_) + tabtexName[i];
  
  ofstream flatex, tabsys, tabsum, tabpol, sumpol, glbsum;
  flatex.open (tabtexName[0]);
  tabsys.open (tabtexName[1]);
  tabsum.open (tabtexName[2]);
  tabpol.open (tabtexName[3]);
  sumpol.open (tabtexName[4]);
  glbsum.open (tabtexName[5]);

  flatex << std::setw(4)          //occupy n characters
	 << std::setprecision(2); //with after decimal point 
  flatex.fill('0');
  sprintf(ss,"%%Y(%dS) y:%d\n",ipeak,irap);
  flatex << ss;
  tabsys << ss;
  tabsum << ss;
  tabpol << ss;
  sumpol << ss;
  glbsum << ss;

  flatex << "\\begin{tabular}{c";
  tabsys << "\\begin{tabular}{c|c";
  tabsum << "\\begin{tabular}{c|c|c|c|c|c";
  tabpol << "\\begin{tabular}{c|c|c|c|c|c";
  sumpol << "\\begin{tabular}{c c c c c |  c c c c";
  for (int isys=0; isys<nsys+6; isys++)
    flatex << "|c";
  for (int isys=0; isys<nsys; isys++) {
    if(isys>nsystotable-1) continue;
    //if(isys==5) continue; //exclude systematic pt resolution
    tabsys << "|c";
  }
  flatex << "} \\hline\n";
  tabsys << "} \\hline\n";
  tabsum << "} \\hline\n";
  tabpol << "} \\hline\n";
  sumpol << "} \\hline\n";

  //print table (final, after processing)
  for (int ibin=0; ibin<nbin; ibin++) {//row

    if(ibin==0) {
      printf("pt-bin \t   pt-av.  central %4s %4s", "stat-hi","stat-lo");
      if(!rapmode_)
	sprintf(ss,"\\pt (\\GeVc) ");
      else 
	sprintf(ss,"$|y|$");
      tabsum << ss;
      sumpol << ss;
      if(!rapmode_) 
	sprintf(ss,"$p_{\\rm T}$ "); //no space in syst tables to include GeVc
      else
	sprintf(ss,"$|y|$");
      tabsys << ss;
      flatex << ss;
      tabpol << ss;

      //sprintf(ss," & $\\sigma(\\Upsilon)$ & \tstat. &");
      sprintf(ss," & ");
      tabsys << ss;
      sprintf(ss," & $\\sigma(\\Upsilon{\\rm(%dS)})$ & \t%4s &", ipeak, "stat.");
      tabsum << ss;
      sprintf(ss," & $\\sigma$ & \t%4s &", "stat.$\/ \\sigma$");
      sumpol << ss;
      sprintf(ss," & $\\sigma$ & \t%4s &", "stat.");
      flatex << ss;
      for(int isys=0;isys<nsys;isys++) { 
	if(hilo.Contains(sysname[isys]))
	  printf("%15s", sysname[isys].Data());
	else
	  printf("%11s", sysname[isys].Data());
	sprintf(ss,"%11s & ", sysnameT[isys].Data()); 
	flatex<<ss;
	if(isys>nsystotable-1) continue;
	tabsys<<ss;
      }
      printf("%11s %11s %11s\n", "other","tot.syst","tot.error");
      sprintf(ss,"%11s & %11s & %11s & %11s\n", "ext.","$\\sum_{\\rm{syst.}}$","${\\cal{L}}$","$\\Delta\\sigma$ "); 
      flatex<<ss;
      //sprintf(ss,"%11s", "${\\cal{L}}$"); 
      sprintf(ss,"%11s", "add."); 
      tabsys<<ss;
      sprintf(ss,"%11s & %11s & %11s\n", "$\\sum_{\\rm{syst.}}$","lumi.","$\\Delta\\sigma$"); 
      tabsum<<ss;
      sprintf(ss,"%11s & %11s\n", "$\\sum_{\\rm{syst.}}\/\\sigma$","$\\Delta\\sigma\/\\sigma$"); 
      sumpol<<ss;
      //cout << "npol" << npol_ << endl; return;
      TString polnameTT[npol] = {"Unpolarized","Helicity transverse","Helicity longitudinal", "Collins-Soper transverse","Collins-Soper longitudinal"};
      for(int i=0; i<npol_; i++) {
	sprintf(ss,"& %11s ", polnameTT[i].Data());
	tabpol<<ss;
	sprintf(ss,"& %11s ", polnameT[i].Data());
	if(i)
	  sumpol<<ss;
      }
      sprintf(ss,"\\\\ \\hline\n"); 
      flatex << ss;
      tabsys << ss;
      tabsum << ss;
      tabpol << ss;
      sumpol << ss;
    } 

    if(ibin==0) {
      printf("pt-integrat\t ");
      //sprintf(ss," $0 - \\infty$\t");
      if(rapmode_)
	sprintf(ss,"$%4.1f-%4.1f$\t",pte[0],pte[nbin-1]);
      else
	sprintf(ss,"$%4.0f-%4.0f$\t",pte[0],pte[nbin-1]);
      tabsys << ss;
      tabsum << ss;
      sumpol << ss;
      if(rapmode_)
	sprintf(ss,"%4.1f:%4.1f\t",pte[0],pte[nbin-1]);
      else
	sprintf(ss,"%4.0f:%4.0f\t",pte[0],pte[nbin-1]);
      flatex << ss;
      tabpol << ss;
    } else {
      printf("[%4.1f-%4.1f]",pte[ibin-1],pte[ibin]);
      if(rapmode_)
	sprintf(ss,"$%4.1f-%4.1f$\t",pte[ibin-1],pte[ibin]);
      else
	sprintf(ss,"$%4.0f-%4.0f$\t",pte[ibin-1],pte[ibin]);
      tabsys << ss;
      tabsum << ss;
      sumpol << ss;
      if(rapmode_)
	sprintf(ss,"%3.1f:%3.1f\t",pte[ibin-1],pte[ibin]);
      else
	sprintf(ss,"%2.0f:%2.0f\t",pte[ibin-1],pte[ibin]);
      flatex << ss;
      tabpol << ss;
    }

    norm = xsec[ibin]/100.;
    
    //printf("%8.2f %8.1f%% %4.1f%%", results[ibin][0][0], results[ibin][1][0]/norm, results[ibin][1][1]/norm);
    //sprintf(ss,"%4.2f & \t${}_{-%4.1f}^{+%4.1f}$ &", results[ibin][0][0], results[ibin][1][0]/norm, results[ibin][1][1]/norm);
    printf("%8.2f %8.1f%% %4.1f%%", results[ibin][0][0], results[ibin][1][0]/norm, results[ibin][1][1]/norm);
    sprintf(ss,"%4.2f & \t$\\pm{%4.1f}$ &", results[ibin][0][0], results[ibin][1][0]/norm);
    flatex << "&" << ss; 
    if(!ibin)
      sprintf(ss,"%4.2f & \t${%4.1f}$ &", results[ibin][0][0], results[ibin][1][0]/norm);
    else
      sprintf(ss,"%4.2f & \t${%4.0f}$ &", results[ibin][0][0], results[ibin][1][0]/norm);
    sumpol << "&" << ss; 
    tabsys << "&"; 
    sprintf(ss,"%5.3f & \t $\\pm$ %5.3f &",  results[ibin][0][0], results[ibin][1][0]*binwpt[ibin]);
    tabsum << "&" << ss; 

    for(int isys=0;isys<nsys;isys++) {    
      //for(int isys=0;isys<nsystotable;isys++) {    
      //if(isys==4) continue; //exclude systematic pt resolution ATTENTION: hardcoded!!!! (see also above)
      printf(" %+8.1f%%",results[ibin][isys+2][lo]/norm);
      if(!hilo.Contains(sysname[isys])) {
	sprintf(ss," $\\pm%4.1f$ &",results[ibin][isys+2][lo]/norm);
      } else {
	sprintf(ss,"${}_{-%4.1f}^{+%4.1f}$ &",fabs(results[ibin][isys+2][lo]/norm),results[ibin][isys+2][hi]/norm);
      }
      flatex<<ss;
      if(isys>nsystotable-1) continue;
      if(!hilo.Contains(sysname[isys])) {
	sprintf(ss," $%4.1f$ &",results[ibin][isys+2][lo]/norm);
      } else {
	sprintf(ss,"$%4.1f\\,(%4.1f)$ &",fabs(results[ibin][isys+2][hi]/norm),fabs(results[ibin][isys+2][lo])/norm);
      }
      tabsys<<ss;
    }

    //debug:
    //printline(xsec,nbin,"\n\nxsec");
    //printline(external,nbin,"abs");
    //printline(external,nbin,"rel",xsec);
    //return;

    //merge leftover syst columns into last (add.)
    double lastcl = pow(external[ibin],2); //exclude luminosity!!! in syst tables
    //double lastch = pow(external[ibin],2);
    //printf("\nxxxxxxxxx merging soource: %d %d \t",nsystotable,nsys);
    for(int isys=nsystotable;isys<nsys;isys++) {    
      lastcl += pow(results[ibin][isys+2][lo],2);
      //lastch += pow(results[ibin][isys+2][hi],2);
      //printf("%d (%s,%8.6f) ",isys,sysnameT[isys].Data(),results[ibin][isys+2][lo]);
    }
    //printf("\n");
    //return;
    lastcl = sqrt(lastcl)/norm;
    //lastch = sqrt(lastch)/norm;
    //printf("lastc: ibin:%d %8.6f ext:%8.6f\n",ibin,lastcl,external[ibin]/norm);
    sprintf(ss,"$%4.1f$ \\\\ \n", lastcl);
    //sprintf(ss,"$\\pm%4.1f$ \\\\ \\hline\n", lastcl);
    //results[ibin][nsys+2][lo]/norm);//other
    //  note this is done and needed above :results[ibin][nsys+2][k]=othersys[ibin];
    //  ie, the lumi is included in the calculations 
    tabsys<<ss;

    printf("%+8.1f%% "        ,results[ibin][nsys+2][lo]/norm); //other
    printf("%+8.1f%% %+4.1f%%",results[ibin][nsys+3][lo]/norm,results[ibin][nsys+3][hi]/norm); //tot.sys
    printf("%+8.1f%% %+4.1f%%",results[ibin][nsys+4][lo]/norm,results[ibin][nsys+4][hi]/norm); //tot.err
    cout << endl;
    sprintf(ss,"$\\pm$%4.1f & ${}_{-%4.1f}^{+%4.1f}$ & $\\pm%4.0f$ & ${}_{-%4.1f}^{+%4.1f}$ \\\\ \\hline\n",
	    results[ibin][nsys+2][lo]/norm, //other
	    sqrt(pow(results[ibin][nsys+3][lo],2)-pow(xsec[ibin]*lumirelerr,2))/norm, //tot.sys w/o lumi (lo)
	    sqrt(pow(results[ibin][nsys+3][hi],2)-pow(xsec[ibin]*lumirelerr,2))/norm, //tot.sys w/o lumi (hi)
	    //results[ibin][nsys+3][lo]/norm,results[ibin][nsys+3][hi]/norm,//tot.sys
	    lumirelerr*100, 
	    results[ibin][nsys+4][lo]/norm,results[ibin][nsys+4][hi]/norm); //tot.err
    flatex<<ss;

    double afac = results[ibin][0][0]/xsec[ibin];  /* affects integrated fit, when giving absolute unc values,
						      integtated value extracted from sum instead of fit;
						      also factor of two as fit includes all (affects mode1)
						   */

    sprintf(ss,"${}_{-%5.3f}^{+%5.3f}$ & $\\pm$ %5.3f & ${}_{-%5.3f}^{+%5.3f}$ \\\\ \\hline\n",
	    //results[ibin][nsys+3][lo],results[ibin][nsys+3][hi],//tot.sys
	    //sqrt(pow(results[ibin][nsys+3][lo],2)-pow(xsec[ibin]*lumirelerr,2))*binwpt[ibin], //tot.sys w/o lumi (lo)
	    //sqrt(pow(results[ibin][nsys+3][hi],2)-pow(xsec[ibin]*lumirelerr,2))*binwpt[ibin], //tot.sys w/o lumi (hi)
	    //xsec[ibin]*lumirelerr*binwpt[ibin], //lumi
	    sqrt(pow(results[ibin][nsys+3][lo],2)-pow(xsec[ibin]*lumirelerr,2))*binwpt[ibin]*afac, //tot.sys w/o lumi (lo)
	    sqrt(pow(results[ibin][nsys+3][hi],2)-pow(xsec[ibin]*lumirelerr,2))*binwpt[ibin]*afac, //tot.sys w/o lumi (hi)
	    xsec[ibin]*lumirelerr*binwpt[ibin]*afac, //lumi
	    results[ibin][nsys+4][lo]*binwpt[ibin]*afac,results[ibin][nsys+4][hi]*binwpt[ibin]*afac); //tot.err
    tabsum<<ss;

    sprintf(ss,"$%4.0f\\,(%4.0f)$ & $%4.0f\\,(%4.0f)$",
	    //sprintf(ss,"${}_{-%4.0f}^{+%4.0f}$ & ${}_{-%4.0f}^{+%4.0f}$",
	    //sprintf(ss,"${}_{-%4.0f}^{+%4.0f}$ & ${}_{-%4.0f}^{+%4.0f}$",
	    sqrt(pow(results[ibin][nsys+3][hi]/norm,2)-pow(lumirelerr*100,2)),  //tot.sys w/o systematics! (hi)
	    sqrt(pow(results[ibin][nsys+3][lo]/norm,2)-pow(lumirelerr*100,2)),  //tot.sys w/o systematics! (lo)
	    //results[ibin][nsys+3][lo]/norm,results[ibin][nsys+3][hi]/norm,//tot.sys
	    results[ibin][nsys+4][hi]/norm,results[ibin][nsys+4][lo]/norm); //tot.err    
    sumpol<<ss;

    for(ipol = 0; ipol<npol_; ipol++) {
      sprintf(ss," & $ %6.2f \\pm %5.2f {}_{-%4.2f}^{+%4.2f} \\pm%4.2f$",
	      polarization[ibin][ipol][0]*binwpt[ibin], //central
	      polarization[ibin][ipol][1]*binwpt[ibin], //stat.
	      //syst (Note simplification taken here: from unpolarized directly)
	      sqrt(pow(results[ibin][nsys+3][lo]/norm,2)-pow(lumirelerr*100,2))/100.*polarization[ibin][ipol][0]*binwpt[ibin],  //tot.sys w/o systematics! (lo)
	      sqrt(pow(results[ibin][nsys+3][hi]/norm,2)-pow(lumirelerr*100,2))/100.*polarization[ibin][ipol][0]*binwpt[ibin],  //tot.sys w/o systematics! (hi)

	      lumirelerr*polarization[ibin][ipol][0]*binwpt[ibin]//lumi
	      );
      tabpol<<ss;
      if(ipol) {
	sprintf(ss," & %+5.0f",
		(polarization[ibin][ipol][0]-polarization[ibin][0][0])/polarization[ibin][0][0]*100);
	sumpol<<ss;
      }
    }
   sprintf(ss,"\\\\ \\hline\n");
   tabpol<<ss;
   sprintf(ss,"\\\\ \n");
   sumpol<<ss;
  }

  sprintf(ss,"\\hline\n");
  sumpol<<ss;
  tabsys<<ss;

  int kbin;
  double* rapb = getrapbinedges("nominal",kbin);
  double* ptbb = getptbinedges("nominal",kbin,true);

  //printline(rapb,nbin,"debug rap edges");
  sprintf(ss,"$\\Upsilon$(%dS)\\,, ",ipeak);
  if      (anamode_==0)
    sprintf(ss,"%s $%2.0f<|y|<%2.0f\\;:\\quad$ ",ss, rapb[0],rapb[1]);
  else if (anamode_==1) 
    sprintf(ss,"%s $%2.0f<|y|<%2.0f\\;:\\quad$ ",ss, rapb[irap],rapb[irap+1]);
  else if (anamode_==2) 
    sprintf(ss,"%s $\\pt<%3.0f\\quad\\;:\\quad$ ",ss, ptbb[kbin-1]);
  //cout << ss << " " << nbin << endl; return;
  glbsum << ss;
  sprintf(ss," & $\\sigma = %5.3f \\pm %5.3f\\;,$", xglb, xglbE[0]);
  glbsum << ss;
  sprintf(ss," & $\\sum d\\sigma = %5.3f \\pm %5.3f\\;,$ ", xsum, xsumE[0]);
  glbsum << ss;
  sprintf(ss," & $\\Delta = %5.1f \\%%$ ", fabs(xglb-xsum)/xglb*100);
  glbsum << ss;
  sprintf(ss," \\\\ %% %ds y:%d mode:%d \n",ipeak,irap,anamode_);
  glbsum << ss;

  sprintf(ss,"\\end{tabular}\n");
  flatex << ss;
  tabsys << ss;
  tabsum << ss;
  sumpol << ss;
  tabpol << ss;
  flatex.close();
  tabsys.close();
  tabsum.close();
  tabpol.close();
  sumpol.close();
  glbsum.close();


  //return;//
  //tot unc xsec graph
  double ptl[nbin], pth[nbin];
  double xsecdiff[nbin];
  double sysButLumihi[nbin], sysButLumilo[nbin];
  double ptnull[nbin], ptesys[nbin];

  //double *pte = getptbinedges(peak,nbin_);  // pt/rap bin edges -- ?? this changed value after being defined above
  printline(pte,nbin,"debug pt edges");
  printline(pt,nbin,"debug pt center");
  //return;

  for(int ibin=0; ibin<nbin; ibin++) { 
    
    xsecdiff[ibin]  = xsec[ibin];
    xsecdiff[ibin] /= binw[ibin]; 
    totunclo[ibin] /= binw[ibin]; 
    totunchi[ibin] /= binw[ibin]; 
    xsecEh[ibin]   /= binw[ibin];
    xsecEl[ibin]   /= binw[ibin]; 
    sysButLumilo[ibin]=sqrt(pow(totunclo[ibin],2)-pow(xsecdiff[ibin]*lumirelerr,2));
    sysButLumihi[ibin]=sqrt(pow(totunchi[ibin],2)-pow(xsecdiff[ibin]*lumirelerr,2));
    
    if(rapmode_ && ibin) 
      assert((binw[ibin]-0.4)/0.4<0.000001); 

    //ptl[ibin]=ptb[ibin];
    //pth[ibin]=ptb[ibin];
    if(ibin==0)
      ptl[ibin] = pte[0];
    else {
      ptl[ibin]=pt[ibin]-pte[ibin-1];
      pth[ibin]=pte[ibin]-pt[ibin];
    }
    ptnull[ibin]=0.15;//ptb[ibin];    
    if(rapmode_)
      ptnull[ibin]=0.015;    
    
    //Divide all plots by factor of 2, in case of choice of y-differential, instead of |y|-differential
    if(!absrap && ibin) {
      xsecdiff[ibin]*=0.5;
      totunclo[ibin]*=0.5;
      totunchi[ibin]*=0.5;
      xsecEl  [ibin]*=0.5;
      xsecEh  [ibin]*=0.5;
      sysButLumilo[ibin]*=0.5;
      sysButLumihi[ibin]*=0.5;
    }

    ptesys[ibin]=0.;
  }

//  for(int ibin=0; ibin<nbin; ibin++) { 
//    //pt-diff: already normalized in fitter (both central values and uncertainties!)
//    xsecdiff[ibin] = xsec[ibin]; 
//    //printf("ptbin: %d pt:%4.1f ptl:%4.1f pth%4.1f \n", ibin, pt[ibin],pte[ibin],pte[ibin]);
//    if(rapmode_ && ibin) {
//      //rap-diff: already un-normalize by pt above (in get results/getunceertainties), normalize now by rapidity
//      xsecdiff[ibin] /= (pte[ibin]-pte[ibin-1]); 
//      totunclo[ibin] /= (pte[ibin]-pte[ibin-1]); 
//      totunchi[ibin] /= (pte[ibin]-pte[ibin-1]); 
//      assert(((pte[ibin]-pte[ibin-1])-0.4)/0.4<0.000001); 
//      //ptl[ibin]=pt [ibin]-ptb[ibin-1];
//      //pth[ibin]=ptb[ibin]-pt[ibin];
//    } else {
//      //ptl[ibin]=pt [ibin]-ptb[ibin-1];
//      //pth[ibin]=ptb[ibin]-pt[ibin];
//      //ptl[ibin]=0; pth[ibin]=0;
//    }
//    ptl[ibin]=ptb[ibin];
//    pth[ibin]=ptb[ibin];
//  }


  //return;//
  TGraphAsymmErrors* xsecgr = new TGraphAsymmErrors(nbin, pt, xsecdiff, ptl, pth, totunclo, totunchi);
  //TGraphAsymmErrors* xsecgr    = new TGraphAsymmErrors(nbin, pt, xsecdiff, ptnull, ptnull, totunclo, totunchi);

  TGraphAsymmErrors* xsecgrsta = new TGraphAsymmErrors(nbin, pt, xsecdiff, ptnull, ptnull, xsecEl, xsecEh);

  TGraphAsymmErrors* xsecgrsys = new TGraphAsymmErrors(nbin, pt, xsecdiff, ptesys, ptesys, sysButLumilo, sysButLumihi);

  /*
  //stat only unc xsec graph  
  TString fnvar = getFileName(peak);
  TFile fvar; fvar.Open(fnvar,"READ");
  gDirectory->Cd(fnvar+":/plots");
  TGraphAsymmErrors* ptdiff = (TGraphAsymmErrors*)gROOT->FindObject(TString::Format("xsection_ups%dS_rap%d_pt", ipeak_, irap_));
  fvar.Close();
  */
  
  if(rapmode_) 
    //cross section graph, rapidity differential
    plotxsecgrRap(xsecgr,xsecgrsta,xsecgrsys);
  else 
    //cross section graph, pt differential
    plotxsecgr(xsecgr,xsecgrsta,xsecgrsys);
  
  
  if(pol)
    //cross section polarization graph
    plotpolgr(polariz,polname);
  
}


TString getFileName(TString name) {
  return TString::Format("mode%d/fitres_%ds_upsilonYieldWeighted_%s/xsection.root",anamode_,ipeak_,name.Data());
  //TString pk = TString::Format("_%ds",ipeak_);
  //return "mode" + str(anamode_) + "/fitres" + pk  + "_upsilonYieldWeighted_" + name + "/xsection.root";
}

double* getptbinsRap(TString filevar, const int nbin, int center) {
  double *diff = new double[nbin];
  TString fnvar = getFileName(filevar);
  TFile fvar; fvar.Open(fnvar,"READ");
  gDirectory->Cd(fnvar+":/plots");
  difvar = (TH2F*)gROOT->FindObject(TString::Format("pt_y_bins"));
  if(center==0)
    diff[0]=-1;
  else
    diff[0]=0;
  for(int i=0; i<nbin-1; i++) {
    if       (center==0) {
      diff[i+1] = difvar->GetYaxis()->GetBinCenter(i+1);
    } else if(center==1) {
      diff[i+1] = difvar->GetYaxis()->GetBinCenter(i+1) - difvar->GetYaxis()->GetXbins()->GetArray()[i];
    } else if(center==2) {
      diff[i+1] = difvar->GetYaxis()->GetXbins()->GetArray()[i+1] - difvar->GetYaxis()->GetBinCenter(i+1);
    }
  }
  //diff[i] = UpsRapBinEdges[i];
  //printf("aa %d  %d\n", difvar->GetYaxis()->GetNbins(), nbin);
  //for(int i=0; i<nbin; i++) 
  //printf("zz i:%d %f\n",i,difvar->GetYaxis()->GetXbins()->GetArray()[i]);
  if(verbose) {
    for(int i=0; i<nbin; i++)
      printf("%4.1f - ",diff[i]);
    cout << endl;
  }
  fvar.Close();
  return diff;
}

double* getptbinedges(TString filevar, int &nbin, bool ispt) {
  TString fnvar = getFileName(filevar);
  //TFile* fvar = TFile::Open(fnvar,"READ");//this gives some incosistencies in read out values when closing)
  TFile fvar; fvar.Open(fnvar,"READ");
  gDirectory->Cd(fnvar+":/plots");
  difvar = (TH2F*)gROOT->FindObject(TString::Format("pt_y_bins"));
  double* diff;
  if(rapmode_ && !ispt) {
    nbin = difvar->GetYaxis()->GetNbins() + 1;
    diff = difvar->GetYaxis()->GetXbins()->GetArray();
  } else {
    nbin = difvar->GetXaxis()->GetNbins() + 1;
    diff = difvar->GetXaxis()->GetXbins()->GetArray();
  }
  fvar.Close();
  //for(int i=0;i<nbin;i++)  printf("edges [%d]:%3.1f\t",i,diff[i]); printf("\nnbin:%d\n",nbin);
  return diff;
}


double* getrapbinedges(TString filevar, int &nbin) {
  TString fnvar = getFileName(filevar);
  TFile fvar; fvar.Open(fnvar,"READ");
  gDirectory->Cd(fnvar+":/plots");
  difvar = (TH2F*)gROOT->FindObject(TString::Format("pt_y_bins"));
  double* diff;
  nbin = difvar->GetYaxis()->GetNbins() + 1;
  diff = difvar->GetYaxis()->GetXbins()->GetArray();
  //for(int i=0;i<nbin+1;i++)  printf("%3.1f\t",diff[i]); printf("%f\n",nbin);
  fvar.Close();
  return diff;
}


double* getptbins(TString filevar, int &nbin,int center) {
  if(rapmode_) {
    //get n rap bin
    //tbd:get nbins from filevar (just rerun the fit)
    //TString fnvar = getFileName(filevar);
    TString fnvar = getFileName("nominal");
    TFile fvar; fvar.Open(fnvar,"READ");
    gDirectory->Cd(fnvar+":/plots");
    difvarr = (TH2F*)gROOT->FindObject(TString::Format("pt_y_bins"));
    nbin = difvarr->GetYaxis()->GetNbins() +1;

    double* diff = getptbinsRap(filevar, nbin, center);
    fvar.Close();
    return diff;
  }
  TString fnvar = getFileName(filevar);
  TFile fvar; fvar.Open(fnvar,"READ");
  TGraphAsymmErrors *difvar;
  gDirectory->Cd(fnvar+":/plots");
  difvar = (TGraphAsymmErrors*)gROOT->FindObject(TString::Format("xsection_ups%dS_rap%d_pt", ipeak_, irap_));
  const int nb = difvar->GetN() +1;
  nbin=nb;
  double *cen = difvar->GetX();
  double *lef = difvar->GetEXlow();
  double *rig = difvar->GetEXhigh();
  double *diff = new double[nb];
  if(center==0) 
    diff[0]=-1;
  else
    diff[0] = 0.;
  for(int i=0; i<nb-1; i++) {
    if(center==0) 
      diff[i+1] = cen[i];
    else if(center==1)//lo
      diff[i+1] = lef[i];
    else if(center==2)//hi
      diff[i+1] = rig[i];
  }
  if(verbose) {
    for(int i=0; i<nb; i++)
      printf("%4.1f - ",diff[i]);
    cout << endl;
  }
  fvar.Close();
  return diff;
}

double* getUncertainty(TString filevar, int &nbin, bool hi){
  if(rapmode_) {
    double* diff = getUncertaintyRap(filevar, nbin, hi);
    return diff;
  }

  TString fnvar = getFileName(filevar);
  TFile fvar; fvar.Open(fnvar,"READ");
  TGraphAsymmErrors *difvar, *totvar; 
  gDirectory->Cd(fnvar+":/plots");
  difvar = (TGraphAsymmErrors*)gROOT->FindObject(TString::Format("xsection_ups%dS_rap%d_pt", ipeak_, irap_));
  totvar = (TGraphAsymmErrors*)gROOT->FindObject(TString::Format("xsection_ups%dS_total"   , ipeak_));
  const int nb = difvar->GetN() + 1;
  nbin = nb;
  double *diff = new double[nb];
  if(hi) 
    diff[0] = (totvar->GetEYhigh())[0];
  else
    diff[0] = (totvar->GetEYlow())[0];
  for(int i=0; i<nb-1; i++) {
    if(hi)
      diff[i+1] = (difvar->GetEYhigh())[i];
    else
      diff[i+1] = (difvar->GetEYlow())[i];
  }

  //pt binning un-normalization
  int nn;
  double* ptbb = getptbinedges("nominal",nn,true);
  assert(nn==nbin);
  //printline(ptbb,nn,"debug");
  for(int i=0; i<nbin-1; i++) 
    diff[i+1] *= (ptbb[i+1]-ptbb[i]);

  fvar.Close();
  return diff;
}


double* getUncertaintyRap(TString filevar, const int nbin, bool hi){
  TString fnvar = getFileName(filevar);
  TFile fvar; fvar.Open(fnvar,"READ");
  TGraphAsymmErrors *difvar[nbin], *totvar; 
  gDirectory->Cd(fnvar+":/plots");
  totvar = (TGraphAsymmErrors*)gROOT->FindObject(TString::Format("xsection_ups%dS_total"   , ipeak_));
  for(int i=0; i<nbin; i++) {
    difvar[i] = (TGraphAsymmErrors*)gROOT->FindObject(TString::Format("xsection_ups%dS_rap%d_pt", ipeak_, i));
  }
  //const int nb = difvar->GetN() + 1;
  //nbin = nb;
  double *diff = new double[nbin];
  if(hi) 
    diff[0] = (totvar->GetEYhigh())[0];
  else
    diff[0] = (totvar->GetEYlow())[0];
  for(int i=0; i<nbin-1; i++) {
    if(hi)
      diff[i+1] = (difvar[i]->GetEYhigh())[0];
    else
      diff[i+1] = (difvar[i]->GetEYlow())[0];
  }

  //pt un-normalization
  int nb;
  double* ptbb = getptbinedges("nominal",nb,true);
  assert(anamode_==2 && nb==2);
  //printline(ptbb,nb,"debug");
  for(int i=0; i<nbin-1; i++) 
    diff[i+1] *= (ptbb[1]-ptbb[0]);
  //for(int i=0; i<nbin-1; i++) 
  // diff[i+1] *= ptwidth_;

  if(verbose) {
    printf("dg::%s\t",filevar.Data());
    for(int i=0; i<nbin; i++) 
      printf("(%d) %6.4f\t",i,diff[i]);
    cout << "\t(unc., hi:" << hi << ")" << endl;
  }
  fvar.Close();
  return diff;
}

double* getResultsRap(TString filevar, const int nbin){
  //nbin=nbinn+1;
  TString fnvar = getFileName(filevar);
  TFile fvar; fvar.Open(fnvar,"READ");
  gDirectory->Cd(fnvar+":/plots");
  totvar = (TGraphAsymmErrors*)gROOT->FindObject(TString::Format("xsection_ups%dS_total"   , ipeak_));
  TGraphAsymmErrors *difvar[nbin], *totvar; 
  for(int i=0; i<nbin; i++) {
    difvar[i] = (TGraphAsymmErrors*)gROOT->FindObject(TString::Format("xsection_ups%dS_rap%d_pt", ipeak_, i));
  }
  //const int nb = difvar->GetN() + 1;
  //nbin = nb;`
  double *diff = new double[nbin];
  diff[0] = (totvar->GetY())[0];
  for(int i=0; i<nbin-1; i++) 
    diff[i+1] = (difvar[i]->GetY())[0];

  //pt un-normalization
  int nb;
  double* ptbb = getptbinedges("nominal",nb,true);
  assert(anamode_==2 && nb==2);
  //printline(ptbb,nb,"debug");
  for(int i=0; i<nbin-1; i++) 
    diff[i+1] *= (ptbb[1]-ptbb[0]);

  if(verbose) {
    printf("dg::%s\t",filevar.Data());
    for(int i=0; i<nbin; i++) 
      printf("(%d) %6.4f\t",i,diff[i]);
    cout << endl;
  }
  fvar.Close();
  return diff;
}

double* getResults(TString filevar, int &nbin){
  if(rapmode_) {
    //get n rap bin
    //tbd:get nbins from filevar (just rerun the fit)
    //TString fnvar = getFileName(filevar);
    TString fnvar = getFileName("nominal");
    TFile fvar; fvar.Open(fnvar,"READ");
    gDirectory->Cd(fnvar+":/plots");
    difvarr = (TH2F*)gROOT->FindObject(TString::Format("pt_y_bins"));
    nbin = difvarr->GetYaxis()->GetNbins() +1;

    double* diff = getResultsRap(filevar, nbin);
    fvar.Close();
    return diff;
  }
  TString fnvar = getFileName(filevar);
  TFile fvar; fvar.Open(fnvar,"READ");
  TGraphAsymmErrors *difvar, *totvar; 
  gDirectory->Cd(fnvar+":/plots");
  //cout << gDirectory->pwd() << gDirectory->ls() << endl;
  difvar = (TGraphAsymmErrors*)gROOT->FindObject(TString::Format("xsection_ups%dS_rap%d_pt", ipeak_, irap_));
  totvar = (TGraphAsymmErrors*)gROOT->FindObject(TString::Format("xsection_ups%dS_total"   , ipeak_));
  const int nb = difvar->GetN() + 1;
  nbin = nb;
  double *diff = new double[nb];
  diff[0] = (totvar->GetY())[0];
  for(int i=0; i<nbin-1; i++)
    diff[i+1] = (difvar->GetY())[i];

  //pt binning un-normalization
  int nn;
  double* ptbb = getptbinedges("nominal",nn,true);
  assert(nn==nbin);
  //printline(ptbb,nn,"debug");
  for(int i=0; i<nbin-1; i++) 
    diff[i+1] *= (ptbb[i+1]-ptbb[i]);

  fvar.Close();
  return diff;
}

double* getVariation(TString filevar, TString fileref){

  int nbvar, nbref;
  double *sysv, *sysr; 
  sysv = getResults(filevar, nbvar);
  sysr = getResults(fileref, nbref);
  assert(nbvar==nbref);
  double *diff = new double[nbref];
  for(int i=0; i<nbref; i++) {
    diff[i] = sysv[i] - sysr[i];
  }

  if(verbose) {
    printf("var:%s ref:%s\n",filevar.Data(), fileref.Data()); 
    for(int i=0; i<nbref; i++) {
      printf("\t%d var:%7.5f ref:%7.5f diff:%7.5f\n",i, sysv[i],sysr[i],diff[i]);
    }
  }

  return diff;
}

void printline(double* diff, int nbin, TString n, double* norm) {
  bool donorm = true;
  printf("%-10s",n.Data());
  for(int j=0; j<nbin; j++) {
    if(donorm && norm!=NULL) {
      printf("\t%+5.1f%%", diff[j] / norm[j] * 100);
      //      printf("\t%5.1+f%%", diff[j] / norm[j] * 100);
    } else {
      printf("\t%5.2f",diff[j]); 
    }
  }
  cout << endl;
}

void printptline(TString fvar) {
  int nbin;
  double *ptt = getptbinedges(fvar,nbin);
  //return; //xxxx
  //double *ptb = getptbins(fvar,nbin,0);
  printf("%17s %2.0f-%-2.0f "," ",ptt[0],ptt[nbin-1]);

  for(int j=0; j<nbin-1; j++) {
  if(!rapmode_)
    printf("%2.0f-%-2.0f   ",ptt[j],ptt[j+1]);
  else
    printf("%3.1f-%-3.1f   ",ptt[j],ptt[j+1]);
  }
  cout << endl;
}

void drawrule(int nbin) {
  for(int i=0; i<15; i++) printf("_");
  for(int j=0; j<nbin; j++) {
    //printf("|");
    for(int i=0; i<8; i++) printf("_");
  }
  cout << endl;
}


void plotpolgr(vector <TString>  polariz,  vector <TString>  polname) {
  
  assert((polariz.size())==5);
  assert((polname.size())==5);

  int nbinn;
  double *pt     = getptbins(polariz[0],nbinn,0);
  double *pterr  = getptbins(polariz[0],nbinn,1);
  //for(int j=0; j<nbinn; j++) printf("%d %f\n",j, pt[j]); return;
  //tot unc xsec graph
  const int nbin=nbinn;
  double ptl[nbin], pth[nbin],xerr[nbin];
  for(int j=0; j<nbin; j++) { 
    ptl[j]=0; pth[j]=0;
    xerr[j]=0.;
  }
  
  const int ngr = polariz.size();
  double *pol[ngr];

  TGraphAsymmErrors polgr[ngr];

  bool first = true;
  double *pte = getptbinedges("nominal",nbin);  // pt/rap bin edges
  //printline(pte,nbin,"debug pt edges"); return;
  double binw[nbin];
  for(int ibin=0; ibin<nbin; ibin++) {
    if(ibin==0)  binw[ibin]=1.;
    else         binw[ibin]=pte[ibin]-pte[ibin-1];
  }
  //printline(binw,nbin,"debug binw"); return;

  for(int i=0; i<ngr; i++) {
    int nn=nbin;  
    pol  [i] = getResults(polariz[i],nn); 

    //renomalize by rapidity with
    for(int ibin=0; ibin<nbin; ibin++) { 
      //pt-diff: already normalized in fitter 

      (pol[i])[ibin] /= binw[ibin];
      
      //printf("dddddd %f\n",pte[ibin]-pte[ibin-1]);
      if(rapmode_ && ibin) {	assert(((pte[ibin]-pte[ibin-1])-0.4)/0.4<0.000001);} 
      
    }

    //shift all values -- for avoiding distoting the interpolation lines
    if(diffonly) {
      nbinn = nbin-1;
      //(pol[0])[0]=0; 
      for(int ibin=0; ibin<nbinn; ibin++) { 
	if(first) {
	  pt [ibin]   = pt [ibin+1];
	  ptl[ibin]   = ptl[ibin+1];
	  pth[ibin]   = pth[ibin+1];
	}
	xerr[ibin]   = xerr[ibin+1];
	pol[i][ibin] = pol[i][ibin+1];
      }
      first = false;
    }
    polgr[i] = TGraphAsymmErrors(nbin, pt, pol[i], ptl,pth,xerr,xerr);
  }

  TCanvas b;
  TMultiGraph *mg = new TMultiGraph();
  
  TLegend *leg = new TLegend(0.7,0.5,0.9,0.7);	
  for(int i=0; i<ngr; i++) {
    //polgr[i].SetMarkerSize(0.6);
    polgr[i].SetLineWidth  (0.01);
    polgr[i].SetLineColor  (2+10*i);
    polgr[i].SetMarkerColor(2+10*i);
    polgr[i].SetMarkerStyle(20+i);
    mg->Add(&polgr[i]);
    leg->AddEntry(&polgr[i],polname[i],"pl");
  }

  mg->Draw("apc");//apf/a4
  
  TString ytitle;
  if(rapmode_) {
    ytitle = TString::Format("d#sigma/d|y| * #Beta(#Upsilon(%dS)#rightarrow#mu#mu) (nb)",ipeak_);
    mg->GetXaxis()->SetTitle("|y^{#Upsilon}|");  
  } else {
    ytitle = TString::Format("d#sigma/dp_{T}* #Beta(#Upsilon(%dS)#rightarrow#mu#mu)  (nb/(GeV/c))",ipeak_);
    mg->GetXaxis()->SetTitle("p_{T}^{#Upsilon} (GeV/c)");
    gPad->SetLogy();
  }
  mg->GetYaxis()->SetTitleOffset(1.0); //needed
  mg->GetYaxis()->SetTitle(ytitle);

  
  TPaveText *tp = new TPaveText(0.55,0.75,0.9,0.85,"brNDC");
  tp->SetBorderSize(0);
  tp->AddText("CMS, #sqrt{s} = 7 TeV");
  double* rapedge = getrapbinedges("nominal",nbin);
  TString lablumi = TString::Format("L = %2.0f pb^{-1}", lumiforlabel);
  TString labrap = rapmode_ ? "" : TString::Format(", |y^{#Upsilon}|<%2.0f",rapedge[1]);
  //TString labrap;
  //if(!rapmode_) labrap = TString::Format("|y^{#Upsilon}|<%2.0f",rapedge[irap_+1]);
  //else         labrap = TString::Format("|y^{#Upsilon}|<%2.0f",rapedge[irap_+1]);
  //if(rapedge[irap_]!=0) labrap = TString::Format("%2.0f<",rapedge[irap_]) + labrap;
  tp->AddText(TString::Format("%s %s",lablumi.Data(),labrap.Data())); 

  tp->Draw();
  leg->Draw();

  double maxy = ymax_;
  if(rapmode_)
    if      (ipeak_==1) maxy = 8;
    else if (ipeak_==2) maxy = 4;
    else if (ipeak_==3) maxy = 2;
  mg->SetMaximum(maxy);

  mg->SetMinimum(ymin_);
  if(diffonly)
    mg->GetXaxis()->SetRangeUser(0,ptmax_);
  else
    mg->GetXaxis()->SetRangeUser(-2,ptmax_);

  TString pname; 
  //if(rapmode_) pname = TString::Format("mode%d/polgr_%ds_ydiff", anamode_,ipeak_);
  //else         
  pname = TString::Format("mode%d/polgr_%ds_y%d", anamode_,ipeak_,irap_);

  b.SaveAs(pname+".pdf");
  b.SaveAs(pname+".gif");
  return;
}


void plotxsecgrRap(TGraphAsymmErrors* xsecgr, TGraphAsymmErrors* xsecgrsta, TGraphAsymmErrors* xsecgrsys) {

  gROOT->LoadMacro("setTDRStyle_modified.C");
  setTDRStyle();


  TCanvas a;
  TMultiGraph *mgr = new TMultiGraph();

  mgr->Add(xsecgr);
  //mgr->Add(xsecgrsys);
  //mgr->Add(xsecgrsta);
  mgr->Draw("ap");


  //TString ytitle = TString::Format("d#sigma/d|y| * #Beta(#Upsilon(%dS)#rightarrow#mu#mu) (nb)",ipeak_);
  TString ytitle = TString::Format("d#sigma/dy #times #Beta(#mu#mu)  (nb)");
  TString xtitle("y^{#Upsilon}");
  mgr->GetYaxis()->SetTitleOffset(1.0); //needed
  mgr->GetYaxis()->SetTitle(ytitle);
  mgr->GetXaxis()->SetTitle(xtitle);
  
  double absrapfact = absrap?1:0.5;

  double maxy = ymax_;
  if(rapmode_) {
    if      (ipeak_==1) maxy = 8;
    else if (ipeak_==2) maxy = 4;
    else if (ipeak_==3) maxy = 2;
    maxy *= absrapfact*1.2;
  }
  mgr->SetMaximum(maxy);
  

  mgr->SetMinimum(0);

  if(diffonly)
    mgr->GetXaxis()->SetRangeUser(0,2);
  else
    mgr->GetXaxis()->SetRangeUser(-2,2);

  TPaveText *tp = new TPaveText(0.55,0.72,0.9,0.9,"brNDC");
  //TPaveText *tp = new TPaveText(0.55,0.75,0.9,0.85,"brNDC");
  tp->SetBorderSize(0);
  tp->AddText("CMS, #sqrt{s} = 7 TeV");
  tp->AddText(TString::Format("L =%3.0f pb^{-1}",lumiforlabel));
  tp->AddText("");
  tp->AddText(TString::Format("#Upsilon(%dS)",ipeak_)); 
  tp->Draw();
  
  //TString pname = TString::Format("mode%d/xsec_%ds_rapdiff", anamode_,ipeak_);
  TString pname = TString::Format("mode%d/xsection_%ds_y0", anamode_,ipeak_);
  a.SaveAs(pname + ".gif");
  a.SaveAs(pname + ".pdf");
  TFile o (pname + ".root","recreate");  
  mgr->Write("mg1");  
  o.Close();
  
  TCanvas c;
  TH1F *hh = gPad->DrawFrame(0,0,2,maxy);
  hh->GetYaxis()->SetTitleOffset(1.0); //needed
  hh->SetXTitle(xtitle);
  hh->SetYTitle(ytitle);
  tp->Draw();

  double dummyx[1] = {0.}, dummyy[1] = {0.};
  dummy = new TGraph(1, dummyx, dummyy);

  bool theory = true;

  TLegend *leg; 
  if(!theory) 
  leg = new TLegend(0.15,0.7,0.55,0.9);	
  else
    leg = new TLegend(0.15,0.7,0.57,0.88);	
  //leg = new TLegend(0.15,0.65,0.55,0.9);	
  leg->AddEntry(xsecgr,   "total unc.","ple");
  //leg->AddEntry(xsecgrsys,"stat. #oplus syst. unc.","p");
  //leg->AddEntry(dummy,"global lumi unc. (11%)","e0");
  leg->AddEntry(xsecgrsys,"total unc. except lumi.","p");
  //leg->AddEntry(dummy,"luminosity (11%)","e0");
  leg->AddEntry(xsecgrsta,"statistical unc.","p");


  xsecgrsta->SetLineColor  (38);//(kBlue);
  xsecgrsta->SetMarkerColor(38);//(kBlue);
  xsecgrsta->SetFillColor  (38);//(kBlue);
  xsecgrsta->SetMarkerStyle(21);
  xsecgrsta->SetMarkerSize(0.8);
  //xsecgrsta->SetFillStyle(0);


  xsecgrsys->SetLineColor  (kRed);
  xsecgrsys->SetMarkerColor(kRed);
  xsecgrsys->SetMarkerStyle(25);
  xsecgrsys->SetMarkerSize(0.8);
  xsecgrsys->SetFillStyle(0);
  
  xsecgr   ->Draw("e0");
  xsecgrsta->Draw("e2");
  xsecgrsys->Draw("[]");

  double pythia_norm[3]; //1s/2s/3s
  pythia_norm[0] = 7.37/15.17;
  pythia_norm[1] = 2.16/6.47;
  pythia_norm[2] = 1.02/1.7;

  if(theory) {
    TFile* th = new TFile(TString::Format("../GenCrossSection/genXSupsilon%dSrap.root",ipeak_),"read"); // open the file
    TH1F* pythia = (TH1F*)gDirectory->Get("genY"); // get the hist
    //pythia->Rebin(2); pythia->Scale(0.5);
    pythia->Scale(pythia_norm[ipeak_-1]*absrapfact);
    pythia->SetLineColor(kBlue);
    pythia->SetMarkerColor(kBlue);
    xsecgrsys->SetMarkerSize(1.2);
    xsecgrsys->SetLineWidth(2);
    pythia->SetLineWidth(2);
    pythia->SetLineStyle(2);
    pythia->SetMarkerStyle(26);

    //TLegend *tleg; 
    //tleg = new TLegend(0.15,0.2,0.55,0.27);	
    //tleg->AddEntry(pythia,"Pythia (normalized)","l");
    //tleg->Draw();    

    leg->AddEntry(pythia,"PYTHIA (normalized)","l");
    pythia->Draw("same HIST C");
  }

  leg->Draw();


  c.SaveAs(pname + "_detailed.gif");
  c.SaveAs(pname + "_detailed.pdf");

  th->Close();

  return;
}

void plotxsecgr(TGraphAsymmErrors* xsecgr, TGraphAsymmErrors* xsecgrsta, TGraphAsymmErrors* xsecgrsys) { 

  gROOT->LoadMacro("setTDRStyle_modified.C");
  setTDRStyle();

  TCanvas a;
  TMultiGraph *mgr = new TMultiGraph();

  //ptdiff->SetLineColor(2);
  //ptdiff->SetLineWidth(1);
  //xsecgr->SetLineWidth(3);
  //xsecgr->SetLineColor(kBlue);
  //xsecgr->SetMarkerSize(0.7);
  //xsecgr->SetMarkerStyle(4);
  mgr->Add(xsecgr);

  //if(!diffonly)
  //mgr->Add(ptdiff);
  mgr->Draw("ap");
  gPad->SetLogy();
  TString ytitle = TString::Format("d#sigma/dp_{T} #times  #Beta(#mu#mu)  (nb/(GeV/c))",ipeak_);
  TString xtitle("p_{T}^{#Upsilon} (GeV/c)");
  mgr->GetYaxis()->SetTitleOffset(1.0); //needed
  mgr->GetYaxis()->SetTitle(ytitle);
  mgr->GetXaxis()->SetTitle(xtitle);

  mgr->SetMaximum(ymax_);
  mgr->SetMinimum(ymin_);
  if(diffonly)
    mgr->GetXaxis()->SetRangeUser(0,ptmax_);
  else
    mgr->GetXaxis()->SetRangeUser(-2,ptmax_);

  TPaveText *tp = new TPaveText(0.55,0.65,0.9,0.85,"brNDC");
  //TPaveText *tp = new TPaveText(0.55,0.75,0.9,0.85,"brNDC");
  //tp->SetFillColor(0);
  tp->SetBorderSize(0);
  //tp->SetTextAlign(22);
  //tp->SetTextColor(12);
  //tp->SetTextSize(0.04);
  tp->AddText("CMS,  #sqrt{s} = 7 TeV");
  //tp->AddText("L = 3 pb^{-1}, |y^{#Upsilon}|<2");
  int nbin;
  double* rapedge = getrapbinedges("nominal",nbin);
  TString lablumi = TString::Format("L = %2.0f pb^{-1}", lumiforlabel);
  TString labrap = TString::Format("|y^{#Upsilon}|<%2.0f",rapedge[irap_+1]);
  if(rapedge[irap_]!=0) labrap = TString::Format("%2.0f<",rapedge[irap_]) + labrap;
  tp->AddText(TString::Format("%s, %s",lablumi.Data(),labrap.Data())); 
  tp->AddText("");
  tp->AddText(TString::Format("#Upsilon(%dS)",ipeak_)); 
  tp->Draw();
  

  TString pname = TString::Format("mode%d/xsection_%ds_y%d", anamode_,ipeak_,irap_);
  a.SaveAs(pname + ".gif");
  a.SaveAs(pname + ".pdf");

  //keep order: tot, totbutlumi (denoted syst), sta 
  mgr->Add(xsecgrsys);
  mgr->Add(xsecgrsta);

  TFile o (pname + ".root","recreate");  
  //a.SaveAs(TString::Format("mode%d/xsection_%ds_y%d.gif", anamode_,ipeak_,irap_));
  //a.SaveAs(TString::Format("mode%d/xsection_%ds_y%d.pdf", anamode_,ipeak_,irap_));
  //TFile o(TString::Format("mode%d/xsection_%ds_y%d.root", anamode_,ipeak_,irap_),"recreate");  
  mgr->Write("mg1");  
  o.Close();

  TCanvas c;
  TH1F *hh = gPad->DrawFrame(0,ymin_,ptmax_,ymax_);
  gPad->SetLogy();
  hh->GetYaxis()->SetTitleOffset(1.0); //needed
  hh->SetXTitle(xtitle);
  hh->SetYTitle(ytitle);
  tp->Draw();

  double dummyx[1] = {0.}, dummyy[1] = {0.};
  dummy = new TGraph(1, dummyx, dummyy);

  TLegend *leg = new TLegend(0.15,0.2,0.55,0.4);	
  leg->AddEntry(xsecgr,   "measurement","pl");
  leg->AddEntry(xsecgrsta,"statistical unc.","p");
  leg->AddEntry(xsecgrsys,"stat. #oplus syst. unc.","p");
  leg->AddEntry(dummy,"global lumi unc. (11%)","e0");
  leg->Draw();

  xsecgrsta->SetLineColor  (38);
  xsecgrsta->SetMarkerColor(38);
  xsecgrsta->SetFillColor  (38);
  xsecgrsta->SetMarkerStyle(21);
  xsecgrsta->SetMarkerSize(1);
  //xsecgrsta->SetFillStyle(0);

  xsecgrsys->SetLineColor  (kRed);
  xsecgrsys->SetMarkerColor(kRed);
  xsecgrsys->SetMarkerStyle(25);
  xsecgrsys->SetMarkerSize(0.8);
  xsecgrsys->SetFillStyle(0);

  //xsecgrsys->SetLineWidth(2);
  //xsecgr   ->SetLineWidth(2);

  xsecgr   ->Draw("e0");
  xsecgrsta->Draw("e2");
  //xsecgrsys->Draw("e2");
  xsecgrsys->Draw("[]");

  c.SaveAs(pname + "_detailed.gif");
  c.SaveAs(pname + "_detailed.pdf");

  return;
}
