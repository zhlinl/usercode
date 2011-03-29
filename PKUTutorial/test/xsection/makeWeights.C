#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TString.h"
#include "TTree.h"
#include "TLorentzVector.h"


/*
function:
- computes, appends weight column to inout three
notes:
- missing: tracking efficiency
- note double vs single muon trigger
*/

Float_t get(TH2F* h, Float_t pt, Float_t eta);
bool passMuCut(double muEta, double muPt);

TString accHistoName     = "acc_pt35_eta16_pt25_eta24_trk_unpol"; //(pt,y)
TString effTrckHistoName = "h_hits2d"; //(eta,pt)
TString effTrigHistoName = "TH2F_tagAndProbe_Trigger_TM_pt_abseta_fit_eff"; //(pt,eta)
TString effMuidHistoName = "TH2F_tagAndProbe_TM_pt_abseta_fit_eff"; //(Pt,eta)

TString treeDir = "yieldUpsilonTree";
TString treeName = "probe_tree";

TString inFile  ="data.root";
TString accFile ="acceptance.root";
TString trkFile ="effMap_track.root";
TString idFile  ="effMaps_JPsi_DATA_MuonID.root";
TString trigFile="effMaps_JPsi_DATA_Trigger.root";		 
TString outFile ="upsilonYieldWeighted.root";
TString mapDir ="maps/";

void makeWeights(
		 int dosys   = 0, 
		 bool sys_hi = 0,
		 int ipeak   = 0, //0,1,2 = 1s/2s/3s
		 int pol     = 0
		 ) {

  
  ///                                  1         2       3          4          5          6         7        8        9        10      11       12      13        14          15
  enum systType           {nominal=0, AccSta , EtrkSta , EmuidSta , EtrigSta , EtrecoSta, ptscale , ptreso , ptspec , vtxpos , nofsr , tnpmc , mctrue , tnpmcUps , linear812 , other , thelast};
  const int nsys = thelast;
  TString sysN[nsys] = {"nominal","Acc",   "Etrk",   "Emuid",   "Etrig",   "Etreco",  "ptscale","ptreso","ptspec","vtxpos","nofsr","tnpmc","mctrue","tnpmcUps","linear812","other"        };
  bool syst[nsys] = {false};
  
  //vertex probability (using jpsi)
  double VtxProbEff   = 0.9916;
  double VtxProbEff_e = 0.0009;
  //track quality (using tnp)
  double trkQualEff   = 0.9866;
  double trkQualEff_e = 0.0005;
  //highest vertex prob candidate poer event (from Upsilon data)
  double bestCandEff   = 0.998;
  double bestCandEff_e = 0.002;
  //mass scale correction
  double ms_a0   = 3.8e-4;
  double ms_a0_e = 1.9e-4;
  double ms_a1   = 0.;
  double ms_a1_e = 0.;
  double ms_a2   = 3.0e-4;
  double ms_a2_e = 0.7e-4;
  double ms_a3   = 0.;
  double ms_a3_e = 0.;
  
  TString polN[5] = {"unpol","helT","helL","csT","csL"      }; 

  outFile.ReplaceAll(".root","");

  if(pol) {
    accHistoName.ReplaceAll("unpol",polN[pol]);
    sysN[0].Append(TString::Format("_%s",polN[pol].Data()));
  }

  if(ipeak) {
    accHistoName.Append(TString::Format("_%ds",ipeak+1)); 
    accFile.ReplaceAll(".root",TString::Format("_%ds.root",ipeak+1)); 
    sysN[0].Append(TString::Format("_%ds",ipeak+1)); 
  }
  
  syst[dosys] = true;

  outFile.Append("_"+sysN[dosys]);
  //for now systematics computed only from 1S, unpolarized
  if(ipeak || pol) { 
    assert(dosys==0);
    outFile.ReplaceAll("_nominal","");
  }


  
  if(syst[AccSta]) { //Acc sta
    if(sys_hi) {
      accFile          = "low_acceptance.root";
      outFile.Append("Hi");
    } else {
      accFile          = "hi_acceptance.root";
      outFile.Append("Lo");
    }
  } else if (syst[EtrkSta]) {//eff track sta: include track quality, vertex prob
    if(sys_hi) {
      effTrckHistoName = "h_effminus";
      VtxProbEff      -= VtxProbEff_e;
      trkQualEff      -= trkQualEff_e;
      bestCandEff     -= bestCandEff_e;
      outFile.Append("Hi");
    } else {
      effTrckHistoName = "h_effplus";
      VtxProbEff      += VtxProbEff_e;
      trkQualEff      += trkQualEff_e;
      bestCandEff     += bestCandEff_e;
      outFile.Append("Lo");
    }
  } else if (syst[EmuidSta]) { //eff muid sta
    if(sys_hi) {effMuidHistoName.ReplaceAll("TH2F","TH2F_low"  ); outFile.Append("Hi");} 
    else       {effMuidHistoName.ReplaceAll("TH2F","TH2F_high" ); outFile.Append("Lo");}
  } else if (syst[EtrigSta]) { //eff trig sta
    if(sys_hi) {effTrigHistoName.ReplaceAll("TH2F","TH2F_low"  ); outFile.Append("Hi");} 
    else       {effTrigHistoName.ReplaceAll("TH2F","TH2F_high" ); outFile.Append("Lo");}
  } else if (syst[EtrecoSta]) { //eff muid+trig sta (coherent/conservative)
    if(sys_hi) {effTrigHistoName.ReplaceAll("TH2F","TH2F_low"  ); 
                effMuidHistoName.ReplaceAll("TH2F","TH2F_low"  ); outFile.Append("Hi");} 
    else       {effTrigHistoName.ReplaceAll("TH2F","TH2F_high" ); 
                effMuidHistoName.ReplaceAll("TH2F","TH2F_high" ); outFile.Append("Lo");}
  } else if(syst[ptscale]) { //Acc pt bias
    accFile = "acceptance_syst_unpol.root";
    if (sys_hi) {accHistoName+= "_ptup";      outFile.Append("Hi");}
    else        {accHistoName+= "_ptdown";    outFile.Append("Lo");}
  } else if(syst[ptreso]) {
    accFile = "acceptance_syst_unpol.root";
    if(sys_hi)  {accHistoName+= "_sigmaup";   outFile.Append("Hi");}
    else        {accHistoName+= "_sigmadown"; outFile.Append("Lo");}
  } else if(syst[ptspec]) { // Acc PTspectrum
    accFile = "acceptance_syst_ptspec.root";
    accHistoName = "acc_pt35_eta16_pt25_eta24_trk_unpol_pt";// outFile.Append("_ptspec");
  } else if(syst[vtxpos]) { // Acc Vertex
    accFile = "acceptance_syst_vtxpos.root";
    accHistoName = "acc_pt35_eta16_pt25_eta24_trk_unpol_vtx";// outFile.Append("_vtxpos");
  } else if(syst[nofsr]) { // Acc FSR
    accFile = "acceptance_syst_nofsr.root";
    accHistoName += "_nofsr";
  } else if(syst[tnpmc]) { // Tnp MC
    idFile  = "effMaps_JPsi_MC_MuonID.root ";
    trigFile= "effMaps_JPsi_MC_Trigger.root";
  } else if(syst[mctrue]) { // tnp mc truth (jpsi)
    idFile  = "effMaps_JPsi_MC_MuonID_MCTruth.root";
    trigFile= "effMaps_JPsi_MC_Trigger_MCTruth.root";
    effMuidHistoName.ReplaceAll("fit_eff","mcTrue_cnt_eff");
    effTrigHistoName.ReplaceAll("fit_eff","mcTrue_cnt_eff");
    //effMuidHistoName = "Th2f_histoMuFromTk_TMI_pt_eta_mcTrue_cnt_eff";
    //effTrigHistoName = "TH2F_high_histoTrigger_L1DoubleMuOpenpt_eta_mcTrue_cnt_eff";
  } else if(syst[tnpmcUps]) {
    idFile  = "effMaps_AllUpsilons_MC_MuonID.root";
    trigFile= "effMaps_AllUpsilons_MC_Trigger.root";
    effMuidHistoName.ReplaceAll("fit_eff","mcTrue_cnt_eff");
    effTrigHistoName.ReplaceAll("fit_eff","mcTrue_cnt_eff");
  } else if(syst[other]) {
    //mass mpm
    if(sys_hi) { ms_a0  -= ms_a0_e;  ms_a1 -= ms_a1_e;  ms_a2 -= ms_a2_e;  ms_a3   -= ms_a3_e; outFile.Append("Hi");}
    else       { ms_a0  += ms_a0_e;  ms_a1 += ms_a1_e;  ms_a2 += ms_a2_e;  ms_a3   += ms_a3_e; outFile.Append("Lo");}
  }


  outFile.Append(".root");

  //printf("id:%s\n",idFile.Data());

  accFile  = mapDir + accFile;
  trigFile = mapDir + trigFile;
  idFile   = mapDir + idFile;
  trkFile  = mapDir + trkFile;

  cout << "\tpeak:" << ipeak +1 << "S" 
       << "\tpol:"  << polN[pol]
       << "\n\tsystematic: " << sysN[dosys]
       << "\n\tproduced tree: " << outFile
       << "\n\taccFile :" << accFile  << "\n\t\taccHistoName    :" << accHistoName     
       << "\n\ttrigFile:" << trigFile << "\n\t\teffTrigHistoName:" << effTrigHistoName 
       << "\n\tidFile  :" << idFile   << "\n\t\teffMuidHistoName:" << effMuidHistoName 
       << "\n\ttrkFile :" << trkFile  << "\n\t\teffTrckHistoName:" << effTrckHistoName 
       << endl;

  // open file with the yield
  TFile infile(inFile);
  gDirectory->Cd(treeDir);
  TTree* inTree = (TTree*)gROOT->FindObject(treeName);
  if(!inTree){
    cout<<"Could not access yield tree!"<<endl;
    return;
  }

  Float_t invariantMass, upsPt, upsRapidity, muPlusPt, muPlusEta, muMinusPt, muMinusEta, muMinusPhi, muPlusPhi, muPlusCharge;
  //Int_t trigger;
  inTree->SetBranchAddress("invariantMass",&invariantMass);
  inTree->SetBranchAddress("upsPt"        ,&upsPt);
  inTree->SetBranchAddress("upsRapidity"  ,&upsRapidity);
  inTree->SetBranchAddress("muPlusPt"     ,&muPlusPt);
  inTree->SetBranchAddress("muPlusEta"    ,&muPlusEta);
  inTree->SetBranchAddress("muMinusPt"    ,&muMinusPt);
  inTree->SetBranchAddress("muMinusEta"   ,&muMinusEta);
  inTree->SetBranchAddress("muPlusPhi"    ,&muPlusPhi);
  inTree->SetBranchAddress("muMinusPhi"   ,&muMinusPhi);
  inTree->SetBranchAddress("muPlusCharge" ,&muPlusCharge);

  // open file with acceptance
  std::cout<<"step 1"<<std::endl;

  TFile accfile(accFile);
  TH2F* acc = (TH2F*)gROOT->FindObject(accHistoName);
  if(!acc){
    cout<<"Could not access acceptance histogram!"<<endl;
    return;
  }
  // open file with id efficiency
  TFile idfile(idFile);
  TH2F* id = (TH2F*)gROOT->FindObject(effMuidHistoName);
  if(!id){
    cout<<"Could not access id efficiency histogram!"<<endl;
    return;
  }
  // open file with trigger efficiency
  TFile trigfile(trigFile);
  TH2F* trig = (TH2F*)gROOT->FindObject(effTrigHistoName);
  if(!trig){
    cout<<"Could not access trigger efficiency histogram!"<<endl;
    return;
  }
  // open file with trigger efficiency
  TFile trkfile(trkFile);
  TH2F* trk = (TH2F*)gROOT->FindObject(effTrckHistoName);
  if(!trk){
    cout<<"Could not access tracking efficiency histogram!"<<endl;
    return;
  }
  
  // create file for the yield and weights
  TFile outfile(outFile,"recreate");
  //gDirectory->mkdir("upsilonYield")->cd();

  TTree *outTree = inTree->CloneTree(0);
  //TTree* outTree = (TTree*)inTree->CopyTree("8 < invariantMass < 12");

  Float_t weight, weightTrig, weightAcc, weightMuid, weightTrk;

  //outTree->Show();
  outTree->Branch("weight",     &weight,     "weight/F");
  outTree->Branch("weightTrig", &weightTrig, "weightTrig/F");
  outTree->Branch("weightAcc",  &weightAcc,  "weightAcc/F");
  outTree->Branch("weightMuid", &weightMuid, "weightMuid/F");
  outTree->Branch("weightTrk",  &weightTrk,  "weightTrk/F");

  //counters
  int cnt_all(0), cnt_mass(0), cnt_mucut(0), cnt_trig(0), cnt_wAcc(0), cnt_wTrg(0), cnt_wTrk(0), cnt_wTag(0);
  double sum_wei(0.),sum_wacc(0.),sum_wtrg(0.), sum_wtrk(0.), sum_wtag(0.);
  double sum_eff(0.),sum_acc (0.),sum_trg (0.), sum_trk (0.), sum_tag (0.);
  int cnt_cowboys(0);
  int passEvt(0);

  bool pass = true;
  
  //loop through the candidates, filter, and calculate weights
  Int_t nentries = inTree->GetEntries();
  cout << "nentries:" << nentries << endl;

  for (Int_t i=0;i<nentries; i++) {

    if(i%10000==0) std::cout<< "=>" << i<<std::endl;
    inTree->GetEntry(i);
        
    TString dump("");
    double wAcc(0.), wTag(0.), wTrg(0.), wTrk(0.);
    
    //note: needed to swap here
    //      as acceptance and efficiency histograms have different axis order
    wAcc = get(acc ,upsPt   ,fabs(upsRapidity));
    wTag = get(id  ,muPlusPt,fabs(muPlusEta))  * get(id  ,muMinusPt,fabs(muMinusEta));
    wTrg = get(trig,muPlusPt,fabs(muPlusEta))  * get(trig,muMinusPt,fabs(muMinusEta));
    wTrk = get(trk ,fabs(muPlusEta), muPlusPt) * get(trk ,fabs(muMinusEta), muMinusPt);
    //wTrg = ( 1. - (1.-get(trig,muPlusPt,fabs(muPlusEta))) * (1.-get(trig,muMinusPt,fabs(muMinusEta))) );
    //    printf("== %5.2f  %5.2f  %5.2f  %5.2f \n", 	   wAcc,wTag,wTrk,wTrg);

    cnt_all++;
    pass = true;
    
    if(invariantMass > 8 && invariantMass < 14)
      cnt_mass++;
    else 
      pass &= false;
    
    if(passMuCut(muPlusEta,muPlusPt) && passMuCut(muMinusEta,muMinusPt))
      cnt_mucut++;
    else 
      pass &= false;
    
    pass &= (fabs(upsRapidity)<2.0);

    pass &= (fabs(upsPt)<30.0); //impose this cut for consistency throughout downstream


    //mass correction
    TLorentzVector m1, m2, ups;
    m1.SetPtEtaPhiM((1.0 + ms_a0 + ms_a1*fabs(muPlusEta ) + ms_a2*muPlusEta *muPlusEta  + ms_a3*muPlusPt )*muPlusPt ,muPlusEta ,muPlusPhi , 0.105658);
    m2.SetPtEtaPhiM((1.0 + ms_a0 + ms_a1*fabs(muMinusEta) + ms_a2*muMinusEta*muMinusEta + ms_a3*muMinusPt)*muMinusPt,muMinusEta,muMinusPhi, 0.105658);
    ups = m1+m2;
    invariantMass = ups.M();
    upsPt = ups.Pt();
    //upsRapidity = ups.Rapidity(); //do not correct rapidity as we already cut on it


    if(!pass) continue;
    
    if(!wAcc)
      pass &= false;
    else cnt_wAcc++; 
    if(!wTag)
      pass &= false;
    else cnt_wTag++; 
    if(!wTrk)
      pass &= false;
    else cnt_wTrk++; 
    if(!wTrg)
      pass &= false;
    else cnt_wTrg++;
    
    double wTrkQualEff = 1./pow(trkQualEff,2);
    double wUpsGlbEff  = 1./VtxProbEff/bestCandEff;
    
    weight = wAcc * wTrk * wTag * wTrg * wTrkQualEff * wUpsGlbEff;
    

    //double bend = muPlusCharge*(muPlusPhi-muMinusPhi);
    //bool cowboy = (bend<-3.1415 | (bend>0 & bend<3.1415)); 
    //if(!cowboy) pass &= false; else cnt_cowboy++;

    dump = TString::Format("m:%5.2f upspt:%5.2f upsrap:%5.2f pt+:%5.2f pt-:%5.2f eta+:%5.2f eta-:%5.2f w:%5.2f (acc:%5.2f trk:%5.2f muid:%5.2f trg:%5.2f upsGlb:%5.2f)",
			   invariantMass,upsPt, upsRapidity, muPlusPt, muMinusPt, muPlusEta, muMinusEta,weight,wAcc, wTrk, wTag, wTrg, wUpsGlbEff);
			   
    if(!weight)
      cout << "ZERO WEI: " << dump << endl;


    //    if(upsPt>20 && fabs(invariantMass - 9.43)<0.3) printf("upspt: %5.3f mass: %5.3f\n", upsPt,invariantMass);


    if(!pass) cout << dump << endl;
    
    //    pass = true; weight=1; wAcc=1; wTrg=1; wTag=1; wTrk=1;
    
    if(!pass)
      continue;

    if(!weight) 
      cout << "skip event "<< cnt_all << " (w="<<weight<<") " << dump << " wacc:" << wAcc << " wtag:" << wTag << " wtrig:" << wTrg << "\n";
    
    assert(weight!=0);
 
   // weights
    weight     = 1./weight;
    weightAcc  = 1./wAcc;
    weightTrig = 1./wTrg;
    weightMuid = 1./wTag;
    weightTrk  = 1./wTrk;

    // average weights
    sum_wei  += weight;
    sum_wacc += weightAcc ;
    sum_wtrg += weightTrig;
    sum_wtag += weightMuid;
    sum_wtrk += weightTrk ;

    // average efficiencies
    sum_eff += 1./weight;
    sum_acc += wAcc;
    sum_trg += wTrg;
    sum_tag += wTag;
    sum_trk += wTrk;
    
    passEvt++;

    //weight = 1.;
    outTree->Fill();
  }


  // TBD: check for events with mutiple candidates!
  
  // print stats
  printf("\nSTATS|\tall:%d\tmass:%d\tmucut:%d\ttrig:%d\twAcc:%d\twTag:%d\twTrk:%d\twTrg:%d [cowboys:%d]\n",
	 cnt_all, cnt_mass, cnt_mucut, cnt_trig, 
	 cnt_wAcc, cnt_wTag, cnt_wTrk, cnt_wTrg, 
	 cnt_cowboys);

  printf("\tsystematics:%s-%s \taverage weight: %7.4f\n",sysN[dosys].Data(),(sys_hi?"hi":"lo"),sum_wei/passEvt);

  // average weights
  sum_wei  /= passEvt;
  sum_wacc /= passEvt;
  sum_wtrg /= passEvt;
  sum_wtag /= passEvt;
  sum_wtrk /= passEvt;

  // average efficiencies
  sum_eff /= passEvt;
  sum_acc /= passEvt;
  sum_trg /= passEvt;
  sum_tag /= passEvt;
  sum_trk /= passEvt;

  printf("average all-eff:%5.2f\t acc:%5.2f\t trg:%5.2f\t muid:%5.2f\t trk:%5.2f\n",
	 sum_eff,sum_acc,sum_trg,sum_tag,sum_trk);
  printf("average all-wei:%5.2f\twacc:%5.2f\twtrg:%5.2f\twmuid:%5.2f\twtrk:%5.2f\n ",
	 sum_wei,sum_wacc,sum_wtrg,sum_wtag,sum_wtrk);
  
  cout << "Produced tree: " << outFile << endl;


  //close everything
  outTree->Write();
  outfile.Close();
  infile.Close();
  idfile.Close();
  trigfile.Close();
}

Float_t get(TH2F* h, Float_t pt, Float_t eta){
  //Int_t FindBin(Double_t x, Double_t y = 0, Double_t z = 0)
  int i = h->FindBin(pt,eta);
  return h->GetBinContent(i);
}

//check muon cuts
bool passMuCut(double muEta, double muPt) {
  // Muon pt > 3.5 for |eta|<1.6, pt>2.5 for 1.6<|eta|<2.4
  if(fabs(muEta)>2.4)
    return false;
  if(fabs(muEta)<1.6 && muPt<3.5)
    return false;
  if(muPt<2.5)
    return false;
  return true;
}
