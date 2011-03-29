divide(const char * file="acceptance_basehistos_1.root",const char * outname="acceptance.root") {
  TFile* f1 = new TFile(file);


  TH2F* genPtRap1 = (TH2F*) f1->Get("genUps1");
  TH2F* genPtRap2 = (TH2F*) f1->Get("genUps2");
  TH2F* genPtRap3 = (TH2F*) f1->Get("genUps3");
  TH2F* genPtRap4 = (TH2F*) f1->Get("genUps4");
  TH2F* genPtRap5 = (TH2F*) f1->Get("genUps5");

  cout << "entries = " << genPtRap1->GetEntries() << endl;


  TH2F* recoPtRap1 = (TH2F*) f1->Get("recoUps1");
  TH2F* recoPtRap2 = (TH2F*) f1->Get("recoUps2");
  TH2F* recoPtRap3 = (TH2F*) f1->Get("recoUps3");
  TH2F* recoPtRap4 = (TH2F*) f1->Get("recoUps4");
  TH2F* recoPtRap5 = (TH2F*) f1->Get("recoUps5");


  TFile out(outname,"recreate");
  TH2F* acc1 = (TH2F*)genPtRap1->Clone("acc_pt35_eta16_pt25_eta24_trk_unpol");
  acc1->Sumw2();
  acc1->Divide(recoPtRap1,genPtRap1,1,1,"B");

  TH2F* acc2 = (TH2F*)genPtRap2->Clone("acc_pt35_eta16_pt25_eta24_trk_helT");
  acc2->Sumw2();
  acc2->Divide(recoPtRap2,genPtRap2,1,1,"B");

  TH2F* acc3 = (TH2F*)genPtRap3->Clone("acc_pt35_eta16_pt25_eta24_trk_helL");
  acc3->Sumw2();
  acc3->Divide(recoPtRap3,genPtRap3,1,1,"B");

  TH2F* acc4 = (TH2F*)genPtRap4->Clone("acc_pt35_eta16_pt25_eta24_trk_csT");
  acc4->Sumw2();
  acc4->Divide(recoPtRap4,genPtRap4,1,1,"B");

  TH2F* acc5 = (TH2F*)genPtRap5->Clone("acc_pt35_eta16_pt25_eta24_trk_csL");
  acc5->Sumw2();
  acc5->Divide(recoPtRap5,genPtRap5,1,1,"B");

  acc1->Write();
  acc2->Write();
  acc3->Write();
  acc4->Write();
  acc5->Write();

  out.Close();
}
