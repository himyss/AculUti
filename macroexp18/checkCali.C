void checkCali() {
  TFile *f = new TFile("/home/muzalevsky/work/exp1803/data/siCal/SQY_R_full_calibrated_spectra.root");
  const Int_t nhists = 16;
  TH1F *h[nhists],*hsumm;
  hsumm = new TH1F("hsumm","summ cal spec",4095,0,10);
  TString hname,cName;
  //TTree *t = (TTree*)f->Get("AnalysisxTree");
  for(Int_t i=0;i<nhists;i++) {
    hname.Form("HistSQY_R[%d]Efull",i);
    h[i] = (TH1F*)f->Get(hname.Data());
    //cout  << h[i]->GetXaxis()->GetXmax()<< " " << i<< endl;
  }
  TCanvas *c1 = new TCanvas("c1"," summ cal spec",1000,1000);
  Float_t content;
  for(Int_t i=0;i<4095;i++) {
    content = 0.;
    for(Int_t j=0;j<nhists;j++) {
      if(j==0) continue;
      content += h[j]->GetBinContent(i);    
    }
    hsumm->SetBinContent(i,content);
  }

  TF1* g1 = new TF1("g1", "gaus", 4.2, 4.8);
  g1->SetParLimits(0,1.,150.);
  g1->SetParLimits(1,4.,5.11);
  g1->SetParLimits(2,0.1,3.);

  TF1* g2 = new TF1("g2", "gaus", 5., 5.5);
  g2->SetParLimits(0,1.,150.);
  g2->SetParLimits(1,5.12,5.5);
  g2->SetParLimits(2,0.1,3.);

  TF1* g3 = new TF1("g3", "gaus", 5.65, 5.95);
  g3->SetParLimits(0,1.,150.);
  g3->SetParLimits(1,5.65,5.95);
  g3->SetParLimits(2,0.1,3.);

  TF1* g4 = new TF1("g4", "gaus", 7.25, 7.65);
  g4->SetParLimits(0,1.,150.);
  g4->SetParLimits(1,7.25,7.65);
  g4->SetParLimits(2,0.1,3.);

  c1->cd();
  hsumm->Draw();
  hsumm->Fit("g1","R+");
  hsumm->Fit("g2","R+");
  hsumm->Fit("g3","R+");
  hsumm->Fit("g4","R+");
return;
  TCanvas *c[4];
  for(Int_t i=0;i<4;i++){
    cName.Form("c%d",i+1);
    c[i] = new TCanvas(cName.Data(),"calibrated spectra",1000,1000);
    c[i]->Divide(2,2);
  }

  /*c1->cd(1);
  h[0]->Draw();
  h[0]->Rebin(8);
  h[0]->GetXaxis()->SetRangeUser(4.,8.);
  h[0]->Fit("g1","R+");
  h[0]->Fit("g2","R+");
  h[0]->Fit("g3","R+");
  h[0]->Fit("g4","R+");*/

	ofstream outcalfile;
	outcalfile.open("/home/muzalevsky/work/exp1803/data/dataCal/check_cali.txt");
	if (!outcalfile.is_open()) {
		cout <<"Output file was not opened" << endl;
		return;
	}

  Int_t count=-1;
  Int_t nPad;
  for(Int_t i=0;i<nhists;i++) {
    count = i/4;
    nPad = (i%4)+1;
    c[count]->cd(nPad);

    h[i]->Draw();
    h[i]->Rebin(8);
    h[i]->Fit("g1","R+");
    h[i]->Fit("g2","R+");
    h[i]->Fit("g3","R+");
    h[i]->Fit("g4","R+");
    h[i]->GetXaxis()->SetRangeUser(4.,8.);
    c[count]->Update();
    outcalfile << g1->GetParameter(1) << "\t" << g2->GetParameter(1) << "\t" << g3->GetParameter(1) << "\t" << g4->GetParameter(1) << endl; 
  }

  return;
}
