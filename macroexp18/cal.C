void cal() {
	Double_t r1,r2,mrl1,mrl2,mrh1,mrh2,mrl3,mrh3;
// choosing file
//night
	TFile *f = new TFile("/home/muzalevsky/work/calib/si_20_04.root");
	TTree *t = (TTree*)f->Get("AnalysisxTree");	
	r1 = 500; r2 = 1000; // its a gange of fitting
	mrl1 = 637.5; mrh1 = 662.5; // limits for means of gauses
	mrl2 = 742.5; mrh2 = 762.5; 
	mrl3 = 812.5; mrh3 = 837.5; 


//umamplifight.root file
/*	TFile *f = new TFile("unamplified_Co.root");	
	r1 = 510; r2 = 605; 
	mrl1 = 515; mrh1 = 535; 
	mrl2 = 575; mrh2 = 595; */
//////

	//cout<<h->GetNbinsX()<<endl;

	/*TH1F *sq0 = new TH1F("sq0","signals of the 0 strip",400,0,900);

 	TF1 *g1 = new TF1("g1","gaus",mrl1,mrh1); 
	g1->SetParLimits(1,mrl1,mrh1);
	g1->SetParLimits(2,1.,15.);
  	TF1 *g2 = new TF1("g2","gaus",mrl2,mrh2); // range of fit
	g2->SetParLimits(1,mrl2,mrh2);
	g2->SetParLimits(2,1.,15.);
  	TF1 *g3 = new TF1("g3","gaus",mrl3,mrh3); // range of fit
	g3->SetParLimits(1,mrl3,mrh3);
	g3->SetParLimits(2,1.,15.);

	t->Draw("NeEvent.SQ20[4] >> sq0","","");

	sq0->Fit("g1","R");
	sq0->Fit("g2","R+");
	sq0->Fit("g3","R+");
*/
	Double_t x[4];
	Double_t y[4],y1[4];
	x[0] = 206.225924;
	x[1] = 240.106738;
	x[2] = 264.471587;
	x[3] = 344.189266;


/*cout  << 648.432473 << " " << g1->GetParameter(1) << endl; 
cout  << 753.847543 << " " << g2->GetParameter(1) << endl; 
cout  << 823.561856 << " " << g3->GetParameter(1) << endl; 

	//y[0] = 4.87064;
	///y[1] = 5.59033;
	//y[2] = 6.11469;
	//y[3] = 0;

	y[0] = 4.09804;
	y[1] = 4.86926;
	y[2] = 5.42563;
	y[3] = 7.06078;*/

  y1[0] = 4.49274;	//		E1	//in MeV including dead layer of alpha source//энергии альфа источника
  y1[1] = 5.21546;// 				E2	//in MeV//
  y1[2] = 5.74011;//				E3	//in MeV//
  y1[3] = 7.46059;//	

//-0.202695 0.0235211

	TCanvas *c1 = new TCanvas("c1","SQL_X[17]");
  /*c1->Divide(2,1);
  c1->cd(1);
	TGraph *lin = new TGraph(4,x,y);
	TF1 *fit = new TF1("fit","[0]*x + [1]");
	lin->Draw("A*");
	lin->Fit("fit");*/

 // c1->cd(2);
	TGraph *lin1 = new TGraph(4,x,y1);
  TF1 *fit1 = new TF1("fit1","[0]*x + [1]");
	TF1 *fit2 = new TF1("fit2", "pol1", 0, 1000);
	lin1->Draw("A*");
	lin1->Fit("fit2");
 
return;
}
