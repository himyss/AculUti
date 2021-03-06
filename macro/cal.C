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

	TH1F *sq0 = new TH1F("sq0","signals of the 0 strip",400,0,900);

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

	Double_t x[4];
	Double_t y[4];
	x[0] = 626.707305;
	x[1] = 725.447846;
	x[2] = 794.557693;
	x[3] = 60.8185;

cout  << 648.432473 << " " << g1->GetParameter(1) << endl; 
cout  << 753.847543 << " " << g2->GetParameter(1) << endl; 
cout  << 823.561856 << " " << g3->GetParameter(1) << endl; 

	//y[0] = 4.87064;
	///y[1] = 5.59033;
	//y[2] = 6.11469;
	//y[3] = 0;

	y[0] = 4.55772;
	y[1] = 5.29633;
	y[2] = 5.83365;
	y[3] = 0;

	TCanvas *c1 = new TCanvas();
	TGraph *lin = new TGraph(4,x,y);
	TF1 *fit = new TF1("fit","[0]*x + [1]");
	lin->Draw("A*");
	lin->Fit("fit");
/*
	TCanvas *c2 = new TCanvas();
	c2->cd();
	fit->SetRange(-100,900);
	fit->Draw();
	



	Double_t mean1,mean2,amp1,amp2,k,b;
	mean1 = g1->GetParameter("Mean");		
	mean2 = g2->GetParameter("Mean");	
	amp1 = 1173.2;
	amp2 = 1332.5;

	k = (amp1-amp2)/(mean1-mean2); // collibration params
	b = amp1 - k*mean1;

	// from channels to energies
	Double_t eMin,eMax,xMin,xMax;	
	xMin = h->GetXaxis()->GetXmin(); 
	xMax = h->GetXaxis()->GetXmax();

	eMin = xMin*k + b;
	eMax = xMax*k + b;

	TH1F *he = new TH1F("he","energy dist",h->GetNbinsX(),eMin,eMax);

	for(Int_t i = 0; i < h->GetNbinsX(); i++) {
		he->SetBinContent(i,h->GetBinContent(i));
	}
	
	TF1 *ge1 = new TF1("ge1","gaus",r1*k+b,r2*k+b); // range of fit
	ge1->SetParLimits(1,mrl1*k+b,mrh1*k+b);
  	TF1 *ge2 = new TF1("ge2","gaus",r1*k+b,r2*k+b); // range of fit
	ge2->SetParLimits(1,mrl2*k+b,mrh2*k+b);
	he->Fit("ge1","R");
	he->Fit("ge2","R+");

/// results
	Double_t res1,res2;
	res1 = ge1->GetParameter("Sigma")*2.355*100/ge1->GetParameter("Mean");
	cout<< "this is resolution of first peak 1173 keV: " << res1 << endl;
	res2 = ge2->GetParameter("Sigma")*2.355*100/ge2->GetParameter("Mean");
	cout<< "this is resolution of second peak 1333 keV: " << res2 << endl;
*/
return;
}
