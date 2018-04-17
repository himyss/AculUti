draw() {

	TFile *fr = new TFile("/home/muzalevsky/AculUtils/exp1803/SQ300_full_calibrated_spectra.root", "READ");

/*	TH1F *h0 = new TH1F("h0","hist 0", 1000,.0,10.);
	TH1F *h1 = new TH1F("h1","hist 1", 1000,.0,10.);
	TH1F *h2 = new TH1F("h2","hist 2", 1000,.0,10.);
	TH1F *h3 = new TH1F("h3","hist 3", 1000,.0,10.);
	TH1F *h4 = new TH1F("h4","hist 4", 1000,.0,10.);
	TH1F *h5 = new TH1F("h5","hist 5", 1000,.0,10.);
	TH1F *h6 = new TH1F("h6","hist 6", 1000,.0,10.);
	TH1F *h7 = new TH1F("h7","hist 7", 1000,.0,10.);
	TH1F *h8 = new TH1F("h8","hist 8", 1000,.0,10.);
	TH1F *h9 = new TH1F("h9","hist 9", 1000,.0,10.);
	TH1F *h10 = new TH1F("h10","hist 10", 1000,.0,10.);
	TH1F *h11 = new TH1F("h11","hist 11", 1000,.0,10.);
	TH1F *h12 = new TH1F("h12","hist 12", 1000,.0,10.);
	TH1F *h13 = new TH1F("h13","hist 13", 1000,.0,10.);
	TH1F *h14 = new TH1F("h14","hist 14", 1000,.0,10.);
	TH1F *h15 = new TH1F("h15","hist 15", 1000,.0,10.);*/


	TH1F *h0 = (TH1F*)fr->Get("HistSQ300[0]Efull");
	TH1F *h1 = (TH1F*)fr->Get("HistSQ300[1]Efull");
	TH1F *h2 = (TH1F*)fr->Get("HistSQ300[2]Efull");
	TH1F *h3 = (TH1F*)fr->Get("HistSQ300[3]Efull");
	TH1F *h4 = (TH1F*)fr->Get("HistSQ300[4]Efull");
	TH1F *h5 = (TH1F*)fr->Get("HistSQ300[5]Efull");
	TH1F *h6 = (TH1F*)fr->Get("HistSQ300[6]Efull");
	TH1F *h7 = (TH1F*)fr->Get("HistSQ300[7]Efull");
	TH1F *h8 = (TH1F*)fr->Get("HistSQ300[8]Efull");
	TH1F *h9 = (TH1F*)fr->Get("HistSQ300[9]Efull");
	TH1F *h10 = (TH1F*)fr->Get("HistSQ300[10]Efull");
	TH1F *h11 = (TH1F*)fr->Get("HistSQ300[11]Efull");
	TH1F *h12 = (TH1F*)fr->Get("HistSQ300[12]Efull");
	TH1F *h13 = (TH1F*)fr->Get("HistSQ300[13]Efull");
	TH1F *h14 = (TH1F*)fr->Get("HistSQ300[14]Efull");
	TH1F *h15 = (TH1F*)fr->Get("HistSQ300[15]Efull");

	h0->GetXaxis()->SetRangeUser(3,7);
	h1->GetXaxis()->SetRangeUser(3,7);
	h1->SetLineColor(2);
	h0->Draw();
	h1->Draw("same");





}
