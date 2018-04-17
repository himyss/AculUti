#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

void readTree() {

	TFile *fr = new TFile("/home/muzalevsky/AculUtils/exp1804/raw.root");
	TTree *tr = (TTree*)fr->Get("raw");

//	return;

	cout<<">>> Making histograms"<<endl;
	TH1I* hist20[16];
	TString histName;

	TH1I *htest = new TH1I("htest", "htest", 300, 0, 1000);

	/*for(Int_t i=0; i<16;i++){
      		histName.Form("hist20_%d", i);
		tr->SetBranchAddress(histName.Data(), &hist20[i]);
	}*/
	tr->SetBranchAddress("htest", &htest);
	htest->Draw();
/*

	TFile *fw = new TFile("/home/muzalevsky/AculUtils/exp1804/raw.root", "RECREATE");
	TTree *tw = new TTree("raw", "raw hists");

	Float_t SQ20[16];
	tw->Branch("SQ20",SQ20,"SQ20[16]/F");

//	const Long64_t nevents = tr->GetEntries();
	const Long64_t nevents = tr->GetEntries();

	for (Int_t i = 0; i < 15; i++) {
		tr->GetEvent(i);
//		cout << revent->SQX_L[0] << endl;
		SQ20[0] = revent->SQ20[0];
		cout << SQ20[0] << endl;
		tw->Fill();
	}

	tw->Write();
	fw->Close();
*/
}
