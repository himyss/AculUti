#include <string.h>
#include <iostream.h>
#include <fstream.h>

using namespace std;

void SQXL_calib(){ 
				
	string 	input_file_name_cs = "/media/users_NAS/Muzalevsky/si_1000_LR_02_0001.root";		
	Int_t 	max_channel = 300;

	TTree* 		t;
	UShort_t    NeEvent_SQ20[16];
	UShort_t    NeEvent_SQX_L[32];
	TBranch    *b_NeEvent_SQ20;
	TBranch    *b_NeEvent_SQX_L;

	Float_t parSQ1[32], parSQ2[32];
  TString line1;
  ifstream myfile1;
  Int_t count=0;
  myfile1.open("/home/muzalevsky/AculUtils/macro/sqx_l_ec.clb");
  while (! myfile1.eof() ){
    count++;
		line1.ReadLine(myfile1);
    if(count<3) continue;
		if(line1.IsNull()) break;
    sscanf(line1.Data(),"%g %g", parSQ1+count-3,parSQ2+count-3);
    //cout << parSQ1[count-3] << " " << parSQ2[count-3] << endl;
  }
  
	/*parSQ1[0] = -0.3762; 	parSQ2[0] = 0.02216;
	parSQ1[1] = -0.35;	parSQ2[1] = 0.02198;		
	parSQ1[2] = -0.2505;	parSQ2[2] = 0.02134;		
	parSQ1[3] = -0.254;	parSQ2[3] = 0.02215;		
	parSQ1[4] = -0.2537;	parSQ2[4] = 0.02153;		
	parSQ1[5] = -0.2664;	parSQ2[5] = 0.02178;		
	parSQ1[6] = -0.5497;	parSQ2[6] = 0.02013;		
	parSQ1[7] = -0.2586;	parSQ2[7] = 0.02241;		
	parSQ1[8] = -0.2682;	parSQ2[8] = 0.02255;		
	parSQ1[9] = -0.2441;	parSQ2[9] = 0.02177;		
	parSQ1[10] = -0.2519;	parSQ2[10] = 0.0218;		
	parSQ1[11] = -0.2714;	parSQ2[11] = 0.02245;		
	parSQ1[12] = -0.2769;	parSQ2[12] = 0.0223;		
	parSQ1[13] = -0.2572;	parSQ2[13] = 0.0218;		
	parSQ1[14] = -0.2795;	parSQ2[14] = 0.02127;		
	parSQ1[15] = -0.2625;	parSQ2[15] = 0.02182;		
	parSQ1[16] = -0.3148;	parSQ2[16] = 0.02145;		
	parSQ1[17] = -0.3162;	parSQ2[17] = 0.02203;		
	parSQ1[18] = -0.3325;	parSQ2[18] = 0.02177;		
	parSQ1[19] = -0.3102;	parSQ2[19] = 0.02249;		
	parSQ1[20] = -0.3252;	parSQ2[20] = 0.02202;		
	parSQ1[21] = -0.3111;	parSQ2[21] = 0.02191;		
	parSQ1[22] = -0.3158;	parSQ2[22] = 0.02141;		
	parSQ1[23] = -0.2989;	parSQ2[23] = 0.02181;		
	parSQ1[24] = -0.3175;	parSQ2[24] = 0.02124;		
	parSQ1[25] = -0.32;	parSQ2[25] = 0.02203;		
	parSQ1[26] = -0.3286;	parSQ2[26] = 0.02192;		
	parSQ1[27] = -0.3109;	parSQ2[27] = 0.02223;		
	parSQ1[28] = -0.3135;	parSQ2[28] = 0.02134;		
	parSQ1[29] = -0.3222;	parSQ2[29] = 0.0222;		
	parSQ1[30] = -0.3256;	parSQ2[30] = 0.02148;		
	parSQ1[31] = -0.3173;	parSQ2[31] = 0.02251;		*/

	f = new TFile(input_file_name_cs.c_str());
	f->GetObject("AnalysisxTree",t);
	t->SetMakeClass(1);
	//TTree *t = (TTree*)f->Get("AnalysisxTree");
	t->SetBranchAddress("NeEvent.SQ20[16]", NeEvent_SQ20, &b_NeEvent_SQ20);
	t->SetBranchAddress("NeEvent.SQX_L[32]", NeEvent_SQX_L, &b_NeEvent_SQX_L);
	const Long64_t nentries1 = t->GetEntries();
  cout << nentries1 << endl;

	//cout<<">>> Making histograms"<<endl;
	TH1I* hist20[16];
	TString histName;
	TH1I *htest = new TH1I("htest", "htest", 300, 0, 1000);
	Int_t raw20[16];	
	Float_t SQX_L[32];

	for(Int_t i=0; i<16;i++){
      		histName.Form("hist20_%d", i);
		hist20[i] = new TH1I(histName.Data(), histName.Data(), max_channel, 0, 1000);
	}
	
	TFile *fw = new TFile("raw.root", "RECREATE");
	TTree *tw = new TTree("raw", "raw hists");
	//tw->Branch("htest","TH1I",&htest,32000,0);
	//tw->Branch("raw20",&raw20,"raw20/I");
	tw->Branch("SQX_L",&SQX_L,"SQX_L[32]/F");

	/*for(Int_t i=0; i<16;i++){
      		histName.Form("hist20_%d", i);
		cout << histName.Data() << endl;
		tw->Branch(histName.Data(),"TH1I",hist20[i],32000,0);
	}*/
  Int_t maxE = 50000;
  cout<<">>> filling TREE up to "<<maxE<< " event"<<endl;
	//for (Long64_t jentry=0; jentry<nentries1;jentry++) {
	for (Long64_t jentry=0; jentry<maxE;jentry++) {
		t->GetEntry(jentry);
		for(Int_t i=0; i<32; i++){
			//hist20[i]->Fill(NeEvent_SQ20[i]);
			//htest->Fill(NeEvent_SQ20[0]);
			//raw20[i] = NeEvent_SQ20[i];
			SQX_L[i] = NeEvent_SQX_L[i];
		}
		tw->Fill();			
	}
	
	//htest->Draw();
	tw->Write();
	fw->Close();
/*
  TH1F *hpx    = new TH1F("hpx","This is the px distribution",100,-4,4);
   TH2F *hpxpy  = new TH2F("hpxpy","py vs px",40,-4,4,40,-4,4);
   TProfile *hprof  = new TProfile("hprof","Profile of pz versus px",100,-4,4,0,20);
   T->Branch("hpx","TH1F",&hpx,32000,0);
   T->Branch("hpxpy","TH2F",&hpxpy,32000,0)
*/
	//hist20[0]->Draw();
}
